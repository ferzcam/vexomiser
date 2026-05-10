"""
Train INDIGENA on PAVS dataset and evaluate on the test split.

INDIGENA uses TransD KGE embeddings on the UPheno ontology graph
augmented with gene-HP/MP and disease-HP associations.
At inference, genes are scored against each test case's HPO terms
using BMA (Best-Match Average) over phenotype embeddings.

Graph options (cumulative):
  Graph 1 (baseline): UPheno OWL2VecStar projection (HP + MP + UPHENO terms)
  Graph 2: + gene → has_phenotype → MP/HP
             Primary: MGI mouse-knockout phenotypes (human ortholog → MP terms)
             Fallback: PAVS training case HPO terms for genes with no MGI ortholog
  Graph 3: + disease → has_symptom → HP (training diseases only; inductive)
  Graph 4: + gene → associated_with → disease (supervised; training only)

Usage (run from vexomiser root):
    ~/miniforge3/envs/indigena/bin/python eval/indigena_train.py \\
        --upheno-edges ~/Git/indigena/data/upheno_owl2vecstar_edges.tsv \\
        --mgi-gene-phenotypes ~/Git/indigena/data/gene_phenotypes.csv \\
        --hom-file ~/Git/indigena/data/HOM_MouseHumanSequence.rpt \\
        --graph2 --graph3 --graph4 \\
        --track 1 --eval-split test
"""

import mowl
mowl.init_jvm("10g")   # must come before any mowl/JPype imports

from mowl.utils.random import seed_everything
from mowl.projection import Edge
from pykeen.models import TransD
from pykeen.training import SLCWATrainingLoop
from pykeen.training.callbacks import TrainingCallbackHint
import torch as th
from torch.optim import Adam

import os
import logging

import click as ck
import numpy as np
import pandas as pd
from tqdm import tqdm

from metrics import compute_metrics, print_as_tex

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(REPO_ROOT, "data")
RESULTS_DIR = os.path.join(DATA_DIR, "results")
MODELS_DIR = os.path.join(DATA_DIR, "models")

METRIC_KEYS = ["mr", "mrr", "hits@1", "hits@3", "hits@10", "hits@100", "auc"]


# ---------------------------------------------------------------------------
# IRI helpers — must match UPheno OWL2VecStar projection format
# ---------------------------------------------------------------------------

def hp_iri(hp_id: str) -> str:
    """HP:0001263 → http://purl.obolibrary.org/obo/HP_0001263"""
    return "http://purl.obolibrary.org/obo/" + hp_id.replace(":", "_")


def gene_iri(symbol: str) -> str:
    return f"http://mowl.borg/gene/{symbol}"


def disease_iri(disease_id: str) -> str:
    """MONDO:0010088 → http://mowl.borg/MONDO_0010088"""
    return "http://mowl.borg/" + disease_id.replace(":", "_")


# ---------------------------------------------------------------------------
# MGI ortholog + phenotype loading
# ---------------------------------------------------------------------------

def load_human_to_mgi(hom_file: str) -> dict:
    """
    Parse HOM_MouseHumanSequence.rpt and return {human_symbol: mgi_iri}.
    Groups rows by HomoloGene ID; pairs each human symbol with the mouse MGI ID
    in the same group.
    """
    df = pd.read_csv(hom_file, sep="\t", dtype=str)
    df.columns = df.columns.str.strip()

    # Column names from the file header
    # "HomoloGene ID", "Common Organism Name", "NCBI Taxon ID", "Symbol",
    # "EntrezGene ID", "Mouse MGI ID", ...
    # Column name differs between full file ("DB Class Key") and test fixture ("HomoloGene ID")
    group_col = "DB Class Key" if "DB Class Key" in df.columns else "HomoloGene ID"
    groups = df.groupby(group_col)
    human_to_mgi = {}
    for _, grp in groups:
        mouse_rows = grp[grp["NCBI Taxon ID"] == "10090"]
        human_rows = grp[grp["NCBI Taxon ID"] == "9606"]
        if mouse_rows.empty or human_rows.empty:
            continue
        mgi_id = mouse_rows.iloc[0]["Mouse MGI ID"]  # e.g. "MGI:87867"
        if not isinstance(mgi_id, str) or not mgi_id.startswith("MGI:"):
            continue
        mgi_iri = "http://mowl.borg/" + mgi_id.replace(":", "_")  # mowl.borg/MGI_87867
        for human_symbol in human_rows["Symbol"]:
            human_to_mgi[str(human_symbol)] = mgi_iri
    return human_to_mgi


def load_mgi_gene_phenotypes(mgi_pheno_csv: str, graph1_entities: set) -> dict:
    """
    Load pre-computed MGI gene→MP associations (indigena/data/gene_phenotypes.csv).
    Returns {mgi_iri: [mp_iri, ...]} filtered to phenotypes in the UPheno graph.
    """
    df = pd.read_csv(mgi_pheno_csv)
    mgi_to_phenos: dict[str, set] = {}
    for _, row in df.iterrows():
        gene = str(row["Gene"])    # e.g. http://mowl.borg/MGI_99604
        pheno = str(row["Phenotype"])  # e.g. http://purl.obolibrary.org/obo/MP_0000633
        if pheno not in graph1_entities:
            continue
        mgi_to_phenos.setdefault(gene, set()).add(pheno)
    return {g: sorted(p) for g, p in mgi_to_phenos.items()}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def normalize_disease_id(row) -> str:
    """Return MONDO ID if present, else a synthetic ID unique to the case."""
    did = row.get("disease_id", "")
    if pd.notna(did) and str(did).startswith("MONDO:"):
        return str(did)
    return f"SYNTHETIC:{row['case_id']}"


def load_cases(track: int) -> pd.DataFrame:
    path = os.path.join(DATA_DIR, f"track{track}_cases.tsv")
    df = pd.read_csv(path, sep="\t")
    df["hpo_list"] = df["hpo_terms"].apply(
        lambda s: [t.split("|")[0] for t in str(s).split(";")]
        if pd.notna(s) else []
    )
    df["disease_id_norm"] = df.apply(normalize_disease_id, axis=1)
    return df


def load_split_ids(split: str) -> set:
    path = os.path.join(DATA_DIR, "splits", f"{split}.tsv")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Split file not found: {path}. "
            "Run: uv run python eval/generate_splits.py"
        )
    return set(pd.read_csv(path, sep="\t")["case_id"])


# ---------------------------------------------------------------------------
# Graph construction
# ---------------------------------------------------------------------------

def build_graph(upheno_edges_path: str,
                train_cases: pd.DataFrame,
                test_disease_ids: set,
                graph2: bool, graph3: bool, graph4: bool,
                human_to_mgi: dict = None,
                mgi_to_phenos: dict = None):
    """
    Returns (triples, gene2pheno, disease2pheno).

    gene2pheno:    gene_symbol → list of phenotype IRIs (MP from MGI + HP from cases)
    disease2pheno: disease_id_norm → list of HP IRIs (from training cases)

    If human_to_mgi and mgi_to_phenos are provided, MGI mouse-knockout MP
    phenotypes (via human ortholog mapping) are used as the primary gene-phenotype
    source. Training-case HPO terms are added as a fallback for genes with no
    mouse ortholog.
    """
    logger.info("Loading UPheno OWL2VecStar edges (Graph 1)...")
    triples = []
    graph1_entities: set[str] = set()
    with open(upheno_edges_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            src, rel, dst = parts
            triples.append((src, rel, dst))
            graph1_entities.add(src)
            graph1_entities.add(dst)

    logger.info(f"  {len(triples):,} Graph-1 triples, {len(graph1_entities):,} entities")

    # ---- gene2pheno: MGI mouse knockouts (primary) + training cases (fallback) ----
    gene2pheno: dict[str, set] = {}
    use_mgi = human_to_mgi is not None and mgi_to_phenos is not None

    if use_mgi:
        mgi_covered = 0
        for gene in train_cases["gene_symbol"].dropna().unique():
            mgi_iri = human_to_mgi.get(gene)
            if mgi_iri and mgi_iri in mgi_to_phenos:
                gene2pheno[gene] = set(mgi_to_phenos[mgi_iri])
                mgi_covered += 1
        logger.info(f"  MGI phenotypes loaded for {mgi_covered} genes")

    # Fallback: training-case HPO terms for genes not covered by MGI
    fallback = 0
    for _, row in train_cases.iterrows():
        gene = row["gene_symbol"]
        if not isinstance(gene, str):
            continue
        if use_mgi and gene in gene2pheno:
            continue    # already have richer MGI phenotypes
        for hp in row["hpo_list"]:
            iri = hp_iri(hp)
            if iri not in graph1_entities:
                continue
            gene2pheno.setdefault(gene, set()).add(iri)
            fallback += 1
    if fallback:
        logger.info(f"  Training-case HPO fallback: {fallback} edges for genes without MGI data")

    gene2pheno = {g: sorted(p) for g, p in gene2pheno.items()}

    # ---- disease2pheno: aggregate over training cases ----
    disease2pheno: dict[str, set] = {}
    for _, row in train_cases.iterrows():
        disease = row["disease_id_norm"]
        for hp in row["hpo_list"]:
            iri = hp_iri(hp)
            if iri not in graph1_entities:
                continue
            disease2pheno.setdefault(disease, set()).add(iri)
    disease2pheno = {d: sorted(p) for d, p in disease2pheno.items()}

    if graph2:
        logger.info("Adding Graph-2 edges (gene → MP/HP)...")
        n = 0
        for gene, phenos in gene2pheno.items():
            g = gene_iri(gene)
            for p in phenos:
                triples.append((g, "http://mowl.borg/has_phenotype", p))
                n += 1
        logger.info(f"  {n:,} gene-phenotype triples added")

    if graph3:
        logger.info("Adding Graph-3 edges (disease → HP, excluding test diseases)...")
        n = 0
        for disease, phenos in disease2pheno.items():
            if disease in test_disease_ids:
                continue    # inductive: keep test diseases unseen
            d = disease_iri(disease)
            for p in phenos:
                triples.append((d, "http://mowl.borg/has_symptom", p))
                n += 1
        logger.info(f"  {n:,} disease-HP triples added")

    if graph4:
        logger.info("Adding Graph-4 edges (gene → associated_with → disease)...")
        n = 0
        seen = set()
        for _, row in train_cases.iterrows():
            gene = row["gene_symbol"]
            disease = row["disease_id_norm"]
            if disease in test_disease_ids:
                continue
            pair = (gene, disease)
            if pair in seen:
                continue
            seen.add(pair)
            g, d = gene_iri(gene), disease_iri(disease)
            triples.append((g, "http://mowl.borg/associated_with", d))
            n += 1
        logger.info(f"  {n:,} gene-disease triples added")

    # Verify inductive constraint
    trained_entities = {t[0] for t in triples} | {t[2] for t in triples}
    overlapping = test_disease_ids & {disease_iri(d) for d in test_disease_ids} & trained_entities
    if graph3 and overlapping:
        logger.warning(f"{len(overlapping)} test disease IRIs are in the training graph!")

    logger.info(f"Total: {len(triples):,} triples")
    return triples, gene2pheno, disease2pheno


# ---------------------------------------------------------------------------
# BMA scoring (adapted from indigena/evaluation.py → compare_vectorized)
# ---------------------------------------------------------------------------

def compare_bma(all_genes_pheno_vecs: th.Tensor,
                disease_pheno_vecs: th.Tensor,
                gene_pheno_counts: th.Tensor) -> th.Tensor:
    """
    Vectorized BMA between one disease and all genes.

    all_genes_pheno_vecs : (G, max_P, D)  — padded with zeros for shorter genes
    disease_pheno_vecs   : (Q, D)
    gene_pheno_counts    : (G,)            — actual phenotype count per gene

    Returns scores (G,) in [0, 1].
    """
    num_genes, max_phenos, emb_dim = all_genes_pheno_vecs.shape

    # (G*max_P, D) x (D, Q) → (G*max_P, Q)
    sim = th.matmul(all_genes_pheno_vecs.view(-1, emb_dim), disease_pheno_vecs.T)
    # Zero entries come from padding — mask before sigmoid
    sim[sim == 0] = -th.inf
    sim = th.sigmoid(sim).view(num_genes, max_phenos, -1)

    # Gene-centric: for each gene phenotype, best-match across disease phenotypes
    gene_max, _ = sim.max(dim=-1)           # (G, max_P)
    gene_centric = gene_max.sum(dim=-1) / th.clamp(gene_pheno_counts, min=1.0)

    # Disease-centric: for each disease phenotype, best-match across gene phenotypes
    disease_max, _ = sim.max(dim=1)         # (G, Q)
    disease_centric = disease_max.mean(dim=-1)

    return (gene_centric + disease_centric) / 2.0


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def evaluate(model, test_cases: pd.DataFrame, gene2pheno: dict,
             eval_genes: list, triples_factory, out_file: str):
    """
    Score every test case against all eval_genes and write TSV.
    Genes without any known phenotype (not in gene2pheno or empty) get score 0.
    """
    entity_to_id = triples_factory.entity_to_id
    entity_ids = th.tensor(list(entity_to_id.values()))
    entity_embs = model.entity_representations[0](indices=entity_ids.to(next(model.parameters()).device))
    entity_embs = entity_embs.cpu().detach()
    emb_dim = entity_embs.shape[1]

    gene_to_idx = {g: i for i, g in enumerate(eval_genes)}
    n_genes = len(eval_genes)

    # Pre-build padded gene phenotype tensor
    pheno_counts = []
    for gene in eval_genes:
        phenos = gene2pheno.get(gene, [])
        phenos = [p for p in phenos if p in entity_to_id]
        pheno_counts.append(len(phenos))
    max_phenos = max(pheno_counts) if pheno_counts else 1

    all_gene_vecs = th.zeros(n_genes, max_phenos, emb_dim)
    for i, gene in enumerate(eval_genes):
        phenos = gene2pheno.get(gene, [])
        phenos = [p for p in phenos if p in entity_to_id]
        if phenos:
            ids = th.tensor([entity_to_id[p] for p in phenos])
            all_gene_vecs[i, :len(phenos), :] = entity_embs[ids]

    pheno_counts_t = th.tensor(pheno_counts, dtype=th.float32)

    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    skipped = 0
    with open(out_file, "w") as f:
        for _, row in tqdm(test_cases.iterrows(), total=len(test_cases), desc="Evaluating"):
            gene = row["gene_symbol"]
            case_id = row["case_id"]

            if gene not in gene_to_idx:
                skipped += 1
                continue

            gene_idx = gene_to_idx[gene]

            # Disease phenotypes = this case's HPO terms
            disease_phenos = [hp_iri(h) for h in row["hpo_list"] if hp_iri(h) in entity_to_id]
            if not disease_phenos:
                # No known HP terms → write uniform-zero scores (rank = random)
                scores = [0.0] * n_genes
            else:
                d_vecs = entity_embs[th.tensor([entity_to_id[p] for p in disease_phenos])]
                with th.no_grad():
                    scores = compare_bma(all_gene_vecs, d_vecs, pheno_counts_t).tolist()

            f.write(
                f"{gene}\t{case_id}\t{gene_idx}\t"
                + "\t".join(f"{s:.6f}" for s in scores)
                + "\n"
            )

    if skipped:
        logger.warning(f"Skipped {skipped} cases whose causal gene is not in eval_genes.")

    return out_file


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

@ck.command()
@ck.option("--upheno-edges", required=True,
           help="Pre-computed UPheno OWL2VecStar edges TSV "
                "(indigena/data/upheno_owl2vecstar_edges.tsv).")
@ck.option("--mgi-gene-phenotypes", default=None,
           help="MGI gene→MP phenotype CSV (indigena/data/gene_phenotypes.csv). "
                "If given, used as primary gene-phenotype source for Graph 2.")
@ck.option("--hom-file", default=None,
           help="MGI human-mouse ortholog file (HOM_MouseHumanSequence.rpt). "
                "Required when --mgi-gene-phenotypes is set.")
@ck.option("--track", type=ck.Choice(["1", "2"]), default="1", show_default=True)
@ck.option("--eval-split", type=ck.Choice(["val", "test"]), default="test", show_default=True,
           help="Split to evaluate on.")
@ck.option("--graph2", is_flag=True, help="Add gene-phenotype edges (Graph 2).")
@ck.option("--graph3", is_flag=True, help="Add disease-HP edges, inductive (Graph 3).")
@ck.option("--graph4", is_flag=True, help="Add gene-disease edges (Graph 4, supervised).")
@ck.option("--embedding-dim", type=int, default=100, show_default=True)
@ck.option("--batch-size", type=int, default=2048, show_default=True)
@ck.option("--learning-rate", type=float, default=0.001, show_default=True)
@ck.option("--num-epochs", type=int, default=300, show_default=True)
@ck.option("--random-seed", type=int, default=0, show_default=True)
@ck.option("--only-eval", is_flag=True,
           help="Skip training; load existing model checkpoint and evaluate.")
def main(upheno_edges, mgi_gene_phenotypes, hom_file,
         track, eval_split, graph2, graph3, graph4,
         embedding_dim, batch_size, learning_rate, num_epochs, random_seed, only_eval):

    track = int(track)
    if graph4:
        graph3 = True
    if graph3:
        graph2 = True

    seed_everything(random_seed)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(MODELS_DIR, exist_ok=True)

    # ---- Load cases ----
    logger.info("Loading PAVS cases...")
    all_cases = load_cases(track)
    eval_genes = sorted(all_cases["gene_symbol"].dropna().unique().tolist())
    logger.info(f"Gene pool: {len(eval_genes)} genes")

    train_ids = load_split_ids("train")
    val_ids = load_split_ids("val")
    eval_ids = load_split_ids(eval_split)

    train_cases = all_cases[all_cases["case_id"].isin(train_ids)].reset_index(drop=True)
    val_cases = all_cases[all_cases["case_id"].isin(val_ids)].reset_index(drop=True)
    eval_cases = all_cases[all_cases["case_id"].isin(eval_ids)].reset_index(drop=True)
    logger.info(f"Train: {len(train_cases)} | Val: {len(val_cases)} | Eval ({eval_split}): {len(eval_cases)}")

    # Disease IDs in test/val (for inductive exclusion)
    test_disease_ids = set(eval_cases["disease_id_norm"]) | set(val_cases["disease_id_norm"])

    # ---- MGI ortholog + phenotype data ----
    human_to_mgi, mgi_to_phenos = None, None
    mgi_tag = ""
    if mgi_gene_phenotypes and hom_file:
        logger.info("Loading MGI ortholog mapping...")
        human_to_mgi = load_human_to_mgi(hom_file)
        logger.info(f"  {len(human_to_mgi)} human→MGI mappings")
        # gene_phenotypes.csv needs graph1_entities filter; load it lazily inside build_graph
        # Pass file path and let build_graph load after graph1 is built
        mgi_tag = "_mgi"

    # ---- Build graph ----
    graph_tag = ("4" if graph4 else "3" if graph3 else "2" if graph2 else "1")
    file_id = (
        f"indigena_transd_track{track}_graph{graph_tag}{mgi_tag}"
        f"_seed{random_seed}_dim{embedding_dim}_bs{batch_size}_lr{learning_rate}"
    )
    model_path = os.path.join(MODELS_DIR, f"{file_id}.pt")

    def _build():
        nonlocal mgi_to_phenos
        # Load MGI phenotypes after graph1 is known (needs graph1_entities filter)
        # We pass the file path and load inside build_graph; for simplicity, pre-load here
        # with a placeholder filter (will be refined inside build_graph)
        h2m = human_to_mgi
        m2p = None
        if mgi_gene_phenotypes and hom_file:
            # Temporarily load without entity filter to avoid chicken-and-egg;
            # build_graph will filter to graph1_entities internally
            m2p_raw = pd.read_csv(mgi_gene_phenotypes)
            m2p = {}
            for _, row in m2p_raw.iterrows():
                g, p = str(row["Gene"]), str(row["Phenotype"])
                m2p.setdefault(g, set()).add(p)
            m2p = {g: sorted(ps) for g, ps in m2p.items()}
            mgi_to_phenos = m2p
        return build_graph(
            upheno_edges, train_cases, test_disease_ids,
            graph2, graph3, graph4,
            human_to_mgi=h2m, mgi_to_phenos=m2p
        )

    if not only_eval:
        triples, gene2pheno, disease2pheno = _build()
        logger.info(f"Gene2pheno entries: {len(gene2pheno)} | Disease2pheno: {len(disease2pheno)}")

        mowl_triples = [Edge(s, r, d) for s, r, d in sorted(triples)]
        triples_factory = Edge.as_pykeen(mowl_triples)
        logger.info(f"Triples factory: {triples_factory.num_entities} entities, "
                    f"{triples_factory.num_relations} relations")

        device = "cuda" if th.cuda.is_available() else "cpu"
        logger.info(f"Training on {device}")
        model = TransD(
            triples_factory=triples_factory,
            embedding_dim=embedding_dim,
            relation_dim=embedding_dim,
            random_seed=random_seed,
        ).to(device)

        optimizer = Adam(params=model.get_grad_params(), lr=learning_rate)
        training_loop = SLCWATrainingLoop(
            model=model,
            triples_factory=triples_factory,
            optimizer=optimizer,
        )

        logger.info(f"Training for {num_epochs} epochs...")
        training_loop.train(
            triples_factory=triples_factory,
            num_epochs=num_epochs,
            batch_size=batch_size,
        )

        th.save(model.state_dict(), model_path)
        logger.info(f"Model saved to {model_path}")

    else:
        triples, gene2pheno, disease2pheno = _build()
        mowl_triples = [Edge(s, r, d) for s, r, d in sorted(triples)]
        triples_factory = Edge.as_pykeen(mowl_triples)

        device = "cuda" if th.cuda.is_available() else "cpu"
        model = TransD(
            triples_factory=triples_factory,
            embedding_dim=embedding_dim,
            relation_dim=embedding_dim,
            random_seed=random_seed,
        ).to(device)
        model.load_state_dict(th.load(model_path, weights_only=True))
        logger.info(f"Loaded model from {model_path}")

    # ---- Evaluate ----
    out_file = os.path.join(RESULTS_DIR, f"{file_id}_{eval_split}.tsv")
    logger.info(f"Evaluating on {eval_split} split ({len(eval_cases)} cases)...")
    model.eval()
    evaluate(model, eval_cases, gene2pheno, eval_genes, triples_factory, out_file)

    _, macro = compute_metrics(out_file, verbose=False)
    print(f"\n=== INDIGENA — Track {track}  Graph {graph_tag}{mgi_tag}  {eval_split} ===")
    print("MR & MRR & Hits@1 & Hits@3 & Hits@10 & Hits@100 & AUC")
    print(" & ".join(f"{macro[k]:.3f}" for k in METRIC_KEYS))
    logger.info(f"Results written to {out_file}")


if __name__ == "__main__":
    main()
