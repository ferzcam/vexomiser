"""
Export entity embeddings from a trained INDIGENA checkpoint to a TSV file
that can be loaded by the Java IndigenaEmbeddings class.

The TSV format is:
    <entity_IRI>\t<v1>,<v2>,...,<vD>

Usage:
    ~/miniforge3/envs/indigena/bin/python eval/export_embeddings.py \\
        --model-path data/models/indigena_transd_track1_graph4_seed0_dim100_bs2048_lr0.001.pt \\
        --upheno-edges ~/Git/indigena/data/upheno_owl2vecstar_edges.tsv \\
        --graph2 --graph3 --graph4 \\
        --track 1 \\
        --output data/models/indigena_track1_graph4_embeddings.tsv
"""

import mowl
mowl.init_jvm("10g")

from mowl.projection import Edge
import torch as th
import click as ck
import pandas as pd
import os
import sys
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

sys.path.insert(0, os.path.dirname(__file__))
from indigena_train import (
    load_cases, load_split_ids, load_human_to_mgi,
    build_graph, MODELS_DIR,
)


@ck.command()
@ck.option("--model-path", required=True, help="Path to .pt checkpoint.")
@ck.option("--upheno-edges", required=True,
           help="upheno_owl2vecstar_edges.tsv used during training.")
@ck.option("--mgi-gene-phenotypes", default=None,
           help="gene_phenotypes.csv (needed if model was trained with --mgi-gene-phenotypes).")
@ck.option("--hom-file", default=None,
           help="HOM_MouseHumanSequence.rpt (needed with --mgi-gene-phenotypes).")
@ck.option("--graph2/--no-graph2", default=True, show_default=True)
@ck.option("--graph3/--no-graph3", default=True, show_default=True)
@ck.option("--graph4/--no-graph4", default=True, show_default=True)
@ck.option("--track", type=int, default=1, show_default=True)
@ck.option("--output", default=None,
           help="Output TSV path. Default: model_path with .pt replaced by _embeddings.tsv.")
@ck.option("--hpo-gene-phenotypes", default=None,
           help="genes_to_phenotype.txt — must be supplied if model was trained with --hpo-gene-phenotypes.")
@ck.option("--no-hpo-fallback", is_flag=True,
           help="Must match the flag used during training if --no-hpo-fallback was set.")
def main(model_path, upheno_edges, mgi_gene_phenotypes, hom_file,
         hpo_gene_phenotypes, graph2, graph3, graph4, track, output, no_hpo_fallback):

    if output is None:
        output = model_path.replace(".pt", "_embeddings.tsv")

    # Reconstruct the same graph that was used during training to recover entity_to_id.
    logger.info("Reconstructing training graph to recover entity vocabulary...")
    all_cases = load_cases(track)
    train_ids = load_split_ids("train")
    val_ids = load_split_ids("val")
    eval_ids = load_split_ids("test")

    train_cases = all_cases[all_cases["case_id"].isin(train_ids)].reset_index(drop=True)
    val_cases = all_cases[all_cases["case_id"].isin(val_ids)].reset_index(drop=True)
    eval_cases = all_cases[all_cases["case_id"].isin(eval_ids)].reset_index(drop=True)
    test_disease_ids = set(eval_cases["disease_id_norm"]) | set(val_cases["disease_id_norm"])

    human_to_mgi = None
    mgi_to_phenos = None
    if mgi_gene_phenotypes and hom_file:
        human_to_mgi = load_human_to_mgi(hom_file)
        m2p_raw = pd.read_csv(mgi_gene_phenotypes)
        mgi_to_phenos = {}
        for _, row in m2p_raw.iterrows():
            g, p = str(row["Gene"]), str(row["Phenotype"])
            mgi_to_phenos.setdefault(g, set()).add(p)
        mgi_to_phenos = {g: sorted(ps) for g, ps in mgi_to_phenos.items()}

    triples, _, _ = build_graph(
        upheno_edges, train_cases, test_disease_ids,
        graph2, graph3, graph4,
        human_to_mgi=human_to_mgi, mgi_to_phenos=mgi_to_phenos,
        hpo_gene_phenos_path=hpo_gene_phenotypes,
        no_hpo_fallback=no_hpo_fallback,
    )

    mowl_triples = [Edge(s, r, d) for s, r, d in sorted(triples)]
    triples_factory = Edge.as_pykeen(mowl_triples)
    n_entities = triples_factory.num_entities
    logger.info(f"Graph has {n_entities} entities")

    # Load the embedding weights directly from the state dict — no need to init the full model.
    logger.info(f"Loading checkpoint: {model_path}")
    state = th.load(model_path, map_location="cpu", weights_only=True)
    emb_weight = state["entity_representations.0._embeddings.weight"]  # (N, D)
    logger.info(f"Embedding matrix: {emb_weight.shape}")

    entity_to_id = triples_factory.entity_to_id          # {iri: idx}
    id_to_entity = {v: k for k, v in entity_to_id.items()}

    logger.info(f"Writing embeddings to {output} ...")
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    with open(output, "w") as f:
        for idx in range(n_entities):
            entity = id_to_entity[idx]
            vec = emb_weight[idx].tolist()
            f.write(entity + "\t" + ",".join(f"{v:.6f}" for v in vec) + "\n")

    logger.info(f"Done. Wrote {n_entities} embeddings ({emb_weight.shape[1]}-dim).")


if __name__ == "__main__":
    main()
