"""
Track 2 Exomiser evaluation (phenotype + variant).

Each case has a spiked VCF (GIAB HG001 hg38 background + causal variant).
Runs the full Exomiser analysis pipeline via CLI subprocess in parallel, then
scores phenotype similarity for each method via JPype and multiplies by the
per-gene variant score.

Methods produced:
  hiphive           — Exomiser EXOMISER_GENE_COMBINED_SCORE (variant × HiPhive pheno)
  phenix            — variant_score × PhenIX pheno score
  phive             — variant_score × Phive pheno score
  indigena_hiphive  — variant_score × INDIGENA HiPhive pheno score (all organisms)
  indigena_phenix   — variant_score × INDIGENA PhenIX pheno score (human only)
  indigena_phive    — variant_score × INDIGENA Phive pheno score  (mouse only)

Usage:
    python track2_eval.py \\
        --phenotype-data-dir exomiser-data/2406_phenotype \\
        --app-props exomiser-data/application.properties \\
        --split test \\
        --embeddings data/models/indigena_track1_graph4_embeddings.tsv

    # Reuse existing CLI results (skip the 10-min CLI phase):
    python track2_eval.py ... --skip-cli
"""

import json
import logging
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import click as ck
import numpy as np
import pandas as pd
from tqdm import tqdm

from metrics import compute_metrics, print_as_tex

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(REPO_ROOT, "data")
RESULTS_DIR = os.path.join(DATA_DIR, "results")
CLI_TARGET = os.path.join(REPO_ROOT, "exomiser-cli", "target")

METRIC_KEYS = ["mr", "mrr", "hits@1", "hits@3", "hits@10", "hits@100", "auc"]
TEX_HEADER = "MR & MRR & Hits@1 & Hits@3 & Hits@10 & Hits@100 & AUC"


# ---------------------------------------------------------------------------
# Phenopacket v1 conversion (Exomiser v15 CLI expects phenopackets v1)
# ---------------------------------------------------------------------------

def make_phenopacket_v1(src_path: str) -> dict:
    """Convert a PAVS phenopacket v2 JSON to a minimal v1 dict for Exomiser."""
    with open(src_path) as f:
        src = json.load(f)
    return {
        "id": src["id"],
        "subject": {
            "id": src["subject"]["id"],
            "sex": src["subject"].get("sex", "UNKNOWN"),
        },
        "phenotypicFeatures": [
            {"type": pf["type"], "negated": False}
            for pf in src.get("phenotypicFeatures", [])
            if not pf.get("excluded", False) and pf["type"]["id"].startswith("HP:")
        ],
        "metaData": {
            "created": "2024-01-01T00:00:00.000Z",
            "createdBy": "PAVS",
            "resources": [{
                "id": "hp",
                "name": "human phenotype ontology",
                "url": "http://purl.obolibrary.org/obo/hp.owl",
                "version": "hp/releases/2024-04-04",
                "namespacePrefix": "HP",
                "iriPrefix": "http://purl.obolibrary.org/obo/HP_",
            }],
            "phenopacketSchemaVersion": "1.0",
        },
    }


# ---------------------------------------------------------------------------
# CLI runner (called in worker processes)
# ---------------------------------------------------------------------------

def run_one_case(args):
    """Run Exomiser CLI for a single case. Returns (case_id, tsv_path | None)."""
    case_id, phenopacket_path, vcf_path, app_props_path, jar_path, work_dir = args

    case_out_dir = os.path.join(work_dir, case_id.replace(":", "_"))
    os.makedirs(case_out_dir, exist_ok=True)

    pp_v1 = make_phenopacket_v1(phenopacket_path)
    pp_path = os.path.join(case_out_dir, "phenopacket.json")
    with open(pp_path, "w") as f:
        json.dump(pp_v1, f)

    safe_id = case_id.replace(":", "_")
    expected_tsv = os.path.join(case_out_dir, f"{safe_id}.genes.tsv")

    cmd = [
        "java", "-Xmx4g",
        f"-Dspring.config.location=file:{app_props_path}",
        "-jar", jar_path,
        "analyse",
        "--sample", pp_path,
        "--vcf", vcf_path,
        "--assembly", "GRCh38",
        "--preset", "exome",
        "--output-directory", case_out_dir,
        "--output-filename", safe_id,
        "--output-format", "TSV_GENE",
    ]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=300
        )
        if result.returncode != 0:
            err_path = os.path.join(case_out_dir, "error.txt")
            with open(err_path, "w") as ef:
                ef.write(result.stderr[-2000:])
            return case_id, None
        if not os.path.exists(expected_tsv):
            return case_id, None
        return case_id, expected_tsv
    except subprocess.TimeoutExpired:
        return case_id, None
    except Exception:
        return case_id, None


# ---------------------------------------------------------------------------
# Parse Exomiser TSV_GENE output
# ---------------------------------------------------------------------------

def parse_gene_tsv(tsv_path: str) -> dict:
    """Return {gene_symbol: (combined_score, variant_score)} from a TSV_GENE file."""
    scores = {}
    with open(tsv_path) as f:
        header = None
        for line in f:
            if line.startswith("#"):
                header = line.lstrip("#").strip().split("\t")
                continue
            if header is None or not line.strip():
                continue
            parts = line.strip().split("\t")
            row = dict(zip(header, parts))
            gene = row.get("GENE_SYMBOL", "")
            try:
                combined = float(row.get("EXOMISER_GENE_COMBINED_SCORE", 0))
                variant = float(row.get("EXOMISER_GENE_VARIANT_SCORE", 0))
            except ValueError:
                combined, variant = 0.0, 0.0
            if gene:
                if gene not in scores or combined > scores[gene][0]:
                    scores[gene] = (combined, variant)
    return scores


# ---------------------------------------------------------------------------
# JVM / phenotype scoring
# ---------------------------------------------------------------------------

def start_jvm(phenotype_data_dir: str, jar_path: str):
    import glob
    import jpype
    jar_libs = glob.glob(os.path.join(os.path.dirname(jar_path), "lib", "*.jar"))
    core_jars = glob.glob(os.path.join(
        REPO_ROOT, "exomiser-core", "target", "exomiser-core-*.jar"
    ))
    classpath = jar_libs + [jar_path] + core_jars
    jpype.startJVM(
        jpype.getDefaultJVMPath(),
        "-ea", "-Xmx10g",
        classpath=classpath,
        convertStrings=True,
    )


def build_all_prioritisers(phenotype_data_dir: str, embeddings_path=None):
    """
    Build all phenotype prioritisers. Always includes phenix and phive.
    If embeddings_path is given, also includes INDIGENA variants.

    Returns (prioritisers_dict, ds).
    """
    import jpype.imports  # noqa
    from java.nio.file import Paths as JPaths
    from com.zaxxer.hikari import HikariDataSource
    from org.monarchinitiative.exomiser.core.phenotype.dao import (
        HumanPhenotypeOntologyDao, MousePhenotypeOntologyDao, ZebraFishPhenotypeOntologyDao,
    )
    from org.monarchinitiative.exomiser.core.phenotype.service import OntologyServiceImpl
    from org.monarchinitiative.exomiser.core.phenotype import PhenotypeMatchService
    from org.monarchinitiative.exomiser.core.prioritisers.service import (
        ModelServiceImpl, PriorityService,
    )
    from org.monarchinitiative.exomiser.core.prioritisers.dao import DefaultDiseaseDao
    from org.monarchinitiative.exomiser.core.prioritisers import (
        PriorityFactoryImpl, HiPhiveOptions, HiPhivePriority, PhivePriority,
    )
    from org.monarchinitiative.exomiser.core.prioritisers.util import DataMatrixIO

    phenotype_data_dir = os.path.abspath(phenotype_data_dir)
    db_stem = os.path.basename(phenotype_data_dir.rstrip("/"))
    db_path = os.path.join(phenotype_data_dir, db_stem)
    jdbc_url = (
        f"jdbc:h2:file:{db_path}"
        ";ACCESS_MODE_DATA=r"
        ";INIT=SET SCHEMA EXOMISER"
    )

    ds = HikariDataSource()
    ds.setJdbcUrl(jdbc_url)
    ds.setUsername("sa")
    ds.setPassword("")

    ontology_service = OntologyServiceImpl(
        HumanPhenotypeOntologyDao(ds),
        MousePhenotypeOntologyDao(ds),
        ZebraFishPhenotypeOntologyDao(ds),
    )
    phenotype_match_service = PhenotypeMatchService(ontology_service)
    model_service = ModelServiceImpl(ds)
    disease_dao = DefaultDiseaseDao(ds)
    priority_service = PriorityService(model_service, phenotype_match_service, disease_dao)

    rw_path = JPaths.get(os.path.join(phenotype_data_dir, "rw_string_10.mv"))
    data_matrix = DataMatrixIO.loadOffHeapDataMatrix(rw_path)

    phenix_dir = JPaths.get(os.path.join(phenotype_data_dir, "phenix"))
    factory = PriorityFactoryImpl(priority_service, data_matrix, phenix_dir)

    prioritisers = {
        "phenix": factory.makePhenixPrioritiser(),
        "phive":  factory.makePhivePrioritiser(),
    }

    # Resnik-HiPhive / Resnik-Phive / Resnik-PhenIX (always added)
    from org.monarchinitiative.exomiser.core.phenotype import ResnikModelScorerFactory
    resnik_factory = ResnikModelScorerFactory()
    hiphive_opts = HiPhiveOptions.defaults()
    phenix_opts  = HiPhiveOptions.builder().runParams("human").build()
    prioritisers.update({
        "resnik_hiphive": HiPhivePriority(
            hiphive_opts, data_matrix, priority_service, resnik_factory
        ),
        "resnik_phive":  PhivePriority(priority_service, resnik_factory),
        "resnik_phenix": HiPhivePriority(
            phenix_opts, data_matrix, priority_service, resnik_factory
        ),
    })

    if embeddings_path:
        from org.monarchinitiative.exomiser.core.phenotype import (
            IndigenaEmbeddings, IndigenaModelScorerFactory,
        )
        emb = IndigenaEmbeddings.load(JPaths.get(os.path.abspath(embeddings_path)))
        scorer_factory = IndigenaModelScorerFactory(emb)

        prioritisers.update({
            "indigena_hiphive": HiPhivePriority(
                hiphive_opts, data_matrix, priority_service, scorer_factory
            ),
            "indigena_phive":  PhivePriority(priority_service, scorer_factory),
            "indigena_phenix": HiPhivePriority(
                phenix_opts, data_matrix, priority_service, scorer_factory
            ),
        })

    return prioritisers, ds


def score_pheno(prioritiser, cases, eval_genes, gene_entrez, desc="scoring"):
    """Return {case_id: {gene_symbol: pheno_score}} using any Exomiser prioritiser."""
    import jpype.imports  # noqa
    from java.util import ArrayList
    from java.util.stream import Collectors
    from org.monarchinitiative.exomiser.core.model import Gene

    java_genes = ArrayList()
    entrez_to_sym = {}
    for sym in eval_genes:
        if sym in gene_entrez:
            eid = gene_entrez[sym]
            java_genes.add(Gene(sym, eid))
            entrez_to_sym[eid] = sym

    results = {}
    for _, row in tqdm(cases.iterrows(), total=len(cases), desc=desc):
        case_id = row["case_id"]
        hpo_ids = row["hpo_list"]
        pheno_scores = {g: 0.0 for g in eval_genes}
        if hpo_ids:
            java_hpo = ArrayList()
            for h in hpo_ids:
                java_hpo.add(h)
            for r in prioritiser.prioritise(java_hpo, java_genes).collect(Collectors.toList()):
                sym = entrez_to_sym.get(int(r.geneId()))
                if sym:
                    pheno_scores[sym] = float(r.score())
        results[case_id] = pheno_scores
    return results


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_cases(split: str) -> pd.DataFrame:
    df = pd.read_csv(os.path.join(DATA_DIR, "track2_cases.tsv"), sep="\t")
    df = df[df["has_vcf"] == True].reset_index(drop=True)
    df["hpo_list"] = df["hpo_terms"].apply(
        lambda s: [t.split("|")[0] for t in str(s).split(";")]
        if pd.notna(s) else []
    )
    if split != "all":
        split_path = os.path.join(DATA_DIR, "splits", f"{split}.tsv")
        if not os.path.exists(split_path):
            raise FileNotFoundError(
                f"Split file not found: {split_path}. "
                "Run: python eval/generate_splits.py"
            )
        split_ids = set(pd.read_csv(split_path, sep="\t")["case_id"])
        df = df[df["case_id"].isin(split_ids)].reset_index(drop=True)
        logger.info(f"Restricted to '{split}' split: {len(df)} cases with VCF")
    return df


def load_gene_entrez_map() -> dict:
    path = os.path.join(DATA_DIR, "pavs", "PAVS_cases.tsv")
    df = pd.read_csv(path, sep="\t", usecols=["gene_symbol", "gene_id"])
    mapping = {}
    for _, row in df.iterrows():
        gid = str(row["gene_id"])
        if gid.startswith("NCBIGene:"):
            mapping[row["gene_symbol"]] = int(gid.split(":")[1].split("|")[0])
    return mapping


def emit(f, text):
    print(text)
    f.write(text + "\n")
    f.flush()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

@ck.command()
@ck.option("--phenotype-data-dir", required=True,
           help="Path to Exomiser phenotype data dir (e.g. exomiser-data/2406_phenotype)")
@ck.option("--app-props", required=True,
           help="Path to application.properties with hg38 genome data configured")
@ck.option("--split", type=ck.Choice(["train", "val", "test", "all"]), default="test",
           show_default=True)
@ck.option("--workers", default=16, show_default=True,
           help="Number of parallel Exomiser CLI workers")
@ck.option("--embeddings", default=None,
           help="Path to INDIGENA embeddings TSV. When supplied, also runs "
                "indigena_hiphive / indigena_phenix / indigena_phive.")
@ck.option("--work-dir", default=None,
           help="Directory for temp CLI job files and outputs. "
                "Defaults to data/results/track2_jobs/")
@ck.option("--skip-cli", is_flag=True, default=False,
           help="Skip Exomiser CLI phase and reuse existing TSV files in --work-dir.")
def main(phenotype_data_dir, app_props, split, workers, embeddings, work_dir, skip_cli):
    os.makedirs(RESULTS_DIR, exist_ok=True)

    import glob as _glob
    jar_path = None
    for j in _glob.glob(os.path.join(CLI_TARGET, "exomiser-cli-*.jar")):
        jar_path = j
        break
    if not jar_path:
        raise FileNotFoundError(f"Exomiser CLI JAR not found in {CLI_TARGET}")

    if work_dir is None:
        work_dir = os.path.join(RESULTS_DIR, "track2_jobs")
    os.makedirs(work_dir, exist_ok=True)

    logger.info("Loading Track 2 cases...")
    cases = load_cases(split)
    gene_entrez = load_gene_entrez_map()

    eval_genes = sorted(cases["gene_symbol"].dropna().unique().tolist())
    gene_to_index = {g: i for i, g in enumerate(eval_genes)}
    logger.info(f"{len(cases)} cases, {len(eval_genes)} genes in pool")

    vcf_dir = os.path.join(DATA_DIR, "spiked_vcfs_fixed")
    phenopacket_dir = os.path.join(DATA_DIR, "pavs", "phenopackets")

    split_tag = f"_{split}" if split != "all" else ""

    # ------------------------------------------------------------------
    # Phase 1: run Exomiser CLI in parallel (or reuse existing results)
    # ------------------------------------------------------------------
    cli_results = {}  # case_id -> tsv_path | None

    if skip_cli:
        logger.info("--skip-cli: loading existing TSV results...")
        for _, row in cases.iterrows():
            cid = row["case_id"]
            safe_id = cid.replace(":", "_")
            tsv = os.path.join(work_dir, safe_id, f"{safe_id}.genes.tsv")
            cli_results[cid] = tsv if os.path.exists(tsv) else None
        n_ok = sum(1 for v in cli_results.values() if v is not None)
        logger.info(f"Found {n_ok}/{len(cases)} existing TSV files")
    else:
        task_args = []
        for _, row in cases.iterrows():
            cid = row["case_id"]
            vcf = os.path.join(vcf_dir, f"{cid}.vcf.gz")
            pp = os.path.join(phenopacket_dir, f"{cid}.json")
            if not os.path.exists(vcf):
                logger.warning(f"Missing VCF for {cid}, skipping")
                continue
            if not os.path.exists(pp):
                logger.warning(f"Missing phenopacket for {cid}, skipping")
                continue
            task_args.append((cid, pp, vcf, os.path.abspath(app_props), jar_path, work_dir))

        logger.info(f"Running Exomiser CLI ({workers} workers)...")
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(run_one_case, args): args[0] for args in task_args}
            for fut in tqdm(as_completed(futures), total=len(futures), desc="Exomiser CLI"):
                cid, tsv_path = fut.result()
                cli_results[cid] = tsv_path

        n_ok = sum(1 for v in cli_results.values() if v is not None)
        logger.info(f"CLI complete: {n_ok}/{len(task_args)} cases succeeded")

    # ------------------------------------------------------------------
    # Phase 2: collect variant/combined scores; write hiphive output
    # ------------------------------------------------------------------
    hiphive_out = os.path.join(RESULTS_DIR, f"exomiser_hiphive_track2{split_tag}.tsv")
    variant_scores = {}  # case_id -> {gene_symbol: variant_score}

    with open(hiphive_out, "w") as f:
        for _, row in cases.iterrows():
            cid = row["case_id"]
            causal = row["gene_symbol"]
            if causal not in gene_to_index:
                continue

            gene_combined = {g: 0.0 for g in eval_genes}
            gene_variant = {g: 0.0 for g in eval_genes}

            tsv = cli_results.get(cid)
            if tsv:
                for g, (combined, variant) in parse_gene_tsv(tsv).items():
                    if g in gene_combined:
                        gene_combined[g] = combined
                        gene_variant[g] = variant

            variant_scores[cid] = gene_variant
            scores = [gene_combined[g] for g in eval_genes]
            f.write(
                f"{causal}\t{cid}\t{gene_to_index[causal]}\t"
                + "\t".join(str(s) for s in scores) + "\n"
            )

    # ------------------------------------------------------------------
    # Phase 3: JPype phenotype scoring for phenix, phive, indigena_*
    # ------------------------------------------------------------------
    logger.info("Starting JVM for phenotype scoring...")
    start_jvm(phenotype_data_dir, jar_path)
    prioritisers, ds = build_all_prioritisers(phenotype_data_dir, embeddings)

    output_paths = {"hiphive": hiphive_out}

    for pname, prioritiser in prioritisers.items():
        pheno = score_pheno(prioritiser, cases, eval_genes, gene_entrez, desc=pname)
        out_path = os.path.join(RESULTS_DIR, f"exomiser_{pname}_track2{split_tag}.tsv")
        output_paths[pname] = out_path
        with open(out_path, "w") as f:
            for _, row in cases.iterrows():
                cid = row["case_id"]
                causal = row["gene_symbol"]
                if causal not in gene_to_index:
                    continue
                v_scores = variant_scores.get(cid, {})
                p_scores = pheno.get(cid, {})
                scores = [v_scores.get(g, 0.0) * p_scores.get(g, 0.0) for g in eval_genes]
                f.write(
                    f"{causal}\t{cid}\t{gene_to_index[causal]}\t"
                    + "\t".join(str(s) for s in scores) + "\n"
                )

    ds.close()

    # ------------------------------------------------------------------
    # Metrics
    # ------------------------------------------------------------------
    summary_path = os.path.join(RESULTS_DIR, f"exomiser_track2{split_tag}_summary.txt")
    method_order = [
        "hiphive", "indigena_hiphive", "resnik_hiphive",
        "phenix",  "indigena_phenix",  "resnik_phenix",
        "phive",   "indigena_phive",   "resnik_phive",
    ]
    with open(summary_path, "w") as sf:
        emit(sf, f"# Exomiser Track 2 — split: {split}")
        emit(sf, f"# Cases: {len(cases)}  |  Genes: {len(eval_genes)}")
        emit(sf, "")
        for label in method_order:
            path = output_paths.get(label)
            if path is None:
                continue
            _, macro = compute_metrics(path, verbose=False)
            emit(sf, f"## {label}")
            emit(sf, TEX_HEADER)
            emit(sf, " & ".join(f"{macro[k]:.3f}" for k in METRIC_KEYS))
            emit(sf, "")

    logger.info(f"Summary written to {summary_path}")


if __name__ == "__main__":
    main()
