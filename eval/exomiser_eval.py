"""
Evaluate Exomiser phenotype prioritisers (HiPhive, Phive, PhenIX) on the
PAVS benchmark dataset.

Runs on Track 1 (phenotype-only, 5 076 cases) or Track 2 (phenotype +
genotype, 2 815 cases).  Outputs one TSV per prioritiser with per-case scores
across all eval genes, then computes MR / MRR / Hits@1,3,10,100 / AUC.

Output TSV row format (compatible with metrics.compute_metrics):
    gene_symbol <TAB> case_id <TAB> gene_index <TAB> score_0 <TAB> ... score_N

Usage:
    python exomiser_eval.py \\
        --phenotype-data-dir /path/to/2506_phenotype \\
        --track 1
"""

import os
import glob
import json
import logging

import click as ck
import jpype
import numpy as np
import pandas as pd
from tqdm import tqdm

from metrics import compute_metrics, print_as_tex

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(REPO_ROOT, "data")
RESULTS_DIR = os.path.join(DATA_DIR, "results")

CLI_TARGET = os.path.join(REPO_ROOT, "exomiser-cli", "target")

METRIC_KEYS = ["mr", "mrr", "hits@1", "hits@3", "hits@10", "hits@100", "auc"]
TEX_HEADER = "MR & MRR & Hits@1 & Hits@3 & Hits@10 & Hits@100 & AUC"


# ---------------------------------------------------------------------------
# JVM / Exomiser wiring
# ---------------------------------------------------------------------------

def start_jvm(phenotype_data_dir: str):
    jar_main = glob.glob(os.path.join(CLI_TARGET, "exomiser-cli-*.jar"))
    jar_libs = glob.glob(os.path.join(CLI_TARGET, "lib", "*.jar"))
    if not jar_main:
        raise FileNotFoundError(
            f"Exomiser CLI JAR not found in {CLI_TARGET}. "
            "Run: ./mvnw package -pl exomiser-cli -am -DskipTests"
        )
    classpath = jar_libs + jar_main
    jpype.startJVM(
        jpype.getDefaultJVMPath(),
        "-ea",
        "-Xmx10g",
        classpath=classpath,
        convertStrings=True,
    )


def build_priority_factory(phenotype_data_dir: str):
    import jpype.imports  # noqa: E402 — must be after startJVM
    from java.nio.file import Paths as JPaths
    from com.zaxxer.hikari import HikariDataSource
    from org.monarchinitiative.exomiser.core.phenotype.dao import (
        HumanPhenotypeOntologyDao,
        MousePhenotypeOntologyDao,
        ZebraFishPhenotypeOntologyDao,
    )
    from org.monarchinitiative.exomiser.core.phenotype.service import OntologyServiceImpl
    from org.monarchinitiative.exomiser.core.phenotype import PhenotypeMatchService
    from org.monarchinitiative.exomiser.core.prioritisers.service import (
        ModelServiceImpl,
        PriorityService,
    )
    from org.monarchinitiative.exomiser.core.prioritisers.dao import DefaultDiseaseDao
    from org.monarchinitiative.exomiser.core.prioritisers import (
        PriorityFactoryImpl,
        HiPhiveOptions,
    )
    from org.monarchinitiative.exomiser.core.prioritisers.util import DataMatrixIO

    # The H2 database file has the same stem as the directory name
    # H2 requires an absolute path
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

    human_dao = HumanPhenotypeOntologyDao(ds)
    mouse_dao = MousePhenotypeOntologyDao(ds)
    fish_dao = ZebraFishPhenotypeOntologyDao(ds)

    ontology_service = OntologyServiceImpl(human_dao, mouse_dao, fish_dao)
    phenotype_match_service = PhenotypeMatchService(ontology_service)
    model_service = ModelServiceImpl(ds)
    disease_dao = DefaultDiseaseDao(ds)
    priority_service = PriorityService(
        model_service, phenotype_match_service, disease_dao
    )

    rw_path = JPaths.get(os.path.join(phenotype_data_dir, "rw_string_10.mv"))
    data_matrix = DataMatrixIO.loadOffHeapDataMatrix(rw_path)

    phenix_dir = JPaths.get(os.path.join(phenotype_data_dir, "phenix"))

    factory = PriorityFactoryImpl(priority_service, data_matrix, phenix_dir)
    return factory, ds


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_cases(track: int) -> pd.DataFrame:
    path = os.path.join(DATA_DIR, f"track{track}_cases.tsv")
    df = pd.read_csv(path, sep="\t")
    # hpo_terms column: "HP:0001263|label;HP:0000076|label"
    df["hpo_list"] = df["hpo_terms"].apply(
        lambda s: [term.split("|")[0] for term in str(s).split(";")]
        if pd.notna(s) else []
    )
    return df


def load_gene_entrez_map() -> dict:
    """gene_symbol -> entrez_id (int) from PAVS_cases.tsv."""
    path = os.path.join(DATA_DIR, "pavs", "PAVS_cases.tsv")
    df = pd.read_csv(path, sep="\t", usecols=["gene_symbol", "gene_id"])
    mapping = {}
    for _, row in df.iterrows():
        gene_id = str(row["gene_id"])  # e.g. "NCBIGene:285362"
        if gene_id.startswith("NCBIGene:"):
            mapping[row["gene_symbol"]] = int(gene_id.split(":")[1].split("|")[0])
    return mapping


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def run_prioritiser(pname, prioritiser, cases, eval_genes, gene_to_index,
                    java_genes, entrez_to_eval_idx, out_file):
    from java.util import ArrayList
    from java.util.stream import Collectors

    os.makedirs(os.path.dirname(out_file), exist_ok=True)

    with open(out_file, "w") as f:
        for _, row in tqdm(cases.iterrows(), total=len(cases),
                           desc=pname):
            case_id = row["case_id"]
            causal_gene = row["gene_symbol"]
            hpo_ids = row["hpo_list"]

            if causal_gene not in gene_to_index:
                logger.warning(f"Gene {causal_gene} not in eval set, skipping {case_id}")
                continue

            gene_index = gene_to_index[causal_gene]
            scores = [0.0] * len(eval_genes)

            if hpo_ids:
                java_hpo = ArrayList()
                for h in hpo_ids:
                    java_hpo.add(h)

                result_list = (
                    prioritiser.prioritise(java_hpo, java_genes)
                    .collect(Collectors.toList())
                )
                for r in result_list:
                    gid = int(r.geneId())
                    if gid in entrez_to_eval_idx:
                        scores[entrez_to_eval_idx[gid]] = float(r.score())

            f.write(
                f"{causal_gene}\t{case_id}\t{gene_index}\t"
                + "\t".join(str(s) for s in scores)
                + "\n"
            )

    return out_file


def emit(summary_f, text):
    print(text)
    summary_f.write(text + "\n")
    summary_f.flush()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

@ck.command()
@ck.option("--phenotype-data-dir", required=True,
           help="Path to Exomiser phenotype data directory (e.g. /path/to/2406_phenotype)")
@ck.option("--track", type=ck.Choice(["1", "2"]), default="1", show_default=True,
           help="Evaluation track: 1=phenotype-only, 2=phenotype+genotype")
@ck.option("--split", type=ck.Choice(["train", "val", "test", "all"]), default="all",
           show_default=True,
           help="Restrict evaluation to a data split (requires eval/generate_splits.py to have been run)")
def main(phenotype_data_dir, track, split):
    track = int(track)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    logger.info("Loading PAVS cases...")
    all_cases = load_cases(track)

    # Gene pool is always fixed to all genes across the full dataset
    # so rankings are comparable across splits
    eval_genes = sorted(all_cases["gene_symbol"].dropna().unique().tolist())
    gene_to_index = {g: i for i, g in enumerate(eval_genes)}

    if split != "all":
        split_path = os.path.join(DATA_DIR, "splits", f"{split}.tsv")
        if not os.path.exists(split_path):
            raise FileNotFoundError(
                f"Split file not found: {split_path}. "
                "Run: uv run python eval/generate_splits.py"
            )
        split_ids = set(pd.read_csv(split_path, sep="\t")["case_id"])
        cases = all_cases[all_cases["case_id"].isin(split_ids)].reset_index(drop=True)
        logger.info(f"Restricted to '{split}' split: {len(cases)} cases")
    else:
        cases = all_cases

    gene_entrez = load_gene_entrez_map()
    logger.info(f"Track {track}: {len(cases)} cases, {len(eval_genes)} genes in pool")

    logger.info("Starting JVM...")
    start_jvm(phenotype_data_dir)

    import jpype.imports  # noqa
    from java.util import ArrayList
    from org.monarchinitiative.exomiser.core.model import Gene
    from org.monarchinitiative.exomiser.core.prioritisers import HiPhiveOptions

    logger.info("Building Exomiser priority factory...")
    factory, ds = build_priority_factory(phenotype_data_dir)

    java_genes = ArrayList()
    entrez_to_eval_idx = {}
    missing = []
    for i, symbol in enumerate(eval_genes):
        if symbol not in gene_entrez:
            missing.append(symbol)
            continue
        entrez_id = gene_entrez[symbol]
        java_genes.add(Gene(symbol, entrez_id))
        entrez_to_eval_idx[entrez_id] = i

    if missing:
        logger.warning(f"{len(missing)} genes have no Entrez ID and will score 0: {missing[:5]}...")

    prioritisers = {
        "hiphive": factory.makeHiPhivePrioritiser(HiPhiveOptions.defaults()),
        "phive":   factory.makePhivePrioritiser(),
        "phenix":  factory.makePhenixPrioritiser(),
    }

    # Allow running a subset of prioritisers to resume partial runs
    if os.environ.get("PRIORITISERS"):
        subset = os.environ["PRIORITISERS"].split(",")
        prioritisers = {k: v for k, v in prioritisers.items() if k in subset}

    split_tag = f"_{split}" if split != "all" else ""
    summary_path = os.path.join(RESULTS_DIR, f"exomiser_track{track}{split_tag}_summary.txt")
    with open(summary_path, "w") as summary_f:
        emit(summary_f, f"# Exomiser evaluation — Track {track}  split: {split}")
        emit(summary_f, f"# Cases: {len(cases)}  |  Genes: {len(eval_genes)}")
        emit(summary_f, "")

        for pname, prioritiser in prioritisers.items():
            out_file = os.path.join(RESULTS_DIR, f"exomiser_{pname}_track{track}{split_tag}.tsv")
            logger.info(f"Running {pname}...")
            run_prioritiser(pname, prioritiser, cases, eval_genes,
                            gene_to_index, java_genes, entrez_to_eval_idx,
                            out_file)

            _, macro = compute_metrics(out_file, verbose=False)
            emit(summary_f, f"## {pname}")
            emit(summary_f, TEX_HEADER)
            emit(summary_f, " & ".join(f"{macro[k]:.3f}" for k in METRIC_KEYS))
            emit(summary_f, "")

    logger.info(f"Summary written to {summary_path}")
    ds.close()


if __name__ == "__main__":
    main()
