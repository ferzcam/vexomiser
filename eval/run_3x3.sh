#!/usr/bin/env bash
# Run the 3×3 INDIGENA training matrix on the workstation (2 GPUs).
#
# Three training configs × 5 eval modes each (two 3×3 tables).
# G2 source determines the model name; eval modes are fixed at run time.
#
# Usage:
#   bash eval/run_3x3.sh
#
# Prerequisites (run from vexomiser root):
#   conda activate indigena   (or equivalent)
#   uv run python eval/generate_splits.py   (if splits not yet generated)

set -euo pipefail

PYTHON=~/miniforge3/envs/indigena/bin/python
SCRIPT=eval/indigena_train.py

UPHENO=~/Git/indigena/data/upheno_owl2vecstar_edges.tsv
MGI_CSV=~/Git/indigena/data/gene_phenotypes.csv
HOM=~/Git/indigena/data/HOM_MouseHumanSequence.rpt
HPO_G2P=~/Git/vexomiser/data/genes_to_phenotype.txt
PHENO_HPOA=~/Git/vexomiser/data/phenotype.hpoa

COMMON="--upheno-edges $UPHENO
        --hom-file $HOM
        --phenotype-hpoa $PHENO_HPOA
        --eval-gene-phenotypes $HPO_G2P
        --graph2 --graph3 --graph4
        --track 1 --eval-split test"

echo "========================================================"
echo " Config 1/3: G2 = MGI only  (GPU 0)"
echo "========================================================"
CUDA_VISIBLE_DEVICES=0 $PYTHON $SCRIPT \
    $COMMON \
    --mgi-gene-phenotypes $MGI_CSV \
    --no-hpo-fallback \
    &

echo "========================================================"
echo " Config 2/3: G2 = HP only   (GPU 1)"
echo "========================================================"
CUDA_VISIBLE_DEVICES=1 $PYTHON $SCRIPT \
    $COMMON \
    --hpo-gene-phenotypes $HPO_G2P \
    --no-hpo-fallback \
    &

# Wait for GPU 0 and 1 to free up before running config 3
wait

echo "========================================================"
echo " Config 3/3: G2 = MGI + HP  (GPU 0)"
echo "========================================================"
CUDA_VISIBLE_DEVICES=0 $PYTHON $SCRIPT \
    $COMMON \
    --mgi-gene-phenotypes $MGI_CSV \
    --hpo-gene-phenotypes $HPO_G2P \
    &

wait
echo "All 3 training runs complete."
