#!/usr/bin/env bash
# Run Exomiser integration (Track 1 + Track 2) with all three 3×3 training configs.
#
# Steps:
#   1. Export embeddings for MGI-only and HP-only models (MGI+HP already exported).
#   2. Run Track 1 (exomiser_eval.py) with each embedding set.
#   3. Run Track 2 (track2_eval.py --skip-cli) with each embedding set.
#
# Usage (from vexomiser root):
#   bash eval/run_exomiser_3x3.sh
#
# Prerequisites:
#   - 3×3 training completed (eval/run_3x3.sh has been run)
#   - Track 2 CLI phase already done (spiked VCFs processed)
#   - conda activate indigena (or equivalent)

set -euo pipefail

PYTHON=~/miniforge3/envs/indigena/bin/python
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

UPHENO=~/Git/indigena/data/upheno_owl2vecstar_edges.tsv
MGI_CSV=~/Git/indigena/data/gene_phenotypes.csv
HOM=~/Git/indigena/data/HOM_MouseHumanSequence.rpt
HPO_G2P=$REPO_ROOT/data/genes_to_phenotype.txt

MODELS=$REPO_ROOT/data/models
PHENO_DIR=$REPO_ROOT/exomiser-data/2406_phenotype
APP_PROPS=$REPO_ROOT/exomiser-data/application.properties

# Model checkpoints from 3×3 training
PT_MGI=$MODELS/indigena_transd_track1_graph4_mgi_nofallback_seed0_dim100_bs2048_lr0.001.pt
PT_HPO=$MODELS/indigena_transd_track1_graph4_hpo_seed0_dim100_bs2048_lr0.001.pt
PT_MGI_HPO=$MODELS/indigena_transd_track1_graph4_mgi_hpo_seed0_dim100_bs2048_lr0.001.pt

# Embedding output paths
EMB_MGI=$MODELS/indigena_transd_track1_graph4_mgi_nofallback_seed0_dim100_bs2048_lr0.001_embeddings.tsv
EMB_HPO=$MODELS/indigena_transd_track1_graph4_hpo_seed0_dim100_bs2048_lr0.001_embeddings.tsv
EMB_MGI_HPO=$MODELS/indigena_transd_track1_graph4_mgi_hpo_seed0_dim100_bs2048_lr0.001_embeddings.tsv

# -----------------------------------------------------------------------
# Step 1: Export embeddings
# -----------------------------------------------------------------------
echo "========================================================"
echo " Exporting embeddings: MGI-only"
echo "========================================================"
$PYTHON $REPO_ROOT/eval/export_embeddings.py \
    --model-path $PT_MGI \
    --upheno-edges $UPHENO \
    --mgi-gene-phenotypes $MGI_CSV \
    --hom-file $HOM \
    --no-hpo-fallback \
    --graph2 --graph3 --graph4 \
    --track 1 \
    --output $EMB_MGI

echo "========================================================"
echo " Exporting embeddings: HP-only"
echo "========================================================"
$PYTHON $REPO_ROOT/eval/export_embeddings.py \
    --model-path $PT_HPO \
    --upheno-edges $UPHENO \
    --hpo-gene-phenotypes $HPO_G2P \
    --no-hpo-fallback \
    --graph2 --graph3 --graph4 \
    --track 1 \
    --output $EMB_HPO

echo "========================================================"
echo " Embeddings: MGI+HP already exported, skipping."
echo "========================================================"

# -----------------------------------------------------------------------
# Step 2: Track 1 Exomiser eval — baseline (no embeddings, run once)
# -----------------------------------------------------------------------
echo "========================================================"
echo " Track 1: Exomiser baselines"
echo "========================================================"
$PYTHON $REPO_ROOT/eval/exomiser_eval.py \
    --phenotype-data-dir $PHENO_DIR \
    --track 1 \
    --split test

# -----------------------------------------------------------------------
# Step 3: Track 1 Exomiser eval — INDIGENA variants, one per model
# -----------------------------------------------------------------------
for CONFIG in mgi_only hp_only mgi_hp; do
    case $CONFIG in
        mgi_only) EMB=$EMB_MGI ;;
        hp_only)  EMB=$EMB_HPO ;;
        mgi_hp)   EMB=$EMB_MGI_HPO ;;
    esac

    echo "========================================================"
    echo " Track 1: INDIGENA-Exomiser  config=$CONFIG"
    echo "========================================================"
    PRIORITISERS=indigena_hiphive,indigena_phive,indigena_phenix \
    $PYTHON $REPO_ROOT/eval/exomiser_eval.py \
        --phenotype-data-dir $PHENO_DIR \
        --track 1 \
        --split test \
        --embeddings $EMB \
        --tag $CONFIG
done

# -----------------------------------------------------------------------
# Step 4: Track 2 Exomiser eval — INDIGENA variants, one per model
# (--skip-cli reuses existing CLI TSV outputs from the first Track 2 run)
# -----------------------------------------------------------------------
for CONFIG in mgi_only hp_only mgi_hp; do
    case $CONFIG in
        mgi_only) EMB=$EMB_MGI ;;
        hp_only)  EMB=$EMB_HPO ;;
        mgi_hp)   EMB=$EMB_MGI_HPO ;;
    esac

    echo "========================================================"
    echo " Track 2: INDIGENA-Exomiser  config=$CONFIG"
    echo "========================================================"
    $PYTHON $REPO_ROOT/eval/track2_eval.py \
        --phenotype-data-dir $PHENO_DIR \
        --app-props $APP_PROPS \
        --split test \
        --embeddings $EMB \
        --skip-cli \
        --tag $CONFIG
done

echo "All Exomiser 3×3 runs complete."
