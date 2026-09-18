#!/usr/bin/env bash
# Track-1-only depth-k superclass perturbation (eval-only; reuses trained embeddings).
# Modes: superclass_k1/k2/k3  x  baselines + 3 INDIGENA configs.
# Usage: nohup bash eval/run_superclass_depth_track1.sh > data/results/perturbation_logs/superclass_depth.log 2>&1 &
set -uo pipefail
PYTHON=~/miniforge3/envs/indigena/bin/python
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$REPO_ROOT"
PHENO_DIR=$REPO_ROOT/exomiser-data/2406_phenotype
HPO_FILE=$REPO_ROOT/data/pavs_with_noise_hpo_superclass_depth.tsv
MODELS=$REPO_ROOT/data/models
LOGDIR=$REPO_ROOT/data/results/perturbation_logs; mkdir -p "$LOGDIR"

EMB_mgi=$MODELS/indigena_transd_track1_graph4_mgi_nofallback_seed0_dim100_bs2048_lr0.001_embeddings.tsv
EMB_hp=$MODELS/indigena_transd_track1_graph4_hpo_seed0_dim100_bs2048_lr0.001_embeddings.tsv
EMB_mgi_hp=$MODELS/indigena_transd_track1_graph4_mgi_hpo_seed0_dim100_bs2048_lr0.001_embeddings.tsv

MODES=(superclass_k1 superclass_k2 superclass_k3)
CONFIGS=(mgi hp mgi_hp)
INDIGENA_SET=indigena_hiphive,indigena_phive,indigena_phenix
SPLIT=test
ts(){ date '+%Y-%m-%d %H:%M:%S'; }
run_step(){ local l="$1"; shift; echo "[$(ts)] >>> START $l"; if "$@"; then echo "[$(ts)] <<< OK $l"; else echo "[$(ts)] !!! FAIL $l (exit $?) — continuing"; fi; }

echo "[$(ts)] ===== superclass depth Track1 start (modes=${MODES[*]}) ====="
for MODE in "${MODES[@]}"; do
    COL=hpo_terms_$MODE
    echo "[$(ts)] ===== MODE $MODE (column $COL) ====="
    run_step "T1 $MODE baselines" \
        $PYTHON eval/exomiser_eval.py --phenotype-data-dir "$PHENO_DIR" --track 1 --split "$SPLIT" \
            --hpo-file "$HPO_FILE" --hpo-column "$COL" --tag "$MODE"
    for CFG in "${CONFIGS[@]}"; do
        EMB_VAR="EMB_$CFG"; EMB=${!EMB_VAR}
        run_step "T1 $MODE indigena/$CFG" \
            env PRIORITISERS=$INDIGENA_SET \
            $PYTHON eval/exomiser_eval.py --phenotype-data-dir "$PHENO_DIR" --track 1 --split "$SPLIT" \
                --embeddings "$EMB" --hpo-file "$HPO_FILE" --hpo-column "$COL" --tag "${MODE}_${CFG}"
    done
done
echo "[$(ts)] ===== superclass depth Track1 COMPLETE ====="
