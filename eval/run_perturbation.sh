#!/usr/bin/env bash
# Phenotype-perturbation robustness matrix.
#
# Runs every method on each perturbed HPO set from data/pavs_with_noise_hpo.tsv to
# measure how ranking quality degrades as phenotypes become incomplete / noisy /
# imprecise. Five perturbation modes x both tracks x three INDIGENA configs.
#
# Per mode:
#   Track 1 (phenotype-only, exomiser_eval.py):
#       - baselines (hiphive, phive, phenix)           tag=<mode>
#       - indigena_* for each config                   tag=<mode>_<config>
#   Track 2 (phenotype+genotype, track2_eval.py):
#       - perturbed modes: regenerate phenopackets + fresh Exomiser CLI so the
#         CLI-derived hiphive reflects the perturbed HPO; baseline reuses the cache.
#       - baselines (hiphive, phive, phenix)           tag=<mode>
#       - indigena_* for each config (--skip-cli)      tag=<mode>_<config>
#
# Robust to single-step failures: each step is logged and the matrix continues.
#
# Usage (from vexomiser root):
#   nohup bash eval/run_perturbation.sh > data/results/perturbation_logs/run.log 2>&1 &

set -uo pipefail

PYTHON=~/miniforge3/envs/indigena/bin/python
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

PHENO_DIR=$REPO_ROOT/exomiser-data/2406_phenotype
APP_PROPS=$REPO_ROOT/exomiser-data/application.properties
HPO_FILE=$REPO_ROOT/data/pavs_with_noise_hpo.tsv
BASE_PP=$REPO_ROOT/data/pavs/phenopackets
MODELS=$REPO_ROOT/data/models
RESULTS=$REPO_ROOT/data/results
PP_OUT=$RESULTS/perturbed_phenopackets
LOGDIR=$RESULTS/perturbation_logs
mkdir -p "$PP_OUT" "$LOGDIR"

EMB_mgi=$MODELS/indigena_transd_track1_graph4_mgi_nofallback_seed0_dim100_bs2048_lr0.001_embeddings.tsv
EMB_hp=$MODELS/indigena_transd_track1_graph4_hpo_seed0_dim100_bs2048_lr0.001_embeddings.tsv
EMB_mgi_hp=$MODELS/indigena_transd_track1_graph4_mgi_hpo_seed0_dim100_bs2048_lr0.001_embeddings.tsv

MODES=(baseline incomplete_3 noisy_2 imprecise_2 superclass_all)
CONFIGS=(mgi hp mgi_hp)
INDIGENA_SET=indigena_hiphive,indigena_phive,indigena_phenix
SPLIT=test

ts() { date '+%Y-%m-%d %H:%M:%S'; }

# run_step <label> <command...> : log, run, continue on failure.
run_step() {
    local label="$1"; shift
    echo "[$(ts)] >>> START  $label"
    if "$@"; then
        echo "[$(ts)] <<< OK     $label"
    else
        echo "[$(ts)] !!! FAIL   $label (exit $?) — continuing"
    fi
}

echo "[$(ts)] ===== Perturbation matrix start ====="
echo "[$(ts)] modes=${MODES[*]}  configs=${CONFIGS[*]}  split=$SPLIT"

mode_i=0
for MODE in "${MODES[@]}"; do
    mode_i=$((mode_i + 1))
    COL=hpo_terms_$MODE
    echo "[$(ts)] ===== MODE $mode_i/${#MODES[@]}: $MODE (column $COL) ====="

    # ------------------------------------------------------------------
    # Track 1: phenotype-only
    # ------------------------------------------------------------------
    run_step "T1 $MODE baselines" \
        $PYTHON eval/exomiser_eval.py \
            --phenotype-data-dir "$PHENO_DIR" --track 1 --split "$SPLIT" \
            --hpo-file "$HPO_FILE" --hpo-column "$COL" --tag "$MODE"

    for CFG in "${CONFIGS[@]}"; do
        EMB_VAR="EMB_$CFG"; EMB=${!EMB_VAR}
        run_step "T1 $MODE indigena/$CFG" \
            env PRIORITISERS=$INDIGENA_SET \
            $PYTHON eval/exomiser_eval.py \
                --phenotype-data-dir "$PHENO_DIR" --track 1 --split "$SPLIT" \
                --embeddings "$EMB" \
                --hpo-file "$HPO_FILE" --hpo-column "$COL" --tag "${MODE}_${CFG}"
    done

    # ------------------------------------------------------------------
    # Track 2: phenotype+genotype
    # ------------------------------------------------------------------
    if [ "$MODE" = "baseline" ]; then
        # Baseline phenopackets == cached CLI inputs: reuse the existing cache.
        PP_DIR=$BASE_PP
        WORK=$RESULTS/track2_jobs
        SKIP_CLI=(--skip-cli)
    else
        # Regenerate phenopackets with perturbed HPO, run the CLI fresh.
        PP_DIR=$PP_OUT/$MODE
        WORK=$RESULTS/track2_jobs_$MODE
        SKIP_CLI=()
        run_step "T2 $MODE gen-phenopackets" \
            $PYTHON eval/generate_perturbed_phenopackets.py \
                --phenopacket-dir "$BASE_PP" --hpo-file "$HPO_FILE" \
                --hpo-column "$COL" --out-dir "$PP_DIR"
    fi

    # Pass 0: baselines (no embeddings). Runs/refreshes the CLI for this mode.
    run_step "T2 $MODE baselines (CLI)" \
        $PYTHON eval/track2_eval.py \
            --phenotype-data-dir "$PHENO_DIR" --app-props "$APP_PROPS" \
            --split "$SPLIT" --phenopacket-dir "$PP_DIR" --work-dir "$WORK" \
            "${SKIP_CLI[@]}" \
            --hpo-file "$HPO_FILE" --hpo-column "$COL" --tag "$MODE"

    # Config passes: indigena variants, reusing the CLI outputs just produced.
    for CFG in "${CONFIGS[@]}"; do
        EMB_VAR="EMB_$CFG"; EMB=${!EMB_VAR}
        run_step "T2 $MODE indigena/$CFG" \
            $PYTHON eval/track2_eval.py \
                --phenotype-data-dir "$PHENO_DIR" --app-props "$APP_PROPS" \
                --split "$SPLIT" --phenopacket-dir "$PP_DIR" --work-dir "$WORK" \
                --skip-cli --embeddings "$EMB" \
                --hpo-file "$HPO_FILE" --hpo-column "$COL" --tag "${MODE}_${CFG}"
    done
done

echo "[$(ts)] ===== Perturbation matrix COMPLETE ====="
