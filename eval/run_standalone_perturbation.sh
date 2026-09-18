#!/usr/bin/env bash
# Standalone INDIGENA on perturbed phenotypes (eval-only, reuses trained checkpoints).
#
# For each of the 3 training configs and 5 perturbation modes, loads the model and
# re-scores the test split with the perturbed query HPO. One run emits all 5 eval
# modes (eval_mgi, eval_hp_mg, eval_mgi_hp_mg, eval_hp_pd, eval_mgi_hp_pd).
# Output: data/results/<file_id>_test_<mode>_<eval_mode>.tsv
#
# Usage (from vexomiser root):
#   nohup bash eval/run_standalone_perturbation.sh > data/results/perturbation_logs/standalone.log 2>&1 &

set -uo pipefail

PYTHON=~/miniforge3/envs/indigena/bin/python
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"
SCRIPT=eval/indigena_train.py

UPHENO=~/Git/indigena/data/upheno_owl2vecstar_edges.tsv
MGI_CSV=~/Git/indigena/data/gene_phenotypes.csv
HOM=~/Git/indigena/data/HOM_MouseHumanSequence.rpt
HPO_G2P=$REPO_ROOT/data/genes_to_phenotype.txt
PHENO_HPOA=$REPO_ROOT/data/phenotype.hpoa
HPO_FILE=$REPO_ROOT/data/pavs_with_noise_hpo.tsv

COMMON="--upheno-edges $UPHENO --hom-file $HOM --phenotype-hpoa $PHENO_HPOA
        --eval-gene-phenotypes $HPO_G2P --graph2 --graph3 --graph4
        --track 1 --eval-split test --only-eval"

MODES=(baseline incomplete_3 noisy_2 imprecise_2 superclass_all)
mkdir -p "$REPO_ROOT/data/results/perturbation_logs"

ts() { date '+%Y-%m-%d %H:%M:%S'; }
run_step() {
    local label="$1"; shift
    echo "[$(ts)] >>> START  $label"
    if "$@"; then echo "[$(ts)] <<< OK     $label";
    else echo "[$(ts)] !!! FAIL   $label (exit $?) — continuing"; fi
}

echo "[$(ts)] ===== Standalone INDIGENA perturbation start ====="
for CFG in mgi hp mgi_hp; do
    case $CFG in
        mgi)    CFG_FLAGS="--mgi-gene-phenotypes $MGI_CSV --no-hpo-fallback" ;;
        hp)     CFG_FLAGS="--hpo-gene-phenotypes $HPO_G2P --no-hpo-fallback" ;;
        mgi_hp) CFG_FLAGS="--mgi-gene-phenotypes $MGI_CSV --hpo-gene-phenotypes $HPO_G2P" ;;
    esac
    for MODE in "${MODES[@]}"; do
        run_step "standalone $CFG $MODE" \
            $PYTHON $SCRIPT $COMMON $CFG_FLAGS \
                --hpo-file "$HPO_FILE" --hpo-column "hpo_terms_$MODE" --tag "$MODE"
    done
done
echo "[$(ts)] ===== Standalone INDIGENA perturbation COMPLETE ====="
