#!/usr/bin/env bash
# Run this on the workstation after 2406_hg38.zip has finished downloading.
# It unzips the genome data, then launches the Track 2 evaluation.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA_DIR="$REPO_ROOT/exomiser-data"
LOG_DIR="$REPO_ROOT/data/results"
PYTHON="$HOME/miniforge3/envs/indigena/bin/python"

mkdir -p "$LOG_DIR"

echo "=== Unzipping 2406_hg38.zip ==="
cd "$DATA_DIR"
unzip -q 2406_hg38.zip
echo "Unzip complete. Contents:"
ls -lh 2406_hg38/

echo ""
echo "=== Launching Track 2 evaluation ==="
cd "$REPO_ROOT"

nohup "$PYTHON" eval/track2_eval.py \
    --phenotype-data-dir exomiser-data/2406_phenotype \
    --app-props exomiser-data/application.properties \
    --split test \
    --workers 20 \
    --embeddings data/models/indigena_track1_graph4_embeddings.tsv \
    > "$LOG_DIR/track2_eval.log" 2>&1 &

echo "PID $! — tail the log with:"
echo "  tail -f $LOG_DIR/track2_eval.log"
