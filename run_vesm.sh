#!/bin/bash
set -uo pipefail
cd ~/Git/vexomiser
while ! grep -q DONE_VESM data/results/build_vesm.out 2>/dev/null; do pgrep -f build_vesm_cadd >/dev/null || { echo "build died"; exit 1; }; sleep 15; done
export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
PY=~/miniforge3/envs/indigena/bin/python
EMB=/home/zhapacfp/Git/embedpvp2/data/gda/variant_eval/indigena_mp_embeddings.tsv
rm -rf data/results/track2_jobs_vesm
echo "[$(date +%H:%M)] INDIGENA + VESM (CADD channel, 4 workers)"
EXOMISER_ANALYSIS=exomiser-data/analysis_cadd.yml $PY eval/track2_eval.py \
  --phenotype-data-dir exomiser-data/2406_phenotype --app-props exomiser-data/application_vesm.properties \
  --split X --embeddings $EMB --workers 4 --tag _indigena_mp_vesm --work-dir data/results/track2_jobs_vesm \
  > data/results/setting5b_vesm.log 2>&1
cp -f data/results/exomiser_track2_X_summary.txt data/results/summary_vesm.txt
echo "[$(date +%H:%M)] genes.tsv=$(find data/results/track2_jobs_vesm -name '*.genes.tsv'|wc -l)/284 DONE_VESM_RUN"
