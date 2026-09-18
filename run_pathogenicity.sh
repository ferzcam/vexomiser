#!/bin/bash
set -uo pipefail
cd ~/Git/vexomiser
export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
PY=~/miniforge3/envs/indigena/bin/python
EMB=/home/zhapacfp/Git/embedpvp2/data/gda/variant_eval/indigena_mp_embeddings.tsv
C="--phenotype-data-dir exomiser-data/2406_phenotype --app-props exomiser-data/application.properties --split X --embeddings $EMB --workers 16"
echo "[$(date +%H:%M)] setting 4: INDIGENA + CADD"
EXOMISER_ANALYSIS=exomiser-data/analysis_cadd.yml $PY eval/track2_eval.py $C --tag _indigena_mp_cadd --work-dir data/results/track2_jobs_cadd > data/results/setting4_cadd.log 2>&1
cp -f data/results/exomiser_track2_X_summary.txt data/results/summary_cadd.txt
echo "[$(date +%H:%M)] setting 5a: INDIGENA + AlphaMissense"
EXOMISER_ANALYSIS=exomiser-data/analysis_alphamissense.yml $PY eval/track2_eval.py $C --tag _indigena_mp_am --work-dir data/results/track2_jobs_am > data/results/setting5a_am.log 2>&1
cp -f data/results/exomiser_track2_X_summary.txt data/results/summary_am.txt
echo "[$(date +%H:%M)] DONE"
