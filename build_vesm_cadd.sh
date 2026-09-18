#!/bin/bash
set -uo pipefail
cd ~/Git/vexomiser
MM="$HOME/.local/bin/micromamba run -n embedpvp"
V=/home/zhapacfp/Git/embedpvp2/data/vesm/vesm_3b_hg38_missense.tsv.gz
OUT=exomiser-data/cadd/1.7/vesm_as_cadd.tsv.gz
echo "[$(date +%H:%M)] reformatting VESM -> CADD-format (PHRED=-10log10(1-v))"
$MM bash -c "zcat $V | awk -F'\t' 'BEGIN{OFS=\"\t\"; print \"##VESM-3B-as-CADD\"; print \"#Chrom\",\"Pos\",\"Ref\",\"Alt\",\"RawScore\",\"PHRED\"} {v=\$5+0; if(v>0.9999)v=0.9999; if(v<0)v=0; p=-10*log(1-v)/log(10); c=\$1; sub(/^chr/,\"\",c); print c,\$2,\$3,\$4,v,p}' | bgzip -c > $OUT"
echo "[$(date +%H:%M)] indexing"; $MM tabix -f -s1 -b2 -e2 $OUT
echo "[$(date +%H:%M)] rows=$($MM bash -c \"zcat $OUT | tail -n+3 | wc -l\") contigs=$($MM tabix -l $OUT | wc -l) DONE_VESM"
