#!/bin/bash
# Pre-filter spiked VCFs: remove REF="." (malformed GIAB insertions) that crash HTSJDK
set -euo pipefail

SRC=/home/zhapacfp/Git/vexomiser/data/spiked_vcfs
DST=/home/zhapacfp/Git/vexomiser/data/spiked_vcfs_filtered
mkdir -p "$DST"

filter_one() {
    local src=$1
    local dst=$2/$(basename "$1")
    bcftools view -e 'REF="."' -O z -o "$dst" "$src"
    bcftools index -t "$dst"
}
export -f filter_one

ls "$SRC"/*.vcf.gz | grep -v '\.tbi' | \
    xargs -P 24 -I{} bash -c 'filter_one "$@"' _ {} "$DST"

n=$(ls "$DST"/*.vcf.gz | grep -v '\.tbi' | wc -l)
echo "Done: $n filtered VCFs in $DST"
