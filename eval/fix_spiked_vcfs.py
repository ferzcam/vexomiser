"""
Fix spiked VCFs: the causal variant was spiked with REF="." and ALT=ref_allele
(instead of REF=ref_allele, ALT=pathogenic_allele), which HTSJDK cannot parse.

This script recreates a fixed set of VCFs in data/spiked_vcfs_fixed/ by:
  1. Starting with the filtered GIAB background VCF (no REF="." records).
  2. Appending a properly-formatted causal variant record from track2_cases.tsv.
  3. Sorting, bgzip-compressing, and tabix-indexing the result.

The fixed VCFs can then be used directly with Exomiser.
"""
import os
import subprocess
import tempfile
import pandas as pd
from tqdm import tqdm

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, "data")

FILTERED_DIR = os.path.join(DATA, "spiked_vcfs_filtered")
FIXED_DIR    = os.path.join(DATA, "spiked_vcfs_fixed")
CASES_TSV    = os.path.join(DATA, "track2_cases.tsv")


def vcf_sample_name(vcf_gz: str) -> str:
    out = subprocess.check_output(
        ["bcftools", "view", "-h", vcf_gz], text=True
    )
    for line in out.splitlines():
        if line.startswith("#CHROM"):
            return line.strip().split("\t")[-1]
    raise ValueError(f"No #CHROM line in {vcf_gz}")


def fix_one(row):
    cid    = row["case_id"]
    chrom  = row["chrom_fixed"]
    pos    = int(float(row["pos"]))
    ref    = row["ref"]
    alt    = row["alt"]
    gt     = row["gt"]

    src = os.path.join(FILTERED_DIR, f"{cid}.vcf.gz")
    dst = os.path.join(FIXED_DIR,    f"{cid}.vcf.gz")

    if not os.path.exists(src):
        return cid, "missing_src"

    sample = vcf_sample_name(src)

    # Build mini-VCF reusing the source header (so contig/FORMAT defs are present)
    hdr = subprocess.check_output(
        ["bcftools", "view", "-h", src], text=True
    )
    # Ensure GT FORMAT line is in the header
    fmt_line = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    if "##FORMAT=<ID=GT" not in hdr:
        hdr = hdr.replace("#CHROM\t", fmt_line + "#CHROM\t")

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".vcf", delete=False
    ) as tmp:
        tmp.write(hdr)
        tmp.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t100\tPASS\t.\tGT\t{gt}\n")
        tmp_path = tmp.name

    # Compress + index the mini-VCF using bcftools (bgzip not required)
    tmp_gz = tmp_path + ".gz"
    subprocess.run(
        ["bcftools", "view", "-O", "z", "-o", tmp_gz, tmp_path],
        check=True
    )
    subprocess.run(["bcftools", "index", "-t", tmp_gz], check=True)

    # Merge (GIAB background) + (causal variant)
    subprocess.run(
        ["bcftools", "concat", "--allow-overlaps", "-a",
         "-O", "z", "-o", dst, src, tmp_gz],
        check=True
    )
    subprocess.run(["bcftools", "index", "-t", dst], check=True)

    # Clean up temp files
    for f in [tmp_path, tmp_gz, tmp_gz + ".tbi"]:
        try:
            os.unlink(f)
        except OSError:
            pass

    return cid, "ok"


def main(split: str = "test", workers: int = 24):
    import logging
    from concurrent.futures import ProcessPoolExecutor, as_completed

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    log = logging.getLogger(__name__)

    os.makedirs(FIXED_DIR, exist_ok=True)

    cases = pd.read_csv(CASES_TSV, sep="\t")
    cases = cases[cases["has_vcf"] == True].reset_index(drop=True)

    if split != "all":
        split_path = os.path.join(DATA, "splits", f"{split}.tsv")
        split_ids  = set(pd.read_csv(split_path, sep="\t")["case_id"])
        cases = cases[cases["case_id"].isin(split_ids)].reset_index(drop=True)

    log.info(f"Fixing {len(cases)} spiked VCFs → {FIXED_DIR}")

    rows = [row for _, row in cases.iterrows()]
    ok = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(fix_one, r): r["case_id"] for r in rows}
        for fut in tqdm(as_completed(futs), total=len(futs)):
            cid, status = fut.result()
            if status == "ok":
                ok += 1
            else:
                log.warning(f"{cid}: {status}")

    log.info(f"Done: {ok}/{len(cases)} VCFs fixed in {FIXED_DIR}")


if __name__ == "__main__":
    import sys
    split = sys.argv[1] if len(sys.argv) > 1 else "test"
    main(split=split)
