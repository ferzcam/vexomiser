"""
Fix spiked VCFs: the causal variant was spiked with REF="." and ALT=ref_allele
(instead of REF=ref_allele, ALT=pathogenic_allele), which HTSJDK cannot parse.

Strategy: stream the original spiked VCF, replace the single broken record
(REF=".") with the correct REF/ALT from track2_cases.tsv, output to
data/spiked_vcfs_fixed/. No sort or concat needed — the replacement is in-place
so the order is preserved.
"""
import os
import subprocess
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, "data")

ORIG_DIR  = os.path.join(DATA, "spiked_vcfs")
FIXED_DIR = os.path.join(DATA, "spiked_vcfs_fixed")
CASES_TSV = os.path.join(DATA, "track2_cases.tsv")

AWK_SCRIPT = r"""
BEGIN { OFS="\t" }
/^#/ { print; next }
$4 == "." { $4 = ref; $5 = alt; print; next }
{ print }
"""


def fix_one(row):
    cid   = row["case_id"]
    chrom = row["chrom_fixed"]
    pos   = int(float(row["pos"]))
    ref   = row["ref"]
    alt   = row["alt"]

    src = os.path.join(ORIG_DIR, f"{cid}.vcf.gz")
    dst = os.path.join(FIXED_DIR, f"{cid}.vcf.gz")

    if not os.path.exists(src):
        return cid, "missing_src"

    # Decompress → awk replace → bgzip compress
    # The broken record has REF="." — set it to the correct REF/ALT.
    # There is exactly one REF="." record per VCF (the spiked causal variant).
    awk_cmd = (
        f'BEGIN{{OFS="\\t"}} '
        f'/^#/{{print;next}} '
        f'$4=="."{{$4="{ref}";$5="{alt}";print;next}} '
        f'{{print}}'
    )
    try:
        view = subprocess.Popen(
            ["bcftools", "view", "-O", "v", src],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
        )
        awk = subprocess.Popen(
            ["awk", awk_cmd],
            stdin=view.stdout, stdout=subprocess.PIPE
        )
        view.stdout.close()
        compress = subprocess.Popen(
            ["bcftools", "view", "-O", "z", "-o", dst, "-"],
            stdin=awk.stdout, stderr=subprocess.DEVNULL
        )
        awk.stdout.close()

        view.wait()
        awk.wait()
        rc = compress.wait()

        if rc != 0:
            return cid, f"compress_failed:{rc}"

        subprocess.run(
            ["bcftools", "index", "-t", dst],
            check=True, capture_output=True
        )
        return cid, "ok"
    except Exception as e:
        return cid, f"error:{e}"


def main(split: str = "test", workers: int = 24):
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
    errors = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(fix_one, r): r["case_id"] for r in rows}
        for fut in tqdm(as_completed(futs), total=len(futs)):
            cid, status = fut.result()
            if status == "ok":
                ok += 1
            else:
                errors.append((cid, status))

    if errors:
        log.warning(f"{len(errors)} failures:")
        for cid, st in errors[:10]:
            log.warning(f"  {cid}: {st}")

    log.info(f"Done: {ok}/{len(cases)} VCFs fixed in {FIXED_DIR}")


if __name__ == "__main__":
    import sys
    split = sys.argv[1] if len(sys.argv) > 1 else "test"
    main(split=split)
