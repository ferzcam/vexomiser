#!/usr/bin/env python
"""Generate perturbed-HPO phenopackets for Track 2 Exomiser CLI runs.

Track 2's hiphive score is produced by the Exomiser CLI, which reads the patient's
HPO terms from the phenopacket JSON on disk (not from the in-process hpo_list). To
make hiphive reflect a perturbed HPO set, we copy each baseline phenopacket and
replace its included HP phenotypicFeatures with the terms from a chosen perturbation
column of the noise TSV (e.g. hpo_terms_noisy_2).

Only id / subject / included-HP phenotypicFeatures are consumed downstream
(track2_eval.make_phenopacket_v1), but we copy the whole file for fidelity. Labels
are filled from the noise file's `hpo_terms` column where known, else set to the
term id itself (Exomiser keys on the id, so a non-empty placeholder is harmless).

Usage:
    python eval/generate_perturbed_phenopackets.py \
        --phenopacket-dir data/pavs/phenopackets \
        --hpo-file data/pavs_with_noise_hpo.tsv \
        --hpo-column hpo_terms_noisy_2 \
        --out-dir data/results/perturbed_phenopackets/noisy_2
"""
import glob
import json
import os

import click as ck
import pandas as pd


def build_label_map(noise: pd.DataFrame) -> dict:
    """HP id -> label, harvested from the 'id|label;...' hpo_terms column."""
    labels = {}
    for s in noise["hpo_terms"].dropna():
        for tok in str(s).split(";"):
            if "|" in tok:
                hid, lab = tok.split("|", 1)
                labels[hid] = lab
    return labels


@ck.command()
@ck.option("--phenopacket-dir", required=True,
           help="Source directory of baseline phenopacket JSONs.")
@ck.option("--hpo-file", required=True,
           help="Noise TSV keyed by case_id (e.g. data/pavs_with_noise_hpo.tsv).")
@ck.option("--hpo-column", required=True,
           help="Perturbation column to inject, e.g. 'hpo_terms_noisy_2'.")
@ck.option("--out-dir", required=True,
           help="Output directory for perturbed phenopackets.")
def main(phenopacket_dir, hpo_file, hpo_column, out_dir):
    noise = pd.read_csv(hpo_file, sep="\t", dtype=str)
    label_map = build_label_map(noise)
    col = noise.set_index("case_id")[hpo_column]

    os.makedirs(out_dir, exist_ok=True)
    n_written = n_empty = 0
    for path in glob.glob(os.path.join(phenopacket_dir, "*.json")):
        cid = os.path.basename(path)[:-5]
        with open(path) as fh:
            d = json.load(fh)

        terms = []
        if cid in col.index and pd.notna(col.loc[cid]):
            terms = [t.split("|")[0] for t in str(col.loc[cid]).split(";") if t]

        d["phenotypicFeatures"] = [
            {"type": {"id": t, "label": label_map.get(t, t)}, "excluded": False}
            for t in terms
        ]
        with open(os.path.join(out_dir, os.path.basename(path)), "w") as fh:
            json.dump(d, fh)
        n_written += 1
        if not terms:
            n_empty += 1

    print(f"Wrote {n_written} phenopackets to {out_dir} "
          f"({n_empty} with empty HPO set)")


if __name__ == "__main__":
    main()
