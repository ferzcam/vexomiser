#!/usr/bin/env python
"""Build the phenotype-perturbation robustness comparison.

Parses the result TSVs produced by eval/run_perturbation.sh, computes ranking
metrics for every (track, method, config, perturbation-mode) cell, and emits:
  - data/results/perturbation_analysis.tsv   (long-form, all metrics)
  - data/results/perturbation_report.md      (headline tables + drop-from-baseline)

Only files tagged with one of the five perturbation modes are considered, so the
older 3x3 outputs in the same directory are ignored. For Track 2, baseline methods
(hiphive/phive/phenix) are taken from the canonical no-config files; the config-
tagged duplicates emitted by the per-config indigena passes are skipped.
"""
import glob
import os
import re

import pandas as pd

from metrics import compute_metrics

RESULTS_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "results"
)
MODES = ["baseline", "incomplete_3", "noisy_2", "imprecise_2", "superclass_all"]
CONFIGS = ["mgi", "hp", "mgi_hp"]
BASELINE_METHODS = ["hiphive", "phive", "phenix"]
INDIGENA_METHODS = ["indigena_hiphive", "indigena_phive", "indigena_phenix"]
HEADLINE = ["hits@1", "hits@10", "mrr"]

FNAME_RE = re.compile(r"^exomiser_(?P<method>.+?)_track(?P<track>[12])_test_(?P<rest>.+)$")


def parse_name(stem: str):
    """Return (method, track, mode, config) or None if not a perturbation file."""
    m = FNAME_RE.match(stem)
    if not m:
        return None
    method, track, rest = m["method"], int(m["track"]), m["rest"]
    for mode in MODES:
        if rest == mode:
            config = None
            break
        if rest.startswith(mode + "_") and rest[len(mode) + 1:] in CONFIGS:
            config = rest[len(mode) + 1:]
            break
    else:
        return None  # old 3x3 file (no mode token)

    # Keep baselines only from no-config files; indigena only from config files.
    if method in BASELINE_METHODS and config is not None:
        return None
    if method in INDIGENA_METHODS and config is None:
        return None
    if method not in BASELINE_METHODS + INDIGENA_METHODS:
        return None
    return method, track, mode, config


def main():
    rows = []
    for path in sorted(glob.glob(os.path.join(RESULTS_DIR, "exomiser_*_track*_test_*.tsv"))):
        parsed = parse_name(os.path.basename(path)[:-4])
        if parsed is None:
            continue
        method, track, mode, config = parsed
        try:
            _, macro = compute_metrics(path)
        except Exception as e:
            print(f"WARN: failed on {os.path.basename(path)}: {e}")
            continue
        label = method if config is None else f"{method}/{config}"
        rows.append({
            "track": track, "method": method, "config": config or "-",
            "label": label, "mode": mode,
            **{k: macro[k] for k in ["hits@1", "hits@3", "hits@10", "hits@100", "mrr", "mr", "auc"]},
        })

    df = pd.DataFrame(rows)
    df["mode"] = pd.Categorical(df["mode"], categories=MODES, ordered=True)
    out_tsv = os.path.join(RESULTS_DIR, "perturbation_analysis.tsv")
    df.sort_values(["track", "label", "mode"]).to_csv(out_tsv, sep="\t", index=False)
    print(f"Wrote {out_tsv}  ({len(df)} rows)")

    # Stable row ordering: baselines first, then indigena by config.
    def row_order(label):
        base = label.split("/")[0]
        cfg = label.split("/")[1] if "/" in label else ""
        return (base not in BASELINE_METHODS, BASELINE_METHODS.index(base) if base in BASELINE_METHODS
                else INDIGENA_METHODS.index(base), CONFIGS.index(cfg) if cfg in CONFIGS else -1)

    lines = ["# Phenotype-perturbation robustness\n"]
    lines.append("Macro metrics on the **test** split. Columns are perturbation modes; "
                 "`baseline` is the unperturbed control.\n")
    for track in [1, 2]:
        tdf = df[df["track"] == track]
        if tdf.empty:
            continue
        track_name = "Track 1 (phenotype-only)" if track == 1 else "Track 2 (phenotype+genotype)"
        lines.append(f"\n## {track_name}\n")
        for metric in HEADLINE:
            piv = tdf.pivot_table(index="label", columns="mode", values=metric, observed=False)
            piv = piv.reindex(sorted(piv.index, key=row_order))
            lines.append(f"\n### {metric}\n")
            lines.append("| method | " + " | ".join(MODES) + " |")
            lines.append("|" + "---|" * (len(MODES) + 1))
            for label, r in piv.iterrows():
                cells = " | ".join(f"{r[m]:.3f}" if pd.notna(r[m]) else "—" for m in MODES)
                lines.append(f"| {label} | {cells} |")

        # Drop from baseline for the headline metric (hits@1).
        metric = "hits@1"
        piv = tdf.pivot_table(index="label", columns="mode", values=metric, observed=False)
        piv = piv.reindex(sorted(piv.index, key=row_order))
        lines.append(f"\n### {metric}: absolute drop vs baseline (negative = worse)\n")
        lines.append("| method | " + " | ".join(m for m in MODES if m != "baseline") + " |")
        lines.append("|" + "---|" * len(MODES))
        for label, r in piv.iterrows():
            b = r["baseline"]
            cells = " | ".join(
                f"{r[m] - b:+.3f}" if pd.notna(r[m]) and pd.notna(b) else "—"
                for m in MODES if m != "baseline"
            )
            lines.append(f"| {label} | {cells} |")

    report = os.path.join(RESULTS_DIR, "perturbation_report.md")
    with open(report, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote {report}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
