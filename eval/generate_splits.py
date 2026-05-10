"""
Generate stratified train/val/test splits for the PAVS benchmark.

Split is disease-disjoint: all cases sharing the same disease ID go to the
same partition. Cases without a MONDO disease ID receive a unique synthetic ID
(SYNTHETIC:<case_id>) so each is treated as its own disease group.

Stratified by cohort (Saudi, DDD, Mixed) to preserve population distribution.

Ratio: 0.8 / 0.1 / 0.1

Output: data/splits/{train,val,test}.tsv  (same columns as track1_cases.tsv)
"""

import os
import random
import click as ck
import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(REPO_ROOT, "data")
SPLITS_DIR = os.path.join(DATA_DIR, "splits")


def assign_disease_id(cases: pd.DataFrame) -> pd.DataFrame:
    """Fill missing disease_id with a unique synthetic ID per case."""
    cases = cases.copy()
    missing = cases["disease_id"].isna() | (cases["disease_id"] == "")
    cases.loc[missing, "disease_id"] = "SYNTHETIC:" + cases.loc[missing, "case_id"]
    return cases


@ck.command()
@ck.option("--track", type=ck.Choice(["1", "2"]), default="1", show_default=True)
@ck.option("--seed", default=42, show_default=True)
@ck.option("--train", "train_ratio", default=0.8, show_default=True)
@ck.option("--val", "val_ratio", default=0.1, show_default=True)
def main(track, seed, train_ratio, val_ratio):
    random.seed(seed)

    cases = pd.read_csv(os.path.join(DATA_DIR, f"track{track}_cases.tsv"), sep="\t")
    cases = assign_disease_id(cases)
    print(f"Track {track}: {len(cases)} total cases")
    print(f"  Unique diseases: {cases['disease_id'].nunique()}")

    splits = {"train": [], "val": [], "test": []}

    for cohort, group in cases.groupby("cohort"):
        # Group by disease — all cases of the same disease go to the same split
        diseases = list(group.groupby("disease_id"))
        random.shuffle(diseases)

        n = len(diseases)
        n_train = int(n * train_ratio)
        n_val = int(n * val_ratio)

        def collect(disease_groups):
            return pd.concat([g for _, g in disease_groups], ignore_index=True)

        splits["train"].append(collect(diseases[:n_train]))
        splits["val"].append(collect(diseases[n_train:n_train + n_val]))
        splits["test"].append(collect(diseases[n_train + n_val:]))

        n_test = n - n_train - n_val
        print(f"  {cohort}: {n_train} / {n_val} / {n_test} diseases "
              f"({len(group)} cases)")

    os.makedirs(SPLITS_DIR, exist_ok=True)
    for name, parts in splits.items():
        df = pd.concat(parts).reset_index(drop=True)
        path = os.path.join(SPLITS_DIR, f"{name}.tsv")
        df.to_csv(path, sep="\t", index=False)
        print(f"Saved {len(df)} cases → {path}")


if __name__ == "__main__":
    main()
