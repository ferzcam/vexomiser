#!/usr/bin/env python
"""Add depth-k superclass perturbation columns to the noise TSV.

For each k in --depths, every baseline HPO term is REPLACED by its
exact-depth-k HP ancestors (replace-all union, matching the indigena
hierarchy_test build_abstracted_phenotypes.py semantics). A term with no
depth-k ancestor (walked past the root) is kept as-is. Output columns are
named hpo_terms_superclass_k{k}, bare HP ids joined by ";".

Reuses the precomputed ancestor map from the indigena hierarchy_test repo
(keys/values are OBO URIs, exact-depth lists).

  python eval/add_superclass_depth_columns.py --depths 1 2 3
"""
import argparse, json, os
import pandas as pd

OBO = "http://purl.obolibrary.org/obo/"

def id_to_uri(hp):      # HP:0001263 -> http://.../HP_0001263
    return OBO + hp.replace(":", "_")

def uri_to_id(uri):     # http://.../HP_0001263 -> HP:0001263
    return uri.rsplit("/", 1)[-1].replace("_", ":")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", default="data/pavs_with_noise_hpo.tsv")
    ap.add_argument("--ancestors",
                    default=os.path.expanduser(
                        "~/Git/indigena/hierarchy_test/data/hpo_ancestors.json"))
    ap.add_argument("--base-col", default="hpo_terms_baseline")
    ap.add_argument("--depths", type=int, nargs="+", default=[1, 2, 3])
    a = ap.parse_args()

    anc = json.load(open(a.ancestors))["ancestors_by_depth"]
    df = pd.read_csv(a.tsv, sep="\t", dtype=str, keep_default_na=False)

    def abstract(cell, k):
        terms = [t.split("|")[0] for t in str(cell).split(";") if t]
        out = []
        for p in terms:
            levels = anc.get(id_to_uri(p), [])
            ks = levels[k - 1] if len(levels) >= k else []
            if ks:
                out.extend(uri_to_id(u) for u in ks)
            else:
                out.append(p)            # past root / unknown -> keep
        # dedup, preserve first-seen order
        seen, dedup = set(), []
        for t in out:
            if t not in seen:
                seen.add(t); dedup.append(t)
        return ";".join(dedup)

    for k in a.depths:
        col = f"hpo_terms_superclass_k{k}"
        df[col] = df[a.base_col].map(lambda c: abstract(c, k))
        n_changed = (df[col] != df[a.base_col]).sum()
        mean_n = df[col].map(lambda s: len([x for x in s.split(";") if x])).mean()
        print(f"{col}: changed {n_changed}/{len(df)} rows, mean terms/row={mean_n:.2f}")

    df.to_csv(a.tsv, sep="\t", index=False)
    print(f"wrote {a.tsv} with {len(df.columns)} columns")

if __name__ == "__main__":
    main()
