"""
Compare Resnik BMA standalone vs INDIGENA G4 standalone.

For each test case, computes the rank of the causal gene under each method,
then reports:
  - Table: per-method metrics + rank correlation
  - Plot 1: scatter of per-case ranks (Resnik vs INDIGENA)
  - Plot 2: rank difference distribution
  - Plot 3: CDF of ranks for both methods
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import spearmanr, kendalltau
import torch as th

RESNIK_PATH  = "data/results/semsim_resnik_bma_test.tsv"
INDIGENA_PATH = "data/results/indigena_transd_track1_graph4_seed0_dim100_bs2048_lr0.001_test.tsv"
OUT_PLOT     = "data/results/semsim_vs_indigena.png"
OUT_TABLE    = "data/results/semsim_vs_indigena_table.tsv"


def load_ranks(path):
    """Return (case_ids, causal_genes, ranks) from a results TSV."""
    rows = [l.strip().split('\t') for l in open(path)]
    case_ids, genes, ranks = [], [], []
    n_genes = len(rows[0][3:])
    gene_ids = th.arange(n_genes)

    for row in rows:
        gene      = row[0]
        case_id   = row[1]
        position  = int(row[2])
        scores    = th.tensor([-float(x) for x in row[3:]])

        # tie-break with random permutation (same as metrics.py)
        perm = th.randperm(n_genes)
        updated_pos = th.where(gene_ids[perm] == position)[0].item()
        scores_perm = scores[perm]
        order = th.argsort(scores_perm, descending=False)
        rank  = th.where(order == updated_pos)[0].item() + 1

        case_ids.append(case_id)
        genes.append(gene)
        ranks.append(rank)

    return case_ids, genes, np.array(ranks)


def metrics_from_ranks(ranks, n_genes):
    mr   = ranks.mean()
    mrr  = (1 / ranks).mean()
    h1   = (ranks <= 1).mean()
    h3   = (ranks <= 3).mean()
    h10  = (ranks <= 10).mean()
    h100 = (ranks <= 100).mean()
    return dict(MR=mr, MRR=mrr, **{"Hits@1": h1, "Hits@3": h3,
                                    "Hits@10": h10, "Hits@100": h100})


# ── Load ────────────────────────────────────────────────────────────────────
print("Loading Resnik BMA results…")
case_ids_r, genes_r, ranks_r = load_ranks(RESNIK_PATH)

print("Loading INDIGENA G4 results…")
case_ids_i, genes_i, ranks_i = load_ranks(INDIGENA_PATH)

# Align on case_id (both should be same 518 cases)
df_r = pd.DataFrame({"case_id": case_ids_r, "gene": genes_r, "rank_resnik": ranks_r})
df_i = pd.DataFrame({"case_id": case_ids_i, "gene": genes_i, "rank_indigena": ranks_i})
df   = df_r.merge(df_i, on="case_id")
print(f"Aligned cases: {len(df)}")

n_genes = len(open(RESNIK_PATH).readline().split('\t')) - 3

# ── Table ───────────────────────────────────────────────────────────────────
m_r = metrics_from_ranks(df["rank_resnik"].values,  n_genes)
m_i = metrics_from_ranks(df["rank_indigena"].values, n_genes)

spear, spear_p = spearmanr(df["rank_resnik"], df["rank_indigena"])
kend,  kend_p  = kendalltau(df["rank_resnik"], df["rank_indigena"])

print("\n── Metrics ───────────────────────────────────────────────")
print(f"{'Method':<20} {'MR':>7} {'MRR':>7} {'H@1':>7} {'H@3':>7} {'H@10':>7} {'H@100':>8}")
for label, m in [("Resnik BMA", m_r), ("INDIGENA G4", m_i)]:
    print(f"{label:<20} {m['MR']:>7.1f} {m['MRR']:>7.3f} {m['Hits@1']:>7.3f} "
          f"{m['Hits@3']:>7.3f} {m['Hits@10']:>7.3f} {m['Hits@100']:>8.3f}")

print(f"\nSpearman rank correlation: ρ = {spear:.3f}  (p = {spear_p:.2e})")
print(f"Kendall τ:                 τ = {kend:.3f}  (p = {kend_p:.2e})")

# Write table
tbl = pd.DataFrame({
    "method":   ["Resnik BMA", "INDIGENA G4"],
    "MR":       [m_r["MR"],   m_i["MR"]],
    "MRR":      [m_r["MRR"],  m_i["MRR"]],
    "Hits@1":   [m_r["Hits@1"],  m_i["Hits@1"]],
    "Hits@3":   [m_r["Hits@3"],  m_i["Hits@3"]],
    "Hits@10":  [m_r["Hits@10"], m_i["Hits@10"]],
    "Hits@100": [m_r["Hits@100"],m_i["Hits@100"]],
    "spearman_rho": [spear, ""],
    "kendall_tau":  [kend,  ""],
})
tbl.to_csv(OUT_TABLE, sep='\t', index=False)
print(f"Table saved → {OUT_TABLE}")

# ── Plot ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(15, 5))
gs  = gridspec.GridSpec(1, 3, figure=fig)

# 1. Scatter: rank Resnik vs rank INDIGENA
ax1 = fig.add_subplot(gs[0])
ax1.scatter(df["rank_resnik"], df["rank_indigena"],
            alpha=0.35, s=12, color="steelblue", edgecolors="none")
lim = max(df["rank_resnik"].max(), df["rank_indigena"].max()) + 50
ax1.plot([0, lim], [0, lim], "k--", lw=0.8, alpha=0.5)
ax1.set_xlabel("Rank (Resnik BMA)")
ax1.set_ylabel("Rank (INDIGENA G4)")
ax1.set_title(f"Per-case ranks\nSpearman ρ = {spear:.3f}")
ax1.set_xlim(0, lim); ax1.set_ylim(0, lim)

# 2. Rank difference distribution
ax2 = fig.add_subplot(gs[1])
diff = df["rank_resnik"].values - df["rank_indigena"].values
ax2.hist(diff, bins=50, color="coral", edgecolor="white", linewidth=0.4)
ax2.axvline(0, color="black", lw=1, ls="--")
ax2.set_xlabel("Rank(Resnik) − Rank(INDIGENA)")
ax2.set_ylabel("Cases")
ax2.set_title(f"Rank difference\n(positive = Resnik ranks gene higher)")

# 3. CDF of ranks
ax3 = fig.add_subplot(gs[2])
for ranks, label, color in [
    (df["rank_resnik"].values,  "Resnik BMA",  "steelblue"),
    (df["rank_indigena"].values, "INDIGENA G4", "darkorange"),
]:
    sorted_r = np.sort(ranks)
    cdf = np.arange(1, len(sorted_r) + 1) / len(sorted_r)
    ax3.plot(sorted_r, cdf, label=label, color=color, lw=1.8)
ax3.set_xlabel("Rank of causal gene")
ax3.set_ylabel("Fraction of cases")
ax3.set_title("CDF of causal gene rank")
ax3.legend(fontsize=9)
ax3.set_xlim(0, n_genes)

plt.tight_layout()
plt.savefig(OUT_PLOT, dpi=150, bbox_inches="tight")
print(f"Plot saved → {OUT_PLOT}")
