# Depth-k superclass perturbation — Vexomiser (Track 1, test split)

Mirrors the indigena `hierarchy_test` abstraction study inside the Exomiser integration.
Every query phenotype is **replaced** by its exact-depth-`k` HP ancestors (replace-all,
keep-as-is past root) — identical semantics to `build_abstracted_phenotypes.py`. `k=0` is
the unperturbed baseline. Eval-only (reuses trained embeddings/checkpoints). 518 test cases.

Perturbation columns: `data/pavs_with_noise_hpo_superclass_depth.tsv` (original TSV left pristine).
Build: `eval/add_superclass_depth_columns.py` (reuses `hierarchy_test/data/hpo_ancestors.json`).
Run: `eval/run_superclass_depth_track1.sh` (integrated) + `eval/run_standalone_depth.sh` (standalone).
Report: `data/results/perturbation_report_superclass_depth.md` (existing report untouched).

NB. The legacy `superclass_all` mode (one direct-parent hop) ≈ `k=1`. The collapse only
appears at **k≥2**, which is why the single-level `superclass_all` test showed no collapse.

## A. Standalone INDIGENA (per-disease OMIM scoring) vs IC baselines — MRR

The standalone ranker scores `max_d BMA(patient_HPs, disease_d_HPs)` over OMIM disease
models — the direct analog of the hierarchy_test standalone INDIGENA. IC baselines
(HiPhive, PhenIX) are the Resnik analog (PhenIX = HP-only IC = closest to hierarchy_test Resnik).

| method | k0 | k1 | k2 | k3 | k3 / k0 |
|---|---|---|---|---|---|
| **PhenIX (IC / Resnik analog)** | 0.259 | 0.213 | 0.170 | **0.120** | **46%** |
| HiPhive (IC, best baseline) | 0.274 | 0.227 | 0.148 | 0.071 | 26% |
| Standalone INDIGENA — MGI-only | 0.275 | 0.243 | 0.131 | **0.030** | **11%** |
| Standalone INDIGENA — HP-only | 0.281 | 0.214 | 0.101 | 0.034 | 12% |
| Standalone INDIGENA — MGI+HPO | 0.281 | 0.208 | 0.084 | 0.022 | 8% |

Hits@1 (standalone): MGI-only 0.236→0.201→0.089→**0.004**; HP-only 0.237→0.168→0.050→0.006.

**Reading.** At k0 INDIGENA matches or beats the IC baselines. Under abstraction it
**collapses**: by k3 the best standalone INDIGENA holds 12% of baseline MRR (H@1 → ~0),
while PhenIX retains 46%. INDIGENA is ~3.5× worse than PhenIX at k3 despite being ahead
at k0 — exactly the hierarchy_test pattern, now reproduced in Vexomiser.

## B. INDIGENA as drop-in similarity inside Exomiser — MRR

Same gene–disease scaffold; only the pairwise similarity is swapped. Compare each
INDIGENA-X against its matched IC baseline X.

| method | k0 | k1 | k2 | k3 |
|---|---|---|---|---|
| phenix (IC) | 0.259 | 0.213 | 0.170 | 0.120 |
| hiphive (IC) | 0.274 | 0.227 | 0.148 | 0.071 |
| indigena_phenix/mgi | 0.265 | 0.236 | 0.142 | 0.039 |
| indigena_phenix/hp | 0.272 | 0.208 | 0.109 | 0.038 |
| indigena_hiphive/mgi | 0.266 | 0.234 | 0.139 | 0.036 |
| indigena_hiphive/hp | 0.269 | 0.208 | 0.105 | 0.038 |
| indigena_*/mgi_hp | 0.26 | 0.03 | 0.015 | 0.014 |

**Reading.** The matched INDIGENA variants track their IC baseline at k0–k1 but fall
~3× below PhenIX by k3 — the same collapse, even when INDIGENA only supplies the
pairwise term similarity. The `mgi_hp` config is pathological: it collapses immediately
at k1 (MRR 0.26→0.03), a known artifact of mixing mouse+human gene-phenotype annotations.

## Contrast with indigena/hierarchy_test (human leak-strict Model H)

| | hierarchy_test (standalone, human Model H) | Vexomiser (standalone, this study) |
|---|---|---|
| INDIGENA MRR k0→k3 | 0.765 → 0.049 (6%) | 0.275 → 0.030 (11%) |
| Resnik/IC MRR k0→k3 | 0.844 → 0.561 (66%) | PhenIX 0.259 → 0.120 (46%) |
| Verdict | INDIGENA collapses, Resnik graceful | **same** |

Both studies agree: **INDIGENA's symmetric embedding similarity loses discrimination
under ancestor abstraction, while IC/Resnik degrades gracefully.** The effect is monotonic
in depth and only becomes pronounced at k≥2 — so a one-level (superclass_all) test misses it.
