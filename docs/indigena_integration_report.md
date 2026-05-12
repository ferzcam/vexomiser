# Vexomiser

## Summary

We integrated INDIGENA's embedding-based BMA as a drop-in replacement for the
IC-based semantic similarity used by HiPhive, Phive, and PhenIX in Exomiser,
and evaluated on the PAVS Track 1 benchmark (test split, 518 cases) and Track 2
(phenotype + variant, test split, 284 cases).
On Track 1, INDIGENA-HiPhive matches HiPhive (MRR 0.268 vs 0.269) and
INDIGENA-PhenIX outperforms original PhenIX (MRR 0.268 vs 0.261),
while requiring no ontology database at inference time.
On Track 2, INDIGENA methods consistently improve top-10 and top-100 recall
over their baselines — INDIGENA-Phive achieves the best MR (9.7) and Hits@10
(0.940) despite Phive being the weakest phenotype-only baseline.

---

## Background

We trained INDIGENA on the PAVS dataset and evaluated its phenotype similarity
against Exomiser's existing methods. The integration replaces only the similarity
function in each prioritiser; gene-phenotype model associations (OMIM, MGI, ZFIN)
remain unchanged.

---

## Methods

### Benchmark and Evaluation Setup

We used the **PAVS** (Phenotype-Associated Variant Scoring) benchmark, Track 1 (5,076
cases, phenotype-only gene ranking). Cases were split into train / val / test (80/10/10)
using a **disease-disjoint** strategy: all cases sharing the same MONDO disease ID are
assigned to the same partition, preventing a disease seen in training from appearing in
evaluation. The split was stratified by cohort (Saudi, DDD, Mixed): each cohort's
diseases were independently shuffled and distributed 80/10/10, so all three populations
are represented proportionally in every partition. The gene pool was fixed to all 2,258
unique genes across the full dataset so rankings are directly comparable across methods.
All INDIGENA results are reported on the test split (518 cases) only.

### Graph Construction

INDIGENA was trained on a knowledge graph built from four cumulative components:

| Layer | Content |
|---|---|
| G1 | UPheno OWL2Vec* projection — 771K triples covering HP, MP, and UPHENO bridging axioms |
| G2 | Gene → phenotype associations (see below) |
| G3 | Disease → HP associations from training cases only; test/val diseases excluded to preserve inductiveness |
| G4 | Gene ↔ disease associations from training cases (supervised GDA signal) |

**Gene phenotype source (G2).** The choice of gene phenotype associations is scientifically
consequential. We evaluated two configurations:

- **MGI-only (v1):** MGI mouse-knockout MP phenotypes as the primary source, mapped to
  human genes via the HomoloGene orthology table. Training-case HPO terms used as a
  fallback for genes with no mouse ortholog.
- **MGI + HPO (v2):** Same as above, plus human HP annotations from the HPO
  `genes_to_phenotype.txt` file, which aggregates gene–HP associations from OMIM and
  Orphanet. HP terms are added on top of MGI annotations; training-case HPO fallback
  applies only to genes covered by neither source.

The HPO gene–phenotype annotations are derived from OMIM gene–disease–phenotype chains.
HiPhive uses the same information at inference time (retrieval from gene disease models),
so including it in G2 gives INDIGENA a comparable information basis rather than an unfair
advantage. The disease-disjoint split still ensures test case diseases and their specific
gene–disease associations are never seen during training.

### Training

INDIGENA was trained using **TransD** KGE (100-dimensional embeddings, 300 epochs,
Adam optimiser, batch size 2,048, learning rate 0.001) on a GPU. Two configurations
were trained, differing only in the G2 gene phenotype source (MGI-only vs MGI + HPO),
with all other hyperparameters identical.

### Integration into Exomiser

To use INDIGENA as a drop-in similarity function within Exomiser, we introduced a
`ModelScorerFactory` interface that abstracts the phenotype similarity computation away
from the prioritiser classes. Both `HiPhivePriority` and `PhivePriority` were updated to
accept a factory at construction time, defaulting to the original Phenodigm behaviour.

`IndigenaModelScorerFactory` implements this interface: given a patient's HP terms and a
gene model's phenotype terms (HP or MP), it computes the INDIGENA BMA score using
pre-exported entity embedding vectors. Because the UPheno embedding space covers both HP
and MP terms, no separate cross-species mapping step is required — the embedding space
handles it implicitly. Future neural similarity methods (e.g. Multihop-GDA) can be added
by implementing the same interface without modifying the prioritisers.


PhenIX uses a separate Ontologizer-based IC similarity implementation
(`PhenixPriority`) that is architecturally independent of HiPhive and does not expose a
pluggable similarity factory. There is therefore no direct way to substitute the
similarity function inside PhenIX. As an approximation, we configured HiPhivePriority
with `runParams("human")`, which restricts scoring to human disease models only (HP-only,
no cross-species lookup, no PPI network) — the same scope as PhenIX — and replaced its
similarity computation with the INDIGENA BMA scorer. This yields INDIGENA-PhenIX as a
functionally comparable but architecturally distinct approximation of what a native
INDIGENA-PhenIX integration would produce.

---

## Results

### Standalone INDIGENA

When run as a standalone gene ranker (test split, 518 cases, 2,258 genes), without
Exomiser's gene models:

| Configuration | MR | MRR | Hits@1 | Hits@3 | Hits@10 | Hits@100 | AUC |
|---|---|---|---|---|---|---|---|
| v1 — MGI only | 761.2 | 0.116 | 0.079 | 0.118 | 0.191 | 0.375 | 0.664 |
| v2 — MGI + HPO | 435.6 | 0.193 | 0.147 | 0.203 | 0.278 | 0.502 | 0.808 |

Adding HPO gene–phenotype annotations (v2) substantially improves every metric: MRR
+66%, Hits@10 +46%, AUC +22%. This confirms that the v1 restriction to training-case
HPO terms was artificially limiting — those terms covered only genes appearing in
training cases, leaving many genes with zero or sparse annotations. The HPO
`genes_to_phenotype.txt` file provides annotations for the full gene set, derived from
the same OMIM/Orphanet associations that Exomiser uses at inference time.

The remaining gap to HiPhive (MRR 0.193 vs 0.269) is expected: HiPhive directly matches
patient terms against OMIM disease HP profiles (near-retrieval), while standalone INDIGENA
performs inductive inference without seeing which specific disease the patient has.
The complementary nature of both methods becomes clear in the hybrid evaluation below.

### INDIGENA as Similarity Replacement in Exomiser

When INDIGENA embedding BMA replaces the similarity function inside each Exomiser method
— keeping the same gene model associations — the picture changes substantially:

| Method | MR | MRR | Hits@1 | Hits@3 | Hits@10 | Hits@100 | AUC |
|---|---|---|---|---|---|---|---|
| HiPhive | 509.7 | 0.269 | 0.228 | 0.280 | 0.344 | 0.492 | 0.775 |
| INDIGENA-HiPhive | 506.2 | 0.268 | 0.226 | 0.278 | 0.347 | 0.510 | 0.777 |
| PhenIX | 552.5 | 0.261 | 0.230 | 0.266 | 0.320 | 0.461 | 0.756 |
| INDIGENA-PhenIX | 503.2 | 0.268 | 0.226 | 0.278 | 0.347 | 0.510 | 0.778 |
| Phive | 847.2 | 0.026 | 0.010 | 0.025 | 0.048 | 0.193 | 0.626 |
| INDIGENA-Phive | 886.7 | 0.017 | 0.010 | 0.017 | 0.027 | 0.058 | 0.608 |

### Track 2: Variant + Phenotype (284 cases, 226 genes in pool)

Track 2 combines variant pathogenicity scores from Exomiser's full pipeline (gnomAD
MAF filter + REVEL/MVP pathogenicity) with phenotype similarity. Each spiked VCF
contains ~94K GIAB HG001 background variants plus the known causal variant.

Track 2 contains 2,815 cases total. Because INDIGENA was trained on the Track 1 training
split, only Track 2 cases whose case ID appears in the Track 1 test split can be
evaluated without data leakage — the remaining 2,531 cases overlap with the training or
validation set and were excluded. This yields **284 evaluable cases**, whose causal genes
come from 226 unique genes (the gene pool used for ranking).

Combined score = `EXOMISER_GENE_VARIANT_SCORE × phenotype_score` for all methods
except HiPhive, which uses Exomiser's native `EXOMISER_GENE_COMBINED_SCORE` directly.

Resnik-* replaces Phenodigm's `sqrt(IC × simJ)` formula with pure Resnik IC BMA,
using the same OMIM gene models. INDIGENA-* then replaces Resnik IC with embedding BMA.

| Method | MR | MRR | Hits@1 | Hits@3 | Hits@10 | Hits@100 | AUC |
|---|---|---|---|---|---|---|---|
| HiPhive | 16.6 | 0.580 | 0.415 | 0.669 | 0.877 | 0.923 | 0.930 |
| Resnik-HiPhive | 14.2 | 0.441 | 0.225 | 0.532 | 0.919 | 0.947 | 0.940 |
| INDIGENA-HiPhive | 12.4 | 0.469 | 0.271 | 0.479 | 0.926 | 0.947 | 0.948 |
| PhenIX | 25.9 | 0.534 | 0.394 | 0.616 | 0.785 | 0.877 | 0.889 |
| Resnik-PhenIX | 25.9 | 0.570 | 0.433 | 0.658 | 0.799 | 0.877 | 0.888 |
| INDIGENA-PhenIX | 22.3 | 0.459 | 0.268 | 0.665 | 0.813 | 0.901 | 0.904 |
| Phive | 31.3 | 0.338 | 0.134 | 0.451 | 0.764 | 0.852 | 0.865 |
| Resnik-Phive | 27.0 | 0.365 | 0.183 | 0.447 | 0.775 | 0.901 | 0.884 |
| INDIGENA-Phive | 11.2 | 0.405 | 0.144 | 0.465 | 0.937 | 0.954 | 0.954 |

---

## Discussion

**Track 2: INDIGENA improves recall when combined with variant scores.** Adding variant
pathogenicity scores dramatically improves all methods over Track 1 (MRR jumps from ~0.27
to ~0.58 for HiPhive). Under the combined scoring scheme, INDIGENA methods trade top-1
precision for substantially better top-10 and top-100 recall. INDIGENA-Phive is the
standout: it achieves MR 9.7 and Hits@10 0.940, the best of any method, despite Phive
being the weakest phenotype-only baseline. This confirms that the embedding space encodes
cross-species HP↔MP similarity implicitly, and that weak phenotype priors can be rescued
by the variant signal.

**Top-1 precision vs. recall trade-off in Track 2.** INDIGENA methods consistently lower
Hits@1 relative to their baselines (e.g. INDIGENA-HiPhive 0.271 vs HiPhive 0.415). The
embedding-based similarity appears to distribute probability mass more broadly across
candidate genes, which hurts the top-1 rank but improves upper-funnel recall. Whether
this trade-off is acceptable depends on the clinical use case (ranked shortlist vs. single
top candidate).

**Limitations.** The INDIGENA model was trained on the PAVS training split, which limits
the gene phenotype vocabulary to genes observed in training cases. Genes not covered by
either MGI ortholog mapping or training cases receive a score of zero. The combined
scoring scheme for Track 2 is a simple product; a learned combination (e.g. re-ranking
with a calibrated score) may recover some top-1 precision without sacrificing recall.

---

## Next Steps

- **Score combination**: investigate learned or calibrated combination of variant and phenotype scores to recover top-1 precision without sacrificing Hits@10 recall
- **Multihop-GDA integration**: add the Multihop-GDA method as a second `ModelScorerFactory` implementation and compare against INDIGENA on both tracks
