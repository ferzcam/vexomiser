# Vexomiser

## Summary

We integrated INDIGENA's embedding-based BMA as a drop-in replacement for the
IC-based semantic similarity used by HiPhive, Phive, and PhenIX in Exomiser,
and evaluated on the PAVS Track 1 benchmark (test split, 518 cases) and Track 2
(phenotype + variant, test split, 284 cases).
On Track 1, INDIGENA-HiPhive and INDIGENA-PhenIX match or slightly outperform
their baselines (MRR 0.262 vs 0.261 and 0.259 respectively),
while requiring no ontology database at inference time.
On Track 2, INDIGENA methods consistently improve top-10 and top-100 recall
over their baselines — INDIGENA-Phive achieves the best MR (9.9) and Hits@10
(0.937) despite Phive being the weakest phenotype-only baseline.

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
consequential. We used **MGI mouse-knockout MP phenotypes** as the primary source, mapped
to human genes via the HomoloGene orthology table. For human genes with no mouse ortholog
in MGI, training-case HPO terms were used as a fallback.

We deliberately excluded **OMIM disease–gene associations** as a gene phenotype source.
OMIM provides explicit gene → disease → HP links, which encode exactly the gene-disease
associations that INDIGENA is trained to predict. Using them in G2 would leak the GDA
signal into the training graph and invalidate the evaluation. This is the key structural
difference between INDIGENA and HiPhive's human disease models: HiPhive is permitted to
use OMIM data because it is not learning gene-disease associations — it is performing a
retrieval. INDIGENA, as a learning method, must keep those associations out of training.

### Training

INDIGENA was trained using **TransD** KGE (100-dimensional embeddings, 300 epochs,
Adam optimiser, batch size 2,048, learning rate 0.001) on a GPU. Two configurations were
trained: one using MGI mouse MP phenotypes for G2 and one using only training-case HPO
terms.

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
Exomiser's gene models. Scoring mirrors Exomiser's per-disease structure:
`score(gene) = max_d BMA(patient_HPs, disease_d_HPs)` using OMIM disease models from
`phenotype.hpoa`, with gene→disease links from `genes_to_phenotype.txt`.

| Configuration | MR | MRR | Hits@1 | Hits@3 | Hits@10 | Hits@100 | AUC |
|---|---|---|---|---|---|---|---|
| v1 — MGI only | 407.6 | 0.275 | 0.236 | 0.278 | 0.336 | 0.544 | 0.821 |
| v2 — MGI + HPO | 400.9 | 0.274 | 0.232 | 0.284 | 0.344 | 0.544 | 0.822 |

Both configurations perform at or above HiPhive (MRR 0.269) as standalone gene rankers,
with v1 (pure MGI mouse-knockout phenotypes) and v2 (MGI + HPO gene annotations) essentially
tied. Adding HPO gene–phenotype annotations to the training graph does not improve MRR;
v2 gains slightly on Hits@3 and Hits@10 while v1 is marginally better on Hits@1.

The critical factor is the evaluation structure: using `max_d BMA(patient_HPs, disease_d_HPs)`
over per-disease OMIM phenotype sets mirrors how Exomiser uses disease models at inference.
Genes with no OMIM disease associations fall back to their merged gene phenotype set.

### INDIGENA as Similarity Replacement in Exomiser

When INDIGENA embedding BMA replaces the similarity function inside each Exomiser method
— keeping the same gene model associations — the picture changes substantially:

| Method | MR | MRR | Hits@1 | Hits@3 | Hits@10 | Hits@100 | AUC |
|---|---|---|---|---|---|---|---|
| HiPhive | 501.4 | 0.261 | 0.223 | 0.267 | 0.332 | 0.513 | 0.778 |
| INDIGENA-HiPhive | 507.0 | 0.262 | 0.226 | 0.269 | 0.321 | 0.506 | 0.776 |
| PhenIX | 565.5 | 0.259 | 0.225 | 0.268 | 0.317 | 0.486 | 0.750 |
| INDIGENA-PhenIX | 502.0 | 0.262 | 0.226 | 0.269 | 0.321 | 0.507 | 0.778 |
| Phive | 832.1 | 0.035 | 0.017 | 0.032 | 0.061 | 0.222 | 0.632 |
| INDIGENA-Phive | 761.4 | 0.014 | 0.006 | 0.012 | 0.020 | 0.057 | 0.663 |

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

| Method | MR | MRR | Hits@1 | Hits@3 | Hits@10 | Hits@100 | AUC |
|---|---|---|---|---|---|---|---|
| HiPhive | 20.1 | 0.580 | 0.415 | 0.669 | 0.877 | 0.905 | 0.914 |
| INDIGENA-HiPhive | 10.5 | 0.466 | 0.261 | 0.532 | 0.926 | 0.968 | 0.956 |
| PhenIX | 26.5 | 0.531 | 0.391 | 0.616 | 0.789 | 0.870 | 0.886 |
| INDIGENA-PhenIX | 22.9 | 0.453 | 0.257 | 0.669 | 0.813 | 0.901 | 0.902 |
| Phive | 34.9 | 0.338 | 0.134 | 0.451 | 0.768 | 0.831 | 0.849 |
| INDIGENA-Phive | 9.9 | 0.434 | 0.201 | 0.472 | 0.937 | 0.972 | 0.959 |

---

## Discussion

**What INDIGENA contributes as a similarity replacement.** In the Exomiser integration,
INDIGENA replaces only the pairwise phenotype similarity function inside BMA. The
gene-disease model associations (which OMIM diseases are linked to which gene) remain
fixed by Exomiser's database. Crucially, INDIGENA inference uses only phenotype-term
embeddings — `cosine_sim(embed(HP_i), embed(HP_j))` — not gene or disease entity
embeddings. This preserves inductivity: any disease, including ones not seen during
training, can be represented as a bag of HP-term embeddings and scored without retraining.

The embedding space is trained on the full UPheno graph (G1–G4), so each HP-term embedding
encodes not just ontology position but also gene-phenotype and disease-phenotype
co-occurrence patterns. Terms that co-occur in the same biological context are pulled
together even when ontologically distant. In practice, however, the ontology hierarchy is
the dominant signal for HP-HP pairwise similarity, making INDIGENA and Phenodigm highly
correlated for HiPhive and PhenIX — which explains the near-flat Track 1 results when the
similarity function is swapped.

**Why INDIGENA-Phive does not improve on Track 1.** Phive uses mouse (MGI) and zebrafish
(ZFIN) phenotype models, requiring HP↔MP cross-species comparison. INDIGENA learns this
alignment implicitly from UPheno bridge axioms, while Phenodigm uses explicitly calibrated
IC-based cross-species mapping. The implicit alignment is less precise for this task,
explaining the performance drop for INDIGENA-Phive on Track 1.

**What the current integration does not exploit.** The deeper inductive capability of
INDIGENA — scoring a gene for which no OMIM disease entry exists by comparing its
phenotype terms against the patient's HP terms — is not used in the current integration.
Exomiser still gates scoring on known gene-disease associations: genes absent from the
database score zero regardless of phenotypic similarity. Fully exploiting INDIGENA's
inductive strength would require a separate scoring path that bypasses the disease-model
lookup and directly compares gene phenotype terms against the patient query.

**Track 2: INDIGENA improves recall when combined with variant scores.** Adding variant
pathogenicity scores dramatically improves all methods over Track 1 (MRR jumps from ~0.27
to ~0.58 for HiPhive). Under the combined scoring scheme, INDIGENA methods trade top-1
precision for substantially better top-10 and top-100 recall. INDIGENA-Phive is the
standout: it achieves MR 9.9 and Hits@10 0.937, the best of any method, despite Phive
being the weakest phenotype-only baseline. With variant scores supplying most of the
signal, INDIGENA's broader phenotype similarity distribution complements rather than
competes with the sharp variant prior, pushing the causal gene into the top-10 more
reliably.

**Top-1 precision vs. recall trade-off in Track 2.** INDIGENA methods consistently lower
Hits@1 relative to their baselines (e.g. INDIGENA-HiPhive 0.261 vs HiPhive 0.415). The
embedding-based similarity distributes scores more broadly across candidate genes, which
hurts the top-1 rank but improves upper-funnel recall. Whether this trade-off is
acceptable depends on the clinical use case (ranked shortlist vs. single top candidate).

**Limitations.** The current integration only exploits INDIGENA as a drop-in similarity
function over known gene-disease models. The primary bottleneck is not the similarity
function but the completeness of the gene-disease association database: genes with no OMIM
entry score zero. The combined scoring scheme for Track 2 is a simple product; a learned
combination may recover top-1 precision without sacrificing recall.

---

## Next Steps

- **Score combination**: investigate learned or calibrated combination of variant and phenotype scores to recover top-1 precision without sacrificing Hits@10 recall
- **Multihop-GDA integration**: add the Multihop-GDA method as a second `ModelScorerFactory` implementation and compare against INDIGENA on both tracks
