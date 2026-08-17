# EmbedPVP phenotype scorer — embedding export contract

This document fixes the file formats and the identity rules that the Java
`EmbedpvpModelScorer` depends on. It is written **before** the Java, because the Java
cannot recover any of it at runtime: the scorer receives phenotype terms, not a case
identifier, and the case embedding it needs does not exist until an offline fine-tune has
produced it.

## What the scorer computes

EmbedPVP's phenotype term is **not** a comparison of phenotype term sets. The
authoritative definition is `EmbedPVP/embedpvp/main.py::calculate_graph_embed_similarity`:

    MS(variant) = max over the variant's genes g of  (cosine(e_case, e_g) + 1) / 2

with the final EmbedPVP score `0.6 * MS + 0.4 * CADD_norm`, where `CADD_norm` is CADD
min-max rescaled inside that case's own variant table.

Two consequences that shape the contract:

1. **The similarity is between two *entity* embeddings**, one for the case and one for the
   gene — not between two sets of phenotype terms. INDIGENA's `IndigenaModelScorer`
   compares term sets (symmetric BMA of `sigmoid(dot)`); nothing in that class is reusable
   here beyond the TSV loader.
2. **The case embedding is not derivable from the case's HPO terms by the scorer.** It is
   the output of the inductive fine-tune, which runs offline in Python. It must therefore
   be exported to disk and read back, and the scorer must be told *which* case it is
   scoring.

## Files

Both files use the format `IndigenaEmbeddings` already reads — one entity per line,
tab, then comma-separated float components:

    <key>\t<v1>,<v2>,...,<vD>

### 1. Gene embeddings — `<run>_gene_embeddings.tsv`

Key: the mOWL entity IRI for the gene, `http://mowl.borg/<entrez_gene_id>`.

    http://mowl.borg/79621	-0.148693,0.076955,...

Entrez is the join key because that is what Exomiser carries (`GeneModel.entrezGeneId()`,
`ENTREZ_GENE_ID` in the TSV output) and what EmbedPVP itself keys genes by (it resolves
symbols through `Homo_sapiens.gene_info` to a GeneID before looking up
`http://ontology.com/<GeneID>`). Symbols are not used anywhere in the contract — they are
ambiguous and drift between releases.

The existing INDIGENA export (`eval/export_embeddings.py`) already emits gene rows in
exactly this form, so a plain INDIGENA embedding TSV is a valid gene-embedding file.

### 2. Case embeddings — `<run>_case_embeddings.tsv`

Key: the **case id exactly as it appears in the phenopacket and the split files**, e.g.
`PAVS:A0000003`. Not the mOWL IRI: the IRI mangles the colon (`PAVS_A0000003`) and there
is no reason to make the Java undo an encoding.

    PAVS:A0000003	0.031,-0.204,...

An optional third tab-separated field carries the semicolon-joined, sorted HPO term set
the export was computed from:

    PAVS:A0000003	0.031,-0.204,...	HP:0002378;HP:0003593;HP:0012378

This field is a **wiring check, not an input**. The scorer intersects it with the query
terms Exomiser hands it and logs a warning if the two sets are disjoint while both are
non-empty — the signature of a case embedding attached to the wrong patient. It is not
fatal, because Exomiser legitimately drops query terms that are absent from its own HPO
release, so an exact match is not expected.

Transductive models that already contain a case entity (the July
`indigena_mp_embeddings.tsv` has `http://mowl.borg/PAVS_A0000003` rows) can be converted
to this file mechanically; inductive models must export it from the fine-tuned encoder.

## Which case am I scoring?

The scorer is constructed per prioritisation call with query phenotype terms only, so the
case identity has to arrive out of band. Two supported routes, in precedence order:

1. **Explicit** — `new EmbedpvpModelScorerFactory(embeddings, "PAVS:A0000003")`. Used by
   the JPype harnesses, which loop over cases in one JVM.
2. **System property** — `-Dembedpvp.caseId=PAVS:A0000003`, read when the constructor is
   given a null id. Used by the CLI path, where `run_exomiser_variants.py` already forks
   one JVM per case.

Keying by a hash of the query phenotype set was rejected: Exomiser filters query terms
against its own HPO release, so the set the Java sees is not guaranteed to equal the set
the Python export hashed, and a silent hash miss would score every gene 0.

## Known limitation of the `ModelScorer` extension point

`ModelScorer` scores **models** (a gene's mouse orthologs, its diseases), and the
prioritiser then takes the max over a gene's models. A gene with **no** model is never
passed to the scorer at all — and `PhivePriority` does not give it 0, it substitutes the
constant `NO_MOUSE_MODEL_SCORE = 0.6f` (`PhivePriority.java:63`).

EmbedPVP's MS has no such notion: it is defined for every gene that has an embedding,
annotated or not. So `EmbedpvpModelScorer` reproduces EmbedPVP exactly for genes that reach
it, and the prioritiser overrides it for genes that do not.

Measured on 6,972 (case, gene) pairs from 22 test cases
(`embedpvp2/results/dualtrack/`, 2026-08-17):

* the scorer itself agrees with EmbedPVP's own function to **6.7e-8** (float32 vs float64);
* **1,980 / 6,972 pairs (28.4 %)** never reach the scorer and come out of `PhivePriority`
  as exactly 0.6;
* **119** of those had a genuine non-zero EmbedPVP MS that the 0.6 default discarded.

That 28 % is precisely the no-annotation slice design decision D4 says must not be dropped,
so the prioritiser path cannot be used to reproduce EmbedPVP. Closing the gap needs a
`Prioritiser` implementation that scores every gene directly, not a `ModelScorer`.

The variant-level harness in `embedpvp2/code/eval/combine_variant_scores.py` takes the
other route — an external per-gene score table joined onto Exomiser's per-variant output —
and has no such gap. The Java scorer exists so the method is available *inside* Exomiser;
the harness is what the benchmark numbers come from.
