# vexomiser

Fork of [Exomiser](https://github.com/exomiser/Exomiser) extended with vector/embedding-based gene–disease association (GDA) scoring methods.

## Goal

Integrate two GDA methods into Exomiser's variant prioritisation pipeline:

- **INDIGENA** — similarity-based GDA (inductive, embedding-based)
- **Multihop-GDA** — link-prediction / logical query answering for GDA

Both are developed by Fernando Zhapa-Camacho. This project is the delivery channel for both into a clinical variant prioritisation tool.

## Collaborators

- **Azza Althagafi** — generates data splits, weekly Monday meetings (10:00)
- **Fernando Zhapa-Camacho** — Exomiser fork, placeholders for INDIGENA and Multihop-GDA

## Data

Azza's data splits and README live on Ibex:

```
/ibex/scratch/projects/c2014/exomiser_extension/data
```

Read Azza's README there before implementing the scoring hooks — it describes the expected input/output format.

## Development plan

1. Fork Exomiser upstream into this repo.
2. Identify where variant/gene scoring is plugged in (the `GeneScorer` or equivalent extension point).
3. Add placeholder implementations for INDIGENA and Multihop-GDA scores.
4. Wire up to Azza's data splits for evaluation.

## Org-mode project entry

`~/org-mode/research_projects.org` → `* Exomiser extension`
