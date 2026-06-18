# Issue #393 — AlphaGenome chr22 PoC slide deck

**Issue:** [#393](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/393)
**Type:** slide deck only (no analysis notebook of its own)
**Source experiment:** [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) — the AlphaGenome predicted-normal validation PoC ([`../issue_224_alphagenome_exp1/`](../issue_224_alphagenome_exp1/))

## What this is

A lab-seminar Quarto deck condensing the #224 chr22 proof-of-concept: *does AlphaGenome predict tissue-expressed splicing from the GRCh38 reference alone?* The analysis itself lives in the #224 notebook — this folder holds only the presentation artifact (co-located per the slides convention; migrated here from `research/slides/` under [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455)).

## Contents

- `slides.qmd` — the deck (reveal.js; `csl: ../../slides/nature.csl`, `bibliography: refs.bib`)
- `figures/` — deck figures + `_regenerate_figures.py` (regenerates the P-R / threshold / bootstrap plots)
- `refs.bib` — deck-local bibliography
- `slides.html` + `slides_files/` — rendered output (kept on disk, **not** tracked; `quarto render slides.qmd`)

## Data dependency

`figures/_regenerate_figures.py` reads the #224 AlphaGenome predictions at
`../issue_224_alphagenome_exp1/outputs/chr22_stomach_predicted_junctions.parquet`
(cross-experiment dep, single producer in #224; that parquet is API-derived and not tracked in git — see the #224 README).

## Render

```bash
quarto render slides.qmd --to revealjs   # → slides.html (gitignored)
```
