# Feature documentation — convention

Per-feature slide decks (Quarto reveal.js) summarising **what a shipped pipeline feature does, how it performs on the chr22 test config, and where it fits**. Sibling convention to `research/slides/` — same tooling, different audience.

## When to make a deck

**One deck per feature-shipped Issue** with concrete biological output.

- ✅ Implementation Issues that ship a new Snakemake rule, a new conda env, or a substantive behavior change (e.g. [Issue #204](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/204) VDJdb panel rule)
- ❌ Closure tasks, doc-only PRs, lab-notebook PRs, refactor PRs without observable output change
- ❌ Analytical / experiment-tier work → use [`research/slides/`](../../research/slides/README.md) instead

**Distinction from `research/slides/`:** that folder is for analysis Issues (e.g. AlphaGenome chr22 PoC) that produce a decision-rule outcome. Decks here are for *feature ships* — "we built this rule, here's what it produces on test data, here's where it fits."

**Distinction from `research/experiments/`:** that folder holds per-Issue analytical notebooks (Jupyter + outputs/). Decks here are higher-level visual summaries with a different audience (team reference + portfolio peek), not raw analytical work.

## Layout

```
docs/features/
├── README.md                                       # this file
└── issue_NNN_<short-content-desc>/                 # one folder per feature deck
    ├── slides.qmd                                  # the deck
    ├── data/                                       # frozen chr22 outputs the deck reads from
    │   └── *.tsv                                   # e.g. panel.tsv, panel_qc.tsv
    └── figures/                                    # PNG exports + regen script
        ├── _regenerate_figures.py                  # rebuilds PNGs from data/
        └── *.png
```

**Folder naming:** `issue_NNN_<short-content-desc>` — ties to Issue # AND describes content. Example: `issue_204_vdjdb_panel/`.

**Why frozen `data/`:** the deck must be reproducible after the live pipeline outputs at `results/` are gone. The chr22 test outputs are tiny (~14 KB total for the VDJdb deck); committing them alongside is cheap insurance against the deck breaking when CI gets cleaned up.

## Render

```bash
quarto render docs/features/issue_NNN_<short>/slides.qmd
```

Outputs `slides.html` (gitignored) — opens in any browser. HTML is the primary delivery; for a PDF handout, "Print → Save as PDF" from Chrome. PDF render (beamer) is currently disabled pipeline-wide — see [`research/slides/README.md`](../../research/slides/README.md) "Render" for the LaTeX gotcha and follow-up.

## Regenerate figures

```bash
research/.venv/bin/python docs/features/issue_NNN_<short>/figures/_regenerate_figures.py
```

Reads from `data/*.tsv`, writes `figures/*.png`. The regen script is the canonical source of truth — re-run after editing the chr22 data files.
