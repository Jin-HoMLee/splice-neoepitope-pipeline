# Slide decks — convention

Lab-seminar / external-talk slide decks for the splice-neoepitope pipeline. Built in [Quarto](https://quarto.org/) so figures, math, and citations render to **reveal.js HTML** (screen) + **PDF** (handout) from a single Markdown source.

## When to make a deck

**One deck per experiment Issue.** A deck ships in the same PR cycle as the experiment write-up (notebook + decision section). Don't ship retroactively.

- ✅ Experiment Issues (e.g. AlphaGenome chr22 PoC, GTEx filter integration, head-to-head filter comparison)
- ❌ Closure tasks, doc updates, single-fix PRs, sub-issues of a larger experiment

Audience target: **lab seminar / external talk** (~20-25 min, conference-quality figures, deeper methods than a PI weekly).

## Layout

```
research/slides/
├── README.md                                       # this file
├── _template.qmd                                   # reusable template
└── issue_NNN_<short-content-desc>/                 # one folder per deck
    ├── slides.qmd                                  # the deck
    ├── refs.bib                                    # bibliography (BibTeX)
    └── figures/                                    # PNG exports + regen script
        ├── _regenerate_figures.py                  # rebuilds PNGs from notebook outputs
        ├── *.png
```

**Folder naming:** `issue_NNN_<short-content-desc>` — ties to Issue # AND describes content. Example: `issue_393_alphagenome_chr22_poc/`.

**File naming:** `slides.qmd` (generic; doesn't repeat folder name).

## Why Quarto (not Marp / PowerPoint / Keynote)

| Need | Quarto | Marp | PPT/Keynote |
|---|:-:|:-:|:-:|
| Native LaTeX math (e.g. AP = Σₙ (Rₙ−Rₙ₋₁) Pₙ) | ✅ | ⚠️ | ⚠️ |
| BibTeX citations from Zotero | ✅ | ❌ | ❌ |
| reveal.js HTML + PDF from one source | ✅ | ✅ | ❌ |
| PR-reviewable text diffs | ✅ | ✅ | ❌ |
| Dominant in computational biology | ✅ | ⚠️ | — |

PowerPoint/Keynote are binary-blob-hell for git review. Marp is acceptable but weaker on math and citations.

## Install

**Recommended (macOS dev):**

```bash
brew install --cask quarto
```

**Alternative (conda, no sudo, fits project's conda-everywhere pattern):**

```bash
conda create -n splice-slides -c conda-forge quarto -y
conda activate splice-slides
```

Either path provides the `quarto` binary used by the render commands below. The brew install is system-wide; the conda env requires `conda activate splice-slides` before each render.

## Render

From the deck's folder:

```bash
cd research/slides/issue_NNN_<short-desc>/

# HTML (reveal.js) — primary delivery, for screen
quarto render slides.qmd --to revealjs

# Live preview while editing (auto-reloads on save)
quarto preview slides.qmd
```

Outputs land alongside `slides.qmd` as `slides.html` (+ `slides_files/` for assets). Both are gitignored (regenerable from source).

**PDF rendering (handout / archive):** the canonical Quarto path is `--to pdf` (uses beamer + LaTeX). The pilot deck disabled this — Quarto/pandoc's auto-emitted preamble hits a `\makesavenoteenv{longtable}` interaction with `footnotehyper.sty` that fails LaTeX compile. **Workaround for an immediate handout:** open `slides.html` in Chrome → *File → Print → Save as PDF*. A proper LaTeX fix is tracked as a follow-up; once resolved, re-enable the `pdf:` format block in the deck's YAML (commented out in `_template.qmd`).

## Figures

**The notebook stays canonical.** Slide figures regenerate from the experiment's `research/notebooks/<exp>_outputs/*.parquet` cache — slides aren't a second source of truth for results.

Each deck's `figures/_regenerate_figures.py` reads the cached parquet, recomputes the metrics, and writes individual PNGs (one per slide). Re-run after editing the notebook to keep slides in sync:

```bash
conda activate splice-neoepitope-alphagenome   # or whichever env the notebook uses
python research/slides/issue_NNN_<short-desc>/figures/_regenerate_figures.py
```

This pattern keeps notebook and slides reproducible from the same cached prediction parquet.

## Bibliography

**For the pilot:** hand-roll `refs.bib` (~5-10 entries per deck; small enough that automation isn't worth the wiring).

**Upgrade path (future):** `research/scripts/zotero_export_bib.py` (mirroring the existing `zotero_add.py` pattern) that pulls Zotero collection [Z38GTJNW](../scripts/zotero_add.py) → `refs.bib`. Optional further upgrade: [Better BibTeX](https://retorque.re/zotero-better-bibtex/) Zotero plugin auto-sync.

Cite in the `.qmd` with standard Pandoc syntax: `[@avsec2024alphagenome]` → renders inline with reference list slide at the end.

**Citation style:** Nature CSL is cached locally at [`nature.csl`](nature.csl) and referenced as `csl: ../nature.csl` from each deck's YAML. Committing it avoids the `zotero.org/styles/nature` round-trip on every render (offline-render friendly; no flaky-CDN dependency).

## Template usage

`_template.qmd` is a skeleton with placeholder slides for title, framing, methods, results, caveats, decision, references. Copy it into a new deck folder and edit:

```bash
mkdir -p research/slides/issue_NNN_<short-desc>/figures
cp research/slides/_template.qmd research/slides/issue_NNN_<short-desc>/slides.qmd
# edit slides.qmd in place
```

## Cross-references

- Convention established 2026-05-18 in [Issue #401](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/401).
- Pilot deck: [`issue_393_alphagenome_chr22_poc/`](issue_393_alphagenome_chr22_poc/) (covers [Issue #393](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/393) AlphaGenome chr22 PoC).
- Next-up decks: Exp 3 [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) (filter comparison); GTEx integration [Issue #212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212).
