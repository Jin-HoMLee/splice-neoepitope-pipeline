# Slide decks — convention

Lab-seminar / external-talk slide decks for the splice-neoepitope pipeline. Built in [Quarto](https://quarto.org/) so figures, math, and citations render to **reveal.js HTML** (screen) + **PDF** (handout) from a single Markdown source.

## When to make a deck

**Three deck tiers — one deck per qualifying Issue, per tier.** A deck ships in the same PR cycle as the work it presents; don't ship retroactively.

- **Experiment** — per-Issue analysis with a notebook + `outputs/` (`research/experiments/`). E.g. AlphaGenome chr22 PoC, GTEx filter integration, head-to-head filter comparison.
- **Eval primer** — per-tool evaluation Issue (`research/evals/`). E.g. HERMES, ImmSET, Boltz-2.
- **Research decision** — methods/science call that gates downstream work (`research/decisions/`; tool-agnostic, no notebook/outputs). E.g. the Issue #592 cohort-calibration gate.

❌ Closure tasks, doc updates, single-fix PRs, and sub-issues of a larger experiment/eval/decision do **not** get a deck in any tier.

Audience target: **lab seminar / external talk** (~20-25 min, conference-quality figures, deeper methods than a PI weekly).

## Layout

```
research/slides/
├── README.md                                       # this file
├── _template.qmd                                   # reusable template
├── custom.scss                                     # shared project theme (palette, fonts, footers)
├── nature.csl                                      # citation style
└── issue_NNN_<short-content-desc>/                 # one folder per deck
    ├── slides.qmd                                  # the deck
    ├── refs.bib                                    # bibliography (BibTeX)
    └── figures/                                    # PNG exports + regen script
        ├── _regenerate_figures.py                  # rebuilds PNGs from notebook outputs
        ├── *.png
```

**Folder naming:** `issue_NNN_<short-content-desc>` — ties to Issue # AND describes content. Example: `issue_393_alphagenome_chr22_poc/`.

**File naming:** `slides.qmd` (generic; doesn't repeat folder name).

**Decks live in three sibling trees, not under `research/slides/`:**

- `research/experiments/issue_NNN_<short>/slides.qmd` — co-located with the experiment notebook + outputs. See [`CLAUDE.md` "Slide decks for experiment Issues"](../../CLAUDE.md).
- `research/evals/issue_NNN_<tool>/slides.qmd` — per-tool primer for a tool-evaluation Issue. See [`CLAUDE.md` "Slide decks for eval Issues"](../../CLAUDE.md). Example: [`research/evals/issue_218_hermes/`](../evals/issue_218_hermes/).
- `research/decisions/issue_NNN_<short>/slides.qmd` — research-decision deck for a methods/science call that gates downstream work (tool-agnostic option comparison; no notebook/outputs). See [`CLAUDE.md` "Slide decks for research-decision Issues"](../../CLAUDE.md). Example: [`research/decisions/issue_592/`](../decisions/issue_592/).

`research/slides/` itself holds **only the shared scaffolding** (`_template.qmd`, `custom.scss`, `nature.csl`, this README). Decks reference the scaffolding via `theme: [simple, ../../slides/custom.scss]` and `csl: ../../slides/nature.csl`. (`research/slides/issue_393_alphagenome_chr22_poc/` predates the co-location rule and migrates under [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455).)

## Shared theme (`custom.scss`)

All decks share a single SCSS theme at [`custom.scss`](custom.scss) — palette (deep blue primary `#1e5ba8` + muted gold accent `#b8860b`, derived from existing mermaid-diagram colors), typography (Inter / Source Sans 3 with system fallbacks; JetBrains Mono for code), heading underline accents, callout styling, table styling, footer + slide-number layout. Reference from a deck YAML via:

```yaml
format:
  revealjs:
    theme: [simple, ../../slides/custom.scss]    # for decks under research/evals/, research/experiments/, or research/decisions/
    # theme: [simple, custom.scss]               # for decks inside research/slides/ itself
    footer: "Splice Neoepitope Pipeline · JH M Lee Lab"
```

The theme overrides the `simple` reveal.js base — `simple` provides the structural CSS, `custom.scss` layers the project identity on top. **Don't fork the SCSS per deck** — edit `custom.scss` in place to evolve the shared look, and every deck picks up the change on next render.

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

### Diagrams (block diagrams, integration maps, flowcharts)

Three ways to draw a *non-data* diagram, in increasing effort / control:

- **ASCII / monospace** (a plain code fence, as in [`_template.qmd`](_template.qmd)) — quickest, zero-dependency, text-diffable. Crude (monospace box-drawing); no auto-layout, styling, or data. Fine for a rough placeholder or a sketch you'll replace.
- **mermaid** (fenced ```` ```{mermaid} ```` in the `.qmd`) — declarative, **auto-layout**. Best when the diagram's *structure* keeps changing (add / remove / reorder nodes, re-wire edges) and default styling is fine. Text-diffable, no regenerate step. Limited fine styling; can drift between HTML and PDF; cannot plot data.
- **matplotlib** (a `fig_*()` in `figures/_regenerate_figures.py` → PNG) — imperative, **hand-positioned**. Best for (a) anything with *plotted data*, (b) exact house-style polish / layered annotations on arrows, or (c) a stable "hero" diagram drawn once and frozen. Renders identically everywhere; matches `custom.scss` exactly. Costs layout labor on every structural change; PNGs are opaque in git diffs.

**Rule of thumb: ASCII for a throwaway placeholder, mermaid for evolving topology, matplotlib for visual polish or plotted data.** Match the shared palette either way (`#1e5ba8` blue, `#b8860b` gold). Worked matplotlib examples: the NeoGuider eval deck's combine schematic + integration map ([`research/evals/issue_258_neoguider/figures/_regenerate_figures.py`](../evals/issue_258_neoguider/figures/_regenerate_figures.py)).

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
