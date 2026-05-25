# Issue #225 — Slide deck (co-located form)

**Establishing PR:** [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452)
**Issue:** [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) (sub-issue of parent [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) — Experiment 3)
**Date:** 2026-05-22

## Goal

Add a lab-seminar-quality Quarto slide deck for Issue #225 to PR #452, co-located with the existing notebook at [research/experiments/issue_225_normal_junction_filter_strength/slides.qmd](research/experiments/issue_225_normal_junction_filter_strength/slides.qmd). The deck supports human review of the experiment *now* (alongside the notebook in PR review) and serves as the lab-seminar artifact when [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203)'s decision rule is presented later.

This PR also establishes the slides-co-location convention in CLAUDE.md — the rule that experiment slide decks live next to the notebook they document, rather than under the parallel `research/slides/` tree.

## Scope (in this PR)

1. Deck source at `research/experiments/issue_225_normal_junction_filter_strength/slides.qmd`.
2. Figure regenerator at `research/experiments/issue_225_normal_junction_filter_strength/figures/_regenerate_figures.py` (matches prior-art naming).
3. Two regenerated figures: `figures/pr_curve.png`, `figures/caught_bar.png`.
4. One re-used figure: `outputs/filter_venn_chr22.png` (already in the notebook outputs — no copy, referenced by relative path from the deck).
5. Local hand-rolled `refs.bib` (~5-7 citations; re-uses keys from [research/slides/issue_393_alphagenome_chr22_poc/refs.bib](research/slides/issue_393_alphagenome_chr22_poc/refs.bib) where overlap; adds Snaptron + GTEx + MHCflurry where needed).
6. Rendered `slides.html` committed (see "Divergence from prior art" below for rationale).
7. CLAUDE.md update — document the slides-co-location rule (replaces the `research/slides/issue_NNN/slides.qmd` form in the existing "Slide decks for experiment Issues" section).
8. [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455) body update — remove the now-completed lines from its Migrations + Convention-docs sections (the [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) deck and the CLAUDE.md convention-docs update both ship here instead).
9. [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452) Test plan amendment — add deck-related items.

## Out of scope (deferred to [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455) or later)

- Migrating `issue_393_alphagenome_chr22_poc/` from `research/slides/` to the co-located form.
- Migrating `issue_224_*` / `issue_299_*` notebooks under `research/experiments/`.
- PDF / beamer render (disabled per CLAUDE.md footnotehyper conflict — same as prior-art deck).
- Zotero-export automation (`research/scripts/zotero_export_bib.py`) — hand-rolled `refs.bib` matches the established pilot pattern.

## Deck structure (12 slides incl. references, tight lab-seminar form)

| # | Slide | Content | Figure |
|---|---|---|---|
| 1 | **Title** (r-fit-text) | "Is AlphaGenome worth adding as a third normal-junction filter?" + #225/#203 anchor | none |
| 2 | **Why this matters** | Current pipeline = matched-normal filter; problem = MN not always available; three candidate replacements (MN / GTEx / AG). Today's question: does AG add value when stacked with GTEx? | none |
| 3 | **Method (block diagram)** | Mermaid: tumor + 3 filter sources → 3 caught sets → overlap analysis → decision rule | inline mermaid |
| 4 | **Data scope** | Patient_001, chr22, samples, tumor-set size (1,872), filter-set sizes (MN 1,714 / GTEx 880,769 / AG-predicted-normal @ τ = 448) | small table |
| 5 | **AG threshold via F1 sweep** | Universe-restricted F1 sweep against MN ∩ GENCODE (matches #224 §5). F1=0.300 at τ=3.16. Tight one-paragraph note on universe-restriction (positives=259 / universe=7,731) | `figures/pr_curve.png` |
| 6 | **Caught counts** | Per-source bar chart of tumor caught: MN 91 (4.9%) / GTEx 483 (25.8%) / AG 124 (6.6%) / union 503 (26.9%) | `figures/caught_bar.png` |
| 7 | **The killer slide: overlap** | Venn diagram + "AG ⊆ GTEx at the F1-max threshold" + "Only AG = 0" callout | `../outputs/filter_venn_chr22.png` (relative from `figures/` siblings; deck-level path = `outputs/filter_venn_chr22.png`) |
| 8 | **Decision rule** | Three threshold tests rendered as a table; F1<0.5 trips NO-GO clause directly; 0% AG-unique-vs-GTEx makes the verdict over-determined | decision table |
| 9 | **Decision** (r-fit-text) | 🔴 **NO-GO — treat as tissue prior** | none |
| 10 | **Caveats** | 5 items: Snaptron chr22 proxy vs production #211; Exp 2 deferred to #381 (WGS); chr22 only; pan-tissue not tissue-matched; AG threshold = F1-max from §2(b) | none |
| 11 | **Next steps** | Close #225, update [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) Exp 3 row, Exp 2 → [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) pending WGS | none |
| 12 | **References** | Auto-rendered from `refs.bib` | none |

12 slides total counting references. Within the "tight ~9 main + refs" envelope user picked (Decision/Caveats/Next/Refs are short closer slides).

## Figure regenerator (`figures/_regenerate_figures.py`)

**Pattern (mirrors prior art):**
- Single file under `figures/` with `_regenerate_figures.py` name.
- Reads from `REPO_ROOT / "research/experiments/issue_225_normal_junction_filter_strength/outputs/"` for `filter_overlap_table.tsv` and `chr22_gtex_panel.parquet`.
- Reads from notebook's `AG_PARQUET` path (cross-experiment dep, same fragility as the notebook — same #455 fixup tracks it).
- Reads `MATCHED_NORMAL_TSV` + `GENCODE_GTF` for ground-truth universe (same as notebook §2(b)).
- Emits two PNGs (`pr_curve.png`, `caught_bar.png`) next to itself.
- Idempotent: re-running produces byte-identical PNGs given same inputs.
- Run: `conda activate splice-neoepitope-alphagenome && python research/experiments/issue_225_normal_junction_filter_strength/figures/_regenerate_figures.py`

**What's intentionally NOT regenerated:** `outputs/filter_venn_chr22.png` — already produced by the notebook §5 cell, lives in `outputs/` (not `figures/`) because it's a notebook artifact reused by the deck, not a deck-only figure. Reading directly from `outputs/` keeps the notebook canonical per CLAUDE.md.

## Quarto setup

```yaml
---
title: "Is AlphaGenome worth adding as a third normal-junction filter?"
subtitle: "Comparative filter strength on patient_001 (chr22) · Issue #225"
author: "Jin-Ho Lee"
institute: "Splice Neoepitope Pipeline · JH M Lee Lab"
date: 2026-05-22
date-format: "YYYY-MM-DD"
format:
  revealjs:
    theme: simple
    incremental: false
    slide-number: c/t
    progress: true
    fig-align: center
    fig-cap-location: bottom
    transition: fade
    background-transition: fade
    width: 1280
    height: 720
bibliography: refs.bib
csl: ../../slides/nature.csl
execute:
  echo: false
  warning: false
---
```

**Paths (intentional shortenings vs `research/slides/` form):**
- `outputs/filter_venn_chr22.png` (1 level from deck) — was `../../experiments/issue_NNN/outputs/...` (4 levels) under the old layout
- `figures/pr_curve.png` (1 level)
- `bibliography: refs.bib` (local)
- `csl: ../../slides/nature.csl` (shared scaffolding stays centralized)

**Render command:** `cd research/experiments/issue_225_normal_junction_filter_strength/ && quarto render slides.qmd --to revealjs`

## Bibliography

Hand-rolled `refs.bib` — re-use keys from [research/slides/issue_393_alphagenome_chr22_poc/refs.bib](research/slides/issue_393_alphagenome_chr22_poc/refs.bib) where overlap, add new entries where needed:

| Key | Reuse / new | Cited from |
|---|---|---|
| `cotto2023regtools` | reuse | slide 4 (samples processing) |
| `kim2019hisat2` | reuse | slide 4 |
| `mudge2025gencode` | reuse | slide 5 (universe construction) |
| `pedregosa2011scikit` | reuse | slide 5 (sklearn PR curve) |
| `avsec2026alphagenome` | reuse | slide 6 (AG predictions) |
| `wilks2018snaptron` | **new** | slide 4 (GTEx panel construction) |
| `gtex2020v8` | **new** | slide 4 (GTEx context) |

**Implementation step for new entries:**
1. Pre-check Zotero collection Z38GTJNW for existing DOI matches (per `feedback_zotero_dedup_check`) — query `/items/top` for both DOIs.
2. For any missing entry: `python research/scripts/zotero_add.py <DOI> --tags slice-prediction gtex chr22 --note "..."` (tags space-separated per `feedback_zotero_dedup_check`).
3. Hand-copy the BibTeX entry into the deck's `refs.bib` (matches pilot pattern per [research/slides/README.md](research/slides/README.md)). DOIs: Snaptron `10.1093/bioinformatics/bty025`, GTEx v8 `10.1126/science.aaz1776`.

## CLAUDE.md update

**Existing block (CLAUDE.md "Slide decks for experiment Issues" section):**

> Decks live at `research/slides/issue_NNN_<short-content-desc>/slides.qmd`. See [`research/slides/README.md`](research/slides/README.md) for the full convention…

**Updated block:**

- Replace deck-location sentence: decks now live at `research/experiments/issue_NNN_<short-content-desc>/slides.qmd`, co-located with the notebook + outputs they document. Shared scaffolding (`_template.qmd`, `nature.csl`) stays centralized at `research/slides/`.
- Add one-sentence rationale: notebook + outputs + deck rename / archive / migrate as a unit; deck figure paths shorten from `../../experiments/issue_NNN/outputs/...` to `outputs/...`.
- Note that `research/slides/issue_393_alphagenome_chr22_poc/` predates the co-location rule and gets migrated under [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455).

## [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455) scope update

Strike two items (both descope into this PR):

- ~~CLAUDE.md "Experiment notebooks live under `research/experiments/`" section — slides-co-location rule documentation~~ → handled here
- ~~"A new deck for [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) itself" (mentioned under "Out of scope" as a future follow-up)~~ → handled here

[Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455) remains responsible for:
- `issue_224_*` notebook migration
- `issue_299_*` notebook migration
- `issue_393_alphagenome_chr22_poc/` slides migration to co-located form
- Cross-experiment path fix-ups (notebook + README in #225)
- ~~CLAUDE.md slides-co-location rule docs~~ ← removed

## Divergence from prior art (intentional)

1. **Committing `slides.html` and `slides_files/`.** [research/slides/README.md](research/slides/README.md) and the prior-art deck gitignore both. We commit them for this deck because the deck is added explicitly to support PR review — reviewers should be able to open it without a Quarto install. This divergence is local to [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225); a follow-up decision under [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455) can either harmonize (commit HTML for the migrated #393 deck too) or revert (gitignore here once the co-location convention is settled). Not load-bearing — pure ergonomics.

2. **CSL path shortened from `../nature.csl` to `../../slides/nature.csl`.** Direct consequence of the co-location form — the deck is one level deeper relative to `research/slides/`. The figure paths shorten by 3 levels as the trade-off (the intended win of co-location).

## PR #452 Test plan additions

- [ ] `slides.qmd` renders cleanly to `slides.html` via `quarto render slides.qmd --to revealjs`
- [ ] `figures/_regenerate_figures.py` produces `figures/pr_curve.png` + `figures/caught_bar.png` from notebook outputs
- [ ] All figure refs in `slides.qmd` resolve (Venn from `outputs/`, PR curve + bar from `figures/`)
- [ ] CLAUDE.md slides-co-location rule documented
- [ ] [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455) body updated to descope handled items
- [ ] Lab notebook entry for the deck-creation session

## Acceptance criteria

- `slides.html` opens in a browser and reads end-to-end as the experiment story (no broken refs, no half-finished slides).
- Each non-Venn figure regenerates byte-identically from `figures/_regenerate_figures.py`.
- Lab-seminar-quality polish (audience-appropriate text density; no `TODO` / placeholder slides).
- [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455) scope reflects the descoping.
- CLAUDE.md slide-deck section reads consistently with the new layout.

## Risks

- **Quarto install required to render.** Mitigated: `figures/_regenerate_figures.py` runs without Quarto (matplotlib only); `slides.html` is committed so reviewers don't need Quarto either.
- **Snaptron / GTEx citations not yet in Zotero.** Mitigated: known DOIs above; `zotero_add.py` call is a one-liner; pre-check duplicates per [`feedback_zotero_dedup_check`](.claude/memory/scientist/feedback_zotero_dedup_check.md).
- **Cross-experiment `AG_PARQUET` path coupling** (notebook's `research/notebooks/issue_224_*/...`). Mitigated: same fragility as the notebook itself; [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455) tracks the path fix-up under "Cross-experiment dependency fix-ups".
- **Snaptron HTTPS unavailability during render.** Not applicable — figure regenerator reads the cached `outputs/chr22_gtex_panel.parquet`, not a live API call.
