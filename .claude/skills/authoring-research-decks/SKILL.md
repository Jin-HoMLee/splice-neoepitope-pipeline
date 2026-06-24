---
name: authoring-research-decks
description: >-
  Conventions for authoring research artifacts in this repo - experiment
  notebooks (research/experiments/issue_NNN_<short>/ layout, cross-experiment
  data sharing, size bands) and the three Quarto slide-deck tiers (experiment /
  eval / research-decision). Use when creating or editing a slide deck, an
  experiment notebook, or their outputs/ folder.
---

## Experiment notebooks live under `research/experiments/`

Per-Issue experimental work (analysis notebooks + their cached outputs + a slide deck + a one-page README) lives at `research/experiments/issue_NNN_<short-content-desc>/`. The slide deck (`slides.qmd`) is co-located here rather than under `research/slides/` — see "Slide decks for experiment Issues" below.

Layout per experiment:

```
research/experiments/issue_NNN_<short>/
├── README.md          # one-page: goal, parent issue link, status, outputs index, cross-experiment deps
├── notebook.ipynb     # the analysis
└── outputs/           # cached artifacts (parquet, tsv, png)
```

**Distinguish from `research/notebooks/`:** that folder holds stable per-patient analyses (`patient_001_results.ipynb`, `patient_002_results.ipynb`) — long-lived, manuscript-supporting. The experiments/ folder holds scoped per-Issue work that lands once and then becomes a frozen reference. The convention is established by #225; existing per-Issue work (`issue_224_*`, `issue_299_*`) was migrated from `research/notebooks/` into `research/experiments/` under [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455).

### Cross-experiment data sharing

1. **Default:** each experiment owns its outputs in `<experiment>/outputs/`.
2. **Shared between ≥2 experiments → `research/experiments/_shared/`.** Promote an artifact here only when a 2nd consumer materializes (YAGNI; don't pre-share). Filenames carry provenance (`gtex_panel_chr22_snaptron_v1.parquet`, not `gtex_panel.parquet`).
3. **Cross-experiment read, single consumer → explicit path reference.** Document in the consumer's README under "Cross-experiment deps".
4. **Promoted to production → `resources/`.** When an artifact becomes a stable pipeline input.

**Bad practices (any size):** symlinks across experiments, copying artifacts, one experiment writing into another's `outputs/`.

### Size guidance

- **< 10 MB:** check into git.
- **10–100 MB:** check into git, commit a regenerator script alongside.
- **> 100 MB:** keep out of git; store in `gs://splice-neoepitope-project/experiments/<issue>/`; commit a `data_manifest.yaml` in the experiment folder listing artifact paths + checksums + fetch commands. Pin the manifest schema when the first artifact crosses 100 MB.

## Slide decks for experiment Issues

Every experiment-tier Issue ships a Quarto slide deck alongside the per-patient notebook + manuscript work. Decks live **co-located with the notebook** at `research/experiments/issue_NNN_<short-content-desc>/slides.qmd`. Shared scaffolding (`_template.qmd`, `nature.csl`) stays centralized at `research/slides/` and is referenced via `csl: ../../slides/nature.csl` from each deck. See [`research/slides/README.md`](research/slides/README.md) for the full convention (Quarto rationale, render commands, figure-source pattern, Zotero linkage, install). Scope: **one deck per experiment Issue**, not per sub-issue; no decks for closure tasks, doc updates, or single-fix PRs. Audience: lab seminar / external talk. Figures regenerate from `research/experiments/issue_NNN_<short>/outputs/*.parquet` (notebook outputs) and from a local `figures/_regenerate_figures.py` (deck-only figures) so the **notebook stays canonical**. Rationale for co-location: notebook + outputs + deck rename / archive / migrate as a unit; figure paths shorten from `../../experiments/issue_NNN/outputs/...` to `outputs/...`. (the `issue_393_alphagenome_chr22_poc` deck predated this rule and was migrated to `research/experiments/issue_393_alphagenome_chr22_poc/` under [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455).)

**Render tooling:** Quarto via `brew install --cask quarto` (macOS dev). PDF render (beamer) is currently disabled in decks — Quarto/pandoc's auto-emitted preamble hits a `\makesavenoteenv{longtable}` interaction with `footnotehyper.sty` that fails LaTeX compile. HTML reveal.js is the primary delivery; for a PDF handout, "Print → Save as PDF" from Chrome works on `slides.html`. A proper LaTeX fix is tracked as a follow-up.

## Slide decks for eval Issues

Every tool-evaluation Issue (e.g. [Issue #218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218) HERMES, [Issue #201](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/201) ImmSET, [Issue #188](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/188) Boltz-2) also ships a Quarto deck — a **tool primer** distinct from the consolidating publication deck. Decks live at `research/evals/issue_NNN_<tool>/slides.qmd`, sibling to `research/experiments/`. Shared scaffolding (`_template.qmd`, `nature.csl`) at `research/slides/` is referenced via `csl: ../../slides/nature.csl`.

**Why required even though evals are XS/S-sized:** visual condensation per tool is a different artifact from the consolidating decision deck (e.g. [Issue #432](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/432) TCR-pMHC scorer landscape). The per-tool deck answers *"what is this tool, and how does it work?"* — a comprehension precondition that lab notebook prose alone doesn't satisfy. The consolidating deck answers *"why this (a/b/c) decision across all 5 evals?"*. Both have value.

**Format (~8-10 slides; title + 7-9 content):** (i) the question (should we integrate?), (ii) tool primer (architecture, citation, code/weights status), (iii) why-it-works mechanism, (iv) integration map (block diagram showing where it'd plug in), (v) reasons-to-(a)/decline-(b)/skip-(c), (vi) decision + sub-issue ref, (vii) open scientific questions for the integration (if (a)), (viii) references. [`research/evals/issue_218_hermes/`](research/evals/issue_218_hermes/) is the reference example (9 slides total).

**Scope:** one deck per eval Issue. Established 2026-05-26 alongside the HERMES + ImmSET evals. Audience: lab seminar / external talk + manuscript-figure rehearsal for the consolidating publication Issue.

## Slide decks for research-decision Issues

Every research-decision Issue — a science/methods call that **gates downstream work** but is neither a per-tool evaluation nor a data-analysis experiment — ships a Quarto deck that walks the decision arc: question → options → evidence → decision → caveats. Decks live at `research/decisions/issue_NNN_<short>/slides.qmd`, a fourth sibling next to `research/slides/`, `research/evals/`, and `research/experiments/`. Shared scaffolding at `research/slides/` is referenced via relative paths: `csl: ../../slides/nature.csl` plus the shared reveal.js theme `theme: [simple, ../../slides/custom.scss]`. Like the eval/experiment decks, a decision deck carries self-contained inline YAML front-matter — the `_template.qmd` scaffold is copied and edited per deck, not formally extended. Bibliography is deck-local (`bibliography: refs.bib`). See [`research/slides/README.md`](research/slides/README.md) for the shared convention (Quarto rationale, render commands, Zotero linkage, install).

**Why required, and a distinct tier from eval / experiment decks:** the eval primer ([Issue #218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218) HERMES) answers *"what is this tool, how does it work, should we integrate it?"* — tool-centric, profiling one architecture and where it would plug in. The experiment deck presents *findings* from a per-Issue analysis notebook, with figures regenerated from that notebook's `outputs/*.parquet` so the notebook stays canonical. A decision deck is **tool-agnostic and question-centric**: it compares *options* (e.g. labelled cohorts) against scientific evidence to settle a methods choice, answering *"which option, and why?"*. It has no notebook and no `outputs/` directory — its evidence is a fact-checked literature sweep + primary-source checks, and its graphics are diagrammatic (inline `{mermaid}`) rather than data-derived plots. It presents a forward-looking *verdict* that gates future work, not findings from work already run.

**Scope: not every decision warrants a deck.** One deck per research-decision Issue, only when the decision genuinely gates other work and the options + evidence need visual condensation for a lab audience. Routine calls that an Issue-body decision + lab-notebook entry already carry do **not** get a deck; nor do closure tasks, doc updates, or single-fix PRs. Audience: lab seminar / external talk. Established 2026-06-01 ([Issue #597](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/597)), piloted by the [Issue #592](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/592) deck ([PR #596](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/596)). Convention is still under evaluation (n=1 pilot); revisit folding into the eval tier if a second decision deck does not materialize.

**Format (~8-10 slides; title + 7-9 content):** (i) the question (+ why-now gating callout), (ii) the gate — what blocks the downstream Issue (block diagram + the scientific sub-questions), (iii) the options/landscape (comparison table + pick callout), (iv) the key transfer/applicability question (for vs against), (v) corrections the fact-check caught, (vi) design choices the decision resolves (e.g. hold-out scheme), (vii) the decision verdict + forward links, (viii) caveats to read before citing, (ix) references (`::: {#refs}`). [`research/decisions/issue_592/`](research/decisions/issue_592/) is the reference example (the cohort-calibration gate for [Issue #547](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/547); 10 slides — title + 9 content). It predates this convention as the **pilot** that established the tier, and landed as `issue_592/` (no `_<short>` suffix); it keeps that bare name (out of scope for migration per [Issue #597](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/597)). New decks use `issue_NNN_<short>`.
