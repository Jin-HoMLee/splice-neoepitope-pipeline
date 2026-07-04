# Research artifact authoring conventions

Conventions for authoring research artifacts in this repo - experiment notebooks (`research/experiments/issue_NNN_<short>/` layout, cross-experiment data sharing, size bands) and the three Quarto slide-deck tiers (experiment / eval / research-decision).
Read this when creating or editing a slide deck, an experiment notebook, or their `outputs/` folder.

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

### Epic with a shared living artifact

The flat `issue_NNN_<short>/` layout above assumes **one experiment per Issue**. An **epic** whose sub-issues all mutate a single living artifact (a registry, a curated dataset) breaks that assumption: the directory is named for the *parent* epic, the shared artifact + its docs + its tooling are co-owned by every sub-issue, and each sub-issue that produces a *distinct analysis set* (its own notebook + outputs + deck) needs its own home so those sets don't pile flat at the top. The convention for this case: **shared core at the top level, one analysis subfolder per sub-issue**. Most sub-issues only edit the shared artifact and produce no subfolder; create a subfolder only for a sub-issue that yields its own analysis set. Established by [Issue #914](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/914) for the #680 splice-immunogenicity registry.

**Subfolder naming: reuse the parent pattern recursively — `issue_NNN_<short>/`.** A sub-issue is itself a numbered Issue, so the same `issue_NNN_<short>` rule applies; the parent/child relationship is carried by *nesting*, not by a separate `sub_`/`subissue_` prefix (which would re-encode what the path already says). Number-after-prefix also keeps subfolders sorting the same way the parent dir reads. One naming rule governs the whole `research/experiments/` tree, at every depth.

```
research/experiments/issue_NNN_<epic-short>/
├── registry.tsv  README.md  PROVENANCE.md  ...   # shared living core (every sub-issue mutates)
├── validate_*.py  derive_*.py  ...               # shared tooling
├── issue_734_db_audit/                           # one sub-issue's analysis set
│   ├── notebook.ipynb  outputs/  slides.qmd  refs.bib
├── issue_737_sparsity/                           # another sub-issue's analysis set
└── decoy_negatives/                              # shared data artifact (see carve-out)
```

**Carve-out — data artifacts get descriptive names; `issue_NNN_` is for analysis units.** A folder that *is* a single sub-issue's self-contained analysis (its own notebook/deck/outputs) takes `issue_NNN_<short>/`. A folder holding a **data artifact** — not an analysis — keeps a content-descriptive name, with the data file carrying any issue lineage in *its* name. Example: `decoy_negatives/presented_decoys_681.tsv` is a seed of presented-but-untested peptides whose eventual owner is sibling sub-issue #681 (immunopeptidome mining), but which is **produced, validated, and consumed entirely within this epic today** (materialized under #735, validated by `validate_registry.py`, read by the #737 notebook). Best practice is to **organize by *current* ownership, not aspirational ownership**: the artifact lives with the unit that produces+validates+consumes it now, and is *not* relocated into an `issue_681_*` dir that the not-yet-executed #681 has produced nothing in (YAGNI — don't scaffold a home before the work exists). **Migration trigger:** when #681 actually runs and becomes the authoritative producer of the decoy set, it stands up its own `research/experiments/issue_681_<short>/`, the seed migrates there, and this epic switches to a documented cross-experiment read (see "Cross-experiment data sharing" below). The filename's `_681` already records that future lineage.

A subfolder notebook reads the shared core one level up (`EXP = Path.cwd().parent` when `registry.tsv` isn't local) and writes to its own co-located `outputs/`; a subfolder deck's `csl:`/relative links gain one `../` level vs. a top-level deck.

### Cross-experiment data sharing

1. **Default:** each experiment owns its outputs in `<experiment>/outputs/`.
2. **Shared between ≥2 experiments → `research/experiments/_shared/`.** Promote an artifact here only when a 2nd consumer materializes (YAGNI; don't pre-share). Filenames carry provenance (`gtex_panel_chr22_snaptron_v1.parquet`, not `gtex_panel.parquet`).
3. **Cross-experiment read, single consumer → explicit path reference.** Document in the consumer's README under "Cross-experiment deps".
4. **Promoted to production → `references/` or `indices/`.** When an artifact becomes a stable pipeline input. A downloaded or derived *reference input* (genome / annotation / panel) goes to `references/`; a *pipeline-built alignment index* goes to `indices/`. (The old single `resources/` target was retired in [#63](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/63).)

**Bad practices (any size):** symlinks across experiments, copying artifacts, one experiment writing into another's `outputs/`.

### Size guidance

Two axes decide where an `outputs/` artifact lives: **size** and **regenerability**. Size sets the default; **regenerability is the tie-breaker** - the size band alone does not resolve an artifact that cannot be cheaply and deterministically rebuilt. "Regenerable" here means **offline-regenerable**: rebuildable from committed inputs with no paid API, external quota, or network dependency.

- **< 10 MB:** check into git.
- **10-100 MB:** check into git **only if** offline-regenerable from committed inputs, and commit the regenerator script alongside. If it is **not** offline-regenerable - e.g. derived from a paid or quota-limited external API with no committed regenerator - keep it out of git regardless of the band, and give it a **remote mirror + a `data_manifest.yaml` entry** (paths + checksums + fetch commands) so a fresh clone can fetch it rather than it living only on the producing machine.
- **> 100 MB:** keep out of git; store on the project's Cloudflare R2 bucket under `experiments/<issue>/` (the post-GCP-exit object store, [#854](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/854); credentials + S3-compatible client config are in the root `.env`); commit a `data_manifest.yaml` in the experiment folder listing artifact paths + checksums + fetch commands. Pin the manifest schema when the first artifact crosses 100 MB.

**The two existing experiments, reconciled against this rule:**
- `issue_225_normal_junction_filter_strength` - `chr22_gtex_panel.parquet` + tsv + png are Snaptron-derived and regenerable in-notebook, so **committed** (correct under the rule).
- `issue_224_alphagenome_exp1` - `chr22_stomach_predicted_junctions.parquet` (~17 MB) is AlphaGenome-API-derived with no offline regenerator (the notebook sweep needs live API access + quota), so **kept out of git** (correct). Under the regenerability tie-breaker it should additionally carry a **remote R2 mirror + `data_manifest.yaml` entry** so it is fetchable on a fresh clone; that mirror is not yet in place (follow-up), since today it exists only on the producing machine.

## Slide decks for experiment Issues

Every experiment-tier Issue ships a Quarto slide deck alongside the per-patient notebook + manuscript work. Decks live **co-located with the notebook** at `research/experiments/issue_NNN_<short-content-desc>/slides.qmd`. Shared scaffolding (`_template.qmd`, `nature.csl`) stays centralized at `research/slides/` and is referenced via `csl: ../../slides/nature.csl` from each deck. See [`research/slides/README.md`](../research/slides/README.md) for the full convention (Quarto rationale, render commands, figure-source pattern, Zotero linkage, install). Scope: **one deck per experiment Issue**, not per sub-issue; no decks for closure tasks, doc updates, or single-fix PRs. Audience: lab seminar / external talk. Figures regenerate from `research/experiments/issue_NNN_<short>/outputs/*.parquet` (notebook outputs) and from a local `figures/_regenerate_figures.py` (deck-only figures) so the **notebook stays canonical**. Rationale for co-location: notebook + outputs + deck rename / archive / migrate as a unit; figure paths shorten from `../../experiments/issue_NNN/outputs/...` to `outputs/...`. (the `issue_393_alphagenome_chr22_poc` deck predated this rule and was migrated to `research/experiments/issue_393_alphagenome_chr22_poc/` under [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455).)

**Render tooling:** Quarto via `brew install --cask quarto` (macOS dev). PDF render (beamer) is currently disabled in decks — Quarto/pandoc's auto-emitted preamble hits a `\makesavenoteenv{longtable}` interaction with `footnotehyper.sty` that fails LaTeX compile. HTML reveal.js is the primary delivery; for a PDF handout, "Print → Save as PDF" from Chrome works on `slides.html`. A proper LaTeX fix is tracked as a follow-up.

## Slide decks for eval Issues

Every tool-evaluation Issue (e.g. [Issue #218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218) HERMES, [Issue #201](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/201) ImmSET, [Issue #188](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/188) Boltz-2) also ships a Quarto deck — a **tool primer** distinct from the consolidating publication deck. Decks live at `research/evals/issue_NNN_<tool>/slides.qmd`, sibling to `research/experiments/`. Shared scaffolding (`_template.qmd`, `nature.csl`) at `research/slides/` is referenced via `csl: ../../slides/nature.csl`.

**Why required even though evals are XS/S-sized:** visual condensation per tool is a different artifact from the consolidating decision deck (e.g. [Issue #432](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/432) TCR-pMHC scorer landscape). The per-tool deck answers *"what is this tool, and how does it work?"* — a comprehension precondition that lab notebook prose alone doesn't satisfy. The consolidating deck answers *"why this (a/b/c) decision across all 5 evals?"*. Both have value.

**Format (~8-10 slides; title + 7-9 content):** (i) the question (should we integrate?), (ii) tool primer (architecture, citation, code/weights status), (iii) why-it-works mechanism, (iv) integration map (block diagram showing where it'd plug in), (v) reasons-to-(a)/decline-(b)/skip-(c), (vi) decision + sub-issue ref, (vii) open scientific questions for the integration (if (a)), (viii) references. [`research/evals/issue_218_hermes/`](../research/evals/issue_218_hermes/) is the reference example (9 slides total).

**Scope:** one deck per eval Issue. Established 2026-05-26 alongside the HERMES + ImmSET evals. Audience: lab seminar / external talk + manuscript-figure rehearsal for the consolidating publication Issue.

## Slide decks for research-decision Issues

Every research-decision Issue — a science/methods call that **gates downstream work** but is neither a per-tool evaluation nor a data-analysis experiment — ships a Quarto deck that walks the decision arc: question → options → evidence → decision → caveats. Decks live at `research/decisions/issue_NNN_<short>/slides.qmd`, a fourth sibling next to `research/slides/`, `research/evals/`, and `research/experiments/`. Shared scaffolding at `research/slides/` is referenced via relative paths: `csl: ../../slides/nature.csl` plus the shared reveal.js theme `theme: [simple, ../../slides/custom.scss]`. Like the eval/experiment decks, a decision deck carries self-contained inline YAML front-matter — the `_template.qmd` scaffold is copied and edited per deck, not formally extended. Bibliography is deck-local (`bibliography: refs.bib`). See [`research/slides/README.md`](../research/slides/README.md) for the shared convention (Quarto rationale, render commands, Zotero linkage, install).

**Why required, and a distinct tier from eval / experiment decks:** the eval primer ([Issue #218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218) HERMES) answers *"what is this tool, how does it work, should we integrate it?"* — tool-centric, profiling one architecture and where it would plug in. The experiment deck presents *findings* from a per-Issue analysis notebook, with figures regenerated from that notebook's `outputs/*.parquet` so the notebook stays canonical. A decision deck is **tool-agnostic and question-centric**: it compares *options* (e.g. labelled cohorts) against scientific evidence to settle a methods choice, answering *"which option, and why?"*. It has no notebook and no `outputs/` directory — its evidence is a fact-checked literature sweep + primary-source checks, and its graphics are diagrammatic (inline `{mermaid}`) rather than data-derived plots. It presents a forward-looking *verdict* that gates future work, not findings from work already run.

**Scope: not every decision warrants a deck.** One deck per research-decision Issue, only when the decision genuinely gates other work and the options + evidence need visual condensation for a lab audience. Routine calls that an Issue-body decision + lab-notebook entry already carry do **not** get a deck; nor do closure tasks, doc updates, or single-fix PRs. Audience: lab seminar / external talk. Established 2026-06-01 ([Issue #597](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/597)), piloted by the [Issue #592](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/592) deck ([PR #596](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/596)). Convention is still under evaluation (n=1 pilot); revisit folding into the eval tier if a second decision deck does not materialize.

**Format (~8-10 slides; title + 7-9 content):** (i) the question (+ why-now gating callout), (ii) the gate — what blocks the downstream Issue (block diagram + the scientific sub-questions), (iii) the options/landscape (comparison table + pick callout), (iv) the key transfer/applicability question (for vs against), (v) corrections the fact-check caught, (vi) design choices the decision resolves (e.g. hold-out scheme), (vii) the decision verdict + forward links, (viii) caveats to read before citing, (ix) references (`::: {#refs}`). [`research/decisions/issue_592/`](../research/decisions/issue_592/) is the reference example (the cohort-calibration gate for [Issue #547](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/547); 10 slides — title + 9 content). It predates this convention as the **pilot** that established the tier, and landed as `issue_592/` (no `_<short>` suffix); it keeps that bare name (out of scope for migration per [Issue #597](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/597)). New decks use `issue_NNN_<short>`.
