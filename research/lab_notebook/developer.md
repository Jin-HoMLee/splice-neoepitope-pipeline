# Lab Notebook — Developer

Per-role lab notebook for Developer sessions. Started 2026-05-06 as part of the lab-notebook split (see `lab_notebook.md` freeze note for rationale).

Format and rules unchanged from the unified notebook — see `shared/feedback_lab_notebook.md`. New date sections at the TOP; new time sections at the TOP of the date block; entries are immutable once committed.

---

## 2026-05-08

### 15:46 UTC — Editor: Developer

**Headline:** [PR #301](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/301) (Issue #215 — full filtering audit trail) reviewed and ready to merge after a fix iteration: bot review surfaced 2 functional bugs (empty-pipeline early returns + hardwired `proteome` input in the aggregator) plus 3 minor cosmetic items, all addressed across 4 fix commits with 12 new regression tests; re-review came back "all clear, ready to merge"; CI green.

**Work shipped (morning + early afternoon):**

- [PR #301](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/301) (Issue #215 — full filtering audit trail) — bot review iteration:
  - **Bug 1 (empty-pipeline early returns):** `assemble_contigs.py` and `run_mhcflurry.py` early-return paths touched the canonical output (FASTA / empty TSV) but skipped the new stats TSV — would have cascaded into a missing-input failure on the aggregator any time a patient yielded zero contigs or zero peptides. Factored a `_write_zero_stats` helper into each script (`assemble_contigs.py` uses `pd.DataFrame.to_csv`; `run_mhcflurry.py` uses `csv.writer` for parity with the existing main-path stats writer). 6 new tests across both files. Commits 3fccba1, c3f2042.
  - **Bug 2 (proteome optional):** the aggregator hardwired `proteome=...` as a rule input, but `proteome_filter.smk` only defines its rule when `proteome_filter.enabled: true`. With proteome filter disabled (a supported config), the aggregator failed with a missing input. Fix in commit 6e2cbae: `analysis.smk` aggregator input is now a function (`unpack`) gating the `proteome` key on `_PROTEOME_FILTER_ENABLED_REPORT`; `aggregate_filtering_stats.py` makes `proteome_tsv` optional (default `None`); Snakemake entry uses `getattr(sm.input, "proteome", None)`.
  - **Minor 1 (vocab):** `analysis.smk` docstring "binder" → "presenter" (CLAUDE.md vocabulary; folded into the same commit as bug 2).
  - **Minor 2 (NaN render):** `generate_report.py` `_build_filtering_funnel_html` adds `.fillna(0)` after each `.reindex(...)` so partial-category inputs don't render "NaN" cells. The top-level `df.fillna("")` before the pivot was insufficient — reindex introduces new NaN columns post-pivot. Commit 04069b7.
  - **Minor 3 (column validation):** `_read_per_sample_stats` / `_read_per_patient_stats` now raise a clear `ValueError` listing the missing columns (`_PER_SAMPLE_REQUIRED` / `_PER_PATIENT_REQUIRED` set diff) instead of letting schema drift surface as an opaque `KeyError` at the unified-schema reindex downstream. Folded into commit 6e2cbae.
- 259 pytest tests pass (was 247; +12 new). CI green (pytest + snakemake-dry-run both SUCCESS).

**Issues created today:**

- [Issue #304](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/304) (sensitivity-analysis utility scoping ticket; P2/M; blocked on full read of [Prélot et al. bioRxiv 2025.09.10.674685](https://doi.org/10.1101/2025.09.10.674685) for the 35-parameter matrix). Surfaced from this morning's news briefing.
- [Issue #307](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/307) (`fix(scripts): zotero_add.py — bioRxiv preprints need 'publicationTitle', not 'repository'`; P2/S; bug). Caught when adding the Prélot DOI to Zotero — the script's CrossRef preprint branch writes `repository: "bioRxiv"` (semantically correct per Zotero's preprint itemType) while bioRxiv's own .ris export and most CSL citation styles want `publicationTitle: "bioRxiv"` on a `journalArticle` itemType. [Issue #229](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/229) AC 1 (literal `--dry-run` smoke test against bioRxiv) had been deferred — that gap was a yellow flag I missed.

**Memory updates (broadcast at 14:54 UTC):**

- **Extended Always-in-effect** in `shared/MEMORY.md`: **Created-by attribution** now applies to **comments you author** (`gh issue comment`, `gh pr comment`, follow-up replies), not just issue/PR bodies. Place at the bottom of comments. Standup posts unchanged (already carry attribution via `From: <Role>` header). Updated `shared/feedback_github_workflow.md` with placement guidance.
- **NEW Always-in-effect** in `shared/MEMORY.md` (promoted from Reference): **No bare hash-numbers in GitHub text** — `#N` auto-links, so `AC #1`, `step #3`, `finding #7` all autolink to unrelated artifacts. Drop the `#`. Caught in own [Issue #307](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/307) body and [Issue #229](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/229) follow-up comment — both fixed retroactively.
- New shared memory `developer/shared/feedback_zotero_biorxiv.md`: warn before firing `zotero_add.py` on `10.1101/*` DOIs — the CrossRef path is unreliable for bioRxiv. Retire when [Issue #307](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/307) closes.
- Extended `developer/feedback_morning_routine.md`: news briefing scope now requires a concrete Dev hook (config flag, benchmark, code change, infra/dep update) — not just "pipeline-relevant." Caught after StriMap (Sci-territory TCR-pMHC predictor) was wrongly included in this morning's Phase 1.

**Process notes:**

- Bot-review iteration ran cleanly: 4 fix commits → push → re-fire `@claude review` → "ready to merge" verdict in 4m 19s. The shape worked because the bot's findings were specific (file + line + suggested fix) and verifiable against the codebase.
- Snakemake CLI's `--config 'proteome_filter={"enabled": false}'` override didn't apply for nested keys in this version — tried in zsh and the dry-run still listed `proteome_filter_peptides` in the DAG. Disabled-proteome path is covered at the script level by 2 unit tests; the gap (no snakemake-level demonstration of the disabled path) is flagged in the [PR #301 iteration comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/301#issuecomment-4407641304) — could close with a CI matrix entry if it becomes load-bearing.
- Standup-protocol exercise: posted today's broadcast (14:54 UTC) after the user enacted both rule changes (Created-by extension + bare-hash promotion). PM posted a sibling broadcast at 15:17 UTC adding a 7-day archive cadence for `team_memory_broadcasts.md` — picked up via the `/cerebrum`-on-modified-shared-memory pattern.

---

## 2026-05-06

### 17:04 UTC — Editor: Developer

**Headline:** [PR #288](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/288) (snakevision) merged after `@claude review` (4/5 fixes applied); GPU-quota check on [Issue #285](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/285) confirmed migration off P100 is blocked on a Google quota request (decision: wait until tomorrow before filing); patient_002 cloud retry loop swept after 11 attempts / ~7h16m sustained outage; second small chore [PR #291](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/291) updates stale `TODO(#107)` references (libdeflate cap still in place, re-tested).

**Work shipped (afternoon):**

- [PR #288](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/288) (snakevision adoption) — `@claude review` returned 5 findings; 4 fixed in `c1d6528` (version pin `==1.1.0`, idempotent install guard, schema-reference comment for hardcoded TSV columns, runtime warning + doc caveat for `--clean` + `--config samples_tsv=…` mismatch); 1 skipped per reviewer's own "acceptable for this tool" note (YAML paths-with-spaces). Squash-merged at 16:26 UTC; closed [Issue #284](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/284) automatically.
- [PR #290](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/290) (morning lab-notebook entry) — merged at 16:23 UTC; first Developer entry in the new per-role `research/lab_notebook/` split.
- [PR #291](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/291) (chore: libdeflate-cap TODO refs from `#107` → `#237`) — comment-only update in `workflow/envs/hisat2.yaml` and `CLAUDE.md` "Known Dependency Issues". Annotated with the 2026-05-06 re-test result so future maintainers don't re-deliberate. Pending merge at write time.

**Investigation result — [Issue #237](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/237) (libdeflate-cap re-test):**

- Bioconda's 2026-03 samtools/htslib recipe refactor (samtools 1.23.1) did NOT lift the cap. Cross-platform solver test on `linux-64` (`CONDA_SUBDIR=linux-64 conda env create --dry-run`) with modern pins still reports `libdeflate >=1.20,<1.26.0a0` for htslib builds — conflicts with `regtools 1.0.0`'s `libdeflate >=1.26`. Without pins, solver falls back to ancient `samtools==1.6 / htslib==1.9` (2017-era).
- System-samtools workaround (apt 1.13) stays in place. Issue closes via PR #291 merge.

**[Issue #285](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/285) — P100 contingency:**

- Final tally: **11 launch attempts spanning 10:32 → 17:48 BST = ~7h16m of sustained `ZONE_RESOURCE_POOL_EXHAUSTED`**. First time the error has persisted past ~1 hour (CLAUDE.md previously noted only us-central1 transient instances).
- GPU quota check (option 2 from #285): we have **0 quota for T4 / A100 / L4** in either `europe-west1` or `europe-west4`; only K80, P100 (current), V100, P4 are accessible (`limit=1`). Global cap `GPUS_ALL_REGIONS=1`. Migration off P100 = Google quota request (days–weeks for approval), not a single-session task.
- Decision logged on [Issue #285](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/285): wait until tomorrow morning's status check before filing — single-day outage is weak justification; two-day evidence makes the request self-explanatory.
- EOD cron sweep at 17:57 BST: recurring retry loop `9a52fd5c` deleted; capacity still exhausted at sweep time.

**Process notes:**

- Force-pushed `feat/developer/issue-284-snakevision-adopt` after rebase on top of [PR #290](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/290) merge (mergeState had gone BEHIND). Rebase + `--force-with-lease` is the right move when the docs branch lands first and the feature branch then needs a sync — preserves linear history.
- Reaffirmed memory rule: lab-notebook entries are immutable per session; this afternoon's work warrants its own time-section under today's date, not an edit to the morning entry.

---

### 16:14 UTC — Editor: Developer

**Headline:** snakevision adopted via [PR #288](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/288) (closes [Issue #284](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/284)) for interactive Snakemake DAG visualization on rule changes; patient_002 cloud validation ([PR #278](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/278) merged yesterday) blocked all day by 6h+ P100 capacity exhaustion in `europe-west1-b` ([Issue #285](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/285) tracks contingency).

**Work shipped:**

- [PR #288](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/288) (snakevision adoption) — `setup_local.sh` adds `snakevision==1.1.0` to the `snakemake` env (idempotent install guard); new `scripts/visualize_dag.sh` wrapper runs `snakemake --rulegraph` → snakevision → `dag.svg` → opens it. `--clean` flag renders against a temp symlink workspace so all per-sample rules appear (works around snakemake's `--rulegraph` pruning of rules whose outputs already exist; auto-creates placeholder source FASTQs from the samples TSV). Tuned metro-style defaults (`scale=15 node_radius=10 edge_stroke_width=3`). Documented in `docs/installation.md` §6. 235 pytest + dry-run CI green. `@claude review` surfaced 5 findings; 4 fixed in `c1d6528`, 1 skipped per reviewer's own "acceptable for this tool" note.

**Issues created today:**

- [Issue #284](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/284) (snakevision adoption — closing via [PR #288](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/288))
- [Issue #285](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/285) (P100 europe-west1-b capacity contingency — multi-hour outage tracker; first time the error has persisted past 1 hour, now 6h+)

**Blocked:**

- patient_002 cloud validation: 7 launch attempts spanning 10:32–16:50 BST, all `ZONE_RESOURCE_POOL_EXHAUSTED` for `n1-highmem-8 + nvidia-tesla-p100`. WES-baseline results archived to `gs://splice-neoepitope-project/results/_archive/patient_002_wes_normal_2026-04-25/` before retries. 30-min cron retry loop running until 17:57 BST EOD. Will retry tomorrow morning if no luck by then.

**Memory updates:**

- New: `developer/feedback_brainstorming_scope.md` — for XS/S issues, skip the `superpowers:brainstorming` spec doc + `writing-plans` invocation; full flow only for M+ tasks. Saved after user pushed back on default `docs/superpowers/specs/` path.
- New: `developer/feedback_dag_visualization.md` — when changing Snakemake rules, run `bash scripts/visualize_dag.sh` after dry-run, before opening PR.
- Updated: `developer/feedback_morning_routine.md` — three-phase split (news → status → warm-up), each separated by a user pause (was previously two-phase: news → status+warm-up).

**Process notes:**

- snakemake's `--rulegraph` is rule-centric (no source files shown) and prunes rules whose outputs already exist regardless of `--forceall`. The `--clean` flag in `visualize_dag.sh` works around this by rendering against a fresh symlink workspace.
- Lab notebook split (per-role files under `research/lab_notebook/`) landed today via [PR #287](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/287) (Sci-led, closes [Issue #286](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/286)). This is the first Developer entry in the new file.

---
