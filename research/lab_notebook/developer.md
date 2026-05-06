# Lab Notebook — Developer

Per-role lab notebook for Developer sessions. Started 2026-05-06 as part of the lab-notebook split (see `lab_notebook.md` freeze note for rationale).

Format and rules unchanged from the unified notebook — see `shared/feedback_lab_notebook.md`. New date sections at the TOP; new time sections at the TOP of the date block; entries are immutable once committed.

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
