# Lab Notebook — Developer

Per-role lab notebook for Developer sessions. Started 2026-05-06 as part of the lab-notebook split (see `lab_notebook.md` freeze note for rationale).

Format and rules unchanged from the unified notebook — see `shared/feedback_lab_notebook.md`. New date sections at the TOP; new time sections at the TOP of the date block; entries are immutable once committed.

---

## 2026-05-06

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
