# Lab Notebook — Developer

Per-role lab notebook for Developer sessions. Started 2026-05-06 as part of the lab-notebook split (see `lab_notebook.md` freeze note for rationale).

Format and rules unchanged from the unified notebook — see `shared/feedback_lab_notebook.md`. New date sections at the TOP; new time sections at the TOP of the date block; entries are immutable once committed.

---

## 2026-05-13

### 10:25 UTC — Editor: Developer

**Headline:** [Issue #279](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/279) (HISAT2 `--rna-strandness F` for 10x R2 SE) shipped via [PR #358](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/358) — per-sample `strandness` column added to `samples.tsv`, pure helper module + 11 pytest cases, `alignment.smk` wires it through `params.strandness`. Closes [Issue #279](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/279) on merge.

**Work shipped:**

- **Schema:** new `strandness` column appended to `config/samples/*.tsv` headers (3 files). Values: `unstranded | forward | reverse` — biological direction, abstracted from tool-specific HISAT2 SE/PE syntax. Backward-compat: missing column / empty / unrecognized value → `unstranded` (no flag passed), preserves current pipeline behavior on `patient_001` / `patient_001_test`.
- **Pure helper:** `workflow/scripts/strandness.py` — `get_strandness_flag(strandness, is_paired_end)` maps `(biological direction, SE/PE)` → HISAT2 flag string (`F` / `R` / `FR` / `RF` or `""`). 11 pytest cases cover recognized values × SE/PE, missing/empty (backward compat), case-insensitive normalization, unrecognized strings (graceful fallback to `""`).
- **Rule wiring:** `workflow/rules/alignment.smk` — adds `_get_hisat2_strandness(wildcards)` wrapper that reads the row + delegates to the pure helper. `hisat2_align` gains `params.strandness`; shell conditionally injects `--rna-strandness {value}` only when non-empty (preserves current behavior on rows with no strandness entry).
- **Sample-side change:** `patient_002` PBMC normal (`PBMC_scRNA_Pool1_L002`) → `forward` (10x R2 = sense strand per [PR #350](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/350) slide 8). All other rows default to unstranded.

**Design decision (per-sample column over sample-type-driven):**

User picked per-sample column over the simpler sample-type-driven default. Tradeoff: per-sample column adds samples.tsv schema change + helper module (re-sized Issue XS → S during the prep body update) but is more robust — handles future non-10x R2 normals or stranded tumors without an `if sample_type == "normal"` heuristic that would break for bulk RNA-seq normals. Documented in [Issue #279](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/279) body + [PR #358](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/358) body.

**Process notes:**

- TDD: wrote `test_strandness.py` (11 cases) before `strandness.py`. Tests passed on first implementation — no debugging cycle.
- Skipped local `snakemake --dry-run` per the local-pipeline-run-ownership rule (user-only); rely on CI `pipeline-snakemake-dry-run` for rule-parse validation.
- Skipped DAG visualization — the change adds a `params` entry to an existing rule with no new edges, so DAG topology is unchanged.
- Issue body locked in the design decision (per-sample column) + re-sized XS → S during prep, documenting the scope change pre-code per `feedback_scope_discipline.md`.

---

### 08:49 UTC — Editor: Developer

**Headline:** [PR #350](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/350) review fixes + merge ([Issue #348](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/348) closes). [Issue #352](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/352) spike filed — direct evidence the "Last supporting combo is PyTorch 2.7 + CUDA ≤12.6" claim in [CLAUDE.md](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/blob/main/CLAUDE.md) is provably loose; PyTorch 2.12 cu126 wheels may keep Pascal SM 6.0 dispatch alive on our P100s.

**Work shipped:**

- **[PR #350](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/350) review-fix commit `8565450`** — addressed both items from the bot review:
  - §4 (real fix): Salmon column in the cheat-sheet was using PE codes only (`IU`/`ISF`/`ISR`); a SE reader copying `ISF` into a `salmon quant` invocation would error or silently mis-quantify. Split into Salmon SE + Salmon PE columns (`U`/`SF`/`SR` and `IU`/`ISF`/`ISR`) and added a warning blockquote about the SE-leading-`I` drift. Added Salmon `SF` to the pipeline-case summary line.
  - §3 (clarity): R2 label in `06-10x-chemistry.svg` was potentially confusing because the arrow points R→L (sequencing direction) into a "SENSE strand" block while alignment on the genome is L→R. Split the label into a bold conclusion line + italic clarifier noting the arrow shows sequencing direction (R2 primer reads back into cDNA), not alignment direction. ViewBox bumped 280→300 to fit the second line.
- **[Issue #352](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/352) filed** — spike to verify whether PyTorch 2.12 + CUDA 12.6 wheels dispatch on P100 (SM 6.0) at runtime. Surfaced during Phase 1 morning news after digging into the [CUDA 13.2 release thread](https://dev-discuss.pytorch.org/t/introducing-cuda-13-2-and-deprecating-cuda-12-8-release-2-12/3337):
  - **cu126 wheels exist for torch 2.8 → 2.12** (verified at [download.pytorch.org/whl/cu126/torch/](https://download.pytorch.org/whl/cu126/torch/)).
  - **[PyTorch RFC #178665](https://github.com/pytorch/pytorch/issues/178665) build matrix for cu126 in 2.12 explicitly includes Pascal (6.0):** `Maxwell(5.0), Pascal(6.0), Volta(7.0), Turing(7.5), Ampere(8.0, 8.6), Hopper(9.0)`.
  - Original [Pascal-removal dev-discuss thread](https://dev-discuss.pytorch.org/t/cuda-toolkit-version-and-architecture-support-update-maxwell-and-pascal-architecture-support-removed-in-cuda-12-8-and-12-9-builds/3128) explicitly says "removing Maxwell and Pascal GPU support from **CUDA-12.8 binaries**" — only cu128/cu129, not cu126. The cu126 channel has stayed Pascal-compatible through 2.12.
  - → CLAUDE.md's "Last supporting combo is PyTorch 2.7 + CUDA ≤12.6" line conflates the cu128 cut with all of 2.8+. Spike will confirm whether build-matrix listing translates to runtime kernel dispatch on actual P100 hardware — needed before any `python.yaml` refactor to lift the `torch<2.5` pin.

**Process notes:**

- **Foot-gun caught + recovered:** First version of the updated [PR #350](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/350) body included two literal `@claude` mentions ("validated via @claude bot review"), which would have re-triggered the Claude Code GitHub Action via the issue-comment event subscription. Per the no-`@claude`-in-bodies rule, rephrased to "the bot review" before pushing the lab notebook + merging. Reminder for future PR body edits: scan for literal `@claude` strings before submit.
- Bot review came back fast (5 min from ping yesterday at 21:38 UTC) with verdict "Approve with optional improvements." Solid signal-to-noise — caught a real user-facing bug (Salmon SE codes) that I missed during initial review.
- For [Issue #352](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/352): scope is deliberately XS (single SSH session + single `pip install` + smoke-test on a running P100 VM). Filed as a deferred spike rather than starting work today; user will run the smoke-test when convenient.

---

## 2026-05-12

### 21:19 UTC — Editor: Developer

**Headline:** New `docs/slides/` directory + 9-slide Marp visual primer on RNA-seq strandedness ([Issue #348](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/348), [PR #350](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/350)). Research surfaced two technical errors in the [Issue #279](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/279) body — follow-up correction comment to be posted there.

**Work shipped:**

- [Issue #348](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/348) filed (deck request, P3, XS); [PR #350](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/350) opened with the deck.
- `docs/slides/2026-05-12-rna-strandedness-primer.md` — [Marp](https://marp.app) deck, 9 slides, all SVG schematics inline. Covers: dsDNA setup, mRNA orientation, library prep flow (where strand info survives or is lost), PE FR/RF, SE F/R syntax pitfall, 10x 3' GEX R2 chemistry, HISAT2 `--rna-strandness` + XS tag mechanics, regtools `-s XS` dependency, and a tool-mapping cheat sheet (HISAT2 / htseq / Salmon / featureCounts).
- Sets a new convention for the repo: `docs/slides/YYYY-MM-DD-<topic>.md` for visual primers. First entry in the dir.

**Issue #279 corrections surfaced during research:**

1. **SE syntax:** HISAT2 SE takes a single letter (`F` or `R`); `FR`/`RF` is PE-only. Issue body uses `RF` for SE — copy-paste from PE docs.
2. **Direction:** 10x Chromium 3' GEX v3 R2 reads come from the **sense (coding)** strand ([10x Tech Note CG000376](https://assets.ctfassets.net/an68im79xiti/awNZTarmwqmwxKcvDv9wv/b05500661e36290b8ce59689ff889ea8/CG000376_TechNote_Antisense_Intronic_Reads_SingleCellGeneExpression_RevA.pdf), [scg_lib_structs 10xChromium3v3](https://scg-lib-structs.readthedocs.io/en/latest/ge/10xChromium3v3.html)) → forward-stranded. Issue body claims reverse-stranded.

→ Correct flag for 10x R2 SE alignment is `--rna-strandness F` (not `RF`). Implementation work on [Issue #279](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/279) should use `F`.

**Process notes:**

- Skipped spec doc + writing-plans flow per `feedback_brainstorming_scope.md` (XS docs work). Single `AskUserQuestion` covered scope/audience; second covered how to handle the issue-body discrepancy.
- Caught a wrong date in the initial filename (`2026-05-13` instead of today's UTC date `2026-05-12`) right after the first push — fixed via `git mv` + new commit before review. Reminder for future date-prefixed files: always confirm via `date -u` before naming.
- `marp-cli` not installed locally; trusted markdown syntax + manual SVG review. PR test plan asks reviewer to render in VS Code Marp extension as the validation step.

---

### 10:07 UTC — Editor: Developer

**Headline:** Cohort aggregation tool ships as a standalone research-time CLI ([Issue #84](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/84)) — `research/scripts/aggregate_cohort.py`, not a Snakemake rule. Also: [Issue #200](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/200) Snakemake 9 re-scope comment landed.

**Work shipped:**

- `research/scripts/aggregate_cohort.py` ([Issue #84](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/84) — cohort aggregation) — argparse CLI that concatenates per-patient `report.tsv` files (`patient_id | stage | metric | value | notes` schema) into a single cohort table. Two invocation modes: `--inputs PATH [PATH ...]` for explicit paths, or `--patients ID [ID ...] --results-root DIR` that resolves to `{root}/{ID}/reports/report.tsv`. Output sorted by `(patient_id, stage, metric)` for stable diffs. 7 pytest cases covering happy-path concat, sort order, schema validation, missing-file, and CLI entrypoint.
- [Issue #200](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/200) re-scope [comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/200#issuecomment-4429153617) — Snakemake 8→9 has exactly one breaking change (custom logger plugin API); much lighter than 7→8. Re-scopes the upgrade evaluation lower.
- Morning news_log entry shipped via [PR #335](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/335) (Snakemake 9.20 / PyTorch 2.12 / libdeflate; news_log PR is exempt from lab notebook by convention).

**Design choice on Issue #84:** the original sketch in the issue body assumed `PATIENT_IDS` would be globally available to a Snakemake `aggregate_cohort` rule. But `workflow/rules/common.smk::_read_samples_tsv` enforces single-patient-per-run (`raise ValueError(... must contain exactly one patient_id per run ...)`) — cohort aggregation is inherently a **cross-run** operation. User picked the cleanest reframe: ship as a standalone CLI in `research/scripts/` alongside `zotero_add.py`, not as a Snakemake rule. Pairs naturally with `research/notebooks/results_comparison.ipynb` (already named in `research/README.md`).

**Process notes:**

- Pre-empted closure-audit gap by backfilling a `**Priority rationale:**` line on [Issue #84](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/84) (pre-conventions issue from 2026-04-21). No AC checkboxes to tick — the issue body is prose-style, so `check_ac` finds no unticked boxes and stays clean.
- Skipped the spec-doc + writing-plans flow per `feedback_brainstorming_scope.md` (XS/S Issues only need it for M+ tasks). Single design check via `AskUserQuestion` covered the cross-run reframe; user redirected from any of the three pre-baked options to "standalone CLI in research/, not workflow/", which was the right call.

**Closure-audit smoke test (PR #335 + Issue #200 comment earlier today):** workflow ran on the PR-merge but stayed silent (clean) — no marker comment posted. First real-traffic data point for the [Issue #325](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/325) 1-week trial.

---

## 2026-05-11

### 16:32 UTC — Editor: Developer

**Headline:** [PR #332](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/332) opened against [Issue #325](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/325) — closure-audit critic ships as a 2-week probe of Sakana Fugu's "critic-as-default" lesson on a human-supervised pipeline.

**Work shipped:**

- [PR #332](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/332) ([Issue #325](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/325) — closure-audit critic) — GitHub Actions workflow fires on `pull_request.closed && merged` and `issues.closed`, runs `tools/ci/closure_audit.py` (~180 LOC) which performs the 3 closure-ritual checks (AC checkbox count, lab notebook date+`#<N>`-ref, `**Priority rationale:**` substring) and posts one marker-tagged comment listing gaps. Silent on the clean path; post-only (no edit-in-place). Path exemptions: `research/news_log.md`, `research/glossary.md`, `research/lab_notebook/*.md` (the last one to avoid a circular self-check). 9 focused unit tests cover the gnarly parsing (block-slicing across multiple `## YYYY-MM-DD` headers, deferral-comment handling, exemption filter, role resolution from `role:*` labels). Tests run in a separate `ci-tools-pytest` job in [`tests.yml`](.github/workflows/tests.yml) so a tooling test failure doesn't block real pipeline PRs.

**Process notes:**

- Two course-corrections during drafting: (1) the original implementation plan was 10 tasks + 26 tests, which the user correctly flagged as overkill for a ~180-LOC tool — pivoted to inline execution with 9 focused tests; (2) initial file placement under `workflow/scripts/` mixed CI tooling with Snakemake-pipeline files — moved to new top-level `tools/ci/` folder with a `conftest.py` for sys.path setup.
- Spec was preserved at [docs/superpowers/specs/2026-05-11-post-merge-critic-design.md](docs/superpowers/specs/2026-05-11-post-merge-critic-design.md) — design contract worth keeping even though the bloated plan doc was dropped.
- **Kill criterion baked into the PR body:** explicit re-evaluation at 2026-05-25. If the bot's comments are >50% noise, get ignored, or cost more maintenance than they save, rip it out. Issue #325 stays open during the trial.
- `gh auth refresh -s workflow` needed before the push could land — adding new `.github/workflows/*.yml` files requires the `workflow` OAuth scope that the default gh CLI auth lacks. One-time interactive step.

**Why this is a probe, not infrastructure:** the user already controls closure ritual manually as supervisor of 3 AI roles; PM's morning audit is the existing safety net; this critic is a safety-net-for-the-safety-net. The Sakana Fugu "critic-as-default" lesson is really for agent loops with no human in the loop. Worth knowing whether it adds value here, but no commitment to keep it.

---

## 2026-05-10

### 17:04 UTC — Editor: Developer

**Headline:** Two warm-ups shipped — [PR #319](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/319) ([Issue #309](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/309) — CLAUDE.md "Expected unavailability" subsection for the P100 capacity outages) and [PR #320](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/320) ([Issue #307](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/307) — `zotero_add.py` preprint path now emits `journalArticle`+`publicationTitle` to match bioRxiv's own .ris export, fixing citation-rendering breakage). [Issue #321](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/321) opened as a P3/XS follow-up to fix the local `pytest` collection hang that forced a direct-python test runner during PR #320.

**Work shipped:**

- [PR #319](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/319) ([Issue #309](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/309) — P100 unavailability docs) — single-file edit to `CLAUDE.md`: new `### Expected unavailability — P100 in europe-west1-b` subsection under `## Infrastructure` consolidating the 2026-05-06 → 2026-05-08 outage evidence (~46h sustained `ZONE_RESOURCE_POOL_EXHAUSTED`, 11 launch attempts on 05-06 over ~7h16m), mitigation (overnight retry), forward refs to [Issue #285](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/285) (parent epic) and [Issue #310](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/310) (T4/L4 hybrid fallback, blocked on Google T4 quota grant). Last-verified date `2026-05-08`. CI green; squash-merged at 16:35 UTC.

- [PR #320](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/320) ([Issue #307](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/307) — zotero bioRxiv fix) — `crossref_to_zotero` preprint branch flipped from `itemType=preprint`+`repository` to `itemType=journalArticle`+`publicationTitle`, matching the `TY=JOUR`+`JF=bioRxiv` ground truth from bioRxiv's own .ris export (user pasted a real bioRxiv .ris for Lu et al `10.64898/2025.11.30.691400` to anchor the choice). Dropped legacy `archiveID` (self-referential DOI duplicate). `main()` status-line discriminator switched from `item["itemType"] == "preprint"` to `_is_preprint(data)` since `itemType` no longer separates the two paths. 8 unit tests pass (added `test_biorxiv_preprint_matches_native_ris_export` anchored to .ris ground truth + `test_medrxiv_preprint_uses_medrxiv_as_publication_title` post-review for the medRxiv-on-same-codepath case). Bot review approved with 2 nits — both folded in via `7582fb2` (tighter inline comment + medRxiv test). Squash-merged at 16:54 UTC.

**Smoke-test approach for [Issue #307](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/307) AC:**

The "Real bioRxiv smoke test" AC was satisfied by **dry-run** rather than a real Zotero post. Rationale: dry-run output on Lu et al matched the bioRxiv .ris field-for-field — `itemType=journalArticle`, `publicationTitle=bioRxiv`, 5 authors exact, full abstract present, no legacy `repository`/`archiveID`. CrossRef's date (`2025-12-2`) was *finer* than the RIS's coarse `Y1=2025/01/01` default — i.e. the auto-add path is now strictly better than the manual workaround on at least one field. A real post on this DOI would have duplicated the existing manually-imported entry (the paper is tracked in [Issue #316](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/316)); deferred strict "real post" verification to the next bioRxiv DOI added in normal workflow — first one out of the gate post-merge implicitly serves as the live smoke test.

**Memory cleanup (post-merge):**

- Deleted `cerebrum/.../developer/shared/feedback_zotero_biorxiv.md` (warning was specific to the now-fixed CrossRef path).
- Removed the pointer line from `developer/shared/MEMORY.md` (was line 88).

**Issues opened today:**

- [Issue #321](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/321) (`pytest.ini` to exclude `.snakemake/` + `.venv/`) — caught when local `.venv/bin/python -m pytest` hung indefinitely; even `pytest --version` didn't return. Direct `python -c "from test_x import ..."` invocation runs in <1s, so the script + tests are fine; pytest's collection phase appears to be walking conda envs (Snakemake-built locally) for vendor-package conftest.py files. Workaround used in this PR: bypass pytest, run test functions directly. Permanent fix: a `pytest.ini` with `norecursedirs = .snakemake .venv ...` + `testpaths = workflow/tests research/scripts`. P3/XS, role:developer.

**Standup follow-ups (Pending → Done):**

- PM 2026-05-09 11:46 UTC ask (priority rationales for [Issue #310](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/310), [Issue #309](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/309), [Issue #307](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/307), [Issue #304](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/304)) — all 4 already had rationale lines at issue-creation time; PM's audit grep had missed them. Updated [Issue #310](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/310) (P1→P2, body now matches PM's revised triage) and [Issue #307](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/307) (P2→P1) bodies; left [Issue #309](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/309) and [Issue #304](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/304) as-is (already matched). Posted follow-up flagging the audit-format ergonomics gap; PM acked + flipped to Done.
- PM 2026-05-09 11:09 UTC ask ([Issue #215](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/215) lab-notebook closure-audit) — false positive; the 2026-05-08 15:46 UTC entry was bundled into [PR #301](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/301) itself (`git show e382d21 -- research/lab_notebook/developer.md` confirms +36 lines). Posted follow-up pointing at the existing entry; PM acked + flipped to Done.

**Memory updates:**

- **NEW Always-in-effect** in `developer/MEMORY.md`: *"News_log/standup-archive/memory-broadcast PRs are exempt from lab notebook entries."* The doc itself IS the journal record — duplicating it adds no signal. Self-promoted from the implicit `reference_news_log.md` hint ("log = journal, not deliverable") after user flagged it. Also clarified in `shared/reference_news_log.md` directly with an explicit "no lab notebook entry required" paragraph.
- **NEW Always-in-effect in `shared/MEMORY.md`** (added by user mid-session): *"Role-path only — never canonical shared."* Reach shared memory via `developer/shared/<file>` (the `<role>/shared/` symlink resolves to `../shared/`). Caught after my earlier `reference_news_log.md` edit used the canonical `splice-neoepitope-pipeline/shared/...` path — going forward all Read/Write/Edit operations on shared memory must start with the role dir.
- **Retired** in `developer/shared/MEMORY.md`: the `feedback_zotero_biorxiv.md` warning entry — fix landed via [PR #320](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/320), warning is stale.

**Process notes:**

- The `pytest --version` hang was the only real surprise of the session. Strong candidate for tomorrow's first warm-up if the `pytest.ini` fix is XS-shaped (it is — single config file + a CLAUDE.md note). Worth doing before someone else hits the same wall on a different test file.
- Bot review iteration on [PR #320](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/320) was tight — single round, both nits genuinely useful (the inline comment was over-explaining the `_is_preprint` choice; the medRxiv test closes a real coverage gap). Same shape as [PR #301](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/301)'s iteration on 2026-05-08 — fast, specific, verifiable.
- PM's morning closure-audit produced 2 false positives in this session (rationales-already-present and [Issue #215](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/215)-entry-already-present). Both flagged in standup follow-ups with audit-format suggestions (substring grep on `**Priority rationale:**`; raw GitHub blob for closure audit instead of local clone). PM acked the resolutions.

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
