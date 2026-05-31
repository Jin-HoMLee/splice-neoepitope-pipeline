# Lab Notebook — Developer

Per-role lab notebook for Developer sessions. Started 2026-05-06 as part of the lab-notebook split (see `lab_notebook.md` freeze note for rationale).

Format and rules unchanged from the unified notebook — see `shared/feedback_lab_notebook.md`. New date sections at the TOP; new time sections at the TOP of the date block; entries are immutable once committed.

---

## 2026-05-31 — surface pre-cap presenter total for the strong-presenter truncation notice ([PR #579](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/579) closes [Issue #226](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/226))

**What:** The strong-presenter HTML table is capped at `TOP_CANDIDATES_LIMIT=10` rows, but a reader couldn't tell whether "10" was the true count or a silent truncation. Fix persists two bookkeeping metrics in `report.tsv` under the `mhc_prediction` stage — `top_candidates_total` (pre-cap `len(strong_df)`, captured after the quality gate but before `.head(LIMIT)`) and `top_candidates_capped` — and the renderer [`_build_strong_table_html_from_top_candidates`](../../workflow/scripts/generate_report.py) gained `total`/`capped` params that append `"Showing N of M candidates."` when `total > capped`. Chose to persist the count in the long-format `report.tsv` rather than expand the capped `report_top_candidates.tsv` artefact — the capped TSV stays exactly what its name says.

**Why:** parity gap left by the original cap — the capped artefact discarded the information needed to caption itself honestly. Old, low-priority Issue; **verified still valid before starting** (docstring still referenced #226, cap still 10, no addressing commit/PR).

**What I learned / notes:**
- **Proactively caught a leak the RED run confirmed:** `_load_report_tsv` casts *every* `mhc_prediction` metric to int and `_presenter_counts_html` excluded only `total_predictions` — so the two new bookkeeping keys would have leaked into the per-class presenter table as fake presentation classes. Fixed by promoting the exclusion set to a module constant `_MHC_NON_CLASS_METRICS` (frozenset), used by both the loader threshold-skip and the presenter-counts filter/early-exit guard.
- **Foot-gun, now a memory:** during a mutation-test I layered a temp "break the wiring" edit on top of the *uncommitted* implementation, then ran `git checkout -- generate_report.py` to restore — which reverts to HEAD and wiped the entire uncommitted implementation (not just the temp break). Re-applied all 5 source edits from the captured diff. → role memory `git-checkout-discards-uncommitted-work` (use inverse-Edit / `cp` backup / commit-first instead). The mutation check itself succeeded: the over-cap E2E test correctly failed when the call-site wiring was severed, proving it meaningful.
- **Review:** [`@claude review`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/579#issuecomment-4585025420) → clean **Approve**, no bugs; confirmed each design decision (capture point, both-non-None + strict-`>` guard, the leak fix, graceful degradation on missing/old `report.tsv`). One low-severity note (no dedicated E2E test for "old `report.tsv` loaded but lacking the new metrics") declined — that renderer path is already covered at the unit level (`test_no_truncation_notice_when_total_unknown`), and with the new code both metrics always write together when `pred_df` is non-empty. 10 new tests; full suite 370 passed; all 3 CI checks green.
- **Non-issue surfaced + dismissed:** an adversarial diff review flagged the `mhc_prediction` schema row in `research/lab_notebook.md` (~line 1884) as a stale doc-sync gap. The user clarified that file is **historical/append-only** — false positive in context. → project memory `lab-notebook-md-is-historical`.

## 2026-05-30 — post_gh_pr_create matches `gh pr create` after a `VAR=value` prefix ([PR #575](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/575) closes [Issue #561](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/561))

**What:** `matches_pr_create()` in [post_gh_pr_create.py](../../.claude/hooks/post_gh_pr_create.py) skipped the auto-board when `gh pr create` was preceded by a shell env-assignment (`VAR=x gh pr create`, or a `VAR=x` line then `gh pr create`). The fix skips leading `VAR=value` assignment tokens (regex `^[A-Za-z_][A-Za-z0-9_]*=`) at each command-start position, mirroring the harness's subcommand-aware `Bash(gh *)` matcher. `shlex` (posix mode) strips the quotes (`B="x"` → `B=x`) and consumes newlines as whitespace, so the single fix covers both the same-line and newline forms. +5 unit tests.

**Why:** root cause — `shlex` emits `&&`/`;`/`|` as pure-punctuation separator tokens that reset `at_command_start`, but a bare `VAR=x` assignment token does not, so the following `gh pr create` was never seen as a command start. [PR #560](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/560) hit this and had to be boarded by hand.

**What I learned / notes:**
- **Dogfooded in production:** opened [PR #575](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/575) via `DOGFOOD=1 gh pr create …` — the exact previously-broken form — and the hook auto-boarded it + set Status "Ready for review", confirming the fix end-to-end (fire-log line `pr:575 board-add+status:Ready for review`).
- Bot review approved; per its note, added a one-line comment that an `export VAR=x` builtin prefix is deliberately **not** stripped (rare for `gh pr create`, out of scope).
- CI `ci-tools-pytest` is red on an **unrelated** stale live-test (`test_recheck_milestone.py::TestLiveIntegrationSmoke`): milestones #3/#5 legitimately flipped to `[UPDATE NEEDED]` as the calendar advanced (#5 due 2026-05-31). Non-required check; required checks (`pipeline-pytest`, `pipeline-snakemake-dry-run`) pass. Flagged the live recurrence on [Issue #506](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/506) (already P2; the property-based-refactor fix lives there).

## 2026-05-29

### 21:30 UTC — Editor: Developer

**Headline:** [PR #564-vehicle](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/564) — documented a **confirmed Claude Code hook-loading gotcha** in CLAUDE.md § GitHub Safety Wrappers, resolving the auto-board-miss investigation from the 20:45 entry.

**Root cause (systematic debugging, evidence-first):** the `post_gh_pr_create` auto-board hook missed both [PR #560](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/560) and [PR #562](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/562). A trace probe (unconditional line at `main()` entry + a read-only `gh pr view` tool call) proved the hook **never fired in this session** — even though `settings.json` is valid and the matcher correct. A fresh session ran the identical probe and the hook **fired** (`/tmp/pgpc.log` written). The only difference: this session had **edited `.claude/settings.json` mid-session** (wired `check_gh_issue_develop_parent` during PR #560). Conclusion: editing `settings.json` mid-session silently deactivates ALL hooks until restart — no hot-reload. So the misses were a session-loading artifact, **not** a code defect; [Issue #561](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/561)'s `VAR=` matcher bug stays scoped as a separate latent fix (deferred to a fresh session).

**Why doc-worthy:** it silently disables *existing* safety guards too (the `@claude` blocker was dead in this session after my edit). Anyone wiring a hook hits this. Captured as a ⚠️ callout + a role memory. Two doc-process notes: hooks fire only on Claude tool calls (not human terminal commands), and the `/hooks` menu's "create" screen is just its entry UI, not proof of absence.

---

### 20:45 UTC — Editor: Developer

**Headline:** [PR #562](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/562) (closes [Issue #559](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/559)) — merge-time **stray-closing-keyword gate**, the companion to the 20:13 `gh issue develop` guard. Together they close the [PR #543](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/543) → [parent Issue #538](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/538) incident class from both ends: create-time (don't branch off a parent) + merge-time (don't let a stray `close|fix|resolve` + `#N` in the squash body silently shut an unintended Issue).

**What + why (non-routine):** a closing keyword in a *commit-message body* auto-shuts the Issue on merge (the squash commit inherits it), but `closingIssuesReferences` surfaces only PR-**body** link edges — so the existing pre-merge check passed clean and missed exactly the vector that hit PR #543. New detector [`tools/ci/stray_closers.py`](tools/ci/stray_closers.py): a pure, word-boundary-aware `find_stray_closers()` (so `disclose`/`prefix`/`closing` don't match but `auto-close` does) + a thin `gh`-fetching CLI, wired as gate 4 in `scripts/audit_and_merge.sh`. Fails open on interpreter/gh error — never blocks a legitimate merge.

**The most satisfying verification — the gate catches its own origin story:** ran the CLI against the real incident **PR #543** and it flagged the live `close #538` from the commit body (exit 1), plus a *second* stray `Closes #7` I hadn't known about; clean (exit 0) on the just-merged PR #560. So this gate would have blocked the original incident.

**Review (`@-claude`):** no blocking issues. Accepted two: (1) tightened the separator `\s*` → `[ \t]*` so a keyword→`#N` match can't bridge the `"\n"` that `assemble_squash_text` joins fields with — a cross-field false positive GitHub itself wouldn't act on; (2) promoted the "word between keyword and `#N` = no match" invariant to a first-class `find_stray_closers` assertion. Declined two (both reviewer-acknowledged non-blocking): `owner/repo#N` coverage (single-repo project) and CLAUDE.md prose density.

**Self-discipline note (the irony):** had to scrub my own commit messages + PR body so they carry only the intended `closes #559` — a stray `auto-close #538` adjacency there would have made this PR trip the very gate it adds. The `auto-close #538` strings *inside* the files are safe (the gate scans PR body + commit messages, not the diff). 30 tests; full `tools/ci` suite green. **Also surfaced:** the auto-board hook ([PR #558](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/558)) missed both PR #560 and PR #562 — and #562 had no `VAR=` prefix, so [Issue #561](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/561)'s diagnosis is narrower than the real fault; the hook needs a broader look.

---

### 20:13 UTC — Editor: Developer

**Headline:** [PR #560](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/560) (closes [Issue #549](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/549)) — `PreToolUse` guard refusing `gh issue develop <N>` on parent/epic Issues (`subIssuesSummary.total > 0`). Rung-3 mechanism for the [PR #543](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/543) → [parent Issue #538](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/538) orphaning slip (2026-05-28); sibling of `check_at_claude.py`. Fails open on every uncertain path — only a *confirmed* parent denies, so a `gh`/network hiccup never blocks a legitimate develop.

**Why this entry exists (non-routine) — two verification catches around the hook-matcher boundary:**

#### The bot-review "gap" that wasn't — settled by docs, not assumption
Review flagged the narrow `if: Bash(gh issue develop *)` matcher as possibly missing compound commands (`git fetch && gh issue develop N`), since `develop_args()` handles separators (and `test_after_separator` covers it) but the `if` filter "might be start-anchored." Verified against the Claude Code docs (hooks "`if` Condition" + permissions "Compound commands"): the `if` field uses permission-rule `Bash()` semantics, which are **subcommand-aware** — splits on `&&`/`||`/`;`/`|`/newlines, strips `VAR=value` prefixes, fires if *any* subcommand matches. So the compound case **does** fire; the test path is reachable, not dead code. Declined the change, pushed back with citation. Lesson: a matcher-semantics claim is verifiable from docs — don't widen a matcher defensively on an unverified premise.

#### Mis-diagnosed, then verified: the sibling [PR #558](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/558) auto-board hook has a real `VAR=` bug
The new auto-board hook didn't board PR #560. My first read (told to the user) blamed the harness `if` matcher being start-anchored — **wrong**, per the verification above. Empirically tested `matches_pr_create`: `B="x" gh pr create` and the `VAR=`-newline-`gh pr create` form both return `False`, while plain and `&&`-separated return `True`. Root cause is the hook's *own* matcher: `&&`/`;`/`|` are emitted by `shlex` as pure-punctuation separator tokens that reset command-start, but a bare `VAR=x` assignment token does not — so the trailing `gh pr create` is never recognized. Filed [Issue #561](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/561) (follow-up). Lesson: when correcting a wrong claim, verify the replacement empirically before asserting it too.

**Verification:** 33 tests — pure matcher / selector-extraction / repo-parse + orchestration with monkeypatched `gh` I/O + subprocess fail-open; full `tools/ci` suite 155 green; all 3 CI checks green. **Live smoke:** parent #538 (number + URL form) refused with a self-explanatory message; leaf #549 + a non-develop `gh issue view` allowed silently; fire-log line written only on deny (`.claude/hook_fires.jsonl`, gitignored). **Bot review:** item 1 (matcher) declined with docs citation; items 3 (`--repo=HOST/owner/repo` test gap) + 5 (deny-message wording legibility) fixed in `1b7b4a3`; items 2 (`from __future__` non-issue for a standalone hook) + 4 (`-F num` typed-int) confirmations.

---

### 15:10 UTC — Editor: Developer

**Headline:** [PR #558](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/558) (closes [Issue #550](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/550) — `PostToolUse` hook auto-boards a created PR + sets Status). The **rung-3 mechanism** for the board-Status-flip slip that left [PR #6](https://github.com/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/pull/6) + [PR #543](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/543) at Backlog 2026-05-28 — the rule lived in `feedback_project_board.md` (memory, rung-2) and didn't fire on the create beat. Sibling of `scripts/audit_and_merge.sh` (closure-ritual gate on merge).

**Two design pivots + a self-caught bug — the reason this entry exists (non-routine):**

#### Pivot 1 — draft-awareness (caught in review of my own design)
The first cut flipped *every* created PR to "Ready for review". Wrong: a `gh pr create --draft` is opened mid-In-progress (CI / shareable URL), not for review. Fixed to branch on the PR's authoritative `isDraft`: draft → "In progress", non-draft → "Ready for review". Either way it lands ON the board — board-*absence* was the real failure. (User flagged the always-Ready assumption before wiring.)

#### Pivot 2 — the hook caught its OWN false positive, live
The instant it went live this session it fired on my **review-reply comment** to the bot — that comment's `--body` contained the example `… git push && gh pr create …`, and the original `matches_pr_create` regex searched the *raw command line*, so the `&&`-prefixed occurrence **inside the quoted body** matched. It then parsed #558 from the comment URL and re-flipped it (idempotent + one spurious fire-log line). Exact substring-in-body class the `@-claude` mention guard exists for. Fix: tokenize with `shlex` (posix + `punctuation_chars`) so quoted args stay single tokens; only a real `gh pr create` segment (command-start or after a true shell separator) matches; unbalanced quotes fail safe. Verified `cd foo;gh pr create` → match, `gh pr comment … "… && gh pr create …"` → no match.

#### Self-modification gate
The harness **refused** my first edit to `.claude/hooks/` + committed `settings.json` — correct: wiring a hook that fires for all sessions/roles is agent self-modification, and "ok gogo" didn't specifically authorize it. Got explicit sign-off, then proceeded.

**Verification:** 20 tests (pure matcher/parse/should_track/status_for_draft/extract_output + subprocess fail-open/no-op); full `tools/ci` suite green; CI green. **Live dogfood:** hook flipped #558 itself Backlog → Ready for review and wrote its fire-log line — non-draft path end-to-end. Bot review: one real finding (dead `has_project_flag` → removed) + fire-log key `issue`→`pr` clarity; two non-blocking notes declined with reasons.

---

### 14:02 UTC — Editor: Developer

**Headline:** [PR #556](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/556) (closes [Issue #555](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/555) — closure-audit bot honors a routine-ship lab-notebook skip marker). Closes the **enforcement half** of [Issue #483](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/483): that Issue made routine single-PR-closes-single-Issue notebook entries optional *and* dropped the `audit_and_merge.sh` notebook gate — but never taught the post-merge closure-audit bot the exemption. The bot's only skip was **file-based** (`is_exempt`: all-glossary/notebook diffs), so any routine PR touching code still drew a "Lab notebook entry missing" gap. Rule and enforcement contradicted each other.

**Decision (AC 1):** deterministic **opt-out marker** `<!-- skip-lab-notebook: routine -->` in the PR body, over the two alternatives — heuristic trigger-detection ([Issue #483](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/483) itself called this "too brittle") and decline-and-document. The marker mirrors the existing `❎ … deferred` AC-deferral convention: author-controlled, greppable, no false negatives.

**Work shipped:**

- `skip_lab_notebook(pr_body)` — case-insensitive, whitespace-tolerant regex `<!--\s*skip-lab-notebook\b[^>]*-->`. `fetch_pr` now requests `body`; `audit_pr` short-circuits **only** the notebook-check block (`and not skip_lab_notebook(...)`). AC-checkbox + Priority-rationale checks unaffected.
- `audit_issue` deliberately **untouched** — the marker is PR-body-scoped, and decomposition Issues closed without a PR still require an entry, so no skip should reach that path.
- TDD: 6 tests RED→GREEN (skip-honored: value/bare/whitespace+case; skip-absent: no-marker/empty/unrelated-HTML-comment). 32 `tools/ci` tests green locally + in CI.

**Process notes:**

- **Cross-repo AC (item 3).** The `shared/feedback_lab_notebook.md` cross-reference lives in the separate `claude-personas` repo, whose working tree is mid-flight on the MEMORY-slim epic (#538–#542). Made the edit but left it **uncommitted** per the maintainer's call — to be folded into that epic rather than committed onto a parallel session's dirty branch.
- **Review.** Bot review found no bugs; applied the one doc nit (marker value must be `>`-free, since `[^>]*` stops at the first `>` — c75c68b), declined the suggested `audit_pr` integration test (YAGNI for a 2-line short-circuit; the suite is intentionally lean on live-`gh` paths).
- **Not self-exempting.** This is a meta-decision/workflow-rule session → non-routine → entry required. Did *not* place the new marker on this very PR.

---

### 10:55 UTC — Editor: Developer

**Headline:** [PR #554](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/554) (closes [Issue #495](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/495) — closure_audit accepts a PR # **or** any closing-Issue # in the lab-notebook check) shipped TDD-first: 5 new tests, 96 total. Read-side fix for the [PR #494](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/494) false positive — that entry named only its closing Issue (#484), not the PR, and tripped the bot's PR-number-only check.

**Work shipped:**

- `check_lab_notebook` gains an additive `also_accept: Collection[int] = ()`; a day-block passes on `#number` OR any `#also_accept`. Threaded through `check_lab_notebooks_for_issue` → `collect_notebook_gaps`; `audit_pr` now passes the PR's `closingIssuesReferences` numbers. `audit_issue` is unchanged (no PR number).
- Default `()` preserves prior behavior + the gap message verbatim, so the 10 pre-existing notebook tests stay green as a regression guard.

**Process notes:**

- **This entry is the fix demonstrating itself.** Written before merge, it references the closing Issue (#495), not the PR number — the exact "entry predates the PR" shape that tripped [PR #494](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/494). Pre-fix, the post-merge closure-audit bot would have gap-commented here too; with this PR merged, it accepts the entry via `also_accept=[495]`. Eating our own dog food.
- **Review pushback held up under verification.** The bot's lone nit (drop the `(#495)` doc refs, citing CLAUDE.md) rested on a rule not present anywhere in the repo; the file already cites `see #524` / `see #325` for rationale. Kept the pointer, reformatted to the house `(see #N)` form for consistency — `grep`-verified the citation was absent before pushing back.

---

## 2026-05-28

### 20:38 UTC — Editor: Developer

**Headline:** [PR #545](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/545) (closes [Issue #524](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/524) — closure-audit any-role-satisfies check for multi-role Issues) shipped TDD-first with 7 new tests, then a bot-caught dedup regression got fixed in a follow-up commit + 3 more tests. Origin: [PR #518](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/518)'s false-positive (Scientist's notebook satisfied the reference but the bot only audited Developer's because `resolve_roles` returned `role_labels[0]` after `sorted()`).

**Work shipped:**

- `resolve_roles` returns `list[set[str]]` per-Issue (was: flat dedup `list[str]` that dropped all-but-first-alphabetical role per Issue). Positionally aligned with the input list so callers can zip with the original Issue list.
- New `check_lab_notebooks_for_issue(roles, date, n, notebooks)`: any-role-satisfies semantic — returns `None` if at least one role's notebook references `#N`, otherwise a single combined `(role_label, description)` tuple (`"developer or scientist"` joined for the multi-role case).
- `audit_pr` / `audit_issue` callers refactored through a small `collect_notebook_gaps` helper that dedupes by `frozenset(roles)` — added after [bot review](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/545#issuecomment-4567982400) caught a duplicate-gap regression (a PR closing 2 Issues sharing a role label would print the same gap line twice).
- `_load_notebook` gained `encoding="utf-8"` (defensive, applies to all CI runners regardless of process locale).
- Test count: 18 → 21 (7 new for the fix, 3 for the dedup helper); 11 pre-existing regression tests preserved.

**Process notes:**

- **TDD discipline paid off twice.** Red→green on the initial 7 tests confirmed the helper shape was right before any caller touched it. Then when the bot caught the dedup regression, the same testable-helper pattern made adding `test_collect_notebook_gaps_dedupes_identical_role_sets` trivial — no subprocess mocks needed for `audit_pr`/`audit_issue` orchestration because the loop body is now a pure function.
- **Bot caught a regression I didn't.** The old flat-dedup `resolve_roles` had silently dropped a second function — preventing duplicate gap lines when multiple closing Issues shared a role. Refactoring to per-Issue restored the original bug **but** broke this latent behavior. Easy to miss in unit tests; bot's "what if N closing Issues all have role:developer?" reasoning was the right adversarial-thinking lens.
- **Receiving-code-review on the `[roles] = ...` style note:** bot suggested `roles = resolve_roles([labels])[0]` instead of unpacking. I left the unpack as-is initially (it asserts length), but the dedup refactor made the question moot — `audit_issue` now goes through the same `collect_notebook_gaps` path as `audit_pr`, no more single-element unpacking.

---

### 19:30 UTC — Editor: Developer

**Headline:** [PR #519](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/519) (closes [Issue #514](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/514) — `torch<2.5` pin lift via cu126 pip channel) finished out: cloud verification on `neoepitope-pipeline` P100 yielded a **bit-identical** MHCflurry chr22 output vs the pre-bump baseline (81 rows, all structural columns identical, numeric drift at FP-noise level), `@-claude review` returned LGTM in 1m 58s with no blocking findings. Merge-ready; unblocked by [PR #526](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/526) earlier today (see 18:00 UTC entry).

**Work shipped:**

- [PR #519](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/519) cloud verification on `neoepitope-pipeline` (19:02–19:06 UTC): new conda env hash `9681b0f6b9e5ce183dd6efeadcefa736_` resolves `torch 2.12.0+cu126` + `mhcflurry 2.2.1` cleanly via snakemake `--use-conda`. `_has_gpu()` smoke test returns `tensor([0., 0.], device='cuda:0')` on Tesla P100-PCIE-16GB (driver 550.90.07). `--forcerun run_mhcflurry` against this morning's cached upstream (54 peptides × 6 alleles) produced **81 rows = 81 baseline rows** with all structural columns (peptide, best_allele, presentation_class, n_strong_alleles) **bit-identical**, class distribution unchanged (78 non / 3 weak / 0 strong), and FP-determinism drift only (`ic50_nM` max abs 5.6e-03 nM / 0.0001% rel; `presentation_score` max abs 1.3e-07; `presentation_percentile` / `genotype_presentation_score` / `best_presentation_percentile` zero drift). Pre-bump baseline at `gs://splice-neoepitope-project/results/patient_001_test/predictions/mhc_presentation.tsv` from the [PR #526](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/526) chr22 run earlier today — free baseline, no extra run needed.
- PR body updated with verification record; flipped to ready-for-review; both [PR #519](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/519) + [Issue #514](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/514) project board Statuses moved to "Ready for review".
- `@-claude review` posted 19:10:53 UTC; bot LGTM at 19:11:23 UTC (1m 58s elapsed). Two non-blocking notes — `torch ==2.12.0+cu126` spacing is consistent with sibling `mhcflurry >=2.0`; CLAUDE.md PEP 440 sentence "accurate, if a touch loose" — both flagged no-action per the bot itself.

**Process notes:**

- **`--rerun-triggers mtime` does not treat env-yaml mtime as a rerun trigger when outputs already exist.** First snakemake invocation built the new conda env (correctly detecting the python.yaml content change → new hash dir), then said `Nothing to be done` because output mtimes were newer than input mtimes. `--forcerun run_mhcflurry` was required to actually exercise the new env on existing inputs. Worth knowing for future env-only bumps: rebuilding the conda env is necessary but not sufficient to revalidate downstream output.
- No full chr22 `--forceall` needed to validate the other 11 `python.yaml`-consuming rules: they consume CPU-only deps (pandas/numpy/biopython) unchanged by this PR and covered by the 350-test pytest suite. Only `run_mhcflurry` exercises torch.

---

### 18:00 UTC — Editor: Developer

**Headline:** [Issue #522](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/522) (P100 cloud-run blocker) ready to **close** — three-regression compound failure resolved via image-family pivot, plus a fourth regression discovered + fixed mid-session. [PR #526](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/526) ships 7 commits: original apt-cache 404 fix + 535-keep doc flip (both superseded) + dropped-570 install + smoke-test gate + cu124 image pin with `install-nvidia-driver=True` metadata + FRESH_BOOT-conditional 5-min sleep + bot-review hardening (extended retry, version assertion, fallback hint). chr22 ran end-to-end in ~29 min wall (run #3, 17:13–17:42 UTC). [PR #519](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/519) ([Issue #514](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/514) torch pin lift) was blocked on this all along — now unblocked.

#### Yesterday's plan executed; new regression discovered mid-execution

The pause memo wrote a 6-step plan (image-family probe → VM probe → patch script → flip CLAUDE.md → re-run chr22 → unblock #519). Steps 1–4 went smoothly: `gcloud compute images list` on `deeplearning-platform-release` showed only `common-cu129-*` + `pytorch-2-9-cu129-*` (all nvidia-580) as current families, so I extended the search to deprecated-but-READY images. The probe revealed `pytorch-latest-cu124-v20250327-ubuntu-2204` ships proprietary **550.90.07** via DLVM's `install-driver.sh` on first boot — the script's `machine_type =~ ^n` guard picks the **closed** kernel module variant for n1 machine types, which is exactly what Pascal needs (open variants require GSP firmware Pascal lacks). Probe VM: `Tesla P100-PCIE-16GB, Driver 550.90.07, CUDA 12.4` clean.

But step 5 surfaced a **fourth** regression: chr22 runs #1 and #2 crashed at `setup_vm.sh` step 1/8. Run #1 disconnected during `apt-get update`; run #2 disconnected during a held-SSH `cloud-init status --wait` block I added to "fix" run #1. Root cause: Ubuntu's `unattended-upgrades` runs in the background on first boot, and when it bumps `openssh-server`, sshd restarts — killing any held SSH session (and even racing SSH key handshake on new sessions).

#### "Fix #1 was racy, fix #2 just waited"

The fix progression illustrates the trap with these first-boot races: the held-SSH wait (commit `87961bc`) had the same vulnerability it was meant to mitigate. The actual fix (`9a1c134`) replaced it with a **host-side `sleep 300`**, conditional on `FRESH_BOOT=true` (VM was just created, not warm-restarted). Crude, but reliable: cloud-init + unattended-upgrades only fire on first boot, so warm restarts skip the sleep entirely — zero cost. Run #3 then passed end-to-end; `report.html` + 4 TSVs in `gs://splice-neoepitope-project/results/patient_001_test/reports/`.

Lesson worth keeping: **on Google's DLVM images, never hold a long SSH session during the first ~5 minutes after boot.** Short polls work because each is a fresh connection (the smoke test's 18-attempt loop survived fine even during runs #1/#2). The right gate isn't "wait for cloud-init done" via SSH — it's "wait on the host side for the chaos to settle."

#### Bot review: 3 hardening suggestions, all defensible

claude bot LGTM'd CLAUDE.md and yesterday's pause-entry lab notebook write, no blocking issues on the script. Three optional suggestions:

1. **Extended driver-verify window on fresh boot** (18 → 36 attempts). Run #3 passed with 18; bot is over-cautious. But trivial insurance.
2. **`DRIVER_MAJOR >= 575` invariant assertion**. Bot's stated rationale was incorrect — the existing smoke test already catches open-driver-on-Pascal (`nvidia-smi` returns empty when the driver fails to initialize the GPU). But the explicit assertion still has value as a self-documenting requirement ("script requires driver < 575"). Applied on those grounds rather than the bot's.
3. **Inline comment fallback hint** pointing at CLAUDE.md custom-bake instructions if Google ever deletes the deprecated image.

All three applied as one commit (`4864f35`), no re-verification needed (only error paths + a doc comment).

#### Cost + timing

Total session cloud spend: ~$2.50 (probe VM ~$0.50 + three production VM runs at ~$0.50–1.00 each). About 3.5 hours of session work for ~30 min of actual chr22 wall-clock. The decomposition (image probe → strategy → script change → 3 chr22 attempts → bot review → final hardening) was the right shape for this multi-layered regression — a one-shot patch would have left the unattended-upgrades race undiscovered until the next on-call.

#### What this unblocks

[PR #519](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/519) ([Issue #514](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/514)) — torch pin lift via cu126 pip channel — had cloud verification ACs gated on a working P100 pipeline run. With #522 closed, #519 can collect that evidence and merge.

---

## 2026-05-27

### 22:00 UTC — Editor: Developer

**Headline:** [Issue #522](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/522) (P100 cloud-run blocker) deepened into a three-regression compound failure — Ubuntu repo collapsed `nvidia-headless-570-server` into a 580 wrapper; Ubuntu archive's userspace 535 churned past the image's matched kernel-module version; Google pushed a new image build TODAY (`v20260527`) that ships only 580-open which doesn't support Pascal. [PR #526](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/526) lands two partial fixes (`apt-get update` for the original cache 404 + drop the broken 570 install + add a `nvidia-smi` smoke test as a permanent fail-fast gate) but the full fix is blocked on an image-family decision deferred to tomorrow. Detailed investigation log in the [Issue #522 thread comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/522#issuecomment-4558114563).

#### Three regressions, not one

What I filed as a 1-line apt-cache 404 fix this morning turned into a multi-hour investigation. The compounding:

- **Ubuntu repo level:** `apt-cache depends nvidia-headless-570-server` returns `Depends: nvidia-headless-580-server`. Upstream reduced 570-server to a wrapper for 580. `--no-install-recommends` doesn't bypass hard `Depends:`, so the script's "downgrade to 570" cascades 580 in unavoidably — and 580 drops SM 6.0 (Pascal), killing P100.
- **Ubuntu archive churn:** the image's `linux-modules-nvidia-535-server-<kernel>` package ships kernel module at version `535.288.01`. The userspace `nvidia-utils-535-server` advanced in the archive to `535.309.01`. `apt-get update` pulls the userspace forward without a matching kernel-module package, producing irrecoverable mismatch. `nvidia-kernel-source-535-server` doesn't ship a usable `dkms.conf`, so DKMS can't rebuild.
- **Google image-family churn:** the image family `common-cu129-ubuntu-2204-nvidia-580` was pushed today (`-v20260527`) to ship only 580-open variant (`linux-modules-nvidia-580-server-OPEN-...`). On P100 (`10de:15f8`), `nvidia.ko` loads but immediately fails probe — the open driver is Turing+ only per NVIDIA's policy.

#### "Trust the image" was the wrong frame

Mid-session I committed a strategy shift `18baa10`: drop the 570 downgrade entirely, trust the image's pre-installed driver, fail-fast via `nvidia-smi` smoke test on script start. That worked under yesterday's image build (which still shipped 535 baked in) — but Google's `v20260527` push happened during the session, and the new baseline has no Pascal-compat fallback at all. The "trust the image" framing assumed image-family stability that doesn't exist.

The fail-fast smoke test still has standalone value: it correctly turned an opaque downstream `docker run --gpus all` failure (`nvml error: driver not loaded` deep in a TCRdock rule) into a clear actionable script-level error at VM-start time. Keeping that change tomorrow regardless of driver path.

#### Multiple gcloud client crashes during the session

Three separate `gcloud crashed (ConnectionError): Remote end closed connection without response` failures during VM create + describe + ssh during this session. Once during VM create (the script's request landed server-side and went STAGING → TERMINATED via the cleanup trap race). Once during a multi-line ssh investigation command. Not unique to today's investigation — googlecloudsdk 475.0.0 has had intermittent IAP-tunnel websocket disconnects for months. Worth keeping in mind: retry transient gcloud failures rather than treating them as terminal.

#### Branch-name memory paid off

Filed [Issue #522](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/522)'s branch via `gh issue develop 522 --name 522-fix-driver-apt-cache-stale --base main --checkout` — passing `--name` explicitly with an ASCII-only slug per the morning's `feedback_gh_issue_develop_ascii_only.md` memory. The Issue title contains no Unicode this time so the auto-slugifier would've been fine — but applying the rule reflexively keeps the habit warm for next time.

#### Plan for tomorrow's session

1. Image-family investigation: `gcloud compute images list --filter='family ~ "common-cu126|pytorch-2-(1|2|3)-cu12(1|4)"' --project=deeplearning-platform-release` to find a family that ships Pascal-compatible driver (≥470 and <575).
2. Spin up a probe VM (10 min) and confirm `nvidia-smi` works out-of-the-box on a candidate.
3. Patch `run_cloud_gpu.sh:45` to use the new family; keep the smoke-test gate from `18baa10`.
4. Re-rewrite the CLAUDE.md driver bullet to match the new strategy.
5. Re-run chr22 test as the AC verification.
6. Once [PR #526](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/526) lands → unblock [PR #519](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/519) ([Issue #514](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/514) torch pin lift) which has been blocked on this all along.

VM stopped at session-end to avoid cost accrual.

---

### 18:30 UTC — Editor: Developer

**Headline:** Attempted [Issue #514](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/514) (lift `torch<2.5` pin via cu126 pip channel) — code changes landed clean as [PR #519](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/519) (Draft), but cloud verification surfaced two orthogonal infra blockers — filed as [Issue #521](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/521) (P2, script ZONE default drift; closing via [PR #525](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/525)) and [Issue #522](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/522) (P1, NVIDIA-580 archive 404 blocks all cloud runs). [PR #525](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/525) sweeps 5 stale `west1-b` references missed by [Issue #516](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/516)'s CLAUDE.md-only fix. Replaced rebranched-and-reopened [PR #523](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/523) after the Claude Code Action rejected the auto-generated `→`-containing branch name.

#### Cascade pattern: pin-lift → blocked by infra drift → infra drift → blocked by infra failure

[PR #519](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/519) verification path: `python.yaml` change committed → cloud-run attempt 1 fails at VM-find (script defaulted to pre-migration `west1-b`, hit `GPUS_ALL_REGIONS quota exceeded`) → filed [Issue #521](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/521); attempt 2 with `--zone europe-west4-a` fails at the driver-install step (9× 404 on `nvidia-headless-580-server 580.126.20` from Ubuntu's `europe-west4.gce.archive.ubuntu.com` mirror; image ships 580 baked in, setup overlays 570 per CLAUDE.md, apt cascade wants both, 580 archive is gone) → filed [Issue #522](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/522). [PR #519](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/519) stays Draft until [Issue #522](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/522) lands; [PR #525](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/525) closes [Issue #521](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/521) standalone.

#### Branch-name validator caught the Unicode `→` from `gh issue develop`

[Issue #521](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/521)'s title contains a literal `→` (`west1-b → west4-a`). `gh issue develop --base main --checkout` slugified it as `521-...-west1-b-→-west4-a` — git accepts it, but the Claude Code Action's branch-name validator rejects: `Invalid branch name: "..." Branch names must start with an alphanumeric character and contain only alphanumeric characters, forward slashes, hyphens, underscores, periods, hashes (#), plus signs (+), or commas (,).` `@claude review` failed at the checkout step ~4s in, leaving a stale "is working..." comment. Workaround: `git branch -m <ascii-slug>`, force-push, close [PR #523](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/523) with forward pointer, reopen as [PR #525](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/525). New developer-role feedback memory `feedback_gh_issue_develop_ascii_only.md` indexed in MEMORY.md: pass `--name <ascii-slug>` to `gh issue develop` whenever the Issue title has non-ASCII chars.

#### west1-b sweep: scope expansion was correct, not gold-plating

Initial [Issue #521](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/521) AC named line 41 + a grep audit; first commit-pass hit `scripts/run_cloud_gpu.sh:41` + 2 lines in `docs/google_cloud_guide.md`. Full-repo `git grep` (post-PR-open) revealed two more operational defaults: `CLAUDE.md:207` (copy-paste SSH command in HISAT2 cache section, same drift class as the script default) and `research/glossary.md:73` (GCP definition stating pre-migration zone). Both are "anyone reading this today would be misled" surfaces, not historical narrative. Expanded scope to land them in the same PR — preferable to filing yet another follow-up Issue for an identically-shaped one-line fix. Bot review confirmed the expanded sweep was complete (no operational defaults pointing at west1-b remain).

#### Bot review's "unstaged working-tree revert" finding was a false positive

Bot flagged `git diff HEAD -- CLAUDE.md` showed a reverted commit `6905fb7` in the working tree of its checkout. Verified locally: empty diff, clean `git status -s`. Likely a stale-workspace artifact from the bot's failed run on [PR #523](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/523) carrying into [PR #525](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/525). Replied on the PR thread documenting. The other two bot findings (`west4-b` in docs examples — intentional; T4-quota comment — bot agreed with deferral) needed no code change.

#### Process: PR #519 → #522 → #519 unblock chain

[PR #519](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/519)'s cloud-verification AC line items (`_has_gpu()` smoke test + MHCflurry baseline diff on chr22 test patient) are gated on [Issue #522](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/522) — until any cloud run succeeds, no Pascal verification possible. [Issue #522](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/522)'s P1 priority reflects this multiplier on all GPU-touching work. Plausible fixes outlined in the Issue body (5 options: apt-mark hold, --fix-missing, version-pin to fetchable head, image-family bump, custom-image bake) — needs investigation in a fresh session to pick.

---

### 14:30 UTC — Editor: Developer

**Headline:** [PR #510](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/510) shipped, closing [Issue #508](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/508) — added `scripts/board_open_items.py` (paginated open-board listing with role/Status/Priority/Size filters) and promoted the long-existing board-pagination rule from `shared/feedback_board_queries.md` Reference tier into `shared/MEMORY.md` Always-in-effect. Two-layer fix following the established memory ladder (rung 2 inline + rung 3 mechanism). Bot review caught 3 valid robustness issues (GraphQL error envelope unhandled, `subprocess check=True` swallows `gh` stderr, argparse `choices=` missing on enumerated flags) — all applied as one-liners.

#### Slip postmortem: a documented rule sitting at Reference tier behaves the same as no rule

User asked me to "check the board for next best tasks". I ran `gh api graphql` with `first: 100`, returned 2 open Dev items, recommended next picks. User pushed back. Re-query with `pageInfo.hasNextPage` cursor loop → 24 open Dev items (out of 64 total open, 472 board items). Two compounding pitfalls: `first: N` is a cap (not a filter), and the board sorts Done items first — so an unpaginated query returned 97 Done items on the first page, leaving 3 open-state candidates after client-side filtering. The rule was documented at `shared/feedback_board_queries.md` since 2026-05-20, but indexed under "Reference (read when relevant)" — never loaded at session start. Effective rule presence was zero.

#### Memory ladder applied as designed

Per `shared/feedback_memory_escalation.md` ("on 'you forgot X': inline into Always in effect") + `shared/feedback_mechanism_over_memory.md` (repeat-failure shape → mechanism), the fix split cleanly: (1) **rung 2** — promoted the rule into `shared/MEMORY.md` Always-in-effect with telegraphic body + back-link to `feedback_board_queries.md`; (2) **rung 3** — shipped `scripts/board_open_items.py` so future sessions don't redo the cursor loop. Memory teaches the principle (covers cases the script doesn't fit — mutations, unusual fields), script is the canonical implementation for the most common case. **Complementary, not interim/permanent.**

#### Bot review findings — all 3 valid, all 3 robustness gaps

(1) `data["data"]` indexing into the GraphQL response would `KeyError` on `{"errors": [...], "data": null}` envelopes (auth failure, rate limit, partial failure). Bot suggested `if errs := data.get("errors"): print + exit`. (2) `subprocess.run(..., check=True, capture_output=True)` prints a Python traceback instead of the actual `gh` stderr ("Not Found", "requires authentication") — UX-hostile under failure. Fix: drop `check=True`, branch on `r.returncode`, print `r.stderr` to stderr. (3) `--status`/`--priority`/`--size` filters did silent exact-match; passing `--status "in progress"` (lowercase) silently returned 0 results. Free fix: argparse `choices=list(STATUS_ORDER)` etc. since the enumerations already exist for sort ordering. `--role` intentionally kept free-form (extensible label).

#### Pattern: separate-step author flow for review-trigger lifecycle worked clean

Shared Always-in-effect rule: when posting `@claude review`, author flips BOTH PR Status and linked Issue Status to `In review` in the same step as the comment. Done in parallel — three independent calls (gh pr comment + 2× updateProjectV2ItemFieldValue) — clean. After review fixes land, Status stays `In review` (covers both active review and awaiting-iteration), no flip-back to `Ready for review`. Merged → Done auto-set on both PR and Issue by the project workflow.

---

## 2026-05-26

### 12:20 UTC — Editor: Developer

**Headline:** [Issue #200](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/200) eval (Snakemake 8.x → 9.21.0 upgrade path) closed with verdict **defer-with-trigger** — bundle the bump into [Issue #66](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/66) (Google Batch executor migration) when that lands. Three breaking changes between v8.0 and v9.21 (LoggerPlugin v9.0.0, `--cores` explicit v9.2.0, runtime-minutes-not-seconds v9.20.0); none hit our codebase. Migration guide's "only 1 breaking change" framing is technically understated but the net impact for us matches their "hardly any user affected" claim. CLAUDE.md `--configfile` foot-gun persists in 9.x (argparse `nargs="+"` unchanged on `main`); the existing workaround stays valid post-upgrade. `snakemake-executor-plugin-googlebatch` v0.5.1 declares interface compat with 9.x via `snakemake-interface-executor-plugins ^9.0.0` but its dev-test pin is `snakemake = "^8.30.0"` — that gap is the strongest argument for bundling the bump with [Issue #66](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/66) (validation work happens once, not twice). Full findings: [closing eval comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/200#issuecomment-4543995431).

#### Why defer-with-trigger over upgrade-now or skip

Three legs: (1) zero current pain — none of the 3 breaks touch this codebase (no custom logger, `--cores` always explicit, zero `runtime=N` resource specs anywhere in `workflow/`); (2) zero current benefit — none of the 9.x feature additions (`--replace-workflow-config`, `--report-after-run`, Pixi integration, scheduler plugin interface, storage-plugin checksums, named multi-output caching) address a present pain point; (3) strong natural trigger — [Issue #66](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/66) needs to validate the googlebatch plugin against 9.x as part of its scope anyway, so the validation overhead is already budgeted there. **Upgrade-now** would be pure churn (low-risk but benefit-free); **skip-with-rationale** over-strong given the natural trigger is already on the board.

#### No sub-issue carve — trigger lives on existing [Issue #66](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/66)

Heads-up comment posted ([comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/66#issuecomment-4543996909)) capturing the bundle-with-upgrade recommendation. **Pattern: when a deferred task has a natural fire-trigger on an existing Issue, prefer dropping a heads-up there over filing a fresh tracking Issue** — keeps the dependency graph local to where it'll execute, avoids an Issue that mostly just points to another Issue.

#### Drift from the 2026-05-12 news_log claim — refined here, not corrected there

The 2026-05-12 news_log entry said "8→9 migration is small — only 1 breaking change (custom logger API)". Actual changelog count is **3** breaks, but the upstream migration guide also uses the "1 breaking change" framing, so the news_log tracked the guide rather than the source-of-truth `CHANGELOG.md`. Refining the count here (in this lab notebook entry + the closing eval comment) rather than amending the 2026-05-12 news_log entry — entries are immutable.

#### Closure-ritual sequencing on tick-with-deferral

[Issue #200](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/200) body had 5 Tasks; 3 actioned directly, 2 deferred to the [Issue #66](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/66) migration (chr22 integration test + `environment.yaml` pin bump — both belong in the actual upgrade PR, not the eval). Ticked all 5 with inline annotation on the deferred ones per the closure-ritual workaround (`- [x]` with a link to the carrier Issue). Audit-and-merge script parses by literal `- [ ]` count, not by deferral semantics, so explicit ticks + annotation is the correct pattern.

---

## 2026-05-25

### 19:37 UTC — Editor: Developer

**Headline:** [PR #477](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/477) shipped, closing [Issue #345](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/345) — per-sample raw regtools/STAR quantification renamed `alignment/<sample>/junctions.tsv` → `alignment/<sample>/raw_junctions.tsv` to disambiguate from patient-level `junctions/novel_junctions.tsv`. 10 files (9 in the rename commit + 1 README fix from bot review), pure rename, no behavior change. Local chr22 integration verified end-to-end (HISAT2 path, 4/4 steps, 155 records in `novel_junctions.tsv` = 151 tumor_exclusive + 4 normal_shared, matching the post-#370 baseline). Bot review caught one real miss + flagged two correctly-classified frozen-script references; one follow-up [Issue #478](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/478) filed for the stale `RESULTS.md:54` Scientist-territory reference.

#### Scope discipline: extending beyond the AC for a rename is correct if the goal is grep-cleanliness

[Issue #345](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/345)'s AC listed only `alignment.smk`, `filter_junctions.smk`, and `test_filter_junctions.py`. I extended to the two converter scripts (`bed12_to_junctions.py`, `star_sj_to_junctions.py`) + their unit tests because the Issue's stated motivation is **grep-cleanliness** (`find . -name junctions.tsv` returning mixed-meaning hits). Leaving 24 `junctions.tsv` strings in converter tests + docstrings would defeat that. The PR body called this out explicitly so reviewers could push back if they disagreed — bot didn't. **Pattern: when an Issue's AC is narrower than its stated motivation, extend to satisfy the motivation, and document the extension upfront in the PR body so reviewers can challenge it.**

#### Bot review missed in my own pre-PR grep: live README output-tree diagram

My pre-PR grep used `--include="*.py" --include="*.smk" --include="*.yaml" --include="*.yml" --include="*.sh" --include="Snakefile"` to find live references. README.md was excluded from the file-type filter — I only added `*.md` to the *exclusion* sweep (which filtered OUT frozen plans/specs/notebooks). Net effect: README.md:223 never appeared in any grep. Bot review caught it. **Lesson: when greping for rename impact, the include-filter must cover every file class where the renamed string is meaningful — for path-shape renames that's `*.md` (READMEs, docs), not just code.** No memory rule warranted yet (1st instance of this shape; bot caught it free); user and I agreed to wait for a 2nd instance before escalating.

#### Frozen `_regenerate_figures.py` scripts: pushed back on bot's soft "add a note" suggestion

Bot flagged two frozen-experiment scripts (`research/experiments/issue_225_*/figures/_regenerate_figures.py`, `research/slides/issue_393_*/figures/_regenerate_figures.py`) that hardcode the old path, marked informational-not-blocking, and suggested README notes. Skipped both: per the [CLAUDE.md frozen-experiment convention](CLAUDE.md), `research/experiments/issue_NNN_*` and `research/slides/issue_NNN_*` are point-in-time records. Retroactively rewriting frozen scripts weakens "frozen reference" semantics; accreting "this filename used to be X" stickers in experiment READMEs for every future rename creates a maintenance burden. A regeneration attempt fails loudly with `FileNotFoundError` — clear, recoverable failure mode. Bot agreed in framing (informational). **Pattern: when a "useful but not blocking" suggestion conflicts with an established codebase convention, default to the convention and document the reasoning in the reply — don't accrete stickers.**

#### CLAUDE.md:146-153 left as-is — generic placeholder, not a content reference

Bot's own analysis: CLAUDE.md:146-153 uses `results/.../junctions.tsv` as a generic illustrative path in the `--configfile` argparse gotcha example. The filename is incidental to the pitfall being demonstrated. Updating would be scope creep without payoff (the rename's grep-cleanliness goal targets `find -name`, not markdown content search). Left untouched per bot's recommendation.

#### Local integration run targeted only `filter_junctions` to keep verification scoped

For the chr22 belt-and-suspenders run, I targeted `results/patient_001_test/junctions/novel_junctions.tsv` rather than the full `report.html` to skip mhcflurry/tcrdock (heavy models / Docker — out of scope for a rename validation). The rename only touches the alignment → filter chain; CI dry-run + pytest already cover downstream rule-graph wiring. The full DAG is verified weekly via cloud runs anyway. **Pattern: for pure-rename PRs, scope the local integration run to the directly-affected rule's output, not the terminal `all` target.**

---

### 18:49 UTC — Editor: Developer

**Headline:** [PR #469](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/469) shipped, closing [Issue #63](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/63) — strict nf-core separation of `references/` (user-provided + downloaded data) vs `indices/` (pipeline-built alignment indices) vs `resources/test/` (small committed test fixtures + chr22 local-dev cache). 19 files across config, rules, setup_vm.sh, .gitignore, docs, and CI placeholder paths; idempotent migration block in [scripts/setup_vm.sh](scripts/setup_vm.sh) moves pre-#63 layout into place on existing VMs. Local chr22 end-to-end verified (11/11 rules, exit 0) BEFORE PR opened — the verification table covered the new `references/{human_proteome.fasta,vdjdb/<release>,imgt_germlines}` + `~/.mhcflurry/.download_done` sentinel paths, but the production `indices/{hisat2,star}/` path is exercised post-merge (test config keeps `resources/test/hisat2_index/` per the Issue spec). Bot review surfaced 3 items; applied 2 with [a21df77](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/469/commits/a21df77), deferred 1 with reasoning.

#### `resources/` blanket-ignore prevents the `!resources/test/` exception from ever firing

The CLAUDE.md table I wrote claimed `resources/test/` is "partial (test fixtures may be committed)" and the `.gitignore` comment block said "small committed test fixtures only (resources/test/)" — but the actual rule below was `resources/` (blanket, trailing-slash). Bot caught the inconsistency on first review pass. The non-obvious bit was the fix: my first attempt was `resources/` + `!resources/test/` + chr22 re-ignores below, and `git check-ignore -v` reported the blanket `resources/` still winning even for hypothetical new fixture paths. Verified the gitignore quirk in a scratch repo (`/tmp/test_gi/`): per git's docs, *"any matching file excluded by a previous pattern will become included again. It is not possible to re-include a file if a parent directory of that file is excluded. Git doesn't list excluded directories for performance reasons."* The trailing-slash form excludes the **directory** wholesale, and git short-circuits descent — the negation on a later line can't reach anything inside. Switch to `resources/*` (excludes only contents, leaves directory traversable) + `!resources/test/` (re-includes the dir), then chr22 re-ignores work. Documented in the .gitignore comment inline. **General rule: when a `!subdir/` negation isn't taking effect, check whether the parent rule is `parent/` (blanket, blocks descent) vs `parent/*` (contents-only, descends).** The two look almost identical and behave very differently.

#### `models_cache_dir` config key oversells — chose comment-tighten over rename

Bot's 3rd review item: `config/config.yaml` introduces `mhcflurry.models_cache_dir` and the comment suggests "Override to a project-local path (e.g. `indices/mhcflurry`) to share across clones," implying the actual model cache will move there. In current code ([mhc_affinity.smk:12-15](workflow/rules/mhc_affinity.smk)) the key only controls where the **sentinel** (`.download_done`) lives — the model cache itself stays in MHCflurry's platformdirs default (`~/.local/share/mhcflurry/`) until `MHCFLURRY_DATA_DIR` env-var wiring lands (listed as a follow-up in the PR body). Bot's two options: (a) rename to `sentinel_dir`, (b) tighten the comment. Chose (b) — renaming would force a **second rename** once the follow-up wires up the actual cache co-location, at which point `models_cache_dir` becomes accurate. Two renames for an interim state is worse than one honest comment. Comment now explicitly says "only the sentinel is controlled by this key" with the platformdirs default path called out.

#### Skipped the sentinel-naming-consistency nit on a VM-state cost argument

Bot's 2nd item (non-blocker per bot): VDJdb/IMGT sentinels use `.download.done`; MHCflurry uses `.download_done` (period vs underscore in the middle). Pure style nit. Skipping wasn't a default-deference move — renaming would force VMs with existing `~/.mhcflurry/.download_done` to re-run `mhcflurry-downloads fetch` (the sentinel becomes orphaned; the rule's sentinel target wouldn't match). For VDJdb/IMGT the sentinels are under `references/<release>/` which is gitignored and per-VM, so renaming there is cheap, but unifying just one direction still leaves the inconsistency. The cost/benefit favored skip; if anyone cares enough later, file a follow-up Issue. **General rule worth tagging: cosmetic-naming changes to sentinels with persistent on-disk state have a hidden cost — they invalidate the existing state on every machine that has it. Audit that before agreeing.**

#### Closure-ritual gate caught the "First post-merge VM run" aspirational checkbox

The PR Test plan included `- [ ] First post-merge VM run smoke-tests the migration path (will be tracked in a lab notebook entry on the next pipeline run)`. Closure-ritual gate refuses any `- [ ]`. The line was honest about the deferral but the gate doesn't parse inline "(will be tracked...)" deferral text — convention is tick-with-link or remove. No carrier Issue exists for this, so removed the line + added a blockquote below the Test plan noting "the migration path on the production VMs is verified by the first post-merge pipeline run; result will be captured in the next lab notebook entry." The blockquote captures the intent without leaving an unticked box that floats forever. **Pattern: post-merge verification items don't belong on the merging PR's Test plan — they belong as a follow-up Issue OR as a documented expectation outside the checklist.**

---

## 2026-05-22

### 14:15 UTC — Editor: Developer

**Headline:** [PR #467](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/467) closes [Issue #461](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/461) — drops `from __future__ import annotations` from 6 workflow scripts that were latently broken under Snakemake's `script:` wrapper since [PR #428](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/428) (2 days ago). Confirmed end-to-end with a fresh chr22 full-pipeline run (17/17 rules, exit 0, `report.html` produced). Adds a parametrised regression pytest that simulates Snakemake's preamble injection on every `workflow/scripts/*.py`, catching the bug class CI's `snakemake -n` dry-run is structurally blind to.

#### XS estimate was wrong on the first read; empirical verification mattered

Initial sizing was "30-60 min, drop the future imports, done." First file read killed that estimate: every one of the 6 scripts uses PEP 604 (`X | None`) and generic builtins (`list[X]`) extensively — ~70+ annotation sites total. Dropping the future import naïvely would have demanded an `Optional[X]` / `List[X]` rewrite across all 6 files, pushing this to M. Then `git log -p caa0387` revealed PR #428 added the future imports DELIBERATELY to support a lazy-import optimization (Issue #425) — removing them would also break the optimization. So the fix space looked like: (a) bulk-rewrite to typing-module forms, (b) revert PR #428 + lose the optimization, (c) patch Snakemake upstream, or (d) verify the breakage is real first. Picked (d). Two facts then rescued the scope: project's minimum Python is 3.11 (PEP 604 + generic builtins both native, no rewrite needed), and an AST walk found exactly 24 lazy-module annotations (`pd.DataFrame`) needing stringification — zero `Bio.X` annotations. Bulk-rewrite avoided, lazy-import optimization preserved. **Lesson: when sizing a "drop the deprecated thing" fix, first check (1) the project's minimum target version of the underlying feature, (2) whether the "deprecated thing" was added defensively or load-bearing. Both checks were missed on first sizing.**

#### `pipeline-pytest` now closes a structural blind-spot CI had

The regression test [`workflow/tests/test_no_future_imports_in_scripts.py`](workflow/tests/test_no_future_imports_in_scripts.py) parametrises over every `workflow/scripts/*.py`, constructs a Snakemake-style preamble + source, and `py_compile.compile`s the result. 15 scripts × <1 ms each, 0.05s total. Catches the exact bug class that hit PR #457 and would have hit any future PR adding `from __future__` to a `script:`-invoked file. This is the same shape as the "before merging a new rule, run integration on chr22" rule from PR #457's lab notebook entry — except now CI can do it for free instead of relying on a human-driven local run. Both gates are needed: this regression test catches `__future__`-class breakage; the chr22 integration run catches conda-solver, CLI-arg, and NaN-flow bugs that pure static analysis cannot.

#### Full chr22 pipeline as the final empirical check

Before opening the PR, ran `snakemake --cores 4 --use-conda --configfile config/test_config.yaml` end-to-end against the merged-state branch. 17/17 rules, ~4 min wall-clock (cached conda envs from prior runs), `report.html` + `report.tsv` + `mhc_presentation.tsv` + `peptides.tsv` + `novel_junctions.tsv` all produced. 5 of 6 previously-broken scripts ran successfully (`filter_junctions`, `assemble_contigs`, `run_mhcflurry`, `aggregate_filtering_stats`, `generate_report`); the 6th (`run_tcrdock`) didn't execute because `tcrdock.enabled=false` in test config, but its parse-check passes via the regression test. This is the first end-to-end chr22 run since [PR #428](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/428) merged 2 days ago — explaining why the breakage stayed latent (intermediate runs only targeted `panel.tsv`, never traversing the broken downstream rules).

---

### 13:00 UTC — Editor: Developer

**Headline:** [PR #464](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/464) ships the new `docs/features/` convention + first instance — a 9-slide Quarto reveal.js deck for the [Issue #204](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/204) VDJdb panel rule. Bot review round 1 flagged 3 findings (PR# in docstring, DMF5 framing, mermaid clickability); 2 fixed + 1 confirmed intentional ([commit `63837aa`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/464/commits/63837aa), [reply](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/464#issuecomment-4518537556)). User then asked for real DAG over the conceptual mermaid + bottom-margin overflow fixes on slides 3–5 — landed across 5 commits ([6ef46c3](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/464/commits/6ef46c3) through [dd43604](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/464/commits/dd43604)). Bot re-review came back clean ([comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/464#issuecomment-4518745348)) but only covered the first fixup, not the later commits — flagged but skipped a 3rd pass on user-validated visual polish.

#### `docs/features/` carves out a gap between `research/slides/` and `research/experiments/`

The motivation came from a user-flagged absence: substantive feature PRs (like [#457](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/457) closing [#204](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/204)) leave the codebase richer but produce no compact visual artifact summarizing **what shipped + how it performs on the chr22 test config**. The two sibling conventions don't fit: `research/slides/` is for analysis Issues with a decision-rule outcome; `research/experiments/` is for per-Issue analytical notebooks. Same Quarto reveal.js tooling, but the artifact is **pipeline documentation** — hence under `docs/`, not `research/`. Audience: team reference + portfolio peek. The new convention is documented in [`docs/features/README.md`](docs/features/README.md). Tagging this as a generalizable pattern: **a feature-doc deck is appropriate when a substantive Snakemake rule + supporting infra lands and would benefit from a frozen visual artifact (chr22 outputs + a few figures) alongside the per-patient notebook and manuscript text**. Not every PR — closure tasks, doc updates, single-fix PRs don't qualify (per the README).

#### Real-DAG-over-mermaid is the right default for feature-doc decks

The first iteration of slide 3 ("Where it fits") used a hand-drawn mermaid diagram showing `OptiType → fetch_vdjdb_panel → Issue #205/#206`. User pushed back: "can we maybe integrate the content onto the real DAG?" The right answer turned out to be **strictly more informative**: render the actual snakevision rule-graph (chr22 test config) via `scripts/visualize_dag.sh`, freeze a copy at [`docs/features/issue_204_vdjdb_panel/figures/dag.svg`](docs/features/issue_204_vdjdb_panel/figures/dag.svg), and put the conceptual annotations (inputs / outputs / deferred sub-issues) in a side-by-side sidebar. Two non-obvious gotchas surfaced:

- **The default `rule all` does not pull `fetch_vdjdb_panel` into the DAG**, because the rule is gated on `hla.enabled` AND is not yet a transitive dep of `report.html` ([Issue #206](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/206) wires it). The bare `bash scripts/visualize_dag.sh --clean` produces a DAG *without* the new rule — the rule-graph reflects the consumer graph, not the rule registry. **Fix:** pass an explicit panel target after a `--` terminator (`snakemake --configfile ... --rulegraph --forceall -- <report.html> <panel.tsv>`). The `--` is required to escape argparse's `nargs="+"` swallowing the target as a configfile (the existing CLAUDE.md gotcha). Documented the regen command in [`figures/_regenerate_figures.py`](docs/features/issue_204_vdjdb_panel/figures/_regenerate_figures.py) docstring.
- **Visual highlight via hand-edit > Python post-processor here.** Decided against a regen pipeline that re-applies the orange-ring styling after each snakevision run — the SVG won't be regenerated often (only when pipeline rule structure changes, e.g. after #206 lands), and the hand-edit is recipe-documented in the same docstring. Trade-off favors simplicity at the feature-doc tier; would not generalize this to production output where regeneration must be byte-deterministic.

#### Reveal.js silently clips content past 720px — `quarto render` exit-zero is not a slide-fit check

Slides 3 and 4 both rendered cleanly with no warnings and a healthy `Output created: slides.html` — but the user pointed out that the bottom rows of slide 4's filter table and the lower sidebar bullets on slide 3 were invisible to a viewer. Reveal.js renders to a fixed canvas (1280×720 in this deck) and **silently clips** content past the bottom edge. No scrollbar, no indicator, nothing in the build log. The rendered HTML alone gives the same impression of completion regardless of whether content fits. This is a structural blind-spot of the existing workflow: `quarto render` validates Markdown/YAML syntax, not slide-fit. Captured as a memory rule in [`feedback_slide_overflow.md`](file:///Users/jin-holee/.claude/projects/-Users-jin-holee-dev-GitHub-Jin-HoMLee-splice-neoepitope-pipeline/memory/feedback_slide_overflow.md): after every `quarto render`, visually scan each slide; treat the bottom ~10% as breathing room; prefer `font-size` tightening + content compression over `.scrollable` (which hides the problem). Affected slides 3 + 4 + 5 fixed in [c6834ea](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/464/commits/c6834ea) (initial collapse) and [dd43604](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/464/commits/dd43604) (restore bullets + slide 5 resize per user feedback). The collapse-to-paragraphs first try was wrong — user correctly redirected that "really hard to read" comes before fits-on-screen; the right answer is `font-size: 0.62em` + `line-height: 1.25` with the bullet structure preserved.

#### Bot re-review only diffs against its previous-review SHA, not HEAD

When triggered for round 2, the `@-claude review` action reviewed only the post-fixup-`63837aa` state — not the 4 subsequent commits (DAG swap, node highlight, slide overflow). Its review body literally said "Re-Review (post-fixup `63837aa`)" and cited line numbers from that state. The 4 newer commits (visual polish on `figures/dag.svg`, sidebar layout, slide 5 figure size) **did not appear in its scope**. Hypothesis: the Action SHA-pins to its previous review's base and reports the delta since, rather than diffing the full PR against `main`. Not investigated further — for this PR the unreviewed commits were user-validated visual layout, where bot value-add is minimal anyway. **Practical takeaway: do not assume "second `@-claude review` = full re-review of latest HEAD".** If a PR has accumulated commits past the first review's last-reviewed SHA and you actually need bot coverage on them, expect to either close+reopen the PR or live with the gap. Cheaper rule: don't batch commits between bot reviews if you want them all covered — request review on substantive cycles, not at the end.

---

### 09:20 UTC — Editor: Developer

**Headline:** [PR #457](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/457) shipped closing [Issue #204](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/204) — `fetch_vdjdb_panel` Snakemake rule producing per-allele top-10 paired α/β TCR panel with full chain sequences via stitchr. Bot review (`@-claude review` round 1) flagged 4 issues + 3 minor obs; addressed 3 as-is + 1 with reshaped fix ([commit `1941b48`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/457/commits/1941b48), [reply](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/457#issuecomment-4517153084)). Then the local chr22 end-to-end run on sample SRR9143066 surfaced **4 substantive bugs that unit tests + dry-runs + bot review had all missed** ([commit `cfe1bc4`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/457/commits/cfe1bc4)) — closed in one cluster. Final panel: 20 TCRs across 2 alleles (HLA-A\*02:01, HLA-B\*07:02 both `ok`), HLA-C\*07:02 → `empty` (biological zero, 0 VDJdb entries).

#### Local end-to-end runs catch what dry-runs structurally cannot

The 4 bugs surfaced by the chr22 run — `IMGTgeneDL>=0.7.0` pin (PyPI max 0.6.1), missing `from __future__` shebang, missing paired α/β `.notna()` filter, stitchr `-species HUMAN` CLI typo (should be `-s` or omitted) — share a shape: each requires **actually executing the rule body**. CI's `pipeline-snakemake-dry-run` walks the DAG without building conda envs or running scripts; `pipeline-pytest` exercises unit logic with curated fixtures (paired α/β rows only, no stitchr subprocess). Both stayed green throughout. The bot reviewed source code but didn't execute it either. None of these layers would have caught any of the 4. The pattern is general: **anywhere a pipeline rule shells out to a subprocess with CLI args, or depends on the actual conda-env solver result, the integration path is the only honest gate**. Filing this as a memory rule for next time — "before merging a new rule that invokes subprocesses or external CLI tools, run it on the chr22 test config and inspect the output, not just the dry-run".

#### `from __future__ import annotations` is a Snakemake `script:` foot-gun

The 3rd bug — `SyntaxError: from __future__ imports must occur at the beginning of the file` — traces to Snakemake's `PythonScript.write_script` (`snakemake/script/__init__.py:807`), which unconditionally prepends its `snakemake = pickle.loads(...)` preamble before the source. There's a TODO at line 44 (`PY_PREAMBLE_RE = re.compile(...)` with a comment "use this to find the right place for inserting the preamble") but the regex was added without ever being wired up. So **any** workflow script using `from __future__` will SyntaxError under the `script:` directive. 6 other scripts in [workflow/scripts/](workflow/scripts/) (`aggregate_filtering_stats.py`, `assemble_contigs.py`, `filter_junctions.py`, `generate_report.py`, `run_mhcflurry.py`, `run_tcrdock.py`) all share this layout — they're either latently broken or there's a Snakemake-version subtlety I missed. Filed [Issue #461](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/461) to audit + reproduce. The fix for the one script in PR #457 was dropping the future import (the file used `Optional[str]` from `typing`, never PEP 604 `X | None` or generic builtins, so the import was cosmetic).

#### Bot review #1 spent its budget on style; the real bugs needed execution

Comparing the bot review's 4 substantive issues against the 4 bugs the local run actually caught: zero overlap. The bot flagged `lambda wc:` unused wrappers, an unclosed file handle in a test, a wasteful subprocess pattern (alpha + beta stitched before checking alpha's None), and IMGT cache discovery silent-skip. None of those would prevent the panel from being produced. Meanwhile the panel WAS empty in the first 4 runs because of the bugs the bot missed. This isn't a knock on the bot — static review is fundamentally code-shape review, and the bugs the chr22 run caught are dynamic (conda solver, CLI args, NaN values flowing through pandas). Next time: **don't treat bot-review-clean as merge-ready for a new rule**. The Anthropic Agent SDK billing split that lands 2026-06-15 raises the cost of speculative review rounds; better to spend the bot budget on already-executed code, not source-code-shape passes.

#### Process note — exit code 0 from `... 2>&1 | tee log` masked the first failure

The first chr22 run reported `exit code 0` via the background-task notification, but Snakemake itself had errored. `tee`'s exit code shadowed snakemake's non-zero exit. Lesson: when piping snakemake output for capture, either drop `tee` or `set -o pipefail` before the pipeline. Not worth a memory rule; just a sharp mental tag.

---

## 2026-05-21

### 17:50 UTC — Editor: Developer

**Headline:** [PR #449](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/449) shipped, closing [Issue #447](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/447) (document Python execution contexts; fix stale `.venv/bin/python` plan-file refs). In-tree: [CLAUDE.md](CLAUDE.md) gains a `## Python environments` section (4-env matrix + no-root-`.venv` invariant) and 4 stale refs in [docs/superpowers/plans/2026-05-19-issue-17-star-aligner.md](docs/superpowers/plans/2026-05-19-issue-17-star-aligner.md) corrected. Out-of-tree (in `.claude/memory/`): new `shared/feedback_python_environments.md` Always-in-effect rule, `shared/MEMORY.md` index entry, Developer `MEMORY.md` inline pointer + line 30 stale-ref fix, `feedback_test_before_pr.md` body fix. Bot review (`@claude review`) = **Approve** with two non-blocking observations folded into follow-ups, not this PR.

#### "5 envs" collapsed to "4 envs" under critical review — `python -c` is not a separate env

Initial surfacing brief framed this as 5 functional Python envs in the project (per-rule conda, snakemake conda, `workflow/tests/.venv`, `research/.venv`, ad-hoc `python -c`). On critical re-read before drafting CLAUDE.md, the 5th case wasn't actually distinct — a bare `python -c "..."` invocation just uses whatever interpreter the active env exposes, which in practice means the `snakemake` conda env after `conda activate snakemake`. Calling it a separate env was a count-inflation artifact of the original "I keep getting confused about which Python" symptom (every interpreter touch felt like its own decision). Collapsing to 4 envs in the docs was both honest and more useful: the matrix gets a cleaner "use case → canonical path" mapping, and the `python -c` row in the quick-lookup table folds into the `conda activate snakemake` row. The recursion-on-pitch happened **before** any commits — worth recording because the impulse to ship the original framing was strong; the rewrite was 3 minutes of work that prevented stale doc debt.

#### Out-of-tree memory changes don't show in `git diff`, but ARE the bulk of the work

The PR's `git diff main --stat` shows 2 files changed (24 insertions, 4 deletions). The actual session output was 6 files modified: 4 of them inside `.claude/memory/` (gitignored). For future-me looking at the merged commit and wondering "where did the shared rule go?" — it's out-of-tree by design (each role's memory is a per-clone artifact, not source-tree state). The PR + the Issue body + this lab notebook entry are the durable record that ties in-tree to out-of-tree. The shape worked here because Issue #447's AC list explicitly enumerated both halves, so the closure-ritual gate would catch any miss.

#### Vdjdb plan fix (24 hits) deferred via Issue #204 carrier — the file lives on that branch

AC #6 said "Two plan files corrected." Only the star-aligner one is reachable from main (4 hits → 0). The vdjdb plan file ([docs/superpowers/plans/2026-05-20-vdjdb-panel-plan.md](docs/superpowers/plans/2026-05-20-vdjdb-panel-plan.md), 24 hits) lives exclusively on the `204-...` branch — fixing it here would need a cross-branch merge or a rebase dance, and the file will land on main when [Issue #204](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/204) work merges anyway. Solution: defer-tick AC #6 with a link to a heads-up comment filed on Issue #204 ([comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/204#issuecomment-4509276737)) that explicitly enumerates the canonical replacements (pytest invocations → `workflow/tests/.venv/bin/python`; `python -c` quick checks → bare `python` post `conda activate snakemake`). This is the closure-ritual deferral pattern in action — tick `- [x]` with link to carrier Issue, mechanism passes, work doesn't get forgotten. The grep AC (`grep -rn '\.venv/bin/python' docs/ scripts/ workflow/ .github/ | grep -v 'workflow/tests/\.venv'`) returns empty on the 447 branch because the vdjdb file isn't on main yet — so the AC is honestly green at merge-time even with the deferral.

#### Memory broadcast skipped — shared/MEMORY.md absorption is enough

[team_memory_broadcasts.md](.claude/memory/shared/team_memory_broadcasts.md) is the channel for cross-role rule changes that need explicit attention before next session start. For this Python-envs rule, the change is already indexed in `shared/MEMORY.md` as an Always-in-effect entry — Scientist + PM will absorb it automatically at next `/memory` without a broadcast nudge. Broadcasting would have been noise. AC #4 ticked with the rationale inline rather than removed entirely, so the closed Issue carries an audit trail of the deliberate skip (vs. an oversight). [Issue #451](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/451) filed as the open follow-up — would consolidate `workflow/tests/.venv` + `research/.venv` into a single repo-root `.venv`, making `.venv/bin/python` actually correct and dissolving the whole 4-env distinction; gated on Sci sign-off because it crosses notebook-deps ownership.

---

## 2026-05-20

### 16:01 UTC — Editor: Developer

**Headline:** [PR #428](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/428) shipped, closing [Issue #425](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/425) (lazy-import pandas in 6 pandas-heavy scripts). Module-scope `import pandas as pd` (and `from Bio import SeqIO` in [generate_report.py](workflow/scripts/generate_report.py)) moved into the function bodies that actually call `pd.*`, across [generate_report.py](workflow/scripts/generate_report.py), [filter_junctions.py](workflow/scripts/filter_junctions.py), [assemble_contigs.py](workflow/scripts/assemble_contigs.py), [aggregate_filtering_stats.py](workflow/scripts/aggregate_filtering_stats.py), [run_mhcflurry.py](workflow/scripts/run_mhcflurry.py), [run_tcrdock.py](workflow/scripts/run_tcrdock.py). Cold-import per script dropped from ~0.28–0.39s to ~0.03–0.04s (86–90% reduction); full suite + integration runs in 1.42s wall on the M1 8 GB box. One commit per script (6 commits) for surgical review.

#### Key trick — `from __future__ import annotations` defers type hints to strings

The blocker for naive lazy-import was that every script had `def foo(df: pd.DataFrame) -> pd.DataFrame:` signatures, so `pd` had to exist at function-definition time. Adding `from __future__ import annotations` (PEP 563) at the top of each script defers ALL annotations to strings, so `pd.DataFrame` becomes a literal string never evaluated unless something explicitly inspects `__annotations__`. With that in place, the only places that genuinely need `import pandas as pd` are the function bodies that make direct `pd.*` calls (`pd.read_csv`, `pd.DataFrame(...)`, `pd.notna`, etc.). Functions that only *consume* DataFrames as parameters don't need the import — that distinction halved the number of in-function imports in most scripts.

#### Per-script cost varied wildly — generate_report.py needed 11 imports, others 1–4

The 6 scripts had similar cold-import times (~0.28–0.39s, all pandas-dominated) but very different surface areas for the refactor: `generate_report.py` makes pandas calls in 11 distinct functions (report assembly is pandas-heavy by design), while the alignment-pipeline scripts (`run_mhcflurry`, `run_tcrdock`, `aggregate_filtering_stats`) each needed only 1–4. No clean way to centralise — a `_get_pandas()` helper would just move the cost back to the first call site. Kept it as explicit `import pandas as pd` lines at the top of each pandas-using function so future readers see the pattern uniformly.

#### Tests caught a real bug — 5 missed `pd.` call-sites in generate_report.py

Initial pass on [generate_report.py](workflow/scripts/generate_report.py) added in-function imports for the 6 most obvious functions; pytest immediately surfaced `NameError: name 'pd' is not defined` from 5 other call-sites I hadn't seen on the grep sweep (small one-liner uses of `pd.notna` and `pd.DataFrame(...)`). Patched in the same commit before pushing — the test suite acting as a completeness oracle is exactly the reason to gate each commit on `pytest -q` before moving on. Not a methodology failure; the workflow worked as designed (commit-after-green per [feedback_multi_file_workflow.md](memory/feedback_multi_file_workflow.md)).

#### Scope decision — script-side only, test-side deferred

[Issue #425](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/425) explicitly carved test-file refactor out of scope. Rationale: the dominant cost was the script-side transitive load (the script imports pandas → test imports the script → pytest pays pandas cost during collection). With scripts now lazy-loading, the test files' own module-scope pandas imports pay only their direct 0.28–0.62s once-per-file, which is fine when memory isn't tight. Refactoring test files into pytest fixtures would be churn for marginal benefit. The original [Issue #364](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/364) symptom (5-minute hangs under deliberate memory pressure) is now structurally impossible because no script holds pandas at module scope — the page-cache-thrash amplification surface is gone.

---

### 15:30 UTC — Editor: Developer

**Headline:** [PR #424](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/424) shipped, closing [Issue #364](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/364) (local pytest collection hangs under memory pressure). Two deliverables: a committed per-file cold-import profiler at [workflow/tests/scripts/profile_imports.py](workflow/tests/scripts/profile_imports.py) (AC 1) + a "Running on memory-tight machines" section added to [workflow/tests/README.md](workflow/tests/README.md) (AC 2). Remaining AC 3 (verify suite < 60s on M1 8 GB under deliberate paging) carved out as [Issue #425](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/425) and absorbed alongside the lazy-pandas refactor work.

#### 10-min quick-win triage → diagnostic-only investigation

User opened with "any quick wins for the next 10 min?" — I scanned my open dev queue, no PRs blocked on review, and offered three options ([Issue #345](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/345) naming refactor, [Issue #364](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/364) import profile, memory hygiene). Picked #364 as diagnostic-only — concrete signal output, no merge pressure, directly unblocks future local pytest work. The "diagnostic-only" framing kept the budget bounded: the goal was a comment with findings, not a PR.

#### Hypothesis refinement — the original Issue was overstated

Issue body hypothesised that test files import heavy ML libraries (`mhcflurry`, `torch`, `tensorflow`) at module scope. Profiler showed that's **not the case** — none of those are at module scope anywhere. The dominant module-level import is `pandas` (0.26s cold), pulled in by 6 test files + 6 project scripts. Under memory pressure the OS page cache thrashes, amplifying every fresh `pandas` import into a multi-minute event (matches the original `import pandas` hung > 2.5 min observation). Pytest's assertion rewriter compounds it. So the failure mode is **paging-amplified, not import-count-amplified** — refactor target should be lazy-importing pandas in the heavy scripts, which I carved into [Issue #425](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/425).

#### Scope escalation diagnostic → PR

After posting findings as a comment, user asked "are we closing without changes?" — caught me leaving 3 ACs unticked. Two paths offered: fan out to sub-issues, OR ship a small PR covering ACs 1+2 now. User picked the PR route. Profiler script promoted from `/tmp/` to `workflow/tests/scripts/profile_imports.py`, README updated with the two known workarounds (single-file runs, `-p no:assertion`).

#### Workflow correction — close Issues with PR, carve follow-up

Initially left #364 open with "stays open until #425 ships" for the M1-verify AC. User flagged: prefer carving the remaining work into a new Issue and closing the original with the PR. Saved as a new shared memory ([feedback_close_issue_with_pr.md](.claude/memory/shared/feedback_close_issue_with_pr.md)) and applied retroactively — #425 absorbed AC 3, PR body changed from `Refs Issue #364` to `Closes Issue #364`, the three ACs on #364 ticked with deferred-to links. The closure-ritual gate (`scripts/audit_and_merge.sh`) will see all boxes ticked at merge time.

#### Bot review — two nice-to-haves applied, two nits skipped per bot's own non-blocker calls

Bot finished in 2m 4s with "Good to merge as-is" + three observations (commit `79282cb` covers two):
- **Applied:** docstring note in `truncate_to_imports_and_defs` flagging that `ast.Try` ends truncation early (no module-scope try-blocks today, but future-fragile).
- **Applied:** surface last non-empty stderr line on subprocess failure as `  └─ [file] exit=N: <error>`. Naive `print(r.stderr)` dumped Python's full traceback including the synthetic exec'd source — truncating to the last line gives the actionable signal.
- **Skipped:** reversed `'__main__' == __name__` form (no file uses it) + README "while" phrasing nit (bot self-identified as non-blocker).

#### Methodology note — profiler design choice

Tried `python -X importtime -m pytest --collect-only` first; cumulative self-time summed to 0.21s but wall was 20.66s. The 100× gap is pytest's assertion rewriter (bytecode rewriting per test file via a path that bypasses `-X importtime`). So importtime alone is misleading for diagnosing pytest collection cost — needed the cold-subprocess per-file approach (`profile_imports.py`) to isolate the actual import-time-amplifiable cost surface. Documented in the script's docstring so future-me doesn't fall into the same trap.

---

## 2026-05-19

### 18:00 UTC — Editor: Developer

**Headline:** [PR #418](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/418) (trap-based cleanup for `star_index` temp_annotation.gtf, closes [Issue #412](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/412)) shipped. Surgical 5-line diff in [workflow/rules/alignment.smk](workflow/rules/alignment.smk): replaces an end-of-shell-block manual cleanup conditional with a `trap 'rm -f resources/temp_annotation.gtf' EXIT` registered immediately at the top of the `if .gz` branch. Under `set -euo pipefail`, the previous form skipped cleanup whenever STAR aborted mid-run (e.g. OOM during index build on the hardware-constrained spot VM); the trap fires on any exit path.

#### Scope decision — 10-min triage pass, picked smallest concrete win

User had 10 minutes; offered four options (#412 trap, #345 raw_junctions.tsv rename, #352 PyTorch spike notes, backlog triage). User chose triage. Triage surfaced nothing genuinely obsolete across 33 open dev issues, but flagged 3 today-filed `no-ms` issues (#409, #411, #412) and ranked the next-up XS wins. User then picked #412 — surgical, zero blockers, single rule, AC list verbatim from the issue body.

#### Review iteration — bot review at [comment-4490492192](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/418#issuecomment-4490492192)

Bot finished in 1m 45s with LGTM and one minor observation: the original draft placed the trap *after* `gunzip`, which leaves a gap where a `gunzip` failure (corrupted GTF archive) trips `set -e` before the trap registers → partial temp file leaks. The trailing `if [[ -f ... ]]; then rm` form had the same gap, so it was pre-existing rather than a regression, and bot called it "worth considering as a follow-up — not a blocker."

User asked for my recommendation. Picked **apply now**: the move is 1 line, the deviation from the issue's literal AC text ("after the gunzip line") was easy to fix by editing the issue body, and a follow-up issue + PR + review cycle for the same shell line would cost more than the fix itself. Commit `5c2c316`. Issue #412 AC text revised in-flight to specify "**before** the `gunzip` line" with a 2026-05-19 dated note pointing back to the bot review. Reply at [comment-4490566600](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/418#issuecomment-4490566600).

#### Cloud verification — deferred to Issue #378 again

Same pattern as [PR #410](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/410): the third AC ("verify on successful STAR index build that temp file is cleaned up") needs production VM RAM (>8 GB for STAR index build, exceeds M1 limit per CLAUDE.md). [Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378) (patient_002 re-run) already exercises the STAR path end-to-end; trap behavior on the success path is equivalent to the prior manual-cleanup behavior by construction, so this is verification-of-equivalence rather than functional change.

#### Milestone gap flagged

Issue #412 (along with #409, #411) was filed today with no milestone. PR body flagged `i5 - S3 - STAR Polish & Aligner Verification` as the likely home (siblings: #345/#375/#377/#378), but milestone assignment left to PM.

### 13:56 UTC — Editor: Developer

**Headline:** [PR #410](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/410) (production aligner default flipped HISAT2 → STAR, closes [Issue #17](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/17)) shipped. Minimum-viable scope: the two STAR flags that materially affect novel-junction sensitivity (`--twopassMode Basic`, `--limitSjdbInsertNsj 2000000`) plus removal of three silent-throttle filters (`--outSJfilterReads Unique`, `--outSJfilterCountUniqueMin`, `--outSJfilterCountTotalMin`). Paper-fidelity sensitivity tuning explicitly deferred to [Issue #411](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/411) with per-flag benchmarks — avoids cargo-culting the Nature paper command wholesale without measured deltas on our data.

#### Scope decision — middle ground vs full paper-fidelity

Issue body's plan was a one-line config flip; a 2026-04-21 implementation comment from Jin-Ho added ~14 STAR parameter changes from the Nature paper recipe. Three-option analysis (full migration / middle ground / config-flip-only) → went middle. Rationale: the count-filter removals + `--twopassMode Basic` are **correctness** changes (multi-mappers should count; 2-pass is what gives STAR its novel-junction edge over HISAT2), but the sensitivity flags (`--outFilterMatchNminOverLread 0.33`, `--alignIntronMax 500000`, etc.) are **tunables** that should be justified per-flag with a benchmark, not copy-pasted. Cargo-culting a recipe across pipelines without measured deltas is exactly the failure mode the portfolio-lens memory warns against.

#### Implementation — 4 files

- [workflow/rules/alignment.smk](workflow/rules/alignment.smk) — 5-line diff in the `star_align` shell block (2 added, 3 removed)
- [config/config.yaml](config/config.yaml) — `aligner: hisat2` → `star`
- [config/test_config.yaml](config/test_config.yaml) — explicit `aligner: hisat2` pin added. Missing this would have inherited the new `star` default and broken local M1 chr22 dev (32 GB STAR index won't fit on 8 GB RAM). Caught during spec self-review, promoted from out-of-scope to in-scope mid-design.
- [workflow/tests/test_alignment_star_command.py](workflow/tests/test_alignment_star_command.py) — new pytest snapshot test, **first shell-level rule snapshot test in the codebase**. Renders the `star_align` command via `snakemake --dry-run --printshellcmds` with an `aligner=star` config override + stub inputs in `tmp_path`, asserts substring presence/absence of the 5 target flags. 5 tests run in 0.83s (module-scoped fixture; see review iteration below).

#### Scope creep — pytest venv setup documented

The developer memory's `feedback_test_before_pr.md` references `.venv/bin/python -m pytest` but post-2026-05-14 separate-clones migration meant this clone had no `.venv` set up. Added a [workflow/tests/README.md](workflow/tests/README.md) documenting the canonical pattern: `pyenv local 3.13.5` + `python -m venv workflow/tests/.venv` + install from [workflow/tests/requirements-test.txt](workflow/tests/requirements-test.txt). Tests that invoke `snakemake` via subprocess need the `snakemake` conda env activated for the binary on PATH — README documents both run modes. Test runner (`.venv`) and workflow runner (`snakemake` conda env) deliberately split: keeps the workhorse env pristine.

#### dag.pdf cleanup — in-PR

`scripts/visualize_dag.sh` outputs `dag.svg` (gitignored), not `dag.pdf`. The tracked `dag.pdf` was vestigial from `f4f3858` (months ago, before the SVG script existed); `.gitignore` already listed it at line 238 but the rule only takes effect once untracked. Removed in-PR rather than as a follow-up — kept this scope but saved the follow-up overhead.

#### Review iteration — bot review at [comment-4488238391](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/410#issuecomment-4488238391)

Bot finished in 3m 29s. 2 actionable items + 1 semantic note + 1 out-of-scope pre-existing issue:

- `81e858a` — fixture scope-bump from `function` → `module` + swap `tmp_path` for `tmp_path_factory`. 5 tests now share 1 `snakemake -n` subprocess; wall time 2.53s → 0.83s. Bot's perf nit, accepted as a clean win.
- `f25d2c4` — README pyenv version `3.14.4` → `3.13.5`. Bot flagged that 3.14.4 may not be installable on a fresh clone (despite being a stable release since 2025-10) — 3.13.5 is in every pyenv install list since 2024-10. Comment now also notes any stable Python ≥3.10 works; the pin is illustrative, not mandatory.
- **Semantic note:** bot caught that the commit/spec wording "throttles output below paper baseline" was imprecise for the count-filter removals. The removed `--outSJfilterCountUniqueMin 1 1 1 1` flags were *more permissive* than STAR's `3 1 1 1` defaults for non-canonical junctions (position 1); removing them reverts to **stricter** non-canonical thresholds. Combined effect on non-canonical bucket is a mix (stricter unique threshold vs newly-contributing multi-mappers vs more 2-pass discoveries). Net direction on canonical junctions is still more output; non-canonical is non-trivial and queued for [Issue #411](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/411) benchmark. Logged in [comment-4488440732](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/410#issuecomment-4488440732) so future-grep doesn't get misled by the imprecise commit wording.
- **Pre-existing out-of-scope:** `star_index` rule's `resources/temp_annotation.gtf` lacks a trap-based cleanup — `set -euo pipefail` skips the manual cleanup on STAR abort. Filed as [Issue #412](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/412) (P3 refactor).

#### Cloud verification — deferred to Issue #378 via 3 new ACs

End-to-end STAR verification requires production VM RAM that doesn't exist locally. Appended 3 ACs to [Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378) (already-queued patient_002 re-run for the [PR #402](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/402) strand fix) — single cloud run now verifies both PR #402's strand rescue AND PR #410's new STAR defaults. Saves duplicate cloud cost vs running a separate single-sample smoke test as part of this PR. ACs verify: (1) `SJ.out.tab` contains multi-mapper-only rows (`multi_reads > 0 ∧ unique_reads == 0`); (2) 2-pass log line present; (3) no `limitSjdbInsertNsj exceeded` fatal.

### 11:42 UTC — Editor: Developer

**Headline:** Dropped the `workspace/<role>` camping-branch convention from shared memory. Rule's stated justification ("`main` is a shared ref; one worktree's checkout blocks the others") is dead since the 2026-05-14 migration from git worktrees → three independent clones. Each clone has its own `.git`; cross-clone HEAD contention is impossible by construction.

**Why now:** User noted at session start that workspace branches were a workaround for a constraint that no longer exists. Analysis confirmed: no remaining technical reason; the weaker "identify which role parked here" argument was already covered by the clone path itself.

**Changes — 5 memory files:**

- `memory/shared/MEMORY.md` — removed Always-in-effect rule + Reference link
- `memory/shared/feedback_workspace_branch.md` — deleted entirely
- `memory/shared/feedback_team_structure.md` — "Git worktree layout" section rewrote as "Clone layout", drops `git worktree add ... -b workspace/<role>` setup commands, notes the 2026-05-14 migration
- `memory/shared/feedback_lab_notebook.md` — dropped "rather than from `workspace/<role>`" sub-clause in the multi-session-parent guidance (now just "cut from latest `main`")
- `memory/shared/feedback_no_cd.md` — tightened "worktree" → "clone" wording; cwd-persistence rule itself unchanged

**Broadcast decision:** user declined the cross-role broadcast — removals don't require explicit attention; other roles will absorb the cleanup at next `/memory` re-read since the rule simply no longer appears. The broadcast convention is most valuable when rules are *added* (action required) rather than *removed* (action stops).

---

## 2026-05-18

### 21:59 UTC — Editor: Developer

**Headline:** [PR #402](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/402) (STAR strand=0 motif rescue, closes [Issue #374](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/374)) review-response shipped — second silent-contamination axis on the splice-junction extraction surface, sibling to the HISAT2 BED12 anchor-outer bug fixed in [PR #372](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/372). Inline awk in `alignment.smk` was emitting STAR `SJ.out.tab` records with strand `.` whenever STAR couldn't infer strand directly (col 4 = 0); `bedtools getfasta` (called without `-s` in [`assemble_contigs.py`](workflow/scripts/assemble_contigs.py)) then treated `.` as forward orientation, yielding reverse-orientation flanking sequence for true minus-strand junctions and wrong-frame translation downstream.

**Implementation:** new [`workflow/scripts/star_sj_to_junctions.py`](workflow/scripts/star_sj_to_junctions.py) mirrors the structure of [`bed12_to_junctions.py`](workflow/scripts/bed12_to_junctions.py) introduced in [PR #372](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/372). Uses STAR's col 4 directly when ∈ {1,2}; rescues strand from col 5 (intron motif) when col 4 = 0, with the 6 canonical/semi-canonical motifs mapping to ± and motif=0 (truly non-canonical) dropped rather than emitted as `.`. Drop-vs-`.` policy is the explicit lesson from PR #372: contaminating the candidate set is worse than under-recalling. Output format unchanged (`<chrom>:<donor>:<end>:<strand>\t<reads>`), so DAG stays identical and downstream rules are aligner-agnostic.

**Review iteration:** bot review at [comment-4479665786](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/402#issuecomment-4479665786) surfaced 2 actionable items + 3 verified-correct observations. Both items applied:

- `e14b560` — added 4-case parametric test pinning the "col 4 takes priority over motif" invariant (`strand=1 motif=2 → +`, `strand=2 motif=1 → -`, plus the non-canonical-motif counterparts). Reviewer correctly identified that the priority rule was the core of the rescue logic and had no test coverage — `TestDirectStrand` only exercised the agreeing pairs.
- `d3f02d1` — made the `strand_code == 0` branch in `_resolve_strand` explicit; any other code now returns `None` (drop) rather than falling through to motif lookup. No behavior change for STAR-spec inputs (col 4 ∈ {0,1,2}), just tightens intent. Borderline against "don't validate scenarios that can't happen" — applied for clarity, not defense.

Reply at [comment-4482564368](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/402#issuecomment-4482564368).

**Verification trail:**

- Full pytest suite **282 passed, 25 skipped** on `d3f02d1` (was 278 pre-review; +4 new priority cases, no regressions).
- All 3 CI checks (`pipeline-pytest`, `pipeline-snakemake-dry-run`, `ci-tools-pytest`) green on the post-review SHA.
- DAG sanity: `snakemake --rulegraph` output on `origin/main` vs HEAD is identical when sorted (only edge-emission order differs, which is non-deterministic in snakemake's dot output). Required running `bash scripts/prepare_test_data.sh` first to materialize `resources/test/chr22.fa` — `--rulegraph` still does an input-existence check at DAG-build time.

**Process notes:**

- DAG render attempt initially failed with `MissingInputException` because this clone had no `resources/test/`. `--clean` mode in [`visualize_dag.sh`](scripts/visualize_dag.sh) generates placeholders for sample FASTQs only, not for reference files like `chr22.fa` — the symlink workspace inherits the missing parent. Followed user's call to download via `prepare_test_data.sh` rather than tick the box on inspection alone; the airtight sort-diff result was worth the ~15 min download.
- Spurious `git -c commit.gpgsign=true` slipped into the first commit attempt and failed because no GPG secret key is configured for `jinho.michael.lee@gmail.com`. Recent branch commits are unsigned (`git log --pretty='%G?' → N`), so signing is opt-out by user config — re-ran without the override. Not a memory-worthy slip; just a stray flag.
- Two review-response commits split per-concern (test first, then refactor) to match the existing branch's per-file commit pattern (`921fabf test`, `3dff7bb feat`, `160ba9e refactor`, `d7dffec docs`).

---

## 2026-05-17

### 19:31 UTC — Editor: Developer

**Headline:** [PR #389](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/389) (audit_and_merge.sh — closure-ritual gate) shipped, closing [Issue #357](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/357). Merged via the inaugural self-application of the script (the gate enforcing its own ship). Sister mechanism to the just-closed [Issue #325](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/325) (post-merge critic) — together they form a pre-merge gate + post-merge detective pair for the closure-ritual rule that broke 4× in 10 days. Closes 7th of 7 issues on the [`dev-i1` milestone](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/milestone/23) (due 2026-05-19).

**Implementation:** 105-line bash script (67 non-comment) at `scripts/audit_and_merge.sh`. Awk-based extractor for unticked `- [ ]` lines under named `## <Heading>` sections, called against PR body Test plan + every linked Issue's Acceptance criteria (via `gh pr view --json closingIssuesReferences`). Exits 1 with the gaps printed to stderr on any miss; otherwise forwards to `gh pr merge` with default `--squash --delete-branch` (override via positional args).

**Review iteration:** bot review at [comment-4471860832](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/389#issuecomment-4471860832) surfaced one real-but-impractical bug (`echo "$VAR"` swallows leading `-n`/`-e` flags) + 3 minor observations. Applied the `printf '%s\n'` fix at [22bf663](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/commit/22bf663) as defensive idiom (verified empirically that bash echo only flag-parses when the entire arg matches `-[neE]+`, so the trigger condition isn't reachable from real GitHub-rendered bodies — but the fix is canonical and zero-cost). Pushed back on the 3 minor items with technical reasoning (naming-mismatch is contextually correct under control flow; no `--yes` keeps the irreversible-action manual out; `&&` vs `if/then` is idiomatic bash with `set -e` exemption). Reply at [comment-4472262410](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/389#issuecomment-4472262410).

**Smoke-test trail:** inline awk unit tests for 4 behaviors (unticked extraction, count under section, no-heading → empty, multi-section isolation) all PASS in [/tmp/test_audit_awk.sh](file:///tmp/test_audit_awk.sh) (kept local — not added to `tools/ci/` because the awk lives inline in the bash script; CI surface lives at the `pipeline-pytest` boundary). Real-PR smoke against already-merged [PR #391](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/391) (today's news_log) confirmed the audit-pass path; `gh pr merge` returned "already merged" cleanly.

**Why this matters now:** the closure-ritual rule has the highest broken-count in recent memory ([Issue #280](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/280), [Issue #299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/299), [PR #328](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/328), [Issue #347](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/347)) — declarative rules in memory don't survive the action-distance from session start to `gh pr merge`. Mechanism-over-memory ladder rung 3 (memory → inline Always-in-effect → mechanism) per `shared/feedback_mechanism_over_memory.md`. Companion gate to the existing `@claude`-mention PreToolUse hook ([Issue #360](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/360), closed 2026-05-15).

**Adoption path:** CLAUDE.md "Merge workflow" section names the script as the operational path; bare `gh pr merge` documented as bypassing the gate. shared/MEMORY.md `feedback_closure_ritual.md` pointer line updated to lead with the script. The gate doesn't fire on `gh pr merge` directly — it's opt-in via the wrapper, so the discipline is "run the wrapper" rather than "the wrapper auto-fires on merge". Hooks-based auto-fire on `gh pr merge` is a possible follow-up if wrapper-call discipline slips, but not in scope for v1.

---

### 14:49 UTC — Editor: Developer

**Headline:** [Issue #325](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/325) (post-merge closure-audit critic) closing as KEEP — 7-day trial outcome one day ahead of the 2026-05-18 endpoint. Bot continues running; 4/11 production flags surfaced real persistent signal (3 priority-rationale gaps, 1 lab-notebook gap); kill criterion ("adds noise without catching anything the closure ritual doesn't catch within minutes") not met.

**Why this is the final entry on the trial:** [PR #332](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/332) shipped the v1 lean-experiment implementation 2026-05-11; the issue was reopened same day as a trial-tracker with a 2026-05-18 kill/keep decision endpoint. The [2026-05-13 trial summary](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/325#issuecomment-4442526736) (issue comment) already documented the evidence and recommended KEEP. Today's action finalizes that recommendation and closes the trial-tracker.

**Decision rationale recap (from the 2026-05-13 summary):**

| Flag category | Count | Real signal? |
|---|---|---|
| AC checkboxes | ~10 | Transient by read-time, but acting as workflow nudge — [PR #368](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/368) credited the audit explicitly for prompting a post-merge fix. Not noise. |
| Priority rationale | 3 ([Issue #361](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/361), [Issue #351](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/351), [Issue #347](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/347)) | All 3 still missing — persistent signal that would otherwise have stayed invisible. |
| Lab notebook ref | 4 | 1 persistent ([scientist nb, PR #343](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/343) `#342` ≠ `#343` strictness gap); 3 driven to fix. |

Total persistent signal: 4/11. Snapshot-staleness on the transient flags is not noise — v1 spec explicitly accepted that ("Out of scope: comment edit-in-place after fixes"), and the workflow-nudge effect was the desired secondary behavior.

**Parked follow-ups (both deferred, neither blocking):**

- **Trigger delay 10–15 min** — would silence transient AC-checkbox flags but only worth doing once snapshot-staleness dominates reader experience. Re-evaluate after [Issue #357](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/357) (`audit_and_merge.sh`) lands, since pre-merge enforcement should reduce post-merge noise.
- **Broaden lab-notebook number-match** — accept PR # OR parent-issue # OR sub-issue #. The single persistent lab-notebook gap is below the cost threshold for the relaxation.

**Process notes:**

- Lab notebook entry precedes the issue close per `shared/feedback_lab_notebook.md` rule (entry on main before close).
- ACs all ticked at issue creation since they document the v1 lean-experiment scope, which [PR #332](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/332) delivered — no re-tick needed at trial-close.
- Closes the [`dev-i1 - Dev Tooling Quick Wins`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/milestone/23) milestone pending [Issue #357](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/357) — 2 of 7 milestone issues remained; this one resolves via trial-close, [Issue #357](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/357) remains as real implementation work.

---

## 2026-05-16

### 20:18 UTC — Editor: Developer

**Headline:** [Issue #360](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/360) (pre-flight at-mention grep) shipped via [PR #380](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/380) — a rung-3 mechanism-over-memory PreToolUse hook. Bot review surfaced one real bug + one real coverage gap + one threat-model misread; hardened via `6469096` (try/except fail-open, `gh (pr|issue) edit --body` coverage, redundant guard removed, 12 subprocess tests added under [`tools/ci/`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/blob/feat/dev-360-at-claude-grep/tools/ci/test_check_at_claude.py)). Pushed back on the "Critical" finding with explicit threat-model reasoning (hook protects local Claude Code Bash sessions, not the bot itself).

**Why this work happened:**

The bot-mention rule (`memory/shared/feedback_no_at_claude_mention.md`) was inlined into shared Always-in-effect at MEMORY.md line 30, yet broke twice on the same shape — [PR #359](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/359) (2026-05-13) and [Issue #272](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/272) (2026-05-06 originating). Per `memory/shared/feedback_mechanism_over_memory.md`, when a rule breaks ≥2× on the same shape despite being inlined, escalation goes from memory → mechanism. The trigger lives many tool calls deep from session start — by the time the agent is composing a `gh issue comment --body "..."`, the at-mention rule has fallen out of working attention. A PreToolUse hook intercepts at the only point where intent and action align (the `gh` Bash invocation itself), so it cannot be forgotten regardless of how much downstream context has accumulated.

**Implementation reasoning:**

Option B (PreToolUse hook) over Option A (wrapper script): zero memory-discipline overhead. Option A would have required every agent to remember "use `gh_safe_comment.sh` instead of bare `gh`" — exactly the kind of memory-discipline rule that already failed twice. Option B fires automatically on every `gh` Bash call regardless of which command form the agent reached for. Per-machine config trade-off (`.claude/settings.json` is project-committed, but the harness has to be configured to load it — which it is in every Claude Code session) is acceptable.

The hook is intentionally narrow in scope:

```python
re.search(r"(^|[;&|]\s*)gh\s+(pr|issue)\s+(comment|create|edit)\b", cmd)
```

Matches `gh pr|issue` followed by `comment|create|edit` (the three body-writing subcommands), anchored to start-of-command or after a shell separator (`;&|`) to avoid matching `gh` substrings inside unrelated commands. The downstream regex check then re-filters for the bot-username substring, with a canonical-exception clause for the exact `--body "<bot> review"` form (the legitimate way to trigger the review bot on demand).

**Bot review findings — triage:**

| # | Finding | Triage |
|---|---|---|
| 1 | Critical: hook may not fire in GitHub Actions context | **Pushback.** Threat-model mismatch — hook guards local Claude Code sessions, not the bot. The bot posts via the Action's machinery (GitHub API), not Claude Code's Bash tool, so the hook has nothing to enforce on the bot's surface. The "deleted in working tree" observation also self-contradicts: bot quotes `.claude/settings.json:10`/`:11` content while claiming the file is absent. |
| 2 | Bug: no try/except on `json.load(sys.stdin)` | **Accepted + fixed.** Wrapped in `try/except (JSONDecodeError, ValueError) → return 0`. Fail-open is correct — the regex check is the secondary defense, and a hook traceback on malformed harness input would block every `gh` call without clear remediation. |
| 3 | Design: relative command path | **Pushback.** Standard Claude Code hook form; paths resolve from project root. |
| 4 | Design: over-broad `"if": "Bash(gh *)"` | **Pushback.** Script exits in <10ms for non-gh commands; cost is negligible. Tightening would be marginal. |
| 5 | Minor: redundant `total > 0` clause | **Accepted + fixed.** Preceding `if "<bot>" not in cmd: return 0` already guarantees ≥1 match. |
| 6 | Missing: `gh pr edit --body` / `gh issue edit --body` | **Accepted + fixed.** Real gap — `edit --body` overwrites a body and bypasses the comment/create guard. Extended regex `(comment\|create)` → `(comment\|create\|edit)` plus 2 new deny-path tests. |

**Test approach:**

Added [`tools/ci/test_check_at_claude.py`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/blob/feat/dev-360-at-claude-grep/tools/ci/test_check_at_claude.py) with 12 subprocess tests covering allow / deny / fail-open paths. Tests invoke the hook the same way the harness does (pipe PreToolUse JSON to stdin, read deny decision from stdout) — so any wiring drift (e.g. a future settings.json refactor that breaks the hook command path) is caught by `ci-tools-pytest`. Subprocess form was chosen over import-based testing because the hook lives under `.claude/hooks/` (non-standard sys.path) and the subprocess mirror is more honest about what's actually being tested.

**Process notes:**

- Closure ritual applied in full: all 5 ACs on [Issue #360](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/360) ticked before merge (`gh issue edit 360 --body-file ...`).
- Confirmed PR auto-closes [Issue #360](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/360) via `closingIssuesReferences` (PR title contains `closes #360`) — no manual `gh issue close` needed.
- Entry written post-review, pre-merge, per the Always-in-effect rule "Lab notebook entry comes AFTER review, before merge — not before commit" — reflects the final post-review state including the `6469096` hardening commit.
- Sister mechanism [Issue #357](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/357) (`audit_and_merge.sh` — closure-ritual enforcement at merge time) remains open — same shape, same rung-3 escalation, different rule.

**Open follow-ups:**

- None blocking. The `--body-file` foot-gun (a file containing a stray bot at-mention would bypass the regex check since the command line only carries the path) is a known limitation already present for `comment`/`create` and not flagged by the bot — out of scope here.
- If the hook ever needs to also guard `gh release create --notes "..."`, the regex extends naturally.

---

## 2026-05-15

### 18:40 UTC — Editor: Developer

**Headline:** [Issue #370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370) (P0/critical: regtools BED12 anchor-outer encoding shifted every HISAT2-path junction by 100–150 bp on each side; every neoepitope ever predicted on that path was derived from wrong coords) fixed via [PR #372](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/372) merged (squash `ca7ec79`). Two follow-up issues filed for the deferred ACs: [Issue #377](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/377) (CI canary cross-check vs `regtools junctions annotate`) and [Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378) (patient_002 PoC re-run for empirical verification). Closure-audit ritual applied in full despite the merge slip (see Process notes).

**Root-cause reasoning:**

The buggy code at [`alignment.smk:238-240`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/blob/ca7ec79/workflow/rules/alignment.smk) was:

```awk
awk -F'\t' '{if ($5 > 0) print $1":"$2":"$3":"$6"\t"$5}' {output.bed} > {output.junctions}
```

`regtools junctions extract` emits BED12 where cols 2-3 (`chromStart`/`chromEnd`) are the **anchor outer boundaries** of the spliced-read pile-up — not the intron donor/acceptor. The actual intron coords are inset by the anchor lengths recorded in `blockSizes` (col 11) and `blockStarts` (col 12):

```
chr1  16914  17368  JUNC00000004  177  -  ...  2  141,136  0,318
                                                  ^^^^^^^   ^^^^^
                                                  anchors   block starts
donor    (0-based) = 16914 + 141 = 17055
acceptor (0-based, exclusive) = 16914 + 318 = 17232
```

The pipeline was emitting `chr1:16914:17368:-` — shifted by 141 bp on the left and 136 bp on the right. Every annotated junction missed the GENCODE reference; the funnel collapsed to `annotated = 0, tumor_exclusive = 57708` on patient_002 (biologically impossible — any human RNA-seq sample hits tens of thousands of annotated junctions).

**Why it went undetected for as long as it did:** the funnel-stats arithmetic checked out (`junctions_raw = mean_reads_filtered + annotated + normal_shared + tumor_exclusive` summed correctly), `bedtools getfasta` happily extracted *some* sequence from the wrong coords, MHCflurry scored the resulting peptides without complaint, and the `annotated_discarded` field was only added to the funnel output in 2026-04-25 ([Issue #214](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/214)) — so a value of zero had nowhere to be seen for ~5 months. Arithmetic was sanity-checked; biology was not.

**Approach considered (Issue #370 body, summarized):**

| Option | Description | Rejected because |
|---|---|---|
| (A) Fix HISAT2 awk + use STAR col 6 | Apply blockSizes formula on HISAT2; read STAR's `annotated` flag from `SJ.out.tab` col 6 directly | Selected — smallest diff, no new tool |
| (B) Migrate to `regtools junctions annotate` | Pipe BED12 through annotate; consume `known_junction` flag | Annotate's `end += 1` quirk ([regtools/junctions_annotator.cc:66-83](https://github.com/griffithlab/regtools/blob/master/src/junctions/junctions_annotator.cc#L66-L83)) means the output `end` is +1 vs every other convention in the pipeline; requires a shim AND a `SJ.out.tab → BED12` converter for the STAR path |
| (B') | Use annotate for the flag only, do our own coord math | Eliminates the shim by ignoring regtools' coord output — at which point we keep our own coord code anyway |

(A) won. (B) was filed as a future option in [Issue #377](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/377) (CI canary cross-check) — we get the safety-net signal without committing to migration.

**Implementation:**

Extracted the awk replacement into [`workflow/scripts/bed12_to_junctions.py`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/blob/ca7ec79/workflow/scripts/bed12_to_junctions.py) rather than inlining the corrected math back into the `.smk` shell block. Reasons:

1. **Testability.** The arithmetic is the entire correctness story — needs unit tests with controlled fixtures, not "did the pipeline run?"
2. **The bug was a semantic misread by a human reading the awk.** A Python helper with named variables (`donor_0based`, `acceptor_0based_exclusive`) and a comment explaining the genomic-orientation convention is harder to misread than `$2 + bs[1]` / `$2 + bst[2]`.
3. **Tested coverage is wider than possible inline.** 8 tests: core regression, zero-read skip, minus strand, multiple records, GENCODE-match round-trip (the test that would have caught the original bug), unstranded `.` records, malformed short lines, comment/blank lines.

**Format round-trip verified end-to-end:**

```
bed12_to_junctions emits:  chr22:101:200:+   (1-based donor, 0-based exclusive acceptor)
filter_junctions._parse_junction_id:
    start = int("101") - 1 = 100   # 0-based intron start
    end   = int("200")     = 200   # 0-based exclusive intron end
GENCODE BED reference: chr22  100  200  ref_junc  0  +
                        → match ✓
```

The Claude review traced this trip independently and verified it. STAR path is unaffected — `SJ.out.tab` cols 2-3 are 1-based intron donor/acceptor directly, and the existing STAR awk in [`alignment.smk:347-352`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/blob/ca7ec79/workflow/rules/alignment.smk) is accidentally correct (1-based inclusive `[a,b]` → 0-based half-open `[a-1, b)` requires `-1` on start, no change on end — which is exactly what `_parse_junction_id` does).

**Bot review cycle:**

The [Claude review](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/372#issuecomment-4462167720) verdict was "correct and safe to merge" with three minor non-blocking notes. All three were addressed in two follow-up commits before merge:

- **`_snakemake_main()` is dead code** — removed in `c3b7f54`. The rule invokes the helper via `shell:` + CLI, so the `snakemake` global never exists at runtime and the `try/except NameError` always fell through to `_cli_main()`. Reviewer was right; design-pattern carryover from another script in the codebase.
- **Donor/acceptor naming is genomic, not transcript-oriented** — clarifying comment added in `c3b7f54`. On `+` strand, `donor_0based` is the biological 5' splice site; on `-` strand it's the 3' splice site. GENCODE BED uses the same genomic convention so downstream matching is symmetric — no correctness change, but a confused future reader is now disambiguated by the comment.
- **Missing defensive-guard tests** — added in `74e08e5`. Three new tests: unstranded `.` round-trip, lines with `<12` fields skipped, comment + blank lines skipped. Each guard was already in the code; they just weren't covered.

**Closure-audit ritual:**

Four of five [Issue #370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370) ACs were not met at merge-time. Handled per the `shared/feedback_closure_ritual.md` rule:

| AC | Status | Disposition |
|---|---|---|
| CLAUDE.md gotcha note | ✅ Met | Ticked in issue body |
| `annotated_discarded` plausible fraction on re-run | ❎ Empirical, post-merge | Deferred to [Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378) |
| `tumor_exclusive` drops to single-digit-thousands | ❎ Empirical, post-merge | Deferred to [Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378) |
| Spiked annotated junction unit test (both paths) | ❎ Partial | HISAT2 covered by `test_annotated_junction_matches_gencode_reference`; STAR path was never broken — deferred to [Issue #375](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/375) (STAR col 6 emission) |
| CI canary regtools annotate cross-check | ❎ Out-of-scope per original body | Tracked in [Issue #377](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/377) |

Comment posted on [Issue #370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370#issuecomment-4462324530) makes the deferrals visible — the audit trail is the comment, not silent omission.

PR test plan: ticked the two CI boxes (`pipeline-pytest`, `pipeline-snakemake-dry-run` from run `25924541450`), updated the unit-test box to reflect the new test count (262 / 25 skipped after polish commits, vs the body's pre-polish 259 / 25), and reformatted the post-merge re-run line as `- [x] (post-merge: tracked in #378)` per the closure-ritual pattern for explicitly-post-merge verification.

**Process notes / lessons:**

- **I shipped the work before writing this lab notebook entry — that's the failure mode the `feedback_lab_notebook.md` rule was designed to prevent.** When I drafted the closure-audit plan, I listed the lab notebook entry as the *last* step, after merge + board status updates. User caught it ("why are we doing lab notebook entries post merge now? Is that written in memory?"). The rule was explicit; I had forgotten it because it wasn't surfaced in my MEMORY.md Always-in-effect block. Slip acknowledged, rule re-read, and the rule itself was updated to clarify the ordering: **commit → push → open PR → review(s) → notebook entry → merge** (previously stated as "notebook entry → commit → push → merge" which implied notebook-first; the actual intent is that the entry captures the post-review final state but still ships *in the same PR*). The entry-after-review framing was the user's framing — better than the old one because it lets the entry reflect review-driven changes rather than the pre-review draft state.
- **This entry is retroactive.** Per the rule it should have been bundled into [PR #372](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/372) between the polish commits and merge. Since the PR is now squash-merged at `ca7ec79`, it lands on its own `docs/developer/lab-notebook-2026-05-15-1840` branch. Future me reading this entry should see the "[PR #372 merged]" timestamp + "[notebook entry merged]" timestamp diverge and recognize this as the recorded slip, not a new pattern.
- **The "helper file vs inline" call paid for itself within the same review.** Three of the reviewer's three notes were file-level concerns (dead-code branch, naming comment placement, defensive-guard test coverage) that simply don't exist when the math is inline awk in an `.smk` shell block. Inline awk has no place to put a `Path`-typed function signature, an `argparse` entry point, or a `_VALID` enum guard. Extracting the script unlocked the review feedback.
- **TDD on the helper was tighter than I usually run.** Wrote the failing import test first (`ModuleNotFoundError`), then the minimal `convert_bed12_to_junctions` signature, then the GENCODE-match scenario. The minimum-viable implementation didn't need an `if` branch on `reads <= 0` until that test was added — small, but it kept dead defensive code out.

**Cross-references:**

- [Issue #370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370) (root P0 bug) — closed by [PR #372](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/372)
- [PR #372](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/372) (the fix) — squash-merged `ca7ec79`
- [Issue #377](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/377) (CI canary cross-check follow-up, P2)
- [Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378) (patient_002 re-run follow-up, P1)
- [Issue #374](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/374) (STAR strand=0 motif rescue, filed earlier today, P1)
- [Issue #375](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/375) (STAR col 6 annotated flag emission, filed earlier today, P2)
- [PR #371](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/371) (separately-merged `chore: ignore /memory symlink for personas tooling`) — landed before the [Issue #370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370) work began

---

## 2026-05-13

### 14:40 UTC — Editor: Developer

**Headline:** [PR #358](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/358) ([Issue #279](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/279) HISAT2 strandness) merged (squash `5c485f9`). Bot re-review surfaced two real correctness improvements applied + one bot-suggested Snakemake-7 idiom that broke CI and got reverted. [Issue #364](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/364) filed for an unrelated local pytest collection hang surfaced during verification.

**Bot review cycle:**

- **Item 1 — `srcdir()` over `workflow.basedir`:** bot suggested `sys.path.insert(0, srcdir("../scripts"))` as the more idiomatic Snakemake form. Applied + pushed; CI `pipeline-snakemake-dry-run` failed with `NameError: name 'srcdir' is not defined` at parse time. `srcdir()` was a Snakemake-7 helper that is no longer exposed at the `.smk` module scope in Snakemake 8. Reverted to the `workflow.basedir`-based form. Cost: one CI cycle.
- **Item 2 — raise on typo:** bot flagged that silent fallback to `""` (unstranded) on any unrecognized string was a real correctness risk — a typo like `forwrd` would silently produce wrong alignment with no log signal. Applied: introduced `_VALID = {"", "unstranded", "forward", "reverse"}` + `raise ValueError` on anything else. Backward compat preserved for `None` / empty / whitespace / `unstranded`. Tests flipped to `pytest.raises(ValueError, match=...)`. Stronger coverage than before.
- **Item 3 — extract `get_strandness_from_row(row: dict)`:** bot recommended moving the `is_pe` detection out of the Snakemake glue layer into a second pure function. Applied: `_get_hisat2_strandness` in `alignment.smk` reduced to a row-lookup + delegate. Added 5 row-level test cases (PBMC forward SE, PE unstranded, missing `strandness` key, whitespace `fastq2`, PE reverse) grounded in real `samples.tsv` rows. Total test count: 11 → 20.

**Downstream documentation:**

- **CLAUDE.md** gained a "Snakemake 8 Gotchas" section grouping the existing `--configfile`-flag-collapsing note with the new `srcdir()` unavailability note. Records what was tried (with the failing CI evidence) so a future code-review bot or contributor doesn't re-suggest the same fix.
- **`alignment.smk`** has an inline comment at the import block explaining the `srcdir()` constraint for anyone reading that file directly.
- **PR #358 re-review** (post-revert) approved the work without changes; one non-blocking suggestion — scoping the strandness import + `_get_hisat2_strandness` helper inside the `if aligner == "hisat2":` guard so they don't run on STAR-selected configs — was applied in commit `a4e8ebb` as a pure relocation.

**Issue #364 filed (separate concern):**

Local pytest collection hung > 5 min during verification — even a bare `import pandas` hung > 2.5 min. System diagnostics showed only ~60 MB free pages + load average 3.4 on the 8 GB M1 box. Hypothesis: module-level heavy imports across many test files (pandas, mhcflurry, etc.) compound under memory pressure and tip pytest collection into disk-thrashing. CI runners have more RAM so it's a local-only friction. Issue body proposes a per-file import time profile + selective hoist into pytest fixtures.

**Process notes / lessons:**

- **Bot framework-idiom claims need a local smoke-test before commit.** The `srcdir()` suggestion was plausible, well-cited, and CI-evident-wrong. A 5-second local `snakemake -n` parse check would have caught it before the push. The `receiving-code-review` skill principle ("verify before agreeing") applies even when the reviewer is the bot and the claim is about a documented built-in.
- **Inline a discovered gotcha in two places** — the file that has the workaround (so a code reader sees the rationale) AND the project-wide CLAUDE.md (so anyone writing a *new* `.smk` doesn't re-discover it). Did both this round.
- **Local pytest is not a reliable gate today.** Per the test-before-PR rule I want to run the full suite locally, but the system memory pressure makes that unreliable. Workaround: ran only the touched test file (20 cases, 0.02 s), pushed to let CI cover the rest. Tracked under [Issue #364](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/364) so this isn't ad-hoc forever.
- **One self-violation worth naming.** When I applied the bot's review fixes, I committed all 3 files (helper + tests + smk) in a single commit. The multi-file-feature-workflow rule says one-file-per-commit. Justified here — each file in isolation breaks the import chain (helper ↔ tests ↔ consumer), so split commits would each red CI — but the justification belongs in the commit message and was there.

---

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
