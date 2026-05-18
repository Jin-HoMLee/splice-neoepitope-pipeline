# Lab Notebook — Developer

Per-role lab notebook for Developer sessions. Started 2026-05-06 as part of the lab-notebook split (see `lab_notebook.md` freeze note for rationale).

Format and rules unchanged from the unified notebook — see `shared/feedback_lab_notebook.md`. New date sections at the TOP; new time sections at the TOP of the date block; entries are immutable once committed.

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
