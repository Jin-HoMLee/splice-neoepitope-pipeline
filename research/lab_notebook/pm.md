# Lab Notebook — PM

Per-role lab notebook for PM sessions. Started 2026-05-06 as part of the lab-notebook split (see `lab_notebook.md` freeze note for rationale).

Format and rules unchanged from the unified notebook — see `shared/feedback_lab_notebook.md`. New date sections at the TOP; new time sections at the TOP of the date block; entries are immutable once committed.

---

## 2026-05-18

### 10:36 UTC — Editor: PM

#### Morning routine + dev-i1 milestone close + [Issue #247](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/247) Phase 0 ship — PostToolUse recheck hook

Monday morning routine ran clean: PM news (zero items — territory swept yesterday), closure audit (4/4 closed-issue checks passed; one soft observation on [Issue #357](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/357)'s missing lab notebook entry, folded into a standup celebration broadcast at 08:52 UTC), board recap, triage. Weekly full-board sweep surfaced 4 field-incomplete issues, all triaged auto + 1 by-prompt ([Issue #345](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/345) → i5-S3 / P3 / XS; [Issue #352](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/352) → **new dev-i2 milestone** / P2 / XS; [Issue #364](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/364) → dev-i2 / P3 / S; [Issue #394](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/394) → P3 / XS).

**`dev-i1` milestone closed 7/7 one day early.** Created `dev-i2 - Dev Tooling Quick Wins II` (milestone 25, due 2026-06-18) for Dev's follow-up tooling axis. Posted [`team_standup.md` celebration broadcast at 08:52 UTC](memory/shared/team_standup.md) addressed to Developer (cc: Scientist) covering: 7/7 closed, 0 carry-over, the complete closure-ritual safety net (pre-merge gate [Issue #357](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/357) + post-merge detective [Issue #325](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/325) + [Issue #360](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/360) at-mention guard = 3 rung-3 mechanism-over-memory escalations in a single iteration), and the soft observation that [Issue #357](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/357)'s ship merits its own time section in [`developer.md`](research/lab_notebook/developer.md) (substance is in the [2026-05-17 20:15 UTC standup post](memory/shared/team_standup.md) but per the spirit of the closure-ritual rule the lab notebook should mirror it).

**pm-i1 closing-run survey.** Burndown was 4 closed / 3 open with due Thu 2026-05-21. Opted for 6/7 shape: land [Issue #247](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/247) (P1, fully scoped) + [Issue #243](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/243) (P2, fully scoped); slip [Issue #324](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/324) (per-role model routing — M, still skeletal) to `pm-i2 - PM Self-Improvement Tooling`. Per-role model routing belongs in pm-i2 by topic anyway; the slip is honest, not corner-cutting. Audit-trail comment posted on [Issue #324](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/324) explaining the move ([comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/324#issuecomment-4476433931)).

**[Issue #247](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/247) Phase 0 implementation — recheck atom + dispatch hook + settings.local.json wiring.**

Three artifacts shipped on branch `feat/pm/issue-247-milestone-recheck-hook`:

1. [`scripts/pm/recheck_milestone.py`](scripts/pm/recheck_milestone.py) (~110 lines) — recheck atom. CLI: `--issue N` or `--milestone N`. Sums size-weighted days from open issues (XS=0.5, S=1, M=2.5, L=3.5, XL=5 capacity-days), applies `new_due_on = today + remaining/1.5 × 7 days`, prints standardized output, exits 0 (`[No change]`) or 2 (`[UPDATE NEEDED]` when |delta| > 7 days) per AC spec. Single-query GraphQL with aliases for per-issue Size lookup — efficient enough for any milestone size.

2. [`.claude/hooks/recheck_milestone_dispatch.py`](.claude/hooks/recheck_milestone_dispatch.py) (~120 lines) — PostToolUse hook entry. Reads PostToolUse JSON from stdin, pattern-matches the Bash command against 4 trigger shapes, invokes the recheck script, emits `additionalContext` wrapped in `hookSpecificOutput` JSON. Patterns:
   - `\bgh\s+issue\s+close\s+(\d+)` → recheck issue's milestone (trigger 1, close)
   - `\bgh\s+issue\s+edit\s+(\d+)\b[^|;&]*--milestone\b` → recheck **both** milestones via `/issues/N/events` history lookup (resolves title → number) (trigger 2, move)
   - Size field ID (`PVTSSF_lAHOB17eGc4BSomPzhAHGiA`) present in command + `itemId: "PVTI_..."` → resolve item → issue → milestone, recheck (trigger 3, resize)
   - `gh api` + (`-X PATCH` OR `--method PATCH`) + `/milestones/N` (order-agnostic — initial regex required PATCH after milestones path, broke on real commands) → recheck milestone (trigger 4, due_on edit)

3. [`.claude/settings.local.json`](.claude/settings.local.json) (Phase 0, gitignored per AC spec) — added `hooks.PostToolUse` array with single Bash matcher + `if: "Bash(gh *)"` filter + 30s timeout. Phase 1 (promote to committed [`.claude/settings.json`](.claude/settings.json)) deferred to follow-up Issue.

**Three real findings surfaced by the recheck script** during smoke-testing — drift that would otherwise have stayed invisible:

| Milestone | Current `due_on` | Proposed | Delta | Action taken |
|---|---|---|---|---|
| `pm-i3 - PM Workflow Quick Wins II` | 2026-05-24 | 2026-06-06 | +13 days | ✓ PATCHed to 2026-06-06 (live test AC 6) |
| `i2 - S4 - AlphaGenome Predicted-Normal Filter Validation` | 2026-05-22 | 2026-07-04 | +43 days | ⚠ pending decision (scope vs date) |
| `pm-i2 - PM Self-Improvement Tooling` | 2026-06-11 | 2026-08-03 | +53 days | ⚠ pending decision (likely needs displacement to a pm-i4) |

The pm-i3 PATCH was a real `gh api -X PATCH` invocation — the hook fired live and surfaced the post-PATCH recheck as `additionalContext` confirming `[No change]` (delta now 0). That's AC 6 verified live; the harness reloads `settings.local.json` mid-session (good to know — Claude Code doesn't require a session restart for hook config changes).

**AC status at merge:**

- AC 1–3 (script, formula, hook config): ✓ shipped
- AC 4 (manual close test): pattern verified via dispatch unit test with mocked stdin (Test 1 above shows correct JSON output); live `gh issue close N` trigger deferred — rare in practice since most closes happen via `gh pr merge`'s `closingIssuesReferences` auto-close path, which is NOT caught by the current `\bgh\s+issue\s+close\s+(\d+)` pattern. Adding a `gh pr merge` trigger is a clean Phase 1 extension (parse the PR's linked issues, recheck each one's milestone) — surfaced as a known coverage gap, not landed in Phase 0.
- AC 5 (manual move test): pattern verified via dispatch unit test (Test 5 above — used [Issue #324](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/324)'s real move history; both pm-i1 source + pm-i2 destination rechecked correctly via `/issues/N/events` lookup). Live trigger deferred to next natural move event.
- AC 6 (manual due_on PATCH test): ✓ live trigger verified end-to-end as above.
- AC 7 (update [`memory/feedback_milestones.md`](memory/feedback_milestones.md)): ✓ reframed "Recheck `due_on`" section from memory-led to mechanism-led; also tightened the "Sub-issues closed" bullet (was awkward double-duty as both trigger and rationale snippet) to a cleaner "Issue closes" generalization per user catch.
- AC 8 (this entry): ✓

**Phase 1 follow-up candidates** (to be filed as separate Issues, not landed here):

- Promote hook to committed [`.claude/settings.json`](.claude/settings.json) — makes the recheck available to all 3 roles' worktrees.
- Add `gh pr merge` trigger covering the auto-close-via-`closingIssuesReferences` path.
- Add `gh issue reopen` as a counterpart to close (remaining capacity goes back up).
- **Heredoc false-positive guard** — surfaced live during this Issue's own ship. The `gh issue edit 247 --body-file ...` command (with the body file written via `cat > ... <<EOF` heredoc inline) caused the dispatch script to match the close pattern + due_on PATCH pattern against text inside the heredoc body (the body documents these very triggers). Both rechecks were harmless `[No change]`, but the false-positive is real. Mitigation in this session was to run `gh issue edit --body-file <existing-file>` as a standalone Bash command (no inline heredoc) — the standalone command string is clean. Phase 1 hardening should detect heredoc structure and exclude its content before pattern matching.
- Optional: pytest coverage for the dispatch script (mirroring [`tools/ci/test_check_at_claude.py`](tools/ci/test_check_at_claude.py) pattern from [Issue #360](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/360) — clean precedent).

**Process notes:**

- Lab notebook entry written before merge per closure ritual; PR opens next.
- Plan to ship via [`scripts/audit_and_merge.sh`](scripts/audit_and_merge.sh) — inaugural PM-side use of yesterday's [Issue #357](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/357) closure-ritual gate.
- Three "UPDATE NEEDED" findings on milestones (pm-i3, i2-S4, pm-i2) are independent of [Issue #247](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/247)'s ship — but the script's job is to make them visible, and one (pm-i3) was actioned during the smoke test. The other two need separate triage decisions (scope reduction on i2-S4 with Sci; capacity displacement on pm-i2 — perhaps create pm-i4 to absorb overflow).

---

## 2026-05-17

### 19:12 UTC — Editor: PM

#### [Issue #305](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/305) closure — landscape-doc tightening, not memory edit

Issue scope said tighten [`shared/feedback_multi_role_not_multi_agent.md`](memory/shared/feedback_multi_role_not_multi_agent.md) against Anthropic Managed Agents (2026-05-07) OR comment-defer. Initial proposal was a "Concrete contrast" paragraph inside the memory file. **User correction:** the contrast already lives in [`research/multi_agent_landscape.md`](research/multi_agent_landscape.md) (Frameworks → Anthropic Managed Agents) — the memory's job is rule-enforcement, the landscape doc's job is concrete tracking. The memory's own footer cross-references the landscape doc; that pointer was the signal I missed.

**Landed instead.** Tightened the existing Managed Agents entry in [`research/multi_agent_landscape.md`](research/multi_agent_landscape.md) with per-feature mapping: Multiagent Orchestration ↔ multi-role-not-multi-agent autonomy bar; Dreaming ↔ `/memory` rehydration (auto vs human-initiated); Outcomes ↔ no direct analog. Added the third feature (Outcomes) that the prior entry missed; added 9to5Mac source alongside the Anthropic blog primary. Memory file unchanged.

**Adjacent finding — `/cerebrum` → `/memory` rename incomplete.** Landscape doc still used `/cerebrum` (old skill name); fixed in this edit. ~18 other `/cerebrum` references remain across the personas repo (`shared/team_memory_broadcasts.md` at minimum) — out of scope for [Issue #305](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/305), surfaced for follow-up Issue.

**Process slip recap.** Two user corrections in this issue: (1) initial session-opening recommendation surfaced a Developer PR while in PM session (role-boundary slip, caught with *"but PR 389 is the Developer's work"*); (2) initial proposal placed concrete contrast in the memory file when it belonged in the landscape doc (adjacent-artifact slip, caught with *"why did we write this whole memory file again? We also have the research/multi_agent_landscape.md file"*). Sister failures — both about *where* content/work lives in the artifact map. [`feedback_best_next_issue.md`](memory/shared/feedback_best_next_issue.md) Step 1 already covers slip (1); the duplicate-check rule ([`feedback_memory_duplicate_check.md`](memory/shared/feedback_memory_duplicate_check.md)) could be extended to cover slip (2) by adding an adjacent-artifact scan. Candidate follow-up memory-update; not landed here.

---

## 2026-05-16

### 20:08 UTC — Editor: PM

#### Saturday afternoon morning-routine (news section skipped, Friday cleanup skipped — Saturday)

**Trigger.** User opened with *"good afternoon"*, then *"maybe let's do a bit of morning routine without the news section"*. Routine reduced to: closure audit → board recap → standup → triage. Only closure signal: [Issue #370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370) (the P0 regtools coords bug closed yesterday by [PR #372](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/372)).

#### [Issue #370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370) closure audit — two-layer catch

Two distinct audit loops fired on yesterday's close at 18:38 UTC:

- **Closure-audit bot at 18:38 UTC** caught the priority-rationale + lab-notebook gaps. Dev resolved both within 25 min ([PR #379](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/379) for the lab notebook, body edit for the rationale).
- **My follow-up at 19:43 UTC** caught the AC gap — 4 of 5 ACs unticked. Posted as standup nudge + issue comment.

The bot caught what it can (mechanical presence checks — does the body have a `**Priority rationale:**` line? does the lab notebook have today's date header?); the human-PM caught what only the human-PM can (whether the AC content actually matches the merged work). Two complementary layers — neither alone covers both classes.

**Stale-clone false-alarm sub-finding.** Initial grep on [`research/lab_notebook/developer.md`](research/lab_notebook/developer.md) showed no `## 2026-05-15` header, suggesting the bot's flag was still unaddressed. Wrong — local clone was 2 commits behind `origin/main` ([PR #379](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/379) had merged but I hadn't pulled). Fast-forwarded `workspace/pm` to `origin/main` and re-verified — entry is at the top. **Lesson:** before declaring any file "still missing X" during audit, verify the working tree's git state against `origin/main`, not just the local file contents. Adjacent to [`feedback_read_before_claiming.md`](memory/shared/feedback_read_before_claiming.md), one level up (working-tree freshness, not file-content freshness).

#### Sub-issue retroactive linkage to a closed parent — user correction

Initially proposed retroactively linking [Issue #374](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/374), [Issue #375](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/375), [Issue #377](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/377), [Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378) as sub-issues of the (closed) [Issue #370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370) to give the unmet ACs natural defer targets. User pushed back: *"but Issue 370 is already closed. Is that ok?"*

GitHub mechanics allow it (linkage is structural, not status-bound). The conceptual fit fails. Per [`feedback_parent_sub_issues.md`](memory/shared/feedback_parent_sub_issues.md), **linkage = scope, not dependency**. [Issue #370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370) was filed as a single coherent P0 bug, not an epic with sub-scope. The 4 follow-ups are *related-discovered* (surfaced during the same debugging session) but each has independent scope. Retroactively turning [Issue #370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370) into a parent would muddle the closure narrative: closed parents read as "scope complete"; adding open sub-issues makes the closure look premature.

**Right move instead:** comment-defer the unmet ACs on [Issue #370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370) pointing at the specific follow-up issues — same audit-trail visibility, no closed-parent semantics drift.

**Rule gap surfaced.** [`feedback_parent_sub_issues.md`](memory/shared/feedback_parent_sub_issues.md) covers the scope-vs-dependency principle but doesn't explicitly call out the closed-parent edge case. Candidate for a future memory-update PR (not this entry's scope).

#### Milestone capacity — count vs size-weighted days — second user correction

Proposed putting all 4 follow-ups in [`i3 - S3 - Data Preparation - Aligner & Input Format Improvements`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/milestone/11) on topical-fit grounds, citing "6 open / 10 total — moderate but not bloated." User: *"do we limit the milestone sizes by number of issues or by combined estimated duration (days) now?"*

Re-ran with the correct metric. Per [`feedback_milestones.md`](memory/feedback_milestones.md): iteration budget = ~5d size-weighted (XS≈0.5d, S≈1d, M≈2.5d, L≈3.5d, XL≈5d).

[`i3 - S3 - Data Preparation - Aligner & Input Format Improvements`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/milestone/11) committed work — open + closed counts (total committed gates capacity, not remaining):

| # | Size | Days | State |
|---|---|---|---|
| [Issue #17](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/17) | L | 3.5 | Open (Ready) |
| [Issue #127](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/127) | M | 2.5 | Done |
| [Issue #277](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/277) | XS | 0.5 | Done |
| [Issue #279](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/279) | XS | 0.5 | Done |
| [Issue #297](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/297) | S | 1.0 | Open (Backlog) |
| [Issue #370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370) | M | 2.5 | Done |
| **Total** | | **10.5d** | |

That's **2.1× the 5d budget** before adding any of yesterday's follow-ups. Adding 4 more (S+S+S+M = 5.5d) would push it to **16d → 3.2× budget**. The rule explicitly says "Once total reaches ~5d, the iteration is full — even if all current issues are closed. Don't refill freed capacity."

**Why I got it wrong initially.** Anchored on issue count + topical fit. Both real signals — but count is not the capacity gate. Size-weighted days is the rule. I applied the easier-to-eyeball signal and skipped the rule.

**Triggered propose-and-confirm.** Per [`feedback_ask_for_help.md`](memory/feedback_ask_for_help.md) "capacity within ~10% of cap" trigger (here: 2–3× over cap is well past it), surfaced three options via `AskUserQuestion`: (A) new iteration, (B) absorb + slip [`i3 - S3 - Data Preparation - Aligner & Input Format Improvements`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/milestone/11) due_on by ~3 weeks, (C) split — keep [Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378) in i3-S3, move others. User picked **A**.

#### Created [`i5 - S3 - Data Preparation - STAR Polish & Aligner Verification`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/milestone/24)

**Naming.** Per `i<N> - S<stage> - <Stage Name> - <Arc>`. Existing S3 iterations: i2 (GTEx), i3 (HISAT2/regtools), i4 (nf-core). Next unused = **i5**. Arc captures the two halves — STAR polish ([Issue #374](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/374) + [Issue #375](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/375)) + Aligner Verification ([Issue #377](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/377) + [Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378)).

**Capacity check.** New milestone load = 5.5d → 1.1× budget. Marginally over but well within rounding.

**due_on = 2026-06-12.** Parallel with [`i3 - S3 - Data Preparation - Aligner & Input Format Improvements`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/milestone/11) rather than sequenced after it. Rationale: [Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378) is the empirical AC verification for [Issue #370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370) (which lives in i3 - S3) — strict sequencing would leave the [Issue #370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370) fix unvalidated until i3 - S3 is already closed. P100 capacity uncertainty (CLAUDE.md flags sustained exhaustion in `europe-west1-b`) is the main risk on [Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378) — but that slips the issue, not the milestone.

#### Triage applied — all 4 → i5 - S3 - Data Preparation - STAR Polish & Aligner Verification

| # | Size | Priority change | Scope |
|---|---|---|---|
| [Issue #374](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/374) | S | P1 (unchanged) | STAR strand=0 silent contamination |
| [Issue #375](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/375) | S | P2 (unchanged) | STAR col 6 annotated flag |
| [Issue #377](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/377) | S | — → P2 | CI canary regtools annotate cross-check |
| [Issue #378](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/378) | M | — → P1 | patient_002 PoC re-run with [PR #372](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/372) fix |

All 4 already had `**Priority rationale:**` lines in their bodies at creation — Dev did the rationale work; PM propagated milestone + size + missing board priorities.

#### Two open follow-ups carried into next session

1. **Standup at 19:40 UTC yesterday — [Issue #370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370) ACs Pending.** Dev hasn't responded (~24h elapsed at write time). Re-raise threshold is >1d → re-raise Monday morning if still no response.
2. **Older drift from 2026-05-13** — [Issue #352](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/352), [Issue #357](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/357), [Issue #364](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/364), [Issue #365](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/365) still unmilestoned with no priority/size. Deferred to Monday's weekly full-board sweep.

#### Lessons carried forward

- **Run the right metric before triage.** Topical fit + issue count are useful signals but not the capacity gate. Size-weighted days is the rule.
- **Closed-parent retroactive linkage is conceptually wrong even though mechanically allowed.** Scope ≠ historical association.
- **Two-layer audit works.** Mechanical (bot) + judgment (human-PM) catch complementary classes.
- **Verify git state before declaring "X missing."** Local clone staleness is a real false-positive vector during audit.

---

## 2026-05-13

### 14:48 UTC — Editor: PM

#### [PR #344](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/344) merge — Managed Agents primary-source swap + body cleanup

**Trigger.** After [PR #363](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/363) merge, surveyed remaining open PRs: only [PR #344](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/344) (multi-agent landscape doc, open since 2026-05-12) still in flight. CI all green, 8/8 Test plan ticked, Claude review approved with one polish nit + body cosmetics (missing `**Created by:** PM`, deprecated `🤖 Generated with Claude Code` footer from older template). Offered user three paths: merge as-is, fix nit + clean body, or just body cleanup. User chose **path 2** ("go on 2.") — full polish before merge.

**Managed Agents source swap.** Original entry cited [9to5Mac coverage from 2026-05-07](https://9to5mac.com/2026/05/07/anthropic-updates-claude-managed-agents-with-three-new-features/) — third-party. Claude review's reasoning: *"for a doc framed as a portfolio artifact ('we surveyed the field; we're not cargo-culting'), citing a secondary source for the first Anthropic entry looks weaker than it should."* WebSearch on `claude.com / anthropic.com` returned the primary announcement at [claude.com/blog/new-in-claude-managed-agents](https://claude.com/blog/new-in-claude-managed-agents) (2026-05-06, one day before 9to5Mac coverage). One-line edit, committed as `9d66f17`. Symphony's "no public docs yet" line (Claude review nit 2) intentionally left as-is per the maintenance section's stale-watch convention.

**Body cleanup.** Rewrote PR body: `**Created by:** PM` at top per Always-in-effect rule, dropped the auto-generated 🤖 Claude Code footer (artifact of older PR template — the new template is the minimal one per [`shared/feedback_lab_notebook.md` "PR body" rule](https://github.com/Jin-HoMLee/cerebrum/blob/main/splice-neoepitope-pipeline/shared/feedback_lab_notebook.md)), added two new Test plan boxes for today's commits (Pattern Language 2026 morning add at `1cf8dc1`, Managed Agents swap at `9d66f17`). Body now passes the `**Created by:**` audit check.

**Merge.** Squash-merge after CI green-verified on the new commit. [Issue #337](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/337) auto-closed via `Closes` keyword at 14:43:41 UTC.

**Closure audit caught two more gaps.** Audit bot flagged: (a) AC checkboxes on [#337](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/337) — 8/8 unticked (4 Scope + 4 Acceptance criteria); (b) `## 2026-05-13` block in PM notebook didn't reference `#344` (only had the 14:07 UTC #244 entry). Same root-cause as the morning's #244 audit miss — PR-driven close doesn't auto-tick parent ACs. Fixed: all 8 ACs flipped to `[x]` via `gh issue edit --body-file`; this entry IS the lab-notebook fix for gap (b). The **Cross-repo PR-driven close** rule added to [closure_ritual](https://github.com/Jin-HoMLee/cerebrum/commit/f5bff8b) earlier this session would have caught (a) if I had applied it pre-merge — adjacent gap: **same-repo PR-driven close also has the tick-the-source-issue step**, just less obviously remote. Rule wording is already correct (line 22: *"before merging the closing PR, do the box-ticking on the parent issue"*) — the failure was application, not rule-clarity. Audit is doing its job.

### 14:07 UTC — Editor: PM

#### Post-lunch pickup: [Issue #244](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/244) — codify "flag uncertainty before executing" rule

**Pick rationale.** After lunch, scanned open PM Backlog and proposed two candidates: [#244](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/244) (P1, S, pm-i1, due 2026-05-21) and [#353](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/353) (P2, S, pm-i2, due 2026-05-27). Recommended #244 — closer deadline, P1 over P2, calibration value compounds across every future triage. User picked #244 with "ok gogo".

**Scope.** Spec called for new `pm/feedback_ask_for_help.md` codifying when to propose-and-confirm vs. auto-execute, with 3–5 worked examples from real recent triages, plus cross-ref from `pm/feedback_milestones.md` and index entry in `pm/MEMORY.md`. Cross-repo work: issue in splice-neoepitope-pipeline, deliverable in cerebrum repo.

**Mid-task pause for review.** First Write call was rejected with *"what is this now?"* — I had jumped from "ok gogo" straight to a ~120-line draft without surfacing the file content for review. Recapped scope inline, offered to paste-then-write / trim / pause. User said *"ah ok continue!"* — proceeded as drafted. Lesson lands inside the rule itself: even mid-implementation, when about to commit a substantial artifact for the first time, surface the draft before writing. Adjacent to the rule but one level down (artifact-level, not triage-level).

**Worked examples chosen.** Five, mix of ❌ should-have-flagged and ✅ correctly-handled to anchor calibration both directions: (1) 2026-05-13 #337 milestone misattribution → role-cut ambiguity ❌; (2) 2026-05-02 `<role>-i<N>` axis creation → first-of-pattern ✅; (3) 2026-05-13 #353 vocabulary adoption → scope expansion ❌; (4) 2026-05-13 #346 priority rationale backfill → don't-flag ✅; (5) hypothetical 4.5d→5d capacity cap → capacity pressure ✅. The ❌ cases are real misses from earlier today — concrete, dated, traceable to user corrections in the morning routine.

**Cerebrum PR shape.** Branch `feat/pm/issue-244-flag-uncertainty` cut from `origin/main` on cerebrum repo. 3-file commit (`pm/feedback_ask_for_help.md` new, `pm/feedback_milestones.md` cross-ref added in new `### Flag uncertainty before executing the decision tree` subsection right after the decision-tree priority rule, `pm/MEMORY.md` index entry under Role: PM). Excluded a pre-existing uncommitted edit to `shared/feedback_mechanism_over_memory.md` (in-flight work from another session — not mine to touch). Opened as [cerebrum PR #4](https://github.com/Jin-HoMLee/cerebrum/pull/4) with `Closes Jin-HoMLee/splice-neoepitope-pipeline#244` cross-repo reference + 5-item summary + worked-examples list.

**Claude review trigger.** User: *"ping claude for review pls"*. Posted `@claude review` as canonical trigger phrase per `shared/feedback_no_at_claude_mention.md` (the one legitimate `@claude` use — bot-triggered review on a PR comment, not a body/AC mention). Bot reviewed; user merged via squash.

**Cross-repo closure worked.** Auto-closed: [Issue #244](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/244) flipped to CLOSED / COMPLETED at 14:02:29 UTC; project board Status auto-flipped Backlog → In Progress (manual, at start) → Done (auto, on close). Remote feature branch auto-deleted on merge. Local branch left at `[gone]` marker — working tree had uncommitted edits to defer cleanup until user sweeps.

**Closure audit catch.** Bot flagged two gaps: AC checkboxes 4/4 unticked (the issue body still had `- [ ]` for all four, not auto-flipped by cross-repo PR merge) + no `## 2026-05-13` lab notebook header. Ticked all four ACs via `gh issue edit --body-file`; this entry is the lab-notebook backfill. The cross-repo close mechanism doesn't tick ACs in the source issue's body — worth knowing for future cross-repo work patterns (manual tick-step required even when the PR's `Closes` keyword auto-flips Status).

---

## 2026-05-12

### 14:32 UTC — Editor: PM

#### [Issue #337](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/337) scaffold — multi-agent orchestration landscape doc + PR-body minimalism rule

**Doc scaffold.** Created `research/multi_agent_landscape.md` with three sections: Frameworks (Symphony, Managed Agents, Fugu), Methodology / framing (Trends Report, Mason, Code Agent Orchestra), Our position. Each backfilled entry has source link, one-line summary, rationale tagged **Observe / Pattern-confirmation / Counter-position / Pattern borrow + roster reject / Direct reinforcement / Convergent vocabulary** — six distinct stances surfacing the asymmetry between borrowing the *pattern* (one orchestrator + N specialists) and the *roster* (PM/Sci/Dev) on a per-item basis. "Our position" reproduces the 3 framing rules inline (since the cerebrum repo is gitignored and not portfolio-readable) plus the 3-layer orchestration shape (human → PM → Sci/Dev).

**Cross-link side of the work.** Updated all 3 framing memories (`feedback_multi_role_not_multi_agent`, `feedback_domain_bespoke_roles`, `feedback_cerebrum_vs_project`) with a closing `**Cross-reference:**` block pointing at the landscape doc. Goal: when a future role re-reads a framing memory in isolation, they discover the landscape doc as the curated synthesis — and when they read the landscape doc, the inline rationales bottom out at the framing memories. Reciprocal pointers, no duplicated content.

**Reference + workflow hook.** Registered the landscape doc in `shared/reference_docs_inventory.md` as a stale-watch artifact. Maintenance hook lives in PM's `feedback_morning_routine.md` Step 0 ("How to apply"): when news_log gets a `methodology-signal` entry, scan against the landscape doc and pick (a) backfill, (b) update existing rationale, or (c) skip. No new reference memory file created — keeping rule-proliferation in check (same discipline as yesterday's [16:31 UTC] entry where I folded the news-cap into existing morning-routine files instead of spawning `feedback_news_issue_cap.md`).

**PR-body minimalism rule (mid-session capture).** User asked, *"do you feel the PR body is redundant?"* — looking at [PR #341](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/341)'s Summary section (5 bullets that mirrored the lab-notebook entry's `####` topic headings), yes. The lab notebook IS the artifact; markdown renders in the diff view; the Summary section was content re-serialization. User caught me drafting it as a new file (`feedback_lab_notebook_pr_body.md`) — same rule-proliferation slip I just praised yesterday's entry for avoiding. Reverted to extending `shared/feedback_lab_notebook.md` with a new "PR body — minimal, no content re-serialization" section + template + Why + How to apply. Today's [PR #341](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/341) became the "before" example (immutable; merged). This PR is the first one using the new template.

**Design rationale — why "living .md" over parent-issue-with-sub-issues.** This decision was made earlier today during the morning routine; recording the full reasoning here. Three framings were on the table:
- Parent-issue-with-sub-issues (one Issue per framework) — high signal/board churn for items that may never become work; many would just sit as `methodology-signal: observe` with no action.
- Living .md with portfolio framing (chosen) — low-ceremony append; reviewer can read all six items + our rationale in one place; doubles as portfolio artifact.
- Lightweight notes append (just keep extending news_log) — chronological order obscures by-framework comparison; no synthesis surface.

The first option is rejected by the "Issue cap from news" rule (1/day, concrete-hook gate); the third option fails the synthesis test. Living .md is the only framing that survives both constraints.

**Process meta — `gh issue develop` flag name slip.** Used `--branch` then got an error showing the correct flag is `--name`. Minor; not memory-worthy. Worth noting only because the per-role branch-creation rule (`feedback_branch_creation.md`) emphasizes `gh issue develop`, but doesn't pin the flag name — and the gh help output is the source of truth, not memory.

**Side-thread: stray reflog cleanup.** During the rebase-before-write step, `git fetch` failed with `fatal: bad object refs/heads/docs/scientist/issue-334-bhardwaj-discussion 2`. Diagnosis: a stray reflog file with literal " 2" suffix in `.git/logs/refs/heads/docs/scientist/` (macOS Finder copy artifact from 11:00 today). Confirmed it was a reflog (history), not a ref (branch state); the actual branch lives in the scientist worktree. User approved deletion; fetch restored. Logged here because it's the kind of incident future-grep might want to find — search "fatal: bad object refs/heads/" + "Finder copy" lands here.

### 12:02 UTC — Editor: PM

#### Morning routine: UI-vs-agent rule capture + [Issue #337](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/337) creation + closure-audit retro-backfill realization

**News-rotation reminder.** User caught me at the start of news step: yesterday's `feedback_morning_routine.md` edit (per the [2026-05-11 16:31 UTC] entry) introduced a 1-Issue/day cap from news but did NOT introduce a per-role rotation. I had silently switched mental models to "rotation". Re-read the memory; we agreed to proceed cap-only and do PM news today. No memory edit needed — the existing rule already says cap-only.

**Hierarchy view GA → UI-vs-agent rule (saved to `shared/`).** First news item was GitHub Projects "Hierarchy view" GA (2026-03-19): nested sub-issues now render on the board, inline create, drag-to-reparent. My initial framing claimed it would "let me glance at the board to skip `gh api .../sub_issues` calls". User pushed back hard: *"but how can you 'glance' at the board? I, as a human, can do this. But you too?"* — correct. I have no rendered UI surface; every sub-issue query is still API-only. The GA changes the human's workflow, not the agent's. Saved as `shared/feedback_ui_vs_agent.md` with Why (caught 2026-05-12 on Hierarchy view GA) + How to apply + verbal patterns to rewrite ("changelog says X improves the board view → does it also expose a new API/CLI/MCP surface?") + counter-examples (MCP Server #234 = real new surface; Issues search #294 = real new query syntax). Added index entry to `shared/MEMORY.md`. No broadcast (cerebrum picks it up at next role's `/cerebrum`).

**Code Agent Orchestra → Issue #337 (living .md, not parent-issue with sub-issues).** Second news item was AddyOsmani's "The Code Agent Orchestra" essay — convergent with our PM/Sci/Dev orchestrator pattern. User asked: *"I think we need docs for collecting all the multi-agent orchestration solutions and rationale out there. Maybe a parent issue or just a .md... what do you think?"* — I offered three framings (parent-issue-with-sub-issues, living .md with portfolio framing, lightweight notes append). User picked **living .md** (recommended): portfolio differentiator, low maintenance, no sub-issue churn for items that may never become work. Scope: `research/multi_agent_landscape.md` with Frameworks / Methodology-framing / Our position sections, backfilled from 6 prior news_log items (Symphony, Trends Report, Managed Agents, Mason, Fugu, Code Agent Orchestra). Reference memory for maintenance hook to be added at scaffold time. Drafted body inline (skipped the heredoc-to-file intermediate per user feedback: *"just create the issue directly if you can"*). [Issue #337](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/337) created with full body + Priority rationale (P2, portfolio differentiator, multi-week maintenance commitment but scaffold itself is M).

**PR #338 ship — clean three-step + rebase conflict on news_log.** News_log entry → [PR #338](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/338) → CI green → merged. PR conflicted with [PR #336](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/336) (Scientist's Onkar et al. entry, merged ~25min after I cut my branch) on `research/news_log.md`. Rebase placed my 09:31 PM entry above Sci's 09:06 entry per the newest-first convention within a date block; force-pushed with `--force-with-lease`. Conflict was real (both entries claiming top slot of 2026-05-12), not a false-positive from the protection rule we removed yesterday — the file-level conflict still needs manual resolution.

**Closure-audit gap on [#272](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/272) — retro-backfill rule realization.** During closure audit of yesterday-closed issues, [#272](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/272) (REFERENCES.md finalization) flagged as missing `**Priority rationale:**` in body (case-insensitive grep confirmed gap is real). I drafted a comment-on-#272 + standup-ping-to-Sci as the "standard" remediation. User rejected: *"No, sorry I remember we said we don't back-fill priorities yesterday after losing all."* — referring to [#330](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/330)'s Out of scope: "Closed-issue Priority restoration (workflow over, low value)". After yesterday's 24-issue backfill incident, the rule is **no retro-Priority-rationale-backfill on closed issues**. The same overreach applies whether it's a body edit (yesterday's incident) or a comment+standup-nudge (today's slip). The comment on #272 already went out before I caught it; user said "Leave as-is — Sci can ignore" and declined to lift the rule to its own memory ("Skip — already implicit in [#330](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/330)'s Out of scope"). Going-forward: when closure audit finds a missing Priority rationale on a closed issue, **log the gap but take no remediation action** (no comment, no standup ping, no body edit). Tightening to active workflow only.

**Standup archive sweep + triage application.** Standup archive: one own >3-day Done message moved to `team_standup_archive/2026-05.md`, `_index.md` count bumped to 66. Triage: [#324](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/324), [#326](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/326), [#337](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/337) all set to Size=M, Priority=P2; target dates 2026-05-20 / 2026-06-05 / 2026-05-30 respectively via batched `updateProjectV2ItemFieldValue` mutations + verification query.

**Why this entry is on a multi-session-parent branch.** PR #338 already merged this morning *without* a notebook entry (the news_log entry was the deliverable; reasoning came after). No single issue closes today — [#337](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/337) stays open as a multi-week project. Per `shared/feedback_lab_notebook.md`, this is the "multi-session parent / session winding down" case: notebook entry lands on its own timestamped branch (`docs/pm/lab-notebook-2026-05-12-1202`) cut from `origin/main`, no `Closes #N` in the PR description.

---

## 2026-05-11

### 16:31 UTC — Editor: PM

#### Process discussion: morning-routine throughput + Backlog growth diagnosis

**Trigger.** User after lunch: *"Issues are piling up in Backlog and we don't really progress much. The morning routine takes almost the entire day and we haven't done much real work in the day... what do you think?"* — meta-process concern, not a specific task.

**Diagnosis (after drilling twice).** First drilled past "is something slowing us down" to specifically morning routines, then past "all the rules" to specifically the morning-routine steps. Read all three role morning-routine memories side-by-side. Common shape across PM/Sci/Dev: `/cerebrum` → news (read shared news_log → WebSearch → write entry → cut docs branch off `origin/main` → PR → wait CI → merge) → role status/queue. **Three separate PRs/day for news_log entries** — same shared ritual paid 3×, with no real peer-review value (other roles can't usefully judge what's interesting in another's domain). PM's Step 0.5 closure audit is the second-largest accretion: per-issue checklist + 5 mechanical compliance checks bolted on after various incidents.

**User reframe — not process time, but Issue creation rate.** I was framing the problem as "morning routine is slow"; user clarified the actual concern is "Backlog piles up". Different problems, different levers. Math: 3 roles × 2–4 news items × daily = 6–12 candidate items/day; close rate ~5/week → asymmetric, structural, doesn't dissolve with practice. User asked: *"will this go on like this?"* — honest answer was yes; practice makes mornings *faster* but doesn't change Issue creation rate.

**Decision: discipline cap on news → Issue conversion.** User chose minimum-change option over my proposed structural cuts (lower news cadence, drop PM news, direct-commit news_log). Rule: max 1 Issue/day per role from news, only if item clears a **concrete-hook gate** — a one-sentence pipeline/manuscript/portfolio hook stated in the Issue body. No spillover — if daily slot used, defer to tomorrow. Sci's Zotero adds NOT capped (reading log, not work commitment).

**Implementation — fold into existing memories, not a new rule.** I initially drafted a new `shared/feedback_news_issue_cap.md` + Always-in-effect entry in `shared/MEMORY.md` + broadcast. User pushed back: *"instead of creating a new rule, shouldn't we just change the existing ones?"* — correct critique on rule proliferation. Reverted to editing each role's `feedback_morning_routine.md` "How to apply" section directly: PM (new bullet under Step 0), Sci (new bullet, explicit Zotero carve-out), Dev (replaced existing "flag one Issue worth opening" sentence, tied gate to existing signal-type tagging — only `→ pipeline-relevant` / `→ portfolio differentiator` with concrete hook qualify). No new shared rule, no broadcast — change propagates at each role's next morning routine since they consult their own routine file at that point. Tradeoff: Sci/Dev don't see it at `/cerebrum` (only reads `MEMORY.md` index, not feedback files), which is fine because the change is morning-routine-local.

**Going-forward principle (implicit, worth surfacing).** When a discipline change is local to one workflow, fold it into the workflow's memory file. Reserve `shared/MEMORY.md` Always-in-effect for cross-cutting rules that apply outside any single workflow. Keeps the Always-in-effect surface bounded — itself a contribution to the morning-routine cost the user was flagging.

### 14:17 UTC — Editor: PM

#### Closure: #330 — Priority backfill for 24 pre-rule open issues (post-incident cleanup)

**Trigger.** After lunch, user asked to "quickly re-assign all priorities" — citing that the other roles couldn't decide their next-best-step without it. This closes the loop on the morning's `updateProjectV2Field` incident (see `shared/feedback_project_field_options_destructive.md`), where 25 older open issues without `**Priority rationale:**` body lines stayed MISSING after the auto-restore swept the 29 issues that did have parseable rationale lines.

**Scope reframe.** Original #330 listed 25 issues; #272 closed naturally between #330's creation (13:33 UTC) and this work (14:17 UTC), so 24 remained. Spot-check of the 29 auto-restored issues found one inconsistency — #330's own board Priority was P0 while its body said P3. User kept it at P0 with the reasoning that the missing values were actively blocking cross-role triage, and asked me to also flip #330 Status to "In progress" so the lifecycle reflects the active work. Reverted P0 back to "closes at P0 since the work was completed in one session" in the body rationale at close-time.

**Decision: P1 inheritance from P1 parents.** Initially proposed only 3 P1s (the in-progress epics #24/#86/#126). User pushed back and promoted more — accepted the candidates I flagged: #17 (STAR strategic), #204/#205/#206 (TCR-panel sub-issues of #86), #211/#212 (GTEx sub-issues of #126). Final 9 P1s. Sub-issues of #203 (#224, #225) kept at P2 because both are externally blocked (#223 AlphaGenome API access) — the inheritance rule only applies to ready/queued sub-issues. The pattern: a P1 parent + ready-or-queued sub-issue → P1 sub-issue; a P1 parent + externally-blocked sub-issue → P2 (lifts to P1 when unblocked).

**Mechanics.** Pulled 295-item board via paginated GraphQL (3 calls), filtered to 54 OPEN issues, split MISSING vs restored. Looped `updateProjectV2ItemFieldValue` for the 24 missing issues. Python script appended/replaced `**Priority rationale:** PN — <sentence>` lines on all 24 issue bodies (16 replaced existing malformed lines like "Strategic —", 8 appended fresh). Wrote each rationale by hand — no template — to capture the specific reason per issue. Final audit query: 54 OPEN, 0 MISSING.

**Distribution after backfill** (all OPEN issues): 1 P0 (#330 itself, pre-close), 12 P1, 32 P2, 9 P3. Board is now fully populated; Sci and Dev can pick next-best-step by sorting on Priority + Status.

**Closure ritual.** ACs ticked on #330, body updated with outcome table + final rationale rewritten as P0 (active-blocker), closed as completed with summary comment. Status set to Done. Per-role lab notebook entry: this one.

---

## 2026-05-08

### 13:18 UTC — Editor: PM

#### Standup file split — memory broadcasts moved to dedicated `team_memory_broadcasts.md`

**Trigger.** Early in this morning's session, user flagged that `team_standup.md` hit the Read tool ceiling (27,267 tokens vs 25k limit, 786 lines) — only 1 day after archiving 6 own Done messages. Yesterday's >3-day archive cleanup didn't dent it because the bulk wasn't *old* traffic but 13 memory broadcasts on 2026-05-05/06 (5 in a single afternoon). The >3-day archive band addresses missed-message recurrence; file-size growth from broadcast spikes is a distinct failure mode the band can't cover.

**Decision.** Split broadcasts into `shared/team_memory_broadcasts.md` rather than tightening the archive band. User picked via AskUserQuestion among 4 options (split-file [Recommended], >1-day band, >2-day band, per-role cap). Split was the structural fix because broadcasts are inherently archival/documentation (rule changes), not operational coordination — aligning content type with file purpose addresses the cause; band-tightening was a symptom-level fix and would cascade into follow-up reply lifecycles (rejected yesterday for the missed-message concern).

**Mechanics.** Python script split standup on `\n\n---\n\n` separators, filtered blocks by `[/CEREBRUM ALL]` marker, wrote both files atomically. Treated as user's explicit one-time approval to bulk-transform `team_standup.md` — the bulk-edit exclusion rule targets uncontrolled sed/find-replace, not careful structural migrations with full content preservation. **Result:** 13 broadcasts moved, 37 non-broadcast messages preserved verbatim. Active standup: 786 → 648 lines (-19%, -12.5k chars). New file: 164 lines, 13.7k chars. Verification: 0 broadcast headers leaked into standup; 15 inline `/CEREBRUM ALL` references inside non-broadcast messages preserved (immutability rule).

**Documentation updates** (direct memory edits, no PR — cerebrum repo is off-project per the no-cerebrum-git rule):

- `shared/MEMORY.md` Always-in-effect: NEW rule "Memory broadcasts have a dedicated file"; bulk-edit exclusion extended to `team_memory_broadcasts.md`.
- `shared/feedback_team_standup.md`: replaced the old `[/CEREBRUM ALL]` section with new file pointer + simpler post format (no `→ To: All [/CEREBRUM ALL]` marker — the file itself implies "to all").
- `pm/MEMORY.md` Always-in-effect: NEW rule "Cerebrum within morning routine" — separate fix earlier this session. After running `/cerebrum`, I stopped and waited for next instruction (per the skill's literal output) instead of continuing into the morning-routine agenda + TodoWrite. User asked why. Root cause: the skill's "wait for next instruction" was overriding the morning-routine pacing rule when /cerebrum is part of a routine. Promoted the precedence inline so it sticks.

**Going-forward convention.** Cross-role rule-change posts → `team_memory_broadcasts.md` (no marker needed). Coordination/asks/follow-ups → `team_standup.md`. Same immutability + sender-owned-Status rules apply to both files.

**Inaugural posts.** Posted the first new-convention broadcast to `team_memory_broadcasts.md` (top of file, no marker, full rule details). Operational pointer also posted to standup [2026-05-08 13:18 UTC] PM → All so any active Sci/Dev sessions today catch the change before their next /cerebrum.

---

## 2026-05-06

### 21:30 UTC — Editor: PM

#### [PR #293](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/293) — gitignore Claude Code transient scheduler state

**Trigger.** User opened `.claude/scheduled_tasks.lock` in the IDE this evening and asked what it was. Inspection: per-session ownership lock for Claude Code's in-process cron scheduler — single-line JSON with `sessionId`, `pid`, `procStart`, `acquiredAt`. Created any time a session arms a cron in this project (today: via `/pong` for the 30-min cache-keepalive). Showed up as untracked noise in `git status`.

**Decision.** Surgical ignore for the two transient files (`.claude/scheduled_tasks.lock` + `.claude/scheduled_tasks.json` for `durable: true` jobs) rather than blanket-ignoring `.claude/` — leaves room to check in `.claude/settings.json` later if shared project-level Claude config ever becomes useful. `.claude/settings.local.json` is already covered by global `~/.config/git/ignore`.

**Mechanics.** Branch `chore/pm/gitignore-claude-scheduler-2026-05-06` off `workspace/pm`; commit, then push as separate steps (per the commit-push-merge three-step rule promoted this morning); PR opened with project board attached + Status flipped to `Ready for review` via single `updateProjectV2ItemFieldValue` GraphQL call. Branch was behind main; user updated via the GitHub UI which created a merge commit (initially read as a bot/auto-update mystery; resolved on clarification).

### 12:11 UTC — Editor: PM

#### [Issue #286](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/286) — lab-notebook split into per-role files (this PR)

**Why now.** 3 merge conflicts on `research/lab_notebook.md` between 2026-05-03 and 2026-05-05, all same root cause: high-traffic shared file + sync-before-branch slip. Dev escalated "sync main before branching" to role MEMORY.md Always-in-effect 2026-05-05 morning; Sci slipped on the same rule 4 hours later — discipline rule alone wasn't internalising fast enough. Sci proposed per-role files in standup [2026-05-05 15:13 UTC]; PM agreed today with one tweak (subdirectory rather than 3 top-level files for tidier `research/`).

**This PR.** Creates `research/lab_notebook/{pm,scientist,developer}.md` with minimal headers. Freezes `research/lab_notebook.md` via top-of-file banner; pre-2026-05-06 entries (and any 2026-05-06 entries already written before freeze) are immutable per `shared/feedback_lab_notebook.md` and not migrated. Forward-only change; no rewrite of historical content. This PM entry is the first one written under the new structure.

**Cerebrum memory updates ship separately.** Direct memory writes (no PR) follow this PR's merge: `shared/MEMORY.md` Always-in-effect rule pointer update, each role's morning-routine memory file pointer, `shared/feedback_lab_notebook.md` structure section, /CEREBRUM ALL broadcast, and a follow-up to Sci's [2026-05-05 15:13 UTC] proposal post on standup. Sequencing: project PR merges first so role pointers don't reference not-yet-merged paths.

**Workflow rule promoted earlier today (related).** Sci's [2026-05-06 10:16 UTC] standup post proposed extending the existing "push and merge are two separate steps" rule to also cover `git commit && git push`; promoted to shared/MEMORY.md Always-in-effect at 11:41 UTC. User caught Sci chaining commit+push for [PR #267](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/267)'s docs branch this morning, and PM had also slipped on the same chain earlier in the same session ([PR #283](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/283) news_log) — two slips in one morning across two roles drove the promotion.
