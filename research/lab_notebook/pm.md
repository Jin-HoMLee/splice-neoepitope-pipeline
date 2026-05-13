# Lab Notebook — PM

Per-role lab notebook for PM sessions. Started 2026-05-06 as part of the lab-notebook split (see `lab_notebook.md` freeze note for rationale).

Format and rules unchanged from the unified notebook — see `shared/feedback_lab_notebook.md`. New date sections at the TOP; new time sections at the TOP of the date block; entries are immutable once committed.

---

## 2026-05-13

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
