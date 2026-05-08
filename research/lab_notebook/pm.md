# Lab Notebook — PM

Per-role lab notebook for PM sessions. Started 2026-05-06 as part of the lab-notebook split (see `lab_notebook.md` freeze note for rationale).

Format and rules unchanged from the unified notebook — see `shared/feedback_lab_notebook.md`. New date sections at the TOP; new time sections at the TOP of the date block; entries are immutable once committed.

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
