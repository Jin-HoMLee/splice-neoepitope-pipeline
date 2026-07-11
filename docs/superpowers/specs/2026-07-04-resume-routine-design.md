# Light resume session-start routine - design

**Date:** 2026-07-04
**Author:** PM
**Status:** Draft (for review)
**Scope:** Shared (all roles: PM / Scientist / Developer / Memory Manager)

> **Amendment 2026-07-11 ([Issue #1114](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1114)).** The `/coordination` command this spec names (items 4 and 5.2, section 6) has been **renamed to `/inbox`** and its ping scan rewired board-wide; the old command's role-scoped scan was structurally blind to cross-role pings. Read every `/coordination` below as `/inbox`. The design itself is unchanged - only the command's name and the correctness of its step-4 scan. The body is left as written, per the convention that a dated spec is a point-in-time record rather than a living document.

## 1. Problem

A session opened with *"Good morning!"* reliably triggers the full morning routine - that habit works.
But a session opened with a *resume* greeting (*"Hey, I'm back"*, *"let's continue"*, *"I'm back, let's go"*) has no defined behavior.
The agent either flails (no clear routine) or over-applies the full morning ceremony - as happened on 2026-07-04 20:29 local, when *"Hey :) I'm back. Let's gooo!"* started rolling toward all five Kanban beats roughly three hours after the same day's morning cadence had already run.

Re-running the morning cadence on a same-day resume is wrong on two counts:

- It re-runs cadences that are canonically once-a-day-or-less (see premise below), which is redundant work.
- It overloads a fresh session's context with a wall of routine right when the user wants to get to work.

The gap is a missing, named, lightweight *re-orientation* routine for non-morning session starts.

## 2. Premise validation (web-grounded, 2026-07-04)

The split between a per-session re-orientation pass and the once-daily Kanban cadences is grounded on two axes:

- **Kanban cadence model.** The Daily Stand-up is the only canonically daily / per-session cadence (short, board-walk, blockers-first). Replenishment is weekly/biweekly and the Service Delivery Review biweekly-to-quarterly - explicitly not per-session. Re-running those on a resume is a cadence-frequency mismatch the method warns against. ([teachingagile.com](https://teachingagile.com/kanban/introduction/kanban-cadences), [djaa.com](https://djaa.com/kanban-cadences/))
- **Context-switching / re-orientation cost.** Resuming after stepping away requires the brain to pause, reorient, and recall where it left off (~23 min to fully refocus after an interruption). A compact recall step on resume is evidence-based, not ceremony. ([atlassian.com](https://www.atlassian.com/work-management/project-management/context-switching), [toggl.com](https://toggl.com/blog/context-switching))

Conclusion: keep the cadences to their once-daily rhythm (morning routine); add a per-session re-orientation routine (resume routine).

## 3. Design

### 3.1 Trigger - greeting classification (no marker, no state, no nudge)

The routine is selected purely by the opening greeting. No last-session marker, no automatic time-of-day detection, no "did the cadence run today?" nudge - the user's morning-greeting habit is already the reliable signal, and a nudge solves a problem that does not exist.

| Opening | Routine |
|---|---|
| Morning greeting - *"Good morning!"*, *"morning"*, *"gm"* | **Full morning routine** (unchanged) |
| Resume greeting - *"Hey, I'm back"*, *"let's continue"*, *"I'm back, let's go"*, *"let's gooo"*, *"back at it"*, or any session-opening greeting that is not a morning greeting and not a direct task request | **Light resume routine** (new) |
| Direct task request with no greeting - *"fix X"*, *"what's the status of #N"* | No routine render; run the silent session-start spine (Step -1) and go straight to the task |

Rule of thumb: **only an explicit morning greeting triggers the morning routine; any other session-opening greeting triggers the resume routine.**
Genuine ambiguity resolves to the *lighter* routine (resume) - never auto-escalate to the morning ceremony.
The user can always ask for the full routine explicitly (*"morning routine"* / *"good morning"*).

### 3.2 The resume routine - a five-item compact pass

Items 1-2 already fire at every session start under the shared spine (`shared/MEMORY.md` "Always run at every fresh session start"); the resume routine adds items 3-5 and packages all five as one named, compact pass.

1. **Memory check** - `/memory-check` (role + shared `MEMORY.md` re-read; surfaces any new Always-in-effect rule).
2. **Field recall** - read the latest `episodes/` entry + `post-it.md` (re-situate; prune done/stale post-its).
3. **What changed since I stepped away** - a two-to-three line board glance of recent closures / merges / moves. Not the full Service Delivery Review; a "what moved" surface, not a delivery-health audit.
4. **Anything addressed to me** - the `/coordination` scan (board pings `**To:** <role>` via `scan_addressed_comments.py` + open Team Coordination Discussions).
5. **Your next pull** - one line: the top item in flight (`In progress`) or the next `Ready` item, then hand the floor back.

### 3.3 Rendering - one compact message

Unlike the multi-beat morning ceremony (agenda in Message 1, then beat-by-beat with a task list), the resume routine renders as **one compact message**, no TodoWrite ceremony:

- Top line: brief greeting + the memory-check result folded in (per the morning-routine pacing rule, the memory line is the top line, never a standalone stop).
- Items 3-5 as short bullets (a couple of lines each at most).
- Close by handing the floor back with the next-pull as a suggestion, not an auto-start.

Board detail beyond the two-to-three line glance is provided only on request - the default stays light to honor the "do not overload the new session" constraint.

### 3.4 Explicitly out of scope (stays in the morning routine)

The resume routine does **not** run any of these - they remain once-daily in the morning routine:

The Daily Stand-up (named in section 2 as the only canonically per-session cadence) is deliberately absent from this list, and that is not an omission: the five-item pass above - specifically the what-changed glance, the `/coordination` scan, and the next-pull - *is* the compact per-session analog of the Stand-up. Only the heavier once-daily cadences below stay in the morning routine.

- Service Delivery Review (closure audit, milestone health, roadmap-overdue sweep, flow-health aging sweep).
- Replenishment (intake triage `No Status -> Backlog`; commitment `Backlog -> Ready`).
- Signals (external / methodology scan).
- Friday cleanup (branch hygiene + stale-Issue self-review).

If any of these genuinely need attention mid-day, they are pulled deliberately (the user asks, or a specific finding warrants it) - not re-run wholesale on every resume.

### 3.5 Relationship to the morning routine

The two routines are siblings that share the same Step -1 spine (memory check + field recall - authoritatively defined in `shared/MEMORY.md` under "Always run at every fresh session start"; other spec docs that paraphrase that spine more loosely defer to the memory definition).
The morning routine is the daily-cadence ceremony; the resume routine is the per-session re-orientation.
Both are memory-defined and greeting-triggered - consistent with how the morning routine is already delivered (no skill/command; the definition lives in memory and the agent executes it on the greeting).

## 4. Delivery - files changed

The routine is memory-wired, mirroring the morning routine. All personas-repo commits are the Memory Manager's (PM authors, MM commits).

1. **New shared memory file** `shared/feedback_resume_routine.md` - the canonical definition (this spec's section 3, in memory-file form: trigger table, five-item pass, rendering, out-of-scope).
2. **`shared/MEMORY.md`:**
   - Add an Always-in-effect bullet stating the greeting-classification rule (morning greeting -> morning routine; resume greeting -> resume routine).
   - Under "Always run at every fresh session start", note that a resume greeting continues the spine into the resume routine (items 3-5).
   - Add a Reference index line for `feedback_resume_routine.md`.
3. **Spec doc** (this file) lands in the project repo at `docs/superpowers/specs/` (PM commits the project-repo side via PR).

The greeting-classification rule and the resume-routine pointer live in **`shared/MEMORY.md` only** - not duplicated into each role's `MEMORY.md`.
Resuming is universal (every role does it the same way), shared `MEMORY.md` is loaded by all roles, and a shared-only home avoids editing other roles' dirs (against personas governance).
This is a deliberate simplification over an earlier draft that added a per-role pointer bullet.

No scripts, hooks, or new state files are introduced. The routine reuses existing tooling (`/memory-check`, `/coordination`, `scan_addressed_comments.py`, `board_open_items.py`).

## 5. Edge cases

- **Ambiguous greeting** (session-opening but not clearly morning) -> resume routine (lighter default); never auto-escalate.
- **No greeting, direct task** -> silent spine only, straight to the task (no routine render).
- **User explicitly asks for the morning routine mid-day** -> run it; the greeting-classification default is overridable by an explicit request.
- **First session of the day opened with a resume greeting** (user forgot to say "good morning") -> resume routine, no nudge. The user reliably says "good morning" when they want the cadence; honoring that means not second-guessing a non-morning greeting. (This is the deliberate no-nudge decision - accepted trade-off.)

## 6. Verification

- Dry-run the greeting classification against the trigger table (each row maps to the intended routine).
- Confirm the resume routine renders as one compact message with items 1-5 and no morning-cadence beats.
- Confirm the morning routine is unchanged (a *"good morning"* still runs all its cadence beats, which are day- and role-driven, not a fixed count).
- Confirm items 1-2 are not double-run when the spine already executed them in the same session start.

## 7. Follow-ups / rollout

- Tracking filed: parent Issue #1026 (process/flow work, milestone-free) covers the wiring + verify; spec sub-issue #1027 carries this doc and is closed by its PR.
- Spec doc PR (project repo) is not standalone - it closes the spec sub-issue #1027; the parent #1026 stays open until the memory wiring lands.
- Memory wiring (section 4 items 1-3) is committed by the Memory Manager on the personas repo.
- The wiring is shared-only (`shared/feedback_resume_routine.md` + `shared/MEMORY.md`); no per-role `MEMORY.md` edits, so no other-role dirs are touched.

**Created by:** PM
