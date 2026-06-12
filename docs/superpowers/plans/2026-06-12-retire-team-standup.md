# Retire `team_standup.md` → board + Discussions — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Retire the flat-markdown `team_standup.md` async-message file; route work coordination through board (project 9) Issue comments and Issue-less async messages through a GitHub Discussions "Team Coordination" category; rewrite every memory rule + morning-routine step that bakes in the standup-file convention.

**Architecture:** Three edit surfaces, committed separately. (1) **Pipeline repo** on this branch `design/pm/issue-569-retire-standup` (normal PM flow): the design spec [already landed], this plan, and the `.claude/commands/standup.md` command rewrite — this PR `Closes #569`. (2) **Personas repo** `claude-personas-splice-neoepitope-pipeline` (MM-committed per `shared/feedback_personas_governance.md`): all `shared/` + role-memory rewrites, the morning-routine two-channel check, the incidental-prose sweep, and the file retirement to `_retired/`. (3) **GitHub live actions**: drain the 14 live Pending messages to their channel and file the two carve-off follow-up Issues. `#569` closes only after all three are done (the pipeline PR merges last).

**Tech Stack:** Markdown memory files; `gh api graphql` for Discussions; `gh issue comment` / `gh issue create` for board routing; `git mv` for the archival retirement.

---

## Orientation — read before touching anything

**Two repos, two clones:**
- **Pipeline repo** = the project repo. Your cwd `splice-neoepitope-pipeline-pm` is a clone of it. Branch `design/pm/issue-569-retire-standup` is already checked out with the spec committed (`e5d0d17`). `#569` lives here (project board 9).
- **Personas repo** = `/Users/jin-holee/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline`. Holds all the memory files (`shared/`, `pm/`, `scientist/`, `developer/`, `memory_manager/`). Currently clean on `main`. Its Issues are `claude-personas#N` — **distinct** from pipeline `#569`.

**Commit model (`shared/feedback_personas_governance.md`, self-edit model adopted 2026-06-12):**
- Any role **edits** `shared/` + own-role dir directly in the personas working tree. **MM is the sole committer/pusher + reviewer.** So: this plan's executor *makes* the personas-repo edits; **MM commits them** via a board-tracked PR.
- The normal handoff bullet ("Memory edits — for MM to commit") posts to standup — **but standup is what we're retiring.** So this migration's handoff routes to a **board comment on `#569`** instead (work-item → board, per the very rule we're landing). That comment must carry the **before→after** for every *contradict/remove* edit (tier-2 "louder" flag — we are removing existing `shared/` rules, not just adding).

**Cross-repo `#569` mechanics:**
- The personas-repo PR **cannot** `Closes #569` (cross-repo keywords don't auto-close). Reference it in prose as `Jin-HoMLee/splice-neoepitope-pipeline#569` with **neutral phrasing** ("part of …", not "closes …") per `shared/feedback_hash_numbers.md`.
- The **pipeline-repo PR** (this branch) carries `Closes #569` and merges **last**, after the personas PR is merged + messages drained.

**Source-of-truth facts captured at planning time (2026-06-12):**
- Discussions category **Team Coordination** = `DIC_kwDORwn9EM4C-Jo6`, slug `team-coordination`; repo id `R_kgDORwn9EA`. P1 (category exists) is **done**.
- One pre-existing open Discussion in the category: `#571` (2026-05-30, Sci-lane carve question) — a seed post; verify/close during the drain (Task 11), don't treat as live.
- The active `team_standup.md` holds **14 live Pending** messages spanning 2026-06-04 → 2026-06-12 (51 message headers total).
- Full standup-reference map: **49 files**. Bucket A (4 core protocol), B (5 morning-routine), C (27 incidental prose — split C1 live-instructional / C2 historical), D (15 frozen — **exclude**), E (3 role MEMORY).

**Verification idiom (this plan's substitute for TDD):** each edit task ends with a **grep gate** — a command that must return the expected post-edit state. "Red" = the stale string still present; "green" = gone / replaced. Run the gate after the edit; if red, fix before moving on.

---

## File Structure

**Pipeline repo (this branch — PM-committed):**
- `docs/superpowers/specs/2026-06-12-retire-team-standup-design.md` — the spec (already landed)
- `docs/superpowers/plans/2026-06-12-retire-team-standup.md` — this plan
- `.claude/commands/standup.md` — `/standup` command, rewritten to scan the two new channels

**Personas repo (MM-committed PR):**
- `shared/feedback_team_coordination.md` — **NEW** keystone protocol (replaces `feedback_team_standup.md`)
- `shared/feedback_team_standup.md` — **removed** (`git rm`; content superseded by the new file)
- `shared/feedback_standup_two_halves.md` — **removed** (its mechanic — re-raise Pending / archive Done — vanishes under open=live/close=resolved)
- `shared/feedback_standby_trigger.md` — rewritten (stand-by now watches board + Discussions)
- `shared/MEMORY.md` — retire 4 standup Always-in-effect rules + repoint 3 index links
- `shared/feedback_morning_routine.md` — Phase 2 two-channel check + Step −1 actionability split
- `pm/feedback_morning_routine.md`, `scientist/feedback_morning_routine.md`, `scientist/feedback_morning_routine_format.md`, `developer/feedback_morning_routine.md` — role-scoped two-channel check
- `pm/MEMORY.md`, `scientist/MEMORY.md`, `developer/MEMORY.md` — update standup-derived inline rules / index links
- `shared/feedback_personas_governance.md` — handoff mechanic (the one Bucket-C file kept; rest of the sweep is **carved** to a follow-up Issue — Task 9 banner + Task 13 Step 3)
- `shared/team_standup.md` + `shared/team_standup_archive/` → `shared/_retired/` (move, after drain)
- `shared/_retired/README.md` — **NEW** pointer ("coordination moved to board + Discussions")

**GitHub (live):**
- 14 board comments / Discussions (the drained messages)
- 2 carve-off Issues (dispatch digest; `blockedBy`-wiring nudge)

---

## Phase 0 — The keystone: new coordination protocol

### Task 1: Author `shared/feedback_team_coordination.md`

Everything downstream points here, so write it first.

**Files:**
- Create: `<personas>/shared/feedback_team_coordination.md`

- [ ] **Step 1: Write the new protocol file**

```markdown
---
name: team-coordination-board-and-discussions
description: How roles coordinate after the standup file's retirement (#569) — work-item messages go to board Issue comments; Issue-less async messages go to the Team Coordination Discussions category. Routing rule + Discussions convention + graduated immutability.
metadata:
  type: feedback
---

`team_standup.md` (a flat-markdown async message board) was retired 2026-06-12 ([Jin-HoMLee/splice-neoepitope-pipeline#569](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/569)). Coordination now runs on two GitHub-native channels. This file is the protocol; design rationale in `docs/superpowers/specs/2026-06-12-retire-team-standup-design.md` (project repo).

## The routing rule — one test

**Is the message about a specific tracked work item?**

| Message kind | Destination | Mechanism |
|---|---|---|
| Tied to a specific Issue (status, handoff, "decide X here", blocker, review ask) | **The board** | A comment on that Issue. Body opens `**From:** <Role> → **To:** <Role(s)>`; metadata via Status / Priority / Size / `blockedBy` / `role:*` labels |
| Not tied to any Issue (FYI, heads-up, re-raise, team notice, thinking-out-loud) | **Discussions → Team Coordination** | A Discussion, light template (below) |

**Default when unsure:** if you can name the Issue it's about, it's a board comment. There is no GitHub-native cross-role @mention (all roles are the same GitHub user) — delivery is the **morning-routine scan** (`feedback_morning_routine.md` Phase 2), not a push notification. `From:/To:` attribution + the `role:*` label carry addressing.

## Discussions "Team Coordination" convention

Category: **Team Coordination** (`DIC_kwDORwn9EM4C-Jo6`, slug `team-coordination`) on the pipeline repo.

- **One Discussion per topic**; title prefix `[FYI]` / `[REQUEST]` / `[DECISION]`.
- Body opens `**From:** <Role> → **To:** <Role(s)>` (carries the old standup attribution).
- Replies = **threaded comments**.
- **Open = live / needs attention; close = resolved.** The open/closed state *is* the status — no `**Status:**` field, no Pending/Done flips, no archive cadence.
- **The raiser closes** the thread once satisfied (or the addressee closes with a confirming reply).

**Prefix taxonomy** (action-expectation axis — the reader triages without opening):

| Prefix | Means | Example |
|---|---|---|
| `[FYI]` | informational, no response needed (absorbs "heads-up") | "deleted orphaned `workspace/<role>` branches" |
| `[REQUEST]` | asks for input/something, not yet a tracked Issue | "anyone object to dropping the X env before I file it?" |
| `[DECISION]` | a team/process call not tied to a single Issue | "should we adopt convention X across roles?" |

A **re-raise** is a posting-reason (reopen the Discussion / add a comment), not a prefix. `[ACTION]` is **omitted on purpose** — anything actionable on *tracked work* routes to the board as an Issue comment; an `[ACTION]` tag in Discussions would be a smell ("this should've been an Issue"). Harden into a `.github/DISCUSSION_TEMPLATE/team-coordination.yml` only if drift appears (YAGNI).

## Graduated immutability (native edit history)

The old "never amend a post" rule existed because flat-markdown edits erased history untraceably. GitHub removes that failure mode: every comment/Discussion edit carries an "edited" marker + full per-edit history viewable by anyone with read access. So the rule relaxes to a **threshold, not a prohibition**:

- **Minor correction** (typo, formatting, a fix moments after posting) → **edit is fine**; the "edited" marker preserves transparency.
- **Substantive change** (alters meaning, adds material context, or lands after someone may have acted) → **post a follow-up reply**, don't silently rewrite.

The Issue **event timeline is append-only** regardless. Same edit-small / reply-big threshold on the Discussions side.

## Reading patterns

- **Board:** `python3 scripts/board_open_items.py --role <role> --status "In progress"` (and `"Ready for review"`, `"In review"`), plus recent comment activity on your `role:*` Issues. (Project board 9 sorts Done-first — never an unpaginated `projectV2 { items(first: N) }` query; see `feedback_board_queries.md`.)
- **Discussions (open = needs attention):**

```bash
gh api graphql -f query='
{
  repository(owner: "Jin-HoMLee", name: "splice-neoepitope-pipeline") {
    discussions(first: 30, categoryId: "DIC_kwDORwn9EM4C-Jo6", states: [OPEN]) {
      totalCount
      nodes { number title createdAt author { login } }
    }
  }
}'
```

## Posting

- **Board comment:** `gh issue comment <N> --body "**From:** <Role> → **To:** <Role(s)>

<message>

**Created by:** <Role>"` (the `Created by:` footer per `feedback_github_workflow.md`; the `From:/To:` header already attributes, but the explicit footer keeps the artifact-attribution rule uniform).
- **New Discussion:** `gh api graphql` `createDiscussion` mutation with `repositoryId: "R_kgDORwn9EA"`, `categoryId: "DIC_kwDORwn9EM4C-Jo6"`, the `[PREFIX]` title, and the `From:/To:` body. Reply via `addDiscussionComment`; resolve via `closeDiscussion`.
```

- [ ] **Step 2: Grep gate**

Run: `grep -c "Open = live" <personas>/shared/feedback_team_coordination.md`
Expected: `1` (file written, convention present).

- [ ] **Step 3: Stage (do not commit — MM commits in Task 15)**

The file is on disk; it is now readable by every role's next `/memory`. No commit here.

---

## Phase 1 — Retire / rewrite the core protocol files (Bucket A)

### Task 2: Remove the superseded protocol files

**Files:**
- Remove: `<personas>/shared/feedback_team_standup.md`
- Remove: `<personas>/shared/feedback_standup_two_halves.md`

- [ ] **Step 1: Confirm the new file fully absorbs both**

`feedback_team_standup.md` (format, immutability, archive cadence, reading patterns) is superseded by `feedback_team_coordination.md`. `feedback_standup_two_halves.md` (re-raise Pending + archive Done — the dual hygiene sweep) describes a mechanic that **no longer exists** (open/closed state replaces Pending/Done; no archive). Both retire.

- [ ] **Step 2: Remove the files**

```bash
git -C <personas> rm shared/feedback_team_standup.md shared/feedback_standup_two_halves.md
```

(MM runs the actual `git rm` at commit time; the executor may `rm` to keep the working tree consistent for grep gates, and flag it. Either way the files leave the tree.)

- [ ] **Step 3: Grep gate — find every inbound link to repoint**

Run: `grep -rn "feedback_team_standup.md\|feedback_standup_two_halves.md" <personas> --include="*.md" | grep -v "_retired\|team_standup_archive\|\.silent\|drafts/"`
Expected: a list of live references (in `MEMORY.md` files, role memories, `feedback_personas_governance.md`, etc.). **Every hit is a repoint target** for Tasks 5–6 and the Task 10 sweep — note them. The gate is "green" only after those tasks land and this grep returns no live (non-frozen) hits.

### Task 3: Rewrite `shared/feedback_standby_trigger.md`

**Files:**
- Modify: `<personas>/shared/feedback_standby_trigger.md`

- [ ] **Step 1: Repoint the watch target from the standup file to the two channels**

The trigger phrases ("keep an eye out", "stay on stand-by", "babysit the channel") and the `/loop 4m /standup` mechanism survive — but `/standup` now scans **board pings + open Discussions** (Task 14), not the retired file. Edit the file so:
- The title/description and body say the stand-by watches **board comments addressed to your role + open Team Coordination Discussions**, via `/loop 4m /standup`.
- Remove "babysit the standup / watch the standup / monitor for new messages [in the file]" framing; keep "babysit the channel / the inbox" (now = board + Discussions).
- Keep the multi-cron collision caveat verbatim (still true).
- Replace the line `**Why:** User established this pattern 2026-05-05.` tail so it references the channel, not the file.

- [ ] **Step 2: Grep gate**

Run: `grep -in "team_standup\|the standup file\|watch the standup" <personas>/shared/feedback_standby_trigger.md`
Expected: no hits (all file-specific framing gone; `/standup` command name may remain since the command still exists, rewritten).

---

## Phase 2 — MEMORY.md rules (Bucket A core + Bucket E role)

### Task 4: Rewrite the standup Always-in-effect rules in `shared/MEMORY.md`

**Files:**
- Modify: `<personas>/shared/MEMORY.md` (lines ~17, 31, 32, 34 = rules; 53, 57, 58 = index links)

- [ ] **Step 1: Replace the 3 standup-mechanic Always-in-effect rules with one board/Discussions rule**

Remove these three bullets (the flat-file mechanics that no longer exist):
- "**Standup file order:** New messages … go at the TOP of `team_standup.md` …"
- "**Never amend any standup post:** … immutable … bulk-edit must exclude `team_standup.md` …"
- "**Archive standup messages, never delete:** … `team_standup_archive/<YYYY-MM>.md` …"

Replace with a single Always-in-effect bullet:

```markdown
- **Team coordination — board comment or Discussion, never a file:** Work-item messages (tied to a specific Issue) → a comment on that Issue, body opening `**From:** <Role> → **To:** <Role(s)>`. Issue-less async messages (FYI / heads-up / re-raise / team notice) → a GitHub **Discussion** in **Team Coordination** (`DIC_kwDORwn9EM4C-Jo6`), title-prefixed `[FYI]` / `[REQUEST]` / `[DECISION]`. **Open = live / needs attention; close = resolved** — no Pending/Done flips, no archive cadence. Edit-small / reply-big (native edit history makes minor in-place edits safe; substantive changes get a follow-up reply). Full rule: `feedback_team_coordination.md`. (Replaced the retired `team_standup.md` file convention — [#569](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/569), 2026-06-12.)
```

- [ ] **Step 2: Keep the stand-by-trigger bullet but repoint it**

The "Stand-by trigger phrases" bullet (line ~34) stays — edit its tail so `/loop 4m /standup` watches "board pings + open Team Coordination Discussions" and the `src` is `feedback_standby_trigger.md` (unchanged target, rewritten content).

- [ ] **Step 3: Repoint the compact-trigger bullet**

The "durable-capture first" bullet (line ~35) lists "board / lab-notebook / standup" as capture sinks. Change "standup" → "board comment / Discussion".

- [ ] **Step 4: Fix the index links (Shared section, lines ~53, 57, 58)**

- Remove the `[Team standup](team_standup.md)` index line.
- Remove the `[Standup two halves](feedback_standup_two_halves.md)` index line.
- Replace `[Team standup protocol](feedback_team_standup.md)` → `[Team coordination](feedback_team_coordination.md) — board comments + Team Coordination Discussions; routing rule + graduated immutability`.

- [ ] **Step 5: Repoint the bare-path-convention example (line ~15)**

That long bullet ends with "Standup posts: write `team_standup.md` (same dir as shared/MEMORY.md). Past immutable content (archived standup posts …) stays as-is regardless of path form." Change the live instruction to "Coordination posts: board comment or Discussion (no file). Past immutable content (the retired `team_standup.md` + archives under `_retired/`, retired broadcasts under `.silent/.remote/`) stays as-is regardless of path form."

- [ ] **Step 6: Grep gate**

Run: `grep -in "team_standup\|standup post\|standup message\|standup file" <personas>/shared/MEMORY.md`
Expected: hits only where they name the **retired-file as history** (e.g. the `_retired/` path mention). No live "post to standup" / "amend standup" / index-link instruction remains.

### Task 5: Update role MEMORY.md standup-derived rules (Bucket E)

**Files:**
- Modify: `<personas>/shared/feedback_personas_governance.md` (~lines 16, 17, 23, 26, 30 — the coupled exception kept from the carved sweep)
- Modify: `<personas>/scientist/MEMORY.md` (line ~22)
- Modify: `<personas>/developer/MEMORY.md` (lines ~17, 19 — verify)
- Modify: `<personas>/pm/MEMORY.md` (line ~37)
- Modify: `<personas>/pm/project_memory_md_slimming.md` (lines ~63, 73, 112)

- [ ] **Step 1b: `feedback_personas_governance.md` — the handoff mechanic (coupled exception)**

This is the one Bucket-C file kept in the main PR (the rest of the sweep is carved — Task 9 banner). The migration's own handoff (Task 15) routes the "Memory edits — for MM to commit" flag to a **board comment**, which only makes sense if the governance rule says so. Edit:
- The **edit-commit matrix** row `| \`team_standup.md\` + archives | append-only (all roles) | MM only |` → change to `| \`_retired/\` standup archives | frozen history — no edits/appends | MM only |` (the file is retired; nothing is appended anymore).
- The **handoff instruction** "add a **'Memory edits — for MM to commit'** bullet to your standup / lab-notebook entry" → "post a **'Memory edits — for MM to commit'** comment on the tracking Issue (board) / add it to your lab-notebook entry" (the board is the new handoff sink; standup is gone).
- "file a request (issue / standup)" (the "another role's memory" row, ~line 16/30) → "file a request (issue / `[REQUEST]` Discussion)".
- Leave the historical incident references intact (e.g. the 2026-05-29 stranding lineage).

- [ ] **Step 1: scientist/MEMORY.md — the "flip your post" inline rule**

The bullet *"Standup follow-up replies must end with the 'flip your post' prompt …"* describes the Pending→Done close-prompt, which no longer exists (Discussions/Issues close natively). **Remove the bullet entirely** (its mechanic is gone); if a residual coordination-etiquette note is wanted, replace with: *"Coordination replies: close the Discussion (or the board thread) when satisfied; the raiser owns the close. Full rule: `shared/feedback_team_coordination.md`."*

- [ ] **Step 2: developer/MEMORY.md — standup-archive exemption + attribution**

- The "**Standup-archive PRs are exempt from lab notebook entries**" bullet: the standup-archive PR flow is going away (no more archiving). Edit to note the exemption is historical and the active rule is "coordination lives on the board + Discussions, not in PR'd files." Keep the broadcasts/news-log retirement refs (still valid history).
- The "**`Created by:` on every GitHub artifact**" bullet mentions "Standup follow-ups already carry attribution via the `### [...] From: <Role>` header — exempt." Change to "Discussion posts + board coordination comments carry `**From:** <Role>` in the body — still add the `Created by:` footer per the uniform rule."

- [ ] **Step 3: pm/MEMORY.md + pm/project_memory_md_slimming.md**

- `pm/MEMORY.md` line ~37: the closure-audit index entry mentions "false-positive standup pings". Change "standup pings" → "board-comment / Discussion pings".
- `pm/project_memory_md_slimming.md`: this is an **audit-findings doc** describing MEMORY.md slimming. Its standup mentions (line ~63 "never edit `team_standup.md`", ~73 "auto-fire `/loop 4m /standup`", ~112 "`PreToolUse` Edit matcher on `team_standup.md`") are **describing the rules as they were at audit time**. Treat as **C2 historical** — leave unless a mention is a live instruction. (The PreToolUse-matcher note at ~112 is a forward proposal; add a one-line "superseded by #569 retirement" annotation rather than rewriting.)

- [ ] **Step 4: Grep gate**

Run: `grep -in "standup" <personas>/pm/MEMORY.md <personas>/scientist/MEMORY.md <personas>/developer/MEMORY.md`
Expected: no live-instruction hits (only historical annotations, if any).

---

## Phase 3 — Morning routine (Bucket B)

### Task 6: `shared/feedback_morning_routine.md` — two-channel check + Step −1 split

**Files:**
- Modify: `<personas>/shared/feedback_morning_routine.md` (Phase 2 status block ~lines 52–60; Step −1 personas-scan; the "Read before claiming" line ~84)

- [ ] **Step 1: Replace the standup read-step (lines ~53)**

Replace the bullet `- \`team_standup.md\` — messages \`To: <your role>\` (pending + new since last check)` with a **two-channel** block:

```markdown
- **Board pings** — comments `**To:** <your role>` on your `role:*` Issues + recent comment activity (`python3 scripts/board_open_items.py --role <role> --status "In progress"`, then `"Ready for review"`, `"In review"`).
- **Open Discussions** — the Team Coordination category (open = needs attention):
  `gh api graphql -f query='{ repository(owner:"Jin-HoMLee", name:"splice-neoepitope-pipeline"){ discussions(first:30, categoryId:"DIC_kwDORwn9EM4C-Jo6", states:[OPEN]){ totalCount nodes{ number title author{login} } } } }'`
```

- [ ] **Step 2: Delete the standup-hygiene two-halves block (lines ~55–60)**

The "Standup hygiene (two halves — run both …)" block (re-raise Pending >1d, archive Done >3d, bulk-edit exclusions) describes mechanics that no longer exist. **Delete the whole block.** Replace with one line:

```markdown
**Coordination hygiene:** close any Discussion / board thread of yours that's resolved (open = live). No archive cadence — GitHub retains closed threads.
```

- [ ] **Step 3: Step −1 scope-by-actionability split**

Find the Step −1 personas-repo `git status` scan. Add the actionability split (spec §5):

```markdown
**Step −1 personas-repo scan — scope by who can act:**
- **MM** = full-repo `git -C <personas> status` scan **+ drain** (the stranding-catch belongs with the sole committer).
- **PM / Sci / Dev** = flag only **own-dir + `shared/` edits made this session** — not a whole-repo sweep (you can't commit others' stranded edits; re-surfacing them is the `feedback_hook_scope_by_actionability.md` anti-pattern).
- **Suppress-already-flagged** — don't re-surface an edit already handed off in an open board thread / Discussion (same shape as the #639 News-repetition fix).
```

- [ ] **Step 4: Fix the "Read before claiming" example (line ~84)**

`"no Pending" / "no broadcasts" / "standup is empty"` → `"no open Discussions" / "no board pings"`. Update the `Read team_standup.md` instruction to the two-channel queries above.

- [ ] **Step 5: Grep gate**

Run: `grep -in "team_standup\|standup hygiene\|Pending >1\|Done >3" <personas>/shared/feedback_morning_routine.md`
Expected: no hits.

### Task 7: `pm/feedback_morning_routine.md` — PM-scoped two-channel check

**Files:**
- Modify: `<personas>/pm/feedback_morning_routine.md` (lines ~8, 64, 73, 84, 99, 120, 123, 127, 187)

- [ ] **Step 1: Swap every live "standup ping/ask/pointer/handoff" → board-comment / Discussion**

Apply the substitution rubric (Task 10) to each PM-routine reference, all of which are **live instructions**:
- ~8 "standup hygiene shape" → "coordination-hygiene shape (close resolved threads)"
- ~64 "post a standup ask going forward" → "post a board comment / Discussion ask going forward"
- ~73 "Always pair the comment with a standup pointer" → "the board comment **is** the durable+discoverable record now (closed Issues still surface via `board_open_items.py --role`); for a fully Issue-less nudge use a Discussion" — **note:** this rule's original point (closed issues drop off Dev/Sci workflows, so the standup pointer adds discoverability) partly dissolves because the morning routine's board scan already covers role-labelled items. Rewrite to: "leave the comment on the Issue (durable); if the Issue is closed and the nudge is time-sensitive, open a `[REQUEST]` Discussion."
- ~84 "assignment + standup ping" → "assignment + board comment"
- ~99 "queue for next standup" → "queue as a board comment / `[REQUEST]` Discussion"
- ~120, ~123 "push a standup handoff ping" → "post a board comment handoff ping (on the Issue)"
- ~127 Phase 2 "standup messages To: PM + hygiene two-halves" → "board pings + open Discussions To: PM; close resolved threads"
- ~187 "1 standup item to flag" (top-line summary example) → "1 coordination item to flag"

- [ ] **Step 2: Grep gate**

Run: `grep -in "standup" <personas>/pm/feedback_morning_routine.md`
Expected: no hits.

### Task 8: Scientist + Developer morning-routine files

**Files:**
- Modify: `<personas>/scientist/feedback_morning_routine.md`
- Modify: `<personas>/scientist/feedback_morning_routine_format.md` (lines ~ the 2 standup refs)
- Modify: `<personas>/developer/feedback_morning_routine.md`

- [ ] **Step 1: Apply the same two-channel swap, role-scoped**

For each file, replace any "read standup / standup To: <role> / standup hygiene" step with the two-channel check (board pings `To: <role>` + open Discussions) and the "close resolved threads" hygiene line. Mirror Task 7's rubric.

- [ ] **Step 2: Grep gate**

Run: `grep -in "standup" <personas>/scientist/feedback_morning_routine.md <personas>/scientist/feedback_morning_routine_format.md <personas>/developer/feedback_morning_routine.md`
Expected: no live-instruction hits.

---

## Phase 4 — Incidental-prose sweep (Bucket C)

### Task 9: Sweep the C1 live-instructional references

> **⚠️ CARVED to a follow-up Issue (decision 2026-06-12).** This whole task is deferred out of the main migration PR to keep its review surface tight — filed as the 3rd carve-off Issue in Task 13 Step 3, which links back to this section for the file list + rubric. **One exception stays in the main PR:** `feedback_personas_governance.md` (handled in Task 5 Step 1b below), because the migration's own handoff (Task 15) depends on its routing being updated — leaving it would make the migration contradict its own governance. During the deferred window the corpus carries dangling "post to standup" references in the ~13 remaining C1 files; that's the accepted tradeoff, tracked by the sweep Issue. **Subagents executing this plan: SKIP Task 9** — it is reference material for the follow-up Issue, not main-PR work.

**Files (C1 — live instructions; ~14 files):**
- `<personas>/shared/feedback_best_next_issue.md` (~41, 45, 47, 53 — "Check standup messages")
- `<personas>/shared/feedback_milestone_closure_routing.md` (~23, 25, 37, 38, 42, 44, 88–92 — routing-decision via standup)
- `<personas>/shared/feedback_personas_governance.md` (~16, 17, 23, 26, 30 — "add a bullet to your standup", "file a request (issue / standup)", the `team_standup.md` row in the matrix)
- `<personas>/shared/feedback_compact_trigger.md` (~19, 50 — capture sinks "board / lab-notebook / standup")
- `<personas>/shared/feedback_board_hygiene.md` (~92 — "summary report posted to standup")
- `<personas>/shared/feedback_pm_skeletal_only.md` (~22 — "surface in standup")
- `<personas>/shared/feedback_scope_discipline.md` (~9, 60 — "post a standup heads-up to PM")
- `<personas>/shared/feedback_refetch_before_mutate.md` (~11 — "standup nudge to clean up")
- `<personas>/shared/feedback_mechanism_over_memory.md` (~18 — "scan recent standup nudges")
- `<personas>/shared/feedback_ask_user_question.md` (~28 — "What to do about this PM standup message?" example)
- `<personas>/shared/feedback_team_structure.md` (~25, 32, 37, 152 — coordination-channel refs)
- `<personas>/shared/feedback_multi_role_not_multi_agent.md` (~13 — "human-readable artifacts (standup file, sub-issues)")
- `<personas>/shared/feedback_outcome_routing.md` (~16 — standup nudge as a routing channel; the ~104 incident narrative is C2, leave)
- `<personas>/pm/feedback_ask_for_help.md` (~69, 71 — "log the decision in board recap or standup")
- `<personas>/pm/feedback_milestones.md` (~208 — "post a `[MILESTONE ACHIEVED]` … in `team_standup.md`")

- [ ] **Step 1: Apply the substitution rubric, file by file**

**Decision rule for each hit — is the sentence a live instruction or a historical record?**
- **Live instruction** (tells a *future* action): rewrite per the table below.
- **Historical / incident record** (records a *past* event — "caught on standup 2026-05-06", "yesterday's `team_standup.md` miss"): **leave unchanged** — it's a dated fact, not an instruction. (These are the C2 hits; do not over-edit history.)

| Old phrasing | New phrasing |
|---|---|
| "post a standup message / ask / heads-up / ping" | "post a board comment (on the Issue) / a `[FYI]`/`[REQUEST]` Discussion" |
| "standup nudge" | "board-comment nudge" |
| "check / read standup" | "check board pings + open Discussions" |
| "log to standup" | "record in a board comment / Discussion" |
| "surface in standup" | "surface on the board (Issue comment) or in a Discussion" |
| "board recap or standup" | "board recap or a Discussion" |
| "human-readable artifacts (standup file, sub-issues)" | "human-readable artifacts (board comments, Discussions, sub-issues)" |
| `[MILESTONE ACHIEVED]` "in `team_standup.md` addressed to all roles" | `[MILESTONE ACHIEVED]` "as an `[FYI]` Discussion (To: All) — or an Announcements-category post" |

**Two files need more than a phrase swap:**
- `feedback_milestone_closure_routing.md`: its whole routing-decision loop runs through standup ("surface via standup → respond on standup with the routing pick"). Rewrite the loop to use a **`[DECISION]` Discussion** as the async decision-log (the milestone-routing call is Issue-less by nature) + a board comment on each affected Issue as the notification. Keep the incident narrative (~88–92) as historical.
- `feedback_personas_governance.md`: the **edit-commit matrix** has a `team_standup.md + archives` row and the "add a bullet to your standup / lab-notebook entry" handoff. Update: the handoff bullet now posts to a **board comment on the tracking Issue** (or lab-notebook); replace the `team_standup.md` matrix row with a note that the file is retired (`_retired/`, append-only history). The "file a request (issue / standup)" → "file a request (issue / `[REQUEST]` Discussion)".

- [ ] **Step 2: Verify the C2 historical hits were left intact**

C2 (leave unchanged): `feedback_project_file_paths.md` (~15 "yesterday's `team_standup.md` miss"), `feedback_hash_numbers.md` (~84 "standup follow-ups in 2026-05"), `feedback_read_before_claiming.md` (~27, 38 — *these are live* "before claiming standup is empty, Read it" → actually C1: repoint to the two-channel claim; **reclassify to C1** and swap), `feedback_project_vs_meta.md` (~9, 11 — uses "the standup-file pattern" as a *meta-example of a project-specific convention*; leave, it's illustrative), `feedback_notebook_terminology.md` (~12, 21 — "don't say 'notebook' in standup messages" → swap "standup messages" → "board comments / Discussions"; C1), `feedback_portfolio_lens.md` (~26 — "discoverable async coordination via the standup file" → "via the board + Discussions"; C1), `feedback_retiring_rules_with_artifacts.md` (~15 — "the standup is the right channel" → "a board comment / Discussion"; C1), `feedback_milestone_closure_routing.md` incident (~88–92, leave), `feedback_github_workflow.md` (~176 — "Follow-up replies in standup … carry attribution via `From: <Role>`" → "Coordination comments / Discussion posts carry `From: <Role>`"; C1).

> **Note on reclassification:** several files the categorization first marked ambiguous resolve to **C1** on the live-vs-historical test above. When in doubt, apply the test sentence-by-sentence, not file-by-file.

- [ ] **Step 3: Grep gate — no live "standup" instruction survives**

Run:
```bash
grep -rin "standup" <personas> --include="*.md" \
  | grep -v "_retired/\|team_standup_archive/\|\.silent/\|/drafts/\|docs/superpowers/plans/\|docs/superpowers/specs/" \
  | grep -vi "retired\|historical\|2026-05\|caught\|incident\|miss\b"
```
Expected: **empty** (every remaining `standup` mention is either in frozen history, or annotated as a retired/historical fact). Any hit is a missed live instruction — fix it.

---

## Phase 5 — Drain the live Pending messages (P2)

### Task 10: Route the 14 live Pending messages to their channel

**Files:**
- Read: `<personas>/shared/team_standup.md` (the active file)
- Live actions: `gh issue comment` / `gh api graphql createDiscussion`

- [ ] **Step 1: Enumerate the live Pending set**

```bash
grep -nB1 -A10 "Status:\*\* Pending" <personas>/shared/team_standup.md
```
Expect 14 messages dated 2026-06-04 → 2026-06-12 (From: MM/Sci/Dev/PM). For each, read the full message body.

- [ ] **Step 2: Triage each by the routing rule (§1)**

For each Pending message:
- **Names a specific Issue/PR** (status, handoff, blocker, review ask) → **resolve or re-home on the board**: if already actioned, no migration needed (it was waiting on a Pending→Done flip that no longer exists — just confirm Done in your notes). If still open, post the content as a **comment on that Issue** (`**From:**/**To:**` header + `Created by:` footer).
- **Issue-less** (FYI / team notice / cross-role question with no Issue) → open a **Discussion** (`[FYI]`/`[REQUEST]`/`[DECISION]`).
- **Stale** (the 2026-06-04/05 messages may already be resolved): verify against fresh `gh issue view` state before re-homing — **don't migrate noise**. If resolved, drop it.

- [ ] **Step 3: Verify the pre-existing seed Discussion #571**

`gh api graphql` read Discussion #571 (2026-05-30 Sci-lane carve). If its question is long-resolved, **close it** (`closeDiscussion`) with a confirming comment. It predates the convention — bring it into open=live/close=resolved compliance.

- [ ] **Step 4: Gate — no live Pending left to strand**

After routing, every former Pending message is either (a) resolved (its ask satisfied — confirmed via `gh issue view`), or (b) re-homed to a board comment / Discussion. Record the disposition of all 14 in the `#569` handoff comment (Task 15) so MM can verify before the file moves.

---

## Phase 6 — Retire the file (P4)

### Task 11: Move `team_standup.md` + archive to `_retired/`

**Files:**
- Move: `<personas>/shared/team_standup.md` → `<personas>/shared/_retired/team_standup.md`
- Move: `<personas>/shared/team_standup_archive/` → `<personas>/shared/_retired/team_standup_archive/`
- Create: `<personas>/shared/_retired/README.md`

- [ ] **Step 1: Confirm the drain is complete (Task 10 gate green) before moving**

Do **not** move the file until every live Pending message is dispositioned. The move is the point of no return for in-tree discoverability.

- [ ] **Step 2: `git mv` (move, never delete — archive rule)**

```bash
mkdir -p <personas>/shared/_retired
git -C <personas> mv shared/team_standup.md shared/_retired/team_standup.md
git -C <personas> mv shared/team_standup_archive shared/_retired/team_standup_archive
```
(MM runs these at commit time; executor mirrors in the working tree for grep gates.)

- [ ] **Step 3: Write the pointer README**

Create `<personas>/shared/_retired/README.md`:

```markdown
# Retired — team_standup.md

`team_standup.md` (flat-markdown async message board) and its monthly archive were **retired 2026-06-12** ([Jin-HoMLee/splice-neoepitope-pipeline#569](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/569)).

**Coordination moved to two GitHub-native channels:**
- **Work-item messages** → a comment on the relevant Issue (project board 9).
- **Issue-less async messages** → a GitHub **Discussion** in the **Team Coordination** category (`[FYI]` / `[REQUEST]` / `[DECISION]`).

Protocol: `shared/feedback_team_coordination.md`. Design: `docs/superpowers/specs/2026-06-12-retire-team-standup-design.md` (project repo).

The files here are **frozen history** — grep-able for "did we discuss X?", never edited or appended. The old `.silent/.remote/` memory-broadcast archive is a separate retirement ([#482](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/482)).
```

- [ ] **Step 4: Grep gate — no live path reference to the old location**

Run: `grep -rn "shared/team_standup.md\|team_standup_archive" <personas> --include="*.md" | grep -v "_retired/"`
Expected: no hits **except** historical/annotated ones (frozen drafts, plans). Live code/instructions must point to `_retired/` or to the new channels.

---

## Phase 7 — `/standup` command + carve-off Issues + verification

### Task 12: Rewrite the `/standup` command (pipeline repo, this branch)

**Files:**
- Modify: `splice-neoepitope-pipeline-pm/.claude/commands/standup.md`

- [ ] **Step 1: Replace the one-liner**

Current: `Apply the [Team standup protocol](feedback_team_standup.md) — Full standup rules: confirm with user, re-raise own Pending >1 day, archive own Done >3 days (move, never delete)`

New:

```markdown
Scan the two team-coordination channels and surface anything addressed to or relevant to your role; confirm with the user before acting on any request.

1. **Board pings** — comments `**To:** <your role>` on your `role:*` Issues + recent comment activity: `python3 scripts/board_open_items.py --role <role> --status "In progress"` (then `"Ready for review"`, `"In review"`).
2. **Open Discussions** — Team Coordination category (open = needs attention):
   `gh api graphql -f query='{ repository(owner:"Jin-HoMLee", name:"splice-neoepitope-pipeline"){ discussions(first:30, categoryId:"DIC_kwDORwn9EM4C-Jo6", states:[OPEN]){ totalCount nodes{ number title author{login} } } } }'`

Bring any pending item to the user, summarize what it asks, and act only after confirmation. Close a Discussion / board thread when its ask is satisfied (open = live). Full protocol: `shared/feedback_team_coordination.md`.
```

- [ ] **Step 2: Grep gate**

Run: `grep -in "feedback_team_standup\|Pending >1\|archive own Done" splice-neoepitope-pipeline-pm/.claude/commands/standup.md`
Expected: no hits.

- [ ] **Step 3: Commit (pipeline repo — PM flow, this is *not* MM-gated)**

```bash
git -C splice-neoepitope-pipeline-pm add .claude/commands/standup.md
git -C splice-neoepitope-pipeline-pm commit -m "feat(coordination): /standup scans board + Discussions; retire standup-file (#569)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```
(The plan doc + spec are already committed on this branch from the planning session; this commit adds only the command.)

### Task 13: File the three carve-off Issues (P5)

- [ ] **Step 1: Dispatch-digest follow-up Issue**

`gh issue create` (pipeline repo) for the board-query dispatch digest (the periodic "here's what's open per role" summary the standup used to provide informally). Gate it on this migration landing. Include `role:pm`, a `**Priority rationale:**` line, and reference `Jin-HoMLee/splice-neoepitope-pipeline#569` in prose (neutral, non-closing). Apply `arc:board-governance`.

- [ ] **Step 2: `blockedBy`-wiring nudge Issue**

`gh issue create` for the out-of-scope spin-off: actually wire `blockedBy` on the dependency graph (currently unused — 0 open `is:blocked` at pivot). `role:pm`, priority rationale, `arc:board-governance`.

- [ ] **Step 3: Incidental-prose sweep Issue (the carved Phase 4)**

`gh issue create` for the deferred Bucket-C1 sweep (~13 files still referencing the retired standup as a live channel). Body: link this plan's **Task 9** for the exact file list + the substitution rubric + the C1/C2 live-vs-historical test; note `feedback_personas_governance.md` was already handled in the main PR. `role:pm`, a `**Priority rationale:**` line (P2 — corpus-consistency cleanup, non-blocking; the new protocol file already exists so the dangling refs have a valid target), `arc:board-governance`. Reference `Jin-HoMLee/splice-neoepitope-pipeline#569` in prose (neutral, non-closing).

- [ ] **Step 4: Note all three Issue numbers in the `#569` handoff comment (Task 15).**

### Task 14: Final cross-surface verification

- [ ] **Step 1: Corpus consistency gate (personas repo) — scoped to the main-PR surface**

Run the Task 9 Step 3 grep across the personas repo. Because the incidental sweep is **carved**, the expected output is **not empty** — it should contain only the ~13 deferred Bucket-C1 files (tracked by the Task 13 Step 3 Issue). **There must be ZERO hits** in the core/morning-routine/governance/MEMORY surface — i.e. `feedback_team_coordination.md`, `feedback_standby_trigger.md`, `feedback_personas_governance.md`, `shared/MEMORY.md`, all `feedback_morning_routine*.md`, and the role `MEMORY.md` files. Confirm the remaining hits are exactly the deferred set and nothing from the main-PR files.

- [ ] **Step 2: Link-resolution gate**

Run: `grep -rn "feedback_team_standup.md\|feedback_standup_two_halves.md" <personas> --include="*.md" | grep -v "_retired/\|team_standup_archive/\|/drafts/\|docs/superpowers/"`
Expected: empty (every inbound link repointed to `feedback_team_coordination.md` or removed).

- [ ] **Step 3: New-file presence gate**

Run: `test -f <personas>/shared/feedback_team_coordination.md && test -f <personas>/shared/_retired/README.md && test -f <personas>/shared/_retired/team_standup.md && echo OK`
Expected: `OK`.

- [ ] **Step 4: Discussions channel smoke test**

Run the open-Discussions query (Task 1 Step 2 of the protocol). Expected: returns the drained messages (if any were Issue-less) + confirms #571 is closed.

### Task 15: Hand off to MM for the personas-repo commit

- [ ] **Step 1: Post the "Memory edits — for MM to commit" handoff as a board comment on `#569`**

Because standup (the usual handoff sink) is being retired, the handoff routes to a **board comment on `#569`**. The comment must list **every** touched personas-repo file and carry the **before→after** for each *contradict/remove* edit (tier-2 "louder" flag per `feedback_personas_governance.md` — we removed `feedback_team_standup.md`, `feedback_standup_two_halves.md`, and 3 MEMORY.md rules). Include the disposition of all 14 drained Pending messages (Task 10 Step 4) and the 2 carve-off Issue numbers.

```bash
gh issue comment 569 --body "**From:** PM → **To:** Memory Manager

**Memory edits — for MM to commit** (personas repo, #569 standup retirement):

<file list + before→after for the removes + drained-message dispositions + carve-off Issue #s>

**Created by:** PM"
```

- [ ] **Step 2: MM creates the personas-repo PR**

Per `feedback_personas_governance.md` (board-tracked Issue → `gh issue develop` → PR → merge), MM stages the specific files (never `git add -A`), commits, opens a PR, and — since this is a routine-adjacent governance change touching all roles' memory — requests an `@claude review` before merge (the routine-PR-bot-review pilot rule). The PR body references `Jin-HoMLee/splice-neoepitope-pipeline#569` in **neutral** phrasing (cross-repo, non-closing).

- [ ] **Step 3: After the personas PR merges + messages drained, merge the pipeline PR**

The pipeline-repo PR (spec + plan + `/standup` command) carries `Closes #569` and merges **last** via `bash scripts/audit_and_merge.sh <PR>` (closure-ritual gate). Its merge closes `#569` and drives board auto-Done.

---

## Self-Review (run by the plan author)

**Spec coverage (spec §1–§7):**
- §1 routing rule → Task 1 (protocol table) ✓
- §2 Discussions convention → Task 1 ✓
- §3 graduated immutability → Task 1 + Task 4 Step 1 ✓
- §4 morning-routine rewrite (4 files) → Tasks 6, 7, 8 ✓
- §5 Step −1 actionability split → Task 6 Step 3 ✓
- §6 retirement sequencing (P3 land → drain → P4 move → P5 carve) → Phases 1–4 (land) → Phase 5 (drain) → Phase 6 (move) → Task 13 (carve); ordering enforced by Task 11 Step 1 gate ✓
- §7 deliverable & ownership → Orientation + Task 12 (pipeline) + Task 15 (personas/MM) ✓

**Coverage beyond the spec (discovered during planning — flagged):**
- The spec's "Files affected" named the core files but not the **27 incidental-prose files** — leaving them would dangle references to a deleted file. **Decision (2026-06-12): carved to a follow-up Issue** (Task 13 Step 3) to keep the main PR's review surface tight, with one coupled exception kept (`feedback_personas_governance.md`, Task 5 Step 1b) because the migration's own handoff depends on it. Task 9 retains the file list + rubric as reference material for the follow-up. Accepted tradeoff: a temporary window of dangling "post to standup" refs in the ~13 deferred files, tracked by the Issue.
- The spec under-specified the **`/standup` command** (pipeline repo) and **standby-trigger** (personas repo) — both break on file retirement. Added Task 12 + Task 3. This makes the migration **3-surface**, not the "single MM PR" the spec loosely implied — corrected in Orientation.
- Cross-repo `#569` closing mechanics (personas PR can't `Closes`) — corrected in Orientation + Task 15.

**Placeholder scan:** the only intentional `<...>` placeholders are `<personas>` (the personas-repo abs path, defined in Orientation), `<role>`, and the Task 15 comment body (filled at execution from the actual file list). No TODO/TBD.

**Type/name consistency:** category id `DIC_kwDORwn9EM4C-Jo6`, repo id `R_kgDORwn9EA`, new file `feedback_team_coordination.md`, retired dir `_retired/` — used identically throughout.

---

## Connects to

- [Issue #569](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/569) — this migration (pipeline board)
- `docs/superpowers/specs/2026-06-12-retire-team-standup-design.md` — the approved design
- `shared/feedback_personas_governance.md` — the self-edit / MM-commit model the personas PR follows
- [Issue #265](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/265) — multi-role→multi-agent transition ladder (reshaped Rung 2→3)
