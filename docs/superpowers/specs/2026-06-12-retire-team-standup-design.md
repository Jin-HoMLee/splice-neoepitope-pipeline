# Retire `team_standup.md` → board + Discussions — Design

**Issue:** [#569](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/569) (`role:pm`, `arc:board-governance`, `arc-phase:active`)
**Date:** 2026-06-12
**Author:** PM
**Status:** Design — pending implementation plan

## Problem

`team_standup.md` (a flat-markdown async message board in the personas repo's `shared/`) has been the inter-role coordination medium. A 2026-05-30 PM data check showed coordination **already migrated** off it on its own:

- Standup posting collapsed ~10× (posts/week: W18 `33` → W19 `50` → W20 `5` → W21 `16` → W22 `4` partial).
- Issue creation held steady-to-rising (created/week: W19 `27` → W21 `31` → W22 `52`; 30 Issues saw comment activity in the trailing 7 days).

The board now **supersedes** the file for *work* coordination (Status · Priority · Size · Target · milestone, native sub-issues + `blockedBy`, `role:*` labels, an immutable event timeline, `gh`-addressable). The only residual job the file does that the board can't is host **non-work-item async messages** — FYIs, heads-ups, re-raises not tied to any Issue. That is the GitHub **Discussions** job.

Lineage: the multi-role → multi-agent transition ladder ([Issue #265](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/265)) — this is the reshaped Rung 2→3 step (no markdown parser needed once the board *is* the bus).

## Goal

Retire `team_standup.md`. Route **work coordination** through Issues + board; route **non-work-item messages** through a GitHub Discussions "Team Coordination" category. Update the memory rules + morning routine that bake in the standup-file convention.

## Design

### 1. The routing rule

Every message that today goes to the standup file routes by a single test: **is it about a specific tracked work item?**

| Message kind | Destination | Mechanism |
|---|---|---|
| Tied to a specific Issue (status, handoff, "decide X here", blocker, review ask) | **The board** | A comment on that Issue, `From:/To:` attribution; metadata via Status / Priority / Size / `blockedBy` / `role:*` |
| Not tied to any Issue (FYI, heads-up, re-raise, thinking-out-loud, team notice) | **Discussions → Team Coordination** | A Discussion, light template (§2) |

**Default when unsure:** if you can name the Issue it's about, it's a board comment. The standup file's job split cleanly along this seam — the board already absorbed work-coordination (the ~10× posting collapse); Discussions absorbs the Issue-less residue.

There is no GitHub-native cross-role @mention (all roles are the same GitHub user). Delivery is the **morning-routine scan** (§4), not a notification — `From:/To:` attribution + the board's `role:*` label carry addressing.

### 2. Discussions "Team Coordination" convention

The category already exists (`DIC_kwDORwn9EM4C-Jo6`, created 2026-05-30 — P1 done).

- **One Discussion per topic**; title prefix `[FYI]` / `[REQUEST]` / `[DECISION]`.
- Body opens `**From:** <Role> → **To:** <Role(s)>` (carries the standup attribution).
- Replies = **threaded comments**.
- **Open = live / needs attention; close = resolved** — the open/closed state *is* the status; no `**Status:**` field.
- **The raiser closes** the thread once satisfied (or the addressee closes with a confirming reply) — GitHub norm.

**Prefix taxonomy** (web-checked 2026-06-12): tags sort on the proven **action-expectation axis** (what's expected of the reader), so a reader triages without opening the message — one consistent axis, not message-genre:

| Prefix | Means | Discussions example |
|---|---|---|
| `[FYI]` | informational, no response needed (absorbs "heads-up" — synonyms) | "deleted orphaned `workspace/<role>` branches" |
| `[REQUEST]` | asks for input/something, not yet a tracked Issue, less time-bound | "anyone object to dropping the X env before I file it?" |
| `[DECISION]` | a team/process call **not tied to a single Issue** | "should we adopt convention X across roles?" |

Two deliberate choices: (a) **fold "heads-up" into `[FYI]`** and drop "re-raise" as a prefix — a re-raise is a *posting-reason* (reopen the Discussion / note in body), not an action-type; the underlying message is still FYI/REQUEST/DECISION. (b) **`[ACTION]` is omitted on purpose** — anything requiring action on *tracked work* routes to the **board** as an Issue comment, so an `[ACTION]` tag in Discussions would be a smell ("this should've been an Issue"); its absence *reinforces* the routing seam (§1).

Rationale for the light template overall: lightweight-structured beats both free-form (an everything-channel buries signal) and strict schema (200-person ceremony that re-imports the manual bookkeeping the migration sheds); proportional to a 4-role team. If drift appears, harden into a real `.github/DISCUSSION_TEMPLATE/team-coordination.yml` later (YAGNI until then).

### 3. Never-amend → graduated immutability (native)

The standup file's hard "never amend a post, always follow up" rule existed because flat-markdown edits **erased history untraceably** — anyone who read the old version was silently desynced. GitHub removes that failure mode: every comment/Discussion edit carries an "edited" marker + full per-edit history (user + timestamp) viewable by anyone with read access.

So the rule relaxes to the platform norm — a **threshold, not a prohibition**:

- **Minor correction** (typo, formatting, a fix moments after posting) → **edit is fine**; the "edited" marker preserves transparency and keeps the thread readable.
- **Substantive change** (alters meaning, adds material context, or lands after someone may have acted) → **post a follow-up reply**, don't silently rewrite.

Work side: the Issue **event timeline is append-only** regardless. Discussions side: the same edit-small / reply-big threshold.

### 4. Morning-routine rewrite (4 files)

The routine's "read `team_standup.md`" step becomes a **two-channel check**:

- **`shared/feedback_morning_routine.md`** (Phase 2, shared mechanics) — replace the file-read with: (a) board scan (already present — `gh issue list` by `role:*` + Status + recent comment activity), (b) open Discussions in Team Coordination (`gh api .../discussions` filtered to the category; open = needs attention).
- **The 3 role files** (`pm/`, `scientist/`, `developer/` `feedback_morning_routine.md`) — each names the standup-file read in its phase list; swap for the same two-channel check, role-scoped.

No new "since last session" anchor is needed — open/closed state carries it.

### 5. Step −1 scope-by-actionability split

The whole-repo `git status` scan every role runs at Step −1 currently surfaces *all* uncommitted personas-repo state — but PM/Sci/Dev can't commit, so cross-role stranded edits re-surface as recurring noise (the `feedback_hook_scope_by_actionability.md` anti-pattern). Split by **who can act**:

- **MM** = full-repo scan **+ drain** (the stranding-catch belongs with the sole committer).
- **PM / Sci / Dev** = flag only **own-dir edits made this session** — no whole-repo scan.
- **Suppress-already-flagged** — don't re-surface an edit already handed off in an open thread (same shape as the [Issue #639](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/639) News-repetition fix).
- Note: the Discussions migration *shrinks* this scan anyway — coordination posts stop touching the working tree, leaving only real memory edits.

### 6. Retirement sequencing (the live-handoffs problem)

The file can't be deleted today: at design time it holds ~7 live Pending handoffs (including this session's own governance-change handoff). Ordered, inside the single MM-committed PR off #569:

1. **P3 lands** — the rewritten rules go live (board + Discussions are now canonical).
2. **Drain the live messages (P2 generalized)** — each open Pending post routed to its home: work-items → comments on their Issue; Issue-less FYIs → Discussions. The 3 originally-named messages (2026-05-29/30) are stale (~2 weeks) — verify Done/obsolete, don't re-migrate noise.
3. **P4 retire** — `git mv team_standup.md` + `team_standup_archive/` → `_retired/` (move, **never delete** — archive rule); leave a pointer `README.md` ("coordination moved to board + Discussions; see #569").
4. **P5 carve-off** — board-query dispatch digest → a separate follow-up Issue (gated on P3 landing). Plus the out-of-scope spin-off: a small "actually wire `blockedBy`" nudge Issue (the dependency graph is currently unused — 0 open `is:blocked` at pivot).

### 7. Deliverable & ownership

- **This session's output:** this design spec (pipeline repo, `docs/superpowers/specs/`), on branch `design/pm/issue-569-retire-standup`.
- **Execution:** the P3 rewrites (§4, §5) + P4 file move (§6) touch the **personas repo's** `shared/` — they land via **one PR**, authored against this spec, **MM-committed** (MM is sole committer of the personas repo since 2026-06-03, personas#23; cross-role `shared/` + all-roles morning routine = exactly MM's review surface).
- **Repo split:** the spec doc + #569 itself live in the **pipeline repo** (normal PM commit flow). The memory rewrites + file retirement live in the **personas repo** (MM-committed).
- **P1:** tick done — the category already exists.

## Files affected

**Pipeline repo (this branch):**
- `docs/superpowers/specs/2026-06-12-retire-team-standup-design.md` (this file)

**Personas repo (MM-committed, separate PR per §7):**
- `shared/MEMORY.md` — retire standup-file Always-in-effect rules (file order, never-amend-file, archive-Done-3d); replace with board-as-bus + Discussions conventions
- `shared/feedback_team_standup.md` → rewrite/retire as `feedback_team_coordination.md`
- `shared/feedback_morning_routine.md` — Phase 2 two-channel check + Step −1 actionability split
- `pm/`, `scientist/`, `developer/` `feedback_morning_routine.md` — role-scoped two-channel check
- `shared/team_standup.md` + `team_standup_archive/` → `_retired/` + pointer `README.md`

## Out of scope

- **P5 dispatch digest** and the **`blockedBy`-wiring nudge** are carve-off follow-up Issues, not part of this migration.
- Hardening the Discussions convention into a `.github/DISCUSSION_TEMPLATE/` form — deferred until drift is observed (YAGNI).

## Connects to

- [Issue #265](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/265) — transition ladder (reshaped Rung 2→3 step)
- [Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) — Memory Manager role owns the governed memory commits P3 requires
- `shared/feedback_personas_governance.md` — governs the P3 `shared/` edits
