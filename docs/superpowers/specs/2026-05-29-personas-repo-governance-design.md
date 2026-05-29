# Personas-Repo Governance — Edit Boundaries + MM-Owned Commit Lifecycle

> **Status:** approved design (2026-05-29) · **Issue:** [#567](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/567) · **Extends:** [#527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) (Memory Manager rollout)

## Problem

The `claude-personas-splice-neoepitope-pipeline` repo (the memory repo backing the PM / Scientist / Developer / Memory-Manager workflow on `splice-neoepitope-pipeline`) has no enforced governance for *who may edit what* and *who lands changes*. The result, observed 2026-05-29: morning-routine memory edits from multiple role sessions piled up **uncommitted** on a checked-out feature branch (`infra/add-to-project-workflow`), behind `main`, at risk of loss — 10 modified + 4 untracked files with no clear committer. Recovery took a full backup → stash → fast-forward → re-apply → commit cycle.

Root cause is a **governance interregnum**:

- The only *live* rule (`shared/MEMORY.md`, 2026-05-26) is **"Personas-repo git state is not your responsibility"**: roles may edit memory files in-session, but the commit/push lifecycle happens *"outside the session, managed by the user."* In practice nobody committed, so edits stranded.
- The **Memory Manager (MM)** model from [#527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) — roles edit, MM commits — was *designed* (Subs 1–2) but **not onboarded** (Subs 3–9 unshipped; `memory_manager/` dir, personas `CLAUDE.md`, and `.claude/settings.json` do not exist).

So neither model is actually operating, and any role session edits anything with no committer of record.

## Goals

- One coherent, **end-to-end** policy: who may edit which areas **and** who commits/pushes.
- **Prevent stranding** — uncommitted edits must surface immediately, not accumulate invisibly.
- Realize it **without blocking on the stalled MM onboarding**, while still committing to the MM model long-term.
- Prefer **mechanism over memory** where a rule has already proven fragile (this gap recurred despite documented rules).

## Non-goals

- Re-litigating whether MM should exist (decided: yes — it is the committer).
- Memory *content* quality / slimming (tracked separately under the memory-slim epic [#538](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/538)).
- Project-repo (`splice-neoepitope-pipeline`) git workflow — unchanged.

## Decisions (brainstormed + approved 2026-05-29)

1. **Scope:** end-to-end (edit boundaries AND commit/push lifecycle).
2. **Commit owner:** dedicated Memory Manager — roles edit, MM commits/pushes.
3. **`shared/` edits:** MM-curated — roles *propose*, MM is the one who edits + commits `shared/`.
4. **Realization:** phased (rules + scan now; mechanism enforcement once MM onboards).

## §1 — Steady-state model (once MM is onboarded)

| Zone | Edit rights | Non-owner change request | Commit / push |
|---|---|---|---|
| Own `<role>/` dir | that role edits freely | n/a | **MM only** |
| `shared/` (rules all sessions load) | **MM only** | role *proposes* via role-labeled issue / standup note → MM edits | **MM only** |
| Another role's `<role>/` dir | nobody | file a request (issue/standup), e.g. the [#12 dead-link handoff](https://github.com/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/issues/12) | **MM only** |
| `team_standup.md` + archives | append-only (existing rule, unchanged) | n/a | **MM only** |

**Key consequence (approved):** a role's only *direct* write is to its **own dir**. A role **cannot directly edit `shared/`** even inline — it proposes and MM makes the edit. This is the strict reading of decision 3 and the largest day-to-day friction point; accepted deliberately so cross-cutting rules get a second set of eyes (MM's) before they reach every session.

**Own-dir hand-off:** after editing its own dir, a role leaves the existing **"Memory edits — for MM to commit"** bullet in its lab-notebook/standup entry, listing touched files, so the next MM session lands them.

**Session model assumption (confirmed):** PM / Scientist / Developer / MM each run as *distinct* Claude sessions (separate `claude` invocations, distinct cwd / role identity). Phase 2's per-role permission split depends on this.

## §2 — Phase 1: interim (now — MM not yet onboarded)

1. **Land the policy as rules.** Rewrite the live `shared/MEMORY.md` rule *"Personas-repo git state is not your responsibility"* into the §1 model. Because this is itself a `shared/` edit and MM does not yet exist, **PM-as-caretaker makes this edit under a documented bootstrap exception** (this design authorizes it).
2. **Anti-stranding mechanism that works today.** A mandatory **session-start `git -C <personas> status` scan** in each role's morning routine — surfaces uncommitted / stranded state on the first message instead of letting it accumulate. This *is* [#527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) **Sub 6**; it needs no MM and directly prevents the 2026-05-29 incident class.
3. **Interim committer = PM-as-caretaker**, committing/pushing in-session **with the user's push gate**, replacing the failed "user commits outside the session" assumption. Caretaker work is journaled in `pm.md` tagged `(MM caretaker)` per the existing convention.

## §3 — Phase 2: enforcement (after [#527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) Subs 3–9 onboard MM)

1. **Per-role `.claude/settings.json` on the personas repo** — `git commit` / `git push` (and personas-write `git -C`) **denied** in PM / Sci / Dev sessions, **allowed** only in the MM session. Mechanically enforces "MM commits."
2. **PreToolUse hook** — blocks `Write` / `Edit` whose target path is under `shared/**` or another role's `<role>/**` **unless the session is MM**. Each role keeps free write to its own dir; MM may write anywhere. Mechanically enforces the §1 edit boundaries.
3. **PM-caretaker commit role ends**; MM assumes it. The git-status scan (§2.2) stays permanently.

## §4 — Relationship to [#527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527)

This design **refines and extends** the MM rollout:

- §2.2's git-status scan **is** #527 Sub 6 — landing it here satisfies that sub.
- §3 **adds two mechanisms #527 never had.** The original #527 design relied on memory rules alone for "MM commits" + edit boundaries — exactly the adherence model that just proved fragile. Phase 2 makes them deterministic.
- §1 **tightens** #527's edit model: the original was "active roles may edit personas memory (MM commits)"; this restricts `shared/` to MM-only.
- Phase 2 is **gated on MM onboarding** (#527 Subs 3–9), so this Issue and #527 are explicitly coupled.

## Rule changes (Phase 1 deliverables)

| File | Change |
|---|---|
| `shared/MEMORY.md` | Replace Always-in-effect "Personas-repo git state is not your responsibility" with the §1 model (own-dir-edit; `shared/` + cross-role = propose to MM; MM commits; interim = PM-caretaker commits) + the git-status-scan rule. |
| `pm/`, `scientist/`, `developer/feedback_morning_routine.md` | Add the session-start `git -C <personas> status` scan step. |
| (Phase 2) personas `.claude/settings.json` (per-role), new `.claude/hooks/` PreToolUse edit-boundary guard | Created during MM onboarding. |

## Open questions / risks

- **`shared/` friction in the interim.** Until MM is onboarded, `shared/` changes route through PM-caretaker (standing in for MM). This is acceptable but means cross-cutting rule changes are PM-gated in the window.
- **Phase 2 hook path-matching.** The edit-boundary hook must reliably know "the session's role." Mechanism: derive from cwd (each role's clone path) or an explicit env/setting — to be settled in the implementation plan.
- **Coupling risk.** If MM onboarding (#527) stalls indefinitely, Phase 2 never lands and we run on Phase 1 (rules + scan + PM-caretaker) permanently. That is a *safe* degraded state (the scan prevents stranding), not a failure — but full edit-boundary enforcement waits.

## Sequencing

1. **Phase 1** (this Issue): rule rewrite + git-status scan + PM-caretaker commit convention. Independent; lands now.
2. **Phase 2** (folded into #527 Subs 3–9 + new mechanism sub-tasks): permission split + edit-boundary hook. Gated on MM onboarding.
