# Personas-Repo Governance — Phase 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. **Caveat for THIS plan:** it edits the personas repo, and dispatching parallel subagents to edit it would re-create the "too many sessions touching it" problem — execute **inline** (executing-plans), in one PM-caretaker session.

**Goal:** Land Phase 1 of the personas-repo governance ([Issue #567](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/567), sub of [#527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527)): the access/lifecycle policy as memory rules, the session-start git-status scan (= #527 Sub 6), and the PM-caretaker interim commit convention.

**Architecture:** Pure memory-doc edits in the personas repo (`claude-personas-splice-neoepitope-pipeline`): one new shared feedback file holding the full policy, a concise Always-in-effect pointer replacing the obsolete rule in `shared/MEMORY.md`, and a session-start scan step added to each role's morning-routine file. Committed by PM-as-caretaker under the bootstrap exception (MM not yet onboarded), with the user's push gate. The project repo carries only the spec, this plan, and the lab-notebook entry on the #567 branch.

**Tech Stack:** Markdown memory files; `git -C` against the personas repo; `gh` for the project-repo PR. No code, no test harness — "tests" here are `grep` verifications that the right text is present/absent.

**Personas repo path (used throughout):** `~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline` (referenced below as `$PERSONAS`).

---

## File Structure

**Personas repo (`claude-personas-splice-neoepitope-pipeline`):**
- **Create** `shared/feedback_personas_governance.md` — the full governance policy (edit/commit matrix, propose-flows, commit ownership, interim caretaker, the "Memory edits — for MM to commit" flag, Phase-2 roadmap).
- **Modify** `shared/MEMORY.md` — replace the Always-in-effect "Personas-repo git state is not your responsibility" line with a concise new-model pointer; add the new file to the Reference index.
- **Modify** `pm/feedback_morning_routine.md` — add a session-start personas git-status scan step (PM = caretaker committer variant).
- **Modify** `scientist/feedback_morning_routine.md` — add the scan step (surface-and-flag variant).
- **Modify** `developer/feedback_morning_routine.md` — add the scan step (surface-and-flag variant).

**Project repo (`splice-neoepitope-pipeline`, branch `design/pm/issue-567-personas-repo-governance`):**
- **Modify** `research/lab_notebook/pm.md` — non-routine governance-session entry.
- Spec + this plan already committed on the branch.

**Decision baked in:** the policy lives in `shared/` only (auto-loaded by every role), **not** copied into each role's `MEMORY.md`. This diverges from #527 Sub 6's original "per-role MEMORY.md copy" step on purpose — per-role duplication fights the slimming epic ([#538](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/538)); one shared rule covers all roles.

---

## Task 1: Create the governance policy file

**Files:**
- Create: `$PERSONAS/shared/feedback_personas_governance.md`

- [ ] **Step 1: Write the file**

Write to `$PERSONAS/shared/feedback_personas_governance.md`:

```markdown
---
name: personas-repo-governance
description: Who may edit what in the personas memory repo and who commits/pushes. Own-dir edits free; shared/ + cross-role go through MM (propose); MM owns all commits. Interim: PM-caretaker commits with the user push gate + a session-start git-status scan.
metadata:
  type: feedback
---

Governance for the `claude-personas-splice-neoepitope-pipeline` memory repo. Decided 2026-05-29 ([Issue #567](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/567), sub of [Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527)); supersedes the 2026-05-26 "personas-repo git state is not your responsibility" rule. Full design: `docs/superpowers/specs/2026-05-29-personas-repo-governance-design.md` (project repo).

## Edit + commit matrix (steady state — once MM is onboarded)

| Zone | Who may edit | Non-owner change request | Commit / push |
|---|---|---|---|
| Own `<role>/` dir | that role, freely | n/a | MM only |
| `shared/` | MM only | role proposes (role-labeled issue / standup) -> MM edits | MM only |
| Another role's `<role>/` dir | nobody | file a request (issue / standup) | MM only |
| `team_standup.md` + archives | append-only (existing rule) | n/a | MM only |

**Your only direct write is your own `<role>/` dir.** Cross-cutting (`shared/`) and cross-role changes route through MM so they get a second set of eyes before reaching every session.

## Proposing a change you cannot make directly

- **`shared/` rule change** -> file a role-labeled issue OR a standup note describing it; MM makes the edit. (Interim: PM-caretaker makes it.)
- **Another role's memory** (dead link, correction, etc.) -> file a request (issue / standup); the owner role or MM handles it. Example: the dead-link handoff [claude-personas#12](https://github.com/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/issues/12).

## Own-dir edits: flag them for the committer

After editing your own `<role>/` memory in a session, add a **"Memory edits — for MM to commit"** bullet to your lab-notebook / standup entry listing the touched files. That bullet is the structural handoff so the next committer lands them.

## Commit / push ownership

- **Steady state:** MM is the sole committer/pusher for the entire personas repo.
- **Interim (MM not yet onboarded):** **PM-as-caretaker** commits/pushes, in-session, **with the user's push gate** — never chain commit -> push without the user OK. Caretaker work is journaled in `pm.md` tagged `(MM caretaker)`. This replaces the failed "user commits outside the session" assumption that produced the 2026-05-29 stranding incident.

## Anti-stranding: session-start scan

Every PM / Sci / Dev session runs a `git -C <personas-repo> status` scan at session start (see each role's `feedback_morning_routine.md`) and surfaces any uncommitted / stranded state immediately, so edits cannot silently pile up on a stale branch. This is [Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) Sub 6.

## Enforcement roadmap (Phase 2 — gated on MM onboarding)

Once MM is onboarded ([Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) Subs 3–9), two mechanisms make this deterministic (filed as #527 subs): (1) per-role `.claude/settings.json` denying `git commit`/`push` on personas to non-MM sessions; (2) a PreToolUse hook blocking `Write`/`Edit` on `shared/**` and other roles' dirs unless the session is MM. Until then, this file is the convention.
```

- [ ] **Step 2: Verify the file exists and is well-formed**

Run:
```bash
PERSONAS=~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline
test -f "$PERSONAS/shared/feedback_personas_governance.md" && head -4 "$PERSONAS/shared/feedback_personas_governance.md"
```
Expected: file exists; frontmatter `---` / `name: personas-repo-governance` printed.

- [ ] **Step 3: Verify no dead self-links**

Run:
```bash
grep -c "feedback_personas_governance" "$PERSONAS/shared/feedback_personas_governance.md"
```
Expected: `0` (the file does not reference itself).

---

## Task 2: Rewrite the obsolete Always-in-effect rule + index it

**Files:**
- Modify: `$PERSONAS/shared/MEMORY.md` (the "Personas-repo git state is not your responsibility" Always-in-effect line + the Reference index)

- [ ] **Step 1: Replace the obsolete rule line**

In `$PERSONAS/shared/MEMORY.md`, replace the entire Always-in-effect bullet that begins `- **Personas-repo git state is not your responsibility:**` (currently one long line) with:

```markdown
- **Personas-repo edits — own dir only; `shared/` + cross-role go through MM; MM commits:** Edit your own `<role>/` memory freely. **Never directly edit `shared/` or another role's dir** — *propose* instead (`shared/` rule -> role-labeled issue/standup; other-role -> file a request). All commits/pushes are **MM's**. **Interim (MM not yet onboarded):** PM-as-caretaker commits/pushes with the **user's push gate**, after a session-start `git -C <personas> status` scan; flag your own-dir edits with a "Memory edits — for MM to commit" bullet. Supersedes the 2026-05-26 "not your responsibility" rule (which produced the 2026-05-29 stranding incident). Full policy: `shared/feedback_personas_governance.md`; design: [Issue #567](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/567) / [Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527). <!-- src: shared/feedback_personas_governance.md -->
```

- [ ] **Step 2: Verify the old rule is gone and the new one is present**

Run:
```bash
PERSONAS=~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline
grep -c "not your responsibility" "$PERSONAS/shared/MEMORY.md"   # expect 0
grep -c "own dir only.*go through MM" "$PERSONAS/shared/MEMORY.md" # expect 1
```
Expected: first `0`, second `1`.

- [ ] **Step 3: Add the new file to the Reference index**

In `$PERSONAS/shared/MEMORY.md`, find the Reference-section list (the bullets near the bottom linking `feedback_*.md` files) and add:

```markdown
- [Personas-repo governance](feedback_personas_governance.md) — who edits what (own dir only; shared/ + cross-role via MM) + who commits (MM; interim PM-caretaker with user gate) + session-start git-status scan
```

- [ ] **Step 4: Verify the index entry resolves**

Run:
```bash
PERSONAS=~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline
grep -q "(feedback_personas_governance.md)" "$PERSONAS/shared/MEMORY.md" && test -f "$PERSONAS/shared/feedback_personas_governance.md" && echo "OK: index link resolves"
```
Expected: `OK: index link resolves`.

---

## Task 3: Add the session-start git-status scan to each role's morning routine

**Files:**
- Modify: `$PERSONAS/pm/feedback_morning_routine.md`
- Modify: `$PERSONAS/scientist/feedback_morning_routine.md`
- Modify: `$PERSONAS/developer/feedback_morning_routine.md`

- [ ] **Step 1: Add the PM (caretaker) variant**

In `$PERSONAS/pm/feedback_morning_routine.md`, add as a new bullet at the end of the session-start / first-phase step list (the earliest "at session start" steps, before the news/status phases):

```markdown
- **Personas-repo uncommitted-state scan (session start).** Run `git -C ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline status` and surface any uncommitted / stranded state in chat (one line per touched file). As **caretaker committer** (interim, until MM onboards), land legitimate stranded edits with the user's push gate; stage specific files, never `git add -A`. Prevents the silent pile-up that caused the 2026-05-29 incident. Per [Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) Sub 6 / [Issue #567](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/567) governance.
```

- [ ] **Step 2: Add the Scientist (surface-and-flag) variant**

In `$PERSONAS/scientist/feedback_morning_routine.md`, add at the end of the session-start step list:

```markdown
- **Personas-repo uncommitted-state scan (session start).** Run `git -C ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline status` and surface any uncommitted / stranded state in chat. **Do not commit** — surface it and (if it's your own edit) flag it with a "Memory edits — for MM to commit" bullet; PM-caretaker / MM lands it. Per [Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) Sub 6 / [Issue #567](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/567) governance.
```

- [ ] **Step 3: Add the Developer (surface-and-flag) variant**

In `$PERSONAS/developer/feedback_morning_routine.md`, add the **same bullet as Step 2** (Developer is also surface-and-flag, not a committer).

- [ ] **Step 4: Verify all three files carry the scan**

Run:
```bash
PERSONAS=~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline
for r in pm scientist developer; do
  printf "%s: " "$r"; grep -c "Personas-repo uncommitted-state scan" "$PERSONAS/$r/feedback_morning_routine.md"
done
```
Expected: each prints `1`.

---

## Task 4: PM-caretaker commit + push the personas edits (user-gated)

**Files:** all five personas files from Tasks 1–3.

> ⚠️ The working tree may contain **unrelated** stranded edits (e.g. `developer/MEMORY.md`, `developer/hooks-inactive-after-midsession-settings-edit.md`). **Stage only the governance files by exact path** — never `git add -A` / `git add .`. Surface the unrelated edits separately (see Task 6 notes).

- [ ] **Step 1: Stage only the governance files**

Run:
```bash
PERSONAS=~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline
git -C "$PERSONAS" add \
  shared/feedback_personas_governance.md \
  shared/MEMORY.md \
  pm/feedback_morning_routine.md \
  scientist/feedback_morning_routine.md \
  developer/feedback_morning_routine.md
git -C "$PERSONAS" status --short
```
Expected: exactly those 5 paths staged (one `A`, four `M`); any unrelated edits remain unstaged.

- [ ] **Step 2: Commit (do NOT push yet)**

Run:
```bash
PERSONAS=~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline
git -C "$PERSONAS" commit -m "$(cat <<'EOF'
feat(governance): personas-repo edit/commit policy + session-start git-status scan

Phase 1 of project-repo Issue #567 (sub of #527). Lands the access/lifecycle
policy as memory rules:

- shared/feedback_personas_governance.md (new): edit/commit matrix (own-dir
  edits free; shared/ + cross-role go through MM; MM commits), propose-flows,
  PM-caretaker interim convention, "Memory edits — for MM to commit" flag,
  Phase-2 enforcement roadmap.
- shared/MEMORY.md: replace the 2026-05-26 "not your responsibility" Always-in-
  effect rule with the new model + index the new file.
- pm/scientist/developer/feedback_morning_routine.md: session-start
  `git -C <personas> status` scan (= #527 Sub 6); PM caretaker-commits, Sci/Dev
  surface-and-flag.

Committed by PM-as-caretaker under the bootstrap exception (MM not yet onboarded).
EOF
)"
git -C "$PERSONAS" log -1 --stat
```
Expected: commit created; the stat lists exactly the 5 governance files.

- [ ] **Step 3: Surface the diff and push ONLY after the user confirms**

Run (show the diff, then gate):
```bash
PERSONAS=~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline
git -C "$PERSONAS" show --stat HEAD
```
Then ask the user to confirm. On their OK:
```bash
git -C "$PERSONAS" push origin main
```
Expected: push to personas `origin/main` succeeds. Per the commit/push-separate rule, the push is a distinct, user-gated step.

---

## Task 5: Project-repo lab notebook entry

**Files:**
- Modify: `research/lab_notebook/pm.md` (project repo, on branch `design/pm/issue-567-personas-repo-governance`)

- [ ] **Step 1: Add the entry**

Run `date -u +"%H:%M UTC"` for the timestamp, then add a new time-entry at the TOP of the current date block in `research/lab_notebook/pm.md`:

```markdown
### HH:MM UTC — Editor: PM

#### Personas-repo governance Phase 1 — [Issue #567](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/567) (sub of [#527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527))

**Trigger.** User flagged "too many people touching the personas repo" after a stranding incident. Brainstormed → spec → this Phase-1 landing.

**Decisions captured (not in the Issue body).** MM-curated `shared/` was chosen as the *strict* reading (roles cannot inline-edit `shared/`, must propose) over the looser edit-with-flag — accepted the day-to-day friction for the second-set-of-eyes guarantee on cross-cutting rules. #567 was restructured from a parallel parent into a native sub of #527 after the user caught the duplication (scan = #527 Sub 6; Phase-2 mechanisms = new #527 subs). Per-role MEMORY.md copies of the rule were deliberately *not* made (slimming) — the shared rule covers all roles.

**What landed (personas main, PM-caretaker commit).** `shared/feedback_personas_governance.md` (policy), `shared/MEMORY.md` rule rewrite + index, session-start git-status scan in all three role morning routines.

**Followups.** Phase 2 (per-role git-permission split + edit-boundary PreToolUse hook) filed as #527 subs, gated on MM onboarding (Subs 3–9).
```

- [ ] **Step 2: Commit + push the entry on the #567 branch**

Run:
```bash
cd ~/dev/GitHub/Jin-HoMLee/splice-neoepitope-pipeline-pm
git add research/lab_notebook/pm.md
git commit -m "docs(lab-notebook): PM — personas-repo governance Phase 1 (advances #567)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
git push origin design/pm/issue-567-personas-repo-governance
```
Expected: entry committed + pushed on the branch.

---

## Task 6: Open the PR, tick #567 ACs, merge

**Files:** none (GitHub state).

- [ ] **Step 1: Open the PR for #567 (if not already open)**

Run:
```bash
gh pr create --repo Jin-HoMLee/splice-neoepitope-pipeline \
  --project "JH M Lee Lab" \
  --title "design+feat(pm): personas-repo governance Phase 1 (closes #567)" \
  --body "$(cat <<'EOF'
**Created by:** PM

## Summary

Personas-repo governance ([Issue #567](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/567), sub of [#527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527)). This PR carries the spec, the Phase-1 plan, and the lab-notebook entry. The policy content itself landed in the personas repo (PM-caretaker commit; SHA recorded in the lab notebook).

Closes [Issue #567](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/567).

## Test plan

- [x] shared/feedback_personas_governance.md created (policy matrix + propose-flows + caretaker convention + Phase-2 roadmap)
- [x] shared/MEMORY.md "not your responsibility" rule replaced + new file indexed
- [x] session-start git-status scan added to pm/scientist/developer morning routines
- [x] personas edits committed by PM-caretaker + pushed (user-gated)
- [x] spec + plan committed on this branch
EOF
)"
```
Then flip board Status to "Ready for review" (option `8bf9192f`) per the PR-open checklist.

- [ ] **Step 2: Tick all #567 acceptance criteria**

Run `gh issue view 567 --repo Jin-HoMLee/splice-neoepitope-pipeline --json body`, edit each `- [ ]` whose work is done to `- [x]` (design doc ✓, plan ✓, rule rewrite ✓, caretaker convention ✓, Sub 6 scan ✓), and `gh issue edit 567 --body-file`.

- [ ] **Step 3: Merge via the closure-ritual gate**

Run:
```bash
cd ~/dev/GitHub/Jin-HoMLee/splice-neoepitope-pipeline-pm
bash scripts/audit_and_merge.sh <PR_NUMBER>
```
Expected: audit passes (no unticked boxes on PR Test plan or #567 ACs); squash-merge; #567 auto-closes. Then `git checkout main && git pull`.

- [ ] **Step 4: Handle the unrelated stranded Developer edits (separate caretaker action)**

The `developer/MEMORY.md` + `developer/hooks-inactive-after-midsession-settings-edit.md` edits found during this work are **Developer's own-dir content**, not part of #567. Per the new policy, surface them to the user and either (a) land them as a separate PM-caretaker commit attributed to Developer (user-gated push), or (b) flag for the next Developer/MM session. Do NOT fold them into the governance commit.

---

## Self-Review

**1. Spec coverage:**
- §1 model → Task 1 (the matrix + propose-flows in the policy file) + Task 2 (the Always-in-effect summary). ✓
- §2.1 rule rewrite → Task 2. ✓
- §2.2 git-status scan (= #527 Sub 6) → Task 3. ✓
- §2.3 PM-caretaker interim commit → Task 1 (policy section) + Task 4 (the actual caretaker commit + user gate). ✓
- §3 Phase 2 → out of scope for Phase 1 by design; captured as roadmap text (Task 1) + tracked as #527 subs (not built here). ✓
- §"Rule / file changes" table → Tasks 1–3 cover every Phase-1 row. ✓

**2. Placeholder scan:** No TBD/TODO. The only `<PR_NUMBER>` token (Task 6.3) is a runtime value filled from `gh pr create` output — standard, not a content gap.

**3. Type/identifier consistency:** File paths consistent throughout (`shared/feedback_personas_governance.md`, the three `feedback_morning_routine.md` paths). The "Memory edits — for MM to commit" flag string is used identically in Task 1, Task 2, and Task 3. The scan-bullet heading "Personas-repo uncommitted-state scan" matches across Task 3 Steps 1–3 and the verification grep in Step 4.
