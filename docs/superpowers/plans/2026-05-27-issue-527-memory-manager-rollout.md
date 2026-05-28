# Memory Manager Rollout — Implementation Plan (Subs 3-9 of Issue #527)

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bootstrap Memory Manager (MM) as a 4th project role per the design doc at `docs/superpowers/specs/2026-05-27-memory-manager-role-design.md`, completing Subs 3-9 of [parent Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527).

**Architecture:** 6-phase migration — Phase 2 (personas-repo structural setup) and Phases 3-4 (memory file edits) are done by PM sessions using `git -C` against the personas repo as a migration-window exception. Phase 5 (first MM session bootstrap) is user-initiated. Project-repo PRs land lab notebook entries that journal each phase. Personas-repo gets direct commits (no PR convention there).

**Tech Stack:** `gh` CLI, `git` (with `-C` flag for personas-repo migration tasks), GitHub Projects v2 (board #9 cross-repo via `gh api graphql`), Claude Code v2.1.20+ (for `--add-dir` in Phase 5).

**Out of scope:**
- Sub 10 (cross-repo bot blind-spot follow-up — `scripts/audit_and_merge.sh` personas-commit gate). Deferred to its own Dev-tier plan; orthogonal to MM bootstrap completion.
- Any further memory-curation work after Sub 9 (handled by MM in steady state).

**Prerequisite:** Sub 1 (design doc) shipped via [PR #529](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/529) merged 2026-05-27. This plan ships as Sub 2 (its own sub-issue + PR) before Sub 3 (the first build task) starts — same shape as Sub 1.

---

## File Structure

### Files to create (personas repo: `~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/`)

| Path | Responsibility |
|---|---|
| `CLAUDE.md` | Personas-repo project instructions for sessions opened with cwd = personas repo (MM sessions primarily) |
| `.claude/settings.json` | Permissions + hooks scoped to MM session work |
| `memory_manager/MEMORY.md` | MM role memory index with initial Always-in-effect placeholder |
| `memory_manager/shared` (symlink → `../shared`) | Shared-memory access for MM |

### Files to modify (personas repo)

| Path | Change |
|---|---|
| `shared/feedback_team_structure.md` | Add Memory Manager subsection |
| `shared/feedback_domain_bespoke_roles.md` | Add MM workload-separation justification paragraph |
| `shared/feedback_memory_escalation.md` | Rewrite "On you forgot X" rule per Section 3 of design doc |
| `shared/feedback_lab_notebook.md` | Add "Memory edits — for MM to commit" bullet convention + note MM exempt from lab notebook entries |
| `shared/MEMORY.md` | Drop "personas-repo git state is not your responsibility" Always-in-effect line; replace with MM-ownership pointer |
| `pm/MEMORY.md`, `scientist/MEMORY.md`, `developer/MEMORY.md` | Replace personas-repo no-touch reflex with `role:memory_manager` Issue routing |
| `pm/feedback_morning_routine.md`, `scientist/feedback_morning_routine.md`, `developer/feedback_morning_routine.md` | Add `git status` scan on personas repo at session start |

### Files to create (project repo: `~/dev/GitHub/Jin-HoMLee/splice-neoepitope-pipeline-pm/`)

| Path | Responsibility |
|---|---|
| `research/lab_notebook/pm.md` (append) | Per-task lab notebook entries journaling each Sub's PM-side work |

### GitHub state changes

- Enable Issues on `Jin-HoMLee/claude-personas-splice-neoepitope-pipeline` (one-time repo setting via web UI or `gh api`).
- Create label `role:memory_manager` on `Jin-HoMLee/claude-personas-splice-neoepitope-pipeline`.
- Relabel existing curation Issues #248, #326, #346, #353 + file & relabel the MEMORY.md slimming audit Issue.

---

## Task 1: Sub 3 — Personas-repo CLAUDE.md + `.claude/` setup

**Files:**
- Create: `~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/CLAUDE.md`
- Create: `~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/.claude/settings.json`
- Append: `research/lab_notebook/pm.md` (project-repo lab notebook entry)

**Note on workflow:** PM session does this work via `git -C` against the personas repo. This is the migration-window exception (PM normally wouldn't touch personas-repo git, but MM doesn't exist yet to do it).

- [ ] **Step 1.1: File Sub 3 Issue + link to parent #527 as native sub**

Run:

```bash
SUB3_URL=$(gh issue create \
  --title "feat(roles): Sub 3 of #527 — personas-repo CLAUDE.md + .claude/ setup" \
  --milestone "pm-i6 - PM Tooling, Memory & Methodology II" \
  --label "role:pm,role:memory_manager,enhancement" \
  --project "JH M Lee Lab" \
  --body "$(cat <<'EOF'
**Created by:** PM

## Trigger

Sub 3 of [parent Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) (MM rollout epic). Lands the personas-repo workspace configuration so MM sessions opened with cwd = personas repo load conventions + permissions correctly.

## Scope

Add personas-repo root `CLAUDE.md` (project instructions for personas-repo sessions) + `.claude/settings.json` (minimum permissions for `gh` + `git`). Files land via PM `git -C` direct commit to personas-repo main (migration-window exception per design doc Phase 2).

## Acceptance criteria

- [ ] `CLAUDE.md` exists at personas-repo root with sections: repo structure overview, memory file conventions (file-relative paths), role-mapping (which dir = which role), `<role>/shared` symlink mechanism, personas-repo-specific gotchas (no PR convention; direct-to-main).
- [ ] `.claude/settings.json` exists with `gh` + `git` permissions appropriate for MM session work.
- [ ] Both files committed direct to personas-repo `main` via PM `git -C`.
- [ ] Project-repo lab notebook entry journals the work.
- [ ] PR opened, reviewed, merged via `scripts/audit_and_merge.sh`.

## Connects to

- [parent Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) — this is Sub 3.
- Design doc `docs/superpowers/specs/2026-05-27-memory-manager-role-design.md` Phase 2.

**Priority rationale:** P2 — inherited from parent.
EOF
)" 2>&1 | tail -1)

SUB3_NUM=$(echo "$SUB3_URL" | grep -oE '[0-9]+$')

# Link as native sub of parent #527
PARENT_DB_ID=$(gh api repos/Jin-HoMLee/splice-neoepitope-pipeline/issues/527 --jq .id)
SUB3_DB_ID=$(gh api repos/Jin-HoMLee/splice-neoepitope-pipeline/issues/${SUB3_NUM} --jq .id)
echo "{\"sub_issue_id\":${SUB3_DB_ID}}" | gh api repos/Jin-HoMLee/splice-neoepitope-pipeline/issues/527/sub_issues --input -
echo "Sub 3 = #${SUB3_NUM}"
```

Expected: Sub Issue created, returned URL, linked as native sub of #527.

- [ ] **Step 1.2: Sync project-repo main + create branch via `gh issue develop`**

Run:

```bash
cd ~/dev/GitHub/Jin-HoMLee/splice-neoepitope-pipeline-pm
git fetch origin && git pull origin main
gh issue develop ${SUB3_NUM} --checkout --name "feat/pm/issue-${SUB3_NUM}-personas-repo-claude-md"
```

Expected: branch checked out locally.

- [ ] **Step 1.3: Sync personas-repo main (via `git -C`)**

Run:

```bash
git -C ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline fetch origin
git -C ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline pull origin main
```

Expected: personas-repo main up-to-date.

- [ ] **Step 1.4: Draft personas-repo `CLAUDE.md`**

Write to `~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/CLAUDE.md`:

```markdown
# Claude Personas — Splice Neoepitope Pipeline

Memory repo for the multi-role workflow on `Jin-HoMLee/splice-neoepitope-pipeline`. Sessions opened with this repo as cwd are **Memory Manager (MM) sessions**; PM/Sci/Dev sessions operate on the project repo and access this memory via per-role symlinks (`~/.claude/projects/<project-hash>/memory/` → `<role>/`).

## Repo structure

```
claude-personas-splice-neoepitope-pipeline/
├── shared/             ← cross-role memories (loaded by all roles)
├── developer/          ← Developer-role memories
│   └── shared → ../shared
├── pm/                 ← PM-role memories
│   └── shared → ../shared
├── scientist/          ← Scientist-role memories
│   └── shared → ../shared
└── memory_manager/     ← Memory Manager memories (own role)
    └── shared → ../shared
```

Each `<role>/shared` symlink lets every role read shared memories transparently.

## Memory file conventions

- **Paths are file-relative.** `<!-- src: -->` annotations, cross-references in markdown links, and prose mentions all use paths relative to the file containing the reference. Same-dir files: just the filename. Shared files from a role dir: `shared/<file>` (via the symlink). Shared files from inside `shared/`: just the filename.
- **Never prefix with `.claude/memory/...`** — that string silently resolves to `~/.claude/memory/...` and reads the wrong file.

## Role mapping

| Dir | Role | Workspace |
|---|---|---|
| `pm/` | Project Manager | Project-repo clone (board, milestones, triage) |
| `scientist/` | Scientist | Project-repo clone (manuscript, biology) |
| `developer/` | Developer | Project-repo clone (code, pipeline, infra) |
| `memory_manager/` | Memory Manager | **This repo** (curation, git lifecycle) |

## Personas-repo gotchas

- **No PR convention.** Direct commits to `main` (shared "commit, push, merge — three separate steps" rule still applies: stage + commit + surface diff + push on user OK).
- **Git lifecycle is MM-owned.** PM/Sci/Dev can edit personas-repo memory files mid-session (transparently via per-role symlinks) but never commit + push from those sessions. The current edits sit uncommitted until next MM session.
- **Designed for cwd = this repo.** Sessions opened elsewhere can read personas files via per-role symlinks but should not edit + push.

## See also

- Design doc: `docs/superpowers/specs/2026-05-27-memory-manager-role-design.md` in the project repo.
- Parent Issue: [#527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) — MM rollout epic.
- `shared/feedback_team_structure.md` — full role definitions.
```

- [ ] **Step 1.5: Verify `CLAUDE.md` content matches design doc requirements**

Spot-check:
- Repo structure diagram present? ✓
- Memory file conventions section explains file-relative paths + the `.claude/memory/` foot-gun? ✓
- Role mapping table includes all 4 roles? ✓
- Personas-repo gotchas mentions no-PR + MM-owns-git + cwd-design? ✓
- Cross-references to design doc + parent + team_structure.md? ✓

- [ ] **Step 1.6: Draft personas-repo `.claude/settings.json`**

Write to `~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/.claude/settings.json`:

```json
{
  "$schema": "https://json.schemastore.org/claude-code-settings.json",
  "permissions": {
    "allow": [
      "Bash(git status)",
      "Bash(git diff:*)",
      "Bash(git log:*)",
      "Bash(git fetch origin)",
      "Bash(git pull origin main)",
      "Bash(git add:*)",
      "Bash(git commit:*)",
      "Bash(git push)",
      "Bash(git push origin main)",
      "Bash(git restore:*)",
      "Bash(git checkout:*)",
      "Bash(gh issue list:*)",
      "Bash(gh issue view:*)",
      "Bash(gh issue create:*)",
      "Bash(gh issue edit:*)",
      "Bash(gh issue close:*)",
      "Bash(gh issue comment:*)",
      "Bash(gh label list:*)",
      "Bash(gh api repos/Jin-HoMLee/splice-neoepitope-pipeline/issues/*)",
      "Bash(gh api repos/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/issues/*)",
      "Bash(gh api graphql:*)"
    ]
  }
}
```

- [ ] **Step 1.7: Verify `.claude/settings.json` parses as valid JSON**

Run:

```bash
python -c "import json; json.load(open('/Users/jin-holee/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/.claude/settings.json'))" && echo OK
```

Expected: `OK` printed. (Uses bare `python` after `conda activate snakemake` per the Python env rule, or system Python — JSON parsing is stdlib so any 3.x works.)

- [ ] **Step 1.8: Commit both files to personas repo via `git -C`**

Run:

```bash
git -C ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline add CLAUDE.md .claude/settings.json
git -C ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline status
git -C ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline commit -m "$(cat <<'EOF'
feat: personas-repo CLAUDE.md + .claude/ setup for MM sessions

Sub 3 of project-repo Issue #527 (MM rollout epic). Lands personas-repo
workspace configuration so MM sessions opened with cwd = this repo load
conventions + permissions correctly.

CLAUDE.md documents: repo structure, memory file conventions (file-
relative paths), role mapping, personas-repo gotchas (no PR convention,
MM-owned git lifecycle, cwd-design).

.claude/settings.json grants minimum permissions for MM session work
(git operations on personas-repo, gh CLI for both repos via API).

Migration-window exception: PM session committed via git -C since MM
doesn't exist yet. Subsequent personas-repo commits will be MM-only.
EOF
)"
```

Expected: commit created in personas-repo main.

- [ ] **Step 1.9: Push personas-repo to remote (after user confirms diff)**

Run:

```bash
git -C ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline log -1 --stat
echo "Review above diff. Confirm before push: y/N?"
read -r CONFIRM
[ "$CONFIRM" = "y" ] && git -C ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline push origin main
```

Expected: push to personas-repo origin/main after user confirmation.

Per shared `feedback_commit_push_separate.md`: never chain commit and push without user gate. The `read` prompt is the gate.

- [ ] **Step 1.10: Test that personas-repo `CLAUDE.md` loads in a Claude Code session**

Manual check (user action, not agentic):

```bash
cd ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline
claude
# In the spawned session, ask: "What does the personas repo CLAUDE.md say about role mapping?"
# Expected: session quotes the role-mapping table from CLAUDE.md
# Then /exit
```

Expected: CLAUDE.md content is visible to the session at start.

- [ ] **Step 1.11: Write project-repo lab notebook entry**

Run `date -u +"%H:%M UTC"` first to get the timestamp, then append to `research/lab_notebook/pm.md` directly under the existing `## 2026-05-27` header (newest-first), above existing time entries:

```markdown
### HH:MM UTC — Editor: PM

#### Sub 3 of [parent Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) (MM rollout) — personas-repo CLAUDE.md + .claude/ setup

**Trigger.** Sub 3 task on the MM rollout. Personas-repo workspace configuration so MM sessions opened with cwd = personas repo load conventions + permissions correctly.

**What landed.** Two new files in personas-repo main (committed via PM `git -C` per the migration-window exception):
- `CLAUDE.md` at repo root: repo structure, memory file conventions, role mapping (now 4 roles including MM), personas-repo gotchas (no PR, MM-owned git lifecycle, cwd-design).
- `.claude/settings.json`: minimum `gh` + `git` permissions for MM session work.

**Verification.** Opened `claude` in personas-repo cwd; verified CLAUDE.md content loaded into session context. JSON-parsed settings.json successfully.

**Followups.** Sub 4 (`memory_manager/MEMORY.md` initial content) next. Closes [Sub Issue #N](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/N) via [PR #M](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/M).
```

(Replace HH:MM, Sub Issue #N, PR #M placeholders with actual values once known.)

- [ ] **Step 1.12: Commit project-repo lab notebook entry + push branch**

Run:

```bash
cd ~/dev/GitHub/Jin-HoMLee/splice-neoepitope-pipeline-pm
git add research/lab_notebook/pm.md
git commit -m "$(cat <<EOF
docs(lab-notebook): PM 2026-05-27 HH:MM UTC — Sub 3 of #527 (personas-repo CLAUDE.md + .claude/)

Closes #${SUB3_NUM}

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
git push -u origin "feat/pm/issue-${SUB3_NUM}-personas-repo-claude-md"
```

Expected: commit + push to origin.

- [ ] **Step 1.13: Open PR closing Sub 3 + flip board Status**

Run:

```bash
gh pr create --title "feat(pm): Sub 3 of #527 — personas-repo CLAUDE.md + .claude/ setup" \
  --project "JH M Lee Lab" \
  --body "$(cat <<EOF
**Created by:** PM

## Summary

Lab notebook entry for PM 2026-05-27 HH:MM UTC. Lands the personas-repo workspace configuration (CLAUDE.md + .claude/settings.json) so MM sessions opened with cwd = personas repo load conventions + permissions correctly.

Substantive files committed direct to personas-repo main via PM \`git -C\` per the migration-window exception (commit SHA: \`<personas-commit-sha>\`). This PR carries only the project-repo lab notebook journal.

Closes [Issue #${SUB3_NUM}](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/${SUB3_NUM}) (Sub 3 of [parent Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527)).

## Test plan

- [x] CLAUDE.md content matches design doc Section 2 requirements (repo structure, memory conventions, role mapping, gotchas)
- [x] .claude/settings.json parses as valid JSON
- [x] CLAUDE.md verified to load into session opened with personas-repo cwd
- [x] Personas-repo commit SHA recorded in PR body + lab notebook
- [x] Branch cut from origin/main (rebase-before-write satisfied)
EOF
)"
```

Set PR Status to "Ready for review" via GraphQL (replace `<PR_NUM>` after pr create):

```bash
PR_NUM=$(gh pr view --json number --jq .number)
PR_ITEM=$(gh api graphql -f query='{repository(owner:"Jin-HoMLee",name:"splice-neoepitope-pipeline"){pullRequest(number:'${PR_NUM}'){projectItems(first:5){nodes{id project{number}}}}}}' --jq '.data.repository.pullRequest.projectItems.nodes[] | select(.project.number == 9) | .id')
gh api graphql -f query="mutation{updateProjectV2ItemFieldValue(input:{projectId:\"PVT_kwHOB17eGc4BSomP\",itemId:\"${PR_ITEM}\",fieldId:\"PVTSSF_lAHOB17eGc4BSomPzhAHFf8\",value:{singleSelectOptionId:\"8bf9192f\"}}){projectV2Item{id}}}"
```

- [ ] **Step 1.14: Wait for CI + tick last ACs + merge via `audit_and_merge.sh`**

Run:

```bash
gh pr checks ${PR_NUM}
# Expected: ci-tools-pytest, pipeline-pytest, pipeline-snakemake-dry-run all pass

# Tick the last AC on the Sub Issue
gh issue view ${SUB3_NUM} --json body --jq '.body' > /tmp/sub3_body.md
# Manually edit /tmp/sub3_body.md to tick all [ ] boxes to [x]
gh issue edit ${SUB3_NUM} --body-file /tmp/sub3_body.md

# Merge
bash scripts/audit_and_merge.sh ${PR_NUM}
```

Expected: `✓ PR #N merged ...` — Sub 3 auto-closes via `Closes #N` keyword.

- [ ] **Step 1.15: Sync local main**

Run:

```bash
git checkout main
git pull origin main
```

Expected: fast-forward to merged commit.

---

## Task 2: Sub 4 — `memory_manager/MEMORY.md` initial content

**Files:**
- Create: `~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/memory_manager/MEMORY.md`
- Create: `~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/memory_manager/shared` (symlink → `../shared`)
- Append: `research/lab_notebook/pm.md`

- [ ] **Step 2.1: File Sub 4 Issue + link to parent #527 as native sub**

Run:

```bash
SUB4_URL=$(gh issue create \
  --title "feat(roles): Sub 4 of #527 — memory_manager/MEMORY.md initial content" \
  --milestone "pm-i6 - PM Tooling, Memory & Methodology II" \
  --label "role:pm,role:memory_manager,enhancement" \
  --project "JH M Lee Lab" \
  --body "$(cat <<'EOF'
**Created by:** PM

## Trigger

Sub 4 of [parent Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) (MM rollout epic). Creates the `memory_manager/` role directory and its MEMORY.md index, plus the `shared` symlink for cross-role memory access. Required before Phase 5 first MM session bootstrap.

## Scope

- Create `memory_manager/` dir in personas repo.
- Create `memory_manager/MEMORY.md` with frontmatter + initial Always-in-effect placeholder section (MM fleshes out content in Phase 5).
- Create `memory_manager/shared` symlink to `../shared`.
- Lab notebook entry + project-repo PR for journaling.

## Acceptance criteria

- [ ] `memory_manager/MEMORY.md` exists with frontmatter header + Always-in-effect placeholder.
- [ ] `memory_manager/shared` symlink exists and resolves to `../shared`.
- [ ] Files committed direct to personas-repo main via PM `git -C`.
- [ ] Project-repo lab notebook entry journals the work.
- [ ] PR opened, reviewed, merged via `scripts/audit_and_merge.sh`.

## Connects to

- [parent Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) — this is Sub 4.
- Design doc Phase 2.

**Priority rationale:** P2 — inherited from parent.
EOF
)" 2>&1 | tail -1)

SUB4_NUM=$(echo "$SUB4_URL" | grep -oE '[0-9]+$')
SUB4_DB_ID=$(gh api repos/Jin-HoMLee/splice-neoepitope-pipeline/issues/${SUB4_NUM} --jq .id)
echo "{\"sub_issue_id\":${SUB4_DB_ID}}" | gh api repos/Jin-HoMLee/splice-neoepitope-pipeline/issues/527/sub_issues --input -
echo "Sub 4 = #${SUB4_NUM}"
```

- [ ] **Step 2.2: Sync + create branch**

Run:

```bash
cd ~/dev/GitHub/Jin-HoMLee/splice-neoepitope-pipeline-pm
git fetch origin && git pull origin main
gh issue develop ${SUB4_NUM} --checkout --name "feat/pm/issue-${SUB4_NUM}-memory-manager-memory-md"

git -C ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline fetch origin
git -C ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline pull origin main
```

- [ ] **Step 2.3: Create `memory_manager/` dir + symlink + MEMORY.md**

Run:

```bash
PERSONAS=~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline
mkdir -p ${PERSONAS}/memory_manager
ln -s ../shared ${PERSONAS}/memory_manager/shared
```

Write to `${PERSONAS}/memory_manager/MEMORY.md`:

```markdown
# Memory Index — Memory Manager

> Paths in `<!-- src: -->` annotations and index links are **file-relative**. Bare filename = same dir; `shared/<file>` = the shared sibling dir (resolved via the `shared` symlink). **Do not prefix with `.claude/memory/...`** — that string silently resolves to `~/.claude/memory/...` and reads the wrong file.

## Always in effect (no file read required)

- **Workspace is the personas repo.** MM sessions open with cwd = `claude-personas-splice-neoepitope-pipeline/`. Use bare `git ...` commands (no `-C`). Project-repo read access via `--add-dir` at session launch. <!-- src: shared/feedback_team_structure.md -->
- **Personas-repo direct-to-main commits.** No PR convention here; commit + push directly. "Commit, push, merge — three separate steps" rule still applies: stage + commit + surface diff + push on user OK. <!-- src: shared/feedback_github_workflow.md -->
- **Lab notebook is the personas-repo git log.** MM is exempt from `research/lab_notebook/<role>.md` entries — every commit message captures what + why, serving the journaling purpose. <!-- src: shared/feedback_lab_notebook.md ("Entry timing" section, MM exemption clause) -->
- **MM full ownership of memory edits + git lifecycle.** Any role can read/edit personas-repo memory files mid-session, but only MM commits + pushes. Active roles flag uncommitted state via the "Memory edits — for MM to commit" bullet in their lab notebook entry. <!-- src: shared/feedback_memory_escalation.md -->

(MM fleshes out additional Always-in-effect rules in Phase 5 bootstrap session based on early experience.)

## Shared (all sessions)

- [Shared memory index](shared/MEMORY.md) — All cross-role conventions (board rules, lab notebook, GitHub workflow, team structure).

## Role: Memory Manager

(populated in Phase 5 bootstrap — initial drafts for: morning routine for MM sessions; cross-rule consistency check pattern; consolidation methodology; audit script invocation patterns.)
```

- [ ] **Step 2.4: Verify symlink + content**

Run:

```bash
PERSONAS=~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline
ls -la ${PERSONAS}/memory_manager/
# Expected:
# - memory_manager/MEMORY.md (file)
# - memory_manager/shared -> ../shared (symlink)

readlink ${PERSONAS}/memory_manager/shared
# Expected: ../shared

ls ${PERSONAS}/memory_manager/shared/MEMORY.md
# Expected: file exists (resolves to ../shared/MEMORY.md)
```

- [ ] **Step 2.5: Commit + push personas-repo**

Run:

```bash
PERSONAS=~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline
git -C ${PERSONAS} add memory_manager/
git -C ${PERSONAS} status
git -C ${PERSONAS} commit -m "$(cat <<'EOF'
feat: memory_manager/ role dir + initial MEMORY.md

Sub 4 of project-repo Issue #527 (MM rollout epic). Creates the MM
role directory parallel to pm/, scientist/, developer/. Initial MEMORY.md
captures 4 placeholder Always-in-effect rules — workspace, direct-to-main
commit convention, lab notebook exemption, full ownership of memory
git lifecycle.

memory_manager/shared symlinks to ../shared per the existing per-role
symlink pattern.

MM fleshes out further role-specific content in Phase 5 bootstrap.

Migration-window exception: PM session committed via git -C.
EOF
)"
git -C ${PERSONAS} log -1 --stat
echo "Review diff. Confirm push: y/N?"
read -r CONFIRM
[ "$CONFIRM" = "y" ] && git -C ${PERSONAS} push origin main
```

- [ ] **Step 2.6: Write project-repo lab notebook entry**

Run `date -u +"%H:%M UTC"` and append to `research/lab_notebook/pm.md` directly under `## 2026-05-27` header (newest-first), above prior time entries:

```markdown
### HH:MM UTC — Editor: PM

#### Sub 4 of [parent Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) (MM rollout) — memory_manager/ dir + initial MEMORY.md

**Trigger.** Sub 4 of MM rollout. Creates the MM role dir + memory index, parallel to existing pm/, scientist/, developer/.

**What landed.** Personas-repo main now contains:
- `memory_manager/` dir.
- `memory_manager/MEMORY.md` with frontmatter + 4 placeholder Always-in-effect rules (workspace, direct-to-main commits, lab notebook exemption, full git ownership).
- `memory_manager/shared → ../shared` symlink.

**Verification.** `ls -la memory_manager/` shows expected file + symlink; `readlink memory_manager/shared` = `../shared`; reading through the symlink resolves to `shared/MEMORY.md`.

**Followups.** Sub 5 (shared memory updates) next. Closes [Issue #N](url) via [PR #M](url).
```

- [ ] **Step 2.7: Commit + push project-repo, open PR, flip Status, merge**

Run (same shape as Task 1 Steps 1.12-1.15, replacing SUB3 with SUB4 and personas-commit-sha with current value):

```bash
cd ~/dev/GitHub/Jin-HoMLee/splice-neoepitope-pipeline-pm
git add research/lab_notebook/pm.md
git commit -m "$(cat <<EOF
docs(lab-notebook): PM 2026-05-27 HH:MM UTC — Sub 4 of #527 (memory_manager/MEMORY.md initial)

Closes #${SUB4_NUM}

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
git push -u origin "feat/pm/issue-${SUB4_NUM}-memory-manager-memory-md"

# Open PR
gh pr create --title "feat(pm): Sub 4 of #527 — memory_manager/MEMORY.md initial content" \
  --project "JH M Lee Lab" \
  --body "$(cat <<EOF
**Created by:** PM

## Summary

Lab notebook entry for Sub 4 of [parent Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527). Creates the \`memory_manager/\` role dir + initial MEMORY.md + shared symlink in the personas repo (committed direct to personas-repo main via PM \`git -C\`).

Closes [Issue #${SUB4_NUM}](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/${SUB4_NUM}).

## Test plan

- [x] memory_manager/MEMORY.md frontmatter + 4 placeholder Always-in-effect rules present
- [x] memory_manager/shared symlink resolves to ../shared
- [x] Files committed direct to personas-repo main (SHA: \`<personas-commit-sha>\`)
- [x] Branch cut from origin/main
EOF
)"

# Set PR Status → Ready for review (same GraphQL pattern as Task 1)
# ...
# Wait for CI, tick ACs, merge via audit_and_merge.sh
```

Expected: PR opens, CI passes, audit_and_merge.sh confirms ticks + squash-merges, Sub 4 auto-closes.

- [ ] **Step 2.8: Sync local main**

Run:

```bash
git checkout main && git pull origin main
```

---

## Task 3: Sub 5 — Shared memory updates (5 files bundled)

**Files (all in personas repo `shared/`):**
- Modify: `shared/feedback_team_structure.md` — add Memory Manager subsection
- Modify: `shared/feedback_domain_bespoke_roles.md` — add MM workload-separation paragraph
- Modify: `shared/feedback_memory_escalation.md` — rewrite "On you forgot X" rule
- Modify: `shared/feedback_lab_notebook.md` — MM exemption + new "Memory edits — for MM to commit" bullet convention
- Modify: `shared/MEMORY.md` — drop personas-repo no-touch rule; replace with MM-ownership pointer
- Append: `research/lab_notebook/pm.md`

- [ ] **Step 3.1: File Sub 5 Issue + link to parent**

Run (same pattern as Task 1.1, replacing the body's `Sub 3` references with `Sub 5`, title with `Sub 5 of #527 — shared memory updates`).

Body acceptance criteria:

```markdown
- [ ] `shared/feedback_team_structure.md` has a "Memory Manager" subsection with cwd, infrastructure, role-label semantics.
- [ ] `shared/feedback_domain_bespoke_roles.md` has an MM workload-separation paragraph.
- [ ] `shared/feedback_memory_escalation.md` rewritten "On you forgot X" rule per design doc Section 3.
- [ ] `shared/feedback_lab_notebook.md` notes MM exempt from lab notebook entries + adds "Memory edits — for MM to commit" bullet convention.
- [ ] `shared/MEMORY.md` Always-in-effect line "personas-repo git state is not your responsibility" replaced with MM-ownership pointer.
- [ ] Files committed direct to personas-repo main.
- [ ] Project-repo lab notebook entry journals the work.
- [ ] PR opened, reviewed, merged via `scripts/audit_and_merge.sh`.
```

- [ ] **Step 3.2: Sync + create branch**

Same as Task 1.2/2.2 but with `SUB5_NUM`.

- [ ] **Step 3.3: Edit `shared/feedback_team_structure.md` — add MM subsection**

Read the file first to find the right insertion point (after Developer section, before "Why" rationale). Add:

```markdown
**Memory Manager session** — personas-repo curation, memory consistency, git lifecycle
- Cwd at session start: personas repo (`claude-personas-splice-neoepitope-pipeline/`)
- Reads project repo via `--add-dir` at session launch
- Owns: `git add`/`commit`/`push` on the personas repo; the queued memory-curation Issues (`role:memory_manager`); cross-rule consistency sweeps
- Does NOT post to team_standup; does NOT write `research/lab_notebook/<role>.md` (commits are the journal)
- Active roles (PM/Sci/Dev) can EDIT personas-repo memory files mid-session via per-role symlinks; only MM COMMITS them
```

And in the Memory directory layout section, update the diagram to include `memory_manager/`:

```
claude-personas-splice-neoepitope-pipeline/
├── shared/             ← cross-role memories
├── developer/          ← Developer-role memories
│   └── shared → ../shared
├── memory_manager/     ← Memory Manager memories  (NEW)
│   └── shared → ../shared
├── pm/                 ← PM-role memories
│   └── shared → ../shared
└── scientist/          ← Scientist-role memories
    └── shared → ../shared
```

- [ ] **Step 3.4: Verify the edit**

Run:

```bash
grep -A 5 "Memory Manager session" ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/shared/feedback_team_structure.md
# Expected: shows the new subsection

grep -A 3 "memory_manager/" ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/shared/feedback_team_structure.md
# Expected: shows the directory diagram with memory_manager/
```

- [ ] **Step 3.5: Edit `shared/feedback_domain_bespoke_roles.md` — add MM justification**

Add a new section before "Cross-reference" (last section):

```markdown
## Memory Manager — a meta-process role, not a bioinformatics one

The Memory Manager role added 2026-05-27 (via [parent Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527)) is **deliberately not** justified by bioinformatics-domain analogy. Its justification is **workload separation**: 5-7 of the PM Issue queue at the time of role creation were memory-curation work (slimming audits, dedup, audit scripts, rule lifts) that conflated with real project-management work, plus an unowned personas-repo git lifecycle gap.

External framing: "We added a Memory Manager because our memory infrastructure grew structured enough to need a curator; not because other multi-agent setups have one." The role is cross-domain (meta-process — would apply to any multi-role Claude Code setup with growing memory infrastructure), while PM/Sci/Dev remain bioinformatics-bespoke.

Workload validation gate at +4 weeks post-bootstrap (Phase 6 of the design doc) determines whether the role's workload is sustained enough to keep; rollback path is explicit.
```

- [ ] **Step 3.6: Verify**

```bash
grep -A 1 "Memory Manager — a meta-process role" ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/shared/feedback_domain_bespoke_roles.md
# Expected: shows new section header
```

- [ ] **Step 3.7: Rewrite `shared/feedback_memory_escalation.md`**

Read the file first to understand the existing "On you forgot X" rule. Replace the rule body with the new form per design doc Section 3:

```markdown
## On "you forgot X"

Find the existing memory. If it's only behind a link, copy it inline here under Always in effect. If no memory exists, create one AND inline it here. The edit ships at the next Memory Manager session: the active role edits inline (transparently through the per-role symlink); MM commits + pushes the personas repo on next MM session.

For sustained curation needs (multi-rule consolidation, audit scripts, slimming work), file a `role:memory_manager` Issue rather than inlining.

**Migration-window note:** Until the morning-routine `git status` scan on the personas repo is universally adopted (per the role morning-routine updates in `<role>/feedback_morning_routine.md`), active roles flag uncommitted personas-repo edits via the "Memory edits — for MM to commit" bullet in their lab notebook entry (see `shared/feedback_lab_notebook.md`).
```

- [ ] **Step 3.8: Verify the rewrite**

```bash
grep -A 3 'On "you forgot X"' ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/shared/feedback_memory_escalation.md
# Expected: shows new rule text starting with "Find the existing memory..."
```

- [ ] **Step 3.9: Edit `shared/feedback_lab_notebook.md` — MM exemption + new bullet convention**

Find the "When required vs optional — non-routine sessions only" section. Add at the end of the "Optional / skip" subsection:

```markdown
**MM sessions are exempt entirely.** Memory Manager doesn't write to `research/lab_notebook/<role>.md` at all — every personas-repo commit message captures what + why, serving the journaling purpose. The personas-repo git log IS MM's lab notebook.
```

Then find the "Insertion order" section (or the section listing lab notebook entry bullets) and add:

```markdown
## "Memory edits — for MM to commit" bullet convention

When a PM/Sci/Dev session edits personas-repo memory files mid-session (via the per-role symlink), the lab notebook entry must include a bullet under a "Memory edits — for MM to commit" sub-heading listing the touched files:

```markdown
**Memory edits — for MM to commit:**
- `shared/feedback_X.md` — added "Y" rule
- `<role>/MEMORY.md` — inlined the Y rule
```

This is the structural handoff to the next MM session — MM's morning routine includes `git status` on the personas repo and reads these bullets for context.

**Migration-window note:** Until Sub 6 of MM rollout adds `git status` scan to each role's morning routine, this bullet convention is the only handoff signal; afterwards it's still required (the bullet adds context the `git status` alone can't surface — *why* the edit, not just *what*).
```

- [ ] **Step 3.10: Verify**

```bash
grep -B 1 -A 5 "MM sessions are exempt entirely" ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/shared/feedback_lab_notebook.md
# Expected: shows the exemption clause

grep -B 1 -A 5 "Memory edits — for MM to commit" ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/shared/feedback_lab_notebook.md
# Expected: shows the bullet convention section
```

- [ ] **Step 3.11: Edit `shared/MEMORY.md` — replace personas-repo no-touch rule**

Find the Always-in-effect line:

```
- **Personas-repo git state is not your responsibility:** Reading memory files for content is in scope; managing their git lifecycle is not. Never probe personas-repo uncommitted state with `git -C`...
```

Replace with:

```
- **Personas-repo git lifecycle is Memory Manager's responsibility.** PM/Sci/Dev can read AND edit personas-repo memory files mid-session via per-role symlinks (treat as auto-memory dir). Only MM commits + pushes. Active-session edits ship at the next MM session; flag in the lab notebook entry under "Memory edits — for MM to commit" sub-heading. Established 2026-05-27 via [parent Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) (MM rollout); replaces the 2026-05-26 "not your responsibility" rule. Full role definition: `shared/feedback_team_structure.md` (Memory Manager subsection).
```

- [ ] **Step 3.12: Verify the replacement**

```bash
grep -B 1 -A 2 "Personas-repo git lifecycle is Memory Manager" ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/shared/MEMORY.md
# Expected: shows new rule

grep "Personas-repo git state is not your responsibility" ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/shared/MEMORY.md
# Expected: NO output (old rule should be gone)
```

- [ ] **Step 3.13: Commit + push personas-repo bundle**

Run:

```bash
PERSONAS=~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline
git -C ${PERSONAS} add shared/feedback_team_structure.md shared/feedback_domain_bespoke_roles.md shared/feedback_memory_escalation.md shared/feedback_lab_notebook.md shared/MEMORY.md
git -C ${PERSONAS} diff --cached --stat
git -C ${PERSONAS} commit -m "$(cat <<'EOF'
feat: shared memory updates for Memory Manager role

Sub 5 of project-repo Issue #527. Bundles 5 interrelated shared/* edits:

- feedback_team_structure.md: add Memory Manager subsection + update
  memory directory diagram to include memory_manager/.
- feedback_domain_bespoke_roles.md: add MM workload-separation
  justification (not domain analogy); externalised role-roster framing.
- feedback_memory_escalation.md: rewrite "On you forgot X" rule —
  active role edits inline, MM commits later; for sustained curation,
  file role:memory_manager Issue.
- feedback_lab_notebook.md: MM exempt from lab notebook entries
  (personas-repo commits ARE the journal) + new "Memory edits — for
  MM to commit" bullet convention for active-role lab notebook entries.
- MEMORY.md: replace "personas-repo git state is not your
  responsibility" Always-in-effect rule with new MM-ownership pointer.

Migration-window exception: PM committed via git -C.
EOF
)"
git -C ${PERSONAS} log -1 --stat
echo "Review diff. Confirm push: y/N?"
read -r CONFIRM
[ "$CONFIRM" = "y" ] && git -C ${PERSONAS} push origin main
```

- [ ] **Step 3.14: Lab notebook entry + project-repo PR + merge**

Same shape as Task 1 Steps 1.11-1.15 / Task 2 Steps 2.6-2.8, replacing SUB3/SUB4 with SUB5. Lab notebook entry should highlight that this is a 5-file shared-memory bundle and reference the personas commit SHA. PR title: `feat(pm): Sub 5 of #527 — shared memory updates`.

---

## Task 4: Sub 6 — Role MEMORY.md updates + morning routine git status scan

**Files (all in personas repo):**
- Modify: `pm/MEMORY.md` — replace personas-repo no-touch reflex with MM Issue routing
- Modify: `scientist/MEMORY.md` — same
- Modify: `developer/MEMORY.md` — same
- Modify: `pm/feedback_morning_routine.md` — add `git status` scan on personas repo
- Modify: `scientist/feedback_morning_routine.md` — same
- Modify: `developer/feedback_morning_routine.md` — same
- Append: `research/lab_notebook/pm.md`

- [ ] **Step 4.1: File Sub 6 Issue + link to parent**

Run (same pattern as previous tasks). Title: `feat(roles): Sub 6 of #527 — role MEMORY.md updates + morning routine git status scan`.

Body acceptance criteria:

```markdown
- [ ] PM/Sci/Dev MEMORY.md Always-in-effect: drop personas-repo no-touch reflex; replace with `role:memory_manager` Issue routing pointer.
- [ ] PM/Sci/Dev feedback_morning_routine.md: add `git status` on personas repo at session start (surfaces uncommitted memory edits from prior sessions for chat acknowledgment).
- [ ] Files committed direct to personas-repo main.
- [ ] Project-repo lab notebook entry journals the work.
- [ ] PR opened, reviewed, merged via `scripts/audit_and_merge.sh`.
```

- [ ] **Step 4.2: Sync + create branch**

Same as previous tasks with `SUB6_NUM`.

- [ ] **Step 4.3: Edit PM/Sci/Dev MEMORY.md — Always-in-effect updates**

For each role's `MEMORY.md`, find any Always-in-effect rule that says "personas-repo is not your responsibility" (or equivalent) and replace with:

```markdown
- **Personas-repo memory edits.** Read AND edit personas-repo memory files mid-session as natural part of work (per-role symlinks make this transparent). Never commit + push from this session — Memory Manager owns the git lifecycle. For sustained curation work (slimming, dedup, audit), file a `role:memory_manager` Issue. Flag in-session edits via the lab notebook "Memory edits — for MM to commit" bullet. <!-- src: shared/MEMORY.md (Memory Manager rule, 2026-05-27) -->
```

If no such Always-in-effect line exists in a role's MEMORY.md (Sci/Dev may not have one), add the rule above as a new bullet under Always-in-effect.

- [ ] **Step 4.4: Verify each role file**

```bash
PERSONAS=~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline
for role in pm scientist developer; do
  echo "=== ${role}/MEMORY.md ==="
  grep -A 1 "Personas-repo memory edits" ${PERSONAS}/${role}/MEMORY.md
  grep "personas-repo is not your responsibility\|Personas-repo git state is not your responsibility" ${PERSONAS}/${role}/MEMORY.md && echo "WARNING: old rule still present in ${role}/MEMORY.md"
done
# Expected: new rule appears in each; old rule absent
```

- [ ] **Step 4.5: Edit PM/Sci/Dev feedback_morning_routine.md — add git status scan**

For each role's `feedback_morning_routine.md`, find the section describing session-start steps (likely "Phase 1" or "Step −1" or similar). Add at the end of session-start steps:

```markdown
**Personas-repo uncommitted state scan.** At session start, run `git -C ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/ status` and surface any uncommitted state in chat (one-line acknowledgment per touched file). This is the active-role's side of the MM handoff loop: if PM/Sci/Dev session edited personas memory mid-session and forgot to flag, this scan catches it. Established 2026-05-27 via [parent Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) (MM rollout).
```

- [ ] **Step 4.6: Verify each role's morning routine file**

```bash
PERSONAS=~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline
for role in pm scientist developer; do
  echo "=== ${role}/feedback_morning_routine.md ==="
  grep -A 2 "Personas-repo uncommitted state scan" ${PERSONAS}/${role}/feedback_morning_routine.md
done
# Expected: new section appears in each
```

- [ ] **Step 4.7: Commit + push personas-repo (6 files bundle)**

Run:

```bash
PERSONAS=~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline
git -C ${PERSONAS} add pm/MEMORY.md scientist/MEMORY.md developer/MEMORY.md pm/feedback_morning_routine.md scientist/feedback_morning_routine.md developer/feedback_morning_routine.md
git -C ${PERSONAS} diff --cached --stat
git -C ${PERSONAS} commit -m "$(cat <<'EOF'
feat: role MEMORY.md updates + morning-routine git status scan for MM rollout

Sub 6 of project-repo Issue #527. Updates each of PM, Scientist, Developer
role MEMORY.md to replace personas-repo no-touch reflex with the new MM-
ownership rule (active roles can edit personas memory; only MM commits +
pushes). Adds personas-repo `git status` scan to each role's morning
routine so uncommitted memory edits from prior sessions surface at chat
start.

Migration-window exception: PM committed via git -C.
EOF
)"
git -C ${PERSONAS} log -1 --stat
echo "Review diff. Confirm push: y/N?"
read -r CONFIRM
[ "$CONFIRM" = "y" ] && git -C ${PERSONAS} push origin main
```

- [ ] **Step 4.8: Lab notebook entry + project-repo PR + merge**

Same shape as previous tasks, replacing SUB5 with SUB6. Lab notebook entry highlights the 6-file role-memory + morning-routine bundle.

---

## Task 5: Sub 7 — Enable Issues + `role:memory_manager` label on personas repo

**Note:** This is GitHub admin work — no file changes in either repo. Project-repo PR carries only a lab notebook entry.

- [ ] **Step 5.1: File Sub 7 Issue + link to parent**

Run (same pattern). Title: `feat(github): Sub 7 of #527 — enable Issues + role:memory_manager label on personas repo`.

Body acceptance criteria:

```markdown
- [ ] Issues enabled on `Jin-HoMLee/claude-personas-splice-neoepitope-pipeline` (repo setting).
- [ ] Label `role:memory_manager` created on personas repo (color: `8B4F9F` to match project-repo label).
- [ ] Verified by creating a smoke-test Issue on personas repo with the label + immediately closing it (with a comment explaining it's a smoke test).
- [ ] Project-repo lab notebook entry journals the work.
- [ ] PR opened, reviewed, merged via `scripts/audit_and_merge.sh`.
```

- [ ] **Step 5.2: Enable Issues on personas repo**

Run:

```bash
gh api -X PATCH repos/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline -f has_issues=true
```

Expected: API returns 200; `has_issues` now `true`.

- [ ] **Step 5.3: Verify Issues enabled**

```bash
gh api repos/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline --jq '.has_issues'
# Expected: true
```

- [ ] **Step 5.4: Create `role:memory_manager` label on personas repo**

```bash
gh label create "role:memory_manager" \
  --description "Memory Manager work (curation, consistency, lifts)" \
  --color "8B4F9F" \
  --repo Jin-HoMLee/claude-personas-splice-neoepitope-pipeline
```

Expected: label created (no output on success, or label URL).

- [ ] **Step 5.5: Smoke-test by creating + closing a test Issue**

```bash
SMOKE_URL=$(gh issue create --repo Jin-HoMLee/claude-personas-splice-neoepitope-pipeline \
  --title "smoke test: verify role:memory_manager label routing" \
  --label "role:memory_manager" \
  --body "Smoke test for Sub 7 of project-repo Issue #527. Verifies Issues are enabled on personas repo + role:memory_manager label exists. Closing immediately." 2>&1 | tail -1)
SMOKE_NUM=$(echo "$SMOKE_URL" | grep -oE '[0-9]+$')

gh issue close ${SMOKE_NUM} --repo Jin-HoMLee/claude-personas-splice-neoepitope-pipeline \
  --comment "Smoke test passed: Issues enabled + role:memory_manager label routing works. Closing per Sub 7 AC."
```

Expected: Issue creates with label, then closes successfully.

- [ ] **Step 5.6: Lab notebook entry + project-repo PR + merge**

Same shape as previous tasks. Lab notebook entry highlights: Issues now enabled on personas repo, label created, smoke test verified.

---

## Task 6: Sub 8 — Relabel existing curation Issues + file slimming audit Issue

**Files:**
- (no file changes; GitHub state changes)
- Append: `research/lab_notebook/pm.md`

- [ ] **Step 6.1: File Sub 8 Issue + link to parent**

Run (same pattern). Title: `chore(github): Sub 8 of #527 — relabel curation Issues + file MEMORY.md slimming audit Issue`.

Body acceptance criteria:

```markdown
- [ ] Existing Issues #248, #326, #346, #353 relabeled: add `role:memory_manager` (keep `role:pm` as bystander).
- [ ] MEMORY.md slimming audit Issue filed in project repo (was previously only a `project_memory_md_slimming.md` note in PM role memory), labeled `role:memory_manager`, milestone pm-i6.
- [ ] All 5 Issues now appear under `gh issue list --label role:memory_manager`.
- [ ] Project-repo lab notebook entry journals the work.
- [ ] PR opened, reviewed, merged via `scripts/audit_and_merge.sh`.
```

- [ ] **Step 6.2: Relabel #248, #326, #346, #353**

```bash
for n in 248 326 346 353; do
  gh issue edit ${n} --add-label "role:memory_manager"
  echo "Relabeled #${n}"
done
```

Expected: each Issue gets the label added (existing `role:pm` stays).

- [ ] **Step 6.3: Verify relabels**

```bash
gh issue list --label role:memory_manager --state open --json number,title --jq '.[] | "\(.number): \(.title)"'
# Expected: at least 4 Issues listed (#248, #326, #346, #353)
```

- [ ] **Step 6.4: File MEMORY.md slimming audit Issue**

Read `/Users/jin-holee/.claude/projects/-Users-jin-holee-dev-GitHub-Jin-HoMLee-splice-neoepitope-pipeline-pm/memory/project_memory_md_slimming.md` first to extract the audit's scope (29-rule cliff context, per-rule verdicts already drafted, etc.). Then:

```bash
SLIM_URL=$(gh issue create \
  --title "chore(memory): MEMORY.md slimming audit — shared/MEMORY.md 29-rule cliff" \
  --milestone "pm-i6 - PM Tooling, Memory & Methodology II" \
  --label "role:memory_manager" \
  --project "JH M Lee Lab" \
  --body "$(cat <<'EOF'
**Created by:** PM

## Trigger

The PM role memory `project_memory_md_slimming.md` (2026-05-27) handed over an audit: `shared/MEMORY.md` Always-in-effect section has 29 rules vs the documented ~14-rule cliff; PM sessions load 39 rules in total at start. Audit hands over per-rule verdicts (keep / hook / skill / compress / delete) + suggested Issue carve. Non-urgent; was sitting on PM's plate until MM rollout (Issue #527) routed it here.

## Scope

Execute the audit handed over in `pm/project_memory_md_slimming.md`:
- Review each of the 29 Always-in-effect rules in `shared/MEMORY.md` against the per-rule verdict.
- For "keep" verdicts: no action.
- For "compress" verdicts: shorten the rule text without losing meaning; cite the source link.
- For "hook" verdicts: file Dev-tier sub-issue for the relevant `.claude/hooks/` mechanism replacement; remove the rule once hook lands.
- For "skill" verdicts: file Sub-issue to promote the rule into a `superpowers:*` skill; remove the rule once skill lands.
- For "delete" verdicts: remove the rule with rationale in commit message.
- Target: ≤14 rules post-audit, matching the documented cliff.

## Acceptance criteria

- [ ] Each of the 29 rules has a verdict-implementing action committed.
- [ ] `shared/MEMORY.md` Always-in-effect section has ≤14 rules.
- [ ] No information loss: each removed/compressed rule traceable to the replacing mechanism (hook, skill, source-file link).
- [ ] Personas-repo direct commits (MM owns lifecycle).

## Connects to

- [parent Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) — MM rollout this Issue now belongs to.
- PM role memory `project_memory_md_slimming.md` — original audit handover.

**Priority rationale:** P2 — non-urgent; MM picks this up in the validation-gate-period as a representative curation task.
EOF
)" 2>&1 | tail -1)

SLIM_NUM=$(echo "$SLIM_URL" | grep -oE '[0-9]+$')
echo "Slimming audit Issue = #${SLIM_NUM}"
```

- [ ] **Step 6.5: Verify slimming audit Issue now in MM queue**

```bash
gh issue list --label role:memory_manager --state open --json number --jq '. | length'
# Expected: 5 or more (4 existing relabels + 1 new slimming Issue)
```

- [ ] **Step 6.6: Lab notebook entry + project-repo PR + merge**

Same shape as previous tasks. Lab notebook entry highlights: 4 existing Issues relabeled, slimming audit Issue filed (#N), MM queue now ≥5 Issues.

---

## Task 7: Sub 9 — Bootstrap first MM session + 4-week validation gate

**Note:** This task is mostly user-initiated (opening MM session) + observational (4-week validation). The agentic-worker scope here is preparing the bootstrap instructions, validation checklist, and post-validation gate Issue closure. The actual session-opening is user action.

- [ ] **Step 7.1: File Sub 9 Issue + link to parent**

Run (same pattern). Title: `ops: Sub 9 of #527 — bootstrap first MM session + 4-week validation`.

Body acceptance criteria:

```markdown
- [ ] One-time `~/.claude/projects/<personas-hash>/memory/` symlink wired by user.
- [ ] First MM session opened: `cd ~/dev/.../claude-personas-splice-neoepitope-pipeline && claude --add-dir ~/dev/.../splice-neoepitope-pipeline`.
- [ ] MM validation checklist passed: memory dir loads, shared symlink works, --add-dir reaches project repo, CLAUDE.md content visible to session.
- [ ] MM fleshes out `memory_manager/MEMORY.md` Always-in-effect rules beyond the 4 initial placeholders (recorded as a personas-repo commit).
- [ ] MM picks first concrete curation Issue (recommend: slimming audit from Sub 8).
- [ ] First non-bootstrap personas-repo direct commit from MM session.
- [ ] 4-week validation review at +28 days (PM responsibility, observation-only).
```

- [ ] **Step 7.2: Hash discovery instructions for user**

Surface to user (do not run as part of agentic execution):

```
To wire the personas-repo memory dir:

1. Run `claude` once in the personas-repo cwd to let Claude Code auto-create the projects path:

   $ cd ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline
   $ claude
   (in the spawned session, just /exit immediately)

2. Find the new hash dir:

   $ ls -t ~/.claude/projects/ | head -1
   <hash>

3. Symlink the memory dir to the personas-repo's memory_manager/:

   $ ln -s ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/memory_manager ~/.claude/projects/<hash>/memory

4. Verify by re-opening claude in the personas-repo cwd — the session should load memory_manager/MEMORY.md as the role index.
```

This requires user action; not agentic. The agentic worker writes this as instructions in the Sub 9 Issue body and waits for user confirmation that wiring is done.

- [ ] **Step 7.3: First MM session launch instructions**

Once Step 7.2 complete, surface to user:

```
Open first MM session:

   $ cd ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline
   $ claude --add-dir ~/dev/GitHub/Jin-HoMLee/splice-neoepitope-pipeline

Run the validation checklist:
- /memory check loads memory_manager/MEMORY.md
- Read shared/MEMORY.md (transparent through the symlink) — confirm new Always-in-effect MM-ownership rule from Sub 5 is present
- Run `ls ../splice-neoepitope-pipeline/research/lab_notebook/` to confirm --add-dir works (no permissions error)
- Run bare `git status` and confirm cwd = personas repo

If all checks pass, MM session proceeds to first curation task (slimming audit per Sub 8's Issue #N).
```

- [ ] **Step 7.4: MM session fleshes out memory_manager/MEMORY.md (MM action, not PM)**

(Executed by user inside the first MM session, not by this PM-tier agent.) Per design doc Phase 5: MM expands the 4 initial placeholder Always-in-effect rules with additional rules informed by early experience. Persona-repo direct commit.

- [ ] **Step 7.5: Schedule 4-week validation review reminder**

Run:

```bash
# Schedule a one-shot reminder 28 days from today (2026-06-24)
gh issue comment 527 --body "Validation reminder: scheduled for 2026-06-24 (+28 days from MM bootstrap). PM session at that date runs the validation gate per design doc Phase 6 — measure role:memory_manager Issue throughput, MM session count, personas-repo commit cadence, subjective PM cognitive-load delta. Pass → keep MM; fail → rollback per design doc."
```

(Note: this is a comment-based reminder. If `CronCreate` tool is available + appropriate, schedule via cron alternatively — but cron is session-only, so a comment on parent #527 is durable.)

- [ ] **Step 7.6: Close Sub 9 after bootstrap verification**

After user confirms Steps 7.2-7.4 complete + first MM curation work commit landed:

```bash
gh issue close ${SUB9_NUM} --comment "Bootstrap complete. First MM session opened with cwd = personas repo, --add-dir mounted project repo, memory_manager/MEMORY.md loaded, shared symlink resolves correctly. MM picked up first curation Issue (slimming audit #N) and pushed first non-bootstrap commit (personas SHA: \`<sha>\`). 4-week validation reminder scheduled on parent #527."
```

- [ ] **Step 7.7: Lab notebook entry (PM session, observational)**

PM session journals the bootstrap event in `research/lab_notebook/pm.md`:

```markdown
### HH:MM UTC — Editor: PM

#### Sub 9 of [parent Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) (MM rollout) — first MM session bootstrap

**Trigger.** Phase 5 of the MM rollout. User opened first MM session, validated infrastructure, picked first curation Issue.

**What landed.**
- Personas-repo hash symlink wired (`~/.claude/projects/<hash>/memory` → `memory_manager/`).
- First MM session opened: cwd = personas repo, --add-dir mounted project repo.
- Validation checklist passed (memory loads, symlinks resolve, --add-dir works, CLAUDE.md visible).
- MM fleshed out `memory_manager/MEMORY.md` Always-in-effect (added <N> rules beyond initial placeholders).
- First non-bootstrap curation commit: slimming audit progress (Issue #<N> first commit, personas SHA `<sha>`).

**Validation gate scheduled.** 2026-06-24 (+28 days). Per design doc Phase 6: measure throughput + cadence + subjective load delta. Pass → keep MM. Fail → rollback (relabel Issues back to role:pm, retire MM session pattern, keep personas-repo structure as inert dead code).

**Followups.** Sub 10 (cross-repo bot blind-spot follow-up — `scripts/audit_and_merge.sh` personas-commit gate) deferred to Dev-tier plan; orthogonal to this rollout.
```

- [ ] **Step 7.8: Project-repo PR + merge (PM-tier observational record)**

Same shape as previous Subs. PR title: `ops(pm): Sub 9 of #527 — first MM session bootstrap + 4-week validation scheduled`.

---

## Wrap-up: Parent Issue #527 stays Open

After all Subs 3-9 merge + auto-close, the parent Issue #527 still has:
- Sub 10 unchecked (deferred to Dev-tier plan)
- 4-week validation gate pending

Per `feedback_parent_sub_issues.md`, the parent closes via summary comment after ALL subs Done. Since Sub 10 is deferred, parent #527 stays Open until either (a) Sub 10 ships, or (b) validation gate passes and decision is made to close parent without Sub 10 (folding the gate work into a separate pm-i7 Issue).

**Final state at end of this plan:** Subs 3-9 merged; parent Open pending Sub 10 + validation. Memory Manager role operational.

---

## Self-Review

**1. Spec coverage** — each design-doc Phase 2-5 mapped to a Task here:
- Phase 2 (personas-repo structural setup) → Tasks 1 (Sub 3) + 2 (Sub 4) ✓
- Phase 3 (shared memory updates) → Task 3 (Sub 5) ✓
- Phase 4 (role MEMORY.md + Issue relabels) → Tasks 4 (Sub 6) + 6 (Sub 8) ✓
- Phase 4 sub-step "enable Issues + label on personas repo" → Task 5 (Sub 7) ✓
- Phase 5 (first MM session bootstrap + validation reminder) → Task 7 (Sub 9) ✓
- Phase 6 (4-week validation gate) → Step 7.5 schedules; the actual review is a future PM task at +28d (not in this plan's scope by design)
- Sub 10 → explicitly deferred (out-of-scope note at top + Wrap-up section)

**2. Placeholder scan** — searched for: TBD, TODO, "implement later", "add appropriate error handling", "similar to Task N". One conditional remains: `<personas-commit-sha>` and `<sha>` placeholders in PR/lab notebook bodies — these are intentional, populated at execution time from the actual commit (engineer fills in). Documented as such.

**3. Type/identifier consistency:**
- `SUB3_NUM` / `SUB4_NUM` etc. used consistently as the Issue number variable
- `${PERSONAS}` used consistently for the personas-repo path
- Project repo path `~/dev/GitHub/Jin-HoMLee/splice-neoepitope-pipeline-pm/` and personas path `~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/` distinct and unambiguous throughout
- "Sub 3" through "Sub 9" labels match parent Issue #527's sub-issue list ✓

**4. Ambiguity:** Pre-execution clarifications captured:
- Migration-window exception (PM commits via git -C until Phase 5) — stated upfront + repeated in each commit message
- MM exemption from lab notebook (so Sub 9 has no project-repo PR for the MM session's own work; only PM's observational entry) — stated explicitly in Sub 9 task
- 4-week validation is observation-only — explicitly outside execution scope of this plan

---
