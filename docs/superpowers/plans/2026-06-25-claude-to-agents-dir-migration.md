# `.claude/` → `.agents/` Canonical Dir Migration — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `.agents/` the canonical project-local agent-config directory with `.claude/` a committed symlink to it, across all three clones, without breaking the safety-critical PreToolUse guards or any clone's role memory.

**Architecture:** One tracked commit (authored by hand in the PM clone) does the `git mv` rename + the committed `.claude → .agents` symlink + the `settings.json` hook-path rewrite + the `.gitignore` rewrite. A committed, idempotent, self-verifying script (`scripts/migrate_claude_to_agents.sh`) handles the per-clone gitignored stragglers + the git "can't swap a populated dir for a symlink on pull" gotcha in the two receiving clones.

**Tech Stack:** git, bash, Claude Code settings/hooks model.

**Spec:** `docs/superpowers/specs/2026-06-25-claude-to-agents-dir-migration-design.md`
**Issue:** [#861](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/861)
**Branch:** `chore/developer/issue-861-claude-to-agents-dir` (already created; spec already committed as `d11aeb7`)

## Global Constraints

- Canonical config dir is `.agents/`; `.claude/` is a committed symlink (git mode 120000) pointing at it.
- Hook `command` paths in `settings.json` are absolute-canonical: `${CLAUDE_PROJECT_DIR}/.agents/hooks/<file>.py`. The `recheck_dispatch.py --scope shared` arg is preserved verbatim.
- Never delete or lose `settings.local.json` (per-user permissions) or any clone's `memory` symlink. `hook_fires.jsonl`, `hooks/__pycache__/`, `scheduled_tasks.{lock,json}` are droppable (regenerate).
- The receiving-clone script reads each clone's role from its *existing* `memory` symlink target — never hardcode a role.
- All commits must go up through a PR; `main` is protected (`pipeline-pytest`, `pipeline-snakemake-dry-run`, squash-merge). Merge via `scripts/audit_and_merge.sh`.
- ⚠️ Editing `settings.json` silently deactivates all hooks for the *current* session (no hot-reload). Live-hook ACs are verified only in a **fresh session after merge**.
- Em dash is banned in prose; use a plain dash.

---

### Task 1: Author the tracked migration in the PM clone

**Files:**
- Move: `.claude/{commands,hooks,settings.json}` → `.agents/{commands,hooks,settings.json}` (git-tracked rename)
- Create: `.claude` (symlink → `.agents`)
- Modify: `.agents/settings.json` (hook command paths)
- Modify: `.gitignore` (per-clone artifact paths)

**Interfaces:**
- Produces: the canonical `.agents/` tree + `.claude` symlink that Task 2's script and Task 5's receiving clones rely on; the rewritten `.gitignore` that makes the relocated locals invisible to git in every clone.

- [ ] **Step 1: Confirm starting state**

Run:
```bash
git branch --show-current && git status --porcelain && test -d .claude && echo ".claude is a real dir"
```
Expected: branch `chore/developer/issue-861-claude-to-agents-dir`; clean tree; `.claude is a real dir`.

- [ ] **Step 2: `git mv` the three tracked items**

```bash
git mv .claude/commands .agents/commands
git mv .claude/hooks .agents/hooks
git mv .claude/settings.json .agents/settings.json
```

- [ ] **Step 3: Relocate this clone's gitignored locals + re-create the memory symlink**

```bash
MEM_TARGET="$(readlink .claude/memory)"
echo "saved memory target: $MEM_TARGET"
mv .claude/settings.local.json .agents/settings.local.json
mv .claude/hook_fires.jsonl .agents/hook_fires.jsonl 2>/dev/null || true
rm -f .claude/memory
ln -s "$MEM_TARGET" .agents/memory
rm -rf .claude/hooks 2>/dev/null || true   # any stray __pycache__ already moved with hooks; belt-and-suspenders
```

- [ ] **Step 4: Replace the now-empty `.claude` dir with a symlink**

```bash
rmdir .claude && ln -s .agents .claude
```
If `rmdir` fails, `ls -la .claude` to find the leftover untracked file, relocate/remove it, then retry. Do **not** `rm -rf .claude` blindly.

- [ ] **Step 5: Verify the symlink resolves**

Run:
```bash
test -L .claude && echo "symlink: $(readlink .claude)" && test -f .claude/settings.json && echo ".claude/settings.json resolves" && readlink .agents/memory
```
Expected: `symlink: .agents`; `.claude/settings.json resolves`; the saved memory target.

- [ ] **Step 6: Rewrite the hook command paths in `.agents/settings.json`**

Replace the file contents with (only the five `command` values change vs. the original):
```json
{
  "$schema": "https://json.schemastore.org/claude-code-settings.json",
  "hooks": {
    "PreToolUse": [
      {
        "matcher": "Bash",
        "hooks": [
          {
            "type": "command",
            "command": "${CLAUDE_PROJECT_DIR}/.agents/hooks/check_at_claude.py",
            "if": "Bash(gh *)",
            "timeout": 5
          },
          {
            "type": "command",
            "command": "${CLAUDE_PROJECT_DIR}/.agents/hooks/check_gh_issue_develop_parent.py",
            "if": "Bash(gh issue develop *)",
            "timeout": 15
          },
          {
            "type": "command",
            "command": "${CLAUDE_PROJECT_DIR}/.agents/hooks/check_board_query_pagination.py",
            "if": "Bash(gh api *)",
            "timeout": 5
          }
        ]
      }
    ],
    "PostToolUse": [
      {
        "matcher": "Bash",
        "hooks": [
          {
            "type": "command",
            "command": "${CLAUDE_PROJECT_DIR}/.agents/hooks/post_gh_pr_create.py",
            "if": "Bash(gh *)",
            "timeout": 35
          },
          {
            "type": "command",
            "command": "${CLAUDE_PROJECT_DIR}/.agents/hooks/recheck_dispatch.py --scope shared",
            "if": "Bash(gh *)",
            "timeout": 30
          }
        ]
      }
    ]
  }
}
```

- [ ] **Step 7: Verify `settings.json` is valid JSON**

Run:
```bash
conda activate snakemake && python -c "import json; json.load(open('.agents/settings.json')); print('valid json')"
```
Expected: `valid json`.

- [ ] **Step 8: Rewrite the `.gitignore` per-clone artifact block**

Replace the block currently at lines ~312-323 (the `scheduled_tasks` / `settings.local.json` / `memory` / `hook_fires.jsonl` rules) with:
```gitignore
# Claude Code — transient scheduler state (per-session lock + durable cron jobs)
.agents/scheduled_tasks.lock
.agents/scheduled_tasks.json
.claude/scheduled_tasks.lock
.claude/scheduled_tasks.json

# Claude Code — per-user local settings (permissions allow-list, absolute paths)
.agents/settings.local.json
.claude/settings.local.json

# claude-personas role-memory symlink
/.agents/memory
/.claude/memory
.agents/hook_fires.jsonl
.claude/hook_fires.jsonl
```
(Leave `CLAUDE.local.md` and the generic `__pycache__/` rule untouched.)

- [ ] **Step 9: Stage everything and verify the diff shows only intended changes**

Run:
```bash
git add -A && git status --short
```
Expected: renamed tracked files (`R  .claude/commands/... -> .agents/commands/...`, hooks, `settings.json`), new symlink `A  .claude`, modified `.gitignore`. **No** `settings.local.json`, `hook_fires.jsonl`, or `memory` staged (all gitignored). Confirm the symlink mode:
```bash
git ls-files -s .claude
```
Expected: mode `120000` (symlink).

- [ ] **Step 10: Confirm no gitignored local got staged**

Run:
```bash
git diff --cached --name-only | grep -E 'settings\.local\.json|hook_fires\.jsonl|/memory$' && echo "LEAK - unstage it" || echo "no local leaked into the commit"
```
Expected: `no local leaked into the commit`.

- [ ] **Step 11: Commit**

```bash
git commit -m "chore(repo): migrate .claude/ -> .agents/ canonical dir + symlink back (#861)"
```

---

### Task 2: Receiving-clone migration script

**Files:**
- Create: `scripts/migrate_claude_to_agents.sh`

**Interfaces:**
- Consumes: the merged `.agents/` tree + `.claude` symlink + rewritten `.gitignore` from Task 1 (via `git pull`).
- Produces: a migrated receiving clone (base / scientist) with locals relocated under `.agents/` and `.agents/memory` re-created at that clone's own role target.

- [ ] **Step 1: Write the script**

Create `scripts/migrate_claude_to_agents.sh`:
```bash
#!/usr/bin/env bash
# One-shot, idempotent migration of a receiving clone from a real .claude/ dir
# to the canonical .agents/ dir + .claude symlink (Issue #861).
# Run from the clone's repo root, on `main`, after the migration PR has merged.
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

# Idempotency guard: if .claude is already a symlink (or gone), we're done.
if [ ! -d .claude ] || [ -L .claude ]; then
  echo "[migrate] .claude is not a real directory — already migrated. Nothing to do."
  exit 0
fi

# Preconditions: clean tracked tree (gitignored locals are fine).
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "[migrate] ERROR: tracked working tree is dirty. Commit/stash first." >&2
  exit 1
fi

# 1. Save the role-memory symlink target (role-specific; never hardcode).
if [ ! -L .claude/memory ]; then
  echo "[migrate] ERROR: .claude/memory is not a symlink — unexpected; aborting." >&2
  exit 1
fi
MEM_TARGET="$(readlink .claude/memory)"
ROLE="$(basename "$MEM_TARGET")"
echo "[migrate] role detected from memory symlink: $ROLE  (target: $MEM_TARGET)"

# 2. Preserve settings.local.json; drop the regenerable locals so .claude/ holds
#    only tracked files (so git can swap the dir for the incoming symlink).
TMP="$(mktemp -d)"
[ -f .claude/settings.local.json ] && mv .claude/settings.local.json "$TMP/settings.local.json"
rm -f .claude/memory .claude/hook_fires.jsonl .claude/scheduled_tasks.lock .claude/scheduled_tasks.json
rm -rf .claude/hooks/__pycache__

# 3. Pull the migration commit; git now cleanly removes the tracked files from
#    .claude/ and lays down .agents/ + the .claude symlink.
echo "[migrate] pulling origin/main ..."
git pull --ff-only origin main

# 4. Restore the preserved local + re-create the role-memory symlink under .agents/.
[ -f "$TMP/settings.local.json" ] && mv "$TMP/settings.local.json" .agents/settings.local.json
ln -s "$MEM_TARGET" .agents/memory
rmdir "$TMP" 2>/dev/null || true

# 5. Verify.
fail=0
[ -L .claude ] && [ "$(readlink .claude)" = ".agents" ] || { echo "[verify] FAIL: .claude is not a symlink to .agents" >&2; fail=1; }
[ -f .claude/settings.json ] || { echo "[verify] FAIL: .claude/settings.json does not resolve" >&2; fail=1; }
[ -d "$(readlink -f .agents/memory)" ] && [ "$(basename "$(readlink .agents/memory)")" = "$ROLE" ] || { echo "[verify] FAIL: .agents/memory does not resolve to role '$ROLE'" >&2; fail=1; }
ls .agents/hooks/check_at_claude.py >/dev/null 2>&1 || { echo "[verify] FAIL: .agents/hooks not present" >&2; fail=1; }
[ -z "$(git status --porcelain)" ] || { echo "[verify] FAIL: working tree not clean after migration:" >&2; git status --short >&2; fail=1; }

if [ "$fail" -ne 0 ]; then
  echo "[migrate] VERIFICATION FAILED — inspect above. Locals may be in $TMP." >&2
  exit 1
fi
echo "[migrate] OK — $ROLE clone migrated to .agents/ + .claude symlink; tree clean."
```

- [ ] **Step 2: Make it executable**

```bash
chmod +x scripts/migrate_claude_to_agents.sh
```

- [ ] **Step 3: Syntax-check**

Run:
```bash
bash -n scripts/migrate_claude_to_agents.sh && echo "syntax ok"
```
Expected: `syntax ok`.

- [ ] **Step 4: Idempotency smoke test in the PM clone**

The PM clone is already migrated (Task 1), so the guard must short-circuit cleanly.
Run:
```bash
bash scripts/migrate_claude_to_agents.sh
```
Expected: `[migrate] .claude is not a real directory — already migrated. Nothing to do.` and exit 0.

- [ ] **Step 5: Commit**

```bash
git add scripts/migrate_claude_to_agents.sh
git commit -m "chore(scripts): one-shot receiving-clone migrator for .claude -> .agents (#861)"
```

---

### Task 3: Open the PR and request review

**Files:** none (PR metadata only)

- [ ] **Step 1: Push the branch**

```bash
git push -u origin chore/developer/issue-861-claude-to-agents-dir
```

- [ ] **Step 2: Create the PR**

```bash
gh pr create --title "chore(repo): migrate .claude/ -> .agents/ canonical dir + symlink back (#861)" --body "$(cat <<'EOF'
Closes #861.

Makes `.agents/` the canonical agent-config dir; `.claude/` becomes a committed symlink to it (git mode 120000). Directory-level sibling of the AGENTS.md canonicalization (#857). Organizational/naming change, not new functional portability — see the spec for the honest-scope statement.

## What changed
- `git mv` of the tracked subset (`commands/`, `hooks/`, `settings.json`) `.claude/` -> `.agents/`, plus a committed `.claude -> .agents` symlink.
- Hook command paths rewritten to `${CLAUDE_PROJECT_DIR}/.agents/hooks/...` so guard *execution* never depends on symlink resolution (discovery of `.claude/settings.json` still rides the symlink — belt and suspenders).
- `.gitignore` repointed to the `.agents/` artifact paths (the `.claude/` forms kept too; both resolve to the same inode).
- `scripts/migrate_claude_to_agents.sh`: idempotent, self-verifying one-shot migrator the **receiving** clones (base, scientist) run after merge — it relocates the gitignored locals and re-creates each clone's role-memory symlink, working around git's refusal to replace a populated `.claude/` dir with a symlink on pull.

## Per-clone migration runbook
The tracked rename propagates by pull; the local relocation does not. In **each receiving clone** (base, scientist), after this merges:
```bash
git checkout main
bash scripts/migrate_claude_to_agents.sh   # pulls + relocates locals + re-creates memory symlink + verifies
```
The PM clone authored the change and is already migrated.

## Test plan
- [ ] `.agents/` is the real dir; `.claude` is a committed symlink (`git ls-files -s .claude` shows mode 120000).
- [ ] `.agents/settings.json` is valid JSON; hook paths are `${CLAUDE_PROJECT_DIR}/.agents/hooks/...`.
- [ ] No gitignored local (`settings.local.json`, `hook_fires.jsonl`, `memory`) is tracked.
- [ ] `bash -n scripts/migrate_claude_to_agents.sh` passes; idempotency guard short-circuits in an already-migrated clone.
- [ ] Fresh session post-merge (PM clone): `/hooks` shows all three PreToolUse guards; `/doctor` clean; trace-probe confirms a guard fires.
- [ ] Receiving clones (base, scientist) migrated via the script; role memory loads + a guard fires in each.
- [ ] `docs/` / `AGENTS.md` `.claude/` references still resolve (symlink or updated).
EOF
)"
```

- [ ] **Step 2b: Verify CI is green** — wait for `pipeline-pytest`, `pipeline-snakemake-dry-run`, `pipeline-conda-env-solve` to pass (`gh pr checks --watch`).

- [ ] **Step 3: Offer the bot review**

```bash
gh pr comment <PR#> --body "@claude review"
```
Wait for the review; address findings per the receiving-code-review skill.

---

### Task 4: Lab-notebook entry + merge

**Files:**
- Modify: `research/lab_notebook/developer.md` (Issue #861 is `role:developer`; the closure gate keys off the issue's role label, not the driver)

- [ ] **Step 1: Add the lab-notebook entry**

Append a `## 2026-06-25` (merge date) entry to `research/lab_notebook/developer.md` referencing PR + Issue #861, summarizing the canonical-dir migration and the receiving-clone script. Put each sentence on its own line.

- [ ] **Step 2: Commit + push the entry**

```bash
git add research/lab_notebook/developer.md
git commit -m "docs(notebook): developer entry for .claude -> .agents migration (#861)"
git push
```

- [ ] **Step 3: Tick the PR Test plan + Issue #861 Acceptance criteria** that are satisfiable pre-merge (the structural/JSON/script/no-leak boxes). Leave the fresh-session + receiving-clone boxes for Task 5/6, or tick-with-link if deferring per the closure ritual.

- [ ] **Step 4: Merge through the closure gate**

```bash
bash scripts/audit_and_merge.sh <PR#> --squash --delete-branch
```
Expected: clean audit, then merge. If it blocks on an unticked AC box, resolve per the closure ritual (tick with carrier link or move the box) — do not bypass.

---

### Task 5: Migrate the two receiving clones

**Files:** none (runs the committed script in sibling clones)

- [ ] **Step 1: Migrate the base (developer) clone**

```bash
cd /Users/jin-holee/dev/GitHub/Jin-HoMLee/splice-neoepitope-pipeline
git checkout main
bash scripts/migrate_claude_to_agents.sh
```
Expected: ends with `[migrate] OK — developer clone migrated ... tree clean.`

- [ ] **Step 2: Migrate the scientist clone**

```bash
cd /Users/jin-holee/dev/GitHub/Jin-HoMLee/splice-neoepitope-pipeline-scientist
git checkout main
bash scripts/migrate_claude_to_agents.sh
```
Expected: ends with `[migrate] OK — scientist clone migrated ... tree clean.`

- [ ] **Step 3: Verify both clones independently**

For each clone path, run:
```bash
test -L .claude && readlink .claude && basename "$(readlink .agents/memory)" && git status --porcelain && echo "clean"
```
Expected: `.agents`, the correct role (`developer` / `scientist`), no porcelain output, `clean`.

---

### Task 6: Post-merge verification + docs sweep

**Files:**
- Possibly modify: `AGENTS.md`, `docs/remote_routines.md` (only live operational `.claude/` instructions, if any)

- [ ] **Step 1: Sync the PM clone to merged main**

```bash
cd /Users/jin-holee/dev/GitHub/Jin-HoMLee/splice-neoepitope-pipeline-pm
git checkout main && git pull origin main --ff-only && git status --porcelain && echo "PM clone clean on merged main"
```
Expected: `.claude` symlink + `.agents/` intact, clean tree.

- [ ] **Step 2: Docs sweep — find live `.claude/` references**

Run:
```bash
git grep -n "\.claude/" -- AGENTS.md 'docs/*.md' | grep -v 'docs/superpowers/'
```
For each hit, decide: historical record (leave — and it still resolves via the symlink) vs. a forward instruction telling someone to edit/visit a path (update to `.agents/`). Plan/spec docs under `docs/superpowers/` are historical — leave them. Update only live operational prose. Confirm any left-as-`.claude/` path still resolves through the symlink.

- [ ] **Step 3: If any doc edits were made, commit them**

Open a tiny follow-up PR (or fold into a docs commit on a fresh branch) — do not push to `main` directly.

- [ ] **Step 4: Fresh-session live-hook verification (USER / restart required)**

This **cannot** be done in the authoring session (editing `settings.json` killed hooks for it). In a **fresh** Claude Code session in the PM clone:
- `/hooks` — confirm all three PreToolUse guards (`check_at_claude`, `check_gh_issue_develop_parent`, `check_board_query_pagination`) are listed.
- `/doctor` — confirm clean.
- Trace-probe: attempt a blocked pattern (e.g. a `gh issue comment ... --body "@claude foo"` dry intent) and confirm the guard denies it; check `.agents/hook_fires.jsonl` grows.

- [ ] **Step 5: Tick remaining Issue #861 ACs** once Steps 1-4 + Task 5 pass; close-out per the morning-routine closure audit if not already closed by the merge.

---

## Self-Review

**Spec coverage:**
- §1 goal/scope → Task 1 (rename+symlink), honest-scope carried into PR body. ✓
- §2 inventory (tracked vs local; per-clone memory targets) → Task 1 Steps 2-3, Task 2 script role-detection. ✓
- §3 pull gotcha → Task 2 script Steps 2-3 (clear locals → pull → restore). ✓
- §4a originating-by-hand → Task 1. §4b receiving script → Task 2. ✓
- §5 settings.json path rewrite → Task 1 Step 6. ✓
- §6 .gitignore rewrite → Task 1 Step 8. ✓
- §7 rollout + mid-session hook-death caveat → Tasks 3-6, Task 6 Step 4. ✓
- §8 ACs → mapped across PR Test plan (Task 3) + Tasks 5-6. ✓
- §9 risks → mitigations embedded (verify steps, move-aside, idempotency guard, no-leak check). ✓

**Placeholder scan:** `<PR#>` is a runtime value, not a placeholder. No TBD/TODO/"handle edge cases". Script + JSON + commands are complete. ✓

**Type/name consistency:** `scripts/migrate_claude_to_agents.sh`, `.agents/memory`, `${CLAUDE_PROJECT_DIR}/.agents/hooks/...`, role-from-`readlink` — used identically across tasks. ✓
