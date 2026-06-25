#!/usr/bin/env bash
# One-shot, idempotent migration of a receiving clone from a real .claude/ dir
# to the canonical .agents/ dir + .claude symlink (Issue #861).
# Run from the clone's repo root, on `main`, after the migration PR has merged.
#
# BOOTSTRAP GAP (Issue #866): a receiving clone cannot run *its own* copy of this
# script before the migration, because the script is delivered *by* the very
# `git pull` it performs (step 3) - pre-pull it does not exist in that clone.
# Work around it by invoking an ALREADY-MIGRATED clone's copy with the receiving
# clone as the working dir:
#     cd /path/to/receiving-clone && bash /path/to/migrated-clone/scripts/migrate_claude_to_agents.sh
# The `cd "$(git rev-parse --show-toplevel)"` below pins all operations to the
# receiving clone regardless of where the script file itself lives. (A future
# generic migrator should likewise not assume it is already present locally.)
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

# Idempotency guard: if .claude is already a symlink (or gone), we're done.
if [ ! -d .claude ] || [ -L .claude ]; then
  echo "[migrate] .claude is not a real directory - already migrated. Nothing to do."
  exit 0
fi

# Must be on main: step 3 fast-forwards origin/main into the *current* branch,
# which is only correct on main (the runbook checks out main first; this guard
# makes the script self-contained and refuses to pull main into a feature
# branch). Checked after the idempotency guard so an already-migrated clone
# exits 0 cleanly from any branch.
CURRENT_BRANCH="$(git rev-parse --abbrev-ref HEAD)"
if [ "$CURRENT_BRANCH" != "main" ]; then
  echo "[migrate] ERROR: must run on 'main' (currently on '$CURRENT_BRANCH'). Run: git checkout main" >&2
  exit 1
fi

# Preconditions: clean tracked tree (gitignored locals are fine).
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "[migrate] ERROR: tracked working tree is dirty. Commit/stash first." >&2
  exit 1
fi

# 1. Save the role-memory symlink target (role-specific; never hardcode).
if [ ! -L .claude/memory ]; then
  echo "[migrate] ERROR: .claude/memory is not a symlink - unexpected; aborting." >&2
  exit 1
fi
MEM_TARGET="$(readlink .claude/memory)"
ROLE="$(basename "$MEM_TARGET")"
echo "[migrate] role detected from memory symlink: $ROLE  (target: $MEM_TARGET)"

# 2. Preserve settings.local.json; drop the regenerable locals so .claude/ holds
#    only tracked files (so git can swap the dir for the incoming symlink).
TMP="$(mktemp -d)"
restored=0   # flips to 1 once settings.local.json is moved back out of $TMP
recover_on_error() {
  echo "[migrate] ABORTED mid-migration." >&2
  if [ "$restored" -eq 1 ]; then
    echo "[migrate]   - settings.local.json (if any) was already restored to: .agents/settings.local.json" >&2
  else
    echo "[migrate]   - settings.local.json (if it existed) is staged in: $TMP" >&2
  fi
  echo "[migrate]   - re-create the role-memory symlink manually:  ln -s '$MEM_TARGET' .agents/memory   (or .claude/memory if .claude is still a directory)" >&2
}
trap recover_on_error ERR
[ -f .claude/settings.local.json ] && mv .claude/settings.local.json "$TMP/settings.local.json"
rm -f .claude/memory .claude/hook_fires.jsonl .claude/scheduled_tasks.lock .claude/scheduled_tasks.json
rm -rf .claude/hooks/__pycache__

# 3. Pull the migration commit; git now cleanly removes the tracked files from
#    .claude/ and lays down .agents/ + the .claude symlink.
echo "[migrate] pulling origin/main ..."
git pull --ff-only origin main

# 4. Restore the preserved local + re-create the role-memory symlink under .agents/.
[ -f "$TMP/settings.local.json" ] && mv "$TMP/settings.local.json" .agents/settings.local.json
restored=1   # from here on, locals live at .agents/, not $TMP
ln -s "$MEM_TARGET" .agents/memory
rmdir "$TMP" 2>/dev/null || true

# 5. Verify.
fail=0
[ -L .claude ] && [ "$(readlink .claude)" = ".agents" ] || { echo "[verify] FAIL: .claude is not a symlink to .agents" >&2; fail=1; }
[ -f .claude/settings.json ] || { echo "[verify] FAIL: .claude/settings.json does not resolve" >&2; fail=1; }
[ -d .agents/memory ] && [ "$(basename "$(readlink .agents/memory)")" = "$ROLE" ] || { echo "[verify] FAIL: .agents/memory does not resolve to role '$ROLE'" >&2; fail=1; }
ls .agents/hooks/check_at_claude.py >/dev/null 2>&1 || { echo "[verify] FAIL: .agents/hooks not present" >&2; fail=1; }
[ -z "$(git status --porcelain)" ] || { echo "[verify] FAIL: working tree not clean after migration:" >&2; git status --short >&2; fail=1; }

if [ "$fail" -ne 0 ]; then
  trap - ERR   # avoid a second, redundant recovery dump on the exit below
  echo "[migrate] VERIFICATION FAILED - inspect above. Preserved settings.local.json (if any) is at .agents/settings.local.json." >&2
  exit 1
fi
trap - ERR
echo "[migrate] OK - $ROLE clone migrated to .agents/ + .claude symlink; tree clean."
