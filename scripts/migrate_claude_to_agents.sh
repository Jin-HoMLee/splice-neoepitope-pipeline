#!/usr/bin/env bash
# One-shot, idempotent migration of a receiving clone from a real .claude/ dir
# to the canonical .agents/ dir + .claude symlink (Issue #861).
# Run from the clone's repo root, on `main`, after the migration PR has merged.
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

# Idempotency guard: if .claude is already a symlink (or gone), we're done.
if [ ! -d .claude ] || [ -L .claude ]; then
  echo "[migrate] .claude is not a real directory - already migrated. Nothing to do."
  exit 0
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
recover_on_error() {
  echo "[migrate] ABORTED mid-migration." >&2
  echo "[migrate]   - settings.local.json (if it existed) was staged in: $TMP" >&2
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
  echo "[migrate] VERIFICATION FAILED - inspect above. Locals may be in $TMP." >&2
  exit 1
fi
trap - ERR
echo "[migrate] OK - $ROLE clone migrated to .agents/ + .claude symlink; tree clean."
