# Hook Fire-Log Infrastructure — Implementation Plan ([Issue #453](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/453))

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build deterministic fire-log instrumentation + threshold-crossing promotion prompt for the 3 sibling PostToolUse recheck hooks, replacing the "vibes-check" gate on hook graduation from `.claude/settings.local.json` → committed `.claude/settings.json`.

**Architecture:** Centralized writer in `recheck_dispatch.py` — 5 helper functions + 1 wiring refactor + 1 module-level HOOK_CONFIG constant. Scripts under `scripts/pm/` stay untouched (fire predicate sentinel-matches their existing `Status: [No change]` output). One new bash aggregator (`scripts/check_hook_health.sh`) for operator inspection at promotion time.

**Tech Stack:** Python 3 (stdlib only — `json`, `pathlib`, `datetime`, `sys`), pytest (in `workflow/tests/.venv`), bash + `jq` for the aggregator.

**Spec:** [`docs/superpowers/specs/2026-05-28-hook-fire-log-infra-design.md`](../specs/2026-05-28-hook-fire-log-infra-design.md)

---

## Task 1: Scaffold test file + .gitignore

**Files:**
- Create: `workflow/tests/test_recheck_dispatch.py`
- Modify: `.gitignore`

- [ ] **Step 1.1: Create test scaffold**

Write to `workflow/tests/test_recheck_dispatch.py`:

```python
"""Tests for recheck_dispatch.py fire-log helpers (Issue #453)."""
import json
import sys
from pathlib import Path

import pytest

# Make .claude/hooks importable
_HOOKS_DIR = Path(__file__).resolve().parents[2] / ".claude" / "hooks"
sys.path.insert(0, str(_HOOKS_DIR))


def test_placeholder():
    """Removed in Task 2; here to verify test discovery + import path setup."""
    import recheck_dispatch
    assert hasattr(recheck_dispatch, "dispatch")
```

- [ ] **Step 1.2: Append .gitignore line**

Run:
```bash
echo ".claude/hook_fires.jsonl" >> .gitignore
```

Verify:
```bash
tail -3 .gitignore
```
Expected last line: `.claude/hook_fires.jsonl`.

- [ ] **Step 1.3: Verify pytest discovers + placeholder passes**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_recheck_dispatch.py -v
```
Expected: `test_placeholder PASSED`.

- [ ] **Step 1.4: Commit**

```bash
git add workflow/tests/test_recheck_dispatch.py .gitignore
git commit -m "$(cat <<'EOF'
test(pm): scaffold recheck_dispatch test file + gitignore hook_fires.jsonl

Issue #453 Task 1.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: `_log_fire` helper (TDD)

**Files:**
- Modify: `workflow/tests/test_recheck_dispatch.py` (replace placeholder)
- Modify: `.claude/hooks/recheck_dispatch.py` (add LOG_PATH + _log_fire)

- [ ] **Step 2.1: Write failing test**

Replace the `test_placeholder` function in `workflow/tests/test_recheck_dispatch.py` with:

```python
def test_log_fire_jsonl_append(tmp_path, monkeypatch):
    """_log_fire appends one valid JSONL line per call; metadata kwargs preserved."""
    import recheck_dispatch

    log_path = tmp_path / "hook_fires.jsonl"
    monkeypatch.setattr(recheck_dispatch, "LOG_PATH", log_path)

    recheck_dispatch._log_fire("hook_A", issue=1)
    recheck_dispatch._log_fire("hook_A", issue=2, extra="meta")
    recheck_dispatch._log_fire("hook_A", issue=3)

    lines = log_path.read_text().splitlines()
    assert len(lines) == 3
    for line in lines:
        rec = json.loads(line)
        assert rec["hook"] == "hook_A"
        assert "ts" in rec
        assert "issue" in rec
    rec2 = json.loads(lines[1])
    assert rec2.get("extra") == "meta"
    assert rec2.get("issue") == 2
```

- [ ] **Step 2.2: Run to verify fail**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_recheck_dispatch.py::test_log_fire_jsonl_append -v
```
Expected: FAIL with `AttributeError: module 'recheck_dispatch' has no attribute 'LOG_PATH'` (or `_log_fire`).

- [ ] **Step 2.3: Implement `LOG_PATH` + `_log_fire`**

Open `.claude/hooks/recheck_dispatch.py`. Add to imports near the top (after the existing `from pathlib import Path`):

```python
from datetime import datetime, timezone
```

After the existing module-level constants block (around line 30, after `PARENT_STATUS_SCRIPT = ...`), insert:

```python
# ---------------------------------------------------------------------------
# Fire-log infrastructure (Issue #453)
# ---------------------------------------------------------------------------

LOG_PATH = Path(__file__).resolve().parent.parent.parent / ".claude" / "hook_fires.jsonl"


def _log_fire(hook_name: str, issue: int | None = None, **metadata) -> None:
    """Append one JSONL line to LOG_PATH recording a hook fire.

    POSIX line-atomic under O_APPEND for lines <= PIPE_BUF (4096 B). Guards
    against oversized lines by skipping the write rather than risking
    interleaved-byte corruption.
    """
    payload = {
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "hook": hook_name,
        "issue": issue,
        **metadata,
    }
    line = json.dumps(payload, separators=(",", ":")) + "\n"
    if len(line.encode("utf-8")) >= 4096:
        print(f"hook_fires: oversized line (>=4KB) for {hook_name}, skipped", file=sys.stderr)
        return
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(line)
```

- [ ] **Step 2.4: Run to verify pass**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_recheck_dispatch.py::test_log_fire_jsonl_append -v
```
Expected: PASS.

- [ ] **Step 2.5: Commit**

```bash
git add workflow/tests/test_recheck_dispatch.py .claude/hooks/recheck_dispatch.py
git commit -m "$(cat <<'EOF'
feat(pm): _log_fire JSONL writer with oversized-line guard

Issue #453 Task 2. POSIX line-atomic O_APPEND for lines <= PIPE_BUF; guards
oversized lines by stderr-warning + skip rather than risking byte interleaving.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: `_count_fires` helper (TDD)

**Files:**
- Modify: `workflow/tests/test_recheck_dispatch.py` (add test)
- Modify: `.claude/hooks/recheck_dispatch.py` (add _count_fires)

- [ ] **Step 3.1: Write failing test**

Append to `workflow/tests/test_recheck_dispatch.py`:

```python
def test_count_fires_filters_by_hook(tmp_path, monkeypatch):
    """_count_fires returns per-hook counts; unknown hook returns 0; missing log returns 0."""
    import recheck_dispatch

    log_path = tmp_path / "hook_fires.jsonl"
    monkeypatch.setattr(recheck_dispatch, "LOG_PATH", log_path)

    # Missing log → 0 for any hook
    assert recheck_dispatch._count_fires("hook_A") == 0

    recheck_dispatch._log_fire("hook_A", issue=1)
    recheck_dispatch._log_fire("hook_A", issue=2)
    recheck_dispatch._log_fire("hook_B", issue=3)

    assert recheck_dispatch._count_fires("hook_A") == 2
    assert recheck_dispatch._count_fires("hook_B") == 1
    assert recheck_dispatch._count_fires("hook_C") == 0
```

- [ ] **Step 3.2: Run to verify fail**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_recheck_dispatch.py::test_count_fires_filters_by_hook -v
```
Expected: FAIL with `AttributeError: module 'recheck_dispatch' has no attribute '_count_fires'`.

- [ ] **Step 3.3: Implement `_count_fires`**

In `.claude/hooks/recheck_dispatch.py`, append below `_log_fire`:

```python
def _count_fires(hook_name: str) -> int:
    """Count lines in LOG_PATH where hook_name appears in the JSONL record.

    Substring match against the serialized form (`"hook":"<name>"`); cheap
    at K=3 scale, no JSON parsing needed.
    """
    if not LOG_PATH.exists():
        return 0
    needle = f'"hook":"{hook_name}"'
    return sum(1 for line in LOG_PATH.read_text().splitlines() if needle in line)
```

- [ ] **Step 3.4: Run to verify pass**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_recheck_dispatch.py::test_count_fires_filters_by_hook -v
```
Expected: PASS.

- [ ] **Step 3.5: Commit**

```bash
git add workflow/tests/test_recheck_dispatch.py .claude/hooks/recheck_dispatch.py
git commit -m "$(cat <<'EOF'
feat(pm): _count_fires substring-match counter

Issue #453 Task 3. Naive substring match against serialized JSON; cheap at
K=3 scale, no sidecar count file.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: `HOOK_CONFIG` + `_threshold_prompt` (TDD)

**Files:**
- Modify: `workflow/tests/test_recheck_dispatch.py` (add test)
- Modify: `.claude/hooks/recheck_dispatch.py` (add HOOK_CONFIG + _threshold_prompt)

- [ ] **Step 4.1: Write failing test**

Append to `workflow/tests/test_recheck_dispatch.py`:

```python
def test_threshold_prompt_equals_once(tmp_path, monkeypatch):
    """_threshold_prompt fires exactly at count == K, never before or after."""
    import recheck_dispatch

    log_path = tmp_path / "hook_fires.jsonl"
    monkeypatch.setattr(recheck_dispatch, "LOG_PATH", log_path)
    monkeypatch.setattr(recheck_dispatch, "HOOK_CONFIG", {
        "test_hook": {"threshold": 3, "dock": 999},
    })

    recheck_dispatch._log_fire("test_hook", issue=1)
    assert recheck_dispatch._threshold_prompt("test_hook") is None

    recheck_dispatch._log_fire("test_hook", issue=2)
    assert recheck_dispatch._threshold_prompt("test_hook") is None

    recheck_dispatch._log_fire("test_hook", issue=3)
    prompt = recheck_dispatch._threshold_prompt("test_hook")
    assert prompt is not None
    assert "test_hook" in prompt
    assert "3 times" in prompt
    assert "999" in prompt
    assert "check_hook_health.sh" in prompt

    recheck_dispatch._log_fire("test_hook", issue=4)
    # == K, not >=; 4 fires no longer matches.
    assert recheck_dispatch._threshold_prompt("test_hook") is None


def test_threshold_prompt_no_dock(tmp_path, monkeypatch):
    """_threshold_prompt emits '(no dock Issue filed yet)' when dock is None."""
    import recheck_dispatch

    log_path = tmp_path / "hook_fires.jsonl"
    monkeypatch.setattr(recheck_dispatch, "LOG_PATH", log_path)
    monkeypatch.setattr(recheck_dispatch, "HOOK_CONFIG", {
        "dockless_hook": {"threshold": 1, "dock": None},
    })

    recheck_dispatch._log_fire("dockless_hook", issue=1)
    prompt = recheck_dispatch._threshold_prompt("dockless_hook")
    assert prompt is not None
    assert "no dock Issue filed yet" in prompt
```

- [ ] **Step 4.2: Run to verify fail**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_recheck_dispatch.py -v
```
Expected: 2 new tests FAIL with `AttributeError: module 'recheck_dispatch' has no attribute 'HOOK_CONFIG'` (or `_threshold_prompt`).

- [ ] **Step 4.3: Implement `HOOK_CONFIG` + `_threshold_prompt`**

In `.claude/hooks/recheck_dispatch.py`, append below `_count_fires`:

```python
HOOK_CONFIG: dict[str, dict] = {
    "target_sync_check":     {"threshold": 3, "dock": 454},
    "recheck_milestone":     {"threshold": 3, "dock": None},   # dock Issue to be filed when promotion is sought
    "recheck_parent_status": {"threshold": 3, "dock": None},   # dock Issue to be filed when promotion is sought
}


def _threshold_prompt(hook_name: str) -> str | None:
    """Return the 🎯 promotion-review prompt when count == K (exactly); else None.

    Equals-once (not >=) so the prompt fires on exactly the K-th fire and
    never again. No nagging.
    """
    cfg = HOOK_CONFIG.get(hook_name)
    if cfg is None:
        return None
    if _count_fires(hook_name) != cfg["threshold"]:
        return None
    dock = f"Issue #{cfg['dock']}" if cfg["dock"] else "(no dock Issue filed yet)"
    return (
        f"🎯 {hook_name} has now fired {cfg['threshold']} times — review for promotion: "
        f"bash scripts/check_hook_health.sh ; {dock}"
    )
```

- [ ] **Step 4.4: Run to verify pass**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_recheck_dispatch.py -v
```
Expected: all 4 tests PASS.

- [ ] **Step 4.5: Commit**

```bash
git add workflow/tests/test_recheck_dispatch.py .claude/hooks/recheck_dispatch.py
git commit -m "$(cat <<'EOF'
feat(pm): HOOK_CONFIG + _threshold_prompt with equals-once semantics

Issue #453 Task 4. K=3 default per-hook; dock Issue ref per hook
(target_sync_check → #454; the other two → None pending promotion).
Equals-once (not >=) means the 🎯 line fires exactly once per hook
lifetime — no nagging.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: `_is_fire` + `_wrap_warning` + dispatcher wiring

**Files:**
- Modify: `.claude/hooks/recheck_dispatch.py` (add helpers + refactor dispatch() call-sites)

- [ ] **Step 5.1: Add `_is_fire` + `_wrap_warning` helpers**

In `.claude/hooks/recheck_dispatch.py`, append below `_threshold_prompt`:

```python
def _is_fire(hook_name: str, output: str) -> bool:
    """Return True if this output represents a real warning (not 'no change')."""
    if hook_name == "target_sync_check":
        return bool(output)  # function returns None when no fire
    if hook_name in ("recheck_milestone", "recheck_parent_status"):
        # Both scripts emit "Status: [No change]" when no drift is detected.
        return "Status: [No change]" not in output
    return False


def _wrap_warning(hook_name: str, issue: int | None, warning: str, outputs: list[str]) -> None:
    """Append warning to outputs; if it's a real fire, log it + append threshold prompt."""
    outputs.append(warning)  # always emit so user sees recheck context
    if not _is_fire(hook_name, warning):
        return
    _log_fire(hook_name, issue=issue)
    prompt = _threshold_prompt(hook_name)
    if prompt:
        outputs.append(prompt)
```

- [ ] **Step 5.2: Read existing dispatch() to identify call-sites**

Read `.claude/hooks/recheck_dispatch.py` lines 206-270 (the `dispatch()` function body) to inventory all `outputs.append(...)` sites that consume check-helper output. Expected call-sites (verified from current source):

- Line ~212: `outputs.append(f"[milestone recheck — close #{m.group(1)}]\n{run_recheck(...)}")` — wraps recheck_milestone fire
- Line ~215-216: `outputs.append(f"...{run_parent_status_recheck(...)}")` — wraps recheck_parent_status fire
- Line ~226-227, 232-233, 245-246, 263-264: more recheck_milestone fires
- Line ~255-256: another recheck_parent_status fire
- Line ~235-238: `sync_warn = target_sync_check(issue)` + conditional `outputs.append(sync_warn)`

- [ ] **Step 5.3: Refactor each call-site to use `_wrap_warning`**

For each `outputs.append(...)` involving a check's output, replace:

**Before (recheck_milestone close pattern):**
```python
outputs.append(f"[milestone recheck — close #{m.group(1)}]\n{run_recheck('--issue', m.group(1))}")
```

**After:**
```python
issue_n = int(m.group(1))
_wrap_warning(
    "recheck_milestone",
    issue_n,
    f"[milestone recheck — close #{issue_n}]\n{run_recheck('--issue', str(issue_n))}",
    outputs,
)
```

**Before (recheck_parent_status pattern):**
```python
outputs.append(f"{run_parent_status_recheck('--issue', m.group(1))}")
```

**After:**
```python
_wrap_warning(
    "recheck_parent_status",
    int(m.group(1)),
    run_parent_status_recheck('--issue', m.group(1)),
    outputs,
)
```

**Before (target_sync_check pattern around line 235):**
```python
sync_warn = target_sync_check(issue)
if sync_warn:
    outputs.append(sync_warn)
```

**After:**
```python
sync_warn = target_sync_check(issue)
if sync_warn:
    _wrap_warning("target_sync_check", issue, sync_warn, outputs)
```

(Note: keep the `if sync_warn:` guard because `target_sync_check` returns `None` for no-fire; `_wrap_warning` would otherwise append `None` to outputs.)

Apply this refactor to ALL call-sites identified in Step 5.2. Touch each call-site once; no missed sites.

- [ ] **Step 5.4: Smoke-test the refactored dispatcher**

Run a syntax + import check:
```bash
conda activate snakemake
python -c "import sys; sys.path.insert(0, '.claude/hooks'); import recheck_dispatch; print('OK:', dir(recheck_dispatch))" | tr ',' '\n' | grep -E "_log_fire|_count_fires|_threshold_prompt|_is_fire|_wrap_warning|HOOK_CONFIG"
```
Expected: 6 matches (all helpers + HOOK_CONFIG exposed).

Run a sample dispatch invocation against a fake PostToolUse event (the hook reads JSON from stdin):
```bash
echo '{"tool_input": {"command": "gh issue close 999999"}}' | python .claude/hooks/recheck_dispatch.py
```
Expected: exits cleanly; stdout may contain milestone-recheck output if #999999 exists (or empty / error block if not). No Python tracebacks.

- [ ] **Step 5.5: Run full test suite**

Run:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_recheck_dispatch.py -v
```
Expected: all 4 tests PASS (no regression from the refactor).

- [ ] **Step 5.6: Commit**

```bash
git add .claude/hooks/recheck_dispatch.py
git commit -m "$(cat <<'EOF'
feat(pm): _is_fire predicate + _wrap_warning + dispatch() refactor

Issue #453 Task 5. Refactors every check call-site in dispatch() to flow
through _wrap_warning, which logs fires + appends the threshold-crossing
prompt automatically. Scripts in scripts/pm/ untouched (predicate
sentinel-matches "Status: [No change]" for the 2 subprocess scripts).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: `scripts/check_hook_health.sh` aggregator

**Files:**
- Create: `scripts/check_hook_health.sh`

- [ ] **Step 6.1: Write the aggregator script**

Write to `scripts/check_hook_health.sh`:

```bash
#!/usr/bin/env bash
# Summarize fires logged to .claude/hook_fires.jsonl per hook.
# Operator-facing tool for hook-promotion decisions (Issue #453).
set -euo pipefail

LOG_PATH=".claude/hook_fires.jsonl"

# Per-hook dock Issue mapping. Duplicated from HOOK_CONFIG in
# .claude/hooks/recheck_dispatch.py — keep in sync manually (3 entries).
declare -A DOCK_MAP=(
    [target_sync_check]="Issue #454"
    [recheck_milestone]="(no dock Issue filed yet)"
    [recheck_parent_status]="(no dock Issue filed yet)"
)

HOOKS=(target_sync_check recheck_milestone recheck_parent_status)

echo "Hook health summary (${LOG_PATH})"
echo ""

if [[ ! -f "${LOG_PATH}" ]]; then
    echo "No fires logged yet."
    echo ""
    echo "Wired hooks:"
    for hook in "${HOOKS[@]}"; do
        echo "  - ${hook}  (dock: ${DOCK_MAP[$hook]})"
    done
    exit 0
fi

total=0
for hook in "${HOOKS[@]}"; do
    # Filter JSONL lines for this hook (one JSON object per line).
    count=$(jq -c --arg h "${hook}" 'select(.hook == $h)' "${LOG_PATH}" 2>/dev/null | wc -l | tr -d ' ')
    echo "${hook}"
    echo "  Fires:  ${count}"
    if [[ "${count}" -gt 0 ]]; then
        first=$(jq -r --arg h "${hook}" 'select(.hook == $h) | .ts' "${LOG_PATH}" | head -1)
        last=$(jq -r --arg h "${hook}" 'select(.hook == $h) | .ts' "${LOG_PATH}" | tail -1)
        echo "  First:  ${first}"
        echo "  Last:   ${last}"
    fi
    echo "  Dock:   ${DOCK_MAP[$hook]}"
    echo ""
    total=$((total + count))
done

echo "Total: ${total} fires across ${#HOOKS[@]} hooks"
```

- [ ] **Step 6.2: chmod +x**

Run:
```bash
chmod +x scripts/check_hook_health.sh
```

- [ ] **Step 6.3: Manual smoke — empty log case**

Run (after ensuring `.claude/hook_fires.jsonl` does not exist):
```bash
rm -f .claude/hook_fires.jsonl
bash scripts/check_hook_health.sh
```
Expected output:
```
Hook health summary (.claude/hook_fires.jsonl)

No fires logged yet.

Wired hooks:
  - target_sync_check  (dock: Issue #454)
  - recheck_milestone  (dock: (no dock Issue filed yet))
  - recheck_parent_status  (dock: (no dock Issue filed yet))
```

- [ ] **Step 6.4: Manual smoke — populated log case**

Seed a mock log:
```bash
cat > .claude/hook_fires.jsonl <<'EOF'
{"ts":"2026-05-21T10:14:33Z","hook":"target_sync_check","issue":300}
{"ts":"2026-05-25T11:00:00Z","hook":"target_sync_check","issue":301}
{"ts":"2026-05-28T14:25:47Z","hook":"target_sync_check","issue":302}
{"ts":"2026-05-13T09:02:11Z","hook":"recheck_milestone","issue":250}
{"ts":"2026-05-28T15:42:11Z","hook":"recheck_milestone","issue":533}
EOF
bash scripts/check_hook_health.sh
```
Expected output shape:
```
Hook health summary (.claude/hook_fires.jsonl)

target_sync_check
  Fires:  3
  First:  2026-05-21T10:14:33Z
  Last:   2026-05-28T14:25:47Z
  Dock:   Issue #454

recheck_milestone
  Fires:  2
  First:  2026-05-13T09:02:11Z
  Last:   2026-05-28T15:42:11Z
  Dock:   (no dock Issue filed yet)

recheck_parent_status
  Fires:  0
  Dock:   (no dock Issue filed yet)

Total: 5 fires across 3 hooks
```

Then clean up the mock log:
```bash
rm -f .claude/hook_fires.jsonl
```

- [ ] **Step 6.5: Commit**

```bash
git add scripts/check_hook_health.sh
git commit -m "$(cat <<'EOF'
feat(pm): scripts/check_hook_health.sh aggregator for promotion decisions

Issue #453 Task 6. Pure bash + jq; reads .claude/hook_fires.jsonl,
summarizes per-hook fire counts + first/last timestamps + dock Issue
references. Handles "no fires yet" cleanly.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Personas-repo shared memory edit

**Files (in personas repo, NOT project repo):**
- Modify: `~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/shared/feedback_hook_proves_out.md`

**Per shared-rule "personas-repo git state is not your responsibility":** edit the file; do NOT commit/push the personas repo. The user manages personas-repo git lifecycle outside this session.

- [ ] **Step 7.1: Read current state of feedback_hook_proves_out.md**

Resolve the file via the role-clone symlink:
```bash
ls -la .claude/memory/shared/feedback_hook_proves_out.md
```
Expected: symlink target points into the personas repo. Read it:

```bash
cat .claude/memory/shared/feedback_hook_proves_out.md
```

- [ ] **Step 7.2: Add Issue #453 cross-ref + "How to apply" line**

Edit the file in-place (via the symlink — bare relative path, no `git -C`). Append a new section at the bottom (or update existing if already structured around hook lifecycle):

```markdown
## Landed infrastructure ([Issue #453](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/453))

The "proves out" gate is now backed by deterministic instrumentation:

- `.claude/hook_fires.jsonl` (gitignored) — append-only fire log
- `scripts/check_hook_health.sh` — per-hook fire-count aggregator
- Threshold-crossing prompt fires at K=3 (configurable per hook in `HOOK_CONFIG`)

**How to apply:** before any hook-promotion decision (`.claude/settings.local.json` → committed `.claude/settings.json`), run `bash scripts/check_hook_health.sh` from the project root. Use the per-hook fire count + dock-Issue context to inform the (a) promote / (b) keep-local / (c) retire decision.
```

- [ ] **Step 7.3: Verify edit applied**

```bash
tail -15 .claude/memory/shared/feedback_hook_proves_out.md
```
Expected: new "Landed infrastructure" section visible.

(No commit. Personas-repo git is outside session scope.)

---

## Task 8: PR + merge

**Files:** none modified — git plumbing + GitHub interactions.

- [ ] **Step 8.1: Final pre-push test sweep**

Run the full test file one more time:
```bash
workflow/tests/.venv/bin/python -m pytest workflow/tests/test_recheck_dispatch.py -v
```
Expected: all 4 tests PASS.

Verify the gitignore line is present:
```bash
grep "hook_fires.jsonl" .gitignore
```
Expected: `.claude/hook_fires.jsonl`.

Verify the shell script is executable:
```bash
ls -la scripts/check_hook_health.sh
```
Expected: `-rwxr-xr-x` permissions.

- [ ] **Step 8.2: Push branch**

```bash
git push -u origin feat/pm/issue-453-hook-fire-log-infra
```
Expected: branch lands on remote, 7 commits ahead of main (Tasks 1-6 + the spec doc commit from brainstorming).

- [ ] **Step 8.3: Open PR**

```bash
gh pr create \
  --title "feat(pm): hook fire-log infra + threshold-crossing promotion prompt (closes #453)" \
  --body "$(cat <<'EOF'
## Summary

Replaces the vibes-check gate for hook promotion (`.claude/settings.local.json` → committed `.claude/settings.json`) with deterministic instrumentation. Closes [Issue #453](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/453).

## What landed

- **`_log_fire(hook_name, issue, **metadata)`** — POSIX line-atomic JSONL append to `.claude/hook_fires.jsonl` (gitignored). Oversized-line guard prevents corruption if schema grows past 4KB.
- **`_count_fires(hook_name)`** — substring-match counter; cheap at K=3 scale.
- **`HOOK_CONFIG`** — per-hook threshold + dock Issue dict; target_sync_check → #454, the other 2 → None pending promotion.
- **`_threshold_prompt(hook_name)`** — fires the 🎯 promotion-review line exactly once at count == K (not `>=`).
- **`_is_fire(hook_name, output)`** — per-hook fire predicate; matches `Status: [No change]` sentinel for the 2 subprocess scripts; None-vs-str for the inline `target_sync_check`. Scripts in `scripts/pm/` stay untouched.
- **`_wrap_warning(...)`** — wiring helper called from `dispatch()` at every check call-site.
- **`scripts/check_hook_health.sh`** — operator-facing aggregator.
- **4 unit tests** in `workflow/tests/test_recheck_dispatch.py` covering JSONL append, count filtering, threshold equals-once, dockless threshold prompt.
- **`.gitignore`** — appends `.claude/hook_fires.jsonl`.

## Companion change (personas repo, NOT in this PR diff)

- `shared/feedback_hook_proves_out.md` — adds "Landed infrastructure" section citing this Issue + "How to apply" pointing operators at `scripts/check_hook_health.sh` before any promotion decision. Lands via direct commit on the personas repo (managed by user, not this session).

## Design

Spec: [`docs/superpowers/specs/2026-05-28-hook-fire-log-infra-design.md`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/blob/feat/pm/issue-453-hook-fire-log-infra/docs/superpowers/specs/2026-05-28-hook-fire-log-infra-design.md)

## Test plan

- [ ] `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_recheck_dispatch.py -v` — 4/4 pass
- [ ] `bash scripts/check_hook_health.sh` on empty log → prints "No fires logged yet." + wired hooks
- [ ] `bash scripts/check_hook_health.sh` on mock-seeded log → per-hook summary + total
- [ ] Sample dispatcher invocation via stdin JSON → no Python tracebacks
- [ ] CI pipeline checks pass (ci-tools-pytest + pipeline-snakemake-dry-run + pipeline-pytest)

**Created by:** PM
EOF
)" \
  --label "role:pm" \
  --project "JH M Lee Lab"
```

- [ ] **Step 8.4: Flip Status to "In review" on both PR + Issue**

Per shared rule: PR + Issue both flip to `In review` in the same step as requesting review.

Get the PR number from Step 8.3 output, then:
```bash
PR_NUM=<from gh pr create output>
# Find the project item IDs for the PR and the Issue
PR_ITEM_ID=$(gh api graphql -f query="query { repository(owner: \"Jin-HoMLee\", name: \"splice-neoepitope-pipeline\") { pullRequest(number: ${PR_NUM}) { projectItems(first: 5) { nodes { id project { number } } } } } }" --jq ".data.repository.pullRequest.projectItems.nodes[] | select(.project.number == 9) | .id")
ISSUE_ITEM_ID=$(gh api graphql -f query='query { repository(owner: "Jin-HoMLee", name: "splice-neoepitope-pipeline") { issue(number: 453) { projectItems(first: 5) { nodes { id project { number } } } } } }' --jq '.data.repository.issue.projectItems.nodes[] | select(.project.number == 9) | .id')

# Status field + "In review" option ID
STATUS_FIELD_ID="PVTSSF_lAHOB17eGc4BSomPzhAHFf8"
IN_REVIEW_OPTION_ID=$(gh api graphql -f query='query { node(id: "PVT_kwHOB17eGc4BSomP") { ... on ProjectV2 { field(name: "Status") { ... on ProjectV2SingleSelectField { options { id name } } } } } }' --jq '.data.node.field.options[] | select(.name == "In review") | .id')

# Flip both in one mutation
gh api graphql -f query="mutation { setPR: updateProjectV2ItemFieldValue(input: { projectId: \"PVT_kwHOB17eGc4BSomP\", itemId: \"${PR_ITEM_ID}\", fieldId: \"${STATUS_FIELD_ID}\", value: { singleSelectOptionId: \"${IN_REVIEW_OPTION_ID}\" } }) { projectV2Item { id } } setIssue: updateProjectV2ItemFieldValue(input: { projectId: \"PVT_kwHOB17eGc4BSomP\", itemId: \"${ISSUE_ITEM_ID}\", fieldId: \"${STATUS_FIELD_ID}\", value: { singleSelectOptionId: \"${IN_REVIEW_OPTION_ID}\" } }) { projectV2Item { id } } }"
```

- [ ] **Step 8.5: Post @claude review**

```bash
gh pr comment ${PR_NUM} --body "@claude review"
```

Wait for review to land (typical 1-3 min).

- [ ] **Step 8.6: Address review findings (if any)**

When the review lands:
- Read findings via `gh pr view ${PR_NUM} --comments --json comments --jq '.comments[-1].body'`
- For each finding, verify against current code, fix if technically sound (per `superpowers:receiving-code-review` skill), push the fix.
- After every push, CI re-runs. Wait for green before merge.
- If review has no fix-before-merge findings, proceed.

- [ ] **Step 8.7: Tick AC boxes on Issue #453 + PR body**

Before merging, tick all 9 AC boxes on [Issue #453](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/453) (closure-ritual gate). Also tick the Test plan boxes on PR #${PR_NUM} body that are satisfied.

```bash
# Tick Issue #453 ACs
gh issue view 453 --json body --jq '.body' > /tmp/issue453_body.md
sed -i '' 's/^- \[ \] /- [x] /' /tmp/issue453_body.md
gh issue edit 453 --body-file /tmp/issue453_body.md

# Tick PR Test plan
gh pr view ${PR_NUM} --json body --jq '.body' > /tmp/pr_body.md
sed -i '' 's/^- \[ \] /- [x] /' /tmp/pr_body.md
gh pr edit ${PR_NUM} --body-file /tmp/pr_body.md
```

- [ ] **Step 8.8: Merge via closure-ritual gate**

```bash
bash scripts/audit_and_merge.sh ${PR_NUM}
```

Expected: `✓ PR #${PR_NUM} merged (N test-plan boxes ticked, 9 AC boxes ticked + 1/1 priority rationales present).`

- [ ] **Step 8.9: Sync local main + delete merged branch**

```bash
git checkout main
git pull origin main
git branch -D feat/pm/issue-453-hook-fire-log-infra
```

- [ ] **Step 8.10: Verify Issue #453 closed + Status auto-set to Done**

```bash
gh issue view 453 --json state,closedAt --jq .
```
Expected: `"state": "CLOSED"`.

---

## Wrap-up

After all 8 Tasks complete:

- [Issue #453](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/453) closed via PR.
- Fire-log infra live for the 3 sibling recheck hooks.
- `bash scripts/check_hook_health.sh` available for promotion-decision time.
- Threshold-crossing prompt will surface automatically on the 3rd fire of any wired hook.

**Follow-up considerations** (not part of this Issue):
- [Issue #454](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/454) (target-sync promotion) can now be picked up — its prerequisite (this infra) is satisfied.
- [Issue #533](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/533) (target_sync_check extract refactor) can be re-evaluated; this Issue shipping cleanly with the inline function strengthens the YAGNI case for closing #533 as wontfix.

**Lab notebook entry:** apply the "non-routine" criterion at execution time. Default: skip (single PR closes single Issue, body + commits + spec capture everything). Promote to required only if the implementation surfaces a behavioral decision that the spec/commits don't already capture.

---

## Self-review

**Spec coverage:**

| Spec section | Covered by |
|--------------|-----------|
| 3. Architecture (5 helpers + wiring + aggregator + .gitignore) | Tasks 2-6 |
| 4.1 Data shapes (JSONL + HOOK_CONFIG) | Tasks 2, 4 |
| 4.2 Concurrency + counting | Tasks 2, 3 |
| 4.3 Fire predicate + wiring | Task 5 |
| 4.4 Aggregator | Task 6 |
| 4.5 Testing (3 named tests) | Tasks 2, 3, 4 (4 tests — added test_threshold_prompt_no_dock for the dock=None branch beyond the spec's 3) |
| 4.6 Cross-cutting (.gitignore, shared memory) | Tasks 1, 7 |
| AC: `_log_fire` shared utility | Task 2 |
| AC: oversized-line / concurrent-safe append | Task 2 (>= 4KB guard) |
| AC: aggregator no-fires case | Task 6.3 |
| AC: equals-once threshold | Task 4 |
| AC: target_sync_check wired | Task 5.3 (sync_warn → _wrap_warning) |
| AC: milestone-capacity wired | Task 5.3 (run_recheck → _wrap_warning) |
| AC: parent-status wired | Task 5.3 (run_parent_status_recheck → _wrap_warning) |
| AC: unit tests for jsonl append + count + threshold | Tasks 2-4 |
| AC: .gitignore excludes jsonl | Task 1 |
| AC: shared memory cross-ref | Task 7 |

All 9 ACs covered. No spec section orphaned.

**Placeholder scan:** No "TBD", "TODO", "implement later" in the plan body. Every step has the actual content or command.

**Type/identifier consistency:**
- `LOG_PATH` (Path) — defined Task 2, monkeypatched in Tasks 2-4 tests ✓
- `HOOK_CONFIG` (dict[str, dict]) — defined Task 4, monkeypatched in Task 4 tests ✓
- Hook names used as keys: `target_sync_check`, `recheck_milestone`, `recheck_parent_status` — consistent across HOOK_CONFIG (Task 4), _is_fire (Task 5), dispatcher call-sites (Task 5), and aggregator (Task 6) ✓
- `_log_fire(hook_name, issue, **metadata)` signature — matches across Tasks 2 (def), 3 (test calls), 4 (test calls), 5 (`_wrap_warning` call site) ✓
- `_count_fires(hook_name)` signature — matches Tasks 3 (def), 4 (`_threshold_prompt` call) ✓
- `_threshold_prompt(hook_name) → str | None` — matches Tasks 4 (def), 5 (`_wrap_warning` call) ✓
- Variable `outputs` (list[str]) — matches existing `dispatch()` convention; passed to `_wrap_warning` ✓
