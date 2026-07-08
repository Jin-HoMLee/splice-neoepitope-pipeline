"""Tests for `.agents/hooks/check_memory_path_cwd_drift.py` (Issue #1053).

Pure-function unit tests pass an explicit project_root / token so they don't
depend on the real clone location. Subprocess tests pipe PreToolUse JSON to stdin
against the real hook (which derives PROJECT_ROOT from its own path) and exercise
deny / allow / escape-hatch / fail-open end-to-end.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

HOOKS_DIR = Path(__file__).parent.parent.parent / ".agents" / "hooks"
HOOK = HOOKS_DIR / "check_memory_path_cwd_drift.py"
REAL_ROOT = str(HOOKS_DIR.parent.parent)  # the actual project clone root
sys.path.insert(0, str(HOOKS_DIR))

import check_memory_path_cwd_drift as h  # noqa: E402

ROOT = "/proj"


# --- cwd_outside_clone ---


class TestCwdOutsideClone:
    def test_exact_root_inside(self):
        assert h.cwd_outside_clone("/proj", "/proj") is False

    def test_subdir_inside(self):
        assert h.cwd_outside_clone("/proj/workflow", "/proj") is False

    def test_temp_root_outside(self):
        assert h.cwd_outside_clone("/private/tmp/claude/scratch", "/proj") is True

    def test_sibling_clone_outside(self):
        assert h.cwd_outside_clone("/proj-developer", "/proj") is True  # prefix but not subdir

    def test_missing_cwd_not_outside(self):
        assert h.cwd_outside_clone(None, "/proj") is False
        assert h.cwd_outside_clone("", "/proj") is False

    def test_trailing_slash_normalized(self):
        assert h.cwd_outside_clone("/proj/", "/proj") is False


# --- is_relative_memory_token ---


class TestIsRelativeMemoryToken:
    def test_relative_agents_memory(self):
        assert h.is_relative_memory_token(".agents/memory/shared/MEMORY.md") is True

    def test_relative_claude_memory(self):
        assert h.is_relative_memory_token(".claude/memory/foo.md") is True

    def test_dot_slash_prefix(self):
        assert h.is_relative_memory_token("./.agents/memory/x") is True

    def test_absolute_memory_safe(self):
        assert h.is_relative_memory_token("/Users/x/.agents/memory/y") is False

    def test_home_anchored_safe(self):
        assert h.is_relative_memory_token("~/.claude/memory/y") is False

    def test_non_memory_relative(self):
        assert h.is_relative_memory_token("workflow/scripts/x.py") is False

    def test_dynamic_fails_open(self):
        assert h.is_relative_memory_token("$HOME/.agents/memory/x") is False
        assert h.is_relative_memory_token("`pwd`/.agents/memory/x") is False

    def test_empty(self):
        assert h.is_relative_memory_token("") is False

    def test_memory_mid_path_not_at_start(self):
        # a token that merely contains the segment deeper in is not a root-relative ref
        assert h.is_relative_memory_token("foo/.agents/memory/x") is False


# --- command_has_relative_memory_ref ---


class TestCommandScan:
    def test_grep_relative_memory(self):
        assert h.command_has_relative_memory_ref("grep -r foo .agents/memory/shared") == ".agents/memory/shared"

    def test_absolute_memory_command_safe(self):
        assert h.command_has_relative_memory_ref("cat /abs/.agents/memory/x") is None

    def test_no_memory_ref(self):
        assert h.command_has_relative_memory_ref("git status && ls workflow") is None

    def test_quoted_token(self):
        assert h.command_has_relative_memory_ref("cat '.claude/memory/a b.md'") == ".claude/memory/a b.md"

    def test_tokenization_failure_fails_open(self):
        assert h.command_has_relative_memory_ref("cat '.agents/memory/x") is None  # unbalanced quote


# --- offending_ref (per tool) ---


class TestOffendingRef:
    def test_bash(self):
        assert h.offending_ref("Bash", {"command": "cat .agents/memory/x"}) == ".agents/memory/x"

    def test_edit_relative_filepath(self):
        assert h.offending_ref("Edit", {"file_path": ".agents/memory/x"}) == ".agents/memory/x"

    def test_edit_absolute_filepath_safe(self):
        assert h.offending_ref("Edit", {"file_path": "/abs/.agents/memory/x"}) is None

    def test_write_relative(self):
        assert h.offending_ref("Write", {"file_path": ".claude/memory/x"}) == ".claude/memory/x"

    def test_unknown_tool(self):
        assert h.offending_ref("Read", {"file_path": ".agents/memory/x"}) is None


# --- is_allowed ---


class TestIsAllowed:
    def test_set(self):
        assert h.is_allowed({"CLAUDE_ALLOW_CWD_DRIFT": "1"}) is True
        assert h.is_allowed({"CLAUDE_ALLOW_CWD_DRIFT": "true"}) is True

    def test_unset(self):
        assert h.is_allowed({}) is False
        assert h.is_allowed({"CLAUDE_ALLOW_CWD_DRIFT": "0"}) is False


# --- end-to-end via stdin (uses the real project root) ---


def _run(payload: dict, env=None):
    return subprocess.run(
        [sys.executable, str(HOOK)],
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        env=env,
    )


def _bash(cmd: str, cwd: str) -> dict:
    return {"tool_name": "Bash", "tool_input": {"command": cmd}, "cwd": cwd}


DRIFTED = "/private/tmp/claude-scratch"


class TestEndToEnd:
    def test_drift_plus_relative_memory_denied(self):
        proc = _run(_bash("grep foo .agents/memory/shared/MEMORY.md", DRIFTED))
        assert proc.returncode == 0
        out = json.loads(proc.stdout)
        assert out["hookSpecificOutput"]["permissionDecision"] == "deny"

    def test_drift_plus_absolute_memory_allowed(self):
        assert _run(_bash(f"grep foo {REAL_ROOT}/.agents/memory/x", DRIFTED)).stdout.strip() == ""

    def test_drift_plus_non_memory_allowed(self):
        assert _run(_bash("ls workflow", DRIFTED)).stdout.strip() == ""

    def test_no_drift_plus_relative_memory_allowed(self):
        assert _run(_bash("grep foo .agents/memory/x", REAL_ROOT)).stdout.strip() == ""

    def test_edit_relative_memory_from_drift_denied(self):
        payload = {"tool_name": "Edit", "tool_input": {"file_path": ".agents/memory/x"}, "cwd": DRIFTED}
        out = json.loads(_run(payload).stdout)
        assert out["hookSpecificOutput"]["permissionDecision"] == "deny"

    def test_escape_hatch(self):
        env = dict(os.environ, CLAUDE_ALLOW_CWD_DRIFT="1")
        assert _run(_bash("cat .agents/memory/x", DRIFTED), env=env).stdout.strip() == ""

    def test_missing_cwd_fails_open(self):
        payload = {"tool_name": "Bash", "tool_input": {"command": "cat .agents/memory/x"}}
        assert _run(payload).stdout.strip() == ""

    def test_bad_json_fails_open(self):
        proc = subprocess.run(
            [sys.executable, str(HOOK)], input="not json", capture_output=True, text=True
        )
        assert proc.returncode == 0
        assert proc.stdout.strip() == ""

    def test_non_guarded_tool_ignored(self):
        payload = {"tool_name": "Read", "tool_input": {"file_path": ".agents/memory/x"}, "cwd": DRIFTED}
        assert _run(payload).stdout.strip() == ""


def test_hook_is_executable():
    # The harness execs the hook by bare path; a non-executable file EACCESes
    # before the shebang (the blocker caught on PR #1029 / lesson #1032).
    assert os.access(HOOK, os.X_OK), f"{HOOK} must be committed executable (100755)"
