"""Tests for `.agents/hooks/check_commit_push_separation.py` (Issue #626).

Pure-function unit tests import the module directly. Subprocess tests mirror the
harness invocation (pipe PreToolUse JSON to stdin) and exercise the deny / allow
/ escape-hatch / fail-open paths end-to-end.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

HOOKS_DIR = Path(__file__).parent.parent.parent / ".agents" / "hooks"
HOOK = HOOKS_DIR / "check_commit_push_separation.py"
sys.path.insert(0, str(HOOKS_DIR))

import check_commit_push_separation as h  # noqa: E402


# --- chained_violation ---


class TestChainedViolation:
    def test_commit_then_push(self):
        assert h.chained_violation("git commit -m x && git push") == "commit->push"

    def test_push_then_merge(self):
        assert h.chained_violation("git push && gh pr merge 42 --squash") == "push->merge"

    def test_commit_push_merge_all(self):
        # commit->push is detected first (both violations present)
        assert h.chained_violation("git commit -m x && git push && gh pr merge") == "commit->push"

    def test_semicolon_separator(self):
        assert h.chained_violation("git commit -m x ; git push") == "commit->push"

    def test_commit_alone_allowed(self):
        assert h.chained_violation("git commit -m x") is None

    def test_push_alone_allowed(self):
        assert h.chained_violation("git push origin main") is None

    def test_add_then_commit_allowed(self):
        assert h.chained_violation("git add -A && git commit -m x") is None

    def test_fetch_then_push_allowed(self):
        assert h.chained_violation("git fetch && git push") is None

    def test_push_then_commit_not_flagged(self):
        # reverse order is not the dangerous pattern
        assert h.chained_violation("git push && git commit -m x") is None

    def test_quoted_message_not_flagged(self):
        assert h.chained_violation('git commit -m "push and merge later"') is None

    def test_parse_miss_allows(self):
        assert h.chained_violation('git commit -m "x && git push') is None


# --- end-to-end via stdin ---


def _run(payload: dict, env=None):
    return subprocess.run(
        [sys.executable, str(HOOK)],
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        env=env,
    )


def _bash(cmd: str) -> dict:
    return {"tool_name": "Bash", "tool_input": {"command": cmd}}


class TestEndToEnd:
    def test_chained_denied(self):
        proc = _run(_bash("git commit -m x && git push"))
        assert proc.returncode == 0
        out = json.loads(proc.stdout)
        assert out["hookSpecificOutput"]["permissionDecision"] == "deny"

    def test_separate_steps_allowed(self):
        assert _run(_bash("git commit -m x")).stdout.strip() == ""
        assert _run(_bash("git push origin main")).stdout.strip() == ""

    def test_escape_hatch_allows(self):
        env = dict(os.environ, CLAUDE_ALLOW_CHAINED_GIT="1")
        assert _run(_bash("git commit -m x && git push"), env=env).stdout.strip() == ""

    def test_bad_json_fails_open(self):
        proc = subprocess.run(
            [sys.executable, str(HOOK)],
            input="not json",
            capture_output=True,
            text=True,
        )
        assert proc.returncode == 0
        assert proc.stdout.strip() == ""

    def test_non_bash_ignored(self):
        assert _run({"tool_name": "Edit", "tool_input": {"command": "git commit && git push"}}).stdout.strip() == ""


def test_hook_is_executable():
    # The harness execs the hook by bare path; a non-executable file EACCESes
    # before the shebang (the blocker caught on PR #1029).
    assert os.access(HOOK, os.X_OK), f"{HOOK} must be committed executable (100755)"


class TestNewlineForm:
    """Issue #1142: the `&&` chain was caught, the two-line form was not.

    A newline is a command separator in shell but not in `shlex`, so before the
    `normalize_command()` fix `git commit -m x` and `git push` written on two
    lines merged into ONE token block. The commit and push indices then collapsed
    to the same position, `min(commit) < max(push)` was `0 < 0` -> False, and the
    guard silently allowed the exact chain it exists to refuse. Only the `&&`
    spelling was ever enforced.
    """

    def test_newline_separated_commit_then_push_is_caught(self):
        assert h.chained_violation("git commit -m x\ngit push") is not None

    def test_ampersand_form_still_caught(self):
        assert h.chained_violation("git commit -m x && git push") is not None

    def test_commit_alone_still_allowed(self):
        assert h.chained_violation("git commit -m x") is None

    def test_push_alone_still_allowed(self):
        assert h.chained_violation("git push") is None
