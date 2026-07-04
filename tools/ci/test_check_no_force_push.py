"""Tests for `.agents/hooks/check_no_force_push.py` (Issue #626).

Pure-function unit tests import the module directly (the guard does no I/O, only
string inspection). Subprocess tests mirror the harness invocation (pipe
PreToolUse JSON to stdin) and exercise the deny / allow / escape-hatch /
fail-open paths end-to-end.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

HOOKS_DIR = Path(__file__).parent.parent.parent / ".agents" / "hooks"
HOOK = HOOKS_DIR / "check_no_force_push.py"
sys.path.insert(0, str(HOOKS_DIR))

import check_no_force_push as h  # noqa: E402


# --- _is_force_flag / is_force_push ---


class TestForceFlagDetection:
    def test_long_force(self):
        assert h._is_force_flag("--force")

    def test_short_force(self):
        assert h._is_force_flag("-f")

    def test_short_bundle_with_f(self):
        assert h._is_force_flag("-fq")
        assert h._is_force_flag("-qf")

    def test_force_with_lease_bare(self):
        assert h._is_force_flag("--force-with-lease")

    def test_force_with_lease_ref(self):
        assert h._is_force_flag("--force-with-lease=origin/main")

    def test_force_if_includes(self):
        assert h._is_force_flag("--force-if-includes")

    def test_non_force_flags(self):
        for t in ("-u", "--verbose", "-v", "origin", "main", "--set-upstream", "-n", "-q"):
            assert not h._is_force_flag(t)

    def test_double_dash_foo_not_force(self):
        assert not h._is_force_flag("--foo")


class TestIsForcePush:
    def test_plain_force(self):
        assert h.is_force_push(["git", "push", "--force"])

    def test_force_flag_last(self):
        assert h.is_force_push(["git", "push", "origin", "main", "--force"])

    def test_short_f(self):
        assert h.is_force_push(["git", "push", "-f"])

    def test_lease(self):
        assert h.is_force_push(["git", "push", "--force-with-lease"])

    def test_push_no_force_allowed(self):
        assert not h.is_force_push(["git", "push", "origin", "main"])

    def test_non_push_with_force_token(self):
        # a force flag but no `push` subcommand -> not a force push
        assert not h.is_force_push(["git", "commit", "--force"])

    def test_requires_git_token(self):
        assert not h.is_force_push(["push", "--force"])


# --- split_subcommands ---


class TestSplitSubcommands:
    def test_single(self):
        assert h.split_subcommands("git push --force") == [["git", "push", "--force"]]

    def test_chained(self):
        subs = h.split_subcommands("git commit -m x && git push --force")
        assert subs == [["git", "commit", "-m", "x"], ["git", "push", "--force"]]

    def test_quoted_collapses(self):
        # the quoted message is ONE token; no bare push/--force leaks out
        subs = h.split_subcommands('git commit -m "push --force please"')
        assert subs == [["git", "commit", "-m", "push --force please"]]

    def test_unbalanced_quotes_fail_open(self):
        assert h.split_subcommands('git push "--force') is None


# --- command_force_pushes ---


class TestCommandForcePushes:
    def test_standalone_force(self):
        assert h.command_force_pushes("git push --force")

    def test_chained_after_commit(self):
        assert h.command_force_pushes("git commit -m x && git push -f")

    def test_plain_push_allowed(self):
        assert not h.command_force_pushes("git push origin main")

    def test_quoted_message_not_flagged(self):
        assert not h.command_force_pushes('git commit -m "push --force"')

    def test_force_in_other_subcommand_not_attributed(self):
        # --force belongs to a different subcommand than the push
        assert not h.command_force_pushes("git push && ./deploy --force")

    def test_parse_miss_allows(self):
        assert not h.command_force_pushes('git push "--force')


# --- end-to-end via stdin (subprocess) ---


def _run(payload: dict, env=None):
    proc = subprocess.run(
        [sys.executable, str(HOOK)],
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        env=env,
    )
    return proc


def _bash(cmd: str) -> dict:
    return {"tool_name": "Bash", "tool_input": {"command": cmd}}


class TestEndToEnd:
    def test_force_push_denied(self):
        proc = _run(_bash("git push --force origin main"))
        assert proc.returncode == 0
        out = json.loads(proc.stdout)
        assert out["hookSpecificOutput"]["permissionDecision"] == "deny"

    def test_plain_push_allowed(self):
        proc = _run(_bash("git push origin main"))
        assert proc.returncode == 0
        assert proc.stdout.strip() == ""

    def test_escape_hatch_allows(self):
        env = dict(os.environ, CLAUDE_ALLOW_FORCE_PUSH="1")
        proc = _run(_bash("git push --force"), env=env)
        assert proc.stdout.strip() == ""

    def test_bad_json_fails_open(self):
        proc = subprocess.run(
            [sys.executable, str(HOOK)],
            input="not json",
            capture_output=True,
            text=True,
        )
        assert proc.returncode == 0
        assert proc.stdout.strip() == ""

    def test_non_bash_tool_ignored(self):
        proc = _run({"tool_name": "Edit", "tool_input": {"command": "git push --force"}})
        assert proc.stdout.strip() == ""
