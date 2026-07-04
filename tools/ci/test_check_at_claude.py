"""Subprocess tests for `.agents/hooks/check_at_claude.py`.

Mirrors how the Claude Code harness invokes the hook: pipe PreToolUse JSON to
stdin, read deny decision from stdout. Exit code is always 0.
"""

import json
import subprocess
import sys
from pathlib import Path

HOOK = Path(__file__).parent.parent.parent / ".agents" / "hooks" / "check_at_claude.py"


def _run(stdin_payload):
    """Pipe `stdin_payload` (str) to the hook, return (returncode, stdout, stderr)."""
    proc = subprocess.run(
        [sys.executable, str(HOOK)],
        input=stdin_payload,
        capture_output=True,
        text=True,
        timeout=5,
    )
    return proc.returncode, proc.stdout, proc.stderr


def _payload(command: str) -> str:
    return json.dumps({"tool_input": {"command": command}})


def _is_deny(stdout: str) -> bool:
    if not stdout.strip():
        return False
    parsed = json.loads(stdout)
    return (
        parsed.get("hookSpecificOutput", {}).get("permissionDecision") == "deny"
    )


def _deny_reason(stdout: str) -> str:
    """The human-facing `permissionDecisionReason` string from a deny decision."""
    parsed = json.loads(stdout)
    return parsed["hookSpecificOutput"]["permissionDecisionReason"]


class TestAllowPaths:
    def test_non_gh_command_allowed(self):
        rc, out, _ = _run(_payload("ls -la"))
        assert rc == 0 and out.strip() == ""

    def test_gh_status_allowed(self):
        rc, out, _ = _run(_payload("gh pr status"))
        assert rc == 0 and out.strip() == ""

    def test_gh_comment_without_at_claude_allowed(self):
        rc, out, _ = _run(_payload('gh pr comment 1 --body "looks good"'))
        assert rc == 0 and out.strip() == ""

    def test_canonical_review_trigger_double_quoted_allowed(self):
        rc, out, _ = _run(_payload('gh pr comment 1 --body "@claude review"'))
        assert rc == 0 and not _is_deny(out)

    def test_canonical_review_trigger_single_quoted_allowed(self):
        rc, out, _ = _run(_payload("gh pr comment 1 --body '@claude review'"))
        assert rc == 0 and not _is_deny(out)


class TestDenyPaths:
    def test_stray_at_claude_in_comment_body_denied(self):
        rc, out, _ = _run(
            _payload('gh issue comment 1 --body "ping @claude for review"')
        )
        assert rc == 0 and _is_deny(out)

    def test_stray_at_claude_in_issue_create_denied(self):
        rc, out, _ = _run(
            _payload('gh issue create --title "x" --body "as @claude said"')
        )
        assert rc == 0 and _is_deny(out)

    def test_stray_at_claude_in_pr_edit_denied(self):
        # Coverage extension (#6): `gh pr edit --body` was unguarded before.
        rc, out, _ = _run(
            _payload('gh pr edit 1 --body "rewritten body mentioning @claude"')
        )
        assert rc == 0 and _is_deny(out)

    def test_stray_at_claude_in_issue_edit_denied(self):
        # Coverage extension (#6): `gh issue edit --body` was unguarded before.
        rc, out, _ = _run(
            _payload('gh issue edit 1 --body "edited @claude reference"')
        )
        assert rc == 0 and _is_deny(out)

    def test_canonical_plus_stray_denied(self):
        rc, out, _ = _run(
            _payload(
                'gh pr comment 1 --body "@claude review" && '
                'echo "also tagging @claude separately"'
            )
        )
        assert rc == 0 and _is_deny(out)


class TestDenyMessageGuidance:
    """#975: the deny message must surface BOTH escape hatches, so an agent
    denied while legitimately requesting a review can construct the correct
    command from the message alone, without reading the hook source.
    """

    STRAY = 'gh pr comment 1 --body "ping @claude for review"'

    def test_reason_names_canonical_review_trigger(self):
        _, out, _ = _run(_payload(self.STRAY))
        reason = _deny_reason(out)
        # the one allowed form, and how to invoke it
        assert "@claude review" in reason
        assert "--body" in reason

    def test_reason_names_zero_width_workaround(self):
        _, out, _ = _run(_payload(self.STRAY))
        reason = _deny_reason(out)
        assert "@-claude" in reason

    def test_reason_names_guarded_edit_subcommand(self):
        # the guard covers `edit` too; the message should not undersell its scope
        _, out, _ = _run(_payload(self.STRAY))
        assert "edit" in _deny_reason(out)


class TestFailOpen:
    """#2: malformed/empty stdin must fail open, not crash."""

    def test_empty_stdin_allowed(self):
        rc, out, _ = _run("")
        assert rc == 0 and out.strip() == ""

    def test_malformed_json_allowed(self):
        rc, out, _ = _run("{not valid json")
        assert rc == 0 and out.strip() == ""
