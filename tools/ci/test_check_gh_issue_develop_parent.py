"""Tests for `.agents/hooks/check_gh_issue_develop_parent.py` (Issue #549).

Pure-function unit tests import the module directly (no network). Orchestration
tests monkeypatch the two `gh` I/O functions (`resolve_repo`, `sub_issue_total`)
so `main()`'s deny/allow decision is exercised without a live `gh` call. The
subprocess fail-open tests mirror the harness invocation (pipe PreToolUse JSON to
stdin) and only travel paths that return before any `gh` call.

The full live parent-refused / leaf-allowed path is an integration check and is
verified by smoke-testing the hook against real Issues (#538 parent, #549 leaf),
documented in the PR Test plan — not in this suite.
"""

import io
import json
import subprocess
import sys
from pathlib import Path

HOOKS_DIR = Path(__file__).parent.parent.parent / ".agents" / "hooks"
HOOK = HOOKS_DIR / "check_gh_issue_develop_parent.py"
sys.path.insert(0, str(HOOKS_DIR))

import check_gh_issue_develop_parent as h  # noqa: E402


# --- develop_args: command matching ---


class TestDevelopArgs:
    def test_plain_invocation(self):
        assert h.develop_args("gh issue develop 549") == ["549"]

    def test_with_flags(self):
        assert h.develop_args(
            "gh issue develop 549 --name foo --base main"
        ) == ["549", "--name", "foo", "--base", "main"]

    def test_after_separator(self):
        assert h.develop_args("git fetch origin && gh issue develop 538") == ["538"]

    def test_trailing_separator_stops_collection(self):
        # tokens after the `&&` belong to the next command, not develop
        assert h.develop_args("gh issue develop 549 && echo done") == ["549"]

    def test_issue_list_not_matched(self):
        assert h.develop_args("gh issue list") is None

    def test_issue_view_not_matched(self):
        assert h.develop_args("gh issue view 549") is None

    def test_develop_inside_quoted_body_not_matched(self):
        # A PR-comment body discussing the command must NOT match.
        cmd = 'gh pr comment 1 --body "do not run gh issue develop 538"'
        assert h.develop_args(cmd) is None

    def test_separator_inside_quoted_body_not_matched(self):
        cmd = 'gh pr comment 1 --body "example: x && gh issue develop 538"'
        assert h.develop_args(cmd) is None

    def test_unbalanced_quotes_fail_safe(self):
        assert h.develop_args('gh issue develop 549 --name "oops') is None


# --- issue_number: selector extraction ---


class TestIssueNumber:
    def test_bare_number(self):
        assert h.issue_number(["549"]) == 549

    def test_number_after_value_flag(self):
        assert h.issue_number(["--name", "foo", "549"]) == 549

    def test_equals_form_value_flag_self_contained(self):
        assert h.issue_number(["--name=foo", "549"]) == 549

    def test_boolean_flag_not_treated_as_value_flag(self):
        # -c/--checkout takes no value; the number is the next token
        assert h.issue_number(["-c", "549"]) == 549

    def test_mixed_flags(self):
        assert h.issue_number(["--base", "main", "-c", "549"]) == 549

    def test_issue_url(self):
        url = "https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/538"
        assert h.issue_number([url]) == 538

    def test_no_positional_returns_none(self):
        assert h.issue_number(["--list"]) is None

    def test_empty_returns_none(self):
        assert h.issue_number([]) is None

    def test_unrecognizable_positional_returns_none(self):
        assert h.issue_number(["not-a-number"]) is None


# --- repo_from_args: explicit -R/--repo parsing ---


class TestRepoFromArgs:
    def test_short_flag(self):
        assert h.repo_from_args(["549", "-R", "owner/repo"]) == ("owner", "repo")

    def test_long_flag(self):
        assert h.repo_from_args(["--repo", "owner/repo", "549"]) == ("owner", "repo")

    def test_strips_host_segment(self):
        assert h.repo_from_args(
            ["--repo", "github.com/owner/repo", "549"]
        ) == ("owner", "repo")

    def test_equals_form(self):
        assert h.repo_from_args(["--repo=owner/repo"]) == ("owner", "repo")

    def test_equals_form_with_host_prefix(self):
        assert h.repo_from_args(
            ["--repo=github.com/owner/repo"]
        ) == ("owner", "repo")

    def test_absent_returns_none(self):
        assert h.repo_from_args(["549"]) is None


# --- deny_payload: shape ---


def test_deny_payload_shape():
    payload = h.deny_payload(538, 5)
    out = payload["hookSpecificOutput"]
    assert out["hookEventName"] == "PreToolUse"
    assert out["permissionDecision"] == "deny"
    reason = out["permissionDecisionReason"]
    assert "#538" in reason and "5 sub-issue" in reason
    assert "feedback_parent_sub_issues.md" in reason


# --- main(): orchestration (monkeypatched gh I/O, no network) ---


def _run_main(monkeypatch, capsys, cmd, *, repo=("Jin-HoMLee", "repo"), total):
    """Drive main() with stubbed gh I/O; return (rc, parsed_stdout_or_None)."""
    monkeypatch.setattr(h, "resolve_repo", lambda args: repo)
    monkeypatch.setattr(h, "sub_issue_total", lambda o, n, num: total)
    payload = json.dumps({"tool_input": {"command": cmd}})
    monkeypatch.setattr(sys, "stdin", io.StringIO(payload))
    rc = h.main()
    captured = capsys.readouterr().out.strip()
    return rc, (json.loads(captured) if captured else None)


def test_main_denies_parent(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(h, "LOG_PATH", tmp_path / "hook_fires.jsonl")
    rc, out = _run_main(monkeypatch, capsys, "gh issue develop 538", total=5)
    assert rc == 0
    assert out["hookSpecificOutput"]["permissionDecision"] == "deny"
    assert "#538" in out["hookSpecificOutput"]["permissionDecisionReason"]
    # fire-log line written on deny
    logged = (tmp_path / "hook_fires.jsonl").read_text().strip()
    rec = json.loads(logged)
    assert rec["hook"] == "check_gh_issue_develop_parent" and rec["issue"] == 538


def test_main_allows_leaf(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(h, "LOG_PATH", tmp_path / "hook_fires.jsonl")
    rc, out = _run_main(monkeypatch, capsys, "gh issue develop 549", total=0)
    assert rc == 0 and out is None
    # no fire-log on allow
    assert not (tmp_path / "hook_fires.jsonl").exists()


def test_main_fails_open_on_query_error(monkeypatch, capsys):
    # sub_issue_total returns None (gh/network error) → allow
    rc, out = _run_main(monkeypatch, capsys, "gh issue develop 538", total=None)
    assert rc == 0 and out is None


def test_main_fails_open_on_unresolvable_repo(monkeypatch, capsys):
    rc, out = _run_main(monkeypatch, capsys, "gh issue develop 538",
                        repo=None, total=5)
    assert rc == 0 and out is None


# --- subprocess fail-open / no-op (no live gh reached) ---


def _run_subprocess(stdin_payload: str):
    proc = subprocess.run(
        [sys.executable, str(HOOK)],
        input=stdin_payload,
        capture_output=True,
        text=True,
        timeout=10,
    )
    return proc.returncode, proc.stdout, proc.stderr


def _payload(command: str) -> str:
    return json.dumps({"tool_input": {"command": command}})


class TestSubprocessNoOp:
    def test_empty_stdin_fail_open(self):
        rc, out, _ = _run_subprocess("")
        assert rc == 0 and out.strip() == ""

    def test_malformed_json_fail_open(self):
        rc, out, _ = _run_subprocess("{not valid json")
        assert rc == 0 and out.strip() == ""

    def test_non_develop_command_noop(self):
        # returns before any gh call (develop_args is None)
        rc, out, _ = _run_subprocess(_payload("gh issue view 549"))
        assert rc == 0 and out.strip() == ""

    def test_develop_without_number_noop(self):
        # no parseable selector → returns before any gh call
        rc, out, _ = _run_subprocess(_payload("gh issue develop --list"))
        assert rc == 0 and out.strip() == ""
