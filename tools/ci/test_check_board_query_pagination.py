"""Tests for `.claude/hooks/check_board_query_pagination.py` (Issue #717).

Pure-function unit tests import the module directly (no I/O — the guard does no
`gh`/network calls, only string inspection). Subprocess tests mirror the harness
invocation (pipe PreToolUse JSON to stdin) and exercise the full deny / allow /
fail-open paths end-to-end.

The trigger incident (a `first: 100` query reading "Ready is empty" off the
711-item board, 2026-06-12) is documented in the PR Test plan, not re-run here.
"""

import io
import json
import subprocess
import sys
from pathlib import Path

HOOKS_DIR = Path(__file__).parent.parent.parent / ".claude" / "hooks"
HOOK = HOOKS_DIR / "check_board_query_pagination.py"
sys.path.insert(0, str(HOOKS_DIR))

import check_board_query_pagination as h  # noqa: E402


# Representative query fragments (kept on one line where the shape allows).
_UNPAGINATED = (
    "gh api graphql -f query='{ user(login: \"Jin-HoMLee\") { "
    "projectV2(number: 9) { items(first: 100) { nodes { id } } } } }'"
)
_PAGINATED = (
    "gh api graphql -f query='query($after:String){ user(login:\"Jin-HoMLee\"){ "
    "projectV2(number:9){ items(first:100, after:$after){ "
    "pageInfo{ hasNextPage endCursor } nodes{ id } } } } }'"
)
_PER_ISSUE = (
    "gh api graphql -f query='{ repository(owner:\"Jin-HoMLee\",name:\"x\"){ "
    "issue(number:717){ projectItems(first:5){ nodes{ id } } } } }'"
)
_TOTALCOUNT = (
    "gh api graphql -f query='{ user(login:\"Jin-HoMLee\"){ "
    "projectV2(number:9){ items { totalCount } } } }'"
)


# --- api_args: command matching ---


class TestApiArgs:
    def test_plain_invocation(self):
        assert h.api_args("gh api graphql -f query='x'") == "graphql -f query=x"

    def test_after_separator(self):
        # a `gh api` after `&&` is still a command start
        assert h.api_args("git fetch && gh api user --jq .login") == "user --jq .login"

    def test_trailing_separator_stops_collection(self):
        # tokens after the `&&` belong to the next command, not `gh api`
        assert h.api_args("gh api user && echo done") == "user"

    def test_non_api_gh_not_matched(self):
        assert h.api_args("gh issue list --state open") is None

    def test_non_gh_not_matched(self):
        assert h.api_args("python scripts/board_open_items.py --status Ready") is None

    def test_api_inside_quoted_body_not_matched(self):
        # A PR-comment body discussing the command must NOT match.
        cmd = "gh pr comment 1 --body \"don't run gh api graphql projectV2 items(first:100)\""
        assert h.api_args(cmd) is None

    def test_unbalanced_quotes_fail_safe(self):
        assert h.api_args("gh api graphql -f query='oops") is None


# --- is_unpaginated_board_query: the dangerous-shape detector ---


class TestIsUnpaginatedBoardQuery:
    def test_dangerous_shape_true(self):
        assert h.is_unpaginated_board_query(h.api_args(_UNPAGINATED)) is True

    def test_paginated_allowed(self):
        # hasNextPage / endCursor / after: present → safe
        assert h.is_unpaginated_board_query(h.api_args(_PAGINATED)) is False

    def test_after_cursor_alone_allowed(self):
        args = "graphql -f query={ projectV2 { items(first: 50, after: $c) { id } } }"
        assert h.is_unpaginated_board_query(args) is False

    def test_declared_but_unused_after_variable_still_denied(self):
        # `$after:String` is a variable *declaration*, not cursor usage — the
        # `(?<!\$)` lookbehind keeps it from masking a genuinely unpaginated query.
        args = (
            "graphql -f query=query($after:String){ projectV2 { "
            "items(first: 100) { nodes { id } } } }"
        )
        assert h.is_unpaginated_board_query(args) is True

    def test_per_issue_projectitems_allowed(self):
        # no `projectV2` token → not a whole-board scan
        assert h.is_unpaginated_board_query(h.api_args(_PER_ISSUE)) is False

    def test_totalcount_only_allowed(self):
        # `items { totalCount }` has no `items(first:` → not the dangerous shape
        assert h.is_unpaginated_board_query(h.api_args(_TOTALCOUNT)) is False

    def test_whitespace_tolerant_items_first(self):
        args = "graphql -f query={ projectV2 { items( first : 100 ) { id } } }"
        assert h.is_unpaginated_board_query(args) is True

    def test_no_projectv2_allowed(self):
        args = "graphql -f query={ viewer { items(first: 100) { id } } }"
        assert h.is_unpaginated_board_query(args) is False


# --- deny_payload: shape ---


def test_deny_payload_shape():
    out = h.deny_payload()["hookSpecificOutput"]
    assert out["hookEventName"] == "PreToolUse"
    assert out["permissionDecision"] == "deny"
    reason = out["permissionDecisionReason"]
    assert "board_open_items.py" in reason
    assert "feedback_board_queries.md" in reason


# --- main(): orchestration via subprocess (matches harness; no monkeypatch needed) ---


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


def test_main_denies_unpaginated(tmp_path, monkeypatch):
    monkeypatch.setattr(h, "LOG_PATH", tmp_path / "hook_fires.jsonl")
    payload = json.dumps({"tool_input": {"command": _UNPAGINATED}})
    monkeypatch.setattr(sys, "stdin", io.StringIO(payload))
    rc = h.main()
    assert rc == 0
    # fire-log written on deny
    logged = (tmp_path / "hook_fires.jsonl").read_text().strip()
    rec = json.loads(logged)
    assert rec["hook"] == "check_board_query_pagination"


class TestSubprocess:
    def test_denies_unpaginated(self):
        rc, out, _ = _run_subprocess(_payload(_UNPAGINATED))
        assert rc == 0
        decision = json.loads(out)["hookSpecificOutput"]["permissionDecision"]
        assert decision == "deny"

    def test_allows_paginated(self):
        rc, out, _ = _run_subprocess(_payload(_PAGINATED))
        assert rc == 0 and out.strip() == ""

    def test_allows_per_issue(self):
        rc, out, _ = _run_subprocess(_payload(_PER_ISSUE))
        assert rc == 0 and out.strip() == ""

    def test_noop_non_api_command(self):
        rc, out, _ = _run_subprocess(_payload("gh issue list --state open --limit 1000"))
        assert rc == 0 and out.strip() == ""

    def test_empty_stdin_fail_open(self):
        rc, out, _ = _run_subprocess("")
        assert rc == 0 and out.strip() == ""

    def test_malformed_json_fail_open(self):
        rc, out, _ = _run_subprocess("{not valid json")
        assert rc == 0 and out.strip() == ""

    def test_comment_body_discussing_pattern_not_denied(self):
        # documenting the guard in a gh comment must not trip it
        cmd = "gh pr comment 1 --body \"use projectV2 items(first:100) carefully\""
        rc, out, _ = _run_subprocess(_payload(cmd))
        assert rc == 0 and out.strip() == ""
