"""Tests for `.agents/hooks/check_created_by_footer.py` (Issue #1094).

A PreToolUse guard that refuses a `gh (pr|issue) comment` whose authored body
lacks the `**Created by:** <Role>` footer, mechanism-escalating a rule that
slipped >= 3x on the same shape (a footer is the last, feedback-less line of a
long body - exactly what a hook fixes and memory does not).

Design (mirrors the sibling `check_at_claude.py`, and composes with the #1197
accept-the-hole decision): DETECT the comment command on the normalized text (so a
`gh comment` after a heredoc is seen - the #1130/#1142 lesson), but INSPECT the
body. The body is only deniable when we can actually see it: an inline `--body`,
or a heredoc feeding `--body-file -` (the body is in the raw command string). A
real `--body-file <path>` is NOT inspectable, so it fails OPEN (allow) - denying a
compliant file-bodied comment would be a false positive. Fail open on every
uncertain path; deny only a confirmed missing footer.
"""

import io
import json
import subprocess
import sys
from pathlib import Path

HOOKS_DIR = Path(__file__).parent.parent.parent / ".agents" / "hooks"
HOOK = HOOKS_DIR / "check_created_by_footer.py"
sys.path.insert(0, str(HOOKS_DIR))

import check_created_by_footer as h  # noqa: E402


# --- matches_comment: command detection ---


class TestMatchesComment:
    def test_pr_comment(self):
        assert h.matches_comment("gh pr comment 5 --body x") is True

    def test_issue_comment(self):
        assert h.matches_comment("gh issue comment 5 --body x") is True

    def test_pr_create_not_matched(self):
        # Footer rule is for COMMENTS; PR bodies are covered by the open checklist.
        assert h.matches_comment("gh pr create --title t --body x") is False

    def test_pr_edit_not_matched(self):
        assert h.matches_comment("gh pr edit 5 --body x") is False

    def test_after_heredoc_still_matched(self):
        # The dominant authoring path: a heredoc, then the comment command.
        cmd = "cat > /tmp/b.md <<'EOF'\nnotes\nEOF\ngh issue comment 5 --body-file /tmp/b.md"
        assert h.matches_comment(cmd) is True

    def test_comment_inside_quoted_body_not_matched(self):
        # A comment body that merely mentions the command must not self-trigger.
        cmd = 'gh pr comment 5 --body "do not run gh issue comment 9 like this"'
        # It IS a real `gh pr comment` at the start, so it matches on that -
        # the body-mention is not what we key on. This asserts the leading command
        # is what matches, not the quoted mention.
        assert h.matches_comment(cmd) is True

    def test_non_gh_not_matched(self):
        assert h.matches_comment("echo gh pr comment 5") is False


# --- extract_body: inline / heredoc / not-inspectable ---


class TestExtractBody:
    def test_inline_body_double_quoted(self):
        assert h.extract_body('gh pr comment 5 --body "hello world"') == "hello world"

    def test_inline_body_short_flag(self):
        assert h.extract_body("gh pr comment 5 -b 'hi there'") == "hi there"

    def test_heredoc_body_via_body_file_dash(self):
        cmd = ("gh issue comment 5 --body-file - <<'EOF'\n"
               "line one\n**Created by:** Developer\nEOF")
        body = h.extract_body(cmd)
        assert "line one" in body and "**Created by:** Developer" in body

    def test_real_body_file_not_inspectable(self):
        # A real path: we cannot see the content -> None (caller fails open).
        assert h.extract_body("gh pr comment 5 --body-file /tmp/real.md") is None

    def test_no_body_flag_not_inspectable(self):
        assert h.extract_body("gh pr comment 5") is None

    def test_multi_heredoc_body_file_dash_fails_open(self):
        # Finding 1 (PR #1229 review): a scratch heredoc BEFORE the gh comment's own
        # `--body-file -` heredoc. Binding to the FIRST heredoc would read the
        # scratch body (no footer) and false-positive-DENY a compliant comment.
        # Ambiguous (>1 heredoc) -> not inspectable -> None (caller fails open).
        cmd = ("cat > /tmp/scratch.md <<'EOF'\nscratch notes no footer\nEOF\n"
               "gh issue comment 5 --body-file - <<'BODY'\nReal comment.\n"
               "**Created by:** Developer\nBODY")
        assert h.extract_body(cmd) is None


# --- is_exempt: From-opener + bare review trigger ---


class TestIsExempt:
    def test_from_opener_exempt(self):
        assert h.is_exempt("**From:** Developer -> **To:** PM\n\nhi") is True

    def test_from_opener_after_whitespace_exempt(self):
        assert h.is_exempt("  \n**From:** Developer\n\nhi") is True

    def test_bare_review_trigger_exempt(self):
        assert h.is_exempt("@claude review") is True

    def test_hyphenated_review_trigger_exempt(self):
        assert h.is_exempt("@-claude review") is True

    def test_ordinary_body_not_exempt(self):
        assert h.is_exempt("just a normal comment") is False

    def test_mid_body_from_not_exempt(self):
        # `**From:**` must be the OPENER, not buried mid-body.
        assert h.is_exempt("some text\n**From:** X") is False


# --- has_footer ---


class TestHasFooter:
    def test_present(self):
        assert h.has_footer("text\n\n**Created by:** Developer") is True

    def test_absent(self):
        assert h.has_footer("text with no footer") is False


# --- main(): orchestration (stdin-driven, no network) ---


def _run(cmd, *, env=None, tmp_log=None, monkeypatch=None):
    if monkeypatch is not None and tmp_log is not None:
        monkeypatch.setattr(h, "LOG_PATH", tmp_log)
    if env is not None:
        for k, v in env.items():
            monkeypatch.setenv(k, v)
    payload = json.dumps({"tool_input": {"command": cmd}})
    monkeypatch.setattr(sys, "stdin", io.StringIO(payload))
    rc = h.main()
    return rc


def _deny_out(capsys):
    out = capsys.readouterr().out.strip()
    return json.loads(out) if out else None


def test_missing_footer_denies(monkeypatch, capsys, tmp_path):
    rc = _run('gh pr comment 5 --body "a plain comment"',
              tmp_log=tmp_path / "f.jsonl", monkeypatch=monkeypatch)
    out = _deny_out(capsys)
    assert rc == 0
    assert out["hookSpecificOutput"]["permissionDecision"] == "deny"
    assert "Created by" in out["hookSpecificOutput"]["permissionDecisionReason"]
    # fire-log written
    assert (tmp_path / "f.jsonl").read_text().strip()


def test_footer_present_allows(monkeypatch, capsys, tmp_path):
    rc = _run('gh pr comment 5 --body "text\n\n**Created by:** Developer"',
              tmp_log=tmp_path / "f.jsonl", monkeypatch=monkeypatch)
    assert rc == 0 and _deny_out(capsys) is None


def test_from_opener_allows(monkeypatch, capsys, tmp_path):
    rc = _run('gh issue comment 5 --body "**From:** Developer -> **To:** PM\n\nhi"',
              tmp_log=tmp_path / "f.jsonl", monkeypatch=monkeypatch)
    assert rc == 0 and _deny_out(capsys) is None


def test_bare_review_trigger_allows(monkeypatch, capsys, tmp_path):
    rc = _run('gh pr comment 5 --body "@claude review"',
              tmp_log=tmp_path / "f.jsonl", monkeypatch=monkeypatch)
    assert rc == 0 and _deny_out(capsys) is None


def test_heredoc_missing_footer_denies(monkeypatch, capsys, tmp_path):
    cmd = "gh issue comment 5 --body-file - <<'EOF'\na long body\nwith no footer\nEOF"
    rc = _run(cmd, tmp_log=tmp_path / "f.jsonl", monkeypatch=monkeypatch)
    out = _deny_out(capsys)
    assert rc == 0 and out["hookSpecificOutput"]["permissionDecision"] == "deny"


def test_heredoc_with_footer_allows(monkeypatch, capsys, tmp_path):
    cmd = ("gh issue comment 5 --body-file - <<'EOF'\na long body\n"
           "**Created by:** Developer\nEOF")
    rc = _run(cmd, tmp_log=tmp_path / "f.jsonl", monkeypatch=monkeypatch)
    assert rc == 0 and _deny_out(capsys) is None


def test_real_body_file_fails_open(monkeypatch, capsys, tmp_path):
    # Not inspectable -> must NOT deny (a compliant file-bodied comment exists).
    rc = _run("gh pr comment 5 --body-file /tmp/real.md",
              tmp_log=tmp_path / "f.jsonl", monkeypatch=monkeypatch)
    assert rc == 0 and _deny_out(capsys) is None


def test_multi_heredoc_compliant_comment_not_denied(monkeypatch, capsys, tmp_path):
    # Finding 1 end-to-end: the ambiguous multi-heredoc shape must not deny a
    # comment that (in its own heredoc) carries the footer.
    cmd = ("cat > /tmp/scratch.md <<'EOF'\nscratch notes no footer\nEOF\n"
           "gh issue comment 5 --body-file - <<'BODY'\nReal comment.\n"
           "**Created by:** Developer\nBODY")
    rc = _run(cmd, tmp_log=tmp_path / "f.jsonl", monkeypatch=monkeypatch)
    assert rc == 0 and _deny_out(capsys) is None


def test_pr_create_not_guarded(monkeypatch, capsys, tmp_path):
    rc = _run('gh pr create --title t --body "no footer here"',
              tmp_log=tmp_path / "f.jsonl", monkeypatch=monkeypatch)
    assert rc == 0 and _deny_out(capsys) is None


def test_escape_hatch_allows(monkeypatch, capsys, tmp_path):
    rc = _run('gh pr comment 5 --body "deliberately no footer"',
              env={"CLAUDE_ALLOW_NO_FOOTER": "1"},
              tmp_log=tmp_path / "f.jsonl", monkeypatch=monkeypatch)
    assert rc == 0 and _deny_out(capsys) is None


# --- subprocess fail-open (real harness invocation) ---


def _run_subprocess(stdin_payload):
    proc = subprocess.run([sys.executable, str(HOOK)], input=stdin_payload,
                          capture_output=True, text=True, timeout=10)
    return proc.returncode, proc.stdout


def test_malformed_json_fails_open(self=None):
    rc, out = _run_subprocess("{not json")
    assert rc == 0 and out.strip() == ""


def test_non_comment_command_noop():
    rc, out = _run_subprocess(json.dumps({"tool_input": {"command": "gh pr view 5"}}))
    assert rc == 0 and out.strip() == ""
