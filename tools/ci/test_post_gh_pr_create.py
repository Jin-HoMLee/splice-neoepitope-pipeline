"""Tests for `.claude/hooks/post_gh_pr_create.py` (Issue #550).

Pure-function unit tests import the module directly (no network). The fail-open /
no-op subprocess tests mirror the harness invocation (pipe PostToolUse JSON to
stdin) and never reach a live `gh` call because no actionable PR URL is parsed.

The full add+flip path is inherently a live-`gh` integration and is verified by
dogfooding the hook against this Issue's own PR, not in this suite.
"""

import json
import subprocess
import sys
from pathlib import Path

HOOKS_DIR = Path(__file__).parent.parent.parent / ".claude" / "hooks"
HOOK = HOOKS_DIR / "post_gh_pr_create.py"
sys.path.insert(0, str(HOOKS_DIR))

import post_gh_pr_create as h  # noqa: E402


class TestMatchesPrCreate:
    def test_plain_invocation(self):
        assert h.matches_pr_create("gh pr create --title x --body y") is True

    def test_after_separator(self):
        assert h.matches_pr_create("git push -u origin br && gh pr create --fill") is True

    def test_pr_view_not_matched(self):
        assert h.matches_pr_create("gh pr view 1") is False

    def test_pr_list_not_matched(self):
        assert h.matches_pr_create("gh pr list") is False

    def test_echoed_string_not_matched(self):
        # "echo gh pr create" is not an invocation of the command
        assert h.matches_pr_create("echo gh pr create") is False


class TestParsePrUrl:
    def test_basic_url(self):
        assert h.parse_pr_url(
            "https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/556"
        ) == ("Jin-HoMLee", "splice-neoepitope-pipeline", 556)

    def test_url_embedded_in_output(self):
        text = (
            "Creating pull request for branch into main\n"
            "https://github.com/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/pull/7\n"
        )
        assert h.parse_pr_url(text)[2] == 7

    def test_no_url_returns_none(self):
        assert h.parse_pr_url("error: pull request create failed") is None

    def test_takes_last_url(self):
        text = (
            "see https://github.com/Jin-HoMLee/x/pull/1 then "
            "created https://github.com/Jin-HoMLee/x/pull/2"
        )
        assert h.parse_pr_url(text)[2] == 2


def test_has_project_flag():
    assert h.has_project_flag('gh pr create --project "JH M Lee Lab"') is True
    assert h.has_project_flag("gh pr create --title x --body y") is False


def test_should_track():
    assert h.should_track("Jin-HoMLee", "splice-neoepitope-pipeline") is True
    assert h.should_track("Jin-HoMLee", "claude-personas-splice-neoepitope-pipeline") is True
    assert h.should_track("someone-else", "splice-neoepitope-pipeline") is False
    assert h.should_track("Jin-HoMLee", "some-unrelated-repo") is False


def test_status_for_draft():
    # draft PR opened mid-In-progress → In progress; non-draft → Ready for review
    assert h.status_for_draft(True) == ("In progress", "47fc9ee4")
    assert h.status_for_draft(False) == ("Ready for review", "8bf9192f")


def test_extract_output_handles_shapes():
    # bare string
    assert "pull/9" in h.extract_output("done https://github.com/a/b/pull/9")
    # dict with stdout
    assert "pull/9" in h.extract_output({"stdout": "https://github.com/a/b/pull/9\n"})
    # unknown dict shape falls back to stringifying the whole object
    assert "pull/9" in h.extract_output({"weird_key": "https://github.com/a/b/pull/9"})
    # None is safe
    assert h.extract_output(None) == ""


# --- subprocess fail-open / no-op paths (no live gh reached) ---


def _run(stdin_payload: str):
    proc = subprocess.run(
        [sys.executable, str(HOOK)],
        input=stdin_payload,
        capture_output=True,
        text=True,
        timeout=10,
    )
    return proc.returncode, proc.stdout, proc.stderr


def _payload(command: str, output: str = "") -> str:
    return json.dumps({"tool_input": {"command": command}, "tool_response": output})


class TestSubprocessNoOp:
    def test_empty_stdin_fail_open(self):
        rc, out, _ = _run("")
        assert rc == 0 and out.strip() == ""

    def test_malformed_json_fail_open(self):
        rc, out, _ = _run("{not valid json")
        assert rc == 0 and out.strip() == ""

    def test_non_pr_create_command_noop(self):
        rc, out, _ = _run(_payload("gh pr view 1", "anything"))
        assert rc == 0 and out.strip() == ""

    def test_pr_create_without_url_noop(self):
        # matched command but no PR URL in output → no gh call, no output
        rc, out, _ = _run(_payload("gh pr create --fill", "error: validation failed"))
        assert rc == 0 and out.strip() == ""

    def test_untracked_repo_noop(self):
        # a parsed PR URL for an untracked owner/repo must not be boarded
        rc, out, _ = _run(
            _payload("gh pr create --fill", "https://github.com/someone-else/other/pull/3")
        )
        assert rc == 0 and out.strip() == ""
