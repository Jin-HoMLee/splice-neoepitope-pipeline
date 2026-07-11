"""Tests for `.agents/hooks/post_gh_pr_create.py` (Issue #550).

Pure-function unit tests import the module directly (no network). The fail-open /
no-op subprocess tests mirror the harness invocation (pipe PostToolUse JSON to
stdin) and never reach a live `gh` call because no actionable PR URL is parsed.

The full add+flip path is inherently a live-`gh` integration and is verified by
dogfooding the hook against this Issue's own PR, not in this suite.
"""

import io
import json
import subprocess
import sys
from pathlib import Path

HOOKS_DIR = Path(__file__).parent.parent.parent / ".agents" / "hooks"
HOOK = HOOKS_DIR / "post_gh_pr_create.py"
sys.path.insert(0, str(HOOKS_DIR))

import post_gh_pr_create as h  # noqa: E402


class TestMatchesPrCreate:
    def test_plain_invocation(self):
        assert h.matches_pr_create("gh pr create --title x --body y") is True

    def test_after_separator(self):
        assert h.matches_pr_create("git push -u origin br && gh pr create --fill") is True

    def test_var_prefix_same_line_matched(self):
        # `VAR=value gh pr create` — the harness Bash matcher is subcommand-aware
        # and strips the assignment prefix; the hook must mirror that (Issue #561).
        assert h.matches_pr_create('B="x" gh pr create --fill') is True

    def test_var_prefix_newline_matched(self):
        # `VAR=value` then `gh pr create` on the next line — PR #560's actual form,
        # which silently failed to auto-board (Issue #561).
        assert h.matches_pr_create('B="x"\ngh pr create --fill') is True

    def test_multiple_var_prefixes_matched(self):
        # Several leading assignments are all skipped before the command start.
        assert h.matches_pr_create("A=1 B=2 gh pr create --fill") is True

    def test_var_assignment_only_not_matched(self):
        # A bare assignment with no following `gh pr create` must not match.
        assert h.matches_pr_create('B="x"') is False

    def test_var_prefix_before_other_subcommand_not_matched(self):
        # The assignment prefix is stripped, but the command is `gh pr view`, not
        # `gh pr create` — must not match.
        assert h.matches_pr_create('B="x" gh pr view 1') is False

    def test_pr_view_not_matched(self):
        assert h.matches_pr_create("gh pr view 1") is False

    def test_pr_list_not_matched(self):
        assert h.matches_pr_create("gh pr list") is False

    def test_echoed_string_not_matched(self):
        # "echo gh pr create" is not an invocation of the command
        assert h.matches_pr_create("echo gh pr create") is False

    def test_create_inside_quoted_comment_body_not_matched(self):
        # Real false positive (PR #558): a `gh pr comment` whose body discusses
        # `gh pr create` must NOT match — the substring is inside a quoted arg.
        cmd = 'gh pr comment 1 --body "narrow if to Bash(gh pr create *)"'
        assert h.matches_pr_create(cmd) is False

    def test_create_after_literal_separator_in_quoted_body_not_matched(self):
        # The `&&` lives inside the quoted body, so it is not a shell separator.
        cmd = 'gh pr comment 1 --body "example: git push && gh pr create --fill"'
        assert h.matches_pr_create(cmd) is False

    def test_unbalanced_quotes_fail_safe(self):
        # Untokenizable command → fail safe (do not fire).
        assert h.matches_pr_create('gh pr comment 1 --body "oops') is False


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


def test_should_track():
    assert h.should_track("Jin-HoMLee", "splice-neoepitope-pipeline") is True
    assert h.should_track("Jin-HoMLee", "claude-personas-splice-neoepitope-pipeline") is True
    assert h.should_track("someone-else", "splice-neoepitope-pipeline") is False
    assert h.should_track("Jin-HoMLee", "some-unrelated-repo") is False


def test_status_for_draft():
    # draft PR opened mid-In-progress → In progress; non-draft → Ready for review
    assert h.status_for_draft(True) == ("In progress", "47fc9ee4")
    assert h.status_for_draft(False) == ("Ready for review", "8bf9192f")


class TestSkipBotReviewMarker:
    """The trivial-PR opt-out (Issue #1073 AC 2).

    Mirrors the existing `<!-- skip-lab-notebook: routine -->` marker convention
    (`closure_audit._SKIP_LAB_NOTEBOOK`) rather than inventing new vocabulary.
    """

    def test_canonical_marker(self):
        assert h.has_skip_bot_review("body\n<!-- skip-bot-review: trivial -->\n") is True

    def test_case_and_space_tolerant(self):
        assert h.has_skip_bot_review("<!--   SKIP-BOT-REVIEW: trivial   -->") is True

    def test_any_reason_accepted(self):
        # The reason is documentation for the human, not a parsed enum.
        assert h.has_skip_bot_review("<!-- skip-bot-review: docs-only typo -->") is True

    def test_absent(self):
        assert h.has_skip_bot_review("A normal PR body. Closes #1073.") is False

    def test_empty_body(self):
        assert h.has_skip_bot_review("") is False
        assert h.has_skip_bot_review(None) is False

    def test_mention_without_comment_syntax_does_not_match(self):
        # Prose *discussing* the marker must not silently opt the PR out.
        assert h.has_skip_bot_review("We could add skip-bot-review here someday.") is False


class TestShouldRequestReview:
    """AC 1 + AC 2: non-trivial PRs auto-request; drafts and opt-outs do not."""

    def test_non_draft_plain_body_requests(self):
        assert h.should_request_review(is_draft=False, body="Closes #1073.") is True

    def test_draft_does_not_request(self):
        # A draft is not up for review yet - the hook's own existing semantic.
        assert h.should_request_review(is_draft=True, body="Closes #1073.") is False

    def test_skip_marker_does_not_request(self):
        assert h.should_request_review(
            is_draft=False, body="<!-- skip-bot-review: trivial -->"
        ) is False

    def test_draft_and_marker_does_not_request(self):
        assert h.should_request_review(
            is_draft=True, body="<!-- skip-bot-review: trivial -->"
        ) is False


class TestTriggerReconciliation:
    """AC 4: no double-request against the merge-time `bot_review_offer.py` gate.

    The gate is a *detector*: `audit_and_merge.sh` prompts only when the literal
    trigger is absent from the PR's comments. So the reconciliation is not a code
    path - it holds iff the string this hook posts is a string that gate detects.
    Assert exactly that, across the module boundary.

    The falsifier is real: change the hook's trigger to the hyphenated reference
    form `@-claude review` (which does NOT fire the GitHub Action) and this test
    goes red - which is the bug it exists to catch.
    """

    def test_hook_trigger_is_detected_by_the_merge_gate(self):
        sys.path.insert(0, str(Path(__file__).parent))
        import bot_review_offer as gate  # noqa: PLC0415

        posted_comment = {"body": h.REVIEW_TRIGGER}
        assert gate.has_bot_review_offer([posted_comment]) is True

    def test_hook_trigger_is_the_real_action_trigger(self):
        # The hyphenated `@-claude review` reference form must never be posted:
        # it does not fire the Action, so an auto-"requested" PR would merge
        # un-reviewed. Guard the literal.
        assert h.REVIEW_TRIGGER == "@claude review"


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


class TestAutoRequestPath:
    """Matched-pair control over the auto-request decision (Issue #1073).

    Identical inputs, one variable flipped, opposite expected outcomes. If the
    auto-request path were dead, a lone "draft does not request" assertion would
    still pass - the pair is what makes this a check rather than a ritual.

    Also pins the Status consequence: when the hook auto-requests, the card must
    NOT be parked at `Ready for review`. The sibling flip-to-`In review` hook only
    sees Claude's Bash calls, so it cannot see this hook's subprocess comment; if
    this hook set `Ready for review` and posted the trigger, the card would sit in
    the wrong column for the whole review - the exact stranding Issue #996 fixed.
    """

    def _drive(self, monkeypatch, is_draft, body):
        calls = {"comment": [], "flip": [], "status": []}
        monkeypatch.setattr(h, "_pr_view", lambda url: (is_draft, body))
        monkeypatch.setattr(h, "_add_to_board", lambda url: "ITEM_1")
        monkeypatch.setattr(h, "_set_status", lambda item, opt: calls["status"].append(opt))
        monkeypatch.setattr(h, "_request_bot_review", lambda url: calls["comment"].append(url))
        monkeypatch.setattr(h, "_apply_review_request", lambda url: calls["flip"].append(url))
        monkeypatch.setattr(h, "_log_fire", lambda *a, **k: None)
        monkeypatch.setattr(
            sys, "stdin",
            io.StringIO(_payload(
                "gh pr create --fill",
                "https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/1124",
            )),
        )
        assert h.main() == 0
        return calls

    def test_non_draft_requests_review_and_delegates_the_in_review_flip(self, monkeypatch):
        calls = self._drive(monkeypatch, is_draft=False, body="Closes #1073.")
        assert calls["comment"] == [
            "https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/1124"
        ]
        assert len(calls["flip"]) == 1
        # Status delegated to the In-review flip; NOT set to Ready for review here.
        assert calls["status"] == []

    def test_draft_does_not_request_and_stays_in_progress(self, monkeypatch):
        calls = self._drive(monkeypatch, is_draft=True, body="Closes #1073.")
        assert calls["comment"] == []
        assert calls["flip"] == []
        assert calls["status"] == [h.IN_PROGRESS_OPTION]

    def test_skip_marker_does_not_request_and_stays_ready_for_review(self, monkeypatch):
        calls = self._drive(
            monkeypatch, is_draft=False, body="<!-- skip-bot-review: trivial -->"
        )
        assert calls["comment"] == []
        assert calls["flip"] == []
        assert calls["status"] == [h.READY_FOR_REVIEW_OPTION]

    def test_review_request_failure_is_fail_open(self, monkeypatch):
        # A gh hiccup posting the trigger must never break the user's flow: the
        # board add already landed, and the merge-time offer gate still catches
        # the missing review. Degrades to today's behavior, never worse.
        def boom(url):
            raise subprocess.CalledProcessError(1, ["gh"])

        monkeypatch.setattr(h, "_pr_view", lambda url: (False, "Closes #1073."))
        monkeypatch.setattr(h, "_add_to_board", lambda url: "ITEM_1")
        monkeypatch.setattr(h, "_request_bot_review", boom)
        monkeypatch.setattr(
            sys, "stdin",
            io.StringIO(_payload(
                "gh pr create --fill",
                "https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/1124",
            )),
        )
        assert h.main() == 0


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
