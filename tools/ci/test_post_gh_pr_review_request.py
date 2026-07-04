"""Tests for `.agents/hooks/post_gh_pr_review_request.py` (Issue #996).

Pure-function unit tests import the module directly (no network). The subprocess
no-op / fail-open tests pipe PostToolUse JSON to stdin and never reach a live
`gh` call because the command either isn't a review request or resolves to
nothing. The full flip path is a live-`gh` integration verified by dogfooding
the hook against this Issue's own PR, not in this suite.
"""

import json
import subprocess
import sys
from pathlib import Path

HOOKS_DIR = Path(__file__).parent.parent.parent / ".agents" / "hooks"
HOOK = HOOKS_DIR / "post_gh_pr_review_request.py"
sys.path.insert(0, str(HOOKS_DIR))

import post_gh_pr_review_request as h  # noqa: E402


class TestMatchesReviewRequest:
    def test_comment_with_trigger(self):
        assert h.matches_review_request('gh pr comment 996 --body "@claude review"') == "996"

    def test_comment_trigger_case_insensitive(self):
        assert h.matches_review_request('gh pr comment 996 --body "@Claude Review"') == "996"

    def test_comment_trigger_embedded_in_longer_body(self):
        cmd = 'gh pr comment 996 --body "thanks, could you take a look? @claude review please"'
        assert h.matches_review_request(cmd) == "996"

    def test_comment_without_trigger_not_matched(self):
        assert h.matches_review_request('gh pr comment 996 --body "looks good to me"') is None

    def test_hyphenated_reference_form_not_matched(self):
        # `@-claude review` is the guard-dodging reference form, NOT a bot trigger.
        assert h.matches_review_request('gh pr comment 996 --body "posted @-claude review"') is None

    def test_review_subcommand_matched(self):
        assert h.matches_review_request("gh pr review 996 --approve") == "996"

    def test_review_bare_matched(self):
        assert h.matches_review_request("gh pr review 996") == "996"

    def test_url_ref_matched(self):
        cmd = 'gh pr review https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/996'
        assert h.matches_review_request(cmd).endswith("/pull/996")

    def test_after_separator(self):
        assert h.matches_review_request('git push && gh pr review 996 --approve') == "996"

    def test_compound_trigger_in_later_echo_not_matched(self):
        # The trigger lives in a separate `&&`-joined echo, NOT the posted comment
        # body - must not false-positive (segment-bounded scan).
        cmd = 'gh pr comment 996 --body "plain" && echo "@claude review"'
        assert h.matches_review_request(cmd) is None

    def test_compound_triggerless_comment_then_review_matches(self):
        # A triggerless comment must not abort the scan - the later `gh pr review`
        # is the real review request.
        cmd = 'gh pr comment 996 --body "plain" && gh pr review 42'
        assert h.matches_review_request(cmd) == "42"

    def test_var_prefix_matched(self):
        assert h.matches_review_request('GH_TOKEN=x gh pr review 996') == "996"

    def test_flag_before_ref_fails_safe(self):
        # ref not in the first positional slot -> can't parse reliably -> None
        assert h.matches_review_request('gh pr comment --body "@claude review" 996') is None

    def test_other_subcommand_not_matched(self):
        assert h.matches_review_request("gh pr view 996") is None
        assert h.matches_review_request("gh pr merge 996") is None

    def test_comment_carrying_literal_trigger_matches(self):
        # A comment body containing the literal trigger legitimately matches:
        # posting that body WOULD fire the bot, so flipping to In review is
        # correct. The guard keys on command shape + literal trigger, not intent.
        cmd = 'gh pr comment 5 --body "the hook fires on \'@claude review\'"'
        assert h.matches_review_request(cmd) == "5"

    def test_unbalanced_quotes_fail_safe(self):
        assert h.matches_review_request('gh pr comment 996 --body "@claude review') is None

    def test_review_discussed_in_echo_not_matched(self):
        assert h.matches_review_request('echo gh pr review 996') is None


class TestHasTrigger:
    def test_exact(self):
        assert h._has_trigger("@claude review") is True

    def test_case_insensitive(self):
        assert h._has_trigger("@CLAUDE REVIEW") is True

    def test_hyphenated_not_trigger(self):
        assert h._has_trigger("@-claude review") is False

    def test_none_safe(self):
        assert h._has_trigger(None) is False


class TestParsePrUrl:
    def test_basic(self):
        assert h.parse_pr_url(
            "https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/996"
        ) == ("Jin-HoMLee", "splice-neoepitope-pipeline", 996)

    def test_none(self):
        assert h.parse_pr_url("no url here") is None


def test_should_track():
    assert h.should_track("Jin-HoMLee", "splice-neoepitope-pipeline") is True
    assert h.should_track("Jin-HoMLee", "claude-personas-splice-neoepitope-pipeline") is True
    assert h.should_track("someone-else", "splice-neoepitope-pipeline") is False
    assert h.should_track("Jin-HoMLee", "unrelated-repo") is False


class TestShouldFlip:
    def test_ready_for_review_flips(self):
        assert h.should_flip("Ready for review") is True

    def test_ready_flips(self):
        assert h.should_flip("Ready") is True

    def test_in_progress_flips(self):
        # a draft PR under review is a forward move, not a regression
        assert h.should_flip("In progress") is True

    def test_unset_flips(self):
        assert h.should_flip(None) is True

    def test_already_in_review_noop(self):
        assert h.should_flip("In review") is False

    def test_epic_parked_not_touched(self):
        assert h.should_flip("Epic") is False

    def test_done_terminal_not_touched(self):
        assert h.should_flip("Done") is False


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


def _payload(command: str) -> str:
    return json.dumps({"tool_input": {"command": command}})


class TestSubprocessNoOp:
    def test_empty_stdin_fail_open(self):
        rc, out, _ = _run("")
        assert rc == 0 and out.strip() == ""

    def test_malformed_json_fail_open(self):
        rc, out, _ = _run("{not valid json")
        assert rc == 0 and out.strip() == ""

    def test_non_review_command_noop(self):
        rc, out, _ = _run(_payload("gh pr view 996"))
        assert rc == 0 and out.strip() == ""

    def test_comment_without_trigger_noop(self):
        rc, out, _ = _run(_payload('gh pr comment 996 --body "nice work"'))
        assert rc == 0 and out.strip() == ""


# --- orchestration with a mocked gh (in-process, no network) ---


class TestMainFlipPath:
    def test_trigger_flips_linked_issue(self, monkeypatch, capsys):
        import io
        flips = []
        lookups = []
        monkeypatch.setattr(sys, "stdin",
                            io.StringIO(_payload('gh pr comment 996 --body "@claude review"')))
        monkeypatch.setattr(
            h, "_pr_linked_issues",
            lambda ref: ("Jin-HoMLee", "claude-personas-splice-neoepitope-pipeline", [996]))

        def fake_lookup(issue, owner, repo):
            lookups.append((issue, owner, repo))
            return ("ITEM_996", "Ready for review")

        monkeypatch.setattr(h, "_issue_item_and_status", fake_lookup)
        monkeypatch.setattr(h, "_set_status", lambda item_id, opt: flips.append((item_id, opt)))
        monkeypatch.setattr(h, "_log_fire", lambda *a: None)
        assert h.main() == 0
        assert flips == [("ITEM_996", h.IN_REVIEW_OPTION)]
        # the PR's real owner/repo is threaded into the issue lookup (not hardcoded)
        assert lookups == [(996, "Jin-HoMLee", "claude-personas-splice-neoepitope-pipeline")]
        assert "In review" in capsys.readouterr().out

    def test_already_in_review_is_noop(self, monkeypatch, capsys):
        import io
        flips = []
        monkeypatch.setattr(sys, "stdin",
                            io.StringIO(_payload("gh pr review 996 --approve")))
        monkeypatch.setattr(
            h, "_pr_linked_issues",
            lambda ref: ("Jin-HoMLee", "splice-neoepitope-pipeline", [996]))
        monkeypatch.setattr(
            h, "_issue_item_and_status", lambda issue, owner, repo: ("ITEM_996", "In review"))
        monkeypatch.setattr(h, "_set_status", lambda item_id, opt: flips.append(opt))
        assert h.main() == 0
        assert flips == []
        assert capsys.readouterr().out.strip() == ""

    def test_untracked_repo_noop(self, monkeypatch, capsys):
        import io
        flips = []
        monkeypatch.setattr(sys, "stdin",
                            io.StringIO(_payload("gh pr review 996")))
        monkeypatch.setattr(
            h, "_pr_linked_issues", lambda ref: ("someone-else", "other-repo", [1]))
        monkeypatch.setattr(h, "_set_status", lambda item_id, opt: flips.append(opt))
        assert h.main() == 0
        assert flips == []

    def test_fail_open_on_gh_error(self, monkeypatch, capsys):
        import io
        monkeypatch.setattr(sys, "stdin",
                            io.StringIO(_payload("gh pr review 996 --approve")))

        def boom(*a, **k):
            raise subprocess.CalledProcessError(1, ["gh"])

        monkeypatch.setattr(h, "_gh", boom)
        assert h.main() == 0
        assert capsys.readouterr().out.strip() == ""
