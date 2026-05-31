"""Tests for `tools/ci/bot_review_offer.py` (Issue #443).

The detection logic (`has_bot_review_offer`) is pure and exhaustively unit-tested
here. The interactive (a/b/c) prompt + the `gh pr comment` post live in
`scripts/audit_and_merge.sh` (bash) and are verified by a live smoke test against
a real PR (documented in the PR Test plan), not in this suite — mirroring the
`stray_closers.py` split.
"""

import bot_review_offer as bro


def comments(*bodies):
    """Build the `gh pr view --json comments` shape from raw body strings."""
    return [{"body": b} for b in bodies]


class TestRealTriggerDetected:
    def test_canonical_trigger(self):
        assert bro.has_bot_review_offer(comments("@claude review")) is True

    def test_trigger_as_substring_in_sentence(self):
        assert bro.has_bot_review_offer(comments("Running @claude review now, thanks!")) is True

    def test_case_insensitive(self):
        assert bro.has_bot_review_offer(comments("@Claude Review")) is True

    def test_extra_internal_spaces(self):
        assert bro.has_bot_review_offer(comments("@claude   review")) is True

    def test_tab_separator(self):
        # `[ \t]+` covers a tab between mention and "review" — lock the invariant.
        assert bro.has_bot_review_offer(comments("@claude\treview")) is True

    def test_one_of_many_comments(self):
        cs = comments("first comment", "looks good", "@claude review")
        assert bro.has_bot_review_offer(cs) is True


class TestNonTriggers:
    def test_hyphen_reference_form_does_not_count(self):
        # `@-claude review` is the *non-triggering* reference form (the literal
        # substring @claude never appears), so no review was actually requested.
        assert bro.has_bot_review_offer(comments("offer @-claude review after merge")) is False

    def test_mention_without_review_word(self):
        # A bare @claude mention is not a *review* offer.
        assert bro.has_bot_review_offer(comments("@claude can you help with this?")) is False

    def test_newline_does_not_bridge_mention_and_review(self):
        # The trigger must be contiguous; a mention ending a line must NOT
        # cross-match a "review" opening the next line.
        assert bro.has_bot_review_offer(comments("ping @claude\nreview the diff")) is False

    def test_no_trigger_anywhere(self):
        assert bro.has_bot_review_offer(comments("LGTM", "merging now")) is False


class TestTolerantInputs:
    def test_empty_list(self):
        assert bro.has_bot_review_offer([]) is False

    def test_none_list(self):
        assert bro.has_bot_review_offer(None) is False

    def test_none_body_tolerated(self):
        assert bro.has_bot_review_offer([{"body": None}]) is False

    def test_missing_body_key_tolerated(self):
        assert bro.has_bot_review_offer([{"author": {"login": "x"}}]) is False


class TestMainOrchestration:
    """main() with fetch_comments monkeypatched — no network."""

    def _patch(self, monkeypatch, data):
        monkeypatch.setattr(bro, "fetch_comments", lambda pr, repo: data)

    def test_offered_prints_token(self, monkeypatch, capsys):
        self._patch(monkeypatch, comments("@claude review"))
        assert bro.main(["bot_review_offer.py", "443"]) == 0
        assert capsys.readouterr().out.strip() == "OFFERED"

    def test_not_offered_prints_token(self, monkeypatch, capsys):
        self._patch(monkeypatch, comments("LGTM"))
        assert bro.main(["bot_review_offer.py", "443"]) == 0
        assert capsys.readouterr().out.strip() == "NOT_OFFERED"

    def test_gh_error_fails_open_to_offered(self, monkeypatch, capsys):
        def boom(pr, repo):
            raise FileNotFoundError("gh not found")
        monkeypatch.setattr(bro, "fetch_comments", boom)
        # Fail open: a gh hiccup must not disrupt the merge with a spurious prompt.
        assert bro.main(["bot_review_offer.py", "443"]) == 0
        captured = capsys.readouterr()
        assert captured.out.strip() == "OFFERED"
        assert "skipped" in captured.err.lower()

    def test_bad_args_returns_2(self):
        assert bro.main(["bot_review_offer.py"]) == 2
        assert bro.main(["bot_review_offer.py", "notanumber"]) == 2
