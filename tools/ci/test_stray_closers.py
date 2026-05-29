"""Tests for `tools/ci/stray_closers.py` (Issue #559).

The detection logic (`find_stray_closers`) is pure and exhaustively unit-tested
here. The CLI / `gh`-fetching path is a thin integration shim verified by a live
smoke test against a real PR (documented in the PR Test plan), not in this suite.
"""

import stray_closers as sc


def nums(result):
    """Just the flagged issue numbers, for terse assertions."""
    return sorted(n for _kw, n, _line in result)


class TestIncidentCase:
    def test_pr543_to_538_commit_prose(self):
        # The originating incident: a commit-body sentence describing the risk.
        # "auto-close" — the hyphen is a word boundary, so `close #538` matches.
        r = sc.find_stray_closers("would auto-close #538 on merge", set())
        assert nums(r) == [538]
        assert r[0][0].lower() == "close"

    def test_offending_line_is_returned(self):
        text = "line one\nwould auto-close #538 on merge\nline three"
        r = sc.find_stray_closers(text, set())
        assert "auto-close #538" in r[0][2]


class TestIntendedClosersSkipped:
    def test_intended_closer_in_set_not_flagged(self):
        assert sc.find_stray_closers("Closes #549", {549}) == []

    def test_set_accepts_strings(self):
        assert sc.find_stray_closers("Closes #549", {"549"}) == []

    def test_mix_flags_only_stray(self):
        text = "Closes #549\nalso fixes #538 as a side effect"
        assert nums(sc.find_stray_closers(text, {549})) == [538]


class TestNegationIgnored:
    def test_does_not_close_still_flagged(self):
        # GitHub's parser is regex-only; English negation is ignored.
        assert nums(sc.find_stray_closers("this does not close #334", set())) == [334]


class TestWordBoundaries:
    def test_disclose_not_matched(self):
        assert sc.find_stray_closers("disclose #5 in the report", set()) == []

    def test_prefix_not_matched(self):
        # "prefix" ends in "fix" but has no preceding word boundary
        assert sc.find_stray_closers("prefix #5 details", set()) == []

    def test_closing_not_matched(self):
        assert sc.find_stray_closers("closing #5 of the loop", set()) == []

    def test_fixing_not_matched(self):
        assert sc.find_stray_closers("fixing #5 up", set()) == []

    def test_closer_not_matched(self):
        assert sc.find_stray_closers("the closer #5 finished", set()) == []

    def test_hyphen_prefix_is_a_boundary(self):
        assert nums(sc.find_stray_closers("auto-fixes #9", set())) == [9]

    def test_prose_between_keyword_and_ref_not_matched(self):
        # "word between keyword and #N" invariant, e.g. "closed epic Issue #538"
        # — made a first-class assertion per PR #562 review (was only exercised
        # implicitly through main()).
        assert sc.find_stray_closers("closed epic Issue #538", set()) == []

    def test_newline_does_not_bridge_keyword_and_ref(self):
        # assemble_squash_text joins fields with "\n"; a field ending in a bare
        # keyword must NOT cross-match a #N opening the next field (PR #562 review).
        assert sc.find_stray_closers("title ends with fixes\n#538 opens body", set()) == []


class TestKeywordVariants:
    def test_all_variants(self):
        for kw in ("close", "closes", "closed", "fix", "fixes", "fixed",
                   "resolve", "resolves", "resolved"):
            r = sc.find_stray_closers(f"{kw} #7", set())
            assert nums(r) == [7], f"{kw!r} should flag #7"

    def test_case_insensitive(self):
        assert nums(sc.find_stray_closers("FIXES #7", set())) == [7]

    def test_colon_separator(self):
        assert nums(sc.find_stray_closers("Closes: #8", set())) == [8]


class TestMisc:
    def test_dedup_by_number(self):
        # same number referenced twice → reported once
        r = sc.find_stray_closers("fixes #5 and closes #5", set())
        assert nums(r) == [5] and len(r) == 1

    def test_two_distinct_strays(self):
        assert nums(sc.find_stray_closers("fixes #5; resolves #6", set())) == [5, 6]

    def test_keyword_without_hash_not_matched(self):
        assert sc.find_stray_closers("this closes the door", set()) == []

    def test_empty_text(self):
        assert sc.find_stray_closers("", set()) == []

    def test_none_text(self):
        assert sc.find_stray_closers(None, set()) == []


class TestAssembleAndClosingSet:
    def test_assemble_joins_all_fields(self):
        data = {
            "title": "feat: thing",
            "body": "Closes #1",
            "commits": [
                {"messageHeadline": "do x", "messageBody": "auto-close #2"},
                {"messageHeadline": "do y", "messageBody": None},
            ],
        }
        text = sc.assemble_squash_text(data)
        assert "feat: thing" in text and "Closes #1" in text
        assert "auto-close #2" in text and "do y" in text

    def test_assemble_tolerates_missing_fields(self):
        assert sc.assemble_squash_text({}) == "\n"  # title + body both empty
        assert sc.assemble_squash_text({"commits": None}) == "\n"

    def test_closing_set_from(self):
        data = {"closingIssuesReferences": [{"number": 5}, {"number": 9}]}
        assert sc.closing_set_from(data) == {5, 9}

    def test_closing_set_empty(self):
        assert sc.closing_set_from({}) == set()


class TestMainOrchestration:
    """main() with fetch_pr monkeypatched — no network."""

    def _patch(self, monkeypatch, data):
        monkeypatch.setattr(sc, "fetch_pr", lambda pr, repo: data)

    def test_clean_pr_returns_0(self, monkeypatch, capsys):
        self._patch(monkeypatch, {
            "title": "feat: x", "body": "Closes #549",
            "commits": [{"messageHeadline": "x", "messageBody": "closed epic Issue #538"}],
            "closingIssuesReferences": [{"number": 549}],
        })
        # "closed epic Issue #538" is NOT an adjacency (word between) → clean
        assert sc.main(["stray_closers.py", "560"]) == 0
        assert capsys.readouterr().err == ""

    def test_stray_pr_returns_1(self, monkeypatch, capsys):
        self._patch(monkeypatch, {
            "title": "feat: x", "body": "Closes #549",
            "commits": [{"messageHeadline": "x", "messageBody": "would auto-close #538"}],
            "closingIssuesReferences": [{"number": 549}],
        })
        assert sc.main(["stray_closers.py", "560"]) == 1
        assert "#538" in capsys.readouterr().err

    def test_gh_error_fails_open(self, monkeypatch, capsys):
        def boom(pr, repo):
            raise FileNotFoundError("gh not found")
        monkeypatch.setattr(sc, "fetch_pr", boom)
        assert sc.main(["stray_closers.py", "560"]) == 0  # fail open, don't block

    def test_bad_args_returns_2(self, capsys):
        assert sc.main(["stray_closers.py"]) == 2
        assert sc.main(["stray_closers.py", "notanumber"]) == 2
