"""Tests for `scripts/pm/check_closed_recent.py` (Issue #784).

Pure-helper unit tests (window floor, the `>=` boundary, inclusion predicate) plus
a network-free orchestration test with injected `gh` result lists. No live `gh`.
"""

import sys
from datetime import datetime, timezone
from pathlib import Path

SCRIPTS_PM = Path(__file__).parent.parent / "pm"
sys.path.insert(0, str(SCRIPTS_PM))

import check_closed_recent as c  # noqa: E402

UTC = timezone.utc


class TestComputeFloor:
    def test_default_days_1_is_start_of_yesterday(self):
        now = datetime(2026, 7, 4, 20, 30, tzinfo=UTC)
        assert c.compute_floor(now, 1) == datetime(2026, 7, 3, 0, 0, tzinfo=UTC)

    def test_days_0_is_start_of_today(self):
        now = datetime(2026, 7, 4, 20, 30, tzinfo=UTC)
        assert c.compute_floor(now, 0) == datetime(2026, 7, 4, 0, 0, tzinfo=UTC)

    def test_crosses_month_boundary(self):
        now = datetime(2026, 7, 1, 3, 0, tzinfo=UTC)
        assert c.compute_floor(now, 1) == datetime(2026, 6, 30, 0, 0, tzinfo=UTC)

    def test_negative_days_raises(self):
        try:
            c.compute_floor(datetime(2026, 7, 4, tzinfo=UTC), -1)
            assert False, "expected ValueError"
        except ValueError:
            pass


class TestBuildSearch:
    def test_uses_inclusive_boundary(self):
        floor = datetime(2026, 6, 18, 0, 0, tzinfo=UTC)
        assert c.build_search(floor) == "closed:>=2026-06-18"

    def test_is_not_the_exclusive_operator(self):
        # The whole point of #784: the bare `>` dropped the boundary day.
        floor = datetime(2026, 6, 18, 0, 0, tzinfo=UTC)
        assert ">=" in c.build_search(floor)
        assert ">2026" not in c.build_search(floor).replace(">=", "")


class TestIsWithin:
    def test_item_at_boundary_day_start_is_included(self):
        # AC2: an item closed at 00:00:00 of the boundary day must be INCLUDED.
        floor = datetime(2026, 6, 18, 0, 0, tzinfo=UTC)
        assert c.is_within("2026-06-18T00:00:00Z", floor) is True

    def test_item_during_boundary_day_included(self):
        floor = datetime(2026, 6, 18, 0, 0, tzinfo=UTC)
        assert c.is_within("2026-06-18T09:30:00Z", floor) is True

    def test_item_just_before_floor_excluded(self):
        floor = datetime(2026, 6, 18, 0, 0, tzinfo=UTC)
        assert c.is_within("2026-06-17T23:59:59Z", floor) is False

    def test_unparseable_timestamp_kept(self):
        floor = datetime(2026, 6, 18, 0, 0, tzinfo=UTC)
        assert c.is_within(None, floor) is True
        assert c.is_within("", floor) is True


class TestFilterMerged:
    def test_keeps_merged_drops_unmerged(self):
        rows = [
            {"number": 1, "mergedAt": "2026-07-04T09:00:00Z"},
            {"number": 2, "mergedAt": None},   # closed-unmerged -> dropped
            {"number": 3},                       # no mergedAt key -> dropped
        ]
        assert [r["number"] for r in c.filter_merged(rows)] == [1]


class TestCollect:
    FLOOR = datetime(2026, 7, 3, 0, 0, tzinfo=UTC)

    def test_merges_and_sorts_newest_first(self):
        issues = [{"number": 10, "title": "issue ten", "closedAt": "2026-07-03T08:00:00Z"}]
        prs = [{"number": 20, "title": "pr twenty", "mergedAt": "2026-07-04T09:00:00Z"}]
        rows = c.collect(self.FLOOR, issues, prs)
        assert [r["number"] for r in rows] == [20, 10]  # PR is newer -> first
        assert [r["kind"] for r in rows] == ["PR", "Issue"]

    def test_reapplies_boundary_clientside(self):
        # An item before the floor slips through even if the search returned it.
        issues = [{"number": 1, "title": "stale", "closedAt": "2026-07-02T23:00:00Z"}]
        assert c.collect(self.FLOOR, issues, []) == []

    def test_boundary_start_kept(self):
        issues = [{"number": 5, "title": "edge", "closedAt": "2026-07-03T00:00:00Z"}]
        rows = c.collect(self.FLOOR, issues, [])
        assert len(rows) == 1 and rows[0]["number"] == 5

    def test_fail_open_placeholder_sorts_last_not_first(self):
        # A row whose timestamp is unparseable (kept fail-open, when="?") must not
        # jump to the top of a newest-first list on raw-string sort.
        issues = [
            {"number": 1, "title": "real", "closedAt": "2026-07-04T09:00:00Z"},
            {"number": 2, "title": "malformed", "closedAt": "not-a-date"},
        ]
        rows = c.collect(self.FLOOR, issues, [])
        assert [r["number"] for r in rows] == [1, 2]  # real first, placeholder last


class TestRender:
    FLOOR = datetime(2026, 7, 3, 0, 0, tzinfo=UTC)

    def test_empty_says_nothing_closed(self):
        out = c.render([], self.FLOOR)
        assert "nothing closed" in out

    def test_lists_rows_and_count(self):
        rows = [
            {"number": 20, "kind": "PR", "when": "2026-07-04T09:00:00Z", "title": "pr twenty"},
            {"number": 10, "kind": "Issue", "when": "2026-07-03T08:00:00Z", "title": "issue ten"},
        ]
        out = c.render(rows, self.FLOOR)
        assert "#20" in out and "#10" in out
        assert "2 item(s)." in out


class TestBothRepos:
    """Issue #1276 - the recap must cover BOTH repos.

    `REPO` was hard-coded to the project repo, so the entire Memory Manager
    workstream (which lives in the personas repo) was invisible to the
    mechanized recap. On 2026-07-22 a full two-day column of MM closes was
    silently omitted until it was noticed by hand.

    The load-bearing subtlety is NOT the loop over two repos - it is that the
    two repos have COLLIDING issue numbers. A recap that prints a bare `#98`
    is ambiguous between two different Issues, which is the same hazard the
    board tooling already solves with a `pers#` prefix. Rendering must
    disambiguate, or covering both repos trades a missing-data bug for a
    wrong-attribution bug.
    """

    PROJECT = "Jin-HoMLee/splice-neoepitope-pipeline"
    PERSONAS = "Jin-HoMLee/claude-personas-splice-neoepitope-pipeline"

    def test_default_repos_covers_both(self):
        assert c.PROJECT_REPO in c.DEFAULT_REPOS
        assert c.PERSONAS_REPO in c.DEFAULT_REPOS
        assert len(c.DEFAULT_REPOS) == 2

    def test_collect_tags_each_row_with_its_repo(self):
        floor = datetime(2026, 7, 22, 0, 0, tzinfo=UTC)
        rows = c.collect(
            floor,
            [{"number": 5, "title": "project issue", "closedAt": "2026-07-22T10:00:00Z"}],
            [],
            repo=self.PROJECT,
        )
        assert rows[0]["repo"] == self.PROJECT

    def test_rows_from_both_repos_merge_and_sort_newest_first(self):
        floor = datetime(2026, 7, 22, 0, 0, tzinfo=UTC)
        rows = c.collect(floor,
                         [{"number": 1, "title": "older project",
                           "closedAt": "2026-07-22T08:00:00Z"}],
                         [], repo=self.PROJECT)
        rows += c.collect(floor,
                          [{"number": 2, "title": "newer personas",
                            "closedAt": "2026-07-22T20:00:00Z"}],
                          [], repo=self.PERSONAS)
        merged = c.sort_rows(rows)
        assert [r["number"] for r in merged] == [2, 1], merged

    def test_render_disambiguates_colliding_numbers(self):
        """The whole point: same number, two repos, must not render identically.

        Without a repo-distinguishing prefix a reader cannot tell which #98 is
        meant - and #98 is a real live collision on the board today.
        """
        floor = datetime(2026, 7, 22, 0, 0, tzinfo=UTC)
        rows = c.sort_rows(
            c.collect(floor, [{"number": 98, "title": "project ninety-eight",
                               "closedAt": "2026-07-22T09:00:00Z"}], [],
                      repo=self.PROJECT)
            + c.collect(floor, [{"number": 98, "title": "personas ninety-eight",
                                 "closedAt": "2026-07-22T10:00:00Z"}], [],
                        repo=self.PERSONAS)
        )
        out = c.render(rows, floor)
        assert "pers#98" in out, out
        # The project row must NOT also carry the personas prefix.
        assert out.count("pers#98") == 1, out
        assert "#98 " in out.replace("pers#98", ""), out

    def test_render_labels_personas_rows_only(self):
        floor = datetime(2026, 7, 22, 0, 0, tzinfo=UTC)
        rows = c.sort_rows(c.collect(
            floor, [{"number": 7, "title": "project row",
                     "closedAt": "2026-07-22T09:00:00Z"}], [], repo=self.PROJECT))
        out = c.render(rows, floor)
        assert "pers#" not in out, out

    def test_repo_ref_prefix(self):
        assert c.repo_ref(98, self.PERSONAS) == "pers#98"
        assert c.repo_ref(98, self.PROJECT) == "#98"

    def test_unknown_repo_still_renders_a_usable_ref(self):
        """Fail-open: an unexpected repo must not crash a recap."""
        assert "98" in c.repo_ref(98, "Jin-HoMLee/some-other-repo")

    def test_two_repo_aggregation_preserves_the_inclusive_boundary(self):
        """The `>=` boundary must survive the two-repo path, in BOTH repos.

        The original off-by-one (Issue #784) was a bare `>` that excluded the
        whole boundary day. Aggregating across repos adds a new place that bug
        could reappear asymmetrically - e.g. a floor comparison applied to one
        repo's rows but not the other's. Both repos get an item at exactly
        00:00:00Z (the floor instant) and one a second BEFORE it; the boundary
        items must be kept and the pre-floor ones dropped.
        """
        floor = datetime(2026, 7, 22, 0, 0, tzinfo=UTC)
        rows = c.collect(
            floor,
            [{"number": 10, "title": "project at floor", "closedAt": "2026-07-22T00:00:00Z"},
             {"number": 11, "title": "project before floor", "closedAt": "2026-07-21T23:59:59Z"}],
            [], repo=self.PROJECT,
        ) + c.collect(
            floor,
            [{"number": 20, "title": "personas at floor", "closedAt": "2026-07-22T00:00:00Z"},
             {"number": 21, "title": "personas before floor", "closedAt": "2026-07-21T23:59:59Z"}],
            [], repo=self.PERSONAS,
        )
        kept = sorted(r["number"] for r in rows)
        assert kept == [10, 20], kept

    def test_two_repo_aggregation_keeps_merged_prs_from_both(self):
        """PR-side aggregation, since Issues and PRs take different code paths."""
        floor = datetime(2026, 7, 22, 0, 0, tzinfo=UTC)
        rows = c.collect(floor, [],
                         [{"number": 30, "title": "project pr", "mergedAt": "2026-07-22T05:00:00Z"}],
                         repo=self.PROJECT)
        rows += c.collect(floor, [],
                          [{"number": 40, "title": "personas pr", "mergedAt": "2026-07-22T06:00:00Z"}],
                          repo=self.PERSONAS)
        out = c.render(c.sort_rows(rows), floor)
        assert "pers#40" in out and "#30" in out, out
        assert all(r["kind"] == "PR" for r in rows)

    def test_ref_column_fits_a_five_digit_personas_ref(self):
        """Width must not shift the kind column once personas passes 9999.

        board_open_items.py sizes its ref column for `pers#12345`; the two tools
        should stay visually aligned rather than drifting apart at 5 digits.
        """
        floor = datetime(2026, 7, 22, 0, 0, tzinfo=UTC)
        rows = c.collect(floor, [{"number": 12345, "title": "big personas number",
                                  "closedAt": "2026-07-22T09:00:00Z"}], [],
                         repo=self.PERSONAS)
        line = [ln for ln in c.render(rows, floor).splitlines() if "pers#12345" in ln][0]
        assert "pers#12345 Issue" in line, line
