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
