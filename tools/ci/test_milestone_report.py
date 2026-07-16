"""Unit tests for scripts/pm/milestone_report.py — the metrics pure-functions.

Only the metrics layer is covered here (per the design spec §9): the metric
functions are pure (no I/O, no GraphQL, no jinja2/markdown), so they unit-test
cleanly with small synthetic fixtures. The data/aggregation/render layers are
exercised by the pm-i6 pilot dry-run, not by CI unit tests.

Importing this module must NOT pull in jinja2/markdown — those are lazy-imported
inside the render/aggregation functions so this test runs in the bare
``ci-tools-pytest`` env (pytest + pyyaml only).
"""

import ast
import sys
from datetime import datetime, timezone
from pathlib import Path

import pytest

# Make the script importable as a module (mirrors test_recheck_milestone.py)
SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts" / "pm"
sys.path.insert(0, str(SCRIPT_DIR))

import milestone_report as mr


# --- fixtures ---------------------------------------------------------------

def _issue(number, state, created_at, closed_at, roles, state_reason=None, labels=None):
    """Minimal normalized-issue dict carrying only what the metrics read.

    ``state_reason`` mirrors GitHub's ``stateReason`` (COMPLETED / NOT_PLANNED /
    None for open issues); a closed issue with no reason is treated as delivered.

    ``labels`` carries the raw label names; the arrival axis (Issue #811) reads the
    ``unplanned`` marker off it. Defaults to empty, so an issue with no labels is
    committed work - which is what every pre-#811 fixture in this file relies on.
    """
    return {
        "number": number,
        "state": state,
        "created_at": created_at,
        "closed_at": closed_at,
        "roles": roles,
        "state_reason": state_reason,
        "labels": labels if labels is not None else [],
    }


# A small synthetic milestone: 3 closed + 1 still open.
SAMPLE = [
    _issue(1, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-03T00:00:00Z", ["role:scientist"]),
    _issue(2, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-05T00:00:00Z", ["role:developer"]),
    _issue(3, "CLOSED", "2026-06-02T00:00:00Z", "2026-06-08T00:00:00Z", ["role:scientist", "role:developer"]),
    _issue(4, "OPEN", "2026-06-04T00:00:00Z", None, ["role:pm"]),
]


# --- slugify ----------------------------------------------------------------

class TestSlugify:
    def test_spec_example(self):
        assert mr.slugify("pm-i6 - PM Tooling & Memory Methodology II") == "pm-i6-pm-tooling-memory-methodology-ii"

    def test_stage_milestone(self):
        assert mr.slugify("i5 - S3 - Data Preparation") == "i5-s3-data-preparation"

    def test_collapses_and_strips(self):
        assert mr.slugify("  Foo   --  Bar!!  ") == "foo-bar"

    def test_no_leading_trailing_hyphen(self):
        s = mr.slugify("--- weird ---")
        assert not s.startswith("-") and not s.endswith("-")


# --- cycle time -------------------------------------------------------------

class TestLeadTime:
    """These were `TestCycleTime` and asserted `created -> closed` under that name.

    They measure LEAD time, and always did (Issue #1138) - the mislabelling lived in the
    tests as faithfully as in the code, which is why the tests never caught it. Renamed to
    what they actually measure. The SAMPLE fixture carries no `committed_at`, so under the
    corrected definitions its cycle times are all None - see TestCycleVsLeadTime.
    """

    def test_closed_issue_days(self):
        assert mr.lead_time_days(SAMPLE[0]) == 2.0
        assert mr.lead_time_days(SAMPLE[1]) == 4.0

    def test_open_issue_is_none(self):
        assert mr.lead_time_days(SAMPLE[3]) is None

    def test_missing_timestamp_is_none(self):
        assert mr.lead_time_days(_issue(9, "CLOSED", None, "2026-06-03T00:00:00Z", [])) is None

    def test_lead_times_only_closed(self):
        assert sorted(mr.lead_times(SAMPLE)) == [2.0, 4.0, 6.0]

    def test_avg_lead_time(self):
        assert mr.avg_lead_time(SAMPLE) == pytest.approx(4.0)  # (2+4+6)/3

    def test_median_lead_time(self):
        assert mr.median_lead_time(SAMPLE) == pytest.approx(4.0)

    def test_avg_empty_is_none(self):
        assert mr.avg_lead_time([]) is None
        assert mr.median_lead_time([]) is None

    def test_uncommitted_sample_has_no_cycle_time_at_all(self):
        """The corrected contract, stated on the legacy fixture: no commitment act
        recorded -> cycle time is undefined for every one of them."""
        assert mr.cycle_times(SAMPLE) == []
        assert mr.avg_cycle_time(SAMPLE) is None


# --- per-role counts --------------------------------------------------------

class TestPerRoleCounts:
    def test_counts_closed_only_multi_role(self):
        # closed: #1 sci, #2 dev, #3 sci+dev ; open #4 pm excluded
        assert mr.per_role_counts(SAMPLE) == {"role:scientist": 2, "role:developer": 2}

    def test_empty(self):
        assert mr.per_role_counts([]) == {}


# --- close-reason split (delivered vs descoped) -----------------------------

# A milestone mixing COMPLETED + NOT_PLANNED + DUPLICATE closes + one open issue.
# Descoped = NOT_PLANNED (#2) and DUPLICATE (#5); delivered = #1, #3.
MIXED = [
    _issue(1, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-03T00:00:00Z", ["role:scientist"], "COMPLETED"),
    _issue(2, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-05T00:00:00Z", ["role:developer"], "NOT_PLANNED"),
    _issue(3, "CLOSED", "2026-06-02T00:00:00Z", "2026-06-08T00:00:00Z", ["role:scientist", "role:developer"], "COMPLETED"),
    _issue(4, "OPEN", "2026-06-04T00:00:00Z", None, ["role:pm"]),
    _issue(5, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-09T00:00:00Z", ["role:developer"], "DUPLICATE"),
]


class TestCloseReasonSplit:
    def test_delivered_excludes_descoped_and_open(self):
        # NOT_PLANNED (#2), DUPLICATE (#5), and open (#4) all excluded.
        nums = [i["number"] for i in mr.delivered_issues(MIXED)]
        assert nums == [1, 3]

    def test_descoped_includes_not_planned_and_duplicate(self):
        # DUPLICATE is a descope (tracked elsewhere), not a delivery — must not
        # be swept into delivered by a bare "!= NOT_PLANNED" complement.
        nums = sorted(i["number"] for i in mr.descoped_issues(MIXED))
        assert nums == [2, 5]

    def test_closed_still_counts_all_closed(self):
        # closed_issues stays the union (delivered + descoped), unchanged.
        assert sorted(i["number"] for i in mr.closed_issues(MIXED)) == [1, 2, 3, 5]

    def test_missing_reason_on_closed_is_delivered(self):
        # Legacy/None reason on a closed issue must not be dropped from delivered.
        assert mr.delivered_issues(SAMPLE) == mr.closed_issues(SAMPLE)
        assert mr.descoped_issues(SAMPLE) == []

    def test_per_role_counts_excludes_descoped(self):
        # #2 + #5 (developer) are descoped -> excluded; only #1 sci + #3 sci+dev.
        assert mr.per_role_counts(MIXED) == {"role:scientist": 2, "role:developer": 1}


# --- duration & throughput --------------------------------------------------

class TestDurationThroughput:
    def test_duration_earliest_created_to_milestone_close(self):
        # earliest created 2026-06-01 -> milestone closed 2026-06-08 = 7 days
        assert mr.milestone_duration_days(SAMPLE, "2026-06-08T00:00:00Z") == 7.0

    def test_duration_uses_milestone_created_fallback(self):
        # if no issues, fall back to milestone_created_at -> closed
        assert mr.milestone_duration_days([], "2026-06-08T00:00:00Z", "2026-06-01T00:00:00Z") == 7.0

    def test_duration_none_when_unknowable(self):
        assert mr.milestone_duration_days([], None) is None

    def test_throughput_per_week(self):
        # 3 closed over 7 days = 3 per week
        assert mr.throughput_per_week(3, 7.0) == pytest.approx(3.0)

    def test_throughput_zero_duration_guarded(self):
        assert mr.throughput_per_week(3, 0.0) is None
        assert mr.throughput_per_week(3, None) is None


# --- compute_metrics (integration round-trip) -------------------------------

class TestComputeMetrics:
    def test_dict_shape_and_values(self):
        m = mr.compute_metrics(
            SAMPLE,
            {"created_at": "2026-06-01T00:00:00Z", "closed_at": "2026-06-08T00:00:00Z"},
            marker_in_use=True,
        )
        assert m["n_total"] == 4
        assert m["n_closed"] == 3
        assert m["n_delivered"] == 3  # no state_reason -> all delivered
        assert m["n_descoped"] == 0
        assert m["n_carried_forward"] == 1
        assert m["duration_days"] == 7.0
        assert m["throughput_per_week"] == pytest.approx(3.0)
        assert m["avg_lead_time_days"] == pytest.approx(4.0)
        assert m["avg_cycle_time_days"] is None   # fixture has no commitment act
        assert m["median_lead_time_days"] == pytest.approx(4.0)
        assert m["per_role_counts"] == {"role:scientist": 2, "role:developer": 2}

    def test_split_counts_and_delivered_throughput(self):
        # MIXED: 2 delivered (#1,#3) + 2 descoped (#2 NOT_PLANNED, #5 DUPLICATE) + 1 open.
        m = mr.compute_metrics(
            MIXED,
            {"created_at": "2026-06-01T00:00:00Z", "closed_at": "2026-06-08T00:00:00Z"},
            marker_in_use=True,
        )
        assert m["n_closed"] == 4
        assert m["n_delivered"] == 2
        assert m["n_descoped"] == 2
        # throughput keys off *delivered*, not raw closed: 2 over 7 days.
        assert m["throughput_per_week"] == pytest.approx(2.0)
        assert m["per_role_counts"] == {"role:scientist": 2, "role:developer": 1}
        # cycle time is over delivered only: #1=2d, #3=6d -> avg/median 4.0.
        # The descoped #2 (4d) and #5 (8d) must not skew it.
        assert m["avg_lead_time_days"] == pytest.approx(4.0)
        assert m["avg_cycle_time_days"] is None   # fixture has no commitment act
        assert m["median_lead_time_days"] == pytest.approx(4.0)

    def test_partition_invariant(self):
        # delivered + descoped always exactly partitions closed.
        m = mr.compute_metrics(
            MIXED, {"created_at": None, "closed_at": "2026-06-08T00:00:00Z"},
            marker_in_use=True,
        )
        assert m["n_delivered"] + m["n_descoped"] == m["n_closed"]

    def test_empty_milestone_degrades(self):
        m = mr.compute_metrics([], {"created_at": None, "closed_at": None}, marker_in_use=True)
        assert m["n_total"] == m["n_closed"] == m["n_carried_forward"] == 0
        assert m["n_delivered"] == m["n_descoped"] == 0
        assert m["duration_days"] is None
        assert m["throughput_per_week"] is None
        assert m["avg_cycle_time_days"] is None
        assert m["per_role_counts"] == {}


# --- window-mode metrics (weekly SDR, Issue #915) ---------------------------

UNTIL = datetime(2026, 7, 3, tzinfo=timezone.utc)

# Trailing-4-week fixture keyed off UNTIL (2026-07-03). Windows (start, end]:
#   w1 (06-05, 06-12] · w2 (06-12, 06-19] · w3 (06-19, 06-26] · w4 (06-26, 07-03] = reporting week
WINDOW_ISSUES = [
    _issue(1, "CLOSED", "2026-06-29T00:00:00Z", "2026-07-01T00:00:00Z", ["role:pm"], "COMPLETED"),                 # w4 delivered, 2d
    _issue(2, "CLOSED", "2026-06-20T00:00:00Z", "2026-06-27T00:00:00Z", ["role:developer"], "COMPLETED"),          # w4 delivered, 7d
    _issue(4, "CLOSED", "2026-06-30T00:00:00Z", "2026-07-02T00:00:00Z", ["role:pm"], "NOT_PLANNED"),               # w4 descoped
    _issue(3, "CLOSED", "2026-06-18T00:00:00Z", "2026-06-20T00:00:00Z", ["role:scientist"], "COMPLETED"),          # w3 delivered, 2d
    _issue(5, "CLOSED", "2026-06-08T00:00:00Z", "2026-06-10T00:00:00Z", ["role:pm"], "COMPLETED"),                 # w1 delivered, 2d
]


class TestWeekWindows:
    def test_shape_and_chronology(self):
        w = mr.week_windows(UNTIL, 4)
        assert len(w) == 4
        # chronological (oldest first); reporting week (most recent) is last.
        starts = [s for s, _ in w]
        assert starts == sorted(starts)
        assert w[-1][1] == UNTIL                       # reporting week ends at until
        assert (w[-1][1] - w[-1][0]).days == 7         # each window is 7 days

    def test_n_weeks_floor_of_one(self):
        assert len(mr.week_windows(UNTIL, 0)) == 1     # never zero windows


class TestClosedInWindow:
    def test_filters_by_closed_at_half_open(self):
        start = datetime(2026, 6, 26, tzinfo=timezone.utc)
        end = UNTIL
        nums = sorted(i["number"] for i in mr.closed_in_window(WINDOW_ISSUES, start, end))
        assert nums == [1, 2, 4]                        # closed in (06-26, 07-03]

    def test_excludes_boundary_start_includes_end(self):
        # (start, end]: an issue closed exactly at start is excluded; at end included.
        issues = [
            _issue(10, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-26T00:00:00Z", []),  # == start -> out
            _issue(11, "CLOSED", "2026-06-01T00:00:00Z", "2026-07-03T00:00:00Z", []),  # == end -> in
        ]
        nums = [i["number"] for i in mr.closed_in_window(
            issues, datetime(2026, 6, 26, tzinfo=timezone.utc), UNTIL)]
        assert nums == [11]


class TestWeeklySeries:
    def test_trend_buckets_and_medians(self):
        s = mr.weekly_series(WINDOW_ISSUES, UNTIL, 4, marker_in_use=True)
        assert [w["n_delivered"] for w in s] == [1, 0, 1, 2]   # w1, w2(empty), w3, w4
        assert s[1]["median_lead_time_days"] is None          # empty week -> None
        assert s[0]["median_lead_time_days"] == pytest.approx(2.0)   # #5
        assert s[0]["median_cycle_time_days"] is None   # no commitment act
        assert s[3]["median_lead_time_days"] == pytest.approx(4.5)   # #1=2d, #2=7d
        assert s[3]["median_cycle_time_days"] is None   # no commitment act on the fixture
        # descoped #4 must not inflate the reporting-week bucket count.
        assert s[3]["week_end"] == "2026-07-03"


class TestComputeWindowMetrics:
    def test_headline_over_reporting_week_plus_trend(self):
        week_issues = mr.closed_in_window(
            WINDOW_ISSUES, datetime(2026, 6, 26, tzinfo=timezone.utc), UNTIL)
        m = mr.compute_window_metrics(week_issues, WINDOW_ISSUES, UNTIL, 4, marker_in_use=True)
        assert m["n_total"] == 3            # #1, #2, #4 closed in the reporting week
        assert m["n_delivered"] == 2        # #1, #2 (descoped #4 excluded)
        assert m["n_descoped"] == 1         # #4 NOT_PLANNED
        assert m["n_carried_forward"] == 0  # window mode has no carried concept
        assert m["throughput_per_week"] == pytest.approx(2.0)  # 2 delivered / 1 week
        assert m["avg_lead_time_days"] == pytest.approx(4.5)  # (2 + 7) / 2
        assert m["avg_cycle_time_days"] is None   # fixture has no commitment act
        assert m["per_role_counts"] == {"role:pm": 1, "role:developer": 1}  # delivered only
        assert [w["n_delivered"] for w in m["weekly_series"]] == [1, 0, 1, 2]

    def test_zero_ship_week_is_empty_headline(self):
        # No issues closed in the reporting week -> n_total 0 (main() skips on this).
        old = [_issue(20, "CLOSED", "2026-06-08T00:00:00Z", "2026-06-10T00:00:00Z", ["role:pm"], "COMPLETED")]
        week_issues = mr.closed_in_window(
            old, datetime(2026, 6, 26, tzinfo=timezone.utc), UNTIL)
        m = mr.compute_window_metrics(week_issues, old, UNTIL, 4, marker_in_use=True)
        # main() gates the zero-ship skip on n_total == 0 (nothing closed in the
        # reporting week); delivered is 0 and cycle time is undefined.
        assert m["n_total"] == 0
        assert m["n_delivered"] == 0
        assert m["avg_cycle_time_days"] is None


# --- narrative seed digest --------------------------------------------------

class TestFirstProseLine:
    def test_skips_byline_and_timestamp_subheader(self):
        # PM notebook shape: date -> "### HH:MM UTC — Editor: PM" -> "#### <title>"
        body = "\n### 17:53 UTC — Editor: PM\n\n#### Morning routine → governance deep-dive\n\n**Session.** ..."
        assert mr._first_prose_line(body) == "Morning routine → governance deep-dive"

    def test_plain_prose_first_line(self):
        # Developer notebook shape: prose immediately under the date header.
        assert mr._first_prose_line("Refactored the config loader and added tests.") \
            == "Refactored the config loader and added tests."

    def test_skips_html_comment(self):
        assert mr._first_prose_line("<!-- note -->\nReal content here.") == "Real content here."


class TestSeedNarrative:
    ISSUES = [
        {"number": 5, "title": "do a thing", "url": "http://x/5", "state": "CLOSED", "roles": []},
        {"number": 6, "title": "open thing", "url": "http://x/6", "state": "OPEN", "roles": []},
    ]

    def _seed(self):
        return mr.seed_narrative(
            {"title": "m", "created_at": None, "closed_at": None}, self.ISSUES, {}
        )

    def test_issues_emitted_as_links_not_bare_hash(self):
        # Python-Markdown parses a leading "#N" (e.g. "- #6 …") as an <h1>, which
        # blows up the rendered font size. Issues must be emitted as md links so
        # the "#" never sits at list-item-content start.
        md = self._seed()
        assert not any(line.startswith("- #") for line in md.splitlines())
        assert "[#6](http://x/6)" in md   # open -> Carried-forward, as a link

    def test_deliverables_seed_does_not_relist_closed(self):
        # Slimmed seed: closed issues live in the Inventory appendix, not
        # re-listed in the Deliverables narrative (avoids duplicating the table).
        deliverables = self._seed().split("## Carried-forward")[0]
        assert "#5" not in deliverables          # closed #5 not re-listed
        assert "Inventory appendix" in deliverables  # points there instead


# --- arrival axis: committed vs unplanned (Issue #811) -----------------------

# The matched pair the AC asks for: same window, same delivered status, one
# variable flipped (the `unplanned` marker) and opposite expected sides of the
# split. Plus a descoped-and-unplanned issue, which must land on NEITHER side -
# the breakdown keys off *delivered*, so a not-planned close is not shipped work
# no matter how it arrived.
ARRIVAL = [
    _issue(1, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-03T00:00:00Z",
           ["role:pm"], "COMPLETED"),                                    # committed
    _issue(2, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-03T00:00:00Z",
           ["role:pm"], "COMPLETED", ["unplanned"]),                     # unplanned
    _issue(3, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-04T00:00:00Z",
           ["role:dev"], "NOT_PLANNED", ["unplanned"]),                  # descoped: neither
    _issue(4, "OPEN", "2026-06-02T00:00:00Z", None, ["role:pm"]),        # open: neither
]


class TestArrivalAxis:
    def test_matched_pair_lands_on_opposite_sides(self):
        # Identical but for the marker -> opposite sides. This is the control:
        # if the marker were ignored, both would land as committed and this fails.
        b = mr.throughput_breakdown(ARRIVAL, marker_in_use=True)
        assert b["n_committed"] == 1
        assert b["n_unplanned"] == 1

    def test_descoped_unplanned_counts_on_neither_side(self):
        # #3 carries the marker but closed NOT_PLANNED. It is not shipped work,
        # so it must not inflate the unplanned count - otherwise the "unplanned
        # share" would be measuring abandoned work, not absorbed capacity.
        b = mr.throughput_breakdown(ARRIVAL, marker_in_use=True)
        assert b["n_committed"] + b["n_unplanned"] == len(mr.delivered_issues(ARRIVAL)) == 2

    def test_pct_unplanned(self):
        assert mr.throughput_breakdown(ARRIVAL, marker_in_use=True)["pct_unplanned"] == pytest.approx(50.0)

    def test_pct_is_none_not_zero_on_empty_window(self):
        # A zero-ship week must read as "no data", not as a truthful-looking 0%
        # unplanned - the latter would be a number the WIP retune could act on.
        assert mr.throughput_breakdown([], marker_in_use=True)["pct_unplanned"] is None

    def test_unlabelled_issue_is_committed(self):
        # Guards the back-compat assumption every pre-#811 fixture leans on.
        assert mr.is_unplanned(_issue(9, "CLOSED", None, None, [])) is False

    def test_priority_cannot_substitute_for_the_marker(self):
        # The whole premise of #811: urgency and arrival are orthogonal. A P2
        # unplanned fix and a P1 committed item are indistinguishable by band,
        # so the split must key off the marker and nothing else.
        p1_committed = _issue(10, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-02T00:00:00Z",
                              ["role:pm"], "COMPLETED", ["P1"])
        p2_unplanned = _issue(11, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-02T00:00:00Z",
                              ["role:pm"], "COMPLETED", ["P2", "unplanned"])
        b = mr.throughput_breakdown([p1_committed, p2_unplanned], marker_in_use=True)
        assert (b["n_committed"], b["n_unplanned"]) == (1, 1)

    def test_metrics_carry_the_split(self):
        m = mr.compute_metrics(ARRIVAL, {"closed_at": None, "created_at": None}, marker_in_use=True)
        assert m["n_committed"] == 1 and m["n_unplanned"] == 1
        assert m["pct_unplanned"] == pytest.approx(50.0)

    def test_weekly_series_carries_the_split(self):
        until = datetime(2026, 6, 8, tzinfo=timezone.utc)
        series = mr.weekly_series(ARRIVAL, until, 1, marker_in_use=True)
        assert series[0]["n_committed"] == 1
        assert series[0]["n_unplanned"] == 1
        # Conservation: the split must exactly partition the delivered count.
        assert series[0]["n_committed"] + series[0]["n_unplanned"] == series[0]["n_delivered"]


# --- the marker's OWN falsifier (Issue #1180) --------------------------------

# A window of delivered work in which NOT ONE issue carries the marker. From this
# list alone, two completely different worlds are information-theoretically
# identical:
#   (a) the marker IS in use, and we genuinely absorbed no unplanned work -> a real 0%
#   (b) the marker was NEVER applied to anything                          -> 0% is a fiction
# The list cannot tell them apart. Only the injected `marker_in_use` fact can.
NO_MARKER_ANYWHERE = [
    _issue(1, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-03T00:00:00Z",
           ["role:pm"], "COMPLETED"),
    _issue(2, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-03T00:00:00Z",
           ["role:dev"], "COMPLETED"),
]

# Built from codepoints, not typed: this file must CONTAIN no em/en dash (the repo
# guard rejects one on sight) while still asserting their absence elsewhere.
EM_DASH = chr(0x2014)
EN_DASH = chr(0x2013)


class TestMarkerNotInUse:
    """The metric's own falsifier: a 0% that could only ever come back one way.

    Before Issue #1180 the SDR printed a clean, confident `25 / 0 (0% unplanned)` on
    every weekly row - while `gh issue list --label unplanned --state all` returned
    ZERO issues repo-wide. The label had never been applied to anything. The
    measurement shipped; the thing being measured did not.
    """

    def test_matched_pair_real_zero_vs_fabricated_zero(self):
        """THE test. Identical input, one fact flipped, opposite outputs.

        Without this pair the fix is unfalsifiable - which is precisely the disease
        being cured, so shipping the cure without it would be self-refuting.
        """
        real_zero = mr.throughput_breakdown(NO_MARKER_ANYWHERE, marker_in_use=True)
        fabricated = mr.throughput_breakdown(NO_MARKER_ANYWHERE, marker_in_use=False)

        # In use + nobody was unplanned -> a genuine, hard-won zero. Report it.
        assert real_zero["pct_unplanned"] == pytest.approx(0.0)
        assert real_zero["n_committed"] == 2
        assert real_zero["n_unclassifiable"] == 0

        # Not in use -> the SAME list means nothing. Refuse to print a number.
        assert fabricated["pct_unplanned"] is None
        assert fabricated["n_committed"] is None, (
            "'2 committed' is the same false assertion as '0% unplanned', one column over"
        )
        assert fabricated["n_unplanned"] is None
        assert fabricated["n_unclassifiable"] == 2, "the degraded input must be VISIBLE"

    def test_marker_in_use_is_zero_vs_nonzero_on_the_repo_count(self):
        assert mr.is_marker_in_use(0) is False
        assert mr.is_marker_in_use(1) is True

    def test_the_two_none_worlds_stay_distinguishable(self):
        """`pct_unplanned is None` now has two causes; they must not collapse.

        Empty window (no delivered work) and marker-not-in-use both yield None, but
        they are different failures and the report says different things about them -
        so the `marker_in_use` flag has to survive into the metrics dict.
        """
        empty_window = mr.throughput_breakdown([], marker_in_use=True)
        assert empty_window["pct_unplanned"] is None
        assert empty_window["marker_in_use"] is True
        assert empty_window["n_committed"] == 0     # honestly zero: there IS no work

        not_in_use = mr.throughput_breakdown(NO_MARKER_ANYWHERE, marker_in_use=False)
        assert not_in_use["pct_unplanned"] is None
        assert not_in_use["marker_in_use"] is False
        assert not_in_use["n_committed"] is None    # unknowable: there IS work, unclassified

    def test_marker_in_use_is_a_required_keyword(self):
        """A default would let a new call site silently fabricate. It must decide."""
        with pytest.raises(TypeError):
            mr.throughput_breakdown(NO_MARKER_ANYWHERE)

    def test_weekly_series_propagates_unclassifiable(self):
        """The fabricated 0% appeared in EVERY trend row, not just the headline card."""
        until = datetime(2026, 6, 8, tzinfo=timezone.utc)
        series = mr.weekly_series(NO_MARKER_ANYWHERE, until, 1, marker_in_use=False)
        assert series[0]["pct_unplanned"] is None
        assert series[0]["n_committed"] is None
        assert series[0]["n_unclassifiable"] == 2

    def test_console_trend_row_never_leaks_a_none_sentinel(self):
        """The RENDERED row, not just the metrics dict.

        `n_committed=None` is honest in the data, but the console loop f-stringed it
        straight out as the literal `[None+None]` on every row. Asserting on the dict
        (the test above) passes happily while the human-facing output is garbage - so
        this one asserts on what a person actually reads.

        It survived my own live check because that check grepped the output for
        `arrival` - the line I had just fixed. A verification aimed only at what you
        changed cannot show you what you missed.
        """
        until = datetime(2026, 6, 8, tzinfo=timezone.utc)
        series = mr.weekly_series(NO_MARKER_ANYWHERE, until, 1, marker_in_use=False)

        rendered = mr._fmt_arrival(series[0])
        assert "None" not in rendered, "internal sentinel leaked into a human-facing row"
        assert rendered == "unclassifiable"

    def test_console_trend_row_still_shows_the_split_when_in_use(self):
        """Matched pair: same row, marker in use -> the real split, not 'unclassifiable'."""
        until = datetime(2026, 6, 8, tzinfo=timezone.utc)
        series = mr.weekly_series(ARRIVAL, until, 1, marker_in_use=True)
        assert mr._fmt_arrival(series[0]) == "1+1"

    def test_compute_metrics_propagates_unclassifiable(self):
        m = mr.compute_metrics(NO_MARKER_ANYWHERE, {"closed_at": None, "created_at": None},
                               marker_in_use=False)
        assert m["marker_in_use"] is False
        assert m["pct_unplanned"] is None
        assert m["n_unclassifiable"] == 2
        # The delivered headline is unaffected - we know WHAT shipped, not how it arrived.
        assert m["n_delivered"] == 2


class TestCycleVsLeadTime:
    """Cycle time starts at the COMMITMENT POINT; lead time starts at creation.

    The old metric computed created -> closed and called it "cycle time" (Issue #1138).
    That charges uncommitted Backlog dwell to delivery performance: an option can rest in
    Backlog for six weeks (which late-commitment Kanban actively PRESCRIBES), be
    committed, and ship in two days - and the old number called that a 44-day cycle time.
    It penalized the exact behavior the model is built on, and fed that to the WIP retune.
    """

    def test_the_44_day_scenario_the_issue_describes(self):
        """The establishing case, as a number: 6 weeks resting + 2 days working."""
        issue = _issue(1, "CLOSED", "2026-05-01T00:00:00Z", "2026-06-14T00:00:00Z",
                       ["role:pm"], "COMPLETED")
        issue["committed_at"] = "2026-06-12T00:00:00Z"

        assert mr.lead_time_days(issue) == pytest.approx(44.0)   # what it used to report
        assert mr.cycle_time_days(issue) == pytest.approx(2.0)   # what it actually took

    def test_never_committed_cycle_time_is_undefined_not_zero(self):
        """An item that never crossed into Ready has NO cycle time.

        None, never 0.0. A zero would assert an instantaneous delivery; and silently
        dropping it from the mean would bias the average toward exactly the committed
        work - the same "can only come back one way" disease as the #1180 fabricated 0%.
        """
        issue = _issue(1, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-03T00:00:00Z",
                       ["role:pm"], "COMPLETED")   # no committed_at at all
        assert mr.cycle_time_days(issue) is None
        assert mr.lead_time_days(issue) == pytest.approx(2.0), "lead time is still defined"

    def test_never_committed_are_counted_not_dropped(self):
        committed = _issue(1, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-05T00:00:00Z",
                           ["role:pm"], "COMPLETED")
        committed["committed_at"] = "2026-06-03T00:00:00Z"
        barged = _issue(2, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-05T00:00:00Z",
                        ["role:pm"], "COMPLETED")   # never committed

        assert [i["number"] for i in mr.never_committed_issues([committed, barged])] == [2]
        # The cycle-time mean uses only the item that HAS a cycle time...
        assert mr.avg_cycle_time([committed, barged]) == pytest.approx(2.0)
        # ...which is exactly why the excluded one must be counted and surfaced.

    def test_descoped_never_counted_as_never_committed(self):
        """never_committed keys off DELIVERED, so a descoped close is not shipped work."""
        descoped = _issue(3, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-05T00:00:00Z",
                          ["role:pm"], "NOT_PLANNED")
        assert mr.never_committed_issues([descoped]) == []

    def test_open_issue_has_neither(self):
        opened = _issue(4, "OPEN", "2026-06-01T00:00:00Z", None, ["role:pm"])
        opened["committed_at"] = "2026-06-02T00:00:00Z"
        assert mr.cycle_time_days(opened) is None
        assert mr.lead_time_days(opened) is None


class TestCommitmentFetchFailureIsUnknownNotZero:
    """A failed commitment fetch must read as UNKNOWN, never as "nobody committed".

    Caught in review of PR #1185, and it is the sharpest possible instance: on a GraphQL
    error `fetch_commitment_times` returned {}, so every item had no `committed_at`, and
    the report then CONFIDENTLY asserted that the entire delivered population "never
    crossed Backlog -> Ready" - contradicted only by a stderr line the reader never sees.

    One output, three different worlds (nobody committed / we could not read the history /
    the reader is being lied to). That is the exact fabrication #1180 exists to kill,
    reintroduced one layer down, inside the PR that kills it.
    """

    def _delivered(self, n, committed=None, available=True):
        i = _issue(n, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-05T00:00:00Z",
                   ["role:pm"], "COMPLETED")
        if committed:
            i["committed_at"] = committed
        i["commitments_available"] = available
        return i

    def test_matched_pair_genuinely_uncommitted_vs_fetch_failed(self):
        """THE control: identical issue lists, one flag flipped, opposite meanings."""
        genuinely = [self._delivered(1, available=True),
                     self._delivered(2, available=True)]
        fetch_died = [self._delivered(1, available=False),
                      self._delivered(2, available=False)]

        # We CAN see the history, and neither crossed Ready -> a real finding.
        assert len(mr.never_committed_issues(genuinely)) == 2

        # We CANNOT see the history -> we know nothing. Asserting "2 never committed"
        # here would state the opposite of the truth.
        assert mr.never_committed_issues(fetch_died) == []
        assert mr.commitments_available(fetch_died) is False

    def test_metrics_render_never_committed_as_none_when_unavailable(self):
        m = mr.compute_metrics(
            [self._delivered(1, available=False)],
            {"closed_at": None, "created_at": None}, marker_in_use=True)
        assert m["commitments_available"] is False
        assert m["n_never_committed"] is None, "a count here would be a fabricated caveat"
        assert m["avg_cycle_time_days"] is None, "cycle time is unknown too, not zero"

    def test_metrics_report_the_count_when_the_fetch_worked(self):
        """Matched pair to the above: same shape, fetch succeeded -> report the number."""
        m = mr.compute_metrics(
            [self._delivered(1, available=True)],
            {"closed_at": None, "created_at": None}, marker_in_use=True)
        assert m["commitments_available"] is True
        assert m["n_never_committed"] == 1

    def test_availability_defaults_true_for_plain_fixtures(self):
        """Absent flag = available, so pure fixtures behave as before."""
        assert mr.commitments_available([_issue(1, "CLOSED", "a", "b", [])]) is True

    def test_one_unavailable_item_taints_the_whole_report(self):
        """Partial commitment data is deliberately NOT used.

        A half-populated map makes the items from the failed chunk indistinguishable from
        genuinely-uncommitted ones - the same conflation, one chunk smaller.
        """
        mixed = [self._delivered(1, available=True), self._delivered(2, available=False)]
        assert mr.commitments_available(mixed) is False


class TestMarkerSlip:
    """The hand-applied `unplanned` LABEL cross-checked against the OBSERVED arrival
    signal (never crossed Backlog -> Ready). Issue #1188.

    Two independent measures of the same fact - "did this delivered item arrive
    unplanned?":
      - observed:  it has no ``committed_at`` (never entered `Ready`) - recorded by
                   GitHub for free, cannot be forgotten.
      - labelled:  someone applied the ``unplanned`` label at close - hand-maintained.

    A forgotten label is indistinguishable from committed work, so the label biases the
    unplanned share DOWNWARD. Their DISAGREEMENT is the deliverable:
      - ``n_slip``:        observed-unplanned AND unlabelled (the under-report #1144 predicted).
      - ``n_mislabelled``: labelled AND observably committed (the opposite error).
      - ``n_agree``:       observed-unplanned AND labelled.
    """

    def _d(self, n, committed=None, unplanned=False, available=True):
        labels = ["unplanned"] if unplanned else []
        i = _issue(n, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-05T00:00:00Z",
                   ["role:pm"], "COMPLETED", labels)
        if committed:
            i["committed_at"] = committed
        i["commitments_available"] = available
        return i

    def test_matched_pair_agree_vs_slip(self):
        """THE control: two observably-unplanned items, ONLY the label flipped.

        Both never crossed `Ready`, so both are observably unplanned. One carries the
        label (agree), one does not (the slip). If the cross-check ignored the label,
        both would land the same and this fails.
        """
        agree = self._d(1, committed=None, unplanned=True)
        slip = self._d(2, committed=None, unplanned=False)
        m = mr.marker_slip([agree, slip])
        assert m["available"] is True
        assert m["n_observed_unplanned"] == 2   # both never committed
        assert m["n_labelled_unplanned"] == 1   # only #1 carries the label
        assert m["n_slip"] == 1                  # #2: observed-unplanned but unlabelled
        assert m["n_agree"] == 1                 # #1: observed-unplanned and labelled
        assert m["n_mislabelled"] == 0

    def test_mislabelled_is_labelled_but_observably_committed(self):
        """The symmetric error: the label is applied, but the item DID cross the commitment act."""
        mislabelled = self._d(1, committed="2026-06-02T00:00:00Z", unplanned=True)
        clean = self._d(2, committed="2026-06-02T00:00:00Z", unplanned=False)
        m = mr.marker_slip([mislabelled, clean])
        assert m["n_observed_unplanned"] == 0   # both committed
        assert m["n_mislabelled"] == 1          # #1: labelled yet committed
        assert m["n_slip"] == 0
        assert m["n_agree"] == 0

    def test_aggregate_share_agreement_does_not_mean_per_item_fidelity(self):
        """Why the cross-check earns its keep: the two shares can MATCH while the label
        is wrong on every item.

        4 delivered: #1 observed-unplanned+labelled (agree), #2 observed-unplanned+unlabelled
        (slip), #3 committed+labelled (mislabelled), #4 committed+unlabelled (clean). Observed
        share = 2/4 = labelled share = 2/4 - identical in aggregate, yet two of the four items
        are classified WRONGLY by the label. A share comparison alone would call this clean.
        """
        issues = [
            self._d(1, committed=None, unplanned=True),
            self._d(2, committed=None, unplanned=False),
            self._d(3, committed="2026-06-02T00:00:00Z", unplanned=True),
            self._d(4, committed="2026-06-02T00:00:00Z", unplanned=False),
        ]
        m = mr.marker_slip(issues)
        assert m["n_delivered"] == 4
        assert m["observed_pct"] == pytest.approx(50.0)
        assert m["labelled_pct"] == pytest.approx(50.0)
        assert m["n_slip"] == 1        # #2
        assert m["n_mislabelled"] == 1 # #3
        assert m["n_agree"] == 1       # #1

    def test_marker_never_applied_is_total_slip(self):
        """The live 2026-07 finding, as a fixture: observed-unplanned items exist, but NOT
        ONE carries the label -> every one is a slip and the label reports zero."""
        issues = [self._d(1, committed=None, unplanned=False),
                  self._d(2, committed=None, unplanned=False),
                  self._d(3, committed="2026-06-02T00:00:00Z", unplanned=False)]
        m = mr.marker_slip(issues)
        assert m["n_observed_unplanned"] == 2
        assert m["n_labelled_unplanned"] == 0   # the label reports zero...
        assert m["n_slip"] == 2                  # ...while two observably slipped
        assert m["n_agree"] == 0

    def test_unavailable_commitments_is_unknown_never_folded(self):
        """Matched pair to ``test_matched_pair_agree_vs_slip``: same issues, but the
        commitment fetch failed. The observed axis is then unknown, so slip is unknown -
        NOT zero, NOT the whole population (the #1180 / #1185 fabrication class)."""
        died = [self._d(1, committed=None, unplanned=True, available=False),
                self._d(2, committed=None, unplanned=False, available=False)]
        m = mr.marker_slip(died)
        assert m["available"] is False
        assert m["n_observed_unplanned"] is None
        assert m["n_labelled_unplanned"] is None
        assert m["n_slip"] is None
        assert m["n_mislabelled"] is None
        assert m["n_agree"] is None
        assert m["observed_pct"] is None
        assert m["labelled_pct"] is None
        # delivered is still known - we know WHAT shipped, not how it arrived.
        assert m["n_delivered"] == 2

    def test_descoped_never_counted(self):
        """Keys off DELIVERED: a descoped close is not shipped work, however it arrived."""
        descoped = _issue(9, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-05T00:00:00Z",
                          ["role:pm"], "NOT_PLANNED", ["unplanned"])
        descoped["commitments_available"] = True
        m = mr.marker_slip([descoped])
        assert m["n_delivered"] == 0
        assert m["n_observed_unplanned"] == 0
        assert m["n_slip"] == 0

    def test_metrics_carry_the_cross_check(self):
        agree = self._d(1, committed=None, unplanned=True)
        slip = self._d(2, committed=None, unplanned=False)
        m = mr.compute_metrics([agree, slip], {"closed_at": None, "created_at": None},
                               marker_in_use=True)
        assert m["marker_slip"]["n_slip"] == 1
        assert m["marker_slip"]["n_agree"] == 1

    def test_window_metrics_carry_the_cross_check(self):
        agree = self._d(1, committed=None, unplanned=True)
        slip = self._d(2, committed=None, unplanned=False)
        week = [agree, slip]
        m = mr.compute_window_metrics(week, week, UNTIL, 1, marker_in_use=True)
        assert m["marker_slip"]["n_slip"] == 1


class TestNegativeCycleTimeIsDropped:
    def test_committed_after_closed_is_not_a_fast_delivery(self):
        """A close -> reopen -> re-commit can put committed_at AFTER closed_at.

        A negative cycle time is not a fast delivery, it is a nonsense one, and left in it
        would drag the mean below zero. Dropped, not clamped to 0 (a 0 would assert an
        instantaneous delivery - the same lie in the other direction).
        """
        i = _issue(1, "CLOSED", "2026-06-01T00:00:00Z", "2026-06-05T00:00:00Z",
                   ["role:pm"], "COMPLETED")
        i["committed_at"] = "2026-06-09T00:00:00Z"   # after the close
        assert mr.cycle_time_days(i) is None
        assert mr.avg_cycle_time([i]) is None


class TestFirstReadyAt:
    """Extracting the commitment act from the status-change timeline."""

    def _ev(self, when, status, previous=None, project=mr.PROJECT_NUMBER):
        return {"createdAt": when, "status": status, "previousStatus": previous,
                "project": {"number": project}}

    def test_commitment_is_first_entry_INTO_ready_whatever_preceded_it(self):
        """NOT `previousStatus == "Backlog"`.

        An item can be committed straight from intake (No Status -> Ready) - verified live
        on Issue #1162. Keying on the previous status would silently MISS that commitment,
        under-counting commitments and inflating the never-committed population.
        """
        events = [self._ev("2026-07-14T16:05:00Z", "Ready", previous=None)]
        assert mr.first_ready_at(events) == "2026-07-14T16:05:00Z"

    def test_earliest_ready_wins_when_the_item_bounces(self):
        """Ready -> In progress -> Ready: the commitment act is the FIRST crossing."""
        events = [
            self._ev("2026-07-01T00:00:00Z", "Ready", "Backlog"),
            self._ev("2026-07-02T00:00:00Z", "In progress", "Ready"),
            self._ev("2026-07-03T00:00:00Z", "Ready", "In progress"),   # re-entry, not new
        ]
        assert mr.first_ready_at(events) == "2026-07-01T00:00:00Z"

    def test_never_ready_is_none(self):
        events = [
            self._ev("2026-07-01T00:00:00Z", "Backlog", None),
            self._ev("2026-07-02T00:00:00Z", "Done", "Backlog"),        # closed from Backlog
        ]
        assert mr.first_ready_at(events) is None

    def test_another_projects_ready_does_not_count(self):
        events = [self._ev("2026-07-01T00:00:00Z", "Ready", "Backlog", project=999)]
        assert mr.first_ready_at(events) is None, "a Ready on another board is not our commitment"

    def test_empty_timeline_is_none(self):
        assert mr.first_ready_at([]) is None


class TestGeneratedHtmlHasNoEmDash:
    """The report is machine-generated, so the PreToolUse no-em-dash guard cannot see
    it - that guard scans Claude's edits, not a script's output. Left unchecked, every
    report the script ever produces violates the house rule. Issue #1180.

    **The rule is about AUTHORED strings, not the rendered artifact.** The report also
    embeds *data* - Issue titles straight from GitHub, some of which genuinely contain
    an em-dash (e.g. "fix(guards): cross-repo coverage - project-repo guards..."). Those
    must pass through **verbatim**: a report that rewrites a real title is falsifying its
    own source data, which is a far worse sin than a stray dash. So these tests assert
    over what the *script writes*, never over the generated HTML as a whole.
    """

    def test_template_emits_no_em_dash(self):
        template = mr.TEMPLATE_PATH.read_text(encoding="utf-8")
        assert EM_DASH not in template, "em-dash in the report template"
        assert EN_DASH not in template, "en-dash in the report template"

    def test_authored_output_strings_are_not_em_dashes(self):
        """Every em-dash the SCRIPT emits, across all three emitters.

        A template-only check passes while the HTML still carries them, because the
        script emits authored strings from three separate places:
          - the normalizer's Size/Pri/status placeholders (`or "..."` / `else "..."`)
          - the `_fmt` None-placeholder
          - the seeded-narrative boilerplate (comments, TBD suffixes)
        The first live HTML render caught a fourth I had missed by grepping - a Jinja
        `pct` filter - which is precisely why this asserts over the source, not a diff.
        """
        tree = ast.parse(Path(mr.__file__).read_text(encoding="utf-8"))

        # Docstrings are prose ABOUT the code, never emitted into the artifact, so they
        # are exempt (this file's own history: a line-based heuristic kept snagging them
        # and had to be replaced by this AST walk - approximate matching on source text
        # is exactly the class of bug the rest of this PR is about).
        docstrings = set()
        for node in ast.walk(tree):
            if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
                doc = ast.get_docstring(node, clean=False)
                if doc:
                    docstrings.add(doc)

        # A dash inside a `.lstrip(...)`/`.strip(...)` argument is a char-class being
        # stripped OFF INPUT - it emits nothing. Identify those by their AST CONTEXT
        # (an argument to a strip call), not by their value.
        #
        # The first cut skipped them with `text.startswith("#*->")`, which the review
        # correctly called out: that is a value-heuristic, so a future *emitted* string
        # that happened to start with those characters would be silently exempted. Using
        # an approximate match inside a test written to reject approximate matching would
        # have been a poor joke to ship.
        stripped_args = set()
        for node in ast.walk(tree):
            if (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Attribute)
                    and node.func.attr in {"strip", "lstrip", "rstrip"}):
                for arg in node.args:
                    if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                        stripped_args.add(id(arg))

        offenders = []
        for node in ast.walk(tree):
            if not (isinstance(node, ast.Constant) and isinstance(node.value, str)):
                continue
            if id(node) in stripped_args:
                continue
            text = node.value
            if text in docstrings:
                continue
            if EM_DASH in text or EN_DASH in text:
                offenders.append(f"line {node.lineno}: {text[:60]!r}")

        assert not offenders, (
            "the script emits an em/en dash into the generated report: " + str(offenders)
        )
