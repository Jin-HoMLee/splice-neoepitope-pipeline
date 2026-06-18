"""Unit tests for scripts/pm/milestone_report.py — the metrics pure-functions.

Only the metrics layer is covered here (per the design spec §9): the metric
functions are pure (no I/O, no GraphQL, no jinja2/markdown), so they unit-test
cleanly with small synthetic fixtures. The data/aggregation/render layers are
exercised by the pm-i6 pilot dry-run, not by CI unit tests.

Importing this module must NOT pull in jinja2/markdown — those are lazy-imported
inside the render/aggregation functions so this test runs in the bare
``ci-tools-pytest`` env (pytest + pyyaml only).
"""

import sys
from pathlib import Path

import pytest

# Make the script importable as a module (mirrors test_recheck_milestone.py)
SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts" / "pm"
sys.path.insert(0, str(SCRIPT_DIR))

import milestone_report as mr


# --- fixtures ---------------------------------------------------------------

def _issue(number, state, created_at, closed_at, roles):
    """Minimal normalized-issue dict carrying only what the metrics read."""
    return {
        "number": number,
        "state": state,
        "created_at": created_at,
        "closed_at": closed_at,
        "roles": roles,
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

class TestCycleTime:
    def test_closed_issue_days(self):
        assert mr.cycle_time_days(SAMPLE[0]) == 2.0
        assert mr.cycle_time_days(SAMPLE[1]) == 4.0

    def test_open_issue_is_none(self):
        assert mr.cycle_time_days(SAMPLE[3]) is None

    def test_missing_timestamp_is_none(self):
        assert mr.cycle_time_days(_issue(9, "CLOSED", None, "2026-06-03T00:00:00Z", [])) is None

    def test_cycle_times_only_closed(self):
        assert sorted(mr.cycle_times(SAMPLE)) == [2.0, 4.0, 6.0]

    def test_avg_cycle_time(self):
        assert mr.avg_cycle_time(SAMPLE) == pytest.approx(4.0)  # (2+4+6)/3

    def test_median_cycle_time(self):
        assert mr.median_cycle_time(SAMPLE) == pytest.approx(4.0)

    def test_avg_empty_is_none(self):
        assert mr.avg_cycle_time([]) is None
        assert mr.median_cycle_time([]) is None


# --- per-role counts --------------------------------------------------------

class TestPerRoleCounts:
    def test_counts_closed_only_multi_role(self):
        # closed: #1 sci, #2 dev, #3 sci+dev ; open #4 pm excluded
        assert mr.per_role_counts(SAMPLE) == {"role:scientist": 2, "role:developer": 2}

    def test_empty(self):
        assert mr.per_role_counts([]) == {}


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
        )
        assert m["n_total"] == 4
        assert m["n_closed"] == 3
        assert m["n_carried_forward"] == 1
        assert m["duration_days"] == 7.0
        assert m["throughput_per_week"] == pytest.approx(3.0)
        assert m["avg_cycle_time_days"] == pytest.approx(4.0)
        assert m["median_cycle_time_days"] == pytest.approx(4.0)
        assert m["per_role_counts"] == {"role:scientist": 2, "role:developer": 2}

    def test_empty_milestone_degrades(self):
        m = mr.compute_metrics([], {"created_at": None, "closed_at": None})
        assert m["n_total"] == m["n_closed"] == m["n_carried_forward"] == 0
        assert m["duration_days"] is None
        assert m["throughput_per_week"] is None
        assert m["avg_cycle_time_days"] is None
        assert m["per_role_counts"] == {}


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
