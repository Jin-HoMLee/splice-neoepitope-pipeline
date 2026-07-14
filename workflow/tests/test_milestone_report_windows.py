"""Window/bucket boundary semantics for scripts/pm/milestone_report.py (Issue #1099).

The weekly SDR sources its data from GitHub's ``closed:A..B`` search operator, which
is **day-granular and inclusive on both endpoints**. The trend buckets, by contrast,
are ``(start, end]`` datetime intervals. When the two granularities disagree, issues
fall through the gap: they are fetched but land in no bucket at all.

The load-bearing property is **conservation**: every issue the fetch returns must be
counted in exactly one bucket. A count that silently drops issues corrupts the trend
series that CLAUDE.md designates as the input to the WIP-limit retune.

Pure functions only - no network.
"""
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

_PM_DIR = Path(__file__).resolve().parents[2] / "scripts" / "pm"
sys.path.insert(0, str(_PM_DIR))

import milestone_report as mr  # noqa: E402

UTC = timezone.utc


def _issue(number: int, closed_at: str) -> dict:
    """A delivered (COMPLETED) closed issue, in the shape _normalize_issue emits."""
    return {
        "number": number,
        "state": "CLOSED",
        "state_reason": "COMPLETED",
        "created_at": "2026-01-01T00:00:00Z",
        "closed_at": closed_at,
    }


# --- normalize_until ---------------------------------------------------------

def test_normalize_until_anchors_a_mid_day_time_to_end_of_utc_day():
    """A live run's datetime.now() is mid-day; the data it buckets is day-granular."""
    got = mr.normalize_until(datetime(2026, 7, 10, 14, 30, 0, tzinfo=UTC))
    assert got == datetime(2026, 7, 10, 23, 59, 59, 999999, tzinfo=UTC)


def test_normalize_until_is_idempotent():
    once = mr.normalize_until(datetime(2026, 7, 10, 14, 30, 0, tzinfo=UTC))
    assert mr.normalize_until(once) == once


def test_normalize_until_rejects_a_naive_datetime():
    """A naive datetime would be read as *local* time and shift every bucket edge by
    the host's UTC offset - the same class of silent grid drift this helper removes."""
    with pytest.raises(ValueError, match="timezone-aware"):
        mr.normalize_until(datetime(2026, 7, 10, 14, 30, 0))


# --- fetch_span must match what the buckets actually cover -------------------

def test_fetch_span_starts_the_day_after_the_first_bucket_opens():
    """Buckets are half-open at the start: ``(span_start, ...]``.

    ``span_start`` is the last instant of its own day, so that day is NOT covered by
    any bucket. Fetching from ``span_start.date()`` therefore pulls a whole day of
    issues that no bucket can hold - the bottom-edge leak.
    """
    until = mr.normalize_until(datetime(2026, 7, 10, 12, 0, 0, tzinfo=UTC))
    since_iso, until_iso = mr.fetch_span(until, 4)

    windows = mr.week_windows(until, 4)
    span_start = windows[0][0]

    assert since_iso == (span_start + timedelta(days=1)).date().isoformat()
    assert since_iso == "2026-06-13"        # NOT 2026-06-12
    assert until_iso == "2026-07-10"


# --- the load-bearing property: conservation ---------------------------------

@pytest.mark.parametrize("n_weeks", [1, 4, 6])
def test_every_fetched_issue_lands_in_exactly_one_bucket(n_weeks):
    """Conservation. Synthesize an issue closed at *every hour* of the fetch span and
    assert each is bucketed exactly once.

    This is the regression lock for Issue #1099: before the fix, issues closed on the
    fetch's first day (and, on a mid-day anchor, after the last bucket's edge) were
    returned by the fetch and counted by nothing.
    """
    until = mr.normalize_until(datetime(2026, 7, 10, 14, 30, 0, tzinfo=UTC))
    windows = mr.week_windows(until, n_weeks)
    since_iso, until_iso = mr.fetch_span(until, n_weeks)

    # Every hour the day-granular fetch `closed:since..until` would return.
    start = datetime.fromisoformat(since_iso).replace(tzinfo=UTC)
    end = datetime.fromisoformat(until_iso).replace(hour=23, minute=59, second=59, tzinfo=UTC)
    issues, t, n = [], start, 0
    while t <= end:
        issues.append(_issue(n, t.strftime("%Y-%m-%dT%H:%M:%SZ")))
        t += timedelta(hours=1)
        n += 1

    hits = {i["number"]: 0 for i in issues}
    for w_start, w_end in windows:
        for i in mr.closed_in_window(issues, w_start, w_end):
            hits[i["number"]] += 1

    dropped = [n for n, c in hits.items() if c == 0]
    doubled = [n for n, c in hits.items() if c > 1]
    assert not dropped, (
        f"{len(dropped)}/{len(issues)} fetched issues landed in NO bucket "
        f"(first dropped closed_at={issues[dropped[0]]['closed_at']})"
    )
    assert not doubled, f"{len(doubled)} fetched issues were double-counted"


def test_bucket_edge_instant_lands_in_exactly_one_bucket():
    """An issue closed at the exact shared edge of two buckets belongs to the earlier
    one (intervals are ``(start, end]``), never to both and never to neither."""
    until = mr.normalize_until(datetime(2026, 7, 10, 9, 0, 0, tzinfo=UTC))
    windows = mr.week_windows(until, 4)
    edge = windows[0][1]  # boundary shared by bucket 0 (as `end`) and bucket 1 (as `start`)

    issue = _issue(1, edge.strftime("%Y-%m-%dT%H:%M:%S.%fZ"))
    hits = [k for k, (s, e) in enumerate(windows) if mr.closed_in_window([issue], s, e)]
    assert hits == [0], f"edge instant {edge.isoformat()} landed in buckets {hits}, expected [0]"


def test_week_label_names_the_days_the_bucket_actually_covers():
    """``week_start..week_end`` renders as an inclusive day range, so it must be one.

    A reader cross-checking a trend row runs ``closed:<week_start>..<week_end>`` (a
    both-ends-inclusive GitHub search). If the label named the raw half-open ``(start``
    edge, that query would span an extra day and return a bigger number than the row.
    Mislabelling this way is what produced the bogus reference counts in Issue #1099.
    """
    until = mr.normalize_until(datetime(2026, 7, 10, 14, 30, 0, tzinfo=UTC))
    series = mr.weekly_series([], until, 4, marker_in_use=True)

    # Each row covers exactly 7 whole days, inclusive of both labelled endpoints.
    for w in series:
        start = datetime.fromisoformat(w["week_start"])
        end = datetime.fromisoformat(w["week_end"])
        assert (end - start).days == 6, f"{w['week_start']}..{w['week_end']} is not 7 whole days"

    # Consecutive rows must not share a day (that overlap is the double-count trap).
    for earlier, later in zip(series, series[1:]):
        assert earlier["week_end"] < later["week_start"], (
            f"weeks overlap on {earlier['week_end']}: rows would double-count it"
        )

    assert series[0]["week_start"] == "2026-06-13"   # NOT 2026-06-12
    assert series[-1]["week_end"] == "2026-07-10"


def test_weekly_series_totals_match_the_fetched_population():
    """The rendered trend series must sum to the delivered population it was given."""
    until = mr.normalize_until(datetime(2026, 7, 10, 14, 30, 0, tzinfo=UTC))
    since_iso, until_iso = mr.fetch_span(until, 4)

    start = datetime.fromisoformat(since_iso).replace(tzinfo=UTC)
    end = datetime.fromisoformat(until_iso).replace(hour=23, minute=59, second=59, tzinfo=UTC)
    issues, t, n = [], start, 0
    while t <= end:
        issues.append(_issue(n, t.strftime("%Y-%m-%dT%H:%M:%SZ")))
        t += timedelta(hours=6)
        n += 1

    series = mr.weekly_series(issues, until, 4, marker_in_use=True)
    assert sum(w["n_delivered"] for w in series) == len(mr.delivered_issues(issues))
