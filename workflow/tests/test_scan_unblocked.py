"""Tests for scripts/pm/scan_unblocked.py - the freshly-unblocked sweep (Issue #745).

The classification layer is pure and carries the bulk of the coverage; the one
gh-backed helper (pagination) is monkeypatched so nothing touches the network.

The load-bearing tests are the **matched pairs**: identical issue payloads with
exactly one variable flipped (blocker state / clear age / board Status), asserting
OPPOSITE outcomes. A sweep whose tests only ever assert "it fires" cannot fail, and
a check that cannot fail is the thing this repo keeps re-learning not to ship
(`shared/feedback_a_check_must_be_able_to_fail.md`).
"""
import datetime as dt
import sys
from pathlib import Path

import pytest

_PM_DIR = Path(__file__).resolve().parents[2] / "scripts" / "pm"
sys.path.insert(0, str(_PM_DIR))

import scan_unblocked as su  # noqa: E402


NOW = dt.datetime(2026, 7, 14, 12, 0, 0, tzinfo=dt.timezone.utc)


def make_issue(*, number=594, blockers=(), status="Backlog", roles=("role:pm",), title="deck refresh"):
    """Build one GraphQL-shaped issue node.

    `blockers` is a list of (number, state, closed_days_ago | None).
    `status` of None means "not on the board / no Status".
    """
    blocker_nodes = []
    for bnum, bstate, days_ago in blockers:
        closed_at = None
        if days_ago is not None:
            closed_at = (NOW - dt.timedelta(days=days_ago)).isoformat().replace("+00:00", "Z")
        blocker_nodes.append(
            {"number": bnum, "title": f"blocker {bnum}", "state": bstate, "closedAt": closed_at}
        )

    project_items = {"nodes": []}
    if status is not None:
        project_items = {
            "nodes": [{
                "project": {"number": su.PROJECT_NUMBER},
                "fieldValues": {"nodes": [
                    {"field": {"name": "Status"}, "name": status},
                    {"field": {"name": "Priority"}, "name": "P2"},
                ]},
            }]
        }

    return {
        "number": number,
        "title": title,
        "url": f"https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/{number}",
        "labels": {"nodes": [{"name": r} for r in roles] + [{"name": "arc:board-governance"}]},
        "blockedBy": {"nodes": blocker_nodes},
        "projectItems": project_items,
    }


def classify(issue, since_days=su.DEFAULT_SINCE_DAYS):
    return su.classify(issue, now=NOW, since_days=since_days)


# ---------------------------------------------------------------------------
# The establishing case: Issue #594's exact shape must fire.
# ---------------------------------------------------------------------------

def test_594_shape_fires():
    """Both prerequisites closed 3 days ago; dependent still parked in Backlog."""
    issue = make_issue(
        number=594,
        blockers=[(211, "CLOSED", 3), (212, "CLOSED", 3)],
        status="Backlog",
    )
    finding = classify(issue)
    assert finding is not None
    assert finding["number"] == 594
    assert finding["status"] == "Backlog"
    assert [b["number"] for b in finding["cleared_by"]] == [211, 212]
    assert finding["age_days"] == pytest.approx(3.0, abs=0.1)


# ---------------------------------------------------------------------------
# MATCHED PAIR 1 - blocker state. Same issue; flip one blocker to OPEN.
# ---------------------------------------------------------------------------

def test_matched_pair_blocker_state():
    fires = make_issue(blockers=[(211, "CLOSED", 3), (212, "CLOSED", 3)])
    silent = make_issue(blockers=[(211, "CLOSED", 3), (212, "OPEN", None)])

    assert classify(fires) is not None, "all blockers closed -> must fire"
    assert classify(silent) is None, "one blocker still OPEN -> genuinely still blocked"


# ---------------------------------------------------------------------------
# MATCHED PAIR 2 - clear age. Same issue; move the clear outside the window.
# This is the pair that protects late-commitment Kanban: an option resting in
# Backlog with an OLD clear is NOT drift and must not be nagged about forever.
# ---------------------------------------------------------------------------

def test_matched_pair_clear_age():
    fresh = make_issue(blockers=[(211, "CLOSED", 3)])
    stale = make_issue(blockers=[(211, "CLOSED", 40)])

    assert classify(fresh) is not None, "cleared 3d ago -> a state change worth surfacing"
    assert classify(stale) is None, "cleared 40d ago -> a resting option, not a state change"


def test_window_boundary_is_inclusive():
    inside = make_issue(blockers=[(211, "CLOSED", 14)])
    outside = make_issue(blockers=[(211, "CLOSED", 15)])

    assert classify(inside, since_days=14) is not None
    assert classify(outside, since_days=14) is None


def test_since_days_widens_the_window():
    issue = make_issue(blockers=[(211, "CLOSED", 40)])

    assert classify(issue, since_days=14) is None
    assert classify(issue, since_days=60) is not None, "--since-days must actually widen"


# ---------------------------------------------------------------------------
# MATCHED PAIR 3 - board Status. Same issue; someone already pulled it.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("status", ["Backlog", "No Status", None])
def test_uncommitted_statuses_fire(status):
    assert classify(make_issue(blockers=[(211, "CLOSED", 2)], status=status)) is not None


@pytest.mark.parametrize("status", ["Ready", "In progress", "Ready for review", "In review", "Done", "Epic"])
def test_committed_statuses_are_silent(status):
    """Somebody pulled it after the clear - the system worked. Say nothing.

    `Epic` is here on purpose: a parent is parked in Epic by design (Pattern A2)
    and is never itself pulled, so a parent whose blocker clears is not a missed
    commitment.
    """
    assert classify(make_issue(blockers=[(211, "CLOSED", 2)], status=status)) is None


# ---------------------------------------------------------------------------
# The never-blocked case, and the undateable-clear case.
# ---------------------------------------------------------------------------

def test_no_blockers_is_silent():
    """Most of the board. No wired blocker ever -> nothing cleared -> no state change."""
    assert classify(make_issue(blockers=[])) is None


def test_closed_blocker_without_timestamp_is_silent():
    """CLOSED but no closedAt: we cannot date the change, so we cannot honour the
    window. Stay silent rather than invent a timestamp - an unjustifiable finding is
    worse than a miss."""
    assert classify(make_issue(blockers=[(211, "CLOSED", None)])) is None


# ---------------------------------------------------------------------------
# board_status extraction
# ---------------------------------------------------------------------------

def test_board_status_ignores_other_projects():
    nodes = [{
        "project": {"number": 999},
        "fieldValues": {"nodes": [{"field": {"name": "Status"}, "name": "Ready"}]},
    }]
    assert su.board_status(nodes) is None, "a Status on another project must not count"


def test_board_status_returns_none_when_unset():
    nodes = [{"project": {"number": su.PROJECT_NUMBER}, "fieldValues": {"nodes": []}}]
    assert su.board_status(nodes) is None


def test_board_status_reads_our_board():
    nodes = [{
        "project": {"number": su.PROJECT_NUMBER},
        "fieldValues": {"nodes": [
            {"field": {"name": "Priority"}, "name": "P1"},
            {"field": {"name": "Status"}, "name": "Ready"},
        ]},
    }]
    assert su.board_status(nodes) == "Ready"


# ---------------------------------------------------------------------------
# Roles + render
# ---------------------------------------------------------------------------

def test_multi_role_labels_all_surface():
    issue = make_issue(blockers=[(211, "CLOSED", 1)], roles=("role:developer", "role:pm"))
    assert classify(issue)["roles"] == ["role:developer", "role:pm"]


def test_render_empty_states_the_window():
    out = su.render([], since_days=14)
    assert "no findings" in out
    assert "14d" in out


def test_render_names_issue_blocker_and_age():
    finding = classify(make_issue(number=594, blockers=[(211, "CLOSED", 3)]))
    out = su.render([finding], since_days=14)
    assert "#594" in out
    assert "cleared by #211" in out
    assert "Backlog" in out
    assert "role:pm" in out
    assert "never an auto-commit" in out, "advisory framing must survive into the output"


# ---------------------------------------------------------------------------
# Pagination - `first: 100` is a cap, not a filter.
# ---------------------------------------------------------------------------

def test_fetch_paginates(monkeypatch):
    pages = [
        {"data": {"repository": {"issues": {
            "pageInfo": {"hasNextPage": True, "endCursor": "CUR1"},
            "nodes": [make_issue(number=1)],
        }}}},
        {"data": {"repository": {"issues": {
            "pageInfo": {"hasNextPage": False, "endCursor": None},
            "nodes": [make_issue(number=2)],
        }}}},
    ]
    calls = []

    def fake_gh(*args):
        calls.append(args)
        return pages[len(calls) - 1]

    monkeypatch.setattr(su, "gh", fake_gh)

    nodes = su.fetch_open_issues()
    assert [n["number"] for n in nodes] == [1, 2], "second page must be fetched, not dropped"
    assert any("after=CUR1" in a for a in calls[1]), "cursor must be threaded into page 2"


def test_query_uses_true_ceilings_not_samples():
    """`first:` is a cap, not a filter, and GraphQL does not guarantee node order.

    50 is GitHub's documented per-direction blocker ceiling, so it is a TRUE BOUND.
    At a *sampling* cap, an issue could return N all-CLOSED blockers while a still-OPEN
    one sorted past the cap - and classify() would fire a false 'freshly-unblocked'
    finding on an issue that is genuinely still blocked. That is the worst failure this
    sweep can have (it feeds 'never commit a blocked Issue to Ready' the exact input
    that rule exists to prevent), so the bound is pinned here rather than left to review.
    """
    assert "blockedBy(first: 50)" in su._QUERY


# ---------------------------------------------------------------------------
# main() exit-code contract. The module docstring PROMISES fail-open-and-loud;
# a promise with no test is a promise that can only be believed, not checked.
# ---------------------------------------------------------------------------

def test_main_fails_open_and_loud_on_gh_error(monkeypatch, capsys):
    """A sweep that cannot reach GitHub must NOT masquerade as a clean board."""
    def boom(**kwargs):
        raise GhErrorStub("api is down")

    monkeypatch.setattr(su, "sweep", boom)
    monkeypatch.setattr(su, "GhError", GhErrorStub)
    monkeypatch.setattr(sys, "argv", ["scan_unblocked.py"])

    rc = su.main()
    err = capsys.readouterr().err

    assert rc == 2, "a failed sweep must exit 2, never 0 (0 would read as 'no findings')"
    assert "FAILED" in err, "the failure must be loud on stderr, not silent"


class GhErrorStub(Exception):
    """Stand-in for gh_client.GhError so the test needs no network shape."""


def test_main_check_flag_exits_1_on_finding(monkeypatch, capsys):
    monkeypatch.setattr(su, "sweep", lambda **kw: [
        classify(make_issue(number=594, blockers=[(211, "CLOSED", 2)]))
    ])
    monkeypatch.setattr(sys, "argv", ["scan_unblocked.py", "--check"])
    assert su.main() == 1
    assert "#594" in capsys.readouterr().out


def test_main_is_advisory_without_check_flag(monkeypatch, capsys):
    """Matched pair to the above: the SAME finding, without --check, must exit 0.

    Advisory is the house style - the sweep surfaces, it never blocks.
    """
    monkeypatch.setattr(su, "sweep", lambda **kw: [
        classify(make_issue(number=594, blockers=[(211, "CLOSED", 2)]))
    ])
    monkeypatch.setattr(sys, "argv", ["scan_unblocked.py"])
    assert su.main() == 0, "advisory: a finding must not fail the beat unless --check"
    assert "#594" in capsys.readouterr().out


def test_main_clean_board_exits_0(monkeypatch, capsys):
    monkeypatch.setattr(su, "sweep", lambda **kw: [])
    monkeypatch.setattr(sys, "argv", ["scan_unblocked.py", "--check"])
    assert su.main() == 0
    assert "no findings" in capsys.readouterr().out


def test_sweep_sorts_freshest_first(monkeypatch):
    monkeypatch.setattr(su, "fetch_open_issues", lambda: [
        make_issue(number=10, blockers=[(1, "CLOSED", 9)]),
        make_issue(number=20, blockers=[(2, "CLOSED", 1)]),
        make_issue(number=30, blockers=[(3, "OPEN", None)]),   # still blocked -> dropped
    ])
    findings = su.sweep(since_days=14, now=NOW)
    assert [f["number"] for f in findings] == [20, 10]
