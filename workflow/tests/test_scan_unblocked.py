"""Tests for scripts/pm/scan_unblocked.py - the unblocked-but-undecided sweep (Issue #745).

The classification layer is pure and carries the bulk of the coverage; the gh-backed
helpers are monkeypatched so nothing touches the network.

The load-bearing tests are the **matched pairs**: identical payloads with exactly one
variable flipped (blocker state / board Status / ack label), asserting OPPOSITE
outcomes. A sweep whose tests only ever assert "it fires" cannot fail, and a check that
cannot fail is the thing this repo keeps re-learning not to ship
(`shared/feedback_a_check_must_be_able_to_fail.md`).

The single most important test here is `test_findings_never_expire`. The first cut of
this script used a time window, which meant a sweep that did not RUN inside the window
dropped the finding FOREVER. That is silent data loss, and it is the exact failure the
level-triggered redesign exists to remove - so it gets a test that would go red if
anyone ever reintroduces a clock dependency.
"""
import datetime as dt
import sys
from pathlib import Path

import pytest

_PM_DIR = Path(__file__).resolve().parents[2] / "scripts" / "pm"
sys.path.insert(0, str(_PM_DIR))

import scan_unblocked as su  # noqa: E402


NOW = dt.datetime(2026, 7, 14, 12, 0, 0, tzinfo=dt.timezone.utc)


class GhErrorStub(Exception):
    """Stand-in for gh_client.GhError so the tests need no network shape."""


def make_issue(*, number=594, blockers=(), status="Backlog", roles=("role:pm",),
               labels=(), title="deck refresh"):
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

    label_names = list(roles) + list(labels) + ["arc:board-governance"]
    return {
        "number": number,
        "title": title,
        "url": f"https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/{number}",
        "labels": {"nodes": [{"name": n} for n in label_names]},
        "blockedBy": {"nodes": blocker_nodes},
        "projectItems": project_items,
    }


def classify(issue, now=NOW):
    return su.classify(issue, now=now)


# ---------------------------------------------------------------------------
# The establishing case: Issue #594's exact shape must fire.
# ---------------------------------------------------------------------------

def test_594_shape_fires():
    """Both prerequisites closed; dependent still parked in Backlog, undecided."""
    finding = classify(make_issue(
        number=594, blockers=[(211, "CLOSED", 3), (212, "CLOSED", 3)], status="Backlog"))
    assert finding is not None
    assert finding["number"] == 594
    assert finding["status"] == "Backlog"
    assert [b["number"] for b in finding["cleared_by"]] == [211, 212]
    assert finding["waiting_days"] == pytest.approx(3.0, abs=0.1)


# ---------------------------------------------------------------------------
# THE CENTRAL PROPERTY: level-triggered. Nothing expires, ever.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("days_ago", [0, 1, 14, 15, 90, 365, 3650])
def test_findings_never_expire(days_ago):
    """A finding must NOT depend on how long ago the blocker closed, or on when the
    sweep happens to run.

    The first cut expired findings after a 14d window, so a sweep that did not RUN in
    time lost the finding permanently - silent data loss, strictly worse than the
    repeat mention it was avoiding. This test would go red if anyone reintroduced a
    clock dependency into the decision.
    """
    issue = make_issue(blockers=[(211, "CLOSED", days_ago)])
    assert classify(issue) is not None, (
        f"blocker cleared {days_ago}d ago is STILL undecided - a level-triggered loop "
        f"re-reads the world and must report it every run until somebody decides"
    )


def test_classification_is_independent_of_now():
    """Same world, two wildly different clocks -> identical decision."""
    issue = make_issue(blockers=[(211, "CLOSED", 400)])
    a = classify(issue, now=NOW)
    b = classify(issue, now=NOW + dt.timedelta(days=5000))
    assert a is not None and b is not None, "both must fire - the clock decides nothing"
    assert a["number"] == b["number"]


# ---------------------------------------------------------------------------
# MATCHED PAIR 1 - blocker state.
# ---------------------------------------------------------------------------

def test_matched_pair_blocker_state():
    fires = make_issue(blockers=[(211, "CLOSED", 3), (212, "CLOSED", 3)])
    silent = make_issue(blockers=[(211, "CLOSED", 3), (212, "OPEN", None)])

    assert classify(fires) is not None, "all blockers closed -> available -> must fire"
    assert classify(silent) is None, "one blocker still OPEN -> genuinely still blocked"


# ---------------------------------------------------------------------------
# MATCHED PAIR 2 - the ack label. This is what replaces the time window: the
# decline is RECORDED rather than FORGOTTEN.
# ---------------------------------------------------------------------------

def test_matched_pair_ack_label():
    undecided = make_issue(blockers=[(211, "CLOSED", 3)])
    acked = make_issue(blockers=[(211, "CLOSED", 3)], labels=[su.ACK_LABEL])

    assert classify(undecided) is not None, "nobody has decided -> keep surfacing"
    assert classify(acked) is None, "deliberately declined -> recorded -> go quiet"


def test_ack_is_revocable():
    """Removing the label brings the finding back - the decision is revoked, not erased.

    A time window could never do this: once expired, the finding was gone with no way
    to ask for it again.
    """
    acked = make_issue(blockers=[(211, "CLOSED", 999)], labels=[su.ACK_LABEL])
    assert classify(acked) is None

    acked["labels"]["nodes"] = [n for n in acked["labels"]["nodes"]
                               if n["name"] != su.ACK_LABEL]
    assert classify(acked) is not None, "label removed -> back in the sweep"


# ---------------------------------------------------------------------------
# MATCHED PAIR 3 - board Status. Committing it is the OTHER way to clear a finding.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("status", ["Backlog", "No Status", None])
def test_uncommitted_statuses_fire(status):
    assert classify(make_issue(blockers=[(211, "CLOSED", 2)], status=status)) is not None


@pytest.mark.parametrize("status", ["Ready", "In progress", "Ready for review", "In review", "Done", "Epic"])
def test_committed_statuses_are_silent(status):
    """Somebody pulled it after the clear - the system worked. Say nothing.

    `Epic` is here on purpose: a parent is parked in Epic by design (Pattern A2) and is
    never itself pulled, so a parent whose blocker clears is not a missed commitment.
    """
    assert classify(make_issue(blockers=[(211, "CLOSED", 2)], status=status)) is None


# ---------------------------------------------------------------------------
# The never-blocked case, and the undateable clear.
# ---------------------------------------------------------------------------

def test_no_blockers_is_silent():
    """Most of the board. No wired blocker -> nothing to clear -> not our business."""
    assert classify(make_issue(blockers=[])) is None


def test_closed_blocker_without_timestamp_still_fires():
    """A missing closedAt must NOT suppress the finding.

    Under the old window design it did (no date -> cannot honour the window -> silent).
    That was only ever a consequence of the clock dependency. The finding does not
    depend on WHEN it cleared, only THAT it has - so an undateable clear is reported,
    with the date shown as unknown.
    """
    finding = classify(make_issue(blockers=[(211, "CLOSED", None)]))
    assert finding is not None
    assert finding["waiting_days"] is None
    assert finding["cleared_at"] is None


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


def test_ack_label_is_not_reported_as_a_role():
    issue = make_issue(blockers=[(211, "CLOSED", 1)], labels=["infrastructure"])
    assert classify(issue)["roles"] == ["role:pm"]


def test_render_empty_states_the_desired_state():
    out = su.render([])
    assert "no findings" in out
    assert "committed or acknowledged" in out


def test_render_names_both_terminal_actions():
    """The report must tell you how to make a finding STOP - that is what separates a
    to-do list from a nag."""
    finding = classify(make_issue(number=594, blockers=[(211, "CLOSED", 3)]))
    out = su.render([finding])
    assert "#594" in out
    assert "cleared by #211" in out
    assert "COMMIT it" in out
    assert "--ack" in out
    assert "never an auto-commit" in out, "advisory framing must survive into the output"


# ---------------------------------------------------------------------------
# Query shape - `first: N` is a cap, not a filter.
# ---------------------------------------------------------------------------

def test_query_uses_true_ceilings_not_samples():
    """50 is GitHub's documented per-direction blocker ceiling, so it is a TRUE BOUND.

    At a *sampling* cap an issue could return N all-CLOSED blockers while a still-OPEN
    one sorted past the cap - and classify() would fire a false "unblocked" finding on
    an issue that is genuinely still blocked. That is the worst failure this sweep can
    have (it feeds "never commit a blocked Issue to Ready" the exact input that rule
    exists to prevent), so the bound is pinned here rather than left to review.
    """
    assert "blockedBy(first: 50)" in su._QUERY


def test_query_requests_labels_so_ack_is_visible():
    """The ack gate is only as good as the payload it reads."""
    assert "labels(first:" in su._QUERY


# ---------------------------------------------------------------------------
# Pagination
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

    def fake_gh(*args, **kwargs):
        calls.append(args)
        return pages[len(calls) - 1]

    monkeypatch.setattr(su, "gh", fake_gh)

    nodes = su.fetch_open_issues()
    assert [n["number"] for n in nodes] == [1, 2], "second page must be fetched, not dropped"
    assert any("after=CUR1" in a for a in calls[1]), "cursor must be threaded into page 2"


def test_sweep_sorts_longest_undecided_first(monkeypatch):
    monkeypatch.setattr(su, "fetch_open_issues", lambda: [
        make_issue(number=10, blockers=[(1, "CLOSED", 2)]),
        make_issue(number=20, blockers=[(2, "CLOSED", 30)]),
        make_issue(number=30, blockers=[(3, "OPEN", None)]),      # still blocked -> dropped
        make_issue(number=40, blockers=[(4, "CLOSED", 5)], labels=[su.ACK_LABEL]),  # acked
    ])
    findings = su.sweep(now=NOW)
    assert [f["number"] for f in findings] == [20, 10], "longest-undecided first"


# ---------------------------------------------------------------------------
# acknowledge() - the write path
# ---------------------------------------------------------------------------

def test_acknowledge_creates_the_label_before_using_it(monkeypatch):
    """The write path must be SELF-SUFFICIENT.

    `gh issue edit --add-label` does NOT create a missing label - it fails HTTP 422, which
    gh_client treats as deterministic and does not retry. So without an upfront create, the
    FIRST --ack on a fresh clone dies. And --ack is one of the two terminal actions the
    whole level-triggered design rests on: a design whose principle is "every finding has a
    working way to be cleared" cannot gate its clear-path on an undocumented manual step.

    The label exists in our repo today only because it was created BY HAND during testing -
    precisely the invisible prerequisite this test exists to make impossible.
    """
    calls = []
    monkeypatch.setattr(su, "gh", lambda *a, **kw: calls.append(a))

    su.acknowledge(594, "resting on purpose: superseded by the S4 rewrite")

    verbs = [c[0:2] for c in calls]
    assert verbs[0] == ("label", "create"), "label must be ensured BEFORE it is applied"
    assert "--force" in calls[0], "label create must be idempotent (upsert), not a one-shot"
    assert verbs.index(("label", "create")) < verbs.index(("issue", "edit")), \
        "ordering is the whole point: create, then add"

    assert ("issue", "comment") in verbs
    assert ("issue", "edit") in verbs
    edit = next(c for c in calls if c[0:2] == ("issue", "edit"))
    assert su.ACK_LABEL in edit

    body = next(a for c in calls if c[0:2] == ("issue", "comment")
                for a in c if "resting on purpose" in str(a))
    assert "remove the label and it returns" in body.lower(), (
        "the ack must advertise its own revocability - a decision, not a dismissal"
    )


def test_ack_requires_a_reason(monkeypatch, capsys):
    """An unexplained decline is indistinguishable from a dropped finding."""
    monkeypatch.setattr(sys, "argv", ["scan_unblocked.py", "--ack", "594"])
    with pytest.raises(SystemExit):
        su.main()
    assert "--reason" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# main() exit-code contract.
# ---------------------------------------------------------------------------

def test_main_fails_open_and_loud_on_gh_error(monkeypatch, capsys):
    """A sweep that cannot reach GitHub must NOT masquerade as a clean board."""
    def boom(**kwargs):
        raise GhErrorStub("api is down")

    monkeypatch.setattr(su, "sweep", boom)
    monkeypatch.setattr(su, "GhError", GhErrorStub)
    monkeypatch.setattr(sys, "argv", ["scan_unblocked.py"])

    rc = su.main()
    assert rc == 2, "a failed sweep must exit 2, never 0 (0 would read as 'no findings')"
    assert "FAILED" in capsys.readouterr().err, "the failure must be loud on stderr"


def test_main_check_flag_exits_1_on_finding(monkeypatch, capsys):
    monkeypatch.setattr(su, "sweep", lambda **kw: [
        classify(make_issue(number=594, blockers=[(211, "CLOSED", 2)]))
    ])
    monkeypatch.setattr(sys, "argv", ["scan_unblocked.py", "--check"])
    assert su.main() == 1
    assert "#594" in capsys.readouterr().out


def test_main_is_advisory_without_check_flag(monkeypatch, capsys):
    """Matched pair to the above: the SAME finding, without --check, must exit 0."""
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
