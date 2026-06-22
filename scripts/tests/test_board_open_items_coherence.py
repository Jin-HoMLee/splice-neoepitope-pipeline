# scripts/tests/test_board_open_items_coherence.py
#
# Issue #765 — coherence guard: a committed item (Ready/In-progress/… OR
# milestoned) must not carry arc-phase:later (committed = "pull now", later =
# "parked, not now" — a contradiction). Fail-open when arc-phase is absent.
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import board_open_items as b  # noqa: E402


def _norm(status="Backlog", milestone=None, arc_phase=None):
    """Minimal normalized-item shape the coherence check reads."""
    return {"number": 1, "status": status, "milestone": milestone, "arc_phase": arc_phase}


# --- committed-by-status + parked arc → incoherent ---
def test_ready_plus_arc_phase_later_is_incoherent():
    assert b.is_arc_phase_incoherent(_norm(status="Ready", arc_phase="later")) is True


def test_in_progress_plus_later_is_incoherent():
    assert b.is_arc_phase_incoherent(_norm(status="In progress", arc_phase="later")) is True


def test_in_review_plus_later_is_incoherent():
    assert b.is_arc_phase_incoherent(_norm(status="In review", arc_phase="later")) is True


def test_ready_for_review_plus_later_is_incoherent():
    assert b.is_arc_phase_incoherent(_norm(status="Ready for review", arc_phase="later")) is True


# --- off-ladder Epic (un-milestoned parents) is never committed → never flagged ---
# The load-bearing carve-out: epics sit in Epic Status and go un-milestoned by
# design (#690/#776), so they are not "committed" even when parked at later.
# This locks it in — if "Epic" ever leaks into COMMITTED_STATUSES, this fails.
def test_epic_status_unmilestoned_later_is_coherent():
    assert b.is_arc_phase_incoherent(_norm(status="Epic", arc_phase="later")) is False


# --- committed-by-milestone (even while Backlog) + parked arc → incoherent ---
def test_milestoned_backlog_plus_later_is_incoherent():
    assert (
        b.is_arc_phase_incoherent(
            _norm(status="Backlog", milestone="i5 - S3 - Data Preparation", arc_phase="later")
        )
        is True
    )


# --- uncommitted + parked arc → coherent (no contradiction) ---
def test_backlog_unmilestoned_plus_later_is_coherent():
    assert b.is_arc_phase_incoherent(_norm(status="Backlog", arc_phase="later")) is False


def test_no_status_plus_later_is_coherent():
    assert b.is_arc_phase_incoherent(_norm(status="No Status", arc_phase="later")) is False


# --- committed + active/next arc → coherent ---
def test_ready_plus_active_is_coherent():
    assert b.is_arc_phase_incoherent(_norm(status="Ready", arc_phase="active")) is False


def test_ready_plus_next_is_coherent():
    assert b.is_arc_phase_incoherent(_norm(status="Ready", arc_phase="next")) is False


# --- fail-open: no arc-phase label → no opinion ---
def test_no_arc_phase_is_coherent():
    assert b.is_arc_phase_incoherent(_norm(status="Ready", arc_phase=None)) is False
    assert b.is_arc_phase_incoherent(_norm(status="Ready", milestone="i5", arc_phase=None)) is False


# --- normalize now captures the milestone title (the "committed-by-milestone" half) ---
def _item(labels, status="Backlog", milestone_title=None):
    content = {
        "__typename": "Issue", "number": 1, "title": "t", "url": "u", "state": "OPEN",
        "createdAt": None, "updatedAt": None, "closedAt": None,
        "labels": {"nodes": [{"name": n} for n in labels]},
    }
    if milestone_title is not None:
        content["milestone"] = {"title": milestone_title}
    return {
        "content": content,
        "fieldValues": {"nodes": [{"name": status, "field": {"name": "Status"}}]},
    }


def test_normalize_captures_milestone_title():
    n = b.normalize(_item([], milestone_title="i5 - S3 - Data Preparation"))
    assert n["milestone"] == "i5 - S3 - Data Preparation"


def test_normalize_milestone_none_when_absent():
    n = b.normalize(_item([]))
    assert n["milestone"] is None
