# scripts/tests/test_board_open_items_arc.py
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import board_open_items as b  # noqa: E402


def _item(labels):
    return {
        "content": {
            "__typename": "Issue", "number": 1, "title": "t", "url": "u", "state": "OPEN",
            "createdAt": None, "updatedAt": None, "closedAt": None,
            "labels": {"nodes": [{"name": n} for n in labels]},
        },
        "fieldValues": {"nodes": [{"name": "Ready", "field": {"name": "Status"}}]},
    }


def _args(**kw):
    base = dict(role=None, status=None, priority=None, size=None, arc=None, arc_phase=None)
    base.update(kw)
    return argparse.Namespace(**base)


def test_normalize_derives_arc_and_phase():
    n = b.normalize(_item(["role:pm", "arc:board-governance", "arc-phase:active"]))
    assert n["arc"] == "arc:board-governance"
    assert n["arc_phase"] == "active"


def test_normalize_arc_none_when_absent():
    n = b.normalize(_item(["role:pm"]))
    assert n["arc"] is None and n["arc_phase"] is None


def test_arc_phase_label_not_mistaken_for_arc():
    n = b.normalize(_item(["arc-phase:later"]))
    assert n["arc"] is None and n["arc_phase"] == "later"


def test_filter_by_arc_phase():
    it = b.normalize(_item(["arc:scoring-tcr-pmhc", "arc-phase:active"]))
    assert b.matches_filter(it, _args(arc_phase="active"))
    assert not b.matches_filter(it, _args(arc_phase="later"))


def test_filter_by_arc_slug_accepts_short_and_full_form():
    it = b.normalize(_item(["arc:scoring-tcr-pmhc", "arc-phase:active"]))
    assert b.matches_filter(it, _args(arc="scoring-tcr-pmhc"))
    assert b.matches_filter(it, _args(arc="arc:scoring-tcr-pmhc"))
    assert not b.matches_filter(it, _args(arc="cloud-reproducibility"))
