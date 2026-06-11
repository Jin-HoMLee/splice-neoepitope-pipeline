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


def test_normalize_multi_arc_picks_first_and_warns(capsys):
    n = b.normalize(_item(["arc:scoring-tcr-pmhc", "arc:board-governance"]))
    assert n["arc"] == "arc:scoring-tcr-pmhc"
    assert "multiple arc labels" in capsys.readouterr().err


def test_filter_by_arc_phase():
    it = b.normalize(_item(["arc:scoring-tcr-pmhc", "arc-phase:active"]))
    assert b.matches_filter(it, _args(arc_phase="active"))
    assert not b.matches_filter(it, _args(arc_phase="later"))


def test_filter_by_arc_slug_accepts_short_and_full_form():
    it = b.normalize(_item(["arc:scoring-tcr-pmhc", "arc-phase:active"]))
    assert b.matches_filter(it, _args(arc="scoring-tcr-pmhc"))
    assert b.matches_filter(it, _args(arc="arc:scoring-tcr-pmhc"))
    assert not b.matches_filter(it, _args(arc="cloud-reproducibility"))


def test_format_table_omits_arc_columns_by_default():
    it = b.normalize(_item(["role:pm", "arc:board-governance", "arc-phase:active"]))
    out = b.format_table([it])
    assert "Arc" not in out.splitlines()[0]
    assert "board-governance" not in out


def test_format_table_arc_columns_render_slug_and_phase():
    it = b.normalize(_item(["role:pm", "arc:board-governance", "arc-phase:active"]))
    out = b.format_table([it], arc_columns=True)
    header = out.splitlines()[0]
    assert "Arc" in header and "Ph" in header
    # slug is shown without the "arc:" prefix
    assert "board-governance" in out
    assert "arc:board-governance" not in out
    assert "active" in out


def test_format_table_arc_columns_handle_missing_arc():
    it = b.normalize(_item(["role:pm"]))  # no arc / phase labels
    out = b.format_table([it], arc_columns=True)
    assert "Arc" in out.splitlines()[0]
    assert "—" in out  # placeholder, no crash
