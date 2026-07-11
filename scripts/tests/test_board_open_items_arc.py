# scripts/tests/test_board_open_items_arc.py
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import board_open_items as b  # noqa: E402


def _item(labels, children=0):
    content = {
        "__typename": "Issue", "number": 1, "title": "t", "url": "u", "state": "OPEN",
        "createdAt": None, "updatedAt": None, "closedAt": None,
        "labels": {"nodes": [{"name": n} for n in labels]},
    }
    if children:
        content["subIssuesSummary"] = {"total": children}
    return {
        "content": content,
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


def test_normalize_multi_arc_on_a_leaf_picks_first_and_warns(capsys):
    """One arc per LEAF is still the rule - a multi-arc leaf is drift."""
    n = b.normalize(_item(["arc:scoring-tcr-pmhc", "arc:board-governance"]))
    assert n["arc"] == "arc:scoring-tcr-pmhc"
    assert "multiple arc labels" in capsys.readouterr().err


# --- #1103: multi-arc is legal at the parent/initiative tier -----------------
#
# A parent/epic may legitimately span themes (#1036 spans scoring-tcr-pmhc and
# immunogenicity-benchmark). Standard practice permits overlapping labels there
# (Cohn: "theme and epic are labels, not an implied hierarchy"). Warning on it
# pushed toward stripping a true label to satisfy a checker.


def test_normalize_multi_arc_on_a_parent_does_not_warn(capsys):
    n = b.normalize(_item(["arc:scoring-tcr-pmhc", "arc:immunogenicity-benchmark"], children=2))
    assert n["is_parent"] is True
    assert "multiple arc labels" not in capsys.readouterr().err


def test_normalize_multi_arc_parent_still_exposes_all_arcs():
    """`arc` keeps the first (single-value consumers), `arcs` carries the full set."""
    n = b.normalize(_item(["arc:scoring-tcr-pmhc", "arc:immunogenicity-benchmark"], children=2))
    assert n["arc"] == "arc:scoring-tcr-pmhc"
    assert n["arcs"] == ["arc:scoring-tcr-pmhc", "arc:immunogenicity-benchmark"]


def test_normalize_single_arc_leaf_exposes_arcs_too():
    n = b.normalize(_item(["arc:board-governance"]))
    assert n["arcs"] == ["arc:board-governance"]


def test_normalize_no_arc_has_empty_arcs():
    assert b.normalize(_item(["role:pm"]))["arcs"] == []


def test_arc_filter_matches_ANY_arc_of_a_multi_arc_parent():
    """`--arc <slug>` must find a multi-arc parent by *any* of its arcs.

    Matching only the first makes an arc census silently wrong: GitHub returns
    labels in an UNSTABLE order, so `--arc immunogenicity-benchmark` would find
    #1036 or not depending on which label happened to come first. Latent before
    #1103 (multi-arc was drift to be swept); a live wrong answer now that a
    multi-arc parent is a permanent, legitimate state.
    """
    parent = b.normalize(
        _item(["arc:scoring-tcr-pmhc", "arc:immunogenicity-benchmark"], children=2)
    )
    assert b.matches_filter(parent, _args(arc="scoring-tcr-pmhc"))           # first
    assert b.matches_filter(parent, _args(arc="immunogenicity-benchmark"))   # second
    assert b.matches_filter(parent, _args(arc="arc:immunogenicity-benchmark"))
    assert not b.matches_filter(parent, _args(arc="board-governance"))


def test_arc_filter_is_label_order_independent():
    """Reversing label order must not change which arcs match."""
    a = b.normalize(_item(["arc:scoring-tcr-pmhc", "arc:immunogenicity-benchmark"], children=2))
    z = b.normalize(_item(["arc:immunogenicity-benchmark", "arc:scoring-tcr-pmhc"], children=2))
    for slug in ("scoring-tcr-pmhc", "immunogenicity-benchmark"):
        assert b.matches_filter(a, _args(arc=slug)) == b.matches_filter(z, _args(arc=slug))


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


def test_arc_column_marks_a_multi_arc_parent():
    """A row matched by --arc via its 2nd arc must not *display* only its 1st.

    Otherwise `--arc immunogenicity-benchmark` returns #1036 while its Arc cell
    reads `scoring-tcr-pmhc` - a row that looks like it does not match the filter
    it matched. Show a `+N` marker so the cell is honest about the full set.
    """
    parent = _item(["arc:scoring-tcr-pmhc", "arc:immunogenicity-benchmark"], children=2)
    out = b.format_table([b.normalize(parent)], arc_columns=True)
    assert "+1" in out, out


def test_arc_column_unmarked_for_a_single_arc_item():
    out = b.format_table([b.normalize(_item(["arc:board-governance"]))], arc_columns=True)
    assert "+1" not in out
