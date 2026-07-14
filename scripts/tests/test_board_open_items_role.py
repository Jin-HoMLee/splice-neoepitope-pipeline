# scripts/tests/test_board_open_items_role.py
#
# Issue #1153: `--role X` kept only the FIRST `role:` label, so a dual-role Issue
# was invisible in the queue of whichever role was not listed first. Measured on
# the live board 2026-07-14: 12 open Issues hidden from `--role developer`,
# including #1112 - an Issue on which Developer is the Lead/DRI.
#
# This mirrors the arc axis, which already carries the full label set (#1103) for
# exactly this reason. GitHub returns labels in an UNSTABLE order, so first-only
# matching does not merely hide the second role - it hides a *nondeterministic*
# one, which is why this cannot be worked around by "just label it in the right
# order".
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import board_open_items as b  # noqa: E402


def _item(labels):
    return {
        "content": {
            "__typename": "Issue", "number": 1112, "title": "t", "url": "u",
            "state": "OPEN", "createdAt": None, "updatedAt": None, "closedAt": None,
            "labels": {"nodes": [{"name": n} for n in labels]},
        },
        "fieldValues": {"nodes": [{"name": "Ready", "field": {"name": "Status"}}]},
    }


def _args(**kw):
    base = dict(role=None, status=None, priority=None, size=None, arc=None, arc_phase=None)
    base.update(kw)
    return argparse.Namespace(**base)


def test_normalize_collects_every_role_label():
    n = b.normalize(_item(["role:scientist", "role:developer"]))
    assert n["roles"] == ["role:scientist", "role:developer"]


def test_role_key_retains_first_label_for_existing_consumers():
    """`role` keeps its old first-only meaning so the --json shape does not break."""
    n = b.normalize(_item(["role:scientist", "role:developer"]))
    assert n["role"] == "role:scientist"


def test_single_role_item_unchanged():
    n = b.normalize(_item(["role:developer"]))
    assert n["role"] == "role:developer"
    assert n["roles"] == ["role:developer"]


def test_no_role_labels():
    n = b.normalize(_item(["arc:board-governance"]))
    assert n["role"] is None
    assert n["roles"] == []


# --- the matched pair: identical item, one variable flipped, BOTH must match ---
#
# This is the falsifier. Against the pre-fix code the `developer` half goes red
# (the item is filtered out because `role:developer` is not the first label),
# while the `scientist` half passes - so a test that only checked the first role
# would have certified the bug as correct behavior.

def test_dual_role_item_matches_under_its_first_role():
    n = b.normalize(_item(["role:scientist", "role:developer"]))
    assert b.matches_filter(n, _args(role="scientist")) is True


def test_dual_role_item_matches_under_its_second_role():
    n = b.normalize(_item(["role:scientist", "role:developer"]))
    assert b.matches_filter(n, _args(role="developer")) is True


def test_role_filter_still_excludes_an_unrelated_role():
    """The guard must still be able to say no, or it is not a filter."""
    n = b.normalize(_item(["role:scientist", "role:developer"]))
    assert b.matches_filter(n, _args(role="pm")) is False


def test_role_filter_accepts_the_prefixed_form():
    n = b.normalize(_item(["role:scientist", "role:developer"]))
    assert b.matches_filter(n, _args(role="role:developer")) is True


# --- the Role column must SHOW that a row is dual-role -----------------------
#
# Mirrors the arc axis's own rendering pair (test_board_open_items_arc.py:
# test_arc_column_marks_a_multi_arc_parent / ..._unmarked_for_a_single_arc_item).
# Without the marker, a row matched via its SECOND role displays only its first
# and reads as though it did not match the filter that returned it.

def test_role_column_marks_a_dual_role_item():
    n = b.normalize(_item(["role:scientist", "role:developer"]))
    out = b.format_table([n])
    assert "+1" in out


def test_role_column_unmarked_for_a_single_role_item():
    n = b.normalize(_item(["role:developer"]))
    out = b.format_table([n])
    assert "+1" not in out
