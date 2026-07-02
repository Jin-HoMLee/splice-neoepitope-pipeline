"""Tests for scripts/pm/dispatch_digest.py - the per-role board dispatch digest (Issue #721).

The digest rolls up, per role, the committed board work (Ready / In progress /
Ready for review / In review) plus open Team Coordination Discussions and aging
WIP - restoring the at-a-glance visibility the retired team_standup gave (#569).

All logic under test is pure (grouping / aging / rendering) with injected items,
discussions, and `now`, so the tests are deterministic and never touch the network.
"""
import sys
from datetime import datetime, timezone
from pathlib import Path

import pytest

# dispatch_digest lives in scripts/pm/; it self-adds scripts/ to sys.path so it can
# import board_open_items. Mirror both here so the import resolves under pytest.
_SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(_SCRIPTS))
sys.path.insert(0, str(_SCRIPTS / "pm"))

import dispatch_digest as dd  # noqa: E402

NOW = datetime(2026, 7, 2, 12, 0, 0, tzinfo=timezone.utc)


def _item(number, *, role="role:pm", status="Ready", title="t",
          updated="2026-07-02T00:00:00Z"):
    """A minimal normalized item shaped like board_open_items.normalize() output."""
    return {
        "number": number,
        "title": title,
        "url": f"https://github.com/o/r/issues/{number}",
        "kind": "Issue",
        "status": status,
        "role": role,
        "updated_at": updated,
        "priority": None,
        "size": None,
    }


# --------------------------------------------------------------------------- #
# group_by_role_status                                                         #
# --------------------------------------------------------------------------- #

def test_group_buckets_committed_items_by_role_then_status():
    items = [
        _item(1, role="role:pm", status="Ready"),
        _item(2, role="role:pm", status="In progress"),
        _item(3, role="role:developer", status="In review"),
    ]
    g = dd.group_by_role_status(items)
    assert g["pm"]["Ready"][0]["number"] == 1
    assert g["pm"]["In progress"][0]["number"] == 2
    assert g["developer"]["In review"][0]["number"] == 3


def test_group_excludes_uncommitted_statuses():
    """Backlog / No Status / Done are not committed work - dropped from the digest."""
    items = [
        _item(1, status="Backlog"),
        _item(2, status="No Status"),
        _item(3, status="Done"),
        _item(4, status="Ready"),
    ]
    g = dd.group_by_role_status(items)
    committed_numbers = [it["number"] for statuses in g.values()
                         for lst in statuses.values() for it in lst]
    assert committed_numbers == [4]


def test_group_roleless_committed_item_goes_to_unassigned():
    """A committed item with no role: label is a hygiene signal - bucket it, don't drop it."""
    g = dd.group_by_role_status([_item(9, role=None, status="Ready")])
    assert g[dd.UNASSIGNED]["Ready"][0]["number"] == 9


# --------------------------------------------------------------------------- #
# aging_items                                                                  #
# --------------------------------------------------------------------------- #

def test_aging_flags_only_committed_items_idle_at_least_threshold():
    items = [
        _item(1, status="Ready", updated="2026-06-01T00:00:00Z"),   # 31d - aging
        _item(2, status="In progress", updated="2026-07-01T00:00:00Z"),  # 1d - fresh
        _item(3, status="Backlog", updated="2026-01-01T00:00:00Z"),  # old but uncommitted
    ]
    aged = dd.aging_items(items, now=NOW, stale_days=14)
    assert [it["number"] for it in aged] == [1]


def test_aging_sorts_most_dormant_first():
    items = [
        _item(1, status="Ready", updated="2026-06-15T00:00:00Z"),   # 17d
        _item(2, status="Ready", updated="2026-05-01T00:00:00Z"),   # 62d
    ]
    aged = dd.aging_items(items, now=NOW, stale_days=14)
    assert [it["number"] for it in aged] == [2, 1]


# --------------------------------------------------------------------------- #
# render_digest                                                                #
# --------------------------------------------------------------------------- #

def test_render_shows_per_role_counts_and_items():
    items = [
        _item(1, role="role:pm", status="Ready", title="alpha"),
        _item(2, role="role:pm", status="In progress", title="beta"),
    ]
    out = dd.render_digest(items, discussions=[], now=NOW, stale_days=14)
    assert "PM" in out
    assert "#1" in out and "alpha" in out
    assert "#2" in out and "beta" in out


def test_render_lists_open_discussions():
    discussions = [{"number": 910, "title": "MM narrative home", "author": "Jin-HoMLee"}]
    out = dd.render_digest([], discussions=discussions, now=NOW, stale_days=14)
    assert "#910" in out and "MM narrative home" in out


def test_render_surfaces_aging_section_when_present():
    items = [_item(1, role="role:pm", status="Ready", updated="2026-05-01T00:00:00Z")]  # 62d
    out = dd.render_digest(items, discussions=[], now=NOW, stale_days=14)
    assert "Aging" in out
    assert "#1" in out


def test_render_is_deterministic_given_fixed_now():
    items = [_item(1, role="role:pm", status="Ready")]
    a = dd.render_digest(items, discussions=[], now=NOW, stale_days=14)
    b = dd.render_digest(items, discussions=[], now=NOW, stale_days=14)
    assert a == b
