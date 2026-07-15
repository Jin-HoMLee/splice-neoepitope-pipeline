# scripts/tests/test_dispatch_digest_role.py
#
# Issue #1139 - the LAST site of the #1153 first-role-only class.
#
# `board_open_items` (#1153) and `check_roadmap_health` (#1154) were both fixed to
# key off EVERY role label. `dispatch_digest._role_slug` still read the singular
# scalar `item["role"]`, so a multi-role item was bucketed under exactly one role
# and was INVISIBLE to the other. The digest therefore under-reported that role's
# Ready depth, which is what manufactures a phantom `[REPLENISH role]` shortfall
# against `check_ready_queue.sh` - the gate that counts by every role label and so
# disagrees with the digest.
#
# The failure is silent in the worst way: the digest returns a clean, short,
# entirely plausible list. Nothing crashes. An under-count reads exactly like a
# small backlog.
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "pm"))
import dispatch_digest as dd  # noqa: E402


def _item(number, roles, status="Ready"):
    """A board item shaped as board_open_items.normalize() emits it."""
    return {
        "number": number,
        "title": f"item {number}",
        "url": f"https://example.invalid/{number}",
        "status": status,
        # `role` is the first label (what the buggy code read); `roles` is the
        # authoritative full list, exactly as board_open_items sets it.
        "role": roles[0] if roles else None,
        "roles": list(roles),
    }


def test_multirole_item_appears_under_every_role():
    """The blocking one: a dual-role item must be visible to BOTH owners."""
    grouped = dd.group_by_role_status(
        [_item(1, ["role:scientist", "role:developer"])]
    )
    assert "scientist" in grouped, "item vanished from its first role"
    assert "developer" in grouped, (
        "item is invisible to its SECOND role - this is the #1153 under-count: "
        "the digest silently under-reports that role's Ready depth"
    )
    assert grouped["scientist"]["Ready"][0]["number"] == 1
    assert grouped["developer"]["Ready"][0]["number"] == 1


def test_single_role_item_is_unchanged():
    """Matched-pair control: the single-role path must not regress."""
    grouped = dd.group_by_role_status([_item(2, ["role:developer"])])
    assert list(grouped) == ["developer"]
    assert grouped["developer"]["Ready"][0]["number"] == 2


def test_roleless_committed_item_still_surfaces_as_unassigned():
    """A role-less committed item is a hygiene signal, not something to drop."""
    grouped = dd.group_by_role_status([_item(3, [])])
    assert dd.UNASSIGNED in grouped


def test_digest_role_depth_matches_a_by_every_label_count():
    """The regression that motivated the Issue, stated as the gate states it.

    `check_ready_queue.sh` counts an item under EVERY role label it carries. The
    digest must agree, or a role reads as short when it is not.
    """
    items = [
        _item(10, ["role:developer"]),
        _item(11, ["role:scientist", "role:developer"]),
        _item(12, ["role:pm", "role:developer"]),
    ]
    grouped = dd.group_by_role_status(items)
    dev_ready = grouped["developer"]["Ready"]
    assert [it["number"] for it in dev_ready] == [10, 11, 12], (
        "developer owns all three items; a first-label-only bucket sees only one, "
        "under-reporting the lane by 2 and manufacturing a phantom shortfall"
    )


def test_uncommitted_statuses_are_still_dropped():
    """Backlog / Done are not active work; multi-role must not smuggle them in."""
    grouped = dd.group_by_role_status(
        [_item(4, ["role:pm", "role:developer"], status="Backlog")]
    )
    assert grouped == {}
