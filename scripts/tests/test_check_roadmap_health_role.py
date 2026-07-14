# scripts/tests/test_check_roadmap_health_role.py
#
# Issue #1153, second instance. `board_open_items` was fixed to match --role against
# EVERY role label, but `check_roadmap_health` reuses that same `normalize()` output
# and re-implemented its own first-role-only filter - so the identical silent drop
# survived one file over. Caught by the PR #1154 bot review, not by me: I fixed the
# instance and missed the class, which is the failure mode this repo keeps hitting.
#
# The filter is the blocking one (a dual-role overdue Issue vanishes from
# `--role <r>` for whichever role is not listed first). `group_by_role` is the
# softer sibling in the same file: an overdue item owned by two roles must appear
# under BOTH headings, or one of its owners never sees it.
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import check_roadmap_health as crh  # noqa: E402


def _it(roles, number=1112):
    """A normalize()-shaped item: `role` = first label, `roles` = the full set."""
    labels = [f"role:{r}" for r in roles]
    return {"number": number, "role": labels[0] if labels else None, "roles": labels}


# --- the blocking one: the --role filter must not silently drop a dual-role item ---

def test_role_filter_matches_dual_role_under_its_second_role():
    it = _it(["scientist", "developer"])
    assert crh.matches_role(it, "developer") is True


def test_role_filter_matches_dual_role_under_its_first_role():
    it = _it(["scientist", "developer"])
    assert crh.matches_role(it, "scientist") is True


def test_role_filter_still_says_no_to_an_unrelated_role():
    """It must still be able to fail, or it is not a filter."""
    it = _it(["scientist", "developer"])
    assert crh.matches_role(it, "pm") is False


def test_role_filter_accepts_the_prefixed_form():
    it = _it(["scientist", "developer"])
    assert crh.matches_role(it, "role:developer") is True


def test_role_filter_on_a_single_role_item():
    it = _it(["developer"])
    assert crh.matches_role(it, "developer") is True
    assert crh.matches_role(it, "scientist") is False


# --- the softer sibling: grouping must surface a dual-role item to BOTH owners ---

def test_group_by_role_buckets_a_dual_role_item_under_both_roles():
    groups = crh.group_by_role([_it(["scientist", "developer"])])
    assert set(groups) == {"role:scientist", "role:developer"}
    assert groups["role:developer"][0]["number"] == 1112
    assert groups["role:scientist"][0]["number"] == 1112


def test_group_by_role_single_role_unchanged():
    groups = crh.group_by_role([_it(["developer"])])
    assert set(groups) == {"role:developer"}


def test_group_by_role_roleless_item_goes_to_none_bucket():
    groups = crh.group_by_role([{"number": 7, "role": None, "roles": []}])
    assert set(groups) == {"(none)"}


# --- legacy shape: an item with `role` but no `roles` must NOT vanish ---------
#
# Reading `roles` alone would bucket every first-role-only dict into '(none)' and
# make `--role X` match nothing - the same silent-wrong-answer class this Issue is
# about, reintroduced one level down by its own fix. Caught by a pre-existing test
# in workflow/tests/ (a directory my first grep did not search).

def test_legacy_item_without_roles_key_still_groups_under_its_role():
    groups = crh.group_by_role([{"number": 1, "role": "role:pm"}])
    assert set(groups) == {"role:pm"}


def test_legacy_item_without_roles_key_still_matches_the_filter():
    assert crh.matches_role({"number": 1, "role": "role:pm"}, "pm") is True
    assert crh.matches_role({"number": 1, "role": "role:pm"}, "developer") is False
