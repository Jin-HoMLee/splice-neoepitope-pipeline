# scripts/tests/test_check_roadmap_health_epic.py
#
# Issue #1149 - the roadmap-overdue sweep listed any open item with a past Target
# date, with no regard for Status. But a parent/epic is PARKED in the off-ladder
# `Epic` Status and carries neither a milestone nor a Target by convention
# (Issue #690 sub-question A). So a leftover Target on a parked epic surfaced as a
# missed DEADLINE (parent epic #679 read as "3 days overdue" on 2026-07-14) when it
# is really a DATA-HYGIENE finding whose remedy is "clear the date", not "ship it".
#
# Treatment 2 of the two PM offered: route it to its own section rather than drop
# it, so the finding survives under its correct frame instead of being suppressed.
#
# The last test is the falsifier the Issue explicitly asks for: a non-Epic leaf with
# a past Target MUST still be reported. Without it, "no overdue items" would be
# indistinguishable from a filter that suppresses everything.
import os
import sys
from datetime import date

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import check_roadmap_health as crh  # noqa: E402

TODAY = date(2026, 7, 14)


def _item(number, status, target, roles=("role:developer",)):
    return {
        "number": number,
        "title": f"item {number}",
        "url": f"https://example.invalid/{number}",
        "status": status,
        "target_date": target,
        "role": roles[0] if roles else None,
        "roles": list(roles),
    }


class TestEpicParkedParents:
    def test_epic_with_past_target_is_not_a_deadline_finding(self):
        """The #679 shape: a parked epic is not 'overdue', it is mis-dated."""
        items = [_item(679, "Epic", "2026-07-11")]
        assert crh.find_overdue(items, TODAY) == [], (
            "an Epic-parked parent carries no Target by convention, so a leftover "
            "date is a hygiene finding, not a missed deadline"
        )

    def test_epic_with_past_target_is_surfaced_as_a_hygiene_finding(self):
        """Routed, not dropped - the finding survives under the correct frame."""
        items = [_item(679, "Epic", "2026-07-11")]
        flagged = crh.find_parent_targets(items)
        assert [it["number"] for it in flagged] == [679]

    def test_epic_without_a_target_is_flagged_by_neither(self):
        """The correct steady state for a parked epic: no Target at all."""
        items = [_item(680, "Epic", None)]
        assert crh.find_overdue(items, TODAY) == []
        assert crh.find_parent_targets(items) == []

    def test_leaf_with_past_target_is_still_overdue(self):
        """THE FALSIFIER. Flip only the Status: this must still be reported.

        If this passed while the Epic cases also passed, the filter would just be
        suppressing everything, and a green sweep would mean nothing.
        """
        items = [_item(700, "In progress", "2026-07-11")]
        overdue = crh.find_overdue(items, TODAY)
        assert [it["number"] for it in overdue] == [700]
        assert crh.find_parent_targets(items) == [], "a leaf is not a parent-target finding"

    def test_mixed_board_splits_into_the_two_frames(self):
        items = [
            _item(679, "Epic", "2026-07-11"),
            _item(700, "In progress", "2026-07-11"),
            _item(701, "Ready", "2026-07-20"),  # future Target: neither
        ]
        assert [it["number"] for it in crh.find_overdue(items, TODAY)] == [700]
        assert [it["number"] for it in crh.find_parent_targets(items)] == [679]


class TestReport:
    def test_report_names_the_remedy_for_a_parent_target(self):
        report = crh.format_report(
            [], TODAY, parent_targets=[_item(679, "Epic", "2026-07-11")]
        )
        assert "679" in report
        assert "clear" in report.lower(), "the report must name the remedy: clear the date"

    def test_all_clear_report_is_unchanged_when_nothing_is_flagged(self):
        report = crh.format_report([], TODAY, parent_targets=[])
        assert "all clear" in report.lower()
