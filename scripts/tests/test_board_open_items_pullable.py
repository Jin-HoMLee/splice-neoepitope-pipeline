import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import board_open_items as boi  # noqa: E402


def _raw(number, labels=None, blocked=None, start=None):
    """A minimal GraphQL item node as normalize() expects it."""
    field_values = []
    field_values.append({
        "name": "Ready",
        "field": {"name": "Status"},
    })
    if start is not None:
        field_values.append({"date": start, "field": {"name": "Start date"}})
    return {
        "content": {
            "__typename": "Issue",
            "number": number,
            "title": f"#{number}",
            "state": "OPEN",
            "url": f"https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/{number}",
            "body": "",
            "labels": {"nodes": [{"name": n} for n in (labels or [])]},
            "blockedBy": {"nodes": blocked or []},
            "subIssuesSummary": {"total": 0},
        },
        "fieldValues": {"nodes": field_values},
    }


def test_needs_design_label_sets_not_pullable():
    it = boi.normalize(_raw(100, labels=["role:pm", "needs-design"]))
    assert it["not_pullable"] == "needs-design"


def test_clean_issue_is_pullable():
    it = boi.normalize(_raw(101, labels=["role:pm"]))
    assert it["not_pullable"] is None


def test_open_blocker_sets_not_pullable():
    it = boi.normalize(_raw(102, labels=["role:pm"], blocked=[{"number": 5, "state": "OPEN"}]))
    assert it["not_pullable"] == "blocked-by-issue: #5"
