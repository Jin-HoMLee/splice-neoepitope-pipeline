import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts", "pm"))

import pullability  # noqa: E402


def _item(**over):
    base = {"labels": [], "blocked_by": [], "start_date": None}
    base.update(over)
    return base


def test_pullable_when_nothing_gates():
    assert pullability.assess(_item(), today="2026-07-24") is None


def test_open_blocker_gates():
    it = _item(blocked_by=[{"number": 42, "state": "OPEN"}])
    assert pullability.assess(it, today="2026-07-24") == "blocked-by-issue: #42"


def test_closed_blocker_does_not_gate():
    it = _item(blocked_by=[{"number": 42, "state": "CLOSED"}])
    assert pullability.assess(it, today="2026-07-24") is None


def test_needs_design_label_gates():
    assert pullability.assess(_item(labels=["needs-design"]), today="2026-07-24") == "needs-design"


def test_trigger_gated_label_gates():
    assert pullability.assess(_item(labels=["trigger-gated"]), today="2026-07-24") == "trigger-gated"


def test_future_start_date_gates():
    it = _item(start_date="2026-12-01")
    assert pullability.assess(it, today="2026-07-24") == "date-gated: 2026-12-01"


def test_past_start_date_does_not_gate():
    it = _item(start_date="2026-01-01")
    assert pullability.assess(it, today="2026-07-24") is None


def test_ordering_blocker_wins_over_date():
    it = _item(blocked_by=[{"number": 9, "state": "OPEN"}], start_date="2026-12-01")
    assert pullability.assess(it, today="2026-07-24") == "blocked-by-issue: #9"


def test_fails_open_on_empty_item():
    assert pullability.assess({}, today="2026-07-24") is None


def test_fails_open_on_non_dict_blocker_element():
    assert pullability.assess({"blocked_by": [42]}, today="2026-07-24") is None


def test_fails_open_on_non_list_blocked_by():
    assert pullability.assess({"blocked_by": "oops"}, today="2026-07-24") is None
