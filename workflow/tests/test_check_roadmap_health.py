"""Tests for check_roadmap_health.py — the roadmap-overdue (Target-date) sweep.

The overdue logic is pure (injected `today`), so it's deterministic and never
touches the network; main() is driven with a stubbed board fetch (Issue #704).
"""
import json
import sys
from datetime import date
from pathlib import Path

# Make top-level scripts/ importable (standalone CLIs, not workflow/scripts).
_SCRIPTS_DIR = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(_SCRIPTS_DIR))

import board_open_items as boi  # noqa: E402
import check_roadmap_health as crh  # noqa: E402

TODAY = date(2026, 6, 21)


def _raw(number, target, *, role="role:pm", status="Ready"):
    """A minimal ProjectV2 board node (Issue) with a Target-date field value."""
    fvs = [{"name": status, "field": {"name": "Status"}}]
    if target:
        fvs.append({"date": target, "field": {"name": "Target date"}})
    return {
        "content": {
            "__typename": "Issue", "number": number, "title": f"i{number}",
            "state": "OPEN", "url": f"https://example/{number}",
            "createdAt": "2026-05-01T00:00:00Z", "updatedAt": "2026-06-01T00:00:00Z",
            "closedAt": None, "labels": {"nodes": [{"name": role}]},
        },
        "fieldValues": {"nodes": fvs},
    }


# --- find_overdue (pure) ---------------------------------------------------

def test_find_overdue_keeps_only_past_targets():
    items = [
        {"number": 1, "target_date": "2026-06-01", "role": "role:pm"},  # past
        {"number": 2, "target_date": "2026-06-30", "role": "role:pm"},  # future
        {"number": 3, "target_date": None, "role": "role:pm"},          # no target
    ]
    out = crh.find_overdue(items, TODAY)
    assert [i["number"] for i in out] == [1]


def test_find_overdue_today_is_not_overdue():
    # AC: "Target date < today" — the date itself is the deadline, not overdue yet.
    items = [{"number": 1, "target_date": "2026-06-21", "role": "role:pm"}]
    assert crh.find_overdue(items, TODAY) == []


def test_find_overdue_sorted_most_overdue_first():
    items = [
        {"number": 1, "target_date": "2026-06-10", "role": "role:pm"},
        {"number": 2, "target_date": "2026-06-01", "role": "role:pm"},
        {"number": 3, "target_date": "2026-06-15", "role": "role:pm"},
    ]
    out = crh.find_overdue(items, TODAY)
    assert [i["number"] for i in out] == [2, 1, 3]  # oldest Target first


# --- group_by_role (pure) --------------------------------------------------

def test_group_by_role_buckets_and_handles_missing():
    items = [
        {"number": 1, "role": "role:pm"},
        {"number": 2, "role": "role:scientist"},
        {"number": 3, "role": None},
    ]
    g = crh.group_by_role(items)
    assert [i["number"] for i in g["role:pm"]] == [1]
    assert [i["number"] for i in g["role:scientist"]] == [2]
    assert [i["number"] for i in g["(none)"]] == [3]


# --- main() exit-code contract (0 clear / 2 overdue) -----------------------

def test_main_exits_2_when_overdue(monkeypatch, capsys):
    monkeypatch.setattr(boi, "fetch_all_items", lambda: [_raw(1, "2020-01-01")])
    monkeypatch.setattr(sys, "argv", ["check_roadmap_health.py"])
    rc = crh.main()
    assert rc == 2
    assert "i1" in capsys.readouterr().out


def test_main_exits_0_when_clear(monkeypatch, capsys):
    # future Target → no overdue; ancient/future fixtures stay robust to the real clock
    monkeypatch.setattr(boi, "fetch_all_items", lambda: [_raw(1, "2099-01-01")])
    monkeypatch.setattr(sys, "argv", ["check_roadmap_health.py"])
    assert crh.main() == 0


def test_main_role_filter(monkeypatch, capsys):
    raw = [_raw(1, "2020-01-01", role="role:pm"),
           _raw(2, "2020-01-01", role="role:scientist")]
    monkeypatch.setattr(boi, "fetch_all_items", lambda: raw)
    monkeypatch.setattr(sys, "argv", ["check_roadmap_health.py", "--role", "scientist"])
    rc = crh.main()
    out = capsys.readouterr().out
    assert rc == 2
    assert "i2" in out and "i1" not in out


def test_main_json_emits_flat_array(monkeypatch, capsys):
    raw = [_raw(1, "2020-01-01"), _raw(2, "2099-01-01")]  # only #1 overdue
    monkeypatch.setattr(boi, "fetch_all_items", lambda: raw)
    monkeypatch.setattr(sys, "argv", ["check_roadmap_health.py", "--json"])
    rc = crh.main()
    parsed = json.loads(capsys.readouterr().out)
    assert rc == 2
    assert isinstance(parsed, list)
    assert [i["number"] for i in parsed] == [1]
    assert parsed[0]["target_date"] == "2020-01-01"
