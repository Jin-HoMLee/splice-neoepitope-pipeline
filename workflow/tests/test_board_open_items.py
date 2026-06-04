"""Tests for board_open_items.py timestamp threading + recency flags (Issue #642).

Covers:
  - normalize() surfaces issue created/updated/closed timestamps (AC1)
  - --sort-updated / --stale-days recency logic, keyed on Issue.updatedAt (AC2)
  - age helpers + the table Age column

The recency logic is exercised through pure helpers with an injected `now`,
so the tests are deterministic and never touch the network.
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

# Make top-level scripts/ importable (board_open_items.py is a standalone CLI,
# not a Snakemake script:-invoked module under workflow/scripts/).
_SCRIPTS_DIR = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(_SCRIPTS_DIR))

import board_open_items as boi  # noqa: E402

NOW = datetime(2026, 6, 4, 12, 0, 0, tzinfo=timezone.utc)


def _board_item(number, *, updated="2026-06-04T00:00:00Z",
                created="2026-05-01T00:00:00Z", closed=None, role="role:pm",
                status="Ready", priority=None, size=None):
    """A minimal ProjectV2 board node shaped like the GraphQL response."""
    field_values = [{"name": status, "field": {"name": "Status"}}]
    if priority:
        field_values.append({"name": priority, "field": {"name": "Priority"}})
    if size:
        field_values.append({"name": size, "field": {"name": "Size"}})
    return {
        "content": {
            "__typename": "Issue",
            "number": number,
            "title": f"issue {number}",
            "state": "OPEN",
            "url": f"https://example/{number}",
            "createdAt": created,
            "updatedAt": updated,
            "closedAt": closed,
            "labels": {"nodes": [{"name": role}]},
        },
        "fieldValues": {"nodes": field_values},
    }


def _run_main(monkeypatch, capsys, argv, raw_items):
    """Drive main() with a stubbed board fetch; return (exit_code, stdout)."""
    monkeypatch.setattr(boi, "fetch_all_items", lambda: raw_items)
    monkeypatch.setattr(sys, "argv", ["board_open_items.py", *argv])
    rc = boi.main()
    return rc, capsys.readouterr().out


# --- AC1: timestamps surfaced by normalize ---------------------------------

def test_normalize_includes_timestamps():
    out = boi.normalize(_board_item(642, updated="2026-06-03T09:00:00Z",
                                    created="2026-05-15T08:00:00Z"))
    assert out is not None
    assert out["created_at"] == "2026-05-15T08:00:00Z"
    assert out["updated_at"] == "2026-06-03T09:00:00Z"
    assert out["closed_at"] is None


# --- age helpers -----------------------------------------------------------

def test_age_days_basic():
    # 3 days, 12 hours before NOW
    assert boi.age_days("2026-06-01T00:00:00Z", NOW) == 3.5


def test_age_days_handles_missing():
    assert boi.age_days(None, NOW) is None
    assert boi.age_days("", NOW) is None


def test_age_label():
    assert boi.age_label("2026-06-01T00:00:00Z", NOW) == "3d"
    assert boi.age_label(None, NOW) == "—"


# --- AC2: --sort-updated (momentum, most-recent first) ---------------------

def test_apply_recency_sort_updated_desc():
    items = [
        {"number": 1, "updated_at": "2026-06-01T00:00:00Z"},
        {"number": 2, "updated_at": "2026-06-04T00:00:00Z"},
        {"number": 3, "updated_at": "2026-05-20T00:00:00Z"},
    ]
    out = boi.apply_recency(items, sort_updated=True, stale_days=None, now=NOW)
    assert [it["number"] for it in out] == [2, 1, 3]


def test_apply_recency_missing_updated_sorts_last():
    items = [
        {"number": 1, "updated_at": "2026-06-01T00:00:00Z"},
        {"number": 2, "updated_at": None},
        {"number": 3, "updated_at": "2026-06-03T00:00:00Z"},
    ]
    out = boi.apply_recency(items, sort_updated=True, stale_days=None, now=NOW)
    assert [it["number"] for it in out] == [3, 1, 2]


# --- AC2: --stale-days (dormancy, oldest-first, filtered) ------------------

def test_apply_recency_stale_filter_and_asc_sort():
    items = [
        {"number": 1, "updated_at": "2026-06-04T00:00:00Z"},  # 0.5d — fresh
        {"number": 2, "updated_at": "2026-05-01T00:00:00Z"},  # ~34d — stale
        {"number": 3, "updated_at": "2026-05-20T00:00:00Z"},  # ~15d — stale
    ]
    out = boi.apply_recency(items, sort_updated=False, stale_days=14, now=NOW)
    # only the two >=14d items, oldest activity first
    assert [it["number"] for it in out] == [2, 3]


def test_apply_recency_stale_excludes_missing_updated():
    items = [
        {"number": 1, "updated_at": None},
        {"number": 2, "updated_at": "2026-05-01T00:00:00Z"},
    ]
    out = boi.apply_recency(items, sort_updated=False, stale_days=14, now=NOW)
    assert [it["number"] for it in out] == [2]


def test_apply_recency_combined_stale_filter_with_momentum_sort():
    # both flags: stale filter applies, but --sort-updated forces desc order
    items = [
        {"number": 1, "updated_at": "2026-06-04T00:00:00Z"},  # fresh, dropped
        {"number": 2, "updated_at": "2026-05-01T00:00:00Z"},  # ~34d
        {"number": 3, "updated_at": "2026-05-20T00:00:00Z"},  # ~15d
    ]
    out = boi.apply_recency(items, sort_updated=True, stale_days=14, now=NOW)
    assert [it["number"] for it in out] == [3, 2]


# --- AC2: exact stale boundary (inclusive >=) ------------------------------

def test_apply_recency_stale_boundary_is_inclusive():
    # NOW = 2026-06-04T12:00Z; 14 days earlier = 2026-05-21T12:00Z (exactly 14.0d)
    items = [
        {"number": 1, "updated_at": "2026-05-21T12:00:00Z"},  # exactly 14d → kept (>=)
        {"number": 2, "updated_at": "2026-05-21T11:59:00Z"},  # 14d+1m  → kept
        {"number": 3, "updated_at": "2026-05-21T12:01:00Z"},  # 13d23h59m → dropped
    ]
    out = boi.apply_recency(items, sort_updated=False, stale_days=14, now=NOW)
    nums = {it["number"] for it in out}
    assert 1 in nums          # the boundary itself is inclusive — guards >= vs > regressions
    assert 2 in nums
    assert 3 not in nums


# --- normalize() -> apply_recency wiring (AC1 + AC2 in combination) ---------

def test_normalize_output_feeds_apply_recency():
    items = [
        boi.normalize(_board_item(1, updated="2026-05-01T00:00:00Z")),
        boi.normalize(_board_item(2, updated="2026-06-03T00:00:00Z")),
    ]
    out = boi.apply_recency(items, sort_updated=True, stale_days=None, now=NOW)
    # proves normalize's "updated_at" key is the one apply_recency consumes
    assert [it["number"] for it in out] == [2, 1]


# --- empty input -----------------------------------------------------------

def test_apply_recency_and_table_handle_empty():
    assert boi.apply_recency([], sort_updated=True, stale_days=14, now=NOW) == []
    assert boi.format_table([], now=NOW) == "(no items matched)\n"


# --- main(): the check_ready_queue.sh JSON contract + default ordering ------

def test_main_json_emits_flat_array(monkeypatch, capsys):
    """check_ready_queue.sh runs `--status Ready --json | jq length` — the
    output MUST stay a flat top-level array so `jq length` == item count."""
    raw = [_board_item(1), _board_item(2)]
    rc, out = _run_main(monkeypatch, capsys, ["--status", "Ready", "--json"], raw)
    assert rc == 0
    parsed = json.loads(out)
    assert isinstance(parsed, list)      # not an envelope object
    assert len(parsed) == 2              # what `jq length` would report
    assert {it["number"] for it in parsed} == {1, 2}


def test_main_default_no_flags_uses_sort_key(monkeypatch, capsys):
    """No recency flags → original Status/Priority/Size sort_key ordering, not updatedAt."""
    raw = [
        _board_item(1, status="Ready", updated="2026-06-04T00:00:00Z"),        # recent, Status order 3
        _board_item(2, status="In progress", updated="2026-05-01T00:00:00Z"),  # old, Status order 0
    ]
    _, out = _run_main(monkeypatch, capsys, ["--json"], raw)
    # "In progress" sorts above "Ready" regardless of updatedAt
    assert [it["number"] for it in json.loads(out)] == [2, 1]


def test_main_sort_updated_overrides_sort_key(monkeypatch, capsys):
    raw = [
        _board_item(1, status="Ready", updated="2026-06-04T00:00:00Z"),
        _board_item(2, status="In progress", updated="2026-05-01T00:00:00Z"),
    ]
    _, out = _run_main(monkeypatch, capsys, ["--sort-updated", "--json"], raw)
    # momentum: most-recent first, Status ordering ignored
    assert [it["number"] for it in json.loads(out)] == [1, 2]


# --- table Age column ------------------------------------------------------

def test_age_label_clamps_future_to_zero():
    # clock skew: a `now` behind the GitHub timestamp must not render negative
    assert boi.age_label("2026-06-10T00:00:00Z", NOW) == "0d"


def test_format_table_age_column_position_and_missing_render():
    items = [
        boi.normalize(_board_item(1, updated="2026-06-01T00:00:00Z")),  # 3d
        boi.normalize(_board_item(2, updated=None)),                     # missing → —
    ]
    table = boi.format_table(items, now=NOW)
    header = table.splitlines()[0].split()
    assert header.index("Age") == header.index("Sz") + 1   # Age sits right after Sz
    assert header.index("Age") < header.index("Role")
    assert "3d" in table
    assert "—" in table   # missing-timestamp renders in the TABLE, not just the helper
