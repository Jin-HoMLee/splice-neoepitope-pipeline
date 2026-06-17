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

import pytest

# Make top-level scripts/ importable (board_open_items.py is a standalone CLI,
# not a Snakemake script:-invoked module under workflow/scripts/).
_SCRIPTS_DIR = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(_SCRIPTS_DIR))

import board_open_items as boi  # noqa: E402

NOW = datetime(2026, 6, 4, 12, 0, 0, tzinfo=timezone.utc)


def _board_item(number, *, updated="2026-06-04T00:00:00Z",
                created="2026-05-01T00:00:00Z", closed=None, role="role:pm",
                status="Ready", priority=None, size=None, sub_total=None,
                typename="Issue", is_draft=False):
    """A minimal ProjectV2 board node shaped like the GraphQL response.

    `sub_total=None` omits the `subIssuesSummary` block entirely (mimics a node
    the query didn't return one for — e.g. a PR); an int injects
    `subIssuesSummary { total: <n> }` like the Issue fragment does.
    """
    field_values = [{"name": status, "field": {"name": "Status"}}]
    if priority:
        field_values.append({"name": priority, "field": {"name": "Priority"}})
    if size:
        field_values.append({"name": size, "field": {"name": "Size"}})
    content = {
        "__typename": typename,
        "number": number,
        "title": f"issue {number}",
        "state": "OPEN",
        "url": f"https://example/{number}",
        "createdAt": created,
        "updatedAt": updated,
        "closedAt": closed,
        "labels": {"nodes": [{"name": role}]},
    }
    if typename == "PullRequest":
        content["isDraft"] = is_draft
    if sub_total is not None:
        content["subIssuesSummary"] = {"total": sub_total}
    return {
        "content": content,
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
    assert boi.age_days("2026-06-01T00:00:00Z", NOW) == pytest.approx(3.5)


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
    # additive keys reach the JSON output (the PR-body additive-field claim):
    # timestamps (Issue #642) + is_parent (Issue #742)
    assert all(k in parsed[0] for k in ("created_at", "updated_at", "closed_at", "is_parent"))


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


def test_main_stale_days_filters_via_main(monkeypatch, capsys):
    """main() routes --stale-days through apply_recency. Ancient/future fixtures
    keep this robust to the real clock (main() uses datetime.now, not NOW)."""
    raw = [
        _board_item(1, updated="2099-01-01T00:00:00Z"),  # future → age < 14 → dropped
        _board_item(2, updated="2020-01-01T00:00:00Z"),  # ancient → age >= 14 → kept
    ]
    _, out = _run_main(monkeypatch, capsys, ["--stale-days", "14", "--json"], raw)
    assert [it["number"] for it in json.loads(out)] == [2]


def test_main_stale_days_zero_keeps_all_past_items(monkeypatch, capsys):
    """--stale-days 0 (idle >= 0) keeps every past-dated item, oldest-active first —
    it does NOT filter everything out (the documented foot-gun)."""
    raw = [
        _board_item(1, updated="2025-01-01T00:00:00Z"),
        _board_item(2, updated="2020-01-01T00:00:00Z"),
    ]
    _, out = _run_main(monkeypatch, capsys, ["--stale-days", "0", "--json"], raw)
    parsed = json.loads(out)
    assert {it["number"] for it in parsed} == {1, 2}      # both kept
    assert [it["number"] for it in parsed] == [2, 1]       # oldest-active first


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


def test_format_table_now_defaults_to_real_clock():
    # the now=None fallback (datetime.now) renders without error; ancient fixture
    # so the rendered Age is clock-independent
    items = [boi.normalize(_board_item(1, updated="2020-01-01T00:00:00Z"))]
    table = boi.format_table(items)  # no now= → exercises the None branch
    assert "Age" in table.splitlines()[0]
    assert "issue 1" in table


# --- parent-awareness (Issue #742) -----------------------------------------

def test_normalize_is_parent_true_when_subissues():
    n = boi.normalize(_board_item(742, sub_total=3))
    assert n["is_parent"] is True


def test_normalize_is_parent_false_when_zero_subissues():
    n = boi.normalize(_board_item(742, sub_total=0))
    assert n["is_parent"] is False


def test_normalize_is_parent_false_when_summary_absent():
    # a node with no subIssuesSummary block (e.g. the query returned none) → leaf
    n = boi.normalize(_board_item(742, sub_total=None))
    assert n["is_parent"] is False


def test_normalize_pr_is_not_parent():
    # PRs have no sub-issues; the PR fragment carries no subIssuesSummary
    n = boi.normalize(_board_item(742, typename="PullRequest", role="role:pm"))
    assert n["kind"] == "PR"
    assert n["is_parent"] is False


def test_main_json_carries_is_parent(monkeypatch, capsys):
    raw = [_board_item(1, sub_total=2), _board_item(2, sub_total=0)]
    _, out = _run_main(monkeypatch, capsys, ["--json"], raw)
    parsed = {it["number"]: it for it in json.loads(out)}
    assert parsed[1]["is_parent"] is True
    assert parsed[2]["is_parent"] is False


def test_format_table_marks_parent_in_kind():
    parent = boi.normalize(_board_item(1, sub_total=4))
    leaf = boi.normalize(_board_item(2, sub_total=0))
    table = boi.format_table([parent, leaf], now=NOW)
    # parent's Kind cell carries the /P marker, mirroring the /D draft convention
    assert "Issue/P" in table
    # the leaf row stays a plain "Issue" (no stray /P)
    leaf_row = [ln for ln in table.splitlines() if ln.strip().endswith("issue 2")][0]
    assert "/P" not in leaf_row


def test_exclude_parents_filters_out_parents(monkeypatch, capsys):
    raw = [_board_item(1, sub_total=3), _board_item(2, sub_total=0)]
    _, out = _run_main(monkeypatch, capsys, ["--exclude-parents", "--json"], raw)
    nums = [it["number"] for it in json.loads(out)]
    assert nums == [2]  # the parent (#1) is dropped, the leaf (#2) stays


def test_format_table_parent_kind_keeps_column_alignment():
    # "Issue/P" is 7 chars and must not overflow the Kind column and shift the
    # # column out from under its header (the first kind value to exceed the
    # old 5-char budget; "Issue"=5 fit, drafts are "PR/D"=4).
    parent = boi.normalize(_board_item(547, sub_total=4))
    lines = boi.format_table([parent], now=NOW).splitlines()
    header, row = lines[0], lines[2]  # 0=header, 1=separator, 2=data
    assert "Issue/P" in row
    # the issue number sits exactly under the '#' header label
    assert row.index("547") == header.index("#")
