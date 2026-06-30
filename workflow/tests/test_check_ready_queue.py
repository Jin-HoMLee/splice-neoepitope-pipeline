"""Tests for scripts/check_ready_queue.sh - demand-aware per-role Ready trigger.

History:
  - Issue #754 introduced a per-role floor (5) + total cap (18) gate.
  - Issue #902 facet 1 made the floor *demand-aware*: it fires REPLENISH for a
    role only when that role is actually consuming (>= 1 item In progress) AND
    its Ready count is below the floor. A role idle on a quiet/post-ship board
    (0 In progress) no longer nags, even at 0 Ready. Floor 5 -> 3 (one WIP
    slate), cap 18 -> 12 (floor x 3 + 3 headroom).

The script reads the board via board_open_items.py, which hits the network.
For deterministic offline tests it accepts a fixture seam: when
BOARD_ITEMS_JSON_FILE points at a JSON array of {number, status, labels}
items, the script reads from that file instead of calling board_open_items.py.
The script itself filters that array by status (Ready / In progress), so a
fixture controls both the Ready buffer and the In-progress demand signal.

Gate (Issue #902 facet 1):
  - per-role floor 3 for PM / Scientist / Developer (MM excluded)
  - REPLENISH fires only for a role with ready < floor AND in_progress >= 1
  - a role with ready < floor AND in_progress == 0 is "quiet" (holding), not
    a replenish trigger - no needs-attention
  - total cap 12 on the Ready buffer (soft WIP limit), independent of demand
  - exit 2 if any role needs replenish OR total at/over cap; else exit 0
"""
import json
import os
import subprocess
import tempfile
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "check_ready_queue.sh"


def _item(number, *roles, status="Ready"):
    return {"number": number, "status": status,
            "labels": [f"role:{r}" for r in roles]}


def _run(items, *args):
    """Run the script with a JSON fixture injected via BOARD_ITEMS_JSON_FILE."""
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        json.dump(items, f)
        path = f.name
    try:
        env = {**os.environ, "BOARD_ITEMS_JSON_FILE": path}
        return subprocess.run(["bash", str(SCRIPT), *args], env=env,
                              capture_output=True, text=True)
    finally:
        os.unlink(path)


def _spread(role, n, start, status="Ready"):
    return [_item(start + i, role, status=status) for i in range(n)]


def test_consuming_role_below_floor_replenishes():
    # PM has 2 Ready and 1 In progress (consuming) -> REPLENISH pm.
    # Sci/Dev sit at floor (3 Ready) and are silent.
    items = (_spread("pm", 2, 100) + _spread("pm", 1, 150, status="In progress")
             + _spread("scientist", 3, 200) + _spread("developer", 3, 300))
    r = _run(items)
    assert "[REPLENISH pm: 2 < 3]" in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)


def test_quiet_role_below_floor_does_not_replenish():
    # PM has 0 Ready and 0 In progress - the post-ship quiet board. Thin Ready
    # is correct; no REPLENISH nag. Sci/Dev at floor. Healthy, exit 0.
    items = _spread("scientist", 3, 200) + _spread("developer", 3, 300)
    r = _run(items)
    assert "[REPLENISH" not in r.stdout, r.stdout
    assert "[quiet pm: 0 Ready, 0 In progress]" in r.stdout, r.stdout
    assert "healthy" in r.stdout, r.stdout
    assert r.returncode == 0, (r.returncode, r.stdout, r.stderr)


def test_consuming_role_at_zero_ready_replenishes():
    # PM has 0 Ready but 2 In progress (actively consuming, buffer empty) ->
    # this is exactly when to replenish. Sci/Dev at floor.
    items = (_spread("pm", 2, 150, status="In progress")
             + _spread("scientist", 3, 200) + _spread("developer", 3, 300))
    r = _run(items)
    assert "[REPLENISH pm: 0 < 3]" in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)


def test_all_roles_above_floor_healthy():
    # 3 disjoint Ready items per role = 9 total: every floor met, under cap 12.
    items = _spread("pm", 3, 100) + _spread("scientist", 3, 200) + _spread("developer", 3, 300)
    r = _run(items)
    assert "healthy" in r.stdout, r.stdout
    assert "[REPLENISH" not in r.stdout and "[CAP]" not in r.stdout, r.stdout
    assert r.returncode == 0, (r.returncode, r.stdout, r.stderr)


def test_total_at_cap_flags_cap_and_exit_2():
    # Disjoint pm=4/sci=4/dev=4 Ready -> no floor breach, but total=12 == cap.
    items = _spread("pm", 4, 100) + _spread("scientist", 4, 200) + _spread("developer", 4, 300)
    r = _run(items)
    assert "[CAP] Ready at 12 (>= 12)" in r.stdout, r.stdout
    assert "[REPLENISH" not in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)


def test_empty_board_is_quiet_not_replenish():
    # The whole board is empty (post-ship). No demand anywhere -> no nag.
    # This is the key #902 facet-1 reversal of the old floor-5 behavior.
    r = _run([])
    assert "[REPLENISH" not in r.stdout, r.stdout
    assert "healthy" in r.stdout, r.stdout
    assert r.returncode == 0, (r.returncode, r.stdout, r.stderr)


def test_memory_manager_excluded_from_floor():
    # pm/sci/dev met via shared-role items; MM has only 1 Ready + In progress
    # but is excluded from the floor entirely.
    items = [_item(400 + i, "pm", "scientist", "developer") for i in range(3)]
    items += [_item(500, "memory_manager"),
              _item(501, "memory_manager", status="In progress")]
    r = _run(items)
    assert "memory_manager" not in r.stdout, r.stdout
    assert "[REPLENISH" not in r.stdout, r.stdout
    assert r.returncode == 0, (r.returncode, r.stdout, r.stderr)


def test_multiple_consuming_roles_below_floor_each_reported():
    # PM=1 Ready+1 ip, Sci=0 Ready+1 ip -> two replenish lines. Dev at floor.
    items = (_spread("pm", 1, 100) + _spread("pm", 1, 150, status="In progress")
             + _spread("scientist", 1, 250, status="In progress")
             + _spread("developer", 3, 300))
    r = _run(items)
    assert "[REPLENISH pm: 1 < 3]" in r.stdout, r.stdout
    assert "[REPLENISH scientist: 0 < 3]" in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)


def test_consuming_role_below_floor_with_others_quiet():
    # PM consuming + below floor -> REPLENISH. Sci below floor but quiet (0 ip)
    # -> NO replenish. Mixed board: exactly one replenish line.
    items = (_spread("pm", 1, 100) + _spread("pm", 1, 150, status="In progress")
             + _spread("scientist", 1, 200)  # below floor, but no In progress
             + _spread("developer", 3, 300))
    r = _run(items)
    assert "[REPLENISH pm: 1 < 3]" in r.stdout, r.stdout
    assert "[REPLENISH scientist" not in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)


def test_floor_flag_lowers_threshold():
    # PM=2 Ready + 1 ip would breach floor 3, but --floor 2 makes it not below.
    items = (_spread("pm", 2, 100) + _spread("pm", 1, 150, status="In progress")
             + _spread("scientist", 3, 200) + _spread("developer", 3, 300))
    r = _run(items, "--floor", "2")
    assert "[REPLENISH" not in r.stdout, r.stdout
    assert "[CAP]" not in r.stdout, r.stdout  # total Ready 5 < cap 12
    assert r.returncode == 0, (r.returncode, r.stdout, r.stderr)


def test_cap_flag_lowers_threshold():
    # 3/3/3 disjoint Ready = 9; default cap 12 is fine, but --cap 9 trips it.
    items = _spread("pm", 3, 100) + _spread("scientist", 3, 200) + _spread("developer", 3, 300)
    r = _run(items, "--cap", "9")
    assert "[CAP] Ready at 9 (>= 9)" in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)


def test_invalid_floor_exits_1():
    r = _run([], "--floor", "abc")
    assert r.returncode == 1, (r.returncode, r.stdout, r.stderr)
    assert "floor must be" in r.stderr, r.stderr


def test_breakdown_printed_with_ready_and_inprogress():
    # The breakdown is self-documenting: shows both the Ready buffer and the
    # In-progress demand signal so a reader can see why a role did/didn't fire.
    items = (_spread("pm", 2, 100) + _spread("pm", 1, 150, status="In progress")
             + _spread("scientist", 3, 200) + _spread("developer", 3, 300))
    r = _run(items)
    assert "Ready by role:" in r.stdout, r.stdout
    assert "pm=2" in r.stdout, r.stdout
    assert "In progress by role:" in r.stdout, r.stdout
    assert "pm=1" in r.stdout, r.stdout


def test_help_has_no_code_leak():
    r = subprocess.run(["bash", str(SCRIPT), "--help"], capture_output=True, text=True)
    assert r.returncode == 0, (r.returncode, r.stderr)
    assert "set -euo pipefail" not in r.stdout, r.stdout
    assert "SCRIPT_DIR=" not in r.stdout, r.stdout
