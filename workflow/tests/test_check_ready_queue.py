"""Tests for scripts/check_ready_queue.sh — per-role Ready floor + total cap (Issue #754).

The script reads the board via board_open_items.py, which hits the network.
For deterministic offline tests it accepts a fixture seam: when
READY_QUEUE_JSON_FILE points at a JSON array, the script reads Ready items
from that file instead of calling board_open_items.py. We feed crafted
fixtures and assert the per-role [REPLENISH], [CAP], and exit-code behavior.

Gate (Issue #754):
  - per-role floor 5 for PM / Scientist / Developer (MM excluded)
  - total cap 15 on the Ready buffer
  - exit 2 if any role below floor OR total at/over cap; else exit 0
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


def _run(items):
    """Run the script with a JSON fixture injected via READY_QUEUE_JSON_FILE."""
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        json.dump(items, f)
        path = f.name
    try:
        env = {**os.environ, "READY_QUEUE_JSON_FILE": path}
        return subprocess.run(["bash", str(SCRIPT)], env=env,
                              capture_output=True, text=True)
    finally:
        os.unlink(path)


def _spread(role, n, start):
    return [_item(start + i, role) for i in range(n)]


def test_role_below_floor_emits_per_role_replenish_and_exit_2():
    # PM=2 (below floor 5), Sci=5, Dev=5
    items = _spread("pm", 2, 100) + _spread("scientist", 5, 200) + _spread("developer", 5, 300)
    r = _run(items)
    assert "[REPLENISH pm: 2 < 5]" in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)


def test_all_roles_at_floor_under_cap_is_healthy():
    # 5 items each carrying all three roles -> pm=sci=dev=5, total=5 (< cap 15)
    items = [_item(400 + i, "pm", "scientist", "developer") for i in range(5)]
    r = _run(items)
    assert "healthy" in r.stdout, r.stdout
    assert "[REPLENISH" not in r.stdout, r.stdout
    assert r.returncode == 0, (r.returncode, r.stdout, r.stderr)


def test_total_at_cap_flags_cap_and_exit_2():
    # Disjoint: pm=5, sci=5, dev=5 -> no floor breach, but total=15 == cap
    items = _spread("pm", 5, 100) + _spread("scientist", 5, 200) + _spread("developer", 5, 300)
    r = _run(items)
    assert "[CAP] Ready at 15 (>= 15)" in r.stdout, r.stdout
    assert "[REPLENISH" not in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)


def test_memory_manager_excluded_from_floor():
    # pm/sci/dev all met via shared-role items; MM has only 1 (< floor) but is excluded
    items = [_item(400 + i, "pm", "scientist", "developer") for i in range(5)]
    items += [_item(500, "memory_manager")]
    r = _run(items)
    assert "memory_manager" not in r.stdout, r.stdout
    assert "[REPLENISH" not in r.stdout, r.stdout
    assert r.returncode == 0, (r.returncode, r.stdout, r.stderr)


def test_multiple_roles_below_floor_each_reported():
    # PM=1, Sci=0, Dev=5 -> two replenish lines
    items = _spread("pm", 1, 100) + _spread("developer", 5, 300)
    r = _run(items)
    assert "[REPLENISH pm: 1 < 5]" in r.stdout, r.stdout
    assert "[REPLENISH scientist: 0 < 5]" in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)
