"""Tests for scripts/check_ready_queue.sh - proactive per-role floor with
shortfall routing (commit vs groom).

History:
  - Issue #754 introduced a per-role floor (5) + total cap (18) gate.
  - Issue #902 facet 1 (first cut) made the floor demand-aware (gated on
    In progress). That over-corrected: it deleted the stocked-shelf benefit -
    an idle role found an empty Ready queue instead of a curated shortlist.
  - #902 facet 1 (refined 2026-06-30) restores a PROACTIVE floor of 5 (the
    shelf depth the roles liked) and instead makes the *shortfall* interpretable:
      * [REPLENISH role]  : short AND the role has Backlog candidates -> commit
                            the highest-priority DoR-ready ones (never stuff).
      * [GROOMING-GAP role]: short AND NO Backlog candidates -> groom/intake.
    In progress is reported as demand context only; it does NOT gate the floor.

The script reads the board via board_open_items.py, which hits the network. For
deterministic offline tests it accepts a fixture seam: when BOARD_ITEMS_JSON_FILE
points at a JSON array of {number, status, labels} items, the script reads from
that file and filters by status (Ready / In progress / Backlog) itself, so a
fixture controls the buffer, the demand context, and the candidate pool.

Gate:
  - per-role floor 5 for PM / Scientist / Developer (MM excluded)
  - a role with ready < floor -> REPLENISH if it has >= 1 Backlog candidate,
    else GROOMING-GAP; either way it needs attention
  - In progress never affects the floor decision (context only)
  - total cap 18 on the Ready buffer (soft WIP limit)
  - exit 2 if any role is short OR total at/over cap; else exit 0
"""
import json
import os
import subprocess
import tempfile
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "check_ready_queue.sh"


def _item(number, *roles, status="Ready", arc_active=False):
    labels = [f"role:{r}" for r in roles]
    if arc_active:
        labels.append("arc-phase:active")
    return {"number": number, "status": status, "labels": labels}


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


def test_short_with_backlog_candidates_replenishes():
    # PM short (2 < 5) and has Backlog candidates -> REPLENISH. Sci/Dev at floor.
    items = (_spread("pm", 2, 100) + _spread("pm", 5, 150, status="Backlog")
             + _spread("scientist", 5, 200) + _spread("developer", 5, 300))
    r = _run(items)
    assert "[REPLENISH pm: 2 < 5]" in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)


def test_short_with_no_backlog_candidates_is_grooming_gap():
    # PM short (2 < 5) with ZERO pm Backlog candidates -> GROOMING-GAP, not a
    # stuff-junk nag. Sci/Dev at floor.
    items = (_spread("pm", 2, 100)
             + _spread("scientist", 5, 200) + _spread("developer", 5, 300))
    r = _run(items)
    assert "[GROOMING-GAP pm: 2 < 5]" in r.stdout, r.stdout
    assert "[REPLENISH pm" not in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)


def test_floor_is_proactive_not_gated_on_in_progress():
    # PM short (2 < 5) with 0 In progress still fires (proactive shelf) - the
    # key reversal of the first-cut demand gate. Backlog candidates present.
    items = (_spread("pm", 2, 100) + _spread("pm", 4, 150, status="Backlog")
             + _spread("scientist", 5, 200) + _spread("developer", 5, 300))
    r = _run(items)
    assert "[REPLENISH pm: 2 < 5]" in r.stdout, r.stdout
    assert "In progress:" in r.stdout and "pm=0" in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)


def test_in_progress_does_not_force_attention_when_stocked():
    # A role at floor with In-progress work is healthy; In progress is context.
    items = ([_item(400 + i, "pm", "scientist", "developer") for i in range(5)]
             + _spread("pm", 2, 500, status="In progress"))
    r = _run(items)
    assert "healthy" in r.stdout, r.stdout
    assert "[REPLENISH" not in r.stdout and "[GROOMING-GAP" not in r.stdout, r.stdout
    assert r.returncode == 0, (r.returncode, r.stdout, r.stderr)


def test_all_roles_at_floor_healthy():
    # 5 items each carrying all three roles -> pm=sci=dev=5, total 5 < cap 18.
    items = [_item(400 + i, "pm", "scientist", "developer") for i in range(5)]
    r = _run(items)
    assert "healthy" in r.stdout, r.stdout
    assert "[REPLENISH" not in r.stdout and "[GROOMING-GAP" not in r.stdout, r.stdout
    assert r.returncode == 0, (r.returncode, r.stdout, r.stderr)


def test_replenish_prefers_active_arc_backlog_candidates():
    # PM short (2 < 5) with 3 Backlog candidates, 2 of them on an active arc.
    # The REPLENISH nudge flags the active-arc ones to prefer (Issue #931).
    items = (_spread("pm", 2, 100)
             + [_item(150, "pm", status="Backlog", arc_active=True),
                _item(151, "pm", status="Backlog", arc_active=True),
                _item(152, "pm", status="Backlog")]
             + _spread("scientist", 5, 200) + _spread("developer", 5, 300))
    r = _run(items)
    assert "[REPLENISH pm: 2 < 5]" in r.stdout, r.stdout
    assert "3 Backlog candidate(s)" in r.stdout, r.stdout
    assert "2 on an active arc" in r.stdout, r.stdout
    assert "prefer these" in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)


def test_replenish_notes_when_no_active_arc_candidates():
    # PM short with Backlog candidates but NONE on an active arc -> the nudge
    # says so (commit the top DoR-ready pool item / consider promoting an arc).
    items = (_spread("pm", 2, 100) + _spread("pm", 3, 150, status="Backlog")
             + _spread("scientist", 5, 200) + _spread("developer", 5, 300))
    r = _run(items)
    assert "[REPLENISH pm: 2 < 5]" in r.stdout, r.stdout
    assert "none on an active arc" in r.stdout, r.stdout
    assert "prefer these" not in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)


def test_grooming_gap_unaffected_by_arc_awareness():
    # Short with ZERO Backlog candidates stays GROOMING-GAP (arc logic only
    # applies on the REPLENISH branch, which requires >= 1 candidate).
    items = (_spread("pm", 2, 100)
             + _spread("scientist", 5, 200) + _spread("developer", 5, 300))
    r = _run(items)
    assert "[GROOMING-GAP pm: 2 < 5]" in r.stdout, r.stdout
    assert "active arc" not in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)


def test_total_at_cap_flags_cap_and_exit_2():
    # Disjoint pm=6/sci=6/dev=6 Ready -> no floor breach (all >= 5), total 18 == cap.
    items = _spread("pm", 6, 100) + _spread("scientist", 6, 200) + _spread("developer", 6, 300)
    r = _run(items)
    assert "[CAP] Ready at 18 (>= 18)" in r.stdout, r.stdout
    assert "[REPLENISH" not in r.stdout and "[GROOMING-GAP" not in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)


def test_empty_board_grooming_gaps_all_three():
    # Nothing anywhere -> every role short with no candidates -> GROOMING-GAP x3
    # (the board has no committable work at all; intake is the remedy).
    r = _run([])
    for role in ("pm", "scientist", "developer"):
        assert f"[GROOMING-GAP {role}: 0 < 5]" in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)


def test_memory_manager_excluded_from_floor():
    # pm/sci/dev at floor; MM short with a Backlog candidate but excluded.
    items = [_item(400 + i, "pm", "scientist", "developer") for i in range(5)]
    items += [_item(500, "memory_manager"),
              _item(501, "memory_manager", status="Backlog")]
    r = _run(items)
    assert "memory_manager" not in r.stdout, r.stdout
    assert "[REPLENISH" not in r.stdout and "[GROOMING-GAP" not in r.stdout, r.stdout
    assert r.returncode == 0, (r.returncode, r.stdout, r.stderr)


def test_multiple_roles_short_route_independently():
    # PM short with candidates -> REPLENISH; Sci short with none -> GROOMING-GAP;
    # Dev at floor -> silent.
    items = (_spread("pm", 1, 100) + _spread("pm", 3, 150, status="Backlog")
             + _spread("scientist", 0, 250)  # no sci items at all
             + _spread("developer", 5, 300))
    r = _run(items)
    assert "[REPLENISH pm: 1 < 5]" in r.stdout, r.stdout
    assert "[GROOMING-GAP scientist: 0 < 5]" in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)


def test_floor_flag_lowers_threshold():
    # PM=2 Ready would breach floor 5, but --floor 2 makes it not short.
    items = (_spread("pm", 2, 100)
             + _spread("scientist", 2, 200) + _spread("developer", 2, 300))
    r = _run(items, "--floor", "2")
    assert "[REPLENISH" not in r.stdout and "[GROOMING-GAP" not in r.stdout, r.stdout
    assert "[CAP]" not in r.stdout, r.stdout  # total Ready 6 < cap 18
    assert r.returncode == 0, (r.returncode, r.stdout, r.stderr)


def test_cap_flag_lowers_threshold():
    # 5/5/5 disjoint Ready = 15 (every role at floor); default cap 18 is fine,
    # but --cap 9 trips the WIP limit with no floor breach.
    items = _spread("pm", 5, 100) + _spread("scientist", 5, 200) + _spread("developer", 5, 300)
    r = _run(items, "--cap", "9")
    assert "[CAP] Ready at 15 (>= 9)" in r.stdout, r.stdout
    assert "[REPLENISH" not in r.stdout and "[GROOMING-GAP" not in r.stdout, r.stdout
    assert r.returncode == 2, (r.returncode, r.stdout, r.stderr)


def test_invalid_floor_exits_1():
    r = _run([], "--floor", "abc")
    assert r.returncode == 1, (r.returncode, r.stdout, r.stderr)
    assert "floor must be" in r.stderr, r.stderr


def test_breakdown_shows_ready_inprogress_and_backlog():
    # The breakdown is self-documenting across all three relevant columns.
    items = (_spread("pm", 2, 100) + _spread("pm", 4, 150, status="Backlog")
             + _spread("pm", 1, 190, status="In progress")
             + _spread("scientist", 5, 200) + _spread("developer", 5, 300))
    r = _run(items)
    assert "Ready by role:" in r.stdout, r.stdout
    assert "In progress:" in r.stdout, r.stdout
    assert "Backlog candidates:" in r.stdout, r.stdout


def test_help_has_no_code_leak():
    r = subprocess.run(["bash", str(SCRIPT), "--help"], capture_output=True, text=True)
    assert r.returncode == 0, (r.returncode, r.stderr)
    assert "set -euo pipefail" not in r.stdout, r.stdout
    assert "SCRIPT_DIR=" not in r.stdout, r.stdout
