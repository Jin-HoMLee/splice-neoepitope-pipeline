"""Seed-narrative board-carrier forcing-function for the weekly SDR (Issue #1223).

The weekly meta-work SDR is milestone-free, so it has no milestone-close routing
anchor: a Carried-forward or actionable Retrospective item written into its sidecar
would otherwise dead-letter in a read-only HTML artifact. The window-mode seed must
therefore prompt the author to name a board carrier (new Issue / comment on an open
Issue / Discussion) OR mark the item "observation only, no carrier" for each such item.

Milestone mode already has its routing anchor (the milestone close), so it must NOT
carry the forcing prompt - the matched-pair control that proves the prompt is
mode-scoped, not unconditionally emitted.

Pure function - no network, no disk (role-less issues keep `_lab_notebook_seed` inert).
"""
import sys
from pathlib import Path

_PM_DIR = Path(__file__).resolve().parents[2] / "scripts" / "pm"
sys.path.insert(0, str(_PM_DIR))

import milestone_report as mr  # noqa: E402

_MILESTONE = {"title": "Meta-work SDR - week ending 2026-07-21"}
_ISSUES = [
    {"number": 900, "url": "https://example/900", "title": "A closed item",
     "state": "CLOSED", "stateReason": "COMPLETED", "roles": []},
]
_METRICS = {"per_role_counts": {}}


def test_window_mode_forces_a_carrier_decision_per_item():
    out = mr.seed_narrative(_MILESTONE, _ISSUES, _METRICS, mode="window")
    assert "board carrier" in out
    assert "observation only" in out


def test_milestone_mode_omits_the_carrier_forcing_prompt():
    out = mr.seed_narrative(_MILESTONE, _ISSUES, _METRICS, mode="milestone")
    assert "observation only" not in out
