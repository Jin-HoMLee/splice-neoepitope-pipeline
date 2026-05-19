#!/usr/bin/env python3
"""Recheck parent-vs-children Status drift on project board #9.

Implements the audit rule from docs/superpowers/specs/2026-05-19-parent-status-drift-audit-design.md.

Usage:
  scripts/pm/recheck_parent_status.py --issue N      # recheck this issue's parent chain
  scripts/pm/recheck_parent_status.py --all          # audit all parent issues on project #9

Exits 0 if no drift, 2 if drift detected, 1 on error.
"""
from __future__ import annotations

STATUS_LADDER = {
    "Backlog": 0,
    "Ready": 1,
    "In progress": 2,
    "Ready for review": 3,
    "In review": 4,
    "Done": 5,
}


def rank(status: str | None) -> int:
    """Map a Status name to its precedence rank. None/unknown → 0 (Backlog)."""
    return STATUS_LADDER.get(status or "", 0)


LADDER_INVERSE = {v: k for k, v in STATUS_LADDER.items()}


def collective_state(open_children: list[dict]) -> str:
    """Max-rank status across open children. Empty list → 'Done'."""
    if not open_children:
        return "Done"
    max_rank = max(rank(c.get("status")) for c in open_children)
    return LADDER_INVERSE[max_rank]


def classify_drift(parent_status: str | None, open_children: list[dict]) -> str | None:
    """Classify drift for a parent vs its open children.

    Returns one of: 'FORWARD DRIFT', 'BACKWARD DRIFT', 'COMPLETION DRIFT', or None.
    """
    p_rank = rank(parent_status)
    if not open_children:
        # All children closed; parent should be Done
        return None if p_rank == STATUS_LADDER["Done"] else "COMPLETION DRIFT"
    c_rank = rank(collective_state(open_children))
    if p_rank > c_rank:
        return "FORWARD DRIFT"
    if p_rank < c_rank:
        return "BACKWARD DRIFT"
    return None


import json
import subprocess

REPO = "Jin-HoMLee/splice-neoepitope-pipeline"
PROJECT_NUMBER = 9


def gh(*args: str, parse_json: bool = True):
    """Invoke `gh` and parse JSON output (or return raw text)."""
    result = subprocess.run(["gh", *args], capture_output=True, text=True, check=True)
    return json.loads(result.stdout) if parse_json else result.stdout


def parent_issue_number(issue_number: int) -> int | None:
    """Return parent issue number via REST parent_issue_url, or None."""
    data = gh("api", f"repos/{REPO}/issues/{issue_number}")
    url = data.get("parent_issue_url")
    if not url:
        return None
    return int(url.rstrip("/").rsplit("/", 1)[-1])


def open_sub_issues(issue_number: int) -> list[dict]:
    """Return list of open sub-issues (each a dict with 'number' at minimum)."""
    data = gh("api", f"repos/{REPO}/issues/{issue_number}/sub_issues")
    return [c for c in data if c.get("state") == "open"]
