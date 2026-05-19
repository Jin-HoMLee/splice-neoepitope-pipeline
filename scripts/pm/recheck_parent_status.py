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
