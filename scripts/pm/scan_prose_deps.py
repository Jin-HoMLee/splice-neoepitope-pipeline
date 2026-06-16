#!/usr/bin/env python3
"""Reconcile prose-stated Issue dependencies into native GitHub blockedBy edges.

Finds bodies saying "depends on #M" / "blocked on #M" etc., diffs them against
the native blockedBy graph, and reports or wires the missing edges.

Usage:
  scripts/pm/scan_prose_deps.py                 # --report (default): drift table
  scripts/pm/scan_prose_deps.py --issue N       # single-issue scope
  scripts/pm/scan_prose_deps.py --check         # exit 2 if any drift
  scripts/pm/scan_prose_deps.py --apply         # wire all needs-wiring records
  scripts/pm/scan_prose_deps.py --apply --only 745 594   # wire only these dependents

Exits 0 clean / 2 drift-present (--check) / 1 on error.
"""
import argparse
import json
import re
import subprocess
import sys

REPO = "Jin-HoMLee/splice-neoepitope-pipeline"


def gh(*args, parse_json=True):
    """Run a gh command; return parsed JSON (default) or raw stdout text."""
    result = subprocess.run(["gh", *args], capture_output=True, text=True, check=True)
    return json.loads(result.stdout) if parse_json else result.stdout


def fetch_open_issues():
    """All open issues with bodies. Uses the issue list (NOT the project board),
    so the board's Done-first pagination trap does not apply."""
    return gh(
        "issue", "list", "--repo", REPO, "--state", "open", "--limit", "1000",
        "--json", "number,title,body,state",
    )
