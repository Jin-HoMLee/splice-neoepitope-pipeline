#!/usr/bin/env python3
"""Blocking pre-merge cross-repo AC gate for the closure ritual (Issue #665).

The same-repo Acceptance-criteria check in scripts/audit_and_merge.sh keys off
GitHub's native `closingIssuesReferences`, which is only ever populated within a
single repo. So a cross-repo close — the canonical case being a personas-repo PR
that closes a project-repo Issue (the Memory Manager's workflow) — slips the AC
gate entirely: the project Issue's Acceptance criteria are never checked before
the personas PR merges.

This gate closes that gap. It parses the PR body for cross-repo closing
forward-links (`Closes owner/repo#N`, a full issue URL, or the project's
`[Issue #N](url) (closes)` link+keyword form) and audits each target Issue's
Acceptance criteria from *its own* repo. It exits 1 (blocking) on any unticked
AC box, mirroring the same-repo AC check's intent.

Single-sources the parsing + audit in closure_audit (collect_cross_repo_ac_gaps)
so the logic is unit-tested once. Honors the REPO override so it composes with
the sibling gates (stray_closers / lab_notebook_gate / ac_section_lint) when
audit_and_merge.sh runs against a non-default repo (#607); the PR's own repo is
excluded as same-repo. Fails OPEN (exit 0 + warning) on a gh/JSON error, like
those siblings — a hiccup must never block a legitimate merge.

Exit codes:
    0 — no cross-repo closing target has an unticked AC box (or fail-open)
    1 — a cross-repo closing target has unticked Acceptance-criteria boxes
    2 — usage error

Usage:
    python cross_repo_ac_gate.py <PR_NUMBER>
"""

from __future__ import annotations

import json
import os
import subprocess
import sys

import closure_audit as ca


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: cross_repo_ac_gate.py <PR_NUMBER>", file=sys.stderr)
        return 2
    try:
        n = int(sys.argv[1])
    except ValueError:
        print("usage: cross_repo_ac_gate.py <PR_NUMBER>", file=sys.stderr)
        return 2

    repo = os.environ.get("REPO", "Jin-HoMLee/splice-neoepitope-pipeline")
    try:
        gaps = ca.collect_cross_repo_ac_gaps(n, repo=repo)
    except (subprocess.CalledProcessError, json.JSONDecodeError, KeyError, OSError) as e:
        print(f"⚠ cross-repo AC check skipped (gh error: {e})", file=sys.stderr)
        return 0  # fail open — a hiccup must not block a merge

    if not gaps:
        return 0
    for target, msg in gaps:
        print(f"✗ cross-repo close {target} has unticked Acceptance criteria: {msg}",
              file=sys.stderr)
    print("    Tick the cross-repo Issue's AC boxes (or break the closing-keyword↔link "
          "adjacency if it isn't really a close) before merging.", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
