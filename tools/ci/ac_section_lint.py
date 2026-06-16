#!/usr/bin/env python3
"""Non-blocking stray-AC-box lint for the closure-ritual gate (Issue #730).

The closure gate (scripts/audit_and_merge.sh) blocks a merge on unticked boxes
*under* a `## Acceptance criteria` section. But an Issue whose gating boxes live
under a different heading (e.g. `## Plan (phased)`, as Issue #569 did) has *zero*
boxes under "Acceptance criteria", so the blocking AC check sees nothing and
passes silently — the gate's AC enforcement is bypassed by a heading-name
deviation (PR #724).

This lint closes that silent-pass by WARNING (never blocking) when a linked
Issue has unticked `- [ ]` boxes but no `## Acceptance criteria` section, naming
the count + the non-AC heading(s) so the operator can confirm none are gating.
The fix is a convention, not a parser change: gating boxes belong under the one
canonical heading (broadening the parser to read Plan/Tasks/Checklist would
false-block the non-gating boxes those sections routinely carry).

Always exits 0 — this is advisory only. Called by scripts/audit_and_merge.sh
after the blocking AC checks. Fails OPEN (exit 0 + warning) on a gh/JSON error,
mirroring stray_closers.py / lab_notebook_gate.py.

Exit codes:
    0 — always, unless usage is wrong (non-blocking lint; warnings go to stderr)
    2 — usage error

Usage:
    python ac_section_lint.py <PR_NUMBER>
"""

from __future__ import annotations

import json
import os
import subprocess
import sys

import closure_audit as ca


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: ac_section_lint.py <PR_NUMBER>", file=sys.stderr)
        return 2
    try:
        n = int(sys.argv[1])
    except ValueError:
        print("usage: ac_section_lint.py <PR_NUMBER>", file=sys.stderr)
        return 2

    # Honor the REPO override so this gate composes with its sibling gates
    # (stray_closers / lab_notebook_gate) when audit_and_merge.sh runs against a
    # fork (Issue #607). Default matches them — the canonical repo.
    repo = os.environ.get("REPO", "Jin-HoMLee/splice-neoepitope-pipeline")
    try:
        warnings = ca.collect_stray_ac_warnings(n, repo=repo)
    except (subprocess.CalledProcessError, json.JSONDecodeError, KeyError, OSError) as e:
        print(f"⚠ stray-AC-box lint skipped (gh error: {e})", file=sys.stderr)
        return 0  # fail open — advisory only, never block a merge

    for num, msg in warnings:
        print(f"⚠ Issue #{num}: {msg}", file=sys.stderr)
    return 0  # non-blocking: warnings never fail the merge


if __name__ == "__main__":
    sys.exit(main())
