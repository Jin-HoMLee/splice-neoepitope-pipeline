#!/usr/bin/env python3
"""Pre-merge lab-notebook-entry gate (Issue #409).

Thin wrapper around closure_audit.audit_pr_pre_merge. The post-hoc closure-audit
bot (closure_audit.py, run by .github/workflows/closure-audit.yml) flags a
missing lab-notebook entry only AFTER merge, as a marker comment — a cleanup
loop, not a prevention. This gate moves the same check to merge time so a
role-tagged PR cannot merge without its `## <today>` lab-notebook entry. Called
by scripts/audit_and_merge.sh as part of the closure-ritual gate.

Exit codes:
    0 — clean, exempt (only lab_notebook/glossary files touched), skip-marked
        (`<!-- skip-lab-notebook: routine -->` in the PR body), OR a gh/network
        error (fails OPEN with a warning, mirroring stray_closers.py /
        bot_review_offer.py — a hiccup must not block a legitimate merge)
    1 — a role-tagged linked Issue lacks a `## <today>` entry in
        research/lab_notebook/<role>.md referencing the PR/Issue (gaps printed)
    2 — usage error

The entry is read from the working tree, so run this from the PR branch where
the entry was written (the documented closure-ritual flow: write the entry,
then merge).

Usage:
    python lab_notebook_gate.py <PR_NUMBER>
"""

from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime, timezone

import closure_audit as ca


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: lab_notebook_gate.py <PR_NUMBER>", file=sys.stderr)
        return 2
    try:
        n = int(sys.argv[1])
    except ValueError:
        print("usage: lab_notebook_gate.py <PR_NUMBER>", file=sys.stderr)
        return 2

    today = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    try:
        gaps = ca.audit_pr_pre_merge(n, today)
    except (subprocess.CalledProcessError, json.JSONDecodeError, KeyError, OSError) as e:
        print(f"⚠ lab-notebook check skipped (gh error: {e})", file=sys.stderr)
        return 0  # fail open — see module docstring

    if gaps:
        print(f"✗ PR #{n} is missing a lab-notebook entry for {today}:", file=sys.stderr)
        for role, desc in gaps:
            print(f"    {role}: {desc}", file=sys.stderr)
        print(
            f"  Add a '## {today}' block in research/lab_notebook/<role>.md with a "
            f"'### HH:MM UTC — Editor: …' sub-section referencing #{n} "
            f"(or a closing Issue #).",
            file=sys.stderr,
        )
        print(
            "  Bypass: add the entry, or put '<!-- skip-lab-notebook: routine -->' "
            "in the PR body for a routine single-Issue ship.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
