#!/usr/bin/env python3
"""Surface open board items whose project Target date is past (roadmap-overdue sweep).

For the PM morning-routine Service Delivery Review (Phase 1d). The third
`check_*_health` sibling, alongside `check_milestone_health.sh` (milestone
`due_on` clock) and `check_ready_queue.sh` (Ready floor). This one watches the
per-issue **Target date** ProjectV2 field, which can desync from `due_on` (e.g. a
milestone change not paired with a Target re-sync, or an un-milestoned issue
still carrying a stale Target) — so it's the detection backstop for that desync,
NOT a duplicate of the milestone clock.

Language-by-data-source (why `.py`, not `.sh`): the `check_*` family splits by
data source. Plain GitHub REST objects (milestones, issues) → shell + jq (the
two `.sh` checks). ProjectV2 item fields (Status/Priority/Size/Target),
paginated across the whole board → Python, reusing `board_open_items.py`'s
paginator + Target-date extraction rather than re-hand-rolling the cursor loop
and the `fieldValues` union unpacking. See Issue #704.

Usage:
  scripts/check_roadmap_health.py
  scripts/check_roadmap_health.py --role pm
  scripts/check_roadmap_health.py --json | jq '.[].number'

Exit codes:
  0 — all clear (no open item with a past Target date)
  2 — at least one open item is roadmap-overdue
  1 — usage / runtime error
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date, datetime, timezone
from typing import Any

import board_open_items as boi


def _parse_target(s: str | None) -> date | None:
    """Parse a ProjectV2 date field value ('YYYY-MM-DD'); None when absent/blank."""
    if not s:
        return None
    return datetime.strptime(s[:10], "%Y-%m-%d").date()


def find_overdue(items: list[dict[str, Any]], today: date) -> list[dict[str, Any]]:
    """Items whose Target date is strictly before `today`, most-overdue first.

    Items without a Target date are excluded; `today` itself is the deadline, not
    overdue (matches the AC's "Target date < today").
    """
    overdue = []
    for it in items:
        td = _parse_target(it.get("target_date"))
        if td is not None and td < today:
            overdue.append((td, it))
    overdue.sort(key=lambda pair: pair[0])  # oldest Target (most overdue) first
    return [it for _, it in overdue]


def _roles_of(it: dict[str, Any]) -> list[str]:
    """Every `role:` label on an item, tolerating a legacy first-role-only dict.

    Falls back to `[role]` when `roles` is absent rather than treating the item as
    role-less. Reading `roles` alone would silently bucket every legacy-shaped item
    into '(none)' - which is the same silent-wrong-answer class #1153 is about, so
    the fix must not reintroduce it one level down.
    """
    roles = it.get("roles")
    if roles:
        return list(roles)
    role = it.get("role")
    return [role] if role else []


def matches_role(it: dict[str, Any], want: str) -> bool:
    """True if `want` is ANY of the item's role: labels, not merely the first.

    Membership, not equality (#1153). An Issue may carry several `role:` labels
    (the dual-role convention: one role Leads, the other owns scoped ACs), and
    GitHub returns labels in an unstable order - so matching `it["role"]`
    (= labels[0]) drops a dual-role item from one of its owners' views
    nondeterministically. Mirrors `board_open_items.matches_filter`.
    """
    want = want if want.startswith("role:") else f"role:{want}"
    return want in _roles_of(it)


def group_by_role(items: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    """Bucket items by role:* label; a missing role goes to '(none)'.

    An item with two roles appears under BOTH headings (#1153) - an overdue item
    owned jointly by Sci and Dev is overdue for both of them, and bucketing it
    under whichever label GitHub happened to return first means one of its owners
    never sees it in this report.
    """
    groups: dict[str, list[dict[str, Any]]] = {}
    for it in items:
        for role in _roles_of(it) or ["(none)"]:
            groups.setdefault(role, []).append(it)
    return groups


def format_report(overdue: list[dict[str, Any]], today: date) -> str:
    if not overdue:
        return "Roadmap health: all clear — no open item past its Target date.\n"
    lines = [f"Roadmap-overdue (Target date < {today.isoformat()}): {len(overdue)} item(s)\n"]
    for role, items in sorted(group_by_role(overdue).items()):
        lines.append(f"{role} ({len(items)}):")
        for it in items:
            td = _parse_target(it.get("target_date"))
            days = (today - td).days if td else "?"
            lines.append(
                f"  #{it['number']:<5} {td.isoformat() if td else '—'} "
                f"({days}d overdue)  {(it.get('title') or '')[:60]}"
            )
        lines.append("")
    return "\n".join(lines)


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--role", help='Filter by role label (e.g. "pm", "developer", "scientist")')
    p.add_argument("--json", dest="as_json", action="store_true",
                   help="Emit a flat JSON array of overdue items instead of a report")
    args = p.parse_args()

    today = datetime.now(timezone.utc).date()
    try:
        raw = boi.fetch_all_items()
    except SystemExit:
        raise
    except Exception as e:  # network / gh / parse error → runtime error
        print(f"error fetching board: {e}", file=sys.stderr)
        return 1
    normalized = [n for it in raw if (n := boi.normalize(it)) is not None]

    if args.role:
        normalized = [it for it in normalized if matches_role(it, args.role)]

    overdue = find_overdue(normalized, today)

    if args.as_json:
        json.dump(overdue, sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        sys.stdout.write(format_report(overdue, today))

    return 2 if overdue else 0


if __name__ == "__main__":
    sys.exit(main())
