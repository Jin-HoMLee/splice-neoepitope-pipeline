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


# Parents/epics are PARKED in this off-ladder Status and carry neither a milestone
# nor a Target date by convention (Issue #690 sub-question A). A Target on one is
# leftover data, not a commitment.
EPIC_STATUS = "Epic"


def _is_epic_parked(it: dict[str, Any]) -> bool:
    return (it.get("status") or "") == EPIC_STATUS


def find_overdue(items: list[dict[str, Any]], today: date) -> list[dict[str, Any]]:
    """Items whose Target date is strictly before `today`, most-overdue first.

    Items without a Target date are excluded; `today` itself is the deadline, not
    overdue (matches the AC's "Target date < today").

    Epic-parked parents are excluded (Issue #1149). A parked epic carries no Target
    by convention, so a leftover date on one is a DATA-HYGIENE problem ("clear the
    date"), not a missed DEADLINE ("ship it") - reporting it as overdue frames the
    finding wrongly and puts a permanently-noisy line in a daily sweep. They are not
    dropped: find_parent_targets() surfaces them under their correct frame.
    """
    overdue = []
    for it in items:
        if _is_epic_parked(it):
            continue
        td = _parse_target(it.get("target_date"))
        if td is not None and td < today:
            overdue.append((td, it))
    overdue.sort(key=lambda pair: pair[0])  # oldest Target (most overdue) first
    return [it for _, it in overdue]


def find_parent_targets(items: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Epic-parked parents that carry a Target date at all (Issue #1149).

    Not date-compared: ANY Target on a parked epic is the finding, past or future,
    because the convention is that a parent carries none. The remedy is to clear it.
    """
    return [
        it
        for it in items
        if _is_epic_parked(it) and _parse_target(it.get("target_date")) is not None
    ]


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


def _format_parent_targets(parent_targets: list[dict[str, Any]]) -> list[str]:
    """The [PARENT-TARGET] hygiene block (Issue #1149).

    A separate frame on purpose: the remedy is to CLEAR the date, not to ship the
    work, so folding these into the overdue list would prescribe the wrong action.
    """
    lines = [
        f"[PARENT-TARGET] {len(parent_targets)} Epic-parked parent(s) carry a Target date,",
        "which the parent convention forbids (Issue #690-A). Remedy: clear the date.",
        "",
    ]
    for it in parent_targets:
        td = _parse_target(it.get("target_date"))
        lines.append(
            f"  #{it['number']:<5} Target={td.isoformat() if td else '?'}"
            f"  {(it.get('title') or '')[:60]}"
        )
    lines.append("")
    return lines


def format_report(
    overdue: list[dict[str, Any]],
    today: date,
    parent_targets: list[dict[str, Any]] | None = None,
) -> str:
    parent_targets = parent_targets or []
    if not overdue and not parent_targets:
        return "Roadmap health: all clear — no open item past its Target date.\n"

    lines: list[str] = []
    if overdue:
        lines.append(
            f"Roadmap-overdue (Target date < {today.isoformat()}): {len(overdue)} item(s)\n"
        )
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
    if parent_targets:
        lines.extend(_format_parent_targets(parent_targets))
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
    parent_targets = find_parent_targets(normalized)

    if args.as_json:
        json.dump(overdue, sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        sys.stdout.write(format_report(overdue, today, parent_targets))

    return 2 if overdue else 0


if __name__ == "__main__":
    sys.exit(main())
