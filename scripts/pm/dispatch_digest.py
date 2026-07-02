#!/usr/bin/env python3
"""Per-role board dispatch digest (Issue #721).

Rolls up, per role, the *committed* board work (Ready / In progress / Ready for
review / In review) plus open Team Coordination Discussions and aging WIP -
restoring the at-a-glance "what does each role have open/pending" visibility the
retired team_standup gave (#569). The board + Discussions already carry the data;
this just rolls it up into one glanceable digest.

Builds on scripts/board_open_items.py (the paginated ProjectV2 fetch + normalize)
and the open-Discussions query in shared/feedback_team_coordination.md.

Usage:
  scripts/pm/dispatch_digest.py                 # markdown digest to stdout
  scripts/pm/dispatch_digest.py --stale-days 21 # widen the aging-WIP threshold
  scripts/pm/dispatch_digest.py --json          # structured output (scheduled post)
  scripts/pm/dispatch_digest.py --no-discussions # skip the Discussions fetch

Invocation homes: a Warm-up / stand-up beat in the morning routine, or a scheduled
post (see docs/remote_routines.md). It is read-only - no board mutation.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

# board_open_items.py is a top-level module in scripts/; add scripts/ so the import
# resolves whether this runs as a CLI (scripts/pm/) or is imported under pytest.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import board_open_items as boi  # noqa: E402

# The Team Coordination Discussions category on the project repo (constant id -
# same one used in shared/feedback_team_coordination.md / the coordination skill).
DISCUSSION_CATEGORY_ID = "DIC_kwDORwn9EM4C-Jo6"
DISCUSSION_REPO_OWNER = "Jin-HoMLee"
DISCUSSION_REPO_NAME = "splice-neoepitope-pipeline"

# Committed statuses in display order (most-active first). Single-sourced from
# board_open_items.COMMITTED_STATUSES (the Backlog->Ready commitment boundary) so
# membership can never drift from the sibling; STATUS_ORDER supplies the render
# order. Resolves to ["In progress", "In review", "Ready for review", "Ready"].
COMMITTED_STATUSES = sorted(boi.COMMITTED_STATUSES, key=boi.STATUS_ORDER.get)

UNASSIGNED = "(unassigned)"

# Canonical role display order; roles not listed sort after these, alphabetically,
# with UNASSIGNED always last (a role-less committed item is a hygiene signal).
ROLE_ORDER = ["pm", "scientist", "developer", "memory_manager"]
ROLE_DISPLAY = {
    "pm": "PM",
    "scientist": "Scientist",
    "developer": "Developer",
    "memory_manager": "Memory Manager",
    UNASSIGNED: "(unassigned)",
}

DEFAULT_STALE_DAYS = 14


def _role_slug(item: dict[str, Any]) -> str:
    role = item.get("role")
    return role.removeprefix("role:") if role else UNASSIGNED


def group_by_role_status(
    items: list[dict[str, Any]],
) -> dict[str, dict[str, list[dict[str, Any]]]]:
    """Bucket committed items by role slug, then by status.

    Only the four committed statuses are kept; Backlog / No Status / Done are
    dropped (not active work). A role-less committed item lands under UNASSIGNED
    rather than being silently discarded - it is a hygiene signal worth surfacing.
    Insertion order within each bucket follows the input order.
    """
    grouped: dict[str, dict[str, list[dict[str, Any]]]] = {}
    for it in items:
        status = it.get("status")
        if status not in COMMITTED_STATUSES:
            continue
        role = _role_slug(it)
        grouped.setdefault(role, {}).setdefault(status, []).append(it)
    return grouped


def aging_items(
    items: list[dict[str, Any]],
    *,
    now: datetime,
    stale_days: int,
) -> list[dict[str, Any]]:
    """Committed items idle >= stale_days, most-dormant first.

    Keyed on content-level Issue.updatedAt (via board_open_items.age_days), so a
    bare board-field nudge does not reset staleness (Issue #642). Uncommitted
    statuses are excluded - only work someone has committed to can be stalled WIP.
    """
    aged = [
        it
        for it in items
        if it.get("status") in COMMITTED_STATUSES
        and (a := boi.age_days(it.get("updated_at"), now)) is not None
        and a >= stale_days
    ]
    aged.sort(key=lambda it: boi.age_days(it.get("updated_at"), now) or 0.0, reverse=True)
    return aged


def _ordered_roles(grouped: dict[str, Any]) -> list[str]:
    present = set(grouped)
    ordered = [r for r in ROLE_ORDER if r in present]
    extras = sorted(r for r in present if r not in ROLE_ORDER and r != UNASSIGNED)
    tail = [UNASSIGNED] if UNASSIGNED in present else []
    return ordered + extras + tail


def _item_link(it: dict[str, Any]) -> str:
    num = it.get("number")
    url = it.get("url")
    label = f"#{num}"
    linked = f"[{label}]({url})" if url else label
    title = (it.get("title") or "").strip()
    return f"{linked} {title}".rstrip()


def _parse_discussions_response(data: dict[str, Any]) -> list[dict[str, Any]]:
    """Extract discussion rows from a GraphQL response.

    Pure (no I/O) so the fail-soft edges are unit-testable. A null `nodes` list
    coalesces to [] and null node elements (GitHub returns these for inaccessible
    items) are skipped, so neither raises.
    """
    nodes = data["data"]["repository"]["discussions"]["nodes"] or []
    return [
        {
            "number": n.get("number"),
            "title": n.get("title", ""),
            "author": (n.get("author") or {}).get("login", "?"),
            "url": n.get("url"),
        }
        for n in nodes
        if n
    ]


def fetch_discussions() -> list[dict[str, Any]]:
    """Open Team Coordination Discussions as [{number, title, author, url}].

    Fails soft: on any gh/network/parse error, returns [] and warns to stderr, so
    a Discussions outage degrades the digest to board-only rather than aborting it.
    """
    query = """
    query($owner: String!, $name: String!, $cat: ID!) {
      repository(owner: $owner, name: $name) {
        discussions(first: 30, categoryId: $cat, states: [OPEN]) {
          nodes { number title url author { login } }
        }
      }
    }
    """
    cmd = [
        "gh", "api", "graphql",
        "-f", f"query={query}",
        "-F", f"owner={DISCUSSION_REPO_OWNER}",
        "-F", f"name={DISCUSSION_REPO_NAME}",
        "-F", f"cat={DISCUSSION_CATEGORY_ID}",
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"Warning: discussions fetch failed: {r.stderr.strip()}", file=sys.stderr)
            return []
        data = json.loads(r.stdout)
        if errs := data.get("errors"):
            print(f"Warning: discussions GraphQL errors: {errs}", file=sys.stderr)
            return []
        # Parse inside the guard so a malformed shape (null nodes, null element,
        # missing key) fails soft to a board-only digest rather than aborting.
        return _parse_discussions_response(data)
    except (json.JSONDecodeError, KeyError, TypeError, AttributeError) as e:
        print(f"Warning: could not parse discussions response: {e}", file=sys.stderr)
        return []


def render_digest(
    items: list[dict[str, Any]],
    *,
    discussions: list[dict[str, Any]],
    now: datetime,
    stale_days: int,
) -> str:
    """Render the digest as markdown. Pure given (items, discussions, now)."""
    grouped = group_by_role_status(items)
    aged = aging_items(items, now=now, stale_days=stale_days)

    lines: list[str] = []
    lines.append(f"# Board dispatch digest - {now.date().isoformat()}")
    lines.append("")
    lines.append(
        "_Per-role rollup of committed board work (Ready / In progress / "
        "Ready for review / In review), open Team Coordination Discussions, and "
        "aging WIP. Issue #721._"
    )
    lines.append("")

    if not grouped:
        lines.append("_No committed work on the board._")
        lines.append("")
    for role in _ordered_roles(grouped):
        statuses = grouped[role]
        counts = " · ".join(
            f"{s.lower()} {len(statuses.get(s, []))}" for s in COMMITTED_STATUSES
        )
        lines.append(f"## {ROLE_DISPLAY.get(role, role)}  ({counts})")
        for s in COMMITTED_STATUSES:
            bucket = statuses.get(s, [])
            if not bucket:
                continue
            rendered = " · ".join(_item_link(it) for it in bucket)
            lines.append(f"- **{s}:** {rendered}")
        lines.append("")

    lines.append(f"## Aging WIP (committed, idle >= {stale_days}d)")
    if aged:
        for it in aged:
            age = boi.age_label(it.get("updated_at"), now)
            lines.append(
                f"- {_item_link(it)} - {ROLE_DISPLAY.get(_role_slug(it), _role_slug(it))}"
                f" - {it.get('status')} - {age}"
            )
    else:
        lines.append("- (none)")
    lines.append("")

    lines.append(f"## Open Team Coordination Discussions ({len(discussions)})")
    if discussions:
        for d in discussions:
            num = d.get("number")
            url = d.get("url")
            label = f"[#{num}]({url})" if url else f"#{num}"
            lines.append(
                f"- {label} {(d.get('title') or '').strip()} (by {d.get('author', '?')})"
            )
    else:
        lines.append("- (none)")
    lines.append("")

    return "\n".join(lines) + "\n"


def build_digest_json(
    items: list[dict[str, Any]],
    *,
    discussions: list[dict[str, Any]],
    now: datetime,
    stale_days: int,
) -> dict[str, Any]:
    """Structured digest for --json (a scheduled post can render its own view)."""
    grouped = group_by_role_status(items)
    aged = aging_items(items, now=now, stale_days=stale_days)
    return {
        "generated_at": now.isoformat(),
        "stale_days": stale_days,
        "per_role": {
            role: {s: grouped[role].get(s, []) for s in COMMITTED_STATUSES}
            for role in _ordered_roles(grouped)
        },
        "aging_wip": aged,
        "discussions": discussions,
    }


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--stale-days", dest="stale_days", type=int, default=DEFAULT_STALE_DAYS,
        metavar="N", help=f"Aging-WIP threshold in days (default {DEFAULT_STALE_DAYS})",
    )
    p.add_argument(
        "--no-discussions", dest="no_discussions", action="store_true",
        help="Skip the Team Coordination Discussions fetch (board-only digest)",
    )
    p.add_argument(
        "--json", dest="as_json", action="store_true",
        help="Emit a structured JSON digest instead of markdown",
    )
    args = p.parse_args()

    if args.stale_days < 0:
        p.error("--stale-days must be a non-negative integer")

    now = datetime.now(timezone.utc)
    raw = boi.fetch_all_items()
    items = [n for it in raw if (n := boi.normalize(it)) is not None]
    discussions = [] if args.no_discussions else fetch_discussions()

    if args.as_json:
        json.dump(
            build_digest_json(items, discussions=discussions, now=now, stale_days=args.stale_days),
            sys.stdout, indent=2,
        )
        sys.stdout.write("\n")
    else:
        sys.stdout.write(
            render_digest(items, discussions=discussions, now=now, stale_days=args.stale_days)
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
