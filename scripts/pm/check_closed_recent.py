#!/usr/bin/env python3
"""Recap of recently closed Issues + merged PRs (morning-routine Beat 1a/1b).

The Service Delivery Review closure-recap was a hand-typed query each morning:
`gh issue list --state closed --search "closed:>YYYY-MM-DD"`. The `>` operator
excludes the whole boundary day, so `closed:>2026-06-18` returns only items closed
06-19 onward - silently reporting "nothing closed" when items merged *during* the
boundary day. That off-by-one slipped 2026-06-19 (caught by the user). This script
makes the boundary deterministic (`>=`) and covers both closed Issues AND merged
PRs in one place - the closed-side companion to `board_open_items.py`'s open scan.

Usage:

  scripts/pm/check_closed_recent.py                 # default: since the start of yesterday (UTC)
  scripts/pm/check_closed_recent.py --days 3        # since the start of the day 3 days ago
  scripts/pm/check_closed_recent.py --since 2026-07-01

Window floor is a whole day (00:00 UTC): the recap is day-granular, matching how
GitHub's `closed:` search operator works. `--since YYYY-MM-DD` / `--days N`
override the default.

Exit 0 on a completed scan - finding nothing closed is not a failure. A hard `gh`
/ network error on a listing call surfaces loudly (non-zero traceback) rather than
as a false "nothing closed". Issue #784.
"""
import argparse
import json
import subprocess
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

REPO = "Jin-HoMLee/splice-neoepitope-pipeline"
UTC = timezone.utc


# --- pure helpers (unit-tested, no I/O) ---


def compute_floor(now, days=1):
    """Window floor as an aware UTC datetime at 00:00 of the day `days` days ago.

    `days=1` (default) -> the start of yesterday, so the recap includes everything
    closed yesterday and today. The floor is truncated to midnight because the
    `closed:` search operator is day-granular; a sub-day floor would imply a
    precision the query cannot honor.
    """
    if days < 0:
        raise ValueError(f"days must be non-negative, got {days}")
    floor_date = (now.astimezone(UTC) - timedelta(days=days)).date()
    return datetime(floor_date.year, floor_date.month, floor_date.day, tzinfo=UTC)


def floor_date_str(floor):
    """The `YYYY-MM-DD` string used in the search operator."""
    return floor.astimezone(UTC).strftime("%Y-%m-%d")


def build_search(floor):
    """The GitHub search fragment - `>=` (inclusive), never bare `>`.

    Using `>=` includes items closed *during* the boundary day; the bare `>` this
    replaces excluded the whole boundary day, which is the exact bug (Issue #784).
    """
    return f"closed:>={floor_date_str(floor)}"


def parse_ts(ts):
    """Parse a GitHub ISO-8601 timestamp to an aware UTC datetime (or None)."""
    if not ts:
        return None
    try:
        dt = datetime.fromisoformat(ts.replace("Z", "+00:00"))
        return dt if dt.tzinfo else dt.replace(tzinfo=UTC)
    except (ValueError, TypeError):
        return None


def is_within(when_ts, floor):
    """True if `when_ts` (ISO-8601) is at/after `floor`.

    A client-side belt-and-suspenders on top of the `>=` search: pins that an item
    closed at the very start of the boundary day (00:00:00Z) is *included*, so the
    off-by-one cannot silently return. An absent/unparseable timestamp is kept
    (fail-open: better a stray row than a dropped closure).
    """
    dt = parse_ts(when_ts)
    return True if dt is None else dt >= floor


# --- gh I/O ---


def gh_json(*args):
    result = subprocess.run(["gh", *args], capture_output=True, text=True, check=True)
    return json.loads(result.stdout)


def list_closed_issues(search):
    """Closed Issues matching `search` (number, title, closedAt)."""
    return gh_json(
        "issue", "list", "--repo", REPO, "--state", "closed",
        "--search", search, "--limit", "1000",
        "--json", "number,title,closedAt",
    )


def list_merged_prs(search):
    """PRs matching `search`, filtered to the *merged* ones (mergedAt present).

    `--state all --search "closed:>=…"` returns both merged and closed-unmerged
    PRs; only merged PRs belong in a shipped-work recap, so the unmerged ones
    (mergedAt is null) are dropped.
    """
    rows = gh_json(
        "pr", "list", "--repo", REPO, "--state", "all",
        "--search", search, "--limit", "1000",
        "--json", "number,title,mergedAt",
    )
    return [r for r in rows if r.get("mergedAt")]


# --- orchestration ---


def collect(floor, issues, prs):
    """Merge closed Issues + merged PRs into recap rows, newest first.

    `issues` / `prs` are the raw `gh` result lists (injected, so this is testable
    without network). Each is re-filtered client-side via `is_within` so the
    boundary is enforced here too, not only by the search string.
    """
    rows = []
    for it in issues:
        if is_within(it.get("closedAt"), floor):
            rows.append({"number": it["number"], "kind": "Issue",
                         "when": it.get("closedAt") or "?", "title": it.get("title") or ""})
    for pr in prs:
        if is_within(pr.get("mergedAt"), floor):
            rows.append({"number": pr["number"], "kind": "PR",
                         "when": pr.get("mergedAt") or "?", "title": pr.get("title") or ""})
    rows.sort(key=lambda r: r["when"], reverse=True)
    return rows


def render(rows, floor):
    header = f"Closed Issues + merged PRs since {floor_date_str(floor)} (00:00 UTC):"
    if not rows:
        return header + "\n  (nothing closed in window)\n"
    lines = [header]
    for r in rows:
        when = r["when"][:16].replace("T", " ")
        lines.append(f"  #{r['number']:<5} {r['kind']:<5} {when}  {r['title']}")
    lines.append(f"\n  {len(rows)} item(s).")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--since", help="window floor YYYY-MM-DD (overrides --days)")
    group.add_argument("--days", type=int, default=1,
                       help="window = since the start of the day N days ago (default 1 = yesterday)")
    args = parser.parse_args()

    now = datetime.now(UTC)
    if args.since:
        floor = parse_ts(args.since + "T00:00:00Z")
        if floor is None:
            parser.error(f"unparseable --since (want YYYY-MM-DD): {args.since!r}")
    else:
        if args.days < 0:
            parser.error(f"--days must be non-negative, got {args.days}")
        floor = compute_floor(now, args.days)

    search = build_search(floor)
    rows = collect(floor, list_closed_issues(search), list_merged_prs(search))
    print(render(rows, floor), end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
