#!/usr/bin/env python3
"""List open items on the JH M Lee Lab project board (#9) with Status/Priority/Size.

Wraps the paginated ProjectV2 GraphQL query that the board needs:
  - GraphQL `first: N` is a cap, not a filter; we loop on `pageInfo.hasNextPage`.
  - The board sorts Done items first, so unpaginated single-page queries silently
    truncate open work. See `.claude/memory/shared/feedback_board_queries.md`.

For queries that only need label + state (no Status/Priority/Size columns),
prefer `gh issue list --label role:<role> --state open` instead — it's a one-liner
with no pagination needed.

Usage:
  scripts/board_open_items.py
  scripts/board_open_items.py --role developer
  scripts/board_open_items.py --status Ready --priority P1
  scripts/board_open_items.py --json | jq '.[] | select(.size == "S")'
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from typing import Any

OWNER = "Jin-HoMLee"
PROJECT_NUMBER = 9

QUERY = """
query($owner: String!, $number: Int!, $after: String) {
  user(login: $owner) {
    projectV2(number: $number) {
      items(first: 100, after: $after) {
        pageInfo { hasNextPage endCursor }
        nodes {
          content {
            __typename
            ... on Issue {
              number title state url
              labels(first: 20) { nodes { name } }
            }
            ... on PullRequest {
              number title state url isDraft
              labels(first: 20) { nodes { name } }
            }
          }
          fieldValues(first: 20) {
            nodes {
              ... on ProjectV2ItemFieldSingleSelectValue {
                name
                field { ... on ProjectV2SingleSelectField { name } }
              }
            }
          }
        }
      }
    }
  }
}
"""

STATUS_ORDER = {
    "In progress": 0,
    "In review": 1,
    "Ready for review": 2,
    "Ready": 3,
    "Backlog": 4,
    "No Status": 5,
}
PRIORITY_ORDER = {"P0": 0, "P1": 1, "P2": 2, "P3": 3}
SIZE_ORDER = {"XS": 0, "S": 1, "M": 2, "L": 3, "XL": 4}


def fetch_all_items() -> list[dict[str, Any]]:
    items: list[dict[str, Any]] = []
    cursor: str | None = None
    page = 0
    while True:
        page += 1
        cmd = [
            "gh", "api", "graphql",
            "-f", f"query={QUERY}",
            "-F", f"owner={OWNER}",
            "-F", f"number={PROJECT_NUMBER}",
        ]
        if cursor is not None:
            cmd.extend(["-F", f"after={cursor}"])
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stderr, file=sys.stderr)
            sys.exit(r.returncode)
        data = json.loads(r.stdout)
        if errs := data.get("errors"):
            print(f"GraphQL errors: {errs}", file=sys.stderr)
            sys.exit(1)
        proj = data["data"]["user"]["projectV2"]["items"]
        items.extend(proj["nodes"])
        pi = proj["pageInfo"]
        print(
            f"Page {page}: {len(proj['nodes'])} items, hasNext={pi['hasNextPage']}",
            file=sys.stderr,
        )
        if not pi["hasNextPage"]:
            break
        cursor = pi["endCursor"]
    return items


def normalize(item: dict[str, Any]) -> dict[str, Any] | None:
    content = item.get("content") or {}
    if not content:
        return None
    state = content.get("state")
    if state in ("CLOSED", "MERGED"):
        return None
    status = size = priority = None
    for fv in item["fieldValues"]["nodes"]:
        if not fv:
            continue
        field = fv.get("field", {}) or {}
        fname = field.get("name")
        if fname == "Status":
            status = fv.get("name")
        elif fname == "Size":
            size = fv.get("name")
        elif fname == "Priority":
            priority = fv.get("name")
    if status == "Done":
        return None
    labels = [l["name"] for l in content.get("labels", {}).get("nodes", [])]
    role = next((l for l in labels if l.startswith("role:")), None)
    return {
        "number": content.get("number"),
        "title": content.get("title", ""),
        "url": content.get("url", ""),
        "kind": "PR" if content.get("__typename") == "PullRequest" else "Issue",
        "is_draft": content.get("isDraft", False),
        "state": state,
        "status": status or "No Status",
        "priority": priority,
        "size": size,
        "role": role,
        "labels": labels,
    }


def sort_key(it: dict[str, Any]) -> tuple:
    return (
        STATUS_ORDER.get(it["status"], 99),
        PRIORITY_ORDER.get(it["priority"], 99),
        SIZE_ORDER.get(it["size"], 99),
        it["number"] or 0,
    )


def matches_filter(it: dict[str, Any], args: argparse.Namespace) -> bool:
    if args.role:
        want = args.role if args.role.startswith("role:") else f"role:{args.role}"
        if it["role"] != want:
            return False
    if args.status and it["status"] != args.status:
        return False
    if args.priority and it["priority"] != args.priority:
        return False
    if args.size and it["size"] != args.size:
        return False
    return True


def format_table(items: list[dict[str, Any]]) -> str:
    if not items:
        return "(no items matched)\n"
    lines = [
        f"{'Status':<17} {'P':<3} {'Sz':<3} {'Role':<17} {'Kind':<5} {'#':<5} Title",
        f"{'-' * 17} {'-' * 3} {'-' * 3} {'-' * 17} {'-' * 5} {'-' * 5} {'-' * 60}",
    ]
    for it in items:
        role = (it["role"] or "(none)").removeprefix("role:")[:16]
        kind = it["kind"] + ("/D" if it["is_draft"] else "")
        title = (it["title"] or "")[:60]
        lines.append(
            f"{it['status']:<17} {it['priority'] or '?':<3} {it['size'] or '?':<3} "
            f"{role:<17} {kind:<5} {it['number']:<5} {title}"
        )
    return "\n".join(lines) + "\n"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--role", help='Filter by role label (e.g. "developer", "scientist", "pm")')
    p.add_argument("--status", choices=list(STATUS_ORDER), help="Filter by board Status")
    p.add_argument("--priority", choices=list(PRIORITY_ORDER), help="Filter by Priority")
    p.add_argument("--size", choices=list(SIZE_ORDER), help="Filter by Size")
    p.add_argument("--json", dest="as_json", action="store_true", help="Emit JSON array instead of a table")
    args = p.parse_args()

    raw = fetch_all_items()
    normalized = [n for it in raw if (n := normalize(it)) is not None]
    filtered = [it for it in normalized if matches_filter(it, args)]
    filtered.sort(key=sort_key)

    print(
        f"Total board items: {len(raw)}; open: {len(normalized)}; "
        f"after filters: {len(filtered)}",
        file=sys.stderr,
    )

    if args.as_json:
        json.dump(filtered, sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        sys.stdout.write(format_table(filtered))
    return 0


if __name__ == "__main__":
    sys.exit(main())
