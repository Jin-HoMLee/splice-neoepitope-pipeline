#!/usr/bin/env python3
"""Recheck a milestone's due_on against its remaining size-weighted capacity.

Implements the recheck rule from memory/feedback_milestones.md:
  new_due_on = today + remaining-days / 5.0 * 7 days

Usage:
  scripts/pm/recheck_milestone.py --issue N
  scripts/pm/recheck_milestone.py --milestone N

Exits 0 if delta within +-7 days, 2 if UPDATE NEEDED, 1 on error.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import date, datetime, timedelta

REPO = "Jin-HoMLee/splice-neoepitope-pipeline"
PROJECT_NUMBER = 9
THRESHOLD_DAYS = 7
AVAILABILITY_RATE = 5.0  # capacity-days per calendar-week

# Size weight in capacity-days (midpoints of feedback_milestones.md ranges)
SIZE_WEIGHTS = {"XS": 0.5, "S": 1.0, "M": 2.5, "L": 3.5, "XL": 5.0}


def parse_milestone_title(title: str) -> tuple[int, int] | None:
	"""Parse 'i<N> - S<M> - ...' titles. Returns (iteration, stage) or None.

	Role-meta (pm-i*, dev-i*) and legacy (M1, M2) titles return None and fall
	through to pure-capacity behavior.
	"""
	m = re.match(r"^i(\d+)\s*-\s*S(\d+)\s*-\s*", title)
	return (int(m.group(1)), int(m.group(2))) if m else None


def gh(*args: str, parse_json: bool = True) -> object:
    result = subprocess.run(["gh", *args], capture_output=True, text=True, check=True)
    return json.loads(result.stdout) if parse_json else result.stdout


def milestone_for_issue(issue_number: int) -> int | None:
    data = gh("issue", "view", str(issue_number), "--repo", REPO, "--json", "milestone")
    return data["milestone"]["number"] if data.get("milestone") else None


def milestone_meta(milestone_number: int) -> dict:
    return gh("api", f"repos/{REPO}/milestones/{milestone_number}")


def open_issues_in_milestone(milestone_title: str) -> list[int]:
    data = gh(
        "issue", "list",
        "--repo", REPO,
        "--milestone", milestone_title,
        "--state", "open",
        "--limit", "100",
        "--json", "number",
    )
    return [issue["number"] for issue in data]


def sizes_for_issues(issue_numbers: list[int]) -> dict[int, str | None]:
    if not issue_numbers:
        return {}
    owner, name = REPO.split("/")
    aliases = " ".join(
        f'i{n}: issue(number: {n}) {{ '
        f'projectItems(first: 5) {{ nodes {{ project {{ number }} '
        f'fieldValues(first: 20) {{ nodes {{ '
        f'... on ProjectV2ItemFieldSingleSelectValue {{ name field {{ ... on ProjectV2SingleSelectField {{ name }} }} }} '
        f'}} }} }} }} }}'
        for n in issue_numbers
    )
    query = f'query {{ repository(owner: "{owner}", name: "{name}") {{ {aliases} }} }}'
    data = gh("api", "graphql", "-f", f"query={query}")
    repo = data["data"]["repository"]
    sizes: dict[int, str | None] = {}
    for n in issue_numbers:
        node = repo.get(f"i{n}") or {}
        size: str | None = None
        for pi in node.get("projectItems", {}).get("nodes", []):
            if (pi.get("project") or {}).get("number") != PROJECT_NUMBER:
                continue
            for fv in pi.get("fieldValues", {}).get("nodes", []):
                if (fv.get("field") or {}).get("name") == "Size":
                    size = fv.get("name")
        sizes[n] = size
    return sizes


def compute_recheck(milestone_number: int) -> int:
    meta = milestone_meta(milestone_number)
    title = meta["title"]
    current_due_raw = meta.get("due_on")
    current_due_date = (
        datetime.fromisoformat(current_due_raw.replace("Z", "+00:00")).date()
        if current_due_raw else None
    )

    issue_numbers = open_issues_in_milestone(title)
    sizes = sizes_for_issues(issue_numbers)

    print(f"Milestone: {title}")
    print(f"Current due_on: {current_due_date or '—'}")
    print(f"Open issues ({len(issue_numbers)}):")
    remaining = 0.0
    for n in sorted(issue_numbers):
        size = sizes.get(n)
        weight = SIZE_WEIGHTS.get(size or "", 0)
        remaining += weight
        size_disp = f"{size}, ~{weight}d" if size else "no size, 0d"
        print(f"  - #{n} ({size_disp})")
    print(f"Remaining capacity: {remaining}d")

    if remaining == 0:
        unsized_count = sum(1 for n in issue_numbers if sizes.get(n) is None)
        if unsized_count > 0:
            print(f"Proposed due_on: (cannot compute — {unsized_count} open issue(s) missing Size)")
            print("Status: [UNSIZED] — assign Size on the project board, then re-run")
            return 2
        print("Proposed due_on: (no open work)")
        print("Status: [No change] — milestone has no remaining capacity")
        return 0

    calendar_days = int(round(remaining / AVAILABILITY_RATE * 7))
    proposed_due = date.today() + timedelta(days=calendar_days)
    delta = (proposed_due - current_due_date).days if current_due_date else None
    delta_str = f"{delta:+d}" if delta is not None else "n/a"
    print(f"Proposed due_on: {proposed_due} (delta {delta_str} days)")

    if delta is None or abs(delta) <= THRESHOLD_DAYS:
        print("Status: [No change]")
        return 0
    print("Status: [UPDATE NEEDED]")
    return 2


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--issue", type=int, help="Issue number; recheck its milestone")
    group.add_argument("--milestone", type=int, help="Milestone number")
    args = parser.parse_args()

    if args.issue:
        ms = milestone_for_issue(args.issue)
        if ms is None:
            print(f"error: issue #{args.issue} has no milestone", file=sys.stderr)
            return 1
        return compute_recheck(ms)
    return compute_recheck(args.milestone)


if __name__ == "__main__":
    sys.exit(main())
