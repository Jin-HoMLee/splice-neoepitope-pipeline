#!/usr/bin/env python3
"""Recheck parent-vs-children Status drift on project board #9.

Implements the audit rule from docs/superpowers/specs/2026-05-19-parent-status-drift-audit-design.md.

Usage:
  scripts/pm/recheck_parent_status.py --issue N      # recheck this issue's parent chain
  scripts/pm/recheck_parent_status.py --all          # audit all parent issues on project #9

Exits 0 if no drift, 2 if drift detected, 1 on error.
"""
from __future__ import annotations

STATUS_LADDER = {
    "Backlog": 0,
    "Ready": 1,
    "In progress": 2,
    "Ready for review": 3,
    "In review": 4,
    "Done": 5,
}


def rank(status: str | None) -> int:
    """Map a Status name to its precedence rank. None/unknown → 0 (Backlog)."""
    return STATUS_LADDER.get(status or "", 0)


LADDER_INVERSE = {v: k for k, v in STATUS_LADDER.items()}


def collective_state(open_children: list[dict]) -> str:
    """Max-rank status across open children. Empty list → 'Done'."""
    if not open_children:
        return "Done"
    max_rank = max(rank(c.get("status")) for c in open_children)
    return LADDER_INVERSE[max_rank]


def classify_drift(parent_status: str | None, open_children: list[dict]) -> str | None:
    """Classify drift for a parent vs its open children.

    Returns one of: 'FORWARD DRIFT', 'BACKWARD DRIFT', 'COMPLETION DRIFT', or None.
    """
    p_rank = rank(parent_status)
    if not open_children:
        # All children closed; parent should be Done
        return None if p_rank == STATUS_LADDER["Done"] else "COMPLETION DRIFT"
    c_rank = rank(collective_state(open_children))
    if p_rank > c_rank and p_rank >= STATUS_LADDER["In progress"]:
        # Only flag when the parent is making a falsifiable progress claim (In progress
        # or beyond). A Ready parent with Backlog children is the normal post-grooming
        # state and should not be treated as drift.
        return "FORWARD DRIFT"
    if p_rank < c_rank:
        return "BACKWARD DRIFT"
    return None


import json
import subprocess

REPO = "Jin-HoMLee/splice-neoepitope-pipeline"
PROJECT_NUMBER = 9


def gh(*args: str, parse_json: bool = True):
    """Invoke `gh` and parse JSON output (or return raw text)."""
    result = subprocess.run(["gh", *args], capture_output=True, text=True, check=True)
    return json.loads(result.stdout) if parse_json else result.stdout


def parent_issue_number(issue_number: int) -> int | None:
    """Return parent issue number via REST parent_issue_url, or None."""
    data = gh("api", f"repos/{REPO}/issues/{issue_number}")
    url = data.get("parent_issue_url")
    if not url:
        return None
    return int(url.rstrip("/").rsplit("/", 1)[-1])


def open_sub_issues(issue_number: int) -> list[dict]:
    """Return list of open sub-issues (each a dict with 'number' at minimum)."""
    data = gh("api", f"repos/{REPO}/issues/{issue_number}/sub_issues")
    return [c for c in data if c.get("state") == "open"]


def status_for_issue(issue_number: int) -> str | None:
    """Return the project #9 Status name for issue, or None if not on the project."""
    owner, name = REPO.split("/")
    query = (
        f'query {{ repository(owner: "{owner}", name: "{name}") {{ '
        f'issue(number: {issue_number}) {{ '
        f'projectItems(first: 5) {{ nodes {{ '
        f'project {{ number }} '
        f'fieldValues(first: 20) {{ nodes {{ '
        f'... on ProjectV2ItemFieldSingleSelectValue {{ '
        f'name field {{ ... on ProjectV2SingleSelectField {{ name }} }} '
        f'}} }} }} '
        f'}} }} }} }} }}'
    )
    data = gh("api", "graphql", "-f", f"query={query}")
    nodes = (
        data.get("data", {})
        .get("repository", {})
        .get("issue", {})
        .get("projectItems", {})
        .get("nodes", [])
    ) or []
    for pi in nodes:
        if (pi.get("project") or {}).get("number") != PROJECT_NUMBER:
            continue
        for fv in pi.get("fieldValues", {}).get("nodes", []):
            if (fv.get("field") or {}).get("name") == "Status":
                return fv.get("name")
    return None


def audit_parent_chain(issue_number: int) -> list[dict]:
    """Walk up the parent chain from issue_number; audit drift at each level.

    Returns list of records: {issue, status, open_children, collective, drift}.
    Empty list if issue has no parent.
    """
    chain: list[dict] = []
    cursor = parent_issue_number(issue_number)
    seen = {issue_number}
    while cursor is not None and cursor not in seen:
        seen.add(cursor)
        parent_status = status_for_issue(cursor)
        children = open_sub_issues(cursor)
        # Enrich children with their Status
        enriched = [{"number": c["number"], "status": status_for_issue(c["number"])}
                    for c in children]
        chain.append({
            "issue": cursor,
            "status": parent_status,
            "open_children": enriched,
            "collective": collective_state(enriched),
            "drift": classify_drift(parent_status, enriched),
        })
        cursor = parent_issue_number(cursor)
    return chain


def format_record(record: dict) -> str:
    """Render one audit record as a multi-line block matching the spec format."""
    lines: list[str] = []
    lines.append(f"#{record['issue']} — Status: {record['status'] or 'NO STATUS'}")
    children = record["open_children"]
    lines.append(f"  Open sub-issues ({len(children)}):")
    for c in children:
        lines.append(f"    - #{c['number']} ({c['status'] or 'NO STATUS'})")
    lines.append(f"  Collective children state: {record['collective']}")
    drift = record["drift"]
    if drift is None:
        lines.append(f"  Status: [No change]")
    else:
        lines.append(f"  Status: [{drift}]")
    return "\n".join(lines)


import argparse
import sys


def run_issue_mode(issue_number: int) -> int:
    chain = audit_parent_chain(issue_number)
    if not chain:
        print(f"Issue #{issue_number} has no parent — nothing to audit.")
        return 0
    print(f"Parent chain for #{issue_number} (walked {len(chain)} levels):\n")
    drifted = False
    for record in chain:
        print(format_record(record))
        print()
        if record["drift"] is not None:
            drifted = True
    return 2 if drifted else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--issue", type=int, help="Issue number; walk its parent chain")
    group.add_argument("--all", action="store_true", help="Audit all parent issues on project #9")
    args = parser.parse_args(argv)

    if args.issue:
        return run_issue_mode(args.issue)
    return run_all_mode()


def all_parent_issues() -> list[int]:
    """Return numbers of all open issues in the repo that have ≥1 sub-issue."""
    # GitHub REST search: filter open issues, then check sub_issues_summary.total > 0
    # via a per-issue follow-up. Use gh issue list + per-issue REST fetch.
    data = gh("issue", "list", "--repo", REPO, "--state", "open", "--limit", "200",
              "--json", "number")
    parents: list[int] = []
    for issue in data:
        n = issue["number"]
        meta = gh("api", f"repos/{REPO}/issues/{n}")
        if (meta.get("sub_issues_summary") or {}).get("total", 0) > 0:
            parents.append(n)
    return parents


def run_all_mode() -> int:
    parents = all_parent_issues()
    drifted_count = 0
    drift_blocks: list[str] = []
    for p in parents:
        parent_status = status_for_issue(p)
        children = open_sub_issues(p)
        enriched = [{"number": c["number"], "status": status_for_issue(c["number"])}
                    for c in children]
        record = {
            "issue": p,
            "status": parent_status,
            "open_children": enriched,
            "collective": collective_state(enriched),
            "drift": classify_drift(parent_status, enriched),
        }
        if record["drift"] is not None:
            drifted_count += 1
            drift_blocks.append(format_record(record))

    print(f"Audited {len(parents)} parent issues; {drifted_count} drifted.\n")
    for block in drift_blocks:
        print(block)
        print()
    return 2 if drifted_count > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
