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


def find_prior_same_stage(
    iteration: int,
    stage: int,
    all_milestones: list[dict],
) -> dict | None:
    """Highest-iteration prior in the same S-stage chain. Includes closed."""
    candidates = []
    for ms in all_milestones:
        parsed = parse_milestone_title(ms["title"])
        if parsed is None:
            continue
        n, s = parsed
        if s == stage and n < iteration:
            candidates.append((n, ms))
    if not candidates:
        return None
    candidates.sort(key=lambda x: x[0])
    return candidates[-1][1]


def find_open_same_iteration_S5(
    iteration: int,
    all_milestones: list[dict],
) -> dict | None:
    """Find an OPEN i<N> - S5 - ... milestone. Used for paired-S7 gating.

    Loose match: same iteration number is enough. Arc-mismatch (e.g. i4-S7
    'TCR-pMHC Landscape' vs i4-S5 'Google Batch') is a separate data-hygiene
    concern not solvable here.
    """
    for ms in all_milestones:
        parsed = parse_milestone_title(ms["title"])
        if parsed is None:
            continue
        n, s = parsed
        if n == iteration and s == 5 and ms["state"] == "open":
            return ms
    return None


def compute_layered_due_date(
    iteration: int | None,
    stage: int | None,
    capacity_days: float,
    all_milestones: list[dict],
) -> tuple[date, str]:
    """Return (proposed_due, reasoning_note).

    Handles 4 main branches:
      1. No title parse (role-meta etc.) -> pure capacity
      2. S7 paired with open S5 in same iteration -> stack on S5 close
      3. S7 standalone (no paired S5) -> pure capacity
      4. Same-S-stage stacking (closed/undated/normal/overdue prior cases)
    """
    today = date.today()

    if iteration is None or stage is None:
        # Non-S-stage milestone (pm-i*, dev-i*, M1, etc.) — pure capacity
        base = today
        note = ""
    elif stage == 7:
        paired = find_open_same_iteration_S5(iteration, all_milestones)
        if paired and paired.get("due_on"):
            paired_date = date.fromisoformat(paired["due_on"][:10])
            base = max(paired_date, today)
            note = f"(paired-S7: unblocks at M#{paired['number']} close {paired['due_on'][:10]})"
        else:
            # Standalone S7 (e.g. Lit Review i3-S7) — pure capacity
            base = today
            note = "(standalone S7 — no paired open S5)"
    else:
        prior = find_prior_same_stage(iteration, stage, all_milestones)
        if prior is None:
            base = today
            note = "(no prior same-S milestone)"
        elif prior["state"] == "closed":
            base = today
            note = f"(prior M#{prior['number']} closed)"
        elif prior.get("due_on") is None:
            base = today
            note = f"(prior M#{prior['number']} undated — sequencing skipped)"
        else:
            prior_date = date.fromisoformat(prior["due_on"][:10])
            base = max(prior_date, today)
            note = f"(stack after M#{prior['number']} close {prior['due_on'][:10]})"

    calendar_days = int(round(capacity_days / AVAILABILITY_RATE * 7))
    proposed = base + timedelta(days=calendar_days)
    return (proposed, note)


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


def parent_numbers(issue_numbers: list[int]) -> set[int]:
    """Subset of issue_numbers that are parent epics (subIssuesSummary.total > 0).

    Parents carry no Size by convention — size rolls up from their sub-issues
    (shared/feedback_parent_sub_issues.md). compute_recheck excludes them from
    both the capacity sum and the unsized-check, so a milestone holding a parent
    as a roadmap anchor doesn't emit a spurious, un-clearable [UNSIZED] flag
    (Issue #689 — sizing the parent to clear it would itself violate the
    no-size convention).
    """
    if not issue_numbers:
        return set()
    owner, name = REPO.split("/")
    aliases = " ".join(
        f"i{n}: issue(number: {n}) {{ subIssuesSummary {{ total }} }}"
        for n in issue_numbers
    )
    query = f'query {{ repository(owner: "{owner}", name: "{name}") {{ {aliases} }} }}'
    data = gh("api", "graphql", "-f", f"query={query}")
    repo = data["data"]["repository"]
    parents: set[int] = set()
    for n in issue_numbers:
        node = repo.get(f"i{n}") or {}
        total = (node.get("subIssuesSummary") or {}).get("total") or 0
        if total > 0:
            parents.add(n)
    return parents


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
    parents = parent_numbers(issue_numbers)

    print(f"Milestone: {title}")
    print(f"Current due_on: {current_due_date or '—'}")
    print(f"Open issues ({len(issue_numbers)}):")
    remaining = 0.0
    leaf_numbers: list[int] = []
    for n in sorted(issue_numbers):
        # Parent epics are roadmap anchors only — they carry no Size by convention
        # (size rolls up from sub-issues, shared/feedback_parent_sub_issues.md), so
        # exclude them from both the capacity sum and the unsized-check below. Without
        # this a parent-anchored milestone emits an un-clearable [UNSIZED] flag
        # (Issue #689; surfaced 2026-06-11 on pm-i6 with parents #527/#538).
        if n in parents:
            print(f"  - #{n} (parent epic — excluded; size rolls up from sub-issues)")
            continue
        leaf_numbers.append(n)
        size = sizes.get(n)
        weight = SIZE_WEIGHTS.get(size or "", 0)
        remaining += weight
        size_disp = f"{size}, ~{weight}d" if size else "no size, 0d"
        print(f"  - #{n} ({size_disp})")
    print(f"Remaining capacity: {remaining}d")

    # Any unsized open LEAF issue makes the capacity read unreliable: the
    # size-weighted sum silently under-counts true remaining work, which then drives
    # a confident-but-bogus due_on proposal (the 2026-06-03 pm-i6 under-read — one
    # sized issue + two unsized yielded 1.0d and a spurious -28d slip; Issue #618
    # AC2). Flag and bail BEFORE computing a due date, regardless of whether the
    # sized subset happens to be > 0 — previously this guard was nested inside the
    # `remaining == 0` branch and so was bypassed whenever even one issue had a size.
    # Parents are excluded above, so an unsized parent never trips this (Issue #689).
    unsized_count = sum(1 for n in leaf_numbers if sizes.get(n) is None)
    if unsized_count > 0:
        floor = f"; sized subset {remaining}d is a floor" if remaining else ""
        print(f"Proposed due_on: (cannot compute — {unsized_count} open issue(s) missing Size{floor})")
        print("Status: [UNSIZED] — assign Size on the project board, then re-run")
        return 2

    if remaining == 0:
        print("Proposed due_on: (no open work)")
        print("Status: [No change] — milestone has no remaining capacity")
        return 0

    # Sequencing-aware: fetch all milestones once, derive (iteration, stage), apply layered logic
    all_milestones_raw = gh("api", "--paginate", f"repos/{REPO}/milestones?state=all&per_page=100")
    parsed = parse_milestone_title(title)
    iteration, stage = parsed if parsed else (None, None)
    proposed_due, note = compute_layered_due_date(iteration, stage, remaining, all_milestones_raw)

    delta = (proposed_due - current_due_date).days if current_due_date else None
    delta_str = f"{delta:+d}" if delta is not None else "n/a"
    note_suffix = f" {note}" if note else ""
    print(f"Proposed due_on: {proposed_due} (delta {delta_str} days){note_suffix}")

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
