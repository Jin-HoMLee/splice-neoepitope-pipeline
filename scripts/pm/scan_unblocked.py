#!/usr/bin/env python3
"""Freshly-unblocked sweep - surface a dependent when its blocker clears (Issue #745).

## The failure mode

An Issue parked because it was blocked can strand silently when its blocker
**closes**: nothing re-surfaces it, so it sits uncommitted while it is actually
available. The establishing incident is Issue #594 - its prerequisites (#211,
#212) closed on 2026-06-11 and it sat un-pulled for 4 days, caught only by a
manual dig during the morning routine.

## Why the obvious rule is WRONG

The naive detector - *"unblocked AND still in Backlog"* - would fight our own
model. Under late-commitment Kanban (`shared/feedback_board_hygiene.md`), sitting
uncommitted in Backlog is the **normal, correct** resting state of an option, not
drift. Flagging it would emit a permanent nag against items that are exactly where
they belong, and a nag that can only say one thing is a hollow check.

The real signal is the **state change**, not the standing state:

    a blocker closed RECENTLY, and the dependent it was blocking is STILL uncommitted.

That is the #594 shape precisely (the clear went unnoticed for 4 days), and it is
a **replenishment input** - "this is newly available, reconsider it" - never an
instruction to commit. Advisory by house style: it surfaces, it never blocks and
never auto-commits.

## Scope note

This Issue originally also carried a *near-deadline milestone sweep*. That half was
**cut** on 2026-07-14: Issue #902 facet 2 made us an all-Kanban shop, so flow work
commits milestone-free, and only 2 open Issues board-wide still carry a milestone -
the one open milestone with a `due_on` has zero open issues. A deadline sweep would
fire on the empty set, i.e. a check that cannot fail. See the Issue for the numbers.

Usage::

    python3 scripts/pm/scan_unblocked.py                 # sweep, human-readable
    python3 scripts/pm/scan_unblocked.py --since-days 30 # widen the lookback
    python3 scripts/pm/scan_unblocked.py --json          # machine-readable
    python3 scripts/pm/scan_unblocked.py --check         # exit 1 on any finding
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gh_client import gh, GhError  # noqa: E402

REPO = "Jin-HoMLee/splice-neoepitope-pipeline"
OWNER, NAME = REPO.split("/")
PROJECT_NUMBER = 9

# How recently a blocker must have closed for the clear to count as a *state
# change* rather than a standing state. 14d comfortably covers a weekly
# Replenishment cadence (a clear cannot slip through between two sweeps) without
# re-flagging options that have simply been resting for a month.
DEFAULT_SINCE_DAYS = 14

# Board Statuses that mean "not yet committed". A cleared blocker only matters if
# nobody has since pulled the dependent. `Epic` is deliberately absent: a parent is
# parked there by design (Pattern A2) and is never itself pulled, so a parent whose
# blocker clears is not a missed commitment.
UNCOMMITTED_STATUSES = frozenset({None, "", "No Status", "Backlog"})

_QUERY = """
query($owner:String!, $name:String!, $after:String) {
  repository(owner:$owner, name:$name) {
    issues(first: 100, states: OPEN, after: $after) {
      pageInfo { hasNextPage endCursor }
      nodes {
        number
        title
        url
        labels(first: 20) { nodes { name } }
        blockedBy(first: 20) { nodes { number title state closedAt } }
        projectItems(first: 5) {
          nodes {
            project { number }
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
}
"""


# --------------------------------------------------------------------------
# Pure layer - no I/O. Everything that decides anything lives here.
# --------------------------------------------------------------------------

def parse_ts(value):
    """Parse a GitHub ISO-8601 timestamp into an aware datetime, or None."""
    if not value:
        return None
    return dt.datetime.fromisoformat(value.replace("Z", "+00:00"))


def board_status(item_nodes, project_number=PROJECT_NUMBER):
    """Extract the Status single-select value for our board, or None if unset.

    None covers three genuinely different cases that all mean 'uncommitted':
    the issue is not on the board, it is on the board with no Status, or the
    Status field is absent from the payload.
    """
    for node in item_nodes or []:
        if (node.get("project") or {}).get("number") != project_number:
            continue
        for fv in ((node.get("fieldValues") or {}).get("nodes") or []):
            if ((fv.get("field") or {}).get("name")) == "Status":
                return fv.get("name")
    return None


def classify(issue, *, now, since_days):
    """Decide whether one issue is a freshly-unblocked finding. Pure.

    Returns a finding dict, or None. Four ways to be silent, each load-bearing:

    1. **No wired blocker ever.** Nothing cleared, so there is no state change.
    2. **A blocker is still open.** Genuinely still blocked - not our business.
    3. **The last blocker cleared outside the window.** The clear is old news; the
       item is now simply a resting uncommitted option, and flagging it forever is
       the nag this sweep exists NOT to be.
    4. **Already committed** (Ready / In progress / review / Done / Epic). Somebody
       pulled it after the clear - the system worked, say nothing.
    """
    blockers = ((issue.get("blockedBy") or {}).get("nodes")) or []
    if not blockers:
        return None                                            # (1)
    if any(b.get("state") != "CLOSED" for b in blockers):
        return None                                            # (2)

    closed_ats = [ts for ts in (parse_ts(b.get("closedAt")) for b in blockers) if ts]
    if not closed_ats:
        # All blockers report CLOSED but none carries a closedAt. We cannot date the
        # state change, so we cannot honour the window. Stay silent rather than
        # invent a timestamp - a finding we cannot justify is worse than a miss.
        return None
    cleared_at = max(closed_ats)
    age_days = (now - cleared_at).total_seconds() / 86400.0
    if age_days > since_days:
        return None                                            # (3)

    status = board_status((issue.get("projectItems") or {}).get("nodes"))
    if status not in UNCOMMITTED_STATUSES:
        return None                                            # (4)

    roles = sorted(
        label["name"]
        for label in (((issue.get("labels") or {}).get("nodes")) or [])
        if label.get("name", "").startswith("role:")
    )
    return {
        "number": issue["number"],
        "title": issue.get("title", ""),
        "url": issue.get("url", ""),
        "status": status or "No Status",
        "roles": roles,
        "cleared_at": cleared_at.isoformat(),
        "age_days": round(age_days, 1),
        "cleared_by": [
            {"number": b["number"], "title": b.get("title", ""), "closedAt": b.get("closedAt")}
            for b in sorted(blockers, key=lambda b: b["number"])
        ],
    }


def render(findings, *, since_days):
    """Human-readable report. Pure."""
    if not findings:
        return (f"freshly-unblocked sweep: no findings "
                f"(no blocker cleared in the last {since_days}d with its dependent "
                f"still uncommitted).")
    lines = [
        f"freshly-unblocked sweep: {len(findings)} finding(s) - a blocker cleared "
        f"within {since_days}d and the dependent is still uncommitted.",
        "Advisory: this is a Replenishment input (reconsider these), never an auto-commit.",
        "",
    ]
    for f in findings:
        roles = ",".join(f["roles"]) or "role:UNLABELLED"
        lines.append(
            f"  [UNBLOCKED] #{f['number']} ({f['status']}, {roles}) - {f['title'][:60]}"
        )
        for b in f["cleared_by"]:
            lines.append(f"      cleared by #{b['number']} ({b['title'][:50]})")
        lines.append(f"      last blocker closed {f['age_days']}d ago -> {f['url']}")
        lines.append("")
    return "\n".join(lines).rstrip()


# --------------------------------------------------------------------------
# I/O layer
# --------------------------------------------------------------------------

def fetch_open_issues():
    """All open issues with their blockedBy edges + board Status, paginated.

    Paginates on `pageInfo.hasNextPage` - `first: 100` is a cap, not a filter, and a
    single-page read would silently drop everything past the cap
    (`shared/feedback_board_queries.md`).
    """
    nodes, cursor = [], None
    while True:
        args = [
            "api", "graphql",
            "-f", f"query={_QUERY}",
            "-F", f"owner={OWNER}",
            "-F", f"name={NAME}",
        ]
        if cursor:
            args += ["-F", f"after={cursor}"]
        data = gh(*args)
        conn = data["data"]["repository"]["issues"]
        nodes.extend(conn["nodes"])
        page = conn["pageInfo"]
        if not page["hasNextPage"]:
            return nodes
        cursor = page["endCursor"]


def sweep(*, since_days=DEFAULT_SINCE_DAYS, now=None):
    now = now or dt.datetime.now(dt.timezone.utc)
    issues = fetch_open_issues()
    findings = [
        f for f in (classify(i, now=now, since_days=since_days) for i in issues) if f
    ]
    return sorted(findings, key=lambda f: f["age_days"])


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--since-days", type=int, default=DEFAULT_SINCE_DAYS,
                        help=f"lookback window for the blocker's close (default {DEFAULT_SINCE_DAYS})")
    parser.add_argument("--json", action="store_true", help="emit findings as JSON")
    parser.add_argument("--check", action="store_true",
                        help="exit 1 if any finding (for a hygiene gate); default is advisory exit 0")
    args = parser.parse_args()

    try:
        findings = sweep(since_days=args.since_days)
    except GhError as exc:
        # Fail open + loud: a sweep that cannot reach GitHub must not masquerade as
        # a clean board (the silent-green failure this repo keeps re-learning).
        print(f"freshly-unblocked sweep: FAILED to query GitHub - {exc}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(findings, indent=2))
    else:
        print(render(findings, since_days=args.since_days))

    return 1 if (args.check and findings) else 0


if __name__ == "__main__":
    sys.exit(main())
