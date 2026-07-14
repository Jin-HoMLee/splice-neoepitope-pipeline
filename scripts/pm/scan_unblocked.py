#!/usr/bin/env python3
"""Freshly-unblocked sweep - surface a dependent when its blocker clears (Issue #745).

## The failure mode

An Issue parked because it was blocked can strand silently when its blocker
**closes**: nothing re-surfaces it, so it sits uncommitted while it is actually
available. The establishing incident is Issue #594 - its prerequisites (#211,
#212) closed on 2026-06-11 and it sat un-pulled for 4 days, caught only by a
manual dig during the morning routine.

## LEVEL-TRIGGERED, not edge-inferring - and why that distinction is the whole design

The naive detector - *"unblocked AND still in Backlog"* - is wrong: under
late-commitment Kanban (`shared/feedback_board_hygiene.md`), sitting uncommitted in
Backlog is the **normal, correct** resting state of an option, not drift. Flagging it
would nag forever about items that are exactly where they belong.

The **first** cut of this script fixed that with a time window: fire only if a blocker
closed within the last N days. That was an **edge-inferring** design - it reconstructed
"an event happened" from a state snapshot - and it carried a silent-data-loss bug:

    if the sweep did not RUN inside the window, the finding was dropped FOREVER.

Two weeks of not running the morning routine and an unblock event vanishes with no
trace. That is exactly the failure mode the webhook/event literature exists to warn
about, and it is worse than the bug it was written to avoid.

So the window is **gone**. This is now a proper **level-triggered reconciliation loop**
(the Kubernetes/Flux pattern): it does not care about events at all, it looks at the
world *right now* and asks whether it is in the desired state. The desired state is:

    every unblocked Issue has HAD its commitment decision.

and the condition for a finding is a pure statement about the present:

    >=1 wired blocker, ALL of them CLOSED, still uncommitted, and NOT yet acknowledged.

**No time window. Nothing can expire. Missing a run costs nothing** - the next run sees
exactly the same world and reports exactly the same thing.

## Why an acknowledgement label, and why this is not a nag

A level-triggered loop keeps reporting until the world reaches the desired state, so
each finding needs a **terminal action**. There are exactly two, and both are honest:

- **Commit it** -> board Status leaves the uncommitted set -> the finding clears itself.
- **Decline it** (deliberately leave it resting) -> apply the `unblock-ack` label ->
  the finding clears.

That is the difference between a to-do list and a nag: a nag is a warning you cannot
make stop. Here every finding has a definite way to stop, and the way to stop it *is*
the decision we wanted made. An item keeps appearing precisely as long as nobody has
decided about it - which is the point, not a defect.

The output is still **advisory** (house style): it surfaces, it never blocks, and it
never auto-commits. `--ack` records a decline; it does not make one.

## `--check` is a FORCE-DISPOSAL gate, not a passive canary - wire it knowingly

Because findings no longer expire, `--check` stays **exit 1 until every finding is
committed or acked**. Under the old window it would eventually self-clear on its own;
now it will not, by design - each finding *should* get a terminal decision. So `--check`
in a **blocking** context means *"force a decision"*, not *"warn me if something looks
off"*. That is intended, but a consumer expecting eventual auto-green would be surprised,
so do not wire it into a blocking gate without wanting that semantic. The default
(no `--check`) is advisory exit 0, which is the right default and what the morning
routine uses.

**Exit codes:** `0` = ran fine (advisory, or `--check` with nothing to decide); `1` =
ran fine, `--check` gate tripped; `2` = **could not run** (GitHub unreachable / ack
write failed), matching argparse's own exit-2-for-usage-error convention. Note this is
**inverted relative to `scan_prose_deps.py --check`** (which uses 2=drift, 1=infra
error). Both are internally consistent; a routine shelling several PM sweeps must not
assume a uniform mapping. Aligning them is tracked separately rather than silently
changed here.

## Scope note

This Issue originally also carried a *near-deadline milestone sweep*. That half was
**cut** on 2026-07-14: Issue #902 facet 2 made us an all-Kanban shop, so flow work
commits milestone-free, and only 2 open Issues board-wide still carry a milestone -
the one open milestone with a `due_on` has zero open issues. A deadline sweep would
fire on the empty set, i.e. a check that cannot fail. See the Issue for the numbers.

Usage::

    python3 scripts/pm/scan_unblocked.py                    # sweep, human-readable
    python3 scripts/pm/scan_unblocked.py --json             # machine-readable
    python3 scripts/pm/scan_unblocked.py --check            # exit 1 on any finding
    python3 scripts/pm/scan_unblocked.py --ack 594 --reason "resting on purpose: ..."
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

# The acknowledgement marker: "we looked at this cleared blocker and deliberately
# chose to leave the Issue resting." It is what lets this loop be level-triggered
# WITHOUT nagging - the decline is recorded in the world, so the loop can read it
# back and go quiet. A time window would have done the same job by *forgetting*,
# which silently loses any finding the sweep did not happen to run in time to see.
# Remove the label and the finding legitimately returns (the decision is revoked).
ACK_LABEL = "unblock-ack"
ACK_LABEL_COLOR = "0E8A16"
ACK_LABEL_DESCRIPTION = (
    "Unblocked and deliberately left resting; excluded from scan_unblocked.py. "
    "Remove to bring it back."
)

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
        # first: 50 is GitHub's documented per-direction ceiling for blockers, so this
        # is a TRUE BOUND, not a sample (same reasoning + citation as
        # scan_prose_deps.native_blockers). The distinction is load-bearing here: a cap
        # is not a filter, and GraphQL does not guarantee node order, so a *sampling*
        # cap could return N blockers that all read CLOSED while a still-OPEN one sorts
        # past the cap - and classify() would then fire a false "freshly-unblocked"
        # finding on an issue that is genuinely still blocked. That is the worst
        # failure this sweep can have: it would feed "never commit a blocked Issue to
        # Ready" the exact input that rule exists to prevent.
        blockedBy(first: 50) { nodes { number title state closedAt } }
        # Same cap-vs-filter reasoning, lower stakes: if an issue sat on more projects
        # than this and board #9's item sorted last, board_status() would return None,
        # which reads as "uncommitted" and could false-fire on an already-committed
        # item. We run one project, so any bound clears it; 20 is free.
        projectItems(first: 20) {
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


def classify(issue, *, now=None):
    """Decide whether one issue is an unblocked-but-undecided finding. Pure.

    A statement about the world RIGHT NOW - no event history, no clock dependency.
    `now` is used only to *describe* how long a finding has been waiting; it never
    decides anything, so a sweep run today and the same sweep run in a year return
    the same findings. That is what makes this level-triggered.

    Returns a finding dict, or None. Four ways to be silent, each load-bearing:

    1. **No wired blocker ever.** Nothing to clear; most of the board.
    2. **A blocker is still open.** Genuinely still blocked - not our business.
    3. **Already committed** (Ready / In progress / review / Done / Epic). Somebody
       pulled it after the clear - the system worked, say nothing.
    4. **Acknowledged** (`unblock-ack`). Somebody looked and deliberately chose to
       leave it resting. That is a decision, and the decision is what we wanted; the
       loop records it and goes quiet. Remove the label and it legitimately returns.

    Note what is NOT here: a time window. The first cut expired a finding after N days,
    which meant a sweep that did not RUN in time dropped it forever. Silent loss is a
    worse failure than a repeat mention, and a level-triggered loop does not need to
    remember anything - it just re-reads the world.
    """
    blockers = ((issue.get("blockedBy") or {}).get("nodes")) or []
    if not blockers:
        return None                                            # (1)
    if any(b.get("state") != "CLOSED" for b in blockers):
        return None                                            # (2)

    status = board_status((issue.get("projectItems") or {}).get("nodes"))
    if status not in UNCOMMITTED_STATUSES:
        return None                                            # (3)

    labels = [lb.get("name", "") for lb in (((issue.get("labels") or {}).get("nodes")) or [])]
    if ACK_LABEL in labels:
        return None                                            # (4)
    roles = sorted(n for n in labels if n.startswith("role:"))

    # Descriptive only - never a gate. A blocker legitimately may carry no closedAt
    # (deleted/transferred edge cases); that must not suppress a finding, because the
    # finding does not depend on WHEN it cleared, only THAT it has.
    now = now or dt.datetime.now(dt.timezone.utc)
    closed_ats = [ts for ts in (parse_ts(b.get("closedAt")) for b in blockers) if ts]
    cleared_at = max(closed_ats) if closed_ats else None
    waiting_days = (
        round((now - cleared_at).total_seconds() / 86400.0, 1) if cleared_at else None
    )

    return {
        "number": issue["number"],
        "title": issue.get("title", ""),
        "url": issue.get("url", ""),
        "status": status or "No Status",
        "roles": roles,
        "cleared_at": cleared_at.isoformat() if cleared_at else None,
        "waiting_days": waiting_days,
        "cleared_by": [
            {"number": b["number"], "title": b.get("title", ""), "closedAt": b.get("closedAt")}
            for b in sorted(blockers, key=lambda b: b["number"])
        ],
    }


def render(findings):
    """Human-readable report. Pure."""
    if not findings:
        return ("unblocked-but-undecided sweep: no findings "
                "(every unblocked Issue has been committed or acknowledged).")
    lines = [
        f"unblocked-but-undecided sweep: {len(findings)} finding(s) - every blocker is "
        f"CLOSED and the Issue is still uncommitted and unacknowledged.",
        "Advisory: a Replenishment input (reconsider these), never an auto-commit.",
        f"To clear one: COMMIT it (Status leaves Backlog), or DECLINE it "
        f"(--ack <N> --reason ... applies `{ACK_LABEL}`). Both are correct answers.",
        "",
    ]
    for f in findings:
        roles = ",".join(f["roles"]) or "role:UNLABELLED"
        lines.append(
            f"  [UNBLOCKED] #{f['number']} ({f['status']}, {roles}) - {f['title'][:60]}"
        )
        for b in f["cleared_by"]:
            lines.append(f"      cleared by #{b['number']} ({b['title'][:50]})")
        waited = (f"undecided for {f['waiting_days']}d"
                  if f["waiting_days"] is not None else "clear date unknown")
        lines.append(f"      {waited} -> {f['url']}")
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


def sweep(*, now=None):
    now = now or dt.datetime.now(dt.timezone.utc)
    issues = fetch_open_issues()
    findings = [f for f in (classify(i, now=now) for i in issues) if f]
    # Longest-undecided first: the one that has been waiting on a decision the longest
    # is the one most likely to have been forgotten. Undateable clears sort last.
    return sorted(findings, key=lambda f: (f["waiting_days"] is None,
                                           -(f["waiting_days"] or 0)))


def ensure_ack_label():
    """Create the ack label if missing. Idempotent (`--force` upserts).

    The write path MUST be self-sufficient. `gh issue edit --add-label` does not create
    a missing label - it fails HTTP 422, which `gh_client` treats as deterministic and
    does not retry. So without this, the FIRST `--ack` anyone runs on a fresh clone dies,
    and `--ack` is one of the two terminal actions the whole level-triggered design rests
    on. A design whose principle is "every finding must have a working way to be cleared"
    cannot ship with its clear-path gated on an undocumented manual step.

    (The label does exist in our repo today - because I created it by hand while testing,
    which is exactly the kind of invisible prerequisite that makes a thing unreproducible.
    Creating it here puts it in the code instead of in my head.)
    """
    gh("label", "create", ACK_LABEL, "--repo", REPO, "--force",
       "--color", ACK_LABEL_COLOR, "--description", ACK_LABEL_DESCRIPTION,
       parse_json=False)


def acknowledge(number, reason):
    """Record a DELIBERATE decline: label the Issue + leave the reason on the record.

    This does not make the decision - a human does. It writes the decision down so the
    level-triggered loop can read it back and stop asking. The reason is mandatory
    precisely because an unexplained silence is what we are trying to eliminate: a
    finding that vanishes with no rationale is indistinguishable from one that was lost.
    """
    ensure_ack_label()
    body = (
        f"**From:** PM\n\n"
        f"## Acknowledged: unblocked, and deliberately left resting\n\n"
        f"Every wired blocker on this Issue is closed, so it is **available**. "
        f"We looked and chose **not** to commit it now.\n\n"
        f"**Reason:** {reason}\n\n"
        f"Labelled `{ACK_LABEL}`, so the unblocked-but-undecided sweep "
        f"(`scripts/pm/scan_unblocked.py`, Issue #745) stops surfacing it. This is a "
        f"decision, not a dismissal - **remove the label and it returns to the sweep.**\n\n"
        f"**Created by:** PM"
    )
    gh("issue", "comment", str(number), "--repo", REPO, "--body", body, parse_json=False)
    gh("issue", "edit", str(number), "--repo", REPO, "--add-label", ACK_LABEL,
       parse_json=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--json", action="store_true", help="emit findings as JSON")
    parser.add_argument("--check", action="store_true",
                        help="exit 1 if any finding (for a hygiene gate); default is advisory exit 0")
    parser.add_argument("--ack", type=int, metavar="N",
                        help=f"record a deliberate decline on Issue N (applies `{ACK_LABEL}`)")
    parser.add_argument("--reason", help="why it is being left resting (required with --ack)")
    args = parser.parse_args()

    if args.ack is not None:   # `--ack 0` is falsy; `is not None` is the precise guard
        if not args.reason:
            parser.error("--ack requires --reason: an unexplained decline is "
                         "indistinguishable from a dropped finding")
        try:
            acknowledge(args.ack, args.reason)
        except GhError as exc:
            print(f"unblocked sweep: FAILED to acknowledge #{args.ack} - {exc}", file=sys.stderr)
            return 2
        print(f"acknowledged #{args.ack} ({ACK_LABEL}) - it will no longer surface. "
              f"Remove the label to bring it back.")
        return 0

    try:
        findings = sweep()
    except GhError as exc:
        # Fail open + loud: a sweep that cannot reach GitHub must not masquerade as
        # a clean board (the silent-green failure this repo keeps re-learning).
        print(f"unblocked sweep: FAILED to query GitHub - {exc}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(findings, indent=2))
    else:
        print(render(findings))

    return 1 if (args.check and findings) else 0


if __name__ == "__main__":
    sys.exit(main())
