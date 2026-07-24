#!/usr/bin/env python3
"""List open items on the JH M Lee Lab project board (#9) with Status/Priority/Size.

Wraps the paginated ProjectV2 GraphQL query that the board needs:
  - GraphQL `first: N` is a cap, not a filter; we loop on `pageInfo.hasNextPage`.
  - The board sorts Done items first, so unpaginated single-page queries silently
    truncate open work. See `.agents/memory/shared/feedback_board_queries.md`.

For queries that only need label + state (no Status/Priority/Size columns),
prefer `gh issue list --label role:<role> --state open` instead — it's a one-liner
with no pagination needed.

Each item carries its issue/PR created/updated/closed timestamps; `--sort-updated`
(momentum) and `--stale-days N` (dormancy) filter/sort on `Issue.updatedAt` — the
content-level field, NOT the board "Updated" column (which a bare field nudge bumps;
see Issue #642).

**Two-repo aggregation (Issue #999).** Board #9 aggregates items from BOTH the
project repo and the personas repo (`claude-personas-splice-neoepitope-pipeline`),
whose issue/PR numbers COLLIDE (both have a #29, #64, #71, ...). Each item's true
identity is its `url`, never the bare number, so the text table disambiguates by
tagging the personas rows (`pers#71` vs a bare project `71` - the two collide),
and `--json` carries both `url` and an `origin` field (`project` / `personas` /
`other`). When acting on
a listed item, resolve the repo from its `url` before any `gh ... -R` call: a bare
`gh issue view 71` hits whichever repo your cwd defaults to and can silently return
the wrong same-numbered issue.

Usage:
  scripts/board_open_items.py
  scripts/board_open_items.py --role developer
  scripts/board_open_items.py --status Ready --priority P1
  scripts/board_open_items.py --sort-updated --role scientist   # recent momentum
  scripts/board_open_items.py --stale-days 21                    # dormancy sweep
  scripts/board_open_items.py --exclude-parents                  # drop epics (Size/flow sweeps)
  scripts/board_open_items.py --check-coherence                  # committed-but-arc-phase:later drift (#765)
  scripts/board_open_items.py --json | jq '.[] | select(.size == "S")'
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from typing import Any
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / ".agents" / "hooks"))
import graphql_meter  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent / "pm"))
from not_pullable import scan_not_pullable  # noqa: E402
from pullability import assess as assess_pullable  # noqa: E402

OWNER = "Jin-HoMLee"
PROJECT_NUMBER = 9

# Board #9 aggregates two repos with colliding numbers (Issue #999). Personas is
# checked first because its name CONTAINS the project name as a substring.
PERSONAS_REPO = "claude-personas-splice-neoepitope-pipeline"
PROJECT_REPO = "splice-neoepitope-pipeline"
def origin_from_url(url: str) -> str:
    """Classify a board item's origin repo from its URL: project / personas / other."""
    url = url or ""
    if f"/{PERSONAS_REPO}/" in url:
        return "personas"
    if f"/{PROJECT_REPO}/" in url:
        return "project"
    return "other"


def ref_cell(item: dict) -> str:
    """Disambiguated reference for the text table.

    Board #9 aggregates exactly two repos, so tagging ONLY the personas rows fully
    disambiguates them from the collision-partner project rows: a personas item
    renders `pers#123`, everything else (project - the expected default - or an
    unexpected/unknown origin) renders the bare `123`. Asymmetric on purpose: the
    overwhelmingly-common project case stays uncluttered.
    """
    number = item.get("number")
    if item.get("origin") == "personas":
        return f"pers#{number}"
    return str(number)

QUERY = """
query($owner: String!, $number: Int!, $after: String) {
  rateLimit { cost remaining }
  user(login: $owner) {
    projectV2(number: $number) {
      items(first: 100, after: $after) {
        pageInfo { hasNextPage endCursor }
        nodes {
          content {
            __typename
            ... on Issue {
              number title state url body
              createdAt updatedAt closedAt
              subIssuesSummary { total }
              milestone { title }
              labels(first: 20) { nodes { name } }
              blockedBy(first: 50) { nodes { number state } }
            }
            ... on PullRequest {
              number title state url isDraft
              createdAt updatedAt closedAt
              labels(first: 20) { nodes { name } }
            }
          }
          fieldValues(first: 20) {
            nodes {
              ... on ProjectV2ItemFieldSingleSelectValue {
                name
                field { ... on ProjectV2SingleSelectField { name } }
              }
              ... on ProjectV2ItemFieldDateValue {
                date
                field { ... on ProjectV2FieldCommon { name } }
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

# Statuses at/past the late-commitment Backlog→Ready boundary (CLAUDE.md "Board
# status governance"). Backlog / No Status are uncommitted; Done is filtered out
# upstream in normalize().
COMMITTED_STATUSES = frozenset({"Ready", "Ready for review", "In review", "In progress"})


def is_arc_phase_incoherent(it: dict[str, Any]) -> bool:
    """True iff a *committed* item is parked at arc-phase:later — a contradiction.

    Committed = board Status past the Backlog→Ready boundary OR a milestone
    assigned (the commitment signal, set at Backlog→Ready). arc-phase:later means
    "parked, not now", so committed+later is incoherent. Fails open (False) when
    arc-phase is absent — no arc opinion, no flag. Issue #765.
    """
    if it.get("arc_phase") != "later":
        return False
    return it.get("status") in COMMITTED_STATUSES or it.get("milestone") is not None


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
        graphql_meter.log_graphql_spend("board_open_items", data, query_name="board_page")
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
    # Draft items (e.g. the pinned "📌 Active arc slate" reference card, #759)
    # are board furniture, not tracked work — they carry no number/role/status
    # and would otherwise surface as spurious No-Status, role-less rows that
    # pollute intake/hygiene sweeps. Skip them.
    if content.get("__typename") == "DraftIssue":
        return None
    state = content.get("state")
    if state in ("CLOSED", "MERGED"):
        return None
    status = size = priority = target_date = start_date = None
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
        elif fname == "Target date":
            target_date = fv.get("date")
        elif fname == "Start date":
            start_date = fv.get("date")
    if status == "Done":
        return None
    labels = [l["name"] for l in content.get("labels", {}).get("nodes", [])]
    blocked_by = [b for b in (content.get("blockedBy") or {}).get("nodes", []) if b]
    # An Issue may carry several `role:` labels - the dual-role convention gives one
    # role the Lead/DRI and scopes individual ACs to the other. Keep the FULL set, for
    # the same reason `arcs` does (#1103): GitHub returns labels in an unstable order,
    # so filtering on the first one alone hides a dual-role item under a role that
    # varies between calls. #1153: that hid 12 open Issues from `--role developer`,
    # including one where Developer was the Lead.
    role_labels = [l for l in labels if l.startswith("role:")]
    # `role` keeps the first for single-value consumers (the Role column, --json
    # readers); `roles` carries the full set so the filter cannot silently truncate.
    role = role_labels[0] if role_labels else None
    arc_labels = [l for l in labels if l.startswith("arc:")]
    # A parent/epic has >=1 sub-issue. PRs (and the rare node the query returned
    # no summary for) have no subIssuesSummary block -> treated as a leaf. Mirrors
    # `.agents/hooks/check_gh_issue_develop_parent.py` (total > 0 -> parent).
    # Computed BEFORE the multi-arc warning: multi-arc is drift on a leaf but
    # legal on a parent (#1103).
    sub_summary = content.get("subIssuesSummary") or {}
    is_parent = (sub_summary.get("total") or 0) > 0
    # One arc per LEAF; a parent/initiative may legitimately span themes (#1103).
    # #1036 spans scoring-tcr-pmhc + immunogenicity-benchmark, and warning on it
    # pushed toward stripping a true label to satisfy the checker.
    if len(arc_labels) > 1 and not is_parent:
        print(
            f"Warning: issue #{content.get('number')} has multiple arc labels: {arc_labels}",
            file=sys.stderr,
        )
    # `arc` keeps the first for single-value consumers (filters, the Arc column);
    # `arcs` carries the full set so a multi-arc parent is not silently truncated.
    arc = arc_labels[0] if arc_labels else None
    arc_phase = next(
        (l.removeprefix("arc-phase:") for l in labels if l.startswith("arc-phase:")),
        None,
    )
    is_issue = content.get("__typename") != "PullRequest"
    return {
        "arcs": arc_labels,
        "number": content.get("number"),
        "title": content.get("title", ""),
        "url": content.get("url", ""),
        "origin": origin_from_url(content.get("url", "")),
        "kind": "PR" if content.get("__typename") == "PullRequest" else "Issue",
        "is_draft": content.get("isDraft", False),
        "is_parent": is_parent,
        "state": state,
        "status": status or "No Status",
        # Issues only — the PullRequest GraphQL fragment doesn't fetch milestone,
        # so a PR is "committed" via Status alone (it's always in a committed
        # status anyway). None for both real `"milestone": null` and an omitted key.
        "milestone": (content.get("milestone") or {}).get("title"),
        "priority": priority,
        "size": size,
        "target_date": target_date,
        "start_date": start_date,
        "blocked_by": blocked_by,
        "role": role,
        "roles": role_labels,
        "arc": arc,
        "arc_phase": arc_phase,
        "labels": labels,
        # "Not pullable" gate (Issue #1294): a short reason string, or None when
        # the item is pullable. Computed from the one authoritative predicate
        # (pullability.assess) over natively-owned sources - labels, blockedBy,
        # Start date - never the prose scan. Derived, never stored, so it cannot
        # drift from those sources. PRs are always None - forced by the `is_issue`
        # guard below, because a PR is not a Ready-queue candidate anyway (not
        # because the PullRequest fragment omits any of those fields - Start date
        # is a fieldValues entry common to every project item, PR included).
        "not_pullable": (
            assess_pullable(
                {"labels": labels, "blocked_by": blocked_by, "start_date": start_date},
                today=datetime.now(timezone.utc).date().isoformat(),
            )
            if is_issue
            else None
        ),
        "created_at": content.get("createdAt"),
        "updated_at": content.get("updatedAt"),
        "closed_at": content.get("closedAt"),
    }


def sort_key(it: dict[str, Any]) -> tuple:
    return (
        STATUS_ORDER.get(it["status"], 99),
        PRIORITY_ORDER.get(it["priority"], 99),
        SIZE_ORDER.get(it["size"], 99),
        it["number"] or 0,
    )


_MIN_DT = datetime.min.replace(tzinfo=timezone.utc)


def _parse_iso(ts: str | None) -> datetime | None:
    """Parse a GitHub ISO-8601 timestamp (e.g. '2026-06-04T14:15:30Z')."""
    if not ts:
        return None
    return datetime.fromisoformat(ts.replace("Z", "+00:00"))


def age_days(updated_at: str | None, now: datetime) -> float | None:
    """Days since `updated_at`; None when the timestamp is missing/empty."""
    dt = _parse_iso(updated_at)
    if dt is None:
        return None
    return (now - dt).total_seconds() / 86400.0


def age_label(updated_at: str | None, now: datetime) -> str:
    """Compact age-column value, e.g. '3d', or '—' when unknown.

    Clamped at 0 so a server/local clock skew (a `now` slightly behind the
    GitHub `updatedAt`) renders '0d', not a confusing negative.
    """
    a = age_days(updated_at, now)
    return "—" if a is None else f"{max(0, int(a))}d"


def _updated_key(it: dict[str, Any]) -> datetime:
    """Sort key on Issue.updatedAt; missing timestamps sort oldest (epoch min)."""
    return _parse_iso(it.get("updated_at")) or _MIN_DT


def apply_recency(
    items: list[dict[str, Any]],
    *,
    sort_updated: bool,
    stale_days: int | None,
    now: datetime,
) -> list[dict[str, Any]]:
    """Filter/sort open items by last activity (Issue.updatedAt).

    - ``stale_days``: keep only items idle >= N days (dormancy sweep).
    - ``sort_updated``: order most-recently-active first (momentum).

    When both are set the stale filter applies and momentum (descending) ordering
    wins; with only ``stale_days`` the result sorts oldest-active first so the most
    dormant items float to the top. Keyed on content-level ``Issue.updatedAt`` —
    NOT ``ProjectV2Item.updatedAt`` — so a bare board-field nudge does not reset
    staleness (see Issue #642).
    """
    out = list(items)
    if stale_days is not None:
        out = [
            it for it in out
            if (a := age_days(it.get("updated_at"), now)) is not None and a >= stale_days
        ]
    if sort_updated:
        out.sort(key=_updated_key, reverse=True)
    elif stale_days is not None:
        out.sort(key=_updated_key)
    return out


def matches_filter(it: dict[str, Any], args: argparse.Namespace) -> bool:
    if args.role:
        want = args.role if args.role.startswith("role:") else f"role:{args.role}"
        # Membership in the FULL set, not equality with the first role - same reason
        # as the arc filter below (#1153 / #1103). A dual-role Issue must surface in
        # BOTH lanes; matching `it["role"]` alone made it appear under one role and
        # vanish under the other, nondeterministically (GitHub label order is
        # unstable). A Lead could not see the work they were DRI for.
        if want not in it["roles"]:
            return False
    if args.status and it["status"] != args.status:
        return False
    if args.priority and it["priority"] != args.priority:
        return False
    if args.size and it["size"] != args.size:
        return False
    if args.arc:
        want = args.arc if args.arc.startswith("arc:") else f"arc:{args.arc}"
        # Membership in the FULL set, not equality with the first arc. A parent may
        # legitimately carry several arcs (#1103), and GitHub returns labels in an
        # UNSTABLE order - so matching only `it["arc"]` (= arc_labels[0]) would find
        # a multi-arc parent under one of its arcs and silently miss it under the
        # other, nondeterministically. That makes an arc census quietly wrong.
        if want not in it["arcs"]:
            return False
    if args.arc_phase and it["arc_phase"] != args.arc_phase:
        return False
    return True


def format_table(
    items: list[dict[str, Any]],
    now: datetime | None = None,
    arc_columns: bool = False,
) -> str:
    if not items:
        return "(no items matched)\n"
    if now is None:
        now = datetime.now(timezone.utc)
    # Arc/phase columns are opt-in (--arc-columns): the base table is already wide,
    # and arc only matters for an arc-scoped sweep (e.g. after --arc-phase active,
    # where the human table otherwise can't show which arc each issue belongs to —
    # only --json exposes it). Issue #689.
    arc_hdr = f"{'Arc':<18} {'Ph':<7} " if arc_columns else ""
    arc_sep = f"{'-' * 18} {'-' * 7} " if arc_columns else ""
    # Kind column is 8 wide to fit the longest marker "Issue/P" (parent); plain
    # "Issue" / draft "PR/D" are shorter. A narrower budget overflows and shifts
    # every column to its right out from under its header.
    # Ref column (was "#") is 10 wide so a tagged personas ref keeps its column
    # even at 5 digits ("pers#12345") without shifting arc/title; a bare project
    # number is shorter (Issue #999 two-repo disambiguation).
    lines = [
        f"{'Status':<17} {'P':<3} {'Sz':<3} {'Age':<5} {'Role':<17} {'Kind':<8} {'Ref':<10} {arc_hdr}Title",
        f"{'-' * 17} {'-' * 3} {'-' * 3} {'-' * 5} {'-' * 17} {'-' * 8} {'-' * 10} {arc_sep}{'-' * 60}",
    ]
    for it in items:
        # A dual-role Issue (#1153) gets a `+N` marker, mirroring the multi-arc cell
        # below. `--role <r>` matches membership in the FULL set, so without this a
        # row matched via its SECOND role would display only its first and read as
        # though it did not match the filter that returned it.
        role_extra = max(len(it.get("roles") or []) - 1, 0)
        role_suffix = f" +{role_extra}" if role_extra else ""
        role = (it["role"] or "(none)").removeprefix("role:")[: 16 - len(role_suffix)] + role_suffix
        # "/D" = draft PR; "/P" = parent/epic Issue (mutually exclusive — a parent
        # is always an Issue, never a draft).
        kind = it["kind"] + ("/D" if it["is_draft"] else "") + ("/P" if it["is_parent"] else "")
        title = (it["title"] or "")[:60]
        age = age_label(it.get("updated_at"), now)
        if arc_columns:
            # A multi-arc parent (#1103) gets a `+N` marker. `--arc <slug>` matches
            # membership in the FULL set, so without this a row matched via its
            # second arc would display only its first and read as though it did not
            # match the filter that returned it.
            extra = max(len(it.get("arcs") or []) - 1, 0)
            suffix = f" +{extra}" if extra else ""
            arc = (it["arc"] or "—").removeprefix("arc:")[: 17 - len(suffix)] + suffix
            phase = it["arc_phase"] or "—"  # fixed vocab: active/next/later
            arc_cell = f"{arc:<18} {phase:<7} "
        else:
            arc_cell = ""
        lines.append(
            f"{it['status']:<17} {it['priority'] or '?':<3} {it['size'] or '?':<3} "
            f"{age:<5} {role:<17} {kind:<8} {ref_cell(it):<10} {arc_cell}{title}"
        )
    return "\n".join(lines) + "\n"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--role", help='Filter by role label (e.g. "developer", "scientist", "pm")')
    p.add_argument("--status", choices=list(STATUS_ORDER), help="Filter by board Status")
    p.add_argument("--priority", choices=list(PRIORITY_ORDER), help="Filter by Priority")
    p.add_argument("--size", choices=list(SIZE_ORDER), help="Filter by Size")
    p.add_argument("--arc", help='Filter by arc label (e.g. "scoring-tcr-pmhc" or "arc:scoring-tcr-pmhc")')
    p.add_argument("--arc-phase", dest="arc_phase", choices=["active", "next", "later"],
                   help="Filter by arc focus phase")
    p.add_argument("--exclude-parents", dest="exclude_parents", action="store_true",
                   help="Drop parent/epic issues (>=1 sub-issue) from the result. "
                        "Parents carry no Size and mirror a child's Status, so the PM "
                        "triage/flow sweeps flag them as false drift; this filters them "
                        "out at the source (Issue #742).")
    p.add_argument("--check-coherence", dest="check_coherence", action="store_true",
                   help="Coherence sweep (Issue #765): keep only items where a committed "
                        "state (Status past Backlog→Ready, OR a milestone assigned) "
                        "contradicts arc-phase:later (parked). Fails open on items with no "
                        "arc-phase. Forces arc columns so the phase is visible.")
    p.add_argument("--arc-columns", dest="arc_columns", action="store_true",
                   help="Add Arc + phase columns to the table (slug after 'arc:'); "
                        "useful with --arc-phase active to see each issue's arc "
                        "(ignored with --json, which already exposes arc/arc_phase)")
    p.add_argument("--sort-updated", dest="sort_updated", action="store_true",
                   help="Sort by last activity (Issue.updatedAt), most-recent first (momentum)")
    p.add_argument("--stale-days", dest="stale_days", type=int, metavar="N",
                   help="Keep only items idle >= N days, oldest-active first (dormancy sweep; "
                        "N=0 keeps all past-dated items, it does not filter everything out)")
    p.add_argument("--json", dest="as_json", action="store_true", help="Emit JSON array instead of a table")
    args = p.parse_args()

    if args.stale_days is not None and args.stale_days < 0:
        p.error("--stale-days must be a non-negative integer")

    now = datetime.now(timezone.utc)
    raw = fetch_all_items()
    normalized = [n for it in raw if (n := normalize(it)) is not None]
    filtered = [it for it in normalized if matches_filter(it, args)]
    if args.exclude_parents:
        filtered = [it for it in filtered if not it["is_parent"]]
    if args.check_coherence:
        filtered = [it for it in filtered if is_arc_phase_incoherent(it)]
    if args.sort_updated or args.stale_days is not None:
        filtered = apply_recency(
            filtered, sort_updated=args.sort_updated, stale_days=args.stale_days, now=now
        )
    else:
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
        sys.stdout.write(
            format_table(filtered, now=now, arc_columns=args.arc_columns or args.check_coherence)
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
