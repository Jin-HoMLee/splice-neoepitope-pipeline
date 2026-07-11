#!/usr/bin/env python3
"""Recheck parent-vs-children Status drift on project board #9.

Implements the audit rule from docs/superpowers/specs/2026-05-19-parent-status-drift-audit-design.md.

Usage:
  scripts/pm/recheck_parent_status.py --issue N      # recheck this issue's parent chain
  scripts/pm/recheck_parent_status.py --all          # audit all parent issues on project #9

Exits 0 if no drift, 2 if drift detected, 1 on error.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Import the shared hardened gh() wrapper from the sibling module (Issue #1017).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from gh_client import GhError, gh  # noqa: E402

# Single-source the AC-box scan from the closure gates (Issue #1067). Re-implementing
# checkbox parsing here would let this flag and the merge/close gates drift on what
# counts as an AC box - the exact drift the `## Acceptance criteria` convention and
# the #730 stray-box lint exist to prevent.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools" / "ci"))
from closure_audit import scan_ac_boxes  # noqa: E402

REPO = "Jin-HoMLee/splice-neoepitope-pipeline"
PROJECT_NUMBER = 9

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

# The dedicated parked Status for parents/epics under the A2 epic-park model
# (Issue #776). It is deliberately OFF the workflow ladder above — a parent in
# `Epic` is not making a leaf-ladder progress claim, so drift classification must
# not compare it against children's collective rank (Issue #794).
EPIC_STATUS = "Epic"

# Softer flag emitted in place of a bare COMPLETION DRIFT when a parent's only
# remaining children are all closed but ≥1 was closed NOT_PLANNED — that scope
# was deferred, not delivered, so a confident rollup-close suggestion would be
# wrong (Issue #632). format_record() owns the surrounding brackets, so this
# constant carries no brackets of its own.
NOT_PLANNED_REVIEW = (
    "REVIEW: parent has a not-planned child — verify scope was delivered, not deferred"
)

# Issue #1067. Under the A2 park a parent's progress is read off GitHub's native
# sub-issue bar, so a **full bar** reads as "done" at a glance. But body scope
# legitimately exceeds the filed sub-issues (design/docs children close first while
# implementation ACs, or references to sibling Issues, remain), so a full bar with
# unticked body ACs is a real and *invisible* drift class: the bar says complete,
# the body says otherwise. Live cases at filing: #527 (7/7 native, AC gated on a
# then-open sibling), #859 (2/2 native, 4 unticked implementation ACs), #665 (2/2
# native, 2 unticked "proper-fix" ACs). Advisory only, mirroring the #632 softening.
#
# Keyed on ANY unticked box in the parent body, not only boxes under a canonical
# `## Acceptance criteria` heading - a correction forced by running it against the
# live cases before shipping. #859's four unticked boxes sit under `## Sub-issues`,
# so an AC-only flag returns **zero on its own motivating example**: a no-op that
# would have looked green. And that is not an authoring mistake to be linted away:
# a parent body IS a roadmap, so an unticked `## Sub-issues` line is *precisely*
# the "body scope exceeds the filed sub-issues" signal this flag exists to catch.
# The message names the heading, so a human can tell a real gating box from an
# exploratory plan step. (The #730 stray-box lint stays the merge-time authority on
# *where* boxes belong; this is a different question - does the bar tell the truth.)
BODY_AC_REVIEW = (
    "REVIEW: native sub-bar is full but {n} body box(es) still unticked ({where}) - "
    "body scope exceeds the filed sub-issues; verify before closing"
)


def collective_state(open_children: list[dict]) -> str:
    """Max-rank status across open children. Empty list → 'Done'."""
    if not open_children:
        return "Done"
    max_rank = max(rank(c.get("status")) for c in open_children)
    return LADDER_INVERSE[max_rank]


def classify_drift(
    parent_status: str | None,
    open_children: list[dict],
    has_not_planned: bool = False,
    body_unticked_count: int = 0,
    body_unticked_where: str = "",
    has_known_children: bool = True,
) -> str | None:
    """Classify drift for a parent vs its open children.

    Returns one of: 'FORWARD DRIFT', 'BACKWARD DRIFT', 'COMPLETION DRIFT',
    NOT_PLANNED_REVIEW, BODY_AC_REVIEW, or None.

    ``has_not_planned`` is True when the parent has ≥1 sub-issue closed
    NOT_PLANNED. It only matters when every child is closed (``open_children``
    empty): a NOT_PLANNED child carries deferred scope, so "all closed" does
    not imply "complete" — emit the softer verify prompt instead of a confident
    rollup-close suggestion (Issue #632).

    ``body_unticked_count`` is the count of unticked boxes anywhere in the parent
    body, AC and non-AC alike (Issue #1067); ``body_unticked_where`` names their
    headings. It matters only when every child is closed: a full native sub-bar
    reads as "done", but body scope routinely exceeds the filed sub-issues, so any
    unticked body box means the bar is lying. Emits a verify prompt rather than the
    confident rollup-close.

    ``has_known_children`` is False when the sub-issue list came back empty, i.e.
    we cannot actually see any children (a cross-repo parent's sub-issues are
    invisible to this repo's API). "All children closed" and "no children visible"
    are indistinguishable from an empty list, so the body-AC prompt is suppressed
    there rather than false-flagged on a bar we never actually read (AC 3).
    """
    p_rank = rank(parent_status)
    if not open_children:
        # All children closed; parent should be Done. Under A2 a completed parent
        # still closes → Done, so this close-the-parent signal is preserved for an
        # Epic-parked parent too (it is not a leaf-status mirror).
        if p_rank == STATUS_LADDER["Done"]:
            return None
        # A NOT_PLANNED child means the completion claim is unverified; downgrade
        # the bare COMPLETION DRIFT to a verify prompt (Issue #632). It outranks the
        # body-AC prompt: deferred scope is the stronger reason not to trust "done".
        if has_not_planned:
            return NOT_PLANNED_REVIEW
        if body_unticked_count > 0 and has_known_children:
            return BODY_AC_REVIEW.format(
                n=body_unticked_count, where=body_unticked_where
            )
        return "COMPLETION DRIFT"
    if parent_status == EPIC_STATUS:
        # A2 epic-park (#776 / #794): a parent parked in the off-ladder `Epic`
        # Status no longer mirrors its children's collective ladder rank — progress
        # is read off GitHub's native sub-issue bar instead. Suppress the ladder
        # mirror so the park is not fought: `Epic` is unranked (rank 0), so without
        # this guard any active child would spuriously read as BACKWARD DRIFT.
        return None
    c_rank = rank(collective_state(open_children))
    if p_rank > c_rank and p_rank >= STATUS_LADDER["In progress"]:
        # Only flag when the parent is making a falsifiable progress claim (In progress
        # or beyond). A Ready parent with Backlog children is the normal post-grooming
        # state and should not be treated as drift.
        return "FORWARD DRIFT"
    if p_rank < c_rank:
        return "BACKWARD DRIFT"
    return None


def parent_issue_number(issue_number: int) -> int | None:
    """Return parent issue number via REST parent_issue_url, or None."""
    data = gh("api", f"repos/{REPO}/issues/{issue_number}")
    url = data.get("parent_issue_url")
    if not url:
        return None
    return int(url.rstrip("/").rsplit("/", 1)[-1])


def all_sub_issues(issue_number: int) -> list[dict]:
    """Return all sub-issues (open and closed), each carrying 'state' and
    'state_reason' from the REST /sub_issues endpoint."""
    return gh("api", f"repos/{REPO}/issues/{issue_number}/sub_issues")


def open_sub_issues(issue_number: int) -> list[dict]:
    """Return list of open sub-issues (each a dict with 'number' at minimum)."""
    return [c for c in all_sub_issues(issue_number) if c.get("state") == "open"]


def has_not_planned_child(issue_number: int) -> bool:
    """True if any sub-issue was closed with state_reason 'not_planned'.

    open_sub_issues() discards closed children, so this is the only path that
    inspects close reasons — required to tell a deferred (NOT_PLANNED) close
    from a delivered (COMPLETED) one when a parent is all-closed (Issue #632).
    """
    return any(
        c.get("state") == "closed" and c.get("state_reason") == "not_planned"
        for c in all_sub_issues(issue_number)
    )


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


def body_unticked(issue_number: int) -> tuple[int, str]:
    """(count, where) of ALL unticked boxes in the issue body: AC + non-AC alike.

    Single-sources `closure_audit.scan_ac_boxes` rather than re-implementing
    checkbox parsing, so the two agree on what a checkbox is. But it sums BOTH
    partitions, because the question here is not "are the ACs done" - it is "does
    the full native sub-bar tell the truth". A parent body is a roadmap, so an
    unticked `## Sub-issues` line is exactly the body-scope-exceeds-children signal
    (see BODY_AC_REVIEW). `where` names the headings so a human can distinguish a
    real gating box from an exploratory plan step.

    Returns (0, "") on any gh/parse error - advisory only, never break the sweep.
    """
    try:
        data = gh("api", f"repos/{REPO}/issues/{issue_number}")
        scan = scan_ac_boxes(data.get("body") or "")
    except Exception:  # noqa: BLE001 - advisory only; never break the sweep
        return 0, ""
    total = scan.ac_unticked + scan.stray_unticked
    if not total:
        return 0, ""
    where: list[str] = []
    if scan.ac_unticked:
        where.append(f"{scan.ac_unticked} under 'Acceptance criteria'")
    if scan.stray_unticked:
        headings = ", ".join(f"'{h}'" for h in scan.stray_headings)
        where.append(f"{scan.stray_unticked} under {headings}")
    return total, "; ".join(where)


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
        # Only consult close reasons + the body ACs when every child is closed —
        # that's the one case classify_drift uses them, and it saves API calls
        # otherwise.
        all_closed = not enriched
        has_np = all_closed and has_not_planned_child(cursor)
        # `has_known_children` distinguishes "the bar is full" from "we cannot see
        # the bar" (a cross-repo parent's sub-issues are invisible here). Both
        # arrive as an empty open-children list, so without this the body-AC prompt
        # would fire on a parent whose children we never actually read (AC 3).
        known_children = bool(all_sub_issues(cursor)) if all_closed else True
        unticked, where = (
            body_unticked(cursor) if all_closed and known_children else (0, "")
        )
        chain.append({
            "issue": cursor,
            "status": parent_status,
            "open_children": enriched,
            "collective": collective_state(enriched),
            "drift": classify_drift(
                parent_status, enriched,
                has_not_planned=has_np,
                body_unticked_count=unticked,
                body_unticked_where=where,
                has_known_children=known_children,
            ),
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


def run_issue_mode(issue_number: int) -> int:
    # Deliberate asymmetry with run_all_mode (Issue #1017 review): the single-chain
    # --issue path has no per-item isolation - a terminal GhError propagates and
    # fails hard (after the shared gh()'s retries). A one-target lookup should fail
    # loudly rather than half-report; only the --all sweep isolates per parent.
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

    if args.issue is not None:
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
        try:
            meta = gh("api", f"repos/{REPO}/issues/{n}")
        except GhError as e:
            # Per-item isolation (Issue #1017): a terminal gh failure on one issue's
            # parent-probe must not abort the whole listing. Skip it and continue.
            print(f"  [skipped #{n}: gh error probing sub-issues - {e}]", file=sys.stderr)
            continue
        if (meta.get("sub_issues_summary") or {}).get("total", 0) > 0:
            parents.append(n)
    return parents


def run_all_mode() -> int:
    parents = all_parent_issues()
    drifted_count = 0
    skipped_count = 0
    drift_blocks: list[str] = []
    for p in parents:
        try:
            parent_status = status_for_issue(p)
            children = open_sub_issues(p)
            enriched = [{"number": c["number"], "status": status_for_issue(c["number"])}
                        for c in children]
            has_np = not enriched and has_not_planned_child(p)
        except GhError as e:
            # Per-item isolation (Issue #1017): one parent's terminal gh failure
            # skips that parent, it does not abort the remaining audit.
            skipped_count += 1
            print(f"  [skipped #{p}: gh error during audit - {e}]", file=sys.stderr)
            continue
        record = {
            "issue": p,
            "status": parent_status,
            "open_children": enriched,
            "collective": collective_state(enriched),
            "drift": classify_drift(parent_status, enriched, has_not_planned=has_np),
        }
        if record["drift"] is not None:
            drifted_count += 1
            drift_blocks.append(format_record(record))

    suffix = f" ({skipped_count} skipped on gh error)" if skipped_count else ""
    print(f"Audited {len(parents)} parent issues; {drifted_count} drifted{suffix}.\n")
    for block in drift_blocks:
        print(block)
        print()
    # A sweep that skipped every parent on gh errors audited nothing - exit with the
    # error code so a 0 can never be misread as "clean board" when it actually means
    # "blind" (Issue #1017 review: the exact silent-miss class this tool exists to
    # prevent). A partial skip keeps the 0/2 signal, surfaced by the suffix above.
    if parents and skipped_count == len(parents):
        return 1
    return 2 if drifted_count > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
