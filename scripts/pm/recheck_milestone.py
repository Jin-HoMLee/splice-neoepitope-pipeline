#!/usr/bin/env python3
"""Recheck a milestone's due_on against its remaining size-weighted capacity.

Implements the recheck rule from memory/feedback_milestones.md:
  new_due_on = today + remaining-days / 5.0 * 7 days

Usage:
  scripts/pm/recheck_milestone.py --issue N
  scripts/pm/recheck_milestone.py --milestone N
  scripts/pm/recheck_milestone.py --milestone N --moved-issue M

Exits 0 if delta within +-7 days (or no open work), 2 if UPDATE NEEDED / UNSIZED,
3 if the post-move listing is stale (verification pending, Issue #406), 1 on error.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from datetime import date, datetime, timedelta

REPO = "Jin-HoMLee/splice-neoepitope-pipeline"
PROJECT_NUMBER = 9
THRESHOLD_DAYS = 7
AVAILABILITY_RATE = 5.0  # capacity-days per calendar-week

# Post-move eventual-consistency reconciliation (Issue #406). `gh issue edit N
# --milestone X` mutates strongly-consistent issue state, but the
# "issues-by-milestone" listing endpoint (open_issues_in_milestone) lags: right
# after a move the SOURCE milestone can still list #N and the DESTINATION may not
# list it yet, which drove a false-positive [UPDATE NEEDED] on a successful move
# (caught live 2026-05-19). When a moved issue is named we retry the laggy
# listing up to STALE_RETRY_ATTEMPTS times until it agrees with the strong read,
# then bail to STALE_STATUS rather than emit a misleading capacity read.
STALE_RETRY_ATTEMPTS = 2
STALE_RETRY_DELAY_SECONDS = 0.5
STALE_STATUS = "[stale state, verification pending]"
STALE_EXIT_CODE = 3

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


GH_MAX_ATTEMPTS = 4
assert GH_MAX_ATTEMPTS >= 1, "gh() runs the loop at least once; a terminal raise needs a result"
GH_BACKOFF_BASE_SECONDS = 2.0
# Clamp a Retry-After hint: GitHub can legally emit a large value (e.g. 3600s)
# during sustained degradation, which would otherwise stall the nightly live job
# for that whole duration before raising (Issue #711 review).
GH_RETRY_AFTER_CAP_SECONDS = 120.0

# Deterministic client errors that won't fix themselves on retry — a bad request
# stays bad. Everything else (5xx, 403/secondary-rate-limit, network/timeout, an
# empty-stderr crash) is treated as transient and retried. We fail TOWARD retrying:
# the goal is de-flaking the live recheck smoke against transient GitHub-API blips
# (Issue #711 D4), and over-retrying a genuine bad arg only wastes a bounded few
# seconds, whereas under-retrying a transient blip reds the whole job.
_DETERMINISTIC_HTTP_RE = re.compile(r"HTTP (400|401|404|410|422)\b")
_RETRY_AFTER_RE = re.compile(r"retry[- ]after[:\s]+(\d+)", re.IGNORECASE)


def _is_transient_gh_error(stderr: str) -> bool:
    return not _DETERMINISTIC_HTTP_RE.search(stderr or "")


def _retry_after_seconds(stderr: str) -> float | None:
    m = _RETRY_AFTER_RE.search(stderr or "")
    return float(m.group(1)) if m else None


def gh(*args: str, parse_json: bool = True, _runner=subprocess.run, _sleep=time.sleep) -> object:
    """Run ``gh`` with retry + exponential backoff on transient failures (Issue #711 D4).

    GitHub secondary rate limits and transient 5xx / replication-lag errors surface
    as a non-zero ``gh`` exit unrelated to the request. The live recheck smoke makes
    ~35-40 calls per CI run, so a single such blip would otherwise red the shared
    ``ci-tools-pytest`` job. Transient non-zero exits are retried up to
    ``GH_MAX_ATTEMPTS`` with exponential backoff, honoring a ``Retry-After`` hint
    when GitHub provides one. A deterministic 4xx is not retried, and a terminal
    failure raises ``CalledProcessError`` — preserving the ``check=True`` contract
    callers depend on. ``_runner`` / ``_sleep`` are injection seams for tests.
    """
    cmd = ["gh", *args]
    result = None
    for attempt in range(GH_MAX_ATTEMPTS):
        result = _runner(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            return json.loads(result.stdout) if parse_json else result.stdout
        if not _is_transient_gh_error(result.stderr):
            break
        if attempt < GH_MAX_ATTEMPTS - 1:
            delay = _retry_after_seconds(result.stderr)
            if delay is None:
                delay = GH_BACKOFF_BASE_SECONDS * (2 ** attempt)
            _sleep(min(delay, GH_RETRY_AFTER_CAP_SECONDS))
    raise subprocess.CalledProcessError(
        result.returncode, cmd, output=result.stdout, stderr=result.stderr
    )


def milestone_for_issue(issue_number: int) -> int | None:
    data = gh("issue", "view", str(issue_number), "--repo", REPO, "--json", "milestone")
    return data["milestone"]["number"] if data.get("milestone") else None


def milestone_and_state_for_issue(issue_number: int) -> tuple[int | None, str]:
    """Strongly-consistent ``(milestone_number, state)`` for an issue (Issue #406).

    Folds the milestone read the reconciler needs with the issue ``state`` in one
    ``gh issue view`` call. ``state`` is gh's ``"OPEN"``/``"CLOSED"``. Used so the
    reconciler can short-circuit a closed moved issue, which never appears in the
    open-only milestone listing and so could otherwise loop forever on a false
    ``[stale state]`` (PR #985 review, finding 1).
    """
    data = gh(
        "issue", "view", str(issue_number), "--repo", REPO, "--json", "milestone,state"
    )
    ms = data["milestone"]["number"] if data.get("milestone") else None
    return ms, data.get("state", "")


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


def _reconcile_moved_issue(
    milestone_number: int,
    milestone_title: str,
    moved_issue: int,
    issue_numbers: list[int],
    _sleep=time.sleep,
) -> tuple[list[int], bool]:
    """Reconcile the laggy milestone listing against a strongly-consistent read
    of ``moved_issue``'s actual milestone before capacity is computed (Issue #406).

    ``milestone_and_state_for_issue`` is backed by ``gh issue view``, which is
    strongly consistent, so it is the source of truth for whether ``moved_issue``
    belongs to THIS milestone. We compare that against the eventually-consistent
    listing and, on disagreement, re-fetch the listing up to
    ``STALE_RETRY_ATTEMPTS`` times (``STALE_RETRY_DELAY_SECONDS`` apart) until it
    converges.

    A CLOSED moved issue short-circuits to converged: it never appears in the
    open-only listing and does not count toward capacity, so there is nothing to
    reconcile - without this guard a closed issue moved INTO a milestone reads as
    should-be-member yet is always absent, looping forever on a false ``[stale state]``
    (PR #985 review, finding 1).

    Returns ``(issue_numbers, stale)``: the (possibly re-fetched) listing and a
    flag that is True when the listing never caught up to the strong read.
    ``_sleep`` is an injection seam for tests.
    """
    actual_ms, state = milestone_and_state_for_issue(moved_issue)
    if state == "CLOSED":
        return issue_numbers, False
    should_be_member = actual_ms == milestone_number
    for attempt in range(STALE_RETRY_ATTEMPTS + 1):
        if (moved_issue in issue_numbers) == should_be_member:
            return issue_numbers, False
        if attempt < STALE_RETRY_ATTEMPTS:
            _sleep(STALE_RETRY_DELAY_SECONDS)
            issue_numbers = open_issues_in_milestone(milestone_title)
    return issue_numbers, True


def sizes_and_parents_for_issues(
    issue_numbers: list[int],
) -> tuple[dict[int, str | None], set[int]]:
    """One GraphQL round-trip returning ``(sizes, parents)`` for ``issue_numbers``.

    Folds what used to be two separate aliased GraphQL queries (``sizes_for_issues``
    + ``parent_numbers``) into a single call: each ``i{n}`` alias now fetches both
    the project Size field value AND ``subIssuesSummary.total`` (Issue #711). This
    halves the GraphQL surface per recheck — material because the live smoke test
    rechecks every open milestone, and one transient ``gh`` non-zero on any call
    reds the shared ``ci-tools-pytest`` job.

    - ``sizes[n]`` — the issue's Size on PROJECT_NUMBER (``None`` if unset / on a
      different board).
    - ``parents`` — the subset that are parent epics (``subIssuesSummary.total > 0``).
      Parents carry no Size by convention — size rolls up from their sub-issues
      (shared/feedback_parent_sub_issues.md). compute_recheck excludes them from
      both the capacity sum and the unsized-check, so a milestone holding a parent
      as a roadmap anchor doesn't emit a spurious, un-clearable [UNSIZED] flag
      (Issue #689 — sizing the parent to clear it would itself violate the
      no-size convention).
    """
    if not issue_numbers:
        return {}, set()
    owner, name = REPO.split("/")
    aliases = " ".join(
        f'i{n}: issue(number: {n}) {{ '
        f'subIssuesSummary {{ total }} '
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
    parents: set[int] = set()
    for n in issue_numbers:
        node = repo.get(f"i{n}") or {}
        total = (node.get("subIssuesSummary") or {}).get("total") or 0
        if total > 0:
            parents.add(n)
        size: str | None = None
        for pi in node.get("projectItems", {}).get("nodes", []):
            if (pi.get("project") or {}).get("number") != PROJECT_NUMBER:
                continue
            for fv in pi.get("fieldValues", {}).get("nodes", []):
                if (fv.get("field") or {}).get("name") == "Size":
                    size = fv.get("name")
        sizes[n] = size
    return sizes, parents


def compute_recheck(
    milestone_number: int,
    moved_issue: int | None = None,
    _sleep=time.sleep,
) -> int:
    meta = milestone_meta(milestone_number)
    title = meta["title"]
    current_due_raw = meta.get("due_on")
    current_due_date = (
        datetime.fromisoformat(current_due_raw.replace("Z", "+00:00")).date()
        if current_due_raw else None
    )

    issue_numbers = open_issues_in_milestone(title)

    # Post-move reconciliation (Issue #406): only when a moved issue is named
    # (the `gh issue edit --milestone` trigger). If the listing endpoint hasn't
    # caught up to the strongly-consistent move, bail with a [stale state] note
    # instead of computing a misleading capacity read off stale membership.
    if moved_issue is not None:
        issue_numbers, stale = _reconcile_moved_issue(
            milestone_number, title, moved_issue, issue_numbers, _sleep=_sleep
        )
        if stale:
            print(f"Milestone: {title}")
            print(f"Current due_on: {current_due_date or '(none)'}")
            print(
                f"Open issues: (listing endpoint has not caught up to the move "
                f"of #{moved_issue})"
            )
            print("Proposed due_on: (skipped - post-move listing is stale)")
            print(
                f"Status: {STALE_STATUS} - re-run recheck once GitHub's "
                f"issues-by-milestone listing catches up"
            )
            return STALE_EXIT_CODE

    sizes, parents = sizes_and_parents_for_issues(issue_numbers)

    print(f"Milestone: {title}")
    print(f"Current due_on: {current_due_date or '(none)'}")
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
    parser.add_argument(
        "--moved-issue", type=int, default=None,
        help="Issue just moved via `gh issue edit --milestone`; enables post-move "
             "eventual-consistency reconciliation of the listing endpoint (Issue #406)",
    )
    args = parser.parse_args()

    if args.issue:
        ms = milestone_for_issue(args.issue)
        if ms is None:
            print(f"error: issue #{args.issue} has no milestone", file=sys.stderr)
            return 1
        return compute_recheck(ms, moved_issue=args.moved_issue)
    return compute_recheck(args.milestone, moved_issue=args.moved_issue)


if __name__ == "__main__":
    sys.exit(main())
