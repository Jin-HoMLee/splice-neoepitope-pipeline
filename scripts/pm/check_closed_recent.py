#!/usr/bin/env python3
"""Recap of recently closed Issues + merged PRs (morning-routine Beat 1a/1b).

The Service Delivery Review closure-recap was a hand-typed query each morning:
`gh issue list --state closed --search "closed:>YYYY-MM-DD"`. The `>` operator
excludes the whole boundary day, so `closed:>2026-06-18` returns only items closed
06-19 onward - silently reporting "nothing closed" when items merged *during* the
boundary day. That off-by-one slipped 2026-06-19 (caught by the user). This script
makes the boundary deterministic (`>=`) and covers both closed Issues AND merged
PRs in one place - the closed-side companion to `board_open_items.py`'s open scan.

Scans BOTH the project and personas repos by default (Issue #1276), so the
Memory Manager workstream is visible rather than silently missing. Personas rows
render with a `pers#` prefix because issue numbers COLLIDE across the two repos.

Usage:

  scripts/pm/check_closed_recent.py                 # default: since the start of yesterday (UTC), both repos
  scripts/pm/check_closed_recent.py --days 3        # since the start of the day 3 days ago
  scripts/pm/check_closed_recent.py --since 2026-07-01
  scripts/pm/check_closed_recent.py --repo Jin-HoMLee/splice-neoepitope-pipeline   # narrow to one repo

Window floor is a whole day (00:00 UTC): the recap is day-granular, matching how
GitHub's `closed:` search operator works. `--since YYYY-MM-DD` / `--days N`
override the default.

Exit 0 on a completed scan - finding nothing closed is not a failure. A hard `gh`
/ network error on a listing call surfaces loudly (non-zero traceback) rather than
as a false "nothing closed". Issue #784.
"""
import argparse
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gh_client import gh  # noqa: E402

PROJECT_REPO = "Jin-HoMLee/splice-neoepitope-pipeline"
PERSONAS_REPO = "Jin-HoMLee/claude-personas-splice-neoepitope-pipeline"

# The recap covers BOTH repos by default (Issue #1276). `REPO` was hard-coded to
# the project repo, so the entire Memory Manager workstream - which lives in the
# personas repo - was invisible to the mechanized recap, and a two-day column of
# MM closes was silently omitted on 2026-07-22 until someone asked why MM work
# never appeared. A recap that cannot see a whole workstream reports "nothing
# closed" for it, which is the same silent-wrong-answer class this script was
# written (Issue #784) to prevent.
DEFAULT_REPOS = (PROJECT_REPO, PERSONAS_REPO)

# Backwards-compatible alias: some callers/tests referenced the old constant.
REPO = PROJECT_REPO

UTC = timezone.utc


# --- pure helpers (unit-tested, no I/O) ---


def repo_ref(number, repo):
    """Render an item reference that is unambiguous ACROSS the two repos.

    Board 9 aggregates two repos whose issue numbers COLLIDE (both have a #98,
    and they are unrelated Issues). Once the recap spans both repos, a bare
    `#98` silently means one of two things - so covering both repos without
    disambiguating would trade a missing-data bug for a wrong-attribution bug,
    which is worse because it looks correct.

    Mirrors the `pers#` prefix `board_open_items.py` already uses, so the two
    tools read the same way rather than inventing a second convention. The
    personas repo is matched first because its name CONTAINS the project repo's
    name as a substring.
    """
    if repo == PERSONAS_REPO or "claude-personas" in (repo or ""):
        return f"pers#{number}"
    return f"#{number}"


def sort_rows(rows):
    """Newest-first, tolerating an unparseable timestamp.

    Sorts on the parsed timestamp, not the raw string: a fail-open placeholder
    ("?") sorts above digits and would jump a malformed row to the top of a
    newest-first list. Unparseable rows sink to the bottom instead.
    """
    return sorted(
        rows,
        key=lambda r: parse_ts(r["when"]) or datetime.min.replace(tzinfo=UTC),
        reverse=True,
    )


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


def filter_merged(rows):
    """Keep only merged PRs (mergedAt present); drop closed-unmerged (mergedAt null).

    Pulled out of the `gh` I/O layer so the merged-vs-closed selection rule - the
    one real business rule on the PR side - is unit-testable without a live `gh`.
    """
    return [r for r in rows if r.get("mergedAt")]


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


def list_closed_issues(search, repo=PROJECT_REPO):
    """Closed Issues matching `search` (number, title, closedAt).

    Goes through the shared `gh_client.gh` wrapper (Issue #1062), the last of the
    hand-rolled `gh` helpers to converge. The local `gh_json()` was bare
    `check=True` with no retry, so a transient 5xx / secondary-rate-limit /
    replication blip raised - and in a *recap* that reads as "nothing closed",
    which is the silent-wrong-answer this script already exists to prevent
    (Issue #784). The wrapper retries transient failures with exponential backoff
    and honors Retry-After, while still raising GhError on a terminal failure, so
    the check=True contract callers rely on is preserved.
    """
    return gh(
        "issue", "list", "--repo", repo, "--state", "closed",
        "--search", search, "--limit", "1000",
        "--json", "number,title,closedAt",
    )


def list_merged_prs(search, repo=PROJECT_REPO):
    """PRs matching `search`, filtered to the *merged* ones (mergedAt present).

    `--state all --search "closed:>=…"` returns both merged and closed-unmerged
    PRs; only merged PRs belong in a shipped-work recap, so the unmerged ones
    (mergedAt is null) are dropped.
    """
    rows = gh(
        "pr", "list", "--repo", repo, "--state", "all",
        "--search", search, "--limit", "1000",
        "--json", "number,title,mergedAt",
    )
    return filter_merged(rows)


# --- orchestration ---


def collect(floor, issues, prs, repo=PROJECT_REPO):
    """Merge one repo's closed Issues + merged PRs into recap rows.

    `issues` / `prs` are the raw `gh` result lists (injected, so this is testable
    without network). Each is re-filtered client-side via `is_within` so the
    boundary is enforced here too, not only by the search string.

    Single-repo by design: each row is tagged with the `repo` it came from, and
    the caller concatenates across repos then calls `sort_rows` once. Sorting
    here would only order within a repo, which would interleave wrongly the
    moment a second repo's rows were appended.
    """
    rows = []
    for it in issues:
        if is_within(it.get("closedAt"), floor):
            rows.append({"number": it["number"], "kind": "Issue", "repo": repo,
                         "when": it.get("closedAt") or "?", "title": it.get("title") or ""})
    for pr in prs:
        if is_within(pr.get("mergedAt"), floor):
            rows.append({"number": pr["number"], "kind": "PR", "repo": repo,
                         "when": pr.get("mergedAt") or "?", "title": pr.get("title") or ""})
    return sort_rows(rows)


def render(rows, floor):
    header = f"Closed Issues + merged PRs since {floor_date_str(floor)} (00:00 UTC):"
    if not rows:
        return header + "\n  (nothing closed in window)\n"
    lines = [header]
    for r in rows:
        when = r["when"][:16].replace("T", " ")
        ref = repo_ref(r["number"], r.get("repo", PROJECT_REPO))
        lines.append(f"  {ref:<9} {r['kind']:<5} {when}  {r['title']}")
    lines.append(f"\n  {len(rows)} item(s).")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--since", help="window floor YYYY-MM-DD (overrides --days)")
    group.add_argument("--days", type=int, default=1,
                       help="window = since the start of the day N days ago (default 1 = yesterday)")
    parser.add_argument(
        "--repo", action="append", metavar="OWNER/NAME",
        help=("repo to scan; repeatable. Default scans BOTH the project and "
              "personas repos, so the Memory Manager workstream is visible "
              "(Issue #1276). Pass explicitly to narrow."),
    )
    args = parser.parse_args()
    repos = tuple(args.repo) if args.repo else DEFAULT_REPOS

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
    # Fetch per repo, tag the rows, then sort ONCE across the merged set - so the
    # recap reads as a single newest-first timeline rather than two concatenated
    # per-repo blocks.
    rows = []
    for repo in repos:
        rows += collect(floor,
                        list_closed_issues(search, repo),
                        list_merged_prs(search, repo),
                        repo=repo)
    print(render(sort_rows(rows), floor), end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
