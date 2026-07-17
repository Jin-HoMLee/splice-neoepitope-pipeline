#!/usr/bin/env python3
"""Reconcile prose-stated Issue dependencies into native GitHub blockedBy edges.

Finds bodies saying "depends on #M" / "blocked on #M" etc., diffs them against
the native blockedBy graph, and reports or wires the missing edges.

Usage:
  scripts/pm/scan_prose_deps.py                 # --report (default): drift table
  scripts/pm/scan_prose_deps.py --issue N       # single-issue scope
  scripts/pm/scan_prose_deps.py --check         # exit 1 if any drift
  scripts/pm/scan_prose_deps.py --apply         # wire all needs-wiring records
  scripts/pm/scan_prose_deps.py --apply --only 745 594   # wire only these dependents

Exit-code convention (shared with scan_unblocked.py, aligned in Issue #1179):
0 = ran fine, nothing to act on; 1 = ran fine, --check gate tripped (drift
present); 2 = could-not-run (a transient per-issue lookup failure / incomplete
scan, which takes precedence over 1). 2 == could-not-run matches argparse's own
exit-2-for-usage-error, so a routine shelling several PM sweeps can treat the
codes uniformly: 1 = a finding, 2 = did not execute.
"""
import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gh_client import GhError, gh  # noqa: E402

REPO = "Jin-HoMLee/splice-neoepitope-pipeline"
REPO_OWNER, REPO_NAME = REPO.split("/")


class IssueLookupError(Exception):
    """A per-issue `gh` lookup could not be completed (a transient 5xx / network
    blip, a GraphQL `errors` payload, or non-JSON stdout). The caller skips that
    one pair and surfaces it, rather than aborting the whole board scan or
    misreading a failed lookup as a real result.

    Abstract: only the concrete subclasses are ever raised. Each sets a concrete
    `action` report-row label AND registers it in `_LOOKUP_FAILED_ACTIONS` below
    (the single source the `--check` gate + report ordering key off). `action`
    is None here so a subclass that forgets to override fails loudly, rather than
    inheriting a plausible-looking string that would silently escape the gate."""
    action = None


class BlockerLookupError(IssueLookupError):
    """The native blockedBy lookup for one issue could not be completed (#989)."""
    action = "edges-lookup-failed"


class MetaLookupError(IssueLookupError):
    """The blocker-meta (state / is_pr) REST lookup could not be completed (#1012)."""
    action = "meta-lookup-failed"


# Single source of truth for the scan-integrity actions: derived from the classes,
# NOT a string convention. Both the `--check` incomplete-scan gate and the report
# sort order below consume this, so a new lookup type is registered in exactly one
# place and cannot silently escape either (edges: #989; meta + this guard: #1012).
_LOOKUP_FAILED_ACTIONS = (BlockerLookupError.action, MetaLookupError.action)


def fetch_open_issues():
    """All open issues with bodies. Uses the issue list (NOT the project board),
    so the board's Done-first pagination trap does not apply."""
    return gh(
        "issue", "list", "--repo", REPO, "--state", "open", "--limit", "1000",
        "--json", "number,title,body,state",
    )


# Blocker-phrase allowlist — deliberately narrow: only verbs that mean "this is
# blocked until #M". Narrative verbs (informs/consumes/relates to/see/fixes/...)
# are NOT listed, so they never match — that IS the exclude mechanism (a narrow
# allowlist), proven by test_parse_narrative_phrases_excluded.
# Note: blocked-by: allows optional # (e.g. "blocked-by:722" or "blocked-by:#722");
# others require # to avoid false matches like "requires 5 steps".
BLOCKER_RE = re.compile(
    r"(?:"
    r"(?:depends on|blocked by|blocked on|gated on|requires)\s*#"
    r"|blocked-by:?\s*#?"
    r")(\d+)",
    re.IGNORECASE,
)


def _normalize(text):
    """Strip markdown that breaks blocker-phrase -> #N adjacency.

    Unwraps issue-ref markdown links ([Issue #719](url) -> #719) and removes
    emphasis/code markers (* `) so a bolded or linked dependency still matches
    BLOCKER_RE. Underscores are left intact (they appear in identifiers)."""
    text = re.sub(r"\[[^\]]*?#(\d+)[^\]]*?\]\([^)]*\)", r"#\1", text)  # unwrap links first
    text = re.sub(r"[*`]", "", text)                                    # then strip emphasis
    return text


def parse_dependencies(number, body):
    """Return sorted, deduped (dependent, blocker) pairs from one issue body.

    Pure: no I/O. Normalizes markdown (emphasis + issue-ref links) then matches
    the blocker-phrase allowlist with strict adjacency; drops self-references.
    KNOWN LIMITATION: a multi-word narrative gap between the phrase and #N
    (e.g. "depends on the registry from #732") is NOT caught — strict adjacency
    is deliberate to avoid false-matching tool/hardware deps ("requires
    NetMHCpan-4.0"). Such cases are surfaced via the human review gate, not here.
    Two more exotic gaps: a link title with multiple issue refs (e.g. "[#5 and #6](url)")
    keeps only the first #N after link-unwrap; and a blocker phrase embedded inside a
    link title (e.g. "[depends on #722](url)") is stripped by link-unwrap and missed.
    Blocker open/closed and PR-vs-issue filtering happen later, in reconcile()."""
    if not body:
        return []
    blockers = {int(m.group(1)) for m in BLOCKER_RE.finditer(_normalize(body))}
    blockers.discard(number)  # self-reference
    return [(number, b) for b in sorted(blockers)]


def classify(dependent, blocker, *, blocker_meta, existing):
    """Classify one (dependent, blocker) pair. Pure given resolved inputs."""
    if dependent == blocker:
        return "self-ref"  # unreachable via reconcile (parse drops self-refs); defensive
    if blocker in existing:
        return "already-wired"
    if blocker_meta["state"] == "closed":
        return "closed-blocker"  # convention: only wire currently-open blockers
    if blocker_meta["is_pr"]:
        return "un-wireable-pr"  # native deps are issue<->issue
    return "needs-wiring"


def issue_meta(number):
    """{'state': 'open'|'closed', 'is_pr': bool} via the REST issues endpoint
    (which serves both issues and PRs; a PR carries a 'pull_request' key).

    Raises MetaLookupError if the lookup could not be completed (a failed `gh`
    call, or non-JSON stdout), so the caller skips that one pair instead of
    aborting the whole scan on a transient blip - the same resilience #989 gave
    the blockedBy lookup, extended to the blocker-meta fetch (#1012)."""
    try:
        obj = gh("api", f"repos/{REPO}/issues/{number}")
    except (GhError, json.JSONDecodeError) as e:
        raise MetaLookupError(
            f"gh issue-meta lookup failed for #{number}: {getattr(e, 'stderr', None) or e}"
        ) from e
    return {"state": obj["state"], "is_pr": "pull_request" in obj}


def native_blockers(number):
    """Set of issue numbers this issue is already natively blockedBy.

    Reads the native blockedBy edge via a raw `gh api graphql` passthrough,
    NOT the `gh issue view --json blockedBy` client field. This is deliberate
    and load-bearing: the `--json blockedBy` field is validated against gh's
    built-in per-version field list and only exists on gh >= 2.94.0, whereas a
    GraphQL query string is passed to GitHub's server unchanged (the blockedBy
    edge is server-side GA since 2025-08-21), so this read is NOT client-
    version-gated. Do not "modernize" it to `--json blockedBy`: that would
    silently add a gh >= 2.94.0 requirement and break on older clients (e.g.
    the CCR sandbox's gh 2.65.0). GitHub enforces a max of 50 blockers per
    direction per issue, so first:50 is a true ceiling, not a sample.

    Returns an empty set if the issue simply has no blockers (a benign null
    issue node with no errors). Raises BlockerLookupError if the lookup could
    not be completed at all -- a failed `gh` call, or a GraphQL `errors`
    payload with no usable issue node -- so the caller skips that one issue
    instead of crashing the whole scan or misreading it as "no blockers"
    (Issue #989)."""
    q = (
        f'query {{ repository(owner: "{REPO_OWNER}", name: "{REPO_NAME}") {{'
        f"  issue(number: {number}) {{ blockedBy(first: 50) {{ nodes {{ number }} }} }}"
        "}}"
    )
    try:
        data = gh("api", "graphql", "-f", f"query={q}")
    except (GhError, json.JSONDecodeError) as e:
        # CalledProcessError = gh exited non-zero (5xx / network / GraphQL error);
        # JSONDecodeError = gh exited 0 with empty/non-JSON stdout. Both mean the
        # lookup could not be completed, so degrade rather than abort the scan.
        raise BlockerLookupError(
            f"gh graphql blockedBy lookup failed for #{number}: {getattr(e, 'stderr', None) or e}"
        ) from e
    # Use `or {}` at every level: a GraphQL top-level error nulls `data` (key
    # present, value None), so `.get("data", {})` returns None, not the default.
    issue = ((data.get("data") or {}).get("repository") or {}).get("issue")
    # gh can exit 0 with a partial response carrying a top-level `errors` array.
    # Trust it only if the issue node still came back; otherwise the lookup is
    # incomplete and must be treated as a failure, not an empty (no-blocker) result.
    if data.get("errors") and issue is None:
        raise BlockerLookupError(
            f"GraphQL errors on blockedBy lookup for #{number}: {data['errors']}"
        )
    nodes = ((issue or {}).get("blockedBy") or {}).get("nodes") or []
    return {n["number"] for n in nodes}


def reconcile(pairs):
    """Resolve + classify each (dependent, blocker) pair into a record dict.
    Caches per-issue lookups so each blocker/dependent is fetched once."""
    meta_cache = {}
    edges_cache = {}

    def cached(cache, n, fn):
        # Memoize fn(n), caching a raised IssueLookupError too, so a flaky
        # per-issue lookup is attempted once (not once per pair sharing it) and
        # re-raised on each subsequent hit.
        if n not in cache:
            try:
                cache[n] = fn(n)
            except IssueLookupError as e:
                cache[n] = e
        result = cache[n]
        if isinstance(result, IssueLookupError):
            raise result
        return result

    records = []
    for dependent, blocker in pairs:
        try:
            existing = cached(edges_cache, dependent, native_blockers)
            bmeta = cached(meta_cache, blocker, issue_meta)
        except IssueLookupError as e:
            # One flaky per-issue lookup (edge OR blocker-meta) must not abort the
            # whole scan. Skip this pair's classification and surface it in the
            # report so a failed lookup is never mistaken for a real result
            # (edges: Issue #989; blocker-meta: Issue #1012).
            print(f"scan_prose_deps: WARNING - skipping #{dependent} -> #{blocker}: {e}",
                  file=sys.stderr)
            records.append({
                "dependent": dependent,
                "blocker": blocker,
                "state": "?",
                "action": e.action,
            })
            continue
        action = classify(dependent, blocker, blocker_meta=bmeta, existing=existing)
        records.append({
            "dependent": dependent,
            "blocker": blocker,
            "state": bmeta["state"],
            "action": action,
        })
    return records


# *-lookup-failed first: they are scan-integrity warnings (a dependency this run
# could not determine), so they sort above the actionable-drift rows.
_ACTION_ORDER = [*_LOOKUP_FAILED_ACTIONS, "needs-wiring", "un-wireable-pr", "already-wired", "closed-blocker", "self-ref"]
_ACTION_RANK = {a: i for i, a in enumerate(_ACTION_ORDER)}


def render_report(records):
    if not records:
        return "scan_prose_deps: no prose-dependency drift found.\n"
    lines = [f"{'DEPENDENT':>9}  {'BLOCKER':>7}  {'STATE':<7}  ACTION",
             "-" * 44]
    key = lambda r: (_ACTION_RANK.get(r["action"], len(_ACTION_ORDER)), r["dependent"], r["blocker"])
    for r in sorted(records, key=key):
        lines.append(f"#{r['dependent']:<8} #{r['blocker']:<6} {r['state']:<7}  {r['action']}")
    return "\n".join(lines) + "\n"


def wire(records):
    """POST a native blockedBy edge for each record (must be needs-wiring).
    Uses the REST dependencies endpoint, which takes the blocker's numeric DB id."""
    wired = []
    for r in records:
        dependent, blocker = r["dependent"], r["blocker"]
        blocker_db_id = gh("api", f"repos/{REPO}/issues/{blocker}", parse_json=True)["id"]
        gh("api", "--method", "POST",
           f"repos/{REPO}/issues/{dependent}/dependencies/blocked_by",
           "-F", f"issue_id={blocker_db_id}", parse_json=False)
        print(f"  wired: #{dependent} blocked_by #{blocker}")
        wired.append((dependent, blocker))
    return wired


def _scan(issue_number=None):
    """Fetch -> parse -> reconcile. Single-issue if issue_number given.

    Single-issue scope fetches that issue directly (the hot path for the
    DoR / best-next per-candidate check) rather than listing the whole board."""
    if issue_number is not None:
        obj = gh("api", f"repos/{REPO}/issues/{issue_number}")
        issues = [{"number": obj["number"], "body": obj.get("body")}]
    else:
        issues = fetch_open_issues()
    pairs = []
    for i in issues:
        pairs.extend(parse_dependencies(i["number"], i.get("body")))
    return reconcile(pairs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--issue", type=int, help="restrict scan to a single dependent issue")
    parser.add_argument("--report", action="store_true",
                        help="print the drift table (default action; explicit form)")
    parser.add_argument("--check", action="store_true",
                        help="exit 1 if any needs-wiring drift, 2 if any lookup failed (incomplete scan)")
    parser.add_argument("--apply", action="store_true", help="wire the needs-wiring edges")
    parser.add_argument("--only", type=int, nargs="*", default=None,
                        help="with --apply: wire only these dependent issue numbers")
    args = parser.parse_args()
    if args.only is not None and not args.apply:
        parser.error("--only requires --apply")

    records = _scan(args.issue)
    needs = [r for r in records if r["action"] == "needs-wiring"]
    lookup_failed = [r for r in records if r["action"] in _LOOKUP_FAILED_ACTIONS]

    if args.apply:
        subset = needs if args.only is None else [r for r in needs if r["dependent"] in args.only]
        if not subset:
            print("scan_prose_deps: nothing to wire.")
            return 0
        wire(subset)
        return 0

    print(render_report(records), end="")
    if args.check:
        # Exit-code convention, shared with scan_unblocked.py (aligned in Issue
        # #1179): 1 = ran fine, gate tripped (drift found); 2 = could-not-run
        # (incomplete lookup). 2 takes precedence - a partially-failed scan can't
        # be trusted to have found all drift - and 2 == could-not-run matches
        # argparse's own exit-2-for-usage-error, so a routine shelling several PM
        # sweeps can read a 2 from any of them as "this did not execute".
        if lookup_failed:
            return 2
        return 1 if needs else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
