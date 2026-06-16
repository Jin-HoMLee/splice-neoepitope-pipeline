#!/usr/bin/env python3
"""Reconcile prose-stated Issue dependencies into native GitHub blockedBy edges.

Finds bodies saying "depends on #M" / "blocked on #M" etc., diffs them against
the native blockedBy graph, and reports or wires the missing edges.

Usage:
  scripts/pm/scan_prose_deps.py                 # --report (default): drift table
  scripts/pm/scan_prose_deps.py --issue N       # single-issue scope
  scripts/pm/scan_prose_deps.py --check         # exit 2 if any drift
  scripts/pm/scan_prose_deps.py --apply         # wire all needs-wiring records
  scripts/pm/scan_prose_deps.py --apply --only 745 594   # wire only these dependents

Exits 0 clean / 2 drift-present (--check) / 1 on error.
"""
import argparse
import json
import re
import subprocess
import sys

REPO = "Jin-HoMLee/splice-neoepitope-pipeline"


def gh(*args, parse_json=True):
    """Run a gh command; return parsed JSON (default) or raw stdout text."""
    result = subprocess.run(["gh", *args], capture_output=True, text=True, check=True)
    return json.loads(result.stdout) if parse_json else result.stdout


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


def parse_dependencies(number, body):
    """Return sorted, deduped (dependent, blocker) pairs from one issue body.

    Pure: no I/O. Matches only the blocker-phrase allowlist; drops self-references.
    Blocker open/closed and PR-vs-issue filtering happen later, in reconcile().
    """
    if not body:
        return []
    blockers = {int(m.group(1)) for m in BLOCKER_RE.finditer(body)}
    blockers.discard(number)  # self-reference
    return [(number, b) for b in sorted(blockers)]


def classify(dependent, blocker, *, blocker_meta, existing):
    """Classify one (dependent, blocker) pair. Pure given resolved inputs."""
    if dependent == blocker:
        return "self-ref"
    if blocker in existing:
        return "already-wired"
    if blocker_meta["state"] == "closed":
        return "closed-blocker"  # convention: only wire currently-open blockers
    if blocker_meta["is_pr"]:
        return "un-wireable-pr"  # native deps are issue<->issue
    return "needs-wiring"


def issue_meta(number):
    """{'state': 'open'|'closed', 'is_pr': bool} via the REST issues endpoint
    (which serves both issues and PRs; a PR carries a 'pull_request' key)."""
    obj = gh("api", f"repos/{REPO}/issues/{number}")
    return {"state": obj["state"], "is_pr": "pull_request" in obj}


def native_blockers(number):
    """Set of issue numbers this issue is already natively blockedBy."""
    q = (
        'query { repository(owner: "Jin-HoMLee", name: "splice-neoepitope-pipeline") {'
        f"  issue(number: {number}) {{ blockedBy(first: 50) {{ nodes {{ number }} }} }}"
        "}}"
    )
    data = gh("api", "graphql", "-f", f"query={q}")
    nodes = data["data"]["repository"]["issue"]["blockedBy"]["nodes"]
    return {n["number"] for n in nodes}


def reconcile(pairs):
    """Resolve + classify each (dependent, blocker) pair into a record dict.
    Caches per-issue lookups so each blocker/dependent is fetched once."""
    meta_cache = {}
    edges_cache = {}

    def meta(n):
        if n not in meta_cache:
            meta_cache[n] = issue_meta(n)
        return meta_cache[n]

    def edges(n):
        if n not in edges_cache:
            edges_cache[n] = native_blockers(n)
        return edges_cache[n]

    records = []
    for dependent, blocker in pairs:
        bmeta = meta(blocker)
        action = classify(dependent, blocker, blocker_meta=bmeta, existing=edges(dependent))
        records.append({
            "dependent": dependent,
            "blocker": blocker,
            "state": bmeta["state"],
            "action": action,
        })
    return records
