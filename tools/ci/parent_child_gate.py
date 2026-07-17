#!/usr/bin/env python3
"""Blocking pre-merge parent/child gate for the closure ritual (Issue #1155).

Moves the parent guard from branch-time to merge-time. The PR #543 -> Issue #538
harm was a PR *auto-closing a parent while its children were still open*, which
orphans them. That harm lands **at merge**, so we now guard there instead of
banning the safe, reversible act of branching off a parent.

The gate refuses a merge whose PR's `closingIssuesReferences` contains a **parent**
Issue (`subIssuesSummary.total > 0`) with at least one **open** child, and names
the open children. A parent whose children are all closed passes.

Fail posture is **B (fail-closed on the parent check)**, a deliberate, documented
divergence from the sibling gates (stray_closers / lab_notebook_gate /
cross_repo_ac_gate), which fail OPEN as best-effort convenience controls. Here the
guarded harm (orphaning children) is irreversible, and our `feedback_fail_safe_not
_fail_open` memory says a gate must never false-PASS - so when the parent/child
state cannot be determined (any gh/network/JSON error), the gate **blocks** with a
"could not verify" message rather than admitting the merge. A transient false-block
costs a re-run; a false-pass orphans an epic's children. Ratified by Jin-Ho
2026-07-17 (Issue #1155 Notes).

Fail-OPEN only where there is genuinely nothing to check: a PR with no closing
references, or a linked Issue that is a leaf (not a parent). A gh error while
*enumerating* the closing set also fails closed - not being able to look is not
the same as having looked and found nothing, and during a GitHub outage the final
`gh pr merge` cannot succeed anyway, so failing closed there costs nothing extra.

Exit codes:
    0 - no linked closing Issue is a parent with an open child (or nothing to check)
    1 - a parent with an open child, OR an undetermined parent state (fail-closed)
    2 - usage error

Usage:
    python parent_child_gate.py <PR_NUMBER>
"""

import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_DEFAULT = "Jin-HoMLee/splice-neoepitope-pipeline"
LOG_PATH = Path(__file__).resolve().parent.parent.parent / ".agents" / "hook_fires.jsonl"

_IO_ERRORS = (
    subprocess.CalledProcessError,
    subprocess.TimeoutExpired,
    FileNotFoundError,
    json.JSONDecodeError,
    KeyError,
    TypeError,
    ValueError,
)


# --- pure aggregation (unit-tested, no I/O) ---


def gate_decision(resolutions):
    """Given per-Issue resolutions, return (blocked: bool, messages: list[str]).

    Each resolution is a dict with `number` and `status` in
    {`leaf`, `clean`, `open_children`, `undetermined`}; an `open_children`
    resolution also carries `open` (a list of open child numbers). Blocks on any
    `open_children` (a real orphaning risk) or `undetermined` (fail-closed).
    `leaf` and `clean` never block.
    """
    blocked = False
    messages = []
    for r in resolutions:
        status = r["status"]
        num = r["number"]
        if status == "open_children":
            blocked = True
            kids = ", ".join(f"#{c}" for c in r["open"])
            messages.append(
                f"✗ PR closes parent Issue #{num}, which still has open "
                f"child(ren): {kids}. Merging would auto-close the parent and "
                "orphan them. Close the children first, or unlink the parent "
                "from the PR's closing references."
            )
        elif status == "undetermined":
            blocked = True
            messages.append(
                f"✗ Could not verify the parent/child state of linked Issue "
                f"#{num} (gh/network error). Failing closed: the orphaning harm "
                "is irreversible, so re-run when GitHub is reachable rather than "
                "admit an unverified merge."
            )
    return blocked, messages


# --- gh I/O (fail-CLOSED: any error becomes an 'undetermined' block) ---


def _gh(*args, timeout=15):
    return subprocess.run(
        ["gh", *args], capture_output=True, text=True, timeout=timeout, check=True
    )


def _fetch_parent_state(number, repo):
    """Return (sub_total, children) for one Issue via a single GraphQL call.

    `children` is a list of {"number": int, "state": "OPEN"|"CLOSED"}. Raises on
    any gh/JSON error (the caller maps that to an 'undetermined' fail-closed
    resolution).
    """
    owner, name = repo.split("/", 1)
    query = (
        "query($o:String!,$n:String!,$num:Int!){"
        "repository(owner:$o,name:$n){issue(number:$num){"
        "subIssuesSummary{total}"
        "subIssues(first:100){nodes{number state}}}}}"
    )
    res = _gh("api", "graphql", "-f", f"query={query}",
              "-f", f"o={owner}", "-f", f"n={name}", "-F", f"num={number}")
    data = json.loads(res.stdout)
    issue = data["data"]["repository"]["issue"]
    total = int(issue["subIssuesSummary"]["total"])
    nodes = issue["subIssues"]["nodes"]
    children = [{"number": int(c["number"]), "state": c["state"]} for c in nodes]
    return total, children


def resolve(number, repo):
    """Classify one linked closing Issue into a resolution dict.

    Fail-CLOSED: any gh/JSON error while fetching parent/child state yields
    status 'undetermined' (which `gate_decision` treats as a block).
    """
    try:
        total, children = _fetch_parent_state(number, repo)
    except _IO_ERRORS:
        return {"number": number, "status": "undetermined"}
    if total == 0:
        return {"number": number, "status": "leaf"}
    open_kids = [c["number"] for c in children if c["state"] == "OPEN"]
    if open_kids:
        return {"number": number, "status": "open_children", "open": open_kids}
    return {"number": number, "status": "clean"}


def fetch_closing_issues(pr, repo):
    """Issue numbers in the PR's native closingIssuesReferences. Raises on error."""
    res = _gh("pr", "view", str(pr), "--repo", repo,
              "--json", "closingIssuesReferences")
    data = json.loads(res.stdout)
    return [ref["number"] for ref in (data.get("closingIssuesReferences") or [])]


def _log_fire(numbers, repo):
    """Append one fire-log line on a block (Issue #453 infra). Never raises."""
    for number in numbers:
        payload = {
            "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "hook": "parent_child_gate",
            "issue": number,
            "repo": repo,
            "action": "blocked-parent-with-open-child-or-undetermined",
        }
        line = json.dumps(payload, separators=(",", ":")) + "\n"
        try:
            LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
            with open(LOG_PATH, "a", encoding="utf-8") as f:
                f.write(line)
        except OSError:
            pass


# --- CLI ---


def main():
    if len(sys.argv) != 2 or not sys.argv[1].isdigit():
        print("usage: parent_child_gate.py <PR_NUMBER>", file=sys.stderr)
        return 2
    pr = sys.argv[1]
    repo = os.environ.get("REPO", REPO_DEFAULT)

    try:
        linked = fetch_closing_issues(pr, repo)
    except _IO_ERRORS as e:
        # Cannot enumerate the closing set -> cannot prove the merge is safe.
        # Fail CLOSED (posture B): not-looked is not the same as nothing-to-check.
        print(f"✗ Could not fetch PR #{pr} closing issues (gh error: {e}). "
              "Failing closed - re-run when GitHub is reachable.", file=sys.stderr)
        return 1

    if not linked:
        return 0  # nothing to check

    resolutions = [resolve(n, repo) for n in linked]
    blocked, messages = gate_decision(resolutions)
    for m in messages:
        print(m, file=sys.stderr)
    if blocked:
        offenders = [r["number"] for r in resolutions
                     if r["status"] in ("open_children", "undetermined")]
        _log_fire(offenders, repo)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
