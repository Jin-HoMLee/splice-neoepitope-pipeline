#!/usr/bin/env python3
"""Detect whether a bot-review trigger has been offered on a PR (Issue #443).

Background — the bot-review-offer slip
--------------------------------------
The project convention is to proactively offer a `@claude review` after any
non-trivial PR (`shared/feedback_github_workflow.md`); the user prefers being
asked rather than relying on the rule being remembered. Two same-morning slips
(PR #441 + PR #442 merged without the offer, 2026-05-21) crossed the
memory→mechanism escalation threshold, so the offer is enforced at the moment of
action — the merge invocation in `scripts/audit_and_merge.sh`.

This module answers one question deterministically: has the *real* review trigger
(`@claude review`) already been posted in the PR's comments? `audit_and_merge.sh`
shells out here and, if not, prompts the user (offer / skip-trivial / cancel).

What counts as "offered"
------------------------
Only the literal trigger `@claude review` (case-insensitive, contiguous on one
line) — the string that actually triggers the Claude Code GitHub Action. The
hyphenated `@-claude review` reference form does NOT count: it never triggers a
review (the literal substring `@claude` is absent), so treating it as "offered"
would let a real PR merge un-reviewed.

Fails OPEN: a `gh` hiccup prints `OFFERED` (skip the prompt) rather than blocking
the merge with a spurious prompt — mirroring the sibling gates in
`scripts/audit_and_merge.sh`.

See https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/443.
"""
import json
import os
import re
import subprocess
import sys

# The real review trigger: `@claude` + horizontal whitespace + `review`,
# case-insensitive. `[ \t]+` (not `\s`) keeps it contiguous on one line, so a
# bare `@claude` ending a line cannot cross-match a `review` opening the next.
# `@-claude review` does not match — the hyphen breaks the `@claude` literal.
_TRIGGER_RE = re.compile(r"@claude[ \t]+review", re.IGNORECASE)


# --- pure detection (unit-tested, no I/O) ---


def has_bot_review_offer(comments):
    """Return True if any comment body contains the `@claude review` trigger.

    `comments` is the `gh pr view --json comments` shape — a list of dicts with a
    "body" key. Tolerant of None list, None/missing bodies.
    """
    for c in comments or []:
        body = (c or {}).get("body") or ""
        if _TRIGGER_RE.search(body):
            return True
    return False


# --- gh I/O + CLI ---


def fetch_comments(pr, repo):
    res = subprocess.run(
        ["gh", "pr", "view", str(pr), "--repo", repo, "--json", "comments"],
        capture_output=True, text=True, check=True, timeout=30,
    )
    return json.loads(res.stdout).get("comments", [])


def main(argv):
    if len(argv) < 2 or not argv[1].isdigit():
        print("usage: bot_review_offer.py <PR_NUMBER>", file=sys.stderr)
        return 2
    pr = argv[1]
    repo = os.environ.get("REPO", "Jin-HoMLee/splice-neoepitope-pipeline")
    try:
        comments = fetch_comments(pr, repo)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired,
            FileNotFoundError, json.JSONDecodeError) as e:
        # Fail OPEN: a gh hiccup must not disrupt the merge with a spurious prompt.
        print(f"⚠ bot-review-offer check skipped (gh error: {e}).", file=sys.stderr)
        print("OFFERED")
        return 0

    print("OFFERED" if has_bot_review_offer(comments) else "NOT_OFFERED")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
