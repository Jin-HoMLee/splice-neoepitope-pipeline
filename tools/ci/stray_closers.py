#!/usr/bin/env python3
"""Detect stray GitHub closing-keyword references in a PR's would-be squash body.

Background — the `closingIssuesReferences` blind spot
-----------------------------------------------------
A closing keyword (`close[sd]?`, `fix(es|ed)?`, `resolve[sd]?`) immediately
followed by `#N` anywhere in a squash-merge commit body auto-closes Issue N on
the default branch — even when the text merely *describes* the behavior in prose.
A PR's `closingIssuesReferences` API surfaces only PR-**body** link edges, NOT
keywords inside commit-message bodies (which the squash commit inherits), so a
pre-merge `closingIssuesReferences` check passes clean and misses them.

That bit PR #543 → parent epic Issue #538 (2026-05-28): the lab-notebook PR's
commit body said "would auto-close #538 on merge" (the hyphen in "auto-close" is
a word boundary), silently closing the epic before its sub-issues completed.

This module assembles the text that could land in the squash commit (PR title +
body + every commit's message) and flags any closing-keyword → `#N` where N is
NOT in the PR's intended closing set (`closingIssuesReferences`). Used as a
pre-merge gate in `scripts/audit_and_merge.sh`, sibling of the closure-ritual
check there.

Escape workaround (for genuine non-closing references): break the
keyword→`#N` adjacency — neutral phrasing like "related to Issue #N" or put the
keyword far from the token. See shared/feedback_hash_numbers.md (Companion
foot-gun section).

See https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/559.
"""
import json
import os
import re
import subprocess
import sys

# GitHub closing keywords, case-insensitive, anchored on word boundaries so
# "disclose"/"prefix"/"closing"/"fixing"/"closer" do NOT match but "auto-close"
# (hyphen = boundary) does. Separator is HORIZONTAL whitespace + optional colon
# ([ \t], not \s) — matching GitHub's `Closes: #N` / `Closes #N` forms while NOT
# bridging a newline: assemble_squash_text joins fields with "\n", and \s would
# let a field ending in a bare keyword cross-match a #N opening the next field
# (a false positive GitHub itself would not act on). Per PR #562 review.
_CLOSER_RE = re.compile(
    r"\b(close[sd]?|fix(?:es|ed)?|resolve[sd]?)\b[ \t]*:?[ \t]*#(\d+)",
    re.IGNORECASE,
)


# --- pure detection (unit-tested, no I/O) ---


def find_stray_closers(text, closing_set):
    """Return [(keyword, number, line), ...] for closing-keyword → #N refs in
    `text` whose N is NOT in `closing_set` (the PR's intended closing issues).

    Deduplicated by issue number (first offending line kept). `closing_set` may
    hold ints or strings. Returns [] for empty/None text.
    """
    intended = {int(n) for n in closing_set}
    seen = set()
    out = []
    for m in _CLOSER_RE.finditer(text or ""):
        number = int(m.group(2))
        if number in intended or number in seen:
            continue
        seen.add(number)
        line = _line_containing(text, m.start())
        out.append((m.group(1), number, line))
    return out


def _line_containing(text, pos):
    """The full line of `text` containing character offset `pos`, stripped."""
    start = text.rfind("\n", 0, pos) + 1
    end = text.find("\n", pos)
    if end == -1:
        end = len(text)
    return text[start:end].strip()


def assemble_squash_text(pr_data):
    """Join the fields that can land in a squash commit body: PR title + body +
    every commit's headline and body. Tolerant of missing/None fields.
    """
    parts = [pr_data.get("title") or "", pr_data.get("body") or ""]
    for c in pr_data.get("commits") or []:
        parts.append(c.get("messageHeadline") or "")
        parts.append(c.get("messageBody") or "")
    return "\n".join(parts)


def closing_set_from(pr_data):
    """The set of issue numbers the PR legitimately closes (closingIssuesReferences)."""
    return {ref["number"] for ref in (pr_data.get("closingIssuesReferences") or [])}


# --- gh I/O + CLI ---


def fetch_pr(pr, repo):
    res = subprocess.run(
        ["gh", "pr", "view", str(pr), "--repo", repo, "--json",
         "title,body,commits,closingIssuesReferences"],
        capture_output=True, text=True, check=True, timeout=30,
    )
    return json.loads(res.stdout)


def main(argv):
    if len(argv) < 2 or not argv[1].isdigit():
        print("usage: stray_closers.py <PR_NUMBER>", file=sys.stderr)
        return 2
    pr = argv[1]
    repo = os.environ.get("REPO", "Jin-HoMLee/splice-neoepitope-pipeline")
    try:
        data = fetch_pr(pr, repo)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired,
            FileNotFoundError, json.JSONDecodeError) as e:
        # Fail OPEN: a gh hiccup must not block a legitimate merge.
        print(f"⚠ stray-closer check skipped (gh error: {e}).", file=sys.stderr)
        return 0

    stray = find_stray_closers(assemble_squash_text(data), closing_set_from(data))
    if not stray:
        return 0

    print(f"✗ PR #{pr} would-be squash body contains stray closing keyword(s) "
          f"targeting Issue(s) NOT in its closing set:", file=sys.stderr)
    for kw, number, line in stray:
        print(f"    {kw} #{number}  ←  {line}", file=sys.stderr)
    print("    These auto-close the Issue(s) on merge (GitHub ignores negation). "
          "Break the keyword→#N adjacency (e.g. \"related to Issue #N\") in the PR "
          "body AND commit messages, or add the Issue to closingIssuesReferences if "
          "the close is intended. See shared/feedback_hash_numbers.md.", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
