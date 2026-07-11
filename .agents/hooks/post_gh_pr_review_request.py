#!/usr/bin/env python3
"""PostToolUse hook: when a bot review is requested on a PR, advance BOTH the PR's
own board card and its linked Issue's card to "In review".

Sibling of `post_gh_pr_create.py`. That hook sets a fresh PR's linked-issue card
to "Ready for review"; the "Ready for review" -> "In review" hop (when a review
actually starts) was entirely manual and kept being forgotten, stranding cards
in "Ready for review" (e.g. #406, #234). Rung-3 mechanism escalation for a
recurring board-hygiene defect (Issue #996).

Fires on the two local commands that start a review:
  - `gh pr comment <ref> --body "@claude review"`  (the canonical bot trigger)
  - `gh pr review <ref> ...`                         (a direct review action)
It resolves the PR's linked Issue(s) via `closingIssuesReferences` and sets each
one's Status on project #9 to "In review", **and** sets the PR's own card too.

Why both (Issue #1108): `post_gh_pr_create` parks the PR card at "Ready for
review" and nothing ever advanced it, so the PR and its Issue split across two
columns and the PR's column depended on whether a human happened to drag it.
GitHub Docs endorse boarding PRs for review tracking, and no built-in automation
covers a review request (only closed -> Done and merged -> Done are native), so
this hook owns the transition for both cards. The PR flip is NOT gated on having
a linked Issue: a `--no-issue` companion PR still belongs in the review column.

Residual (documented, not solved here): a human reviewing directly on github.com
emits no local command, so a local PostToolUse hook cannot see it - that path
still needs a manual flip (or a future GitHub Action / webhook). The bot-review
path is our dominant one, so the hook covers the common case.

Reads PostToolUse hook JSON on stdin. Fails OPEN on any parse miss, untracked
repo, or `gh` error - a board-automation hook must never break the user's flow.
Idempotent (already-"In review" is a no-op) and it does not fight a parent parked
in "Epic" or a terminal "Done" card (both are skipped). On a real flip it emits
an `additionalContext` confirmation and logs one line to `.agents/hook_fires.jsonl`
(gitignored) per the fire-log infra (Issue #453).
"""
from __future__ import annotations

import json
import re
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _shell_parse  # noqa: E402

# Project #9 ("JH M Lee Lab") - user-level project, IDs stable + repo-independent
# (same values as post_gh_pr_create.py / recheck_dispatch.py).
PROJECT_ID = "PVT_kwHOB17eGc4BSomP"
PROJECT_NUMBER = 9
PROJECT_OWNER = "Jin-HoMLee"
STATUS_FIELD_ID = "PVTSSF_lAHOB17eGc4BSomPzhAHFf8"
IN_REVIEW_OPTION = "df73e18b"

# The canonical bot-review trigger the GitHub Action fires on. Matched
# case-insensitively as a literal substring (the hyphenated `@-claude review`
# reference form is NOT a trigger and deliberately does not match).
REVIEW_TRIGGER = "@claude review"

# Statuses we must not overwrite: already-there (idempotent no-op), a parent
# parked off-ladder ("Epic"), or a terminal card ("Done").
SKIP_STATUSES = frozenset({"In review", "Epic", "Done"})

# Only PRs in these repos feed project #9.
TRACKED_REPOS = {
    "splice-neoepitope-pipeline",
    "claude-personas-splice-neoepitope-pipeline",
}

LOG_PATH = Path(__file__).resolve().parent.parent.parent / ".agents" / "hook_fires.jsonl"
_PR_URL_RE = re.compile(r"https://github\.com/([\w.-]+)/([\w.-]+)/pull/(\d+)")
_PUNCT = set("();<>|&")  # shell punctuation_chars -> standalone separator tokens
_ASSIGNMENT_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*=")


# --- pure helpers (unit-tested) ---


def _tokenize(cmd: str) -> list[str] | None:
    """shlex-tokenize honoring quotes + shell punctuation, or None if unbalanced.

    Normalized first (Issue #1130) - heredoc bodies stripped, unquoted newlines
    turned into separators - so a review request written as
    `cat > c.md <<'EOF' ... EOF; gh pr comment N --body-file c.md` is seen. The
    sibling `post_gh_pr_create.py` had the identical gap, which left it silently
    dead on every heredoc-created PR; single-sourcing the normalizer is what keeps
    the two matchers from drifting apart on it again.
    """
    try:
        lex = shlex.shlex(_shell_parse.normalize_command(cmd),
                          posix=True, punctuation_chars=True)
        lex.whitespace_split = True
        return list(lex)
    except ValueError:
        return None


def _segment(tokens: list[str], start: int) -> list[str]:
    """Tokens from `start` up to (not including) the next pure-punctuation
    separator (`&&`, `||`, `;`, `|`, ...), i.e. the rest of the current command
    only. Bounds a body scan to one command so a trigger in a later `&&`-joined
    segment doesn't leak in.
    """
    seg = []
    for tok in tokens[start:]:
        if tok and all(ch in _PUNCT for ch in tok):
            break
        seg.append(tok)
    return seg


def matches_review_request(cmd: str) -> str | None:
    """Return the PR ref (number or URL, as a string) of a review-request command.

    Returns None when `cmd` is not a review request. Tokenizes with shlex (like
    the sibling's `matches_pr_create`) so a `gh pr comment`/`gh pr review`
    appearing inside a quoted argument - e.g. a body that merely discusses the
    command - does not match. Compound commands are handled: a review-request
    segment after a real shell separator (`&&`, `||`, `;`, `|`) matches.

    Two shapes match:
      - `gh pr comment <ref> ...` whose body contains the `@claude review` trigger
      - `gh pr review <ref> ...`  (any direct review action)

    The PR ref is taken as the first bare positional immediately after the
    subcommand - the strong gh convention (`gh pr comment 996 --body ...`). If the
    token right after the subcommand is a flag, that segment can't be parsed
    reliably and is skipped (we don't mis-grab a flag value). Compound commands
    are scanned segment by segment: a non-matching `gh pr comment` segment does
    NOT abort the scan, so a later `&&`-joined `gh pr review` still matches, and a
    comment's trigger scan is bounded to its own segment so a trigger in a
    separate `echo` can't false-positive. Untokenizable input (unbalanced quotes)
    fails safe (None).
    """
    tokens = _tokenize(cmd)
    if tokens is None:
        return None

    at_command_start = True
    for i, tok in enumerate(tokens):
        if tok and all(ch in _PUNCT for ch in tok):  # pure-punctuation = separator
            at_command_start = True
            continue
        if at_command_start and _ASSIGNMENT_RE.match(tok):
            continue  # leading `VAR=value` prefix -> command starts after it
        if at_command_start and tuple(tokens[i:i + 3]) in (
            ("gh", "pr", "comment"), ("gh", "pr", "review"),
        ):
            subcommand = tokens[i + 2]
            ref = tokens[i + 3] if i + 3 < len(tokens) else None
            if ref is not None and not ref.startswith("-"):
                if subcommand == "review":
                    return ref
                # comment: only a body (within THIS segment) carrying the trigger
                if any(_has_trigger(t) for t in _segment(tokens, i + 3)):
                    return ref
            # no match in this segment -> keep scanning later segments
            at_command_start = False
            continue
        at_command_start = False
    return None


def _has_trigger(text: str) -> bool:
    """True if `text` contains the review trigger (case-insensitive literal)."""
    return REVIEW_TRIGGER in (text or "").lower()


def parse_pr_url(text: str) -> tuple[str, str, int] | None:
    """Return (owner, repo, number) of the LAST PR URL in `text`, or None."""
    matches = _PR_URL_RE.findall(text or "")
    if not matches:
        return None
    owner, repo, number = matches[-1]
    return owner, repo, int(number)


def should_track(owner: str, repo: str) -> bool:
    """True only for PRs that belong on project board #9."""
    return owner == PROJECT_OWNER and repo in TRACKED_REPOS


def should_flip(current_status: str | None) -> bool:
    """True if a card in `current_status` should be advanced to In review.

    Skips already-In-review (idempotent), a parent parked in Epic, and a terminal
    Done card. Every other state (Ready for review, Ready, In progress for a draft
    under review, or an unset status) is a legitimate forward move.
    """
    return current_status not in SKIP_STATUSES


# --- gh I/O (fail-open) ---


def _gh(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["gh", *args], capture_output=True, text=True, timeout=15, check=True
    )


def _pr_linked_issues(ref: str) -> tuple[str, str, int, list[int]] | None:
    """Resolve a PR ref to (owner, repo, pr_number, [linked issue numbers]).

    `gh pr view <ref>` accepts a number (resolved against the cwd repo), a URL, or
    a branch, so this works from any tracked clone. Returns None if the URL can't
    be parsed (repo unresolvable). The PR's own number is returned so its board
    card can be advanced alongside its linked Issues (Issue #1108) -- the linked
    list may legitimately be empty (a `--no-issue` companion PR).
    """
    res = _gh("pr", "view", ref, "--json", "url,closingIssuesReferences")
    data = json.loads(res.stdout)
    parsed = parse_pr_url(data.get("url") or "")
    if parsed is None:
        return None
    owner, repo, pr_number = parsed
    issues = [
        n for n in (
            (r or {}).get("number") for r in data.get("closingIssuesReferences") or []
        ) if isinstance(n, int)
    ]
    return owner, repo, pr_number, issues


def _item_and_status(
    kind: str, number: int, owner: str, repo: str
) -> tuple[str | None, str | None]:
    """Return (project-item id on #9, current Status name) for an Issue or PR.

    `kind` is the GraphQL field name: "issue" or "pullRequest". Both content types
    expose `projectItems`, so one resolver serves the linked Issue *and* the PR's
    own card -- keeping a single source of truth for the field-name matching.

    (None, None) if the item is not on the board. Mirrors recheck_dispatch's
    projectItems resolution. The repo is threaded from the PR (not hardcoded) so a
    tracked personas-repo PR resolves its OWN board items rather than a
    same-numbered pipeline issue (`closingIssuesReferences` is same-repo-only, so
    a personas PR's linked numbers are personas numbers).
    """
    query = (
        'query { repository(owner: "' + owner + '", name: "' + repo + '") '
        '{ ' + kind + '(number: ' + str(number) + ') { projectItems(first: 10) { nodes { '
        'id project { number } fieldValues(first: 20) { nodes { '
        '... on ProjectV2ItemFieldSingleSelectValue { name field { '
        '... on ProjectV2FieldCommon { name } } } } } } } } } }'
    )
    res = _gh("api", "graphql", "-f", f"query={query}")
    data = json.loads(res.stdout)
    node = ((data.get("data") or {}).get("repository") or {}).get(kind) or {}
    for item in (node.get("projectItems") or {}).get("nodes") or []:
        if ((item.get("project") or {}).get("number")) != PROJECT_NUMBER:
            continue
        item_id = item.get("id")
        for fv in (item.get("fieldValues") or {}).get("nodes") or []:
            if fv and (fv.get("field") or {}).get("name") == "Status":
                return item_id, fv.get("name")
        return item_id, None
    return None, None


def _set_status(item_id: str, option_id: str) -> None:
    _gh(
        "api", "graphql", "-f",
        "query=mutation($p:ID!,$i:ID!,$f:ID!,$o:String!){"
        "updateProjectV2ItemFieldValue(input:{projectId:$p,itemId:$i,fieldId:$f,"
        "value:{singleSelectOptionId:$o}}){projectV2Item{id}}}",
        "-f", f"p={PROJECT_ID}", "-f", f"i={item_id}",
        "-f", f"f={STATUS_FIELD_ID}", "-f", f"o={option_id}",
    )


def _log_fire(pr: int | None, issues: list[int], repo: str) -> None:
    """Append one fire-log line (Issue #453 infra). Never raises.

    `pr` is the PR number when its own card was flipped, else None (already
    In review / Done / not on the board).
    """
    payload = {
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "hook": "post_gh_pr_review_request",
        "pr": pr,
        "issues": issues,
        "repo": repo,
        "action": "status:In review",
    }
    line = json.dumps(payload, separators=(",", ":")) + "\n"
    try:
        LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
        with open(LOG_PATH, "a", encoding="utf-8") as f:
            f.write(line)
    except OSError:
        pass


# --- orchestration ---


def apply_review_request(ref: str) -> dict | None:
    """Flip the PR's own card + every linked Issue to `In review`. Fails open.

    Extracted from `main()` so it has two callers and stays single-sourced:

    1. `main()` - a review requested from the shell, which this hook sees as a
       Bash tool-call.
    2. `post_gh_pr_create.py` - which auto-requests the review at PR-open time
       (Issue #1073). That request is a *subprocess* `gh pr comment`, not a Claude
       tool-call, so this hook's PostToolUse matcher never sees it. Without a
       shared entry point the auto-requested PR would sit at `Ready for review`
       for its whole review - re-creating the exact stranding Issue #996 fixed.

    Returns a summary dict of what flipped, or None when nothing did.
    """
    # Accumulate what actually flipped OUTSIDE the try, so a mid-sequence gh error
    # still logs + surfaces the board mutations that already landed. The fire-log is
    # the audit trail the review-debt/health tooling reads, so a real flip that goes
    # unrecorded is a silent board/log divergence (Issue #1108 review, finding 2).
    owner = repo = None
    pr_number = None
    flipped_pr = None
    flipped: list[int] = []
    try:
        resolved = _pr_linked_issues(ref)
        if resolved is None:
            return None
        owner, repo, pr_number, issues = resolved
        if not should_track(owner, repo):
            return None

        # The PR's OWN card first (Issue #1108). Deliberately NOT gated on `issues`:
        # a `--no-issue` companion PR has no linked Issue but still belongs in the
        # review column. Previously the empty-`issues` early-return skipped it.
        pr_item, pr_status = _item_and_status("pullRequest", pr_number, owner, repo)
        if pr_item is not None and should_flip(pr_status):
            _set_status(pr_item, IN_REVIEW_OPTION)
            flipped_pr = pr_number

        for issue in issues:
            item_id, status = _item_and_status("issue", issue, owner, repo)
            if item_id is None or not should_flip(status):
                continue
            _set_status(item_id, IN_REVIEW_OPTION)
            flipped.append(issue)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired,
            json.JSONDecodeError, FileNotFoundError):
        pass  # fail open - fall through to log/surface whatever already flipped

    # `owner is None` means we failed before resolving the repo, so nothing flipped.
    if owner is None or (flipped_pr is None and not flipped):
        return None

    _log_fire(flipped_pr, flipped, f"{owner}/{repo}")
    return {
        "owner": owner,
        "repo": repo,
        "pr_number": pr_number,
        "flipped_pr": flipped_pr,
        "flipped": flipped,
    }


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0  # fail open

    cmd = (payload.get("tool_input") or {}).get("command", "")
    ref = matches_review_request(cmd)
    if ref is None:
        return 0

    result = apply_review_request(ref)
    if result is None:
        return 0

    parts = []
    if result["flipped_pr"] is not None:
        parts.append(f"PR #{result['flipped_pr']}")
    if result["flipped"]:
        parts.append(
            "linked Issue(s) " + ", ".join(f"#{n}" for n in result["flipped"])
        )
    print(json.dumps({
        "hookSpecificOutput": {
            "hookEventName": "PostToolUse",
            "additionalContext": (
                f"post_gh_pr_review_request: review requested on PR "
                f"#{result['pr_number']} ({result['owner']}/{result['repo']}) - "
                f"{' + '.join(parts)} Status -> In review."
            ),
        }
    }))
    return 0


if __name__ == "__main__":
    sys.exit(main())
