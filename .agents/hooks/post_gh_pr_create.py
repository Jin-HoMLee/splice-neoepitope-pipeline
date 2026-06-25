#!/usr/bin/env python3
"""PostToolUse hook: auto-add a freshly-created PR to project board #9 and set
its Status to match the PR's review-readiness.

Automates the deterministic steps of the PR-open 4-step checklist
(`feedback_pr_open_checklist.md`); steps 3 (`**Created by:**` attribution) and 4
(Issue mirror on review request) stay manual — they need author judgment a hook
can't infer. Sibling of `scripts/audit_and_merge.sh` (closure-ritual gate on
`gh pr merge`).

Status is **draft-aware**: a non-draft PR-open is the review request, so Status
→ "Ready for review"; a draft PR is opened mid-In-progress (CI / shareable URL,
not up for review yet), so Status → "In progress". Either way the PR lands ON
the board — the original failure (PRs #6 + #543 stuck at Backlog 2026-05-28
because the rule lived in memory, not a mechanism) was about board-absence.
Draft-ness is read from the PR's authoritative `isDraft` field, not the command
string, so `--draft`, `-d`, and repo defaults are all covered.

Reads PostToolUse hook JSON on stdin, parses the PR URL from the tool output,
then runs the `gh` mutations. Fails OPEN on any parse miss, untracked repo, or
`gh` error — a board-automation hook must never break the user's flow. On a
successful set it emits an `additionalContext` confirmation and logs one line to
`.agents/hook_fires.jsonl` (gitignored) per the fire-log infra (Issue #453).

See Issue #550 for the originating rule.
"""
from __future__ import annotations

import json
import re
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

# Project #9 ("JH M Lee Lab") — a user-level project, so these IDs are stable
# and repo-independent (same values used by recheck_dispatch.py).
PROJECT_ID = "PVT_kwHOB17eGc4BSomP"
PROJECT_NUMBER = 9
PROJECT_OWNER = "Jin-HoMLee"
STATUS_FIELD_ID = "PVTSSF_lAHOB17eGc4BSomPzhAHFf8"
IN_PROGRESS_OPTION = "47fc9ee4"
READY_FOR_REVIEW_OPTION = "8bf9192f"

# Only PRs in these repos feed project #9. A guard against accidentally boarding
# an unrelated PR the session happens to open elsewhere.
TRACKED_REPOS = {
    "splice-neoepitope-pipeline",
    "claude-personas-splice-neoepitope-pipeline",
}

LOG_PATH = Path(__file__).resolve().parent.parent.parent / ".agents" / "hook_fires.jsonl"
_PR_URL_RE = re.compile(r"https://github\.com/([\w.-]+)/([\w.-]+)/pull/(\d+)")
_PR_CREATE_PREFIX = ("gh", "pr", "create")
_PUNCT = set("();<>|&")  # shell punctuation_chars → standalone separator tokens
# A leading `VAR=value` assignment at a command-start position (shlex strips the
# quotes in posix mode, so `B="x"` arrives as the token `B=x`). Mirrors the
# harness's subcommand-aware Bash matcher, which strips such prefixes.
_ASSIGNMENT_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*=")


# --- pure helpers (unit-tested) ---


def matches_pr_create(cmd: str) -> bool:
    """True only if `cmd` actually *invokes* `gh pr create` as a command.

    Tokenizes with shlex (quotes + shell punctuation respected) rather than
    substring-matching the raw command line, so `gh pr create` appearing INSIDE
    a quoted argument — e.g. a `gh pr comment` body that discusses the command —
    does NOT match. That false positive tripped on this hook's own PR #558
    review-reply (`... git push && gh pr create ...` inside the comment body).

    Compound commands are handled: a `gh pr create` segment after a real shell
    separator (`&&`, `||`, `;`, `|`) matches; one buried in quotes does not.
    Leading `VAR=value` environment-assignment prefixes are skipped so the
    command start is found after them (the harness Bash matcher does the same) —
    `B="x" gh pr create …` and its `VAR=…`-newline-`gh pr create` form both match.
    Untokenizable input (unbalanced quotes) fails safe → no match.
    """
    try:
        lex = shlex.shlex(cmd or "", posix=True, punctuation_chars=True)
        lex.whitespace_split = True
        tokens = list(lex)
    except ValueError:
        return False
    at_command_start = True
    for i, tok in enumerate(tokens):
        if tok and all(ch in _PUNCT for ch in tok):  # pure-punctuation = separator
            at_command_start = True
            continue
        if at_command_start and _ASSIGNMENT_RE.match(tok):
            continue  # leading `VAR=value` prefix → command starts after it
        if at_command_start and tuple(tokens[i:i + 3]) == _PR_CREATE_PREFIX:
            return True
        at_command_start = False
    return False


def parse_pr_url(text: str) -> tuple[str, str, int] | None:
    """Return (owner, repo, number) of the LAST PR URL in `text`, or None.

    `gh pr create` prints the new PR's URL to stdout; taking the last match is
    robust to leading log lines that may reference other URLs.
    """
    matches = _PR_URL_RE.findall(text or "")
    if not matches:
        return None
    owner, repo, number = matches[-1]
    return owner, repo, int(number)


def should_track(owner: str, repo: str) -> bool:
    """True only for PRs that belong on project board #9."""
    return owner == PROJECT_OWNER and repo in TRACKED_REPOS


def status_for_draft(is_draft: bool) -> tuple[str, str]:
    """(human label, single-select option id) for the PR's board Status.

    A draft PR is mid-In-progress; a non-draft PR-open is the review request.
    """
    if is_draft:
        return "In progress", IN_PROGRESS_OPTION
    return "Ready for review", READY_FOR_REVIEW_OPTION


def extract_output(tool_response) -> str:
    """Coax a searchable string out of the PostToolUse `tool_response`.

    Bash tool responses are usually a dict with a `stdout` key, but the shape is
    not contractually fixed — fall back to stringifying the whole object so a
    URL anywhere in it is still found.
    """
    if tool_response is None:
        return ""
    if isinstance(tool_response, str):
        return tool_response
    if isinstance(tool_response, dict):
        stdout = tool_response.get("stdout")
        if isinstance(stdout, str) and stdout:
            return stdout
    return json.dumps(tool_response)


# --- gh I/O (fail-open) ---


def _gh(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["gh", *args], capture_output=True, text=True, timeout=15, check=True
    )


def _pr_is_draft(pr_url: str) -> bool:
    res = _gh("pr", "view", pr_url, "--json", "isDraft")
    return bool(json.loads(res.stdout).get("isDraft"))


def _add_to_board(pr_url: str) -> str | None:
    """Idempotently add the PR to project #9; return its project-item id.

    `gh project item-add` no-ops server-side if the PR is already on the board
    (the project dedupes by content) and returns the existing item either way,
    so this is safe to call unconditionally.
    """
    res = _gh(
        "project", "item-add", str(PROJECT_NUMBER),
        "--owner", PROJECT_OWNER, "--url", pr_url, "--format", "json",
    )
    return json.loads(res.stdout).get("id")


def _set_status(item_id: str, option_id: str) -> None:
    _gh(
        "api", "graphql", "-f",
        "query=mutation($p:ID!,$i:ID!,$f:ID!,$o:String!){"
        "updateProjectV2ItemFieldValue(input:{projectId:$p,itemId:$i,fieldId:$f,"
        "value:{singleSelectOptionId:$o}}){projectV2Item{id}}}",
        "-f", f"p={PROJECT_ID}", "-f", f"i={item_id}",
        "-f", f"f={STATUS_FIELD_ID}", "-f", f"o={option_id}",
    )


def _log_fire(number: int, repo: str, status_label: str) -> None:
    """Append one fire-log line (Issue #453 infra). Never raises.

    Uses a `pr` key (not the base schema's `issue`) because this hook fires on a
    PR, not an Issue; `scripts/check_hook_health.sh` aggregates only on `hook`
    and `ts`, so the field name is free.
    """
    payload = {
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "hook": "post_gh_pr_create",
        "pr": number,
        "repo": repo,
        "action": f"board-add+status:{status_label}",
    }
    line = json.dumps(payload, separators=(",", ":")) + "\n"
    try:
        LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
        with open(LOG_PATH, "a", encoding="utf-8") as f:
            f.write(line)
    except OSError:
        pass


# --- orchestration ---


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0  # fail open

    cmd = (payload.get("tool_input") or {}).get("command", "")
    if not matches_pr_create(cmd):
        return 0

    parsed = parse_pr_url(extract_output(payload.get("tool_response")))
    if parsed is None:
        return 0  # create failed / --web / --dry-run — nothing to board
    owner, repo, number = parsed
    if not should_track(owner, repo):
        return 0

    pr_url = f"https://github.com/{owner}/{repo}/pull/{number}"
    try:
        label, option = status_for_draft(_pr_is_draft(pr_url))
        item_id = _add_to_board(pr_url)
        if not item_id:
            return 0
        _set_status(item_id, option)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired,
            json.JSONDecodeError, FileNotFoundError):
        return 0  # fail open — never break the user's flow on a gh hiccup

    _log_fire(number, f"{owner}/{repo}", label)
    print(json.dumps({
        "hookSpecificOutput": {
            "hookEventName": "PostToolUse",
            "additionalContext": (
                f"✅ post_gh_pr_create: PR #{number} ({owner}/{repo}) added to "
                f"board #{PROJECT_NUMBER} + Status → {label}. Remaining manual "
                f"checklist steps: `**Created by:**` body attribution + mirror "
                f"the linked Issue's Status on review request."
            ),
        }
    }))
    return 0


if __name__ == "__main__":
    sys.exit(main())
