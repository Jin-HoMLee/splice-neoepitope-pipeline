#!/usr/bin/env python3
"""Pre-flight hook: refuse advancing a sub-issue's board Status past its parent.

The lab board (#9) mirrors a parent's Status from its sub-issue progress: the
canonical rule (memory/shared/feedback_parent_sub_issues.md) is that the *first
sub-issue to start* bumps the parent to In progress. So moving a sub-issue to an
*advanced* Status (In progress / Ready for review / In review) while its parent is
still in **Backlog** (or has **No Status**) is a discipline slip — the parent
should have been bumped first. It silently desyncs the board: an epic reads
"untouched" while its subs are flying (incident #480; the rule is Reference-tier
in memory and didn't load at the Status-mutation moment, so this is the
mechanism-over-memory escalation — Issue #499, companion to the source-agnostic
GitHub Action in Issue #498).

A board Status write is a `gh api graphql … updateProjectV2ItemFieldValue …`
mutation carrying an **itemId** (`PVTI_…`) and a **singleSelectOptionId**, not an
issue number — so the hook resolves itemId → issue → native parent → parent's
project-#9 Status, then applies the precondition. A single command may carry
several mutations (the PR+Issue lifecycle flip is one aliased dual-mutation); each
is checked independently and any one offending block denies.

Sibling of `check_gh_issue_develop_parent.py` / `check_board_query_pagination.py`.
Reads PreToolUse hook JSON on stdin, prints a deny decision on stdout when the
guard fires, exits 0 silently otherwise (the harness treats no-output as allow).

Fails OPEN on every uncertain path — parse miss, untokenizable command, not a
Status mutation, a non-advanced target (→ Backlog / Ready / Done), a PR item or
standalone issue (no parent), a parent not tracked on board #9, or any `gh`
error/timeout. The only path that denies is a *confirmed* sub-issue whose parent
is on board #9 in Backlog / No Status. A guard must never break the flow on a
hiccup.

See https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/499.
"""
from __future__ import annotations

import json
import re
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

LOG_PATH = Path(__file__).resolve().parent.parent.parent / ".claude" / "hook_fires.jsonl"

_PUNCT = set("();<>|&")  # shell punctuation_chars → standalone separator tokens
_API_PREFIX = ("gh", "api")

# Project #9 ("JH M Lee Lab") Status field + the option IDs of its single-select.
PROJECT_ID = "PVT_kwHOB17eGc4BSomP"
STATUS_FIELD_ID = "PVTSSF_lAHOB17eGc4BSomPzhAHFf8"
BACKLOG = "f75ad846"
IN_PROGRESS = "47fc9ee4"
READY_FOR_REVIEW = "8bf9192f"
IN_REVIEW = "df73e18b"
# Moving a sub to any of these requires the parent to be at least In progress.
ADVANCED_TARGETS = {IN_PROGRESS, READY_FOR_REVIEW, IN_REVIEW}
# Parent Status values that block a sub from advancing (epic not yet "in flight").
# None == the parent is on the board with no Status value set (No Status).
BLOCKING_PARENT_STATUS = {None, BACKLOG}

_STATUS_OPTION_NAMES = {
    BACKLOG: "Backlog",
    "61e4505c": "Ready",
    IN_PROGRESS: "In progress",
    READY_FOR_REVIEW: "Ready for review",
    IN_REVIEW: "In review",
    "98236657": "Done",
    None: "No Status",
}

# Each updateProjectV2ItemFieldValue block carries itemId / fieldId / optionId.
_ITEM_ID_RE = re.compile(r'itemId\s*:\s*\\?["\']([^"\'\\]+)')
_FIELD_ID_RE = re.compile(r'fieldId\s*:\s*\\?["\']([^"\'\\]+)')
_OPTION_ID_RE = re.compile(r'singleSelectOptionId\s*:\s*\\?["\']([^"\'\\]+)')


# --- pure helpers (unit-tested, no I/O) ---


def _is_command_start_gh_api(cmd: str) -> bool:
    """True iff a `gh api` invocation appears at a command start.

    Tokenizes with shlex (quotes + shell punctuation respected) so a literal
    "gh api" buried inside a quoted argument — e.g. a PR-comment body discussing
    this hook — does NOT count. Only a `gh api` at a command start (start of line
    or after a `&&`/`||`/`;`/`|` separator) qualifies. Untokenizable input
    (unbalanced quotes) fails safe → False.
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
        if at_command_start and tuple(tokens[i:i + 2]) == _API_PREFIX:
            return True
        at_command_start = False
    return False


def status_mutations(cmd: str) -> list[tuple[str, str]]:
    """Extract (itemId, optionId) for each Status-field write to an advanced target.

    Returns [] unless `cmd` is a command-start `gh api …` invocation. Splits on
    each `updateProjectV2ItemFieldValue` so an aliased dual-mutation (the PR+Issue
    lifecycle flip) yields one tuple per block. A block is kept only when its
    fieldId is the board Status field AND its target optionId is in ADVANCED_TARGETS
    — writes to Backlog / Ready / Done, or to other fields, carry no parent
    precondition and are screened out here (so the gh resolver is never consulted
    for them).
    """
    if not _is_command_start_gh_api(cmd):
        return []
    out: list[tuple[str, str]] = []
    blocks = cmd.split("updateProjectV2ItemFieldValue")
    for block in blocks[1:]:  # text before the first marker has no mutation
        field = _FIELD_ID_RE.search(block)
        if not field or field.group(1) != STATUS_FIELD_ID:
            continue
        item = _ITEM_ID_RE.search(block)
        option = _OPTION_ID_RE.search(block)
        if not item or not option:
            continue
        if option.group(1) not in ADVANCED_TARGETS:
            continue
        out.append((item.group(1), option.group(1)))
    return out


def should_refuse(parent_status_option_id: str | None) -> bool:
    """True iff the parent's Status blocks a sub from advancing (Backlog / No Status)."""
    return parent_status_option_id in BLOCKING_PARENT_STATUS


def deny_payload(sub_number: int, parent_number: int,
                 parent_status_name: str, target_name: str) -> dict:
    """The PreToolUse deny decision for a sub whose parent isn't yet in flight."""
    return {
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "deny",
            "permissionDecisionReason": (
                f"Refusing to set sub-issue #{sub_number} Status → '{target_name}': "
                f"its parent #{parent_number} is still '{parent_status_name}'. The "
                "first sub-issue to start must bump the parent to 'In progress' "
                f"first (set #{parent_number} → In progress, then retry). Otherwise "
                "the board desyncs — the epic reads untouched while its subs are "
                "advancing. Rule: memory/shared/feedback_parent_sub_issues.md "
                "(parent Status mirrors sub-issue progress)."
            ),
        }
    }


# --- gh I/O (fail-open) ---


def _gh(*args: str, timeout: int = 8) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["gh", *args], capture_output=True, text=True, timeout=timeout, check=True
    )


def resolve_parent_status(item_id: str) -> tuple[int, int, str | None] | None:
    """(sub_number, parent_number, parent_status_optionId) for a board item, or None.

    Resolves the board item → its issue content → the issue's native parent →
    the parent's Status value on project #9. Returns None (→ fail open / allow)
    when the item's content is not an Issue (e.g. a PR), the issue has no parent
    (standalone issue or it *is* a parent), the parent isn't tracked on board #9,
    or any `gh`/parse error occurs. When the parent IS on board #9 but has no
    Status value set, the optionId is None (caller treats that as No Status →
    blocking).
    """
    query = (
        "query($id:ID!){node(id:$id){... on ProjectV2Item{content{... on Issue{"
        "number parent{number projectItems(first:20){nodes{project{id}"
        'fieldValueByName(name:"Status"){... on ProjectV2ItemFieldSingleSelectValue'
        "{optionId}}}}}}}}}}"
    )
    try:
        res = _gh("api", "graphql", "-f", f"query={query}", "-f", f"id={item_id}")
        content = json.loads(res.stdout)["data"]["node"]["content"]
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired,
            FileNotFoundError, json.JSONDecodeError, KeyError, TypeError, ValueError):
        return None
    if not content:  # node had no Issue content (e.g. a PullRequest item)
        return None
    parent = content.get("parent")
    if not parent:  # standalone issue or the item is itself a parent → not a sub
        return None
    sub_number = content.get("number")
    parent_number = parent.get("number")
    nodes = (parent.get("projectItems") or {}).get("nodes") or []
    for node in nodes:
        if (node.get("project") or {}).get("id") == PROJECT_ID:
            fv = node.get("fieldValueByName") or {}
            return sub_number, parent_number, fv.get("optionId")
    return None  # parent not tracked on board #9 → can't determine → fail open


def _log_fire(sub_number: int, parent_number: int, target: str) -> None:
    """Append one fire-log line on a deny (Issue #453 infra). Never raises."""
    payload = {
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "hook": "check_parent_status_precondition",
        "sub": sub_number,
        "parent": parent_number,
        "action": f"refused-sub-advance:target={target}",
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
        data = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0  # fail open

    cmd = (data.get("tool_input") or {}).get("command", "")
    mutations = status_mutations(cmd)
    if not mutations:
        return 0  # not a sub-advancing Status mutation → allow

    for item_id, target in mutations:
        resolved = resolve_parent_status(item_id)
        if resolved is None:
            continue  # not a sub / unresolvable → allow this block
        sub_number, parent_number, parent_status = resolved
        if should_refuse(parent_status):
            _log_fire(sub_number, parent_number, target)
            print(json.dumps(deny_payload(
                sub_number, parent_number,
                _STATUS_OPTION_NAMES.get(parent_status, "Backlog/No Status"),
                _STATUS_OPTION_NAMES.get(target, target))))
            return 0  # deny on the first offending block
    return 0


if __name__ == "__main__":
    sys.exit(main())
