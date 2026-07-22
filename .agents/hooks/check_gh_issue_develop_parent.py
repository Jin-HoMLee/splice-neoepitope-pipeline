#!/usr/bin/env python3
"""Pre-flight hook: WARN on `gh issue develop <N>` when Issue N is a parent/epic
(has ≥1 sub-issue).

Post-Issue #1155 this is ADVISORY (`ask`, warn-and-confirm), not a block. The real
enforcement moved to a merge-time gate (tools/ci/parent_child_gate.py) that refuses
a merge whose PR closes a parent with open children. Rationale: the harm (a PR
auto-closing a parent and orphaning its open children - PR #543 → Issue #538,
2026-05-28) lands at MERGE, so the branch-time ban was over-broad, banning the
safe, reversible act of branching to prevent an unsafe, irreversible one. This
hook now warns and points at the merge-time gate; branching proceeds on confirm.

Sibling of `.agents/hooks/check_at_claude.py` (the @claude mention guard) and
`check_gh_issue_create_repo.py` (the repo-guard, also an `ask`). Reads PreToolUse
hook JSON on stdin, prints an `ask` decision on stdout when the guard fires, exits
0 silently otherwise (the harness treats no-output as allow).

Fails OPEN on any parse miss, missing issue number, untokenizable command, or
`gh` error/timeout — a guard hook must never break the user's flow on a hiccup.
The only path that warns is a *confirmed* parent (subIssuesSummary.total > 0).

See https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/549 (original
block) and https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1155
(downgrade to advisory + merge-time gate).
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
import graphql_meter  # noqa: E402

LOG_PATH = Path(__file__).resolve().parent.parent.parent / ".agents" / "hook_fires.jsonl"

_PUNCT = _shell_parse.PUNCT  # single-sourced separator set (Issue #1130 review)
_DEVELOP_PREFIX = ("gh", "issue", "develop")
# `gh issue develop` flags that consume the following token as their value.
_VALUE_FLAGS = {"-b", "--base", "--branch-repo", "-n", "--name", "-R", "--repo"}
_ISSUE_URL_RE = re.compile(r"/issues/(\d+)")


# --- pure helpers (unit-tested, no I/O) ---


def develop_args(cmd: str) -> list[str] | None:
    """Return the token list AFTER a real `gh issue develop` invocation, or None.

    Tokenizes with shlex (quotes + shell punctuation respected) so a literal
    "gh issue develop" buried inside a quoted argument — e.g. a PR-comment body
    discussing this very hook — does NOT match. Only a `gh issue develop` at a
    command start (start of line or after a `&&`/`||`/`;`/`|` separator) counts.
    Args are collected up to the next shell separator, so a trailing
    `&& other-cmd` doesn't leak its tokens in. Untokenizable input (unbalanced
    quotes) fails safe → None.

    The command is **normalized first** (Issue #1142): heredoc bodies are stripped
    and unquoted newlines become separators. Without that, a `gh issue develop` on
    the second line of a multi-line command sat after an ordinary word rather than
    a separator, so `at_command_start` was False and the parent guard silently
    allowed branching off an epic - the exact auto-close footgun it exists to stop.
    A newline IS a command separator in shell; shlex does not treat it as one.
    Same defect class as `matches_pr_create` (Issue #1130), same shared fix.
    """
    try:
        lex = shlex.shlex(_shell_parse.normalize_command(cmd),
                          posix=True, punctuation_chars=True)
        lex.whitespace_split = True
        tokens = list(lex)
    except ValueError:
        return None
    at_command_start = True
    for i, tok in enumerate(tokens):
        if tok and all(ch in _PUNCT for ch in tok):  # pure-punctuation = separator
            at_command_start = True
            continue
        if at_command_start and tuple(tokens[i:i + 3]) == _DEVELOP_PREFIX:
            rest = tokens[i + 3:]
            args: list[str] = []
            for t in rest:
                if t and all(ch in _PUNCT for ch in t):
                    break  # next shell command begins — stop collecting
                args.append(t)
            return args
        at_command_start = False
    return None


def issue_number(args: list[str]) -> int | None:
    """Extract the Issue selector (number or issue-URL) from `gh issue develop` args.

    Skips value-taking flags and their values; the first positional token is the
    selector. Accepts a bare integer or a `.../issues/N` URL. Returns None when no
    selector is parseable (→ fail open).
    """
    skip_next = False
    for tok in args:
        if skip_next:
            skip_next = False
            continue
        if tok.startswith("-"):
            # `--name foo` consumes the next token; `--name=foo` is self-contained.
            if tok in _VALUE_FLAGS:
                skip_next = True
            continue
        if tok.isdigit():
            return int(tok)
        m = _ISSUE_URL_RE.search(tok)
        if m:
            return int(m.group(1))
        return None  # first positional isn't a recognizable issue selector
    return None


def repo_from_args(args: list[str]) -> tuple[str, str] | None:
    """Parse an explicit `-R`/`--repo [HOST/]OWNER/REPO` from the develop args.

    Returns (owner, name) when present and well-formed, else None (caller then
    falls back to the cwd repo). Strips an optional leading HOST/ segment.
    """
    val = None
    for i, tok in enumerate(args):
        if tok in ("-R", "--repo"):
            val = args[i + 1] if i + 1 < len(args) else None
            break
        if tok.startswith("-R=") or tok.startswith("--repo="):
            val = tok.split("=", 1)[1]
            break
    if not val:
        return None
    parts = val.split("/")
    if len(parts) >= 2:
        return parts[-2], parts[-1]
    return None


def warn_payload(number: int, total: int) -> dict:
    """The PreToolUse advisory (warn-and-confirm) decision for a confirmed parent.

    Post-Issue #1155 this is `ask`, not `deny`: branching off a parent is a safe,
    reversible act, so the branch-time guard is advisory and the real enforcement
    is the merge-time gate (tools/ci/parent_child_gate.py), which blocks a merge
    that would auto-close a parent with open children. The message points there.
    """
    return {
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "ask",
            "permissionDecisionReason": (
                f"Issue #{number} is a parent/epic ({total} sub-issue(s)). "
                "Branching off it is allowed - but the MERGE will be blocked "
                "until all its children are closed, because merging a PR that "
                "closes a parent auto-closes the epic and orphans its open "
                "children (incident: PR #543 accidentally closed epic Issue #538 "
                "on merge, 2026-05-28). The merge-time gate "
                "(tools/ci/parent_child_gate.py) is the real enforcement now "
                "(Issue #1155). Confirm to branch, or branch off a leaf "
                "sub-issue instead. Rule: memory/shared/feedback_parent_sub_issues.md."
            ),
        }
    }


# --- gh I/O (fail-open) ---


def _gh(*args: str, timeout: int = 8) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["gh", *args], capture_output=True, text=True, timeout=timeout, check=True
    )


def resolve_repo(args: list[str]) -> tuple[str, str] | None:
    """(owner, name) for the develop target: explicit `-R` wins, else cwd repo."""
    explicit = repo_from_args(args)
    if explicit:
        return explicit
    try:
        res = _gh("repo", "view", "--json", "nameWithOwner", "-q", ".nameWithOwner",
                  timeout=6)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired,
            FileNotFoundError):
        return None
    owner, _, name = res.stdout.strip().partition("/")
    return (owner, name) if owner and name else None


def sub_issue_total(owner: str, name: str, number: int) -> int | None:
    """subIssuesSummary.total for the issue, or None on any error (→ fail open)."""
    query = (
        "query($o:String!,$n:String!,$num:Int!){"
        "rateLimit{cost remaining}"
        "repository(owner:$o,name:$n){issue(number:$num){"
        "subIssuesSummary{total}}}}"
    )
    try:
        res = _gh("api", "graphql", "-f", f"query={query}",
                  "-f", f"o={owner}", "-f", f"n={name}", "-F", f"num={number}")
        data = json.loads(res.stdout)
        graphql_meter.log_graphql_spend("check_gh_issue_develop_parent", data, query_name="parent_check")
        return int(data["data"]["repository"]["issue"]["subIssuesSummary"]["total"])
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired,
            FileNotFoundError, json.JSONDecodeError, KeyError, TypeError, ValueError):
        return None


def _log_fire(number: int, repo: str, total: int) -> None:
    """Append one fire-log line on a deny (Issue #453 infra). Never raises."""
    payload = {
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "hook": "check_gh_issue_develop_parent",
        "issue": number,
        "repo": repo,
        "action": f"warned-parent-develop:sub_total={total}",
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
    args = develop_args(cmd)
    if args is None:
        return 0  # not a `gh issue develop` invocation

    number = issue_number(args)
    if number is None:
        return 0  # no parseable issue selector → fail open

    repo = resolve_repo(args)
    if repo is None:
        return 0  # can't resolve target repo → fail open
    owner, name = repo

    total = sub_issue_total(owner, name, number)
    if total is None or total == 0:
        return 0  # query failed (fail open) or a confirmed leaf → allow

    _log_fire(number, f"{owner}/{name}", total)
    print(json.dumps(warn_payload(number, total)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
