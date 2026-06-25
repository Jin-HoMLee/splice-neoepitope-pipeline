#!/usr/bin/env python3
"""Pre-flight hook: refuse an unpaginated ProjectV2 board GraphQL query.

The lab project board (#9) holds 700+ items and sorts Done first, so a single-page
`projectV2 { items(first: N) }` query silently returns mostly-closed items and
*hides the open work past position N* — Ready / In-progress items vanish with no
error. That exact truncation produced a wrong "Ready is empty" read on 2026-06-12
(this guard's trigger incident): a `first: 100` query showed only the first 14%
of 711 items, hiding the real Scientist Ready queue. The rule to use the
paginating helper lives in memory (shared/feedback_board_queries.md) but keeps
slipping at the moment someone types a quick query — so this is the
mechanism-over-memory escalation.

The correct path already exists: `scripts/board_open_items.py` loops on
`pageInfo.hasNextPage` and handles the Done-first sort. For a hand-rolled query, a
`pageInfo { hasNextPage endCursor }` + `after:` cursor loop is required. For a
single-issue lookup, query `issue(number:N){ projectItems … }` instead of
scanning the whole board.

Sibling of `check_gh_issue_develop_parent.py` / `check_at_claude.py`. Reads
PreToolUse hook JSON on stdin, prints a deny decision on stdout when the guard
fires, exits 0 silently otherwise (the harness treats no-output as allow). This
guard does NO `gh`/network I/O — it is pure string inspection of the command.

Fires ONLY on a real command-start `gh api …` invocation whose collected args
contain `projectV2` + an `items(first: …)` connection but NO pagination token
(`hasNextPage` / `endCursor` / `after:`). So the two safe patterns pass untouched:
the per-issue `issue(number:N){ projectItems … }` lookup (no `projectV2.items`)
and `board_open_items.py` itself (carries pageInfo/after; also a python
invocation, never inspected). A `gh pr comment` / `gh issue create` whose *body*
merely discusses the pattern does NOT match — only a `gh api` at a command start
counts (shlex command-start tokenization).

Fails OPEN on any parse miss or untokenizable command — a guard must never break
the user's flow on a hiccup.

See https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/717.
"""
from __future__ import annotations

import json
import re
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path

LOG_PATH = Path(__file__).resolve().parent.parent.parent / ".agents" / "hook_fires.jsonl"

_PUNCT = set("();<>|&")  # shell punctuation_chars → standalone separator tokens
_API_PREFIX = ("gh", "api")
# An `items(first: …)` connection, tolerant of whitespace/newlines in the query.
_ITEMS_FIRST_RE = re.compile(r"items\s*\(\s*first\s*:", re.IGNORECASE)
# Any of these reads as a real cursor loop → the query is paginated → allow.
# The `(?<!\$)` lookbehind makes the `after:` arm match cursor *usage*
# (`items(first:N, after: $c)`) but NOT a variable *declaration* (`$after:String`)
# — so a query that declares `$after` yet never wires it into items() is correctly
# still treated as unpaginated. hasNextPage/endCursor cover the common cursor loop.
_PAGINATION_RE = re.compile(r"hasNextPage|endCursor|(?<!\$)\bafter\s*:", re.IGNORECASE)


# --- pure helpers (unit-tested, no I/O) ---


def api_args(cmd: str) -> str | None:
    """Return the joined token string AFTER a command-start `gh api …`, or None.

    Tokenizes with shlex (quotes + shell punctuation respected) so a literal
    "gh api" buried inside a quoted argument — e.g. a PR-comment body discussing
    this guard — does NOT match. Only a `gh api` at a command start (start of line
    or after a `&&`/`||`/`;`/`|` separator) counts. Args are collected up to the
    next shell separator, so a trailing `&& other-cmd` doesn't leak its tokens in.

    The GraphQL query rides inside `-f query='…'`, so its internal `(){};:` are
    quote-protected and arrive as a single token — the surrounding quotes keep the
    punctuation_chars splitter from breaking the query apart. Untokenizable input
    (unbalanced quotes) fails safe → None.
    """
    try:
        lex = shlex.shlex(cmd or "", posix=True, punctuation_chars=True)
        lex.whitespace_split = True
        tokens = list(lex)
    except ValueError:
        return None
    at_command_start = True
    for i, tok in enumerate(tokens):
        if tok and all(ch in _PUNCT for ch in tok):  # pure-punctuation = separator
            at_command_start = True
            continue
        if at_command_start and tuple(tokens[i:i + 2]) == _API_PREFIX:
            rest = tokens[i + 2:]
            args: list[str] = []
            for t in rest:
                if t and all(ch in _PUNCT for ch in t):
                    break  # next shell command begins — stop collecting
                args.append(t)
            return " ".join(args)
        at_command_start = False
    return None


def is_unpaginated_board_query(args: str) -> bool:
    """True iff the collected `gh api` args are a board query that will truncate.

    Dangerous shape: a `projectV2` query with an `items(first: …)` connection and
    NO pagination token. The presence of any pagination token (hasNextPage /
    endCursor / after:) is read as a real cursor loop → safe → allow. The bias is
    deliberately toward allow: a guard should never over-block a query that already
    looks paginated.
    """
    if "projectV2" not in args:
        return False
    if not _ITEMS_FIRST_RE.search(args):
        return False
    if _PAGINATION_RE.search(args):
        return False
    return True


def deny_payload() -> dict:
    """The PreToolUse deny decision for an unpaginated board query."""
    return {
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "deny",
            "permissionDecisionReason": (
                "Refusing an unpaginated ProjectV2 board query: this "
                "`projectV2 { items(first: …) }` has no pagination cursor "
                "(hasNextPage / endCursor / after:). The board (#9) holds 700+ "
                "items sorted Done-first, so a single page silently hides the open "
                "Ready / In-progress work past the cap (it read 'Ready is empty' "
                "wrongly on 2026-06-12). Use `scripts/board_open_items.py` "
                "(--role / --status / --json; it paginates and handles the "
                "Done-first sort), or for a single-issue lookup query "
                "`issue(number:N){ projectItems … }` instead of scanning the whole "
                "board. To roll your own, add a `pageInfo { hasNextPage endCursor "
                "}` + `after:` cursor loop. Rule: "
                "memory/shared/feedback_board_queries.md."
            ),
        }
    }


def _log_fire(snippet: str) -> None:
    """Append one fire-log line on a deny (Issue #453 infra). Never raises."""
    payload = {
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "hook": "check_board_query_pagination",
        "action": "refused-unpaginated-board-query",
        "snippet": snippet[:200],
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
    args = api_args(cmd)
    if args is None:
        return 0  # not a command-start `gh api` invocation
    if not is_unpaginated_board_query(args):
        return 0  # not the dangerous shape → allow

    _log_fire(args)
    print(json.dumps(deny_payload()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
