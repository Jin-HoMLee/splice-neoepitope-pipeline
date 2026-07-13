#!/usr/bin/env python3
"""Pre-flight hook: refuse a single Bash command that chains git commit -> push
or push -> `gh pr merge`.

House rule "Commit, push, merge - three separate steps": each is the last
interception point before a less-reversible state (commit = local, push =
external, merge = permanent). Chaining them (`git commit -m x && git push`,
`git push && gh pr merge`) collapses the user's review gate between steps - the
whole point of keeping them separate is that the user can inspect the result of
one before authorizing the next.

WHY A HOOK, NOT `permissions.deny`: Claude Code splits a compound command on
shell operators and matches each subcommand independently (per the permissions
docs), so a static deny on `git push` cannot distinguish "push chained after
commit" from a standalone "push" - it would block every push. Only a hook that
sees the whole command can flag the *chain* while allowing each step alone.
Sibling of `check_no_force_push.py` / `check_no_cd_outside_cwd.py`.

WHAT FIRES: within ONE Bash command, an ordered pair across separate subcommands -
a `git commit` that precedes a `git push`, or a `git push` that precedes a
`gh pr merge`. Each step ALONE, and unrelated chains (`git add && git commit`,
`git fetch && git push`), pass untouched. shlex quoting collapses a quoted string
to one token, so `git commit -m "push and merge"` carries no bare push/merge
token and does NOT false-positive.

ESCAPE HATCH: set env `CLAUDE_ALLOW_CHAINED_GIT=1` for a deliberate one-liner the
user has authorized.

Reads PreToolUse hook JSON on stdin, prints a deny decision on stdout when the
guard fires, exits 0 silently otherwise. Fails OPEN on any parse miss /
untokenizable command. Fires are appended to `.agents/hook_fires.jsonl`
(Issue #453 infra).

See https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/626.
"""
from __future__ import annotations

import json
import os
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _shell_parse  # noqa: E402

LOG_PATH = Path(__file__).resolve().parent.parent.parent / ".agents" / "hook_fires.jsonl"

_PUNCT = _shell_parse.PUNCT  # single-sourced separator set (Issue #1130 review)
_TRUTHY = {"1", "true", "yes", "on"}


# --- pure helpers (unit-tested, no I/O) ---


def split_subcommands(cmd: str) -> list[list[str]] | None:
    """Tokenize `cmd` and split into subcommand token-lists on shell separators.

    Returns a list of subcommands (each a list of tokens), or None on a
    tokenization failure (unbalanced quotes) -> caller fails open. Quoted strings
    survive as single tokens; `&&`/`||`/`;`/`|`/`&`/`(`/`)` arrive as standalone
    separator tokens.

    The command is **normalized first** (Issue #1142): heredoc bodies are stripped
    and unquoted newlines become separators. A newline is a command separator in
    shell but not in shlex, so without this `git commit -m x` and `git push` written
    on two lines MERGED into one token block, the commit and push indices collapsed
    to the same position, and the guard silently allowed the exact chain it exists
    to refuse. Only the `&&` spelling was ever caught - writing the two commands on
    separate lines walked straight through.
    """
    try:
        lex = shlex.shlex(_shell_parse.normalize_command(cmd),
                          posix=True, punctuation_chars=True)
        lex.whitespace_split = True
        tokens = list(lex)
    except ValueError:
        return None
    subcommands: list[list[str]] = []
    current: list[str] = []
    for tok in tokens:
        if tok and all(ch in _PUNCT for ch in tok):
            if current:
                subcommands.append(current)
                current = []
            continue
        current.append(tok)
    if current:
        subcommands.append(current)
    return subcommands


def is_git_commit(tokens: list[str]) -> bool:
    return "git" in tokens and "commit" in tokens


def is_git_push(tokens: list[str]) -> bool:
    return "git" in tokens and "push" in tokens


def is_gh_pr_merge(tokens: list[str]) -> bool:
    return "gh" in tokens and "pr" in tokens and "merge" in tokens


def chained_violation(cmd: str) -> str | None:
    """Return a label for the dangerous chain, or None.

    `commit->push`  : a `git commit` subcommand precedes a `git push` subcommand.
    `push->merge`   : a `git push` subcommand precedes a `gh pr merge` subcommand.
    Ordering uses subcommand position, so each step alone (a single subcommand)
    can never trigger. False on a parse miss.
    """
    subs = split_subcommands(cmd)
    if subs is None:
        return None
    commit_idx = [i for i, s in enumerate(subs) if is_git_commit(s)]
    push_idx = [i for i, s in enumerate(subs) if is_git_push(s)]
    merge_idx = [i for i, s in enumerate(subs) if is_gh_pr_merge(s)]
    if commit_idx and push_idx and min(commit_idx) < max(push_idx):
        return "commit->push"
    if push_idx and merge_idx and min(push_idx) < max(merge_idx):
        return "push->merge"
    return None


def is_allowed(env=None) -> bool:
    """True when the user has set the CLAUDE_ALLOW_CHAINED_GIT escape hatch."""
    env = os.environ if env is None else env
    return (env.get("CLAUDE_ALLOW_CHAINED_GIT") or "").strip().lower() in _TRUTHY


def deny_payload(label: str) -> dict:
    """The PreToolUse deny decision for a chained git commit/push/merge."""
    pair = {
        "commit->push": "`git commit` and `git push`",
        "push->merge": "`git push` and `gh pr merge`",
    }.get(label, "these git steps")
    return {
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "deny",
            "permissionDecisionReason": (
                f"Refusing to chain {pair} in one command. House rule: commit, "
                "push, merge are three separate steps - each is the last "
                "interception point before a less-reversible state (commit local, "
                "push external, merge permanent), and chaining collapses the review "
                "gate between them. Run them as separate commands so each result "
                "can be inspected first. Deliberate one-liner the user authorized: "
                "set CLAUDE_ALLOW_CHAINED_GIT=1."
            ),
        }
    }


# --- I/O ---


def _log_fire(label: str, snippet: str) -> None:
    """Append one fire-log line on a deny (Issue #453 infra). Never raises."""
    payload = {
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "hook": "check_commit_push_separation",
        "action": "refused-chained-git:" + label,
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

    if data.get("tool_name") != "Bash":
        return 0

    if is_allowed():
        return 0

    cmd = (data.get("tool_input") or {}).get("command", "")
    if not cmd:
        return 0
    label = chained_violation(cmd)
    if label is None:
        return 0

    _log_fire(label, cmd)
    print(json.dumps(deny_payload(label)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
