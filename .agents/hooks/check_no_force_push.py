#!/usr/bin/env python3
"""Pre-flight hook: refuse a `git push` carrying a force flag.

Force-pushing rewrites published history: it can clobber a teammate's commits on
a shared branch and is irreversible on the remote. The house rule is "never
force-push without explicit user confirmation", but an agent cannot ask
mid-tool-call, so this guard denies the push and points at the escape hatch for
the case where the user HAS authorized it.

WHY A HOOK, NOT `permissions.deny`: the official permissions docs call
argument-constrained static denies "fragile" - a pattern like
`Bash(git push --force:*)` misses `git push origin main --force` (flag last),
`-f`, and `--force-with-lease`, and the docs explicitly recommend a PreToolUse
hook for argument filtering. A hook parses the whole command and catches the
force flag in any position. Sibling of the other destructive-command guards
`check_commit_push_separation.py` / `check_no_cd_outside_cwd.py` and of the
existing `check_board_query_pagination.py` / `check_no_emdash.py`.

SUBCOMMAND-SCOPED: the command is split on shell separators
(`&&`/`||`/`;`/`|`/newline) and each subcommand is inspected on its own, so a
force flag in a *different* subcommand is never mis-attributed to the push
(`git push && lint --force` does not fire). shlex quoting collapses a quoted
string to a single token, so `git commit -m "push --force"` carries no bare
`push` / `--force` token and does NOT false-positive.

ESCAPE HATCH: set env `CLAUDE_ALLOW_FORCE_PUSH=1` when the user has explicitly
authorized a force-push (e.g. reflog recovery, a personal branch rewrite).

Reads PreToolUse hook JSON on stdin, prints a deny decision on stdout when the
guard fires, exits 0 silently otherwise (the harness treats no-output as allow).
Fails OPEN on any parse miss / untokenizable command - a guard must never break
the user's flow on a hiccup. Fires are appended to `.agents/hook_fires.jsonl`
(Issue #453 infra).

See https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/626.
"""
from __future__ import annotations

import json
import os
import re
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path

LOG_PATH = Path(__file__).resolve().parent.parent.parent / ".agents" / "hook_fires.jsonl"

_PUNCT = set("();<>|&")  # shell punctuation_chars -> standalone separator tokens
_TRUTHY = {"1", "true", "yes", "on"}

# A single-dash short-flag bundle containing 'f' (e.g. -f, -fq, -qf). Double-dash
# long flags (--force*) are matched separately; the bundle regex never matches
# them because the char after the leading '-' must be a letter, not another '-'.
_SHORT_FORCE_RE = re.compile(r"^-[a-zA-Z]*f[a-zA-Z]*$")


# --- pure helpers (unit-tested, no I/O) ---


def split_subcommands(cmd: str) -> list[list[str]] | None:
    """Tokenize `cmd` and split into subcommand token-lists on shell separators.

    Returns a list of subcommands (each a list of tokens), or None when the
    command cannot be tokenized (unbalanced quotes) -> caller fails open. Uses
    shlex with punctuation_chars so quoted strings survive as single tokens and
    `&&`/`||`/`;`/`|`/`&`/`(`/`)` arrive as standalone separator tokens.
    """
    try:
        lex = shlex.shlex(cmd or "", posix=True, punctuation_chars=True)
        lex.whitespace_split = True
        tokens = list(lex)
    except ValueError:
        return None
    subcommands: list[list[str]] = []
    current: list[str] = []
    for tok in tokens:
        if tok and all(ch in _PUNCT for ch in tok):  # pure-punctuation = separator
            if current:
                subcommands.append(current)
                current = []
            continue
        current.append(tok)
    if current:
        subcommands.append(current)
    return subcommands


def _is_force_flag(tok: str) -> bool:
    """True iff a single token is a git-push force flag (any accepted spelling)."""
    if tok in ("-f", "--force"):
        return True
    if tok.startswith("--force-with-lease") or tok.startswith("--force-if-includes"):
        return True
    return bool(_SHORT_FORCE_RE.match(tok))


def is_force_push(tokens: list[str]) -> bool:
    """True iff a subcommand is `git push` with a force flag.

    Requires the bare tokens `git` and `push` both present (a quoted string like
    "push --force" collapses to one token and so cannot trip this) and at least
    one force flag among the tokens. Requiring `push` as its own token keeps a
    `git log --grep push` style command from matching.
    """
    if "git" not in tokens or "push" not in tokens:
        return False
    return any(_is_force_flag(t) for t in tokens)


def command_force_pushes(cmd: str) -> bool:
    """True iff any subcommand of `cmd` is a force `git push`. False on parse miss."""
    subs = split_subcommands(cmd)
    if subs is None:
        return False
    return any(is_force_push(sub) for sub in subs)


def is_allowed(env=None) -> bool:
    """True when the user has set the CLAUDE_ALLOW_FORCE_PUSH escape hatch."""
    env = os.environ if env is None else env
    return (env.get("CLAUDE_ALLOW_FORCE_PUSH") or "").strip().lower() in _TRUTHY


def deny_payload() -> dict:
    """The PreToolUse deny decision for a force `git push`."""
    return {
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "deny",
            "permissionDecisionReason": (
                "Refusing a force `git push`: it rewrites published history and can "
                "irreversibly clobber commits on a shared branch. House rule: never "
                "force-push without explicit user confirmation. If the user HAS "
                "authorized this (reflog recovery, a personal-branch rewrite), prefer "
                "`--force-with-lease` and set CLAUDE_ALLOW_FORCE_PUSH=1 for this "
                "session. Otherwise confirm with the user first."
            ),
        }
    }


# --- I/O ---


def _log_fire(snippet: str) -> None:
    """Append one fire-log line on a deny (Issue #453 infra). Never raises."""
    payload = {
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "hook": "check_no_force_push",
        "action": "refused-force-push",
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
        return 0  # only guards Bash commands (settings matcher also scopes this)

    if is_allowed():
        return 0

    cmd = (data.get("tool_input") or {}).get("command", "")
    if not cmd or not command_force_pushes(cmd):
        return 0  # not a force push -> allow

    _log_fire(cmd)
    print(json.dumps(deny_payload()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
