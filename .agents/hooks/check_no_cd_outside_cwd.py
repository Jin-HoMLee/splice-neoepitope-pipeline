#!/usr/bin/env python3
"""Pre-flight hook: refuse a `cd` that leaves the project worktree.

House rule "Never cd out of the project worktree": the Bash cwd persists across
tool calls, so a `cd /path/to/another-repo` silently changes the context for
every later command - the classic footgun that makes a git command hit the wrong
repo. The rule's alternative is `git -C <path>` for git on other repos and
absolute paths for file ops. This guard denies a `cd` whose target resolves
outside the worktree; `cd` into a *subdirectory* of the worktree is fine.

WHY A HOOK, NOT `permissions.deny`: a static `Bash(cd *)` deny is all-or-nothing
- it would block `cd subdir/` too. Only a hook can resolve the target path and
decide inside-vs-outside. Sibling of `check_no_force_push.py` /
`check_commit_push_separation.py`.

BOUNDARY: the worktree root (derived from this file's own location, so it is
symlink-safe and needs no env var) PLUS the ephemeral temp roots the agent
legitimately uses for scratch files (`/tmp`, `/private/tmp`, `/var/tmp`,
`/var/folders`, `$TMPDIR`). A `cd` resolving inside any of those is allowed;
anything else (another repo, `$HOME`, a sibling clone) is denied.

MULTI-CLONE HEADS-UP: this project runs as separate role clones + git worktrees,
and the Memory-Manager persona operates from the personas-repo cwd. A cross-tree
hop into a *sibling* clone (not under this worktree, not a temp root) is denied
by design - use `git -C <other-clone> ...` or set CLAUDE_ALLOW_CD for a routine
that genuinely must switch trees.

STATEFUL: chained `cd`s are walked in order against a running virtual cwd, so
`cd workflow && cd ..` (which nets back to the root) does NOT false-positive.

FAILS OPEN (allows) on anything it cannot resolve: `cd -` (previous dir),
`cd $VAR` / command substitution, `cd ~user`, an untokenizable command, or a
missing HOME. A guard must never break a legitimate command on ambiguity.

ESCAPE HATCH: set env `CLAUDE_ALLOW_CD=1` for a deliberate out-of-tree `cd` the
user authorized (prefer `git -C` / absolute paths instead).

Reads PreToolUse hook JSON on stdin, prints a deny decision on stdout when the
guard fires, exits 0 silently otherwise. Fires are appended to
`.agents/hook_fires.jsonl` (Issue #453 infra).

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

PROJECT_ROOT = str(Path(__file__).resolve().parent.parent.parent)
LOG_PATH = Path(PROJECT_ROOT) / ".agents" / "hook_fires.jsonl"

_PUNCT = _shell_parse.PUNCT  # single-sourced separator set (Issue #1130 review)
_TRUTHY = {"1", "true", "yes", "on"}
# Ephemeral roots the agent legitimately uses for scratch files.
_TEMP_ROOTS = ("/tmp", "/private/tmp", "/var/tmp", "/var/folders")
# `cd` builtin options that are NOT the target directory.
_CD_OPTS = {"-L", "-P", "-e", "-@"}


# --- pure helpers (unit-tested, no I/O) ---


def split_subcommands(cmd: str) -> list[list[str]] | None:
    """Tokenize `cmd` and split into subcommand token-lists on shell separators.

    Returns a list of subcommands, or None on a tokenization failure -> caller
    fails open. Quoted strings survive as single tokens.

    The command is **normalized first** (Issue #1142): heredoc bodies are stripped
    and unquoted newlines become separators. A newline is a command separator in
    shell but not in shlex, so without this two `cd`s on two lines MERGED into one
    token block and `cd_target` returned only the FIRST positional - so a benign
    `cd /tmp` on line 1 masked an escaping `cd /etc/secret` on line 2, and the
    guard never examined it.
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


def cd_target(tokens: list[str]):
    """For a `cd` subcommand, return its target: a path string, "" for a bare
    `cd` (goes to $HOME), or None if the subcommand is not a `cd` at all.

    `cd -` returns the literal "-" (previous dir - caller treats as unresolvable).
    Skips `cd` option flags (-L/-P/-e/-@).
    """
    if not tokens or tokens[0] != "cd":
        return None
    for t in tokens[1:]:
        if t == "-":
            return "-"
        if t in _CD_OPTS:
            continue
        if t.startswith("-"):
            continue  # unknown option flag
        return t
    return ""  # bare `cd` -> home


def resolve_cd(target: str, base, env: dict):
    """Resolve a `cd` target to an absolute normalized path, or None if it cannot
    be resolved (dynamic/ambiguous -> caller fails open).

    `base` is the running cwd for relative targets (None if unknown).
    """
    if target == "-":
        return None  # previous dir - unknowable
    if target == "":
        home = env.get("HOME")
        return os.path.normpath(home) if home else None
    if "$" in target or "`" in target:
        return None  # env var / command substitution - unknowable
    if target == "~" or target.startswith("~/"):
        home = env.get("HOME")
        if not home:
            return None
        p = home if target == "~" else os.path.join(home, target[2:])
        return os.path.normpath(p)
    if target.startswith("~"):
        return None  # ~user form - unknowable
    if os.path.isabs(target):
        return os.path.normpath(target)
    if base is None:
        return None  # relative target but running cwd is unknown
    return os.path.normpath(os.path.join(base, target))


def _safe_roots(project_root: str, env: dict) -> set:
    roots = {os.path.normpath(project_root)}
    roots.update(_TEMP_ROOTS)
    tmpdir = env.get("TMPDIR")
    if tmpdir:
        roots.add(os.path.normpath(tmpdir))
    return roots


def is_inside_safe(resolved: str, project_root: str, env: dict) -> bool:
    """True iff `resolved` is the worktree root, a temp root, or a subdir of one."""
    for root in _safe_roots(project_root, env):
        r = root.rstrip("/") or "/"
        if resolved == r or resolved.startswith(r + os.sep):
            return True
    return False


def first_escaping_cd(cmd: str, session_cwd, project_root: str, env: dict):
    """Return the target of the first `cd` that resolves outside the worktree
    (and temp roots), or None if no `cd` escapes.

    Walks chained subcommands with a running virtual cwd so `cd sub && cd ..`
    (net-inside) is not flagged. An unresolvable `cd` makes the running cwd
    unknown (subsequent relative `cd`s then fail open) but never itself denies.
    """
    subs = split_subcommands(cmd)
    if subs is None:
        return None
    running = os.path.normpath(session_cwd) if session_cwd else os.path.normpath(project_root)
    for sub in subs:
        target = cd_target(sub)
        if target is None:
            continue  # not a cd
        resolved = resolve_cd(target, running, env)
        if resolved is None:
            running = None  # unknown from here
            continue
        if not is_inside_safe(resolved, project_root, env):
            return target
        running = resolved
    return None


def is_allowed(env=None) -> bool:
    """True when the user has set the CLAUDE_ALLOW_CD escape hatch."""
    env = os.environ if env is None else env
    return (env.get("CLAUDE_ALLOW_CD") or "").strip().lower() in _TRUTHY


def deny_payload(target: str) -> dict:
    """The PreToolUse deny decision for a `cd` out of the worktree."""
    where = f" (`cd {target}`)" if target else " (bare `cd` -> $HOME)"
    return {
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "deny",
            "permissionDecisionReason": (
                f"Refusing a `cd` out of the project worktree{where}. The Bash cwd "
                "persists across tool calls, so this silently changes context for "
                "every later command (the footgun where a git command hits the "
                "wrong repo). Use `git -C <path> ...` for git on another repo and "
                "absolute paths for file ops - don't move the shell. Deliberate "
                "out-of-tree cd the user authorized: set CLAUDE_ALLOW_CD=1."
            ),
        }
    }


# --- I/O ---


def _log_fire(target: str, snippet: str) -> None:
    """Append one fire-log line on a deny (Issue #453 infra). Never raises."""
    payload = {
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "hook": "check_no_cd_outside_cwd",
        "action": "refused-cd-outside-worktree",
        "target": target[:120],
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
    if not cmd or "cd" not in cmd:
        return 0
    session_cwd = data.get("cwd")
    target = first_escaping_cd(cmd, session_cwd, PROJECT_ROOT, os.environ)
    if target is None:
        return 0

    _log_fire(target, cmd)
    print(json.dumps(deny_payload(target)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
