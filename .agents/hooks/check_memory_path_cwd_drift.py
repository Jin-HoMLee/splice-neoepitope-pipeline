#!/usr/bin/env python3
"""Pre-flight hook: refuse a memory-path operation while the cwd has drifted.

House rule "cwd IS the role's own clone root": the main-session Bash cwd persists
across tool calls, so a `cd` into a scratch/temp dir that is never undone leaves
the shell standing outside the clone. A *later* command that uses a clone-root-
relative memory path (`grep .agents/memory/...`, an Edit of `.agents/memory/x`)
then resolves against the wrong directory - silent wrong-file memory I/O.

This is the narrow residual gap NOT covered by `check_no_cd_outside_cwd.py`,
which deliberately ALLOWS `cd` into temp roots (the scratchpad) - exactly the door
the drift walks through. This guard catches the one high-harm moment: a relative
memory path used while the persisted cwd is outside the clone.

SCOPE (deliberately narrow, so a false-positive deny is unlikely): fires ONLY when
BOTH hold -
  1. the session `cwd` (from the PreToolUse JSON) resolves OUTSIDE the clone
     subtree (not the root, not a subdir), AND
  2. the tool references a RELATIVE path under `.agents/memory/` or
     `.claude/memory/` (Bash: any command token; Edit/Write/MultiEdit: file_path).
An ABSOLUTE memory path is always safe (it does not depend on cwd) -> allowed.
A cwd inside the clone (root or subdir) -> allowed. Non-memory paths -> allowed.

WHY A HOOK, NOT `permissions.deny`: a static rule cannot see the runtime cwd or
resolve absolute-vs-relative. Sibling of `check_no_cd_outside_cwd.py`.

BEST-PRACTICE FRAMING: persistent-shell agent harnesses accept cwd drift as a
known tradeoff and mitigate it by preferring absolute paths / `git -C` and by
making drift LOUD rather than silent. This guard makes the single high-harm drift
moment loud, and its deny message nudges to the real fix.

FAILS OPEN (allows) on: an unresolvable / missing cwd, a dynamic path (`$VAR`,
backticks, command substitution), a tokenization failure, or an unexpected JSON
shape. A guard must never break a legitimate command on ambiguity.

ESCAPE HATCH: set env `CLAUDE_ALLOW_CWD_DRIFT=1` for a deliberate memory op from a
drifted cwd the user authorized (prefer an absolute path / `cd` back instead).

Reads PreToolUse hook JSON on stdin, prints a deny decision on stdout when the
guard fires, exits 0 silently otherwise. Fires are appended to
`.agents/hook_fires.jsonl` (Issue #453 infra).

See https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1053.
"""
from __future__ import annotations

import json
import os
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = str(Path(__file__).resolve().parent.parent.parent)
LOG_PATH = Path(PROJECT_ROOT) / ".agents" / "hook_fires.jsonl"

_TRUTHY = {"1", "true", "yes", "on"}
# The clone-root-relative memory directories a drifted cwd would mis-resolve.
_MEMORY_SEGMENTS = (".agents/memory", ".claude/memory")
_PATH_TOOLS = {"Edit", "MultiEdit", "Write"}


# --- pure helpers (unit-tested, no I/O) ---


def cwd_outside_clone(cwd, project_root: str) -> bool:
    """True iff `cwd` resolves OUTSIDE the clone subtree (not root, not a subdir).

    A missing / empty / unresolvable cwd returns False (cannot tell -> the caller
    fails open and does not fire). Both paths are realpath-normalized so a
    symlinked clone component does not spuriously read as "outside".
    """
    if not cwd:
        return False
    try:
        c = os.path.realpath(cwd)
        r = os.path.realpath(project_root)
    except (OSError, ValueError):
        return False
    if c == r or c.startswith(r + os.sep):
        return False
    return True


def is_relative_memory_token(token: str) -> bool:
    """True iff `token` is a RELATIVE path into a memory dir (so it resolves
    against cwd). Absolute paths and dynamic tokens are not relative memory refs.
    """
    if not token or "$" in token or "`" in token:
        return False  # dynamic -> unresolvable -> fail open
    if token.startswith("/") or token.startswith("~"):
        return False  # absolute / home-anchored -> cwd-independent -> safe
    # A leading `./` is the memory dir itself; a leading `../` is NOT matched -
    # a parent-relative escape (`../.agents/memory`) is a deliberate fail-open gap
    # (not the recorded incident shape; erring toward allow on the deny path).
    stripped = token[2:] if token.startswith("./") else token
    # Anchor at a path boundary so a sibling like `.agents/memory-archive/` does
    # NOT over-match - a false-positive deny is the costly direction here.
    return any(stripped == seg or stripped.startswith(seg + "/") for seg in _MEMORY_SEGMENTS)


def command_has_relative_memory_ref(cmd: str):
    """For a Bash command, return the first relative-memory-path token, or None.

    Returns None on a tokenization failure (caller fails open). Quoted strings
    survive as single tokens.
    """
    try:
        tokens = shlex.split(cmd or "", posix=True)
    except ValueError:
        return None
    for tok in tokens:
        if is_relative_memory_token(tok):
            return tok
    return None


def offending_ref(tool_name: str, tool_input: dict):
    """Return the relative-memory-path reference that would mis-resolve, or None.

    Bash -> scan the command tokens. Edit/Write/MultiEdit -> the `file_path`.
    """
    tool_input = tool_input or {}
    if tool_name == "Bash":
        return command_has_relative_memory_ref(tool_input.get("command", ""))
    if tool_name in _PATH_TOOLS:
        fp = tool_input.get("file_path", "")
        return fp if is_relative_memory_token(fp) else None
    return None


def is_allowed(env=None) -> bool:
    """True when the user has set the CLAUDE_ALLOW_CWD_DRIFT escape hatch."""
    env = os.environ if env is None else env
    return (env.get("CLAUDE_ALLOW_CWD_DRIFT") or "").strip().lower() in _TRUTHY


def deny_payload(ref: str, cwd: str, project_root: str) -> dict:
    """The PreToolUse deny decision for a drifted-cwd memory operation."""
    return {
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "deny",
            "permissionDecisionReason": (
                f"Refusing a relative memory path (`{ref}`) while the shell cwd has "
                f"drifted OUT of the clone root.\n  cwd:   {cwd}\n  clone: {project_root}\n"
                "A clone-root-relative memory path resolves against the wrong "
                "directory from here (silent wrong-file memory I/O). Fix: `cd` back "
                "to the clone root, or use an absolute path (the cwd persists across "
                "tool calls - don't rely on where the shell is standing). Deliberate "
                "drifted memory op the user authorized: set CLAUDE_ALLOW_CWD_DRIFT=1."
            ),
        }
    }


# --- I/O ---


def _log_fire(ref: str, cwd: str, snippet: str) -> None:
    """Append one fire-log line on a deny (Issue #453 infra). Never raises."""
    payload = {
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "hook": "check_memory_path_cwd_drift",
        "action": "refused-memory-path-cwd-drift",
        "ref": ref[:120],
        "cwd": (cwd or "")[:200],
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

    tool_name = data.get("tool_name")
    if tool_name != "Bash" and tool_name not in _PATH_TOOLS:
        return 0

    if is_allowed():
        return 0

    cwd = data.get("cwd")
    if not cwd_outside_clone(cwd, PROJECT_ROOT):
        return 0  # cwd is fine (or unknown) -> nothing to guard

    tool_input = data.get("tool_input") or {}
    ref = offending_ref(tool_name, tool_input)
    if ref is None:
        return 0

    snippet = tool_input.get("command", "") or tool_input.get("file_path", "")
    _log_fire(ref, cwd, snippet)
    print(json.dumps(deny_payload(ref, cwd, PROJECT_ROOT)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
