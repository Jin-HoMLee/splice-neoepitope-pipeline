#!/usr/bin/env python3
"""Pre-flight hook: refuse an Edit/MultiEdit/Write that ADDS an em-dash (U+2014)
or en-dash (U+2013) to Claude's newly-written content.

Jin-Ho's global rule is "never use the em dash, use a plain ASCII hyphen '-'
instead" (~/.claude/CLAUDE.md). It slipped twice in Claude's own additions inside
two sessions (PR #916 em+en, PR #917 em) despite living in memory - the
documented mechanism-over-memory trigger (>= 2x on the same shape). This is the
rung-3 escalation, a sibling of `.agents/hooks/check_at_claude.py` and
`check_gh_issue_develop_parent.py`.

CHARACTER SCOPE DECISION: em (U+2014) AND en (U+2013). Both have slipped and
neither is a plain dash; the en-dash slipped as `S1-S7`. Horizontal ellipsis
(U+2026) and other smart punctuation are deliberately OUT of scope.

DELTA SCOPING: the guard flags only chars it is *net-adding*. For Edit/MultiEdit
it compares each `new_string` against its `old_string`; for Write it compares the
new `content` against the file's current on-disk content (empty for a new file).
So editing a file that legitimately already contains em-dashes (most historical
memory files) never false-positives unless the edit introduces a NEW one. Note
Write is COARSER than Edit/MultiEdit: it compares whole-file counts, so a single
Write that both removes a pre-existing dash and adds a new one nets to zero and
is not flagged (Edit/MultiEdit are span-local and exact). NotebookEdit is
deliberately NOT in GUARDED_TOOLS - notebook cells are out of scope for now.

ESCAPE HATCH: set env `CLAUDE_ALLOW_EMDASH=1` for a genuine verbatim need, or edit
an allowlisted path (the guard's own source/tests, where the chars appear as
fixtures). See `_ALLOWLIST_SUBSTRINGS`.

Reads PreToolUse hook JSON on stdin, prints a deny decision on stdout when the
guard fires, exits 0 silently otherwise (the harness treats no-output as allow).
Fails OPEN on any parse miss / unreadable file / unexpected shape - a guard must
never break legitimate writing on a hiccup. Fires are appended to
`.agents/hook_fires.jsonl` (Issue #453 infra).

See https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/920.
"""
from __future__ import annotations

import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

LOG_PATH = Path(__file__).resolve().parent.parent.parent / ".agents" / "hook_fires.jsonl"

# Named char scope: em-dash + en-dash. Defined via \u escapes so THIS source file
# stays ASCII-clean (and is not itself a fixture the guard would flag).
TARGET_CHARS = {
    "\u2014": "em-dash (U+2014)",
    "\u2013": "en-dash (U+2013)",
}

GUARDED_TOOLS = ("Edit", "MultiEdit", "Write")

# Paths where these chars are legitimate (the guard's own source + test fixtures,
# and the fire-log which may record offending content). Substring match.
_ALLOWLIST_SUBSTRINGS = ("check_no_emdash", ".agents/hook_fires.jsonl")

_TRUTHY = {"1", "true", "yes", "on"}


# --- pure helpers (unit-tested, no I/O) ---


def added_violations(old: str, new: str) -> list[str]:
    """Names of target chars whose count is strictly higher in `new` than `old`.

    Net-added only: a char present in equal measure in both (preserved) or removed
    is not a violation. Empty list == clean.
    """
    out = []
    for ch, name in TARGET_CHARS.items():
        if new.count(ch) > old.count(ch):
            out.append(name)
    return out


def _read_file(path: str) -> str | None:
    """Current on-disk content for a Write target. "" for a new file; None on any
    other read error (caller then fails open)."""
    if not path:
        return ""
    try:
        return Path(path).read_text(encoding="utf-8")
    except FileNotFoundError:
        return ""  # brand-new file: all content is "added"
    except (OSError, UnicodeDecodeError):
        return None  # can't compare reliably -> fail open


def guarded_pairs(tool_name, tool_input, read_file=_read_file):
    """(old, new) string pairs to delta-scan for `tool_name`, or None to fail open.

    Edit -> one (old_string, new_string) pair.
    MultiEdit -> one pair per edit.
    Write -> (current_file_content, content); None if the file is unreadable.
    Any other tool -> None (not guarded).
    """
    if tool_name == "Edit":
        return [(tool_input.get("old_string", "") or "",
                 tool_input.get("new_string", "") or "")]
    if tool_name == "MultiEdit":
        edits = tool_input.get("edits") or []
        pairs = []
        for e in edits:
            if not isinstance(e, dict):
                return None  # unexpected shape -> fail open
            pairs.append((e.get("old_string", "") or "", e.get("new_string", "") or ""))
        return pairs
    if tool_name == "Write":
        old = read_file(tool_input.get("file_path", "") or "")
        if old is None:
            return None  # unreadable existing file -> fail open
        return [(old, tool_input.get("content", "") or "")]
    return None


def is_allowlisted(file_path, env=None) -> bool:
    """True when the write is exempt: env escape hatch set, or an allowlisted path."""
    env = os.environ if env is None else env
    val = (env.get("CLAUDE_ALLOW_EMDASH") or "").strip().lower()
    if val in _TRUTHY:
        return True
    return any(sub in (file_path or "") for sub in _ALLOWLIST_SUBSTRINGS)


def deny_payload(violations: list[str], file_path: str) -> dict:
    """The PreToolUse deny decision for a net-added em/en dash."""
    chars = ", ".join(violations)
    where = f" in {file_path}" if file_path else ""
    return {
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "deny",
            "permissionDecisionReason": (
                f"This edit adds {chars}{where}. House rule: use a plain ASCII "
                "hyphen '-', never an em/en dash (~/.claude/CLAUDE.md). Rewrite the "
                "added text with '-'. The guard scans only the delta you add, so "
                "pre-existing dashes are never flagged. Genuine verbatim need: set "
                "CLAUDE_ALLOW_EMDASH=1, or write to an allowlisted path."
            ),
        }
    }


# --- I/O ---


def _log_fire(tool_name: str, file_path: str, violations: list[str]) -> None:
    """Append one fire-log line on a deny (Issue #453 infra). Never raises."""
    payload = {
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "hook": "check_no_emdash",
        "tool": tool_name,
        "file": file_path,
        "action": "refused-emdash:" + "+".join(violations),
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

    tool_name = data.get("tool_name", "")
    if tool_name not in GUARDED_TOOLS:
        return 0

    tool_input = data.get("tool_input") or {}
    if not isinstance(tool_input, dict):
        return 0  # fail open

    file_path = tool_input.get("file_path", "") or ""
    if is_allowlisted(file_path):
        return 0

    pairs = guarded_pairs(tool_name, tool_input)
    if pairs is None:
        return 0  # fail open

    violations: list[str] = []
    for old, new in pairs:
        for v in added_violations(old, new):
            if v not in violations:
                violations.append(v)
    if not violations:
        return 0

    _log_fire(tool_name, file_path, violations)
    print(json.dumps(deny_payload(violations, file_path)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
