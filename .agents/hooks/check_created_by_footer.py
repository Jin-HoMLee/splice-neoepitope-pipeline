#!/usr/bin/env python3
"""Pre-flight hook: refuse a `gh (pr|issue) comment` whose authored body lacks the
`**Created by:** <Role>` footer (Issue #1094).

The footer rule ("every issue/PR comment you author carries a `**Created by:**
<Role>` line") slipped >= 3x on the same shape in one session, twice AFTER it had
been inlined into Always-in-effect memory - a rung-3 mechanism escalation (memory
tried and demonstrably does not hold). The footer is the last, feedback-less line
of a long body: nothing breaks, nothing warns, so a `PreToolUse` hook fixes it and
memory does not.

Design (mirrors `check_at_claude.py`; composes with the #1197 accept-the-hole
decision):
- **DETECT** the comment command on the NORMALIZED text, so a `gh comment` after a
  heredoc is still seen (the #1130/#1142 lesson - a guard blind to the heredoc
  shape is born dead on the path bodies are actually written).
- **INSPECT** the body, but only DENY when the body is genuinely inspectable: an
  inline `--body`, or a heredoc feeding `--body-file -` (the body is in the raw
  command string). A real `--body-file <path>` is NOT inspectable, so it fails
  OPEN (allow) - denying a compliant file-bodied comment would be a false positive.
- **EXEMPTIONS:** a `**From:** <Role>` opener (coordination messages already carry
  attribution), and the bare bot-review trigger (`@claude review` / `@-claude
  review`), which is a trigger string, not an authored comment.

Fails OPEN on every uncertain path (unparseable command, no visible body, real
`--body-file`, `gh api` forms). Escape hatch: `CLAUDE_ALLOW_NO_FOOTER=1`. Each deny
appends one line to `.agents/hook_fires.jsonl` (Issue #453 fire-log infra).

See https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1094.
"""
from __future__ import annotations

import json
import os
import re
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _shell_parse  # noqa: E402

LOG_PATH = Path(__file__).resolve().parent.parent.parent / ".agents" / "hook_fires.jsonl"

_COMMENT_RE = re.compile(r"(^|[;&|]\s*)gh\s+(pr|issue)\s+comment\b")
_TRIGGER_RE = re.compile(r"^@-?claude review$")
_FOOTER = "**Created by:**"
_FROM_OPENER = "**From:**"


# --- pure helpers (unit-tested, no I/O) ---


def matches_comment(cmd: str) -> bool:
    """True iff `cmd` invokes a `gh (pr|issue) comment` at a command start.

    Detected on the NORMALIZED command so a comment following a heredoc (or on a
    later line) still matches - the shape a long body is actually written with.
    """
    return bool(_COMMENT_RE.search(_shell_parse.normalize_command(cmd)))


def extract_heredoc_bodies(cmd: str) -> list[str]:
    """Every heredoc body in `cmd`, in order (inverse of `strip_heredoc_bodies`).

    Collects the lines between each `<<DELIM` opener and its closing `DELIM`. The
    caller uses the COUNT: a `--body-file -` comment binds to the heredoc only when
    there is exactly one (unambiguous); a second, earlier heredoc (e.g. a scratch
    `cat > f <<EOF` before the `gh comment`) makes the binding ambiguous, so the
    caller fails open rather than read the wrong body (PR #1229 finding 1).
    """
    lines = cmd.splitlines()
    bodies: list[str] = []
    i = 0
    while i < len(lines):
        m = _shell_parse.HEREDOC_RE.search(lines[i])
        i += 1
        if not m:
            continue
        delimiter = m.group(2)
        body: list[str] = []
        while i < len(lines) and lines[i].strip() != delimiter:
            body.append(lines[i])
            i += 1
        bodies.append("\n".join(body))
        i += 1  # skip the closing delimiter line
    return bodies


def extract_body(cmd: str) -> str | None:
    """The comment body if it is INSPECTABLE, else None (caller fails open).

    Inspectable: an inline `--body`/`-b` value, or a heredoc feeding `--body-file -`
    / `-F -` **when it is the only heredoc** (see below). Not inspectable (-> None):
    a real `--body-file <path>` (content is in a file we do not read, per the #1197
    accept-the-hole posture), no body flag, an untokenizable command, or an
    AMBIGUOUS `--body-file -` where more than one heredoc is present - the comment
    body is the heredoc on the `gh` line, but a chained command may open an earlier
    scratch heredoc, so binding to "the first" would read the wrong body and
    false-positive-deny a compliant comment (PR #1229 finding 1). Ambiguity ->
    fail open.
    """
    heredocs = extract_heredoc_bodies(cmd)
    try:
        tokens = shlex.split(_shell_parse.strip_heredoc_bodies(cmd))
    except ValueError:
        return None  # fail open
    body_file = None
    for i, tok in enumerate(tokens):
        if tok in ("--body", "-b") and i + 1 < len(tokens):
            return tokens[i + 1]
        if tok.startswith("--body="):
            return tok.split("=", 1)[1]
        if tok in ("--body-file", "-F") and i + 1 < len(tokens):
            body_file = tokens[i + 1]
        elif tok.startswith("--body-file="):
            body_file = tok.split("=", 1)[1]
    if body_file == "-":
        # Exactly one heredoc -> unambiguously the comment body. Zero (stdin we
        # cannot see) or more than one (ambiguous) -> fail open.
        return heredocs[0] if len(heredocs) == 1 else None
    return None  # real --body-file <path>, or no body flag


def is_exempt(body: str) -> bool:
    """A `**From:**`-opener coordination message, or the bare bot-review trigger."""
    stripped = body.strip()
    if stripped.startswith(_FROM_OPENER):
        return True
    return bool(_TRIGGER_RE.match(stripped))


def has_footer(body: str) -> bool:
    return _FOOTER in body


def deny_payload() -> dict:
    return {
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "deny",
            "permissionDecisionReason": (
                "This `gh (pr|issue) comment` body is missing its "
                "`**Created by:** <Role>` footer. Every authored comment carries it "
                "so a comment's role is unambiguous when multiple personas are "
                "active (rule: shared/feedback_github_workflow.md). Add "
                "`**Created by:** <YourRole>` as the last line. Exempt: a comment "
                "that opens with `**From:** <Role>` (coordination messages already "
                "carry attribution), and the bare `@-claude review` trigger. For a "
                "deliberate no-footer comment, set CLAUDE_ALLOW_NO_FOOTER=1."
            ),
        }
    }


def _log_fire(cmd: str) -> None:
    """Append one fire-log line on a deny (Issue #453 infra). Never raises."""
    payload = {
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "hook": "check_created_by_footer",
        "action": "denied-missing-created-by-footer",
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
    if os.environ.get("CLAUDE_ALLOW_NO_FOOTER") == "1":
        return 0  # explicit escape hatch

    try:
        data = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0  # fail open

    cmd = (data.get("tool_input") or {}).get("command", "")
    if not matches_comment(cmd):
        return 0  # not a `gh (pr|issue) comment`

    body = extract_body(cmd)
    if body is None:
        return 0  # not inspectable -> fail open

    if is_exempt(body) or has_footer(body):
        return 0

    _log_fire(cmd)
    print(json.dumps(deny_payload()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
