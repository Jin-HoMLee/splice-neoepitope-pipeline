#!/usr/bin/env python3
"""Pre-flight hook: refuse `gh (pr|issue) (comment|create|edit)` whose body
contains a literal @claude bot mention, except the canonical `@claude review`
trigger.

See memory/shared/feedback_no_at_claude_mention.md for the originating rule.
Spec: https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/360

Reads PreToolUse hook JSON on stdin, prints a deny decision on stdout when the
guard fires, exits 0 silently otherwise (the harness treats no-output as allow).
"""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _shell_parse  # noqa: E402

_GH_COMMENT_RE = re.compile(r"(^|[;&|]\s*)gh\s+(pr|issue)\s+(comment|create|edit)\b")


def main() -> int:
    try:
        data = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0

    cmd = data.get("tool_input", {}).get("command", "")

    # Detect the command on the NORMALIZED text (Issue #1130): the anchor only
    # accepted a command at string-start or after `;&|`, so a `gh pr comment`
    # following a heredoc - the dominant way a long body is written - never
    # matched, and this guard was silently dead on that path. Normalizing turns an
    # unquoted newline into a separator so the anchor sees the command.
    #
    # Scan for the mention on the RAW command, NOT the normalized one: normalizing
    # strips heredoc bodies, and for this guard the heredoc body is exactly the
    # text we need to inspect. Detect on normalized, inspect on raw.
    if not _GH_COMMENT_RE.search(_shell_parse.normalize_command(cmd)):
        return 0

    if "@claude" not in cmd:
        return 0

    total = cmd.count("@claude")
    canonical = len(re.findall(r'--body\s+"@claude review"', cmd)) + len(
        re.findall(r"--body\s+'@claude review'", cmd)
    )
    if total == canonical:
        return 0

    print(
        json.dumps(
            {
                "hookSpecificOutput": {
                    "hookEventName": "PreToolUse",
                    "permissionDecision": "deny",
                    "permissionDecisionReason": (
                        "Stray @claude mention in `gh (pr|issue) (comment|create|edit)` body. "
                        "The Claude Code GitHub Action triggers on ANY literal @claude "
                        "(incl. inside parens, ACs, code spans, quoted historical refs). "
                        "See memory/shared/feedback_no_at_claude_mention.md. "
                        "To intentionally request a bot review, use exactly "
                        '--body "@claude review" (the canonical trigger is the one '
                        "allowed form). For any OTHER (non-trigger) reference, use the "
                        "zero-width workaround @-claude (literal substring @claude must "
                        "not appear)."
                    ),
                }
            }
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
