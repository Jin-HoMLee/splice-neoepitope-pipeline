#!/usr/bin/env python3
"""Pre-flight hook: refuse `gh (pr|issue) (comment|create)` whose body contains
a literal @claude bot mention, except the canonical `@claude review` trigger.

See memory/shared/feedback_no_at_claude_mention.md for the originating rule.
Spec: https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/360

Reads PreToolUse hook JSON on stdin, prints a deny decision on stdout when the
guard fires, exits 0 silently otherwise (the harness treats no-output as allow).
"""
import json
import re
import sys


def main() -> int:
    try:
        data = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0

    cmd = data.get("tool_input", {}).get("command", "")

    if not re.search(r"(^|[;&|]\s*)gh\s+(pr|issue)\s+(comment|create|edit)\b", cmd):
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
                        "Stray @claude mention in `gh (pr|issue) (comment|create)` body. "
                        "The Claude Code GitHub Action triggers on ANY literal @claude "
                        "(incl. inside parens, ACs, code spans, quoted historical refs). "
                        "See memory/shared/feedback_no_at_claude_mention.md. For "
                        "non-trigger references, use the zero-width workaround @-claude "
                        "(literal substring @claude must not appear)."
                    ),
                }
            }
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
