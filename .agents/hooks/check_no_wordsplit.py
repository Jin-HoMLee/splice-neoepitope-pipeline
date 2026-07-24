#!/usr/bin/env python3
"""Pre-flight hook: refuse a zsh word-split anti-pattern in an ad-hoc Bash command.

The Bash tool runs **zsh**, which - unlike bash - does NOT word-split an unquoted
scalar. So `for X in $var` and `set -- $var` silently collapse the whole value
into ONE item: the loop body runs once over the entire blob, or the positional
params get one giant element. The failure is silent (a loop that processes
nothing looks exactly like a loop that had nothing to process), and it has slipped
five times across two roles despite an Always-in-effect memory rule
(`shared/feedback_zsh_does_not_word_split.md`). Per mechanism-over-memory this is
the rung-3 escalation (Issue #1242): the rule is already inlined and still slips
at action time.

SCOPE - the ad-hoc-command surface only. This guards inline Bash-tool commands,
which `shellcheck` (CI, committed `*.sh` only, Issue #1043) structurally cannot
see - none of the five slips were in a committed script. It flags exactly two
shapes:

    for <name> in $unquoted        # unquoted scalar loop operand
    set -- $unquoted               # unquoted scalar positional reset

HIGH PRECISION is the whole design constraint (an over-eager guard is worse than
none). Everything below PASSES:

    for f in "${arr[@]}"           # quoted array expansion - the correct form
    for f in a b c                 # literal word list
    for f in $(cmd)                # command substitution operand
    for f in *.txt                 # glob
    set -- "$@"                    # quoted positional passthrough
    set -euo pipefail              # shell-option set, no `--`

HOW: the command is normalized (`_shell_parse.normalize_command`, the Issue #1142
newline/heredoc lesson), then **opaque spans are masked** - single/double quotes,
`$(...)`/`$((...))` command-and-arith substitution, and backticks all become
blanks, so any `$var` INSIDE them disappears and cannot fire. The masked command
is split on shell separators AND the loop keywords (`do`/`then`/`else`), so a
`for`/`set` is inspected only when it sits at a real command position - `echo for
f in $x` does not fire because the segment starts with `echo`, not `for`. Only a
bare **named** scalar (`$name` / `${name}`, not `$1`/`$@`/`${arr[@]}`) surviving
in the operand region trips it.

DENY, with an ESCAPE HATCH: set env `CLAUDE_ALLOW_WORDSPLIT=1` for the rare
deliberate single-iteration case (in zsh an unquoted `$var` never splits, so this
is almost always a bug, hence deny over ask - the sibling posture of
`check_no_force_push` / `check_no_emdash`).

BEST-EFFORT, NOT A GUARANTEE (GitHub Safety Wrappers doctrine, Issue #1150): a
command-string matcher has the whole-class blind spots documented there (a
word-split inside a `bash -c '...'` string is masked away and not seen, by
design - that runs bash, where the split may even be intended). The
`.agents/hook_fires.jsonl` log is the detective backstop.

Known limits, deliberately UNDER-fired to preserve precision (PR #1305 review):

- **Parameter expansions with an operator pass unflagged** - `${name:-default}`,
  `${list#pfx}`, `${files//a/b}` word-split in zsh exactly like `$name`, but the
  scalar regex requires `}` immediately after the name (so it cannot match the
  `${arr[@]}` array form), which also excludes the operator forms. Broadening to
  catch them risks re-including arrays, so this is accepted as a false-negative
  rather than a precision regression.
- **A `#` comment is not stripped** - `normalize_command` turns newlines into
  separators but leaves `#` comments in place, so a separator INSIDE a trailing
  comment (`ls  # note; for x in $tmp`) can start a spurious segment. Contrived,
  and comment-stripping is itself error-prone (a `#` inside `${x#foo}` or a quote
  is not a comment), so it is left as a documented under-strip.

Both fail toward NOT firing, matching the precision-first mandate.

Reads PreToolUse hook JSON on stdin, prints a deny decision on stdout when the
guard fires, exits 0 silently otherwise (the harness treats no-output as allow).
Fails OPEN on any parse miss / unexpected shape - a guard must never break the
user's flow on a hiccup.

See https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1242.
"""
from __future__ import annotations

import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _shell_parse  # noqa: E402

LOG_PATH = Path(__file__).resolve().parent.parent.parent / ".agents" / "hook_fires.jsonl"

_TRUTHY = {"1", "true", "yes", "on"}

# A bare NAMED scalar expansion, unquoted: `$name` or `${name}`. Deliberately does
# NOT match `$1`/`$@`/`$*` (special params, out of the named-scalar scope) nor
# `${arr[@]}` (an array subscript - the braced alt requires `}` immediately after
# the name, and the bare alt's negative lookahead rejects a following `[`).
_BARE_SCALAR_RE = re.compile(r"\$\{[A-Za-z_][A-Za-z0-9_]*\}|\$[A-Za-z_][A-Za-z0-9_]*(?![\w\[])")

# Shell separators + the loop/conditional keywords that begin a fresh command
# list. Splitting on these gives command-position anchoring for free: a segment
# starts at a real command boundary, so a `for`/`set` appearing as an *argument*
# (e.g. `echo for f in $x`) lands mid-segment and its anchor fails to match.
# Braces are deliberately NOT split on: a brace GROUP (`{ cmd; }`) is rare in an
# ad-hoc command, and splitting on `{`/`}` would shred a `${var}` expansion (the
# very operand we need to see). A leading brace-group intro is handled by the
# anchors' optional `{ ` prefix instead.
_SEGMENT_SPLIT_RE = re.compile(r"&&|\|\||[;&|()]|\bdo\b|\bthen\b|\belse\b|\bwhile\b|\buntil\b")

# A for-loop header at the start of a segment: `for NAME in <operands>`. An
# optional leading `{ ` absorbs a brace-group intro (`{ for x in $y; }`).
_FOR_RE = re.compile(r"^\s*(?:\{\s+)?for\s+[A-Za-z_][A-Za-z0-9_]*\s+in\b(.*)$", re.DOTALL)
# A positional-param reset at the start of a segment: `set -- <operands>`.
_SET_RE = re.compile(r"^\s*(?:\{\s+)?set\s+--\s+(.*)$", re.DOTALL)


# --- pure helpers (unit-tested, no I/O) ---


def mask_opaque(cmd: str) -> str:
    """Blank every opaque span so a `$var` inside one cannot be seen.

    Opaque = single-quoted, double-quoted, backtick, and `$(...)`/`$((...))`
    substitution spans. Each is replaced char-for-char with spaces (length
    preserved, so nothing shifts). This is what makes the guard high-precision:
    the correct forms carry their variable *inside* quotes (`"${arr[@]}"`,
    `"$@"`) or a substitution (`$(cmd)`), so masking removes them and only a
    genuinely unquoted operand survives to be matched.

    Single left-to-right pass; on entering a construct it blanks through to the
    construct's close without recursing (an inner quote inside `$(...)` is simply
    part of the blanked region). Backslash-escaping is honored inside double
    quotes only (single quotes do not process escapes in shell).
    """
    out = list(cmd)
    n = len(cmd)
    i = 0
    while i < n:
        ch = cmd[i]
        if ch == "'":
            j = i + 1
            while j < n and cmd[j] != "'":
                j += 1
            for k in range(i, min(j + 1, n)):
                out[k] = " "
            i = j + 1
        elif ch == '"':
            j = i + 1
            while j < n and cmd[j] != '"':
                if cmd[j] == "\\":
                    j += 2
                    continue
                j += 1
            for k in range(i, min(j + 1, n)):
                out[k] = " "
            i = j + 1
        elif ch == "`":
            j = i + 1
            while j < n and cmd[j] != "`":
                j += 1
            for k in range(i, min(j + 1, n)):
                out[k] = " "
            i = j + 1
        elif ch == "$" and i + 1 < n and cmd[i + 1] == "(":
            # `$(...)` or `$((...))` - balance parens from the opening `$(`.
            depth = 0
            j = i + 1
            while j < n:
                if cmd[j] == "(":
                    depth += 1
                elif cmd[j] == ")":
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            for k in range(i, min(j + 1, n)):
                out[k] = " "
            i = j + 1
        elif ch == "\\" and i + 1 < n:
            # A backslash-escaped char outside any quote/substitution is a literal:
            # `\$var` is the text "$var", not an expansion. Blank the escape pair so
            # the escaped `$` cannot trip _BARE_SCALAR_RE - an over-fire is the one
            # direction this high-precision guard must avoid (PR #1305 review).
            out[i] = " "
            out[i + 1] = " "
            i += 2
        else:
            i += 1
    return "".join(out)


def _operand_has_bare_scalar(operand: str) -> bool:
    """True iff a (already-masked) operand region contains an unquoted named scalar."""
    return _BARE_SCALAR_RE.search(operand) is not None


def segment_wordsplits(segment: str) -> bool:
    """True iff a single masked segment is a `for … in $var` / `set -- $var` shape."""
    m = _FOR_RE.match(segment)
    if m and _operand_has_bare_scalar(m.group(1)):
        return True
    m = _SET_RE.match(segment)
    if m and _operand_has_bare_scalar(m.group(1)):
        return True
    return False


def command_wordsplits(cmd: str) -> bool:
    """True iff any command-position `for`/`set --` uses an unquoted scalar operand.

    Normalizes (heredoc/newline), masks opaque spans, splits into command-position
    segments, and checks each. False on anything else (fail-open on odd shapes is
    the caller's job; this only ever returns False when it sees no anti-pattern).
    """
    if not cmd:
        return False
    masked = mask_opaque(_shell_parse.normalize_command(cmd))
    for segment in _SEGMENT_SPLIT_RE.split(masked):
        if segment_wordsplits(segment):
            return True
    return False


def is_allowed(env=None) -> bool:
    """True when the user has set the CLAUDE_ALLOW_WORDSPLIT escape hatch."""
    env = os.environ if env is None else env
    return (env.get("CLAUDE_ALLOW_WORDSPLIT") or "").strip().lower() in _TRUTHY


def deny_payload() -> dict:
    """The PreToolUse deny decision for a zsh word-split anti-pattern."""
    return {
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "deny",
            "permissionDecisionReason": (
                "Refusing an unquoted-scalar `for … in $var` / `set -- $var`: the Bash "
                "tool runs zsh, which does NOT word-split an unquoted scalar, so this "
                "iterates ONCE over the whole value (or sets one giant positional param) "
                "- silently, with no error. Use an array for a static list "
                "(`arr=(a b c); for x in \"${arr[@]}\"`), `printf '%s\\n' \"$blob\" | "
                "while read -r x` for a streamed blob, or parse structured data with "
                "jq/python. If you genuinely intend a single-iteration pass, set "
                "CLAUDE_ALLOW_WORDSPLIT=1 for this session."
            ),
        }
    }


# --- I/O ---


def _log_fire(snippet: str) -> None:
    """Append one fire-log line on a deny (Issue #453 infra). Never raises."""
    payload = {
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "hook": "check_no_wordsplit",
        "action": "refused-zsh-wordsplit",
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
        return 0  # only guards Bash commands

    if is_allowed():
        return 0

    cmd = (data.get("tool_input") or {}).get("command", "")
    try:
        fires = command_wordsplits(cmd)
    except Exception:
        return 0  # fail open on any unexpected shape
    if not fires:
        return 0

    _log_fire(cmd)
    print(json.dumps(deny_payload()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
