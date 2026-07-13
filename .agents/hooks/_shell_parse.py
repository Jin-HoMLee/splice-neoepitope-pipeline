"""Shared command-normalization for the `gh`-matching PostToolUse hooks.

Both `post_gh_pr_create.py` and `post_gh_pr_review_request.py` decide whether a
Bash command *invokes* a particular `gh` subcommand by shlex-tokenizing it and
looking for the subcommand at a **command start**. That walk was blind to two
shell constructs, and the blindness was silent (Issue #1130):

1. **Heredoc bodies tokenize into the stream.** `cat > b.md <<'EOF' ... EOF` puts
   every word of the body into the token list, so the token immediately before a
   following `gh` is an ordinary word (the closing delimiter) rather than a
   separator - and the command-start test fails.
2. **A newline is not a separator.** Only pure-punctuation tokens (`;`, `&&`, ...)
   reset the command-start flag, so `cmd1\ngh pr create` never matched either.

Together these mean the hooks **never fired on a heredoc-created PR** - which is
how essentially every PR with a real body is opened. Two shipped automations were
silently dead on the dominant path: the board-add + Status flip (Issue #550 /
Issue #561) and the auto-requested bot review (Issue #1073). No error, no log
line, just nothing happening.

`normalize_command()` fixes both **before** tokenization, so the existing walks
are unchanged. Crucially it preserves the anti-false-positive property those
walks were built for (PR #558: a `gh pr create` appearing *inside a quoted
argument* must NOT match):

- Heredoc **bodies are removed**, not tokenized - so prose inside one cannot
  match, which is strictly safer than the old behavior.
- Newlines are converted to separators **only outside quotes**, so a multi-line
  quoted string stays a single token and a body line beginning "gh pr create ..."
  still cannot masquerade as a command.
"""

from __future__ import annotations

import re

# Shell punctuation that shlex(punctuation_chars=True) emits as standalone
# separator tokens. Shared so the two hooks cannot drift on it.
PUNCT = set("();<>|&")

# A heredoc redirection: `<<EOF`, `<<'EOF'`, `<<"EOF"`, and the `<<-` tab-stripping
# form. The delimiter is captured so the body can be skipped up to it.
_HEREDOC_RE = re.compile(r"<<-?\s*(['\"]?)([A-Za-z_][A-Za-z0-9_]*)\1")


def strip_heredoc_bodies(cmd: str) -> str:
    """Drop heredoc bodies (and their closing delimiter lines), keep everything else.

    The line opening the heredoc is kept - it is a real command - but its body is
    removed rather than tokenized. Removing beats tokenizing here: body prose can
    never be mistaken for a command, so this only ever *reduces* false positives.
    """
    lines = cmd.splitlines()
    out: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        out.append(line)
        match = _HEREDOC_RE.search(line)
        i += 1
        if not match:
            continue
        delimiter = match.group(2)
        # Skip the body, then the closing delimiter line itself. An unterminated
        # heredoc simply consumes the rest, which is the correct read.
        while i < len(lines) and lines[i].strip() != delimiter:
            i += 1
        i += 1  # the delimiter line
    return "\n".join(out)


def newlines_to_separators(cmd: str) -> str:
    """Turn **unquoted** newlines into `;` so they act as command separators.

    Quote-aware on purpose. Replacing newlines blindly (or splitting the command
    by line) would break a multi-line quoted argument into pieces, and a body line
    that happened to begin `gh pr create ...` would then read as a real command -
    reintroducing exactly the PR #558 false positive these matchers exist to
    avoid. Inside quotes, a newline stays a newline and the string stays one token.

    Known limit (PR #1131 review, finding 4): the quote tracker is **not
    backslash-escape aware**, so an odd number of escaped quotes (`\\"`) before an
    in-quote newline can desync the flag and convert that newline to a separator.
    That is the one path here that could *weaken* rather than strengthen the
    false-positive guard. Left as-is deliberately: it needs an escaped quote, then
    a newline, then a literal `gh pr create` at the start of a line, all inside one
    quoted `--body` - and even then the hooks fail open downstream. Recorded rather
    than fixed, so the next reader knows it is a known edge and not an oversight.
    """
    out: list[str] = []
    quote: str | None = None
    for ch in cmd:
        if quote is not None:
            out.append(ch)
            if ch == quote:
                quote = None
            continue
        if ch in "'\"":
            quote = ch
            out.append(ch)
            continue
        out.append(";" if ch == "\n" else ch)
    return "".join(out)


def normalize_command(cmd: str | None) -> str:
    """Make a Bash command safe to walk for a command-start subcommand match."""
    return newlines_to_separators(strip_heredoc_bodies(cmd or ""))
