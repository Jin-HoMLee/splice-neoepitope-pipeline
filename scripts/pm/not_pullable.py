"""Detect body-only "not pullable" gates on a board Issue.

Issue #1248. The Ready-queue floor guard (`scripts/check_ready_queue.sh`) counts
Ready *depth* - "has concrete ACs = pullable" - and is structurally blind to
gates that exist only in an Issue's prose. So a decision-blocked or
trigger-gated queue reads as healthy, and floor-pressure Replenishment
re-commits Issues that look Ready but cannot actually be worked. Issue #841 was
swept into Ready twice on exactly this shape; Issue #929 sat in Ready with two
unsettled decision-fork ACs.

This is a pure function of the body, like the sibling tells it mirrors: no
stored state to drift, any run re-derives it.

Two tells, because the two fixtures fail differently:

  TRIGGER GATE   a prose marker anywhere in the body. Issue #841's "Build only
                 if it slips a 2nd time" lives in `## Notes`, so an AC-scoped
                 scan would miss it entirely.

  NEEDS DESIGN   a decision verb leading an UNTICKED acceptance criterion.
                 Issue #929's "Decide scope:" / "Decide home:". This one IS
                 AC-scoped, because a body that merely narrates a past decision
                 ("we decided to converge on STAR") is not an open fork.

FALSE-POSITIVE DISCIPLINE is the load-bearing half. An Issue that *documents*
the markers must not trip - Issue #1248 itself lists them while proposing this
scan. So the scan is blind to text quoted as code (fenced blocks and inline
spans), which is the same "using a directive vs talking about it" distinction
that Issue #1126 had to draw for the `skip-lab-notebook` marker. An over-eager
guard is worse than none: it gets disabled.

Deliberately NOT a denylist widen. The marker set is small, specific, and
matches phrasings we actually use; growing it on every near-miss is the
non-convergent path the Issue #1150 posture warns about.

KNOWN TENSION, surfaced by the first live run and deliberately NOT resolved here
(Issue #1248 comment thread owns it). The DoR in shared/feedback_board_hygiene.md
says an Issue whose ACs contain an open design decision is not Ready "no matter
how small it looks", and read faithfully that also catches a *research* or *eval*
Issue where deciding IS the deliverable - Issue #659's "Determine whether
ImmunoStruct accepts wild-type-free inputs" is the research question, not a
blocker on starting. This module implements the convention as written rather
than silently narrowing it, because whether research Issues want a carve-out is
a governance call, not an implementation detail. If the answer is yes, the fix
belongs in the DoR first and here second.

KNOWN MISSES, by design. A gate that is not in the body is invisible here:
Issue #876's trigger lives in a conversational verdict and a role's post-it, and
a capability prerequisite phrased "gated on sequence access" (Issue #817) is the
separate capability-prerequisite axis. Recall is deliberately traded for
precision, because this is a linter and a false ALARM gets a guard disabled.
"""
import re
from typing import Optional

# Prose markers that gate an Issue on an external trigger or a date. Matched
# case-insensitively, anywhere in the body outside code. Each one is a phrase
# we actually use, not a guess.
TRIGGER_MARKERS = (
    r"build only if",
    r"build only on",
    r"do not commit before",
    r"trigger-gated",
    r"revival trigger",
    r"revisit when",
)

# A gate is a DIRECTIVE, and a directive opens a clause. Anchoring to a clause
# start is what separates invoking a marker from talking about one, and it is
# the distinction code-stripping alone cannot draw: Issue #1248's prose says
# "chartered to build only on a 2nd trigger event" and "a decision-blocked or
# trigger-gated queue", both mid-clause, while Issue #841's real gate opens a
# sentence ("... incident #1. Build only if it slips a 2nd time"). Same words,
# opposite intent, and position is the signal.
_TRIGGER_RE = re.compile(r"^(?:{})".format("|".join(TRIGGER_MARKERS)), re.IGNORECASE)

# Leading markdown noise on a clause: list bullets, emphasis, blockquote marks.
_CLAUSE_LEAD_RE = re.compile(r"^[\s>]*(?:[-*+]\s+)?[\s*_]*")

# Clause boundaries inside a line: sentence enders plus the colon that opens a
# labelled directive ("Revival trigger: revisit when ...").
_CLAUSE_SPLIT_RE = re.compile(r"(?<=[.!?:])\s+")

# Verbs that, leading an unticked acceptance criterion, mean the scope is a
# decision rather than a deliverable.
_DECISION_AC_RE = re.compile(
    r"^\s*[-*]\s*\[ \]\s*(decide|determine|choose|settle|pick)\b",
    re.IGNORECASE,
)

_FENCED_CODE_RE = re.compile(r"```.*?```", re.DOTALL)
_INLINE_CODE_RE = re.compile(r"`[^`\n]*`")
_AC_HEADING_RE = re.compile(r"^#{1,6}\s*acceptance criteria\s*$", re.IGNORECASE)
_ANY_HEADING_RE = re.compile(r"^#{1,6}\s+\S")


def strip_code(text):
    """Remove fenced blocks and inline spans.

    Quoting a marker as code is how an Issue talks *about* the convention
    instead of invoking it. Stripping rather than flagging keeps the surrounding
    prose scannable, so a body can both document a marker and carry a real one.
    """
    return _INLINE_CODE_RE.sub(" ", _FENCED_CODE_RE.sub(" ", text))


def iter_clauses(text):
    """Yield each clause with its leading markdown noise removed.

    A clause is a line, or a sentence within a line. Markers are tested against
    a clause START, so a marker quoted mid-sentence while describing some other
    Issue's gate does not read as this Issue's own gate.
    """
    for line in text.splitlines():
        for clause in _CLAUSE_SPLIT_RE.split(line):
            yield _CLAUSE_LEAD_RE.sub("", clause)


def acceptance_criteria_lines(body):
    """Return the lines under the canonical `## Acceptance criteria` heading.

    Scoped to the one canonical heading on purpose. Broadening to Plan / Tasks /
    Checklist would sweep in the non-gating boxes those sections routinely carry
    - the same reasoning that keeps `audit_and_merge.sh` gate 2 keyed to a single
    heading, with a separate lint for stray boxes.
    """
    lines = body.splitlines()
    out = []
    in_ac = False
    for line in lines:
        if _AC_HEADING_RE.match(line):
            in_ac = True
            continue
        if in_ac and _ANY_HEADING_RE.match(line):
            break
        if in_ac:
            out.append(line)
    return out


def scan_not_pullable(body):
    # type: (Optional[str]) -> Optional[str]
    """Return a short reason an Issue is not pullable, or None if it is.

    Fails *open* (returns None) on empty or missing input: this feeds a board
    read, and a scan hiccup must never be the thing that makes a healthy queue
    look broken.
    """
    if not body:
        return None

    scannable = strip_code(body)

    for clause in iter_clauses(scannable):
        match = _TRIGGER_RE.match(clause)
        if match:
            return "trigger-gated: body says {!r}".format(match.group(0).lower())

    for line in acceptance_criteria_lines(scannable):
        verb = _DECISION_AC_RE.match(line)
        if verb:
            return "needs-design: unsettled decision AC ({!r})".format(
                verb.group(1).lower()
            )

    return None
