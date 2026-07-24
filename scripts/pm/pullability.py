"""One pullability predicate over natively-owned sources (Issue #1294).

Replaces four scattered "can this be worked now?" checks with a single
function. Each cause is read from exactly one authoritative source and never
copied, so no state can drift. Adding a fifth cause is one branch here.

The prose scanner (not_pullable.py) is NOT consulted here. It is a low-recall
proposer that suggests a label for a human to apply; the label - not the prose -
is what this predicate reads. That is the whole point of Issue #1294: the
decision has one home, and the sources stay where they are natively mastered.
"""
from typing import Optional

# The enumerated reason taxonomy: one authoritative source per cause. Adding a
# fifth cause is one constant plus one branch in assess(), in this one file.
NEEDS_DESIGN_LABEL = "needs-design"
TRIGGER_GATED_LABEL = "trigger-gated"


def _open_blockers(item):
    """Issue numbers of the still-open native blockedBy edges."""
    out = []
    blocked_by = item.get("blocked_by")
    if not isinstance(blocked_by, list):
        return out
    for edge in blocked_by:
        if not isinstance(edge, dict):
            continue
        if (edge.get("state") or "").upper() != "CLOSED":
            num = edge.get("number")
            if num is not None:
                out.append(num)
    return out


def assess(item, *, today=None):
    # type: (dict, Optional[str]) -> Optional[str]
    """Return a short reason the Issue is not pullable now, or None if it is.

    `item` carries already-fetched structured fields:
      - "labels":     list[str] label names
      - "blocked_by": list[{"number", "state"}] native blockedBy edges
      - "start_date": Optional[str] ISO YYYY-MM-DD; a future value gates the Issue

    Sources are checked hardest-first: a hard dependency subsumes a softer gate,
    and a date gate ("not yet") is the weakest. Fails open on any unknown shape.
    """
    open_blockers = _open_blockers(item)
    if open_blockers:
        return "blocked-by-issue: #{}".format(open_blockers[0])

    labels = item.get("labels") or []
    if NEEDS_DESIGN_LABEL in labels:
        return "needs-design"
    if TRIGGER_GATED_LABEL in labels:
        return "trigger-gated"

    start = item.get("start_date")
    if start and today and start > today:
        return "date-gated: {}".format(start)

    return None
