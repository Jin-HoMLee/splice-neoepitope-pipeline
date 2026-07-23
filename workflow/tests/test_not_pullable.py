"""Tests for scripts/pm/not_pullable.py - the body-only "not pullable" scan.

Issue #1248. The Ready-queue floor guard counts Ready *depth* and is blind to
body-only gates, so a decision-blocked or trigger-gated queue reads as healthy
and floor-pressure Replenishment re-commits Issues that cannot actually be
worked. #841 was swept into Ready twice on exactly this shape; #929 sat in Ready
with two `needs-design` decision-fork ACs.

Two distinct tells, because the two fixtures fail differently:

  - TRIGGER GATE  - a prose marker anywhere in the body ("Build only if it
    slips a 2nd time" lives in #841's `## Notes`, not in its ACs).
  - NEEDS DESIGN  - decision verbs inside an UNTICKED acceptance-criterion
    line (#929's "Decide scope:" / "Decide home:").

The false-positive discipline (AC 4) is the load-bearing part. An Issue that
merely *documents* the markers - like #1248 itself, which lists them in
backticks - must NOT trip. That is the same "using a directive vs talking about
it" hazard that bit the `skip-lab-notebook` marker in Issue #1126, and the fix
is the same in spirit: make the scan blind to text that is quoted as code.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "pm"))

from not_pullable import scan_not_pullable  # noqa: E402


# --------------------------------------------------------------------------
# Fixtures modelled on the real Issue bodies named in the #1248 ACs.
# --------------------------------------------------------------------------

# Issue #841 - the trigger gate lives in `## Notes`, in prose, NOT in the ACs.
BODY_841 = """## Problem

The closure-ritual gate verifies a review was **triggered** but never that it
completed.

## Acceptance criteria

- [ ] Gate prints the triggered review's summary/findings to stderr
- [ ] Non-interactive merge blocks when the review is incomplete

## Notes

- This is the **rung-3 mechanism escalation** candidate for the read-before-merge
  rule, currently **incident #1**. Build only if it slips a 2nd time, per the
  mechanism-over-memory ladder - filing now so the escalation is tracked, not to
  build immediately.
"""

# Issue #929 - decision-fork ACs, unticked.
BODY_929 = """## Problem

A semantic critic pass at the merge gate may catch what mechanical gates cannot.

## Acceptance criteria
- [ ] Decide scope: which merges warrant a semantic judge pass vs mechanical-only.
- [ ] Decide home: new PM tooling vs fold into #353.
- [ ] Prototype an isolated-context judge prompt against 1-2 recent merges.
"""

# Issue #1248 itself - it LISTS the markers, in backticks, while proposing the
# scan. It must not trip its own detector.
BODY_1248 = """## Problem

The Ready-queue floor guard counts Ready depth and is blind to body-only gates.

So a decision-blocked or trigger-gated queue reads as healthy, and Replenishment
re-commits Issues that cannot be worked. Issue #841 is chartered to build only on
a 2nd trigger event; it was swept into Ready twice.

## Proposed

Extend the body-scan to also catch trigger-gate markers (`build only if`,
`do not commit before`, `TRIGGER-GATED`, `revival trigger`), then have the floor
guard treat a body-gated Ready item as not-pullable.

## Acceptance criteria

- [ ] The guard detects the greppable not-pullable tells in a candidate's body.
- [ ] A body-gated Ready item is excluded from the pullable-floor count.
"""

BODY_CLEAN = """## Problem

`star_sj_to_junctions.py` reads column 7 only, silently applying a unique-only
policy.

## Acceptance criteria

- [ ] The uniqueness semantic is governed by one explicit config knob.
- [ ] Both aligner paths agree on multimapper handling.
"""


# --------------------------------------------------------------------------
# The two tells fire (AC 1 + AC 3 - #841 and #929 as regression fixtures)
# --------------------------------------------------------------------------

def test_trigger_gate_marker_in_prose_notes_trips():
    """#841: 'Build only if it slips a 2nd time' is a gate, and it is not in an AC."""
    reason = scan_not_pullable(BODY_841)
    assert reason is not None, "#841's trigger gate must be detected"
    assert "trigger" in reason.lower()


def test_needs_design_decision_acs_trip():
    """#929: unticked 'Decide ...' ACs mean the scope is not settled."""
    reason = scan_not_pullable(BODY_929)
    assert reason is not None, "#929's decision-fork ACs must be detected"
    assert "design" in reason.lower() or "decision" in reason.lower()


@pytest.mark.parametrize("verb", ["Decide", "Determine", "Choose"])
def test_each_decision_verb_in_unticked_ac_trips(verb):
    """Every live branch of the decision-verb alternation has direct coverage.

    Bot review on PR #1292 caught that only `decide` was exercised (via
    BODY_929), leaving the other branches able to be dropped from the regex
    while the suite stayed green. `determine` matters most: it is the verb at
    the centre of the #659 governance tension, so it is the branch most likely
    to be re-tuned and was the least protected.
    """
    body = "## Acceptance criteria\n\n- [ ] {} the scope of the work.\n".format(verb)
    assert scan_not_pullable(body) is not None, "verb not detected: {}".format(verb)


@pytest.mark.parametrize("verb", ["Settle", "Pick"])
def test_speculative_verbs_are_deliberately_not_matched(verb):
    """`settle` / `pick` were trimmed from the alternation, on purpose.

    They had no observed instance in any real Issue and are not among the
    decision phrasings the DoR names, so each was unexercised false-positive
    surface on a linter. This pins the trim so a future author does not
    re-add them by reflex without an observed case to justify it.
    """
    body = "## Acceptance criteria\n\n- [ ] {} the scope of the work.\n".format(verb)
    assert scan_not_pullable(body) is None


@pytest.mark.parametrize(
    "marker",
    [
        "Build only if it slips a 2nd time.",
        "Build only on a 2nd trigger event.",
        "Do not commit before ~2026-08-03 - the data does not exist yet.",
        "TRIGGER-GATED: needs a funded GPU before this is real.",
        "Revival trigger: revisit when a durable archive is genuinely needed.",
    ],
)
def test_each_trigger_marker_variant_trips(marker):
    body = f"## Problem\n\nSomething.\n\n## Notes\n\n- {marker}\n"
    assert scan_not_pullable(body) is not None, f"marker not detected: {marker}"


# --------------------------------------------------------------------------
# False-positive discipline (AC 4) - the half that makes this a real check
# --------------------------------------------------------------------------

def test_issue_documenting_the_markers_does_not_trip():
    """#1248 lists the markers in backticks while proposing the scan.

    This is the matched-pair control for the #1126 'using a directive vs
    talking about it' hazard: identical marker text, one quoted as code, and
    the expected outcomes are OPPOSITE. Without this, the scan would flag every
    Issue that discusses the convention - including the one that introduces it.
    """
    assert scan_not_pullable(BODY_1248) is None


@pytest.mark.parametrize(
    "clause, gated",
    [
        # Opens a clause -> this Issue's own gate.
        ("Build only if it slips a 2nd time.", True),
        ("- Build only on a 2nd trigger event.", True),
        ("**Build only if** the trigger fires.", True),
        # Mid-clause -> describing some *other* Issue's gate. Both of these are
        # real sentences from Issue #1248's body.
        ("Issue #841 is chartered to build only on a 2nd trigger event.", False),
        ("So a decision-blocked or trigger-gated queue reads as healthy.", False),
    ],
)
def test_marker_position_decides_gate_vs_description(clause, gated):
    """Matched-pair control: identical marker words, opposite expected outcomes.

    The only variable is whether the marker opens a clause. If this ever passes
    in both directions the anchoring has stopped doing anything, which is the
    failure mode a one-sided test cannot see.
    """
    body = "## Problem\n\nContext.\n\n## Notes\n\n{}\n".format(clause)
    assert (scan_not_pullable(body) is not None) is gated


def test_live_1248_body_shape_does_not_trip():
    """Regression: Issue #1248's real body false-positived before clause anchoring.

    Its prose discusses gates twice mid-sentence AND lists the markers in
    backticks. The Issue that introduces the scan must survive its own scan.
    """
    assert scan_not_pullable(BODY_1248) is None


def test_narrative_mention_of_a_decision_outside_acs_does_not_trip():
    """A body that merely narrates a past decision is not a decision fork."""
    body = """## Problem

We decided to converge on STAR. Jin-Ho will decide the rollout order later, but
that is not gating this work.

## Acceptance criteria

- [ ] Container path runs the chr22 fixture green.
"""
    assert scan_not_pullable(body) is None


def test_ticked_decision_ac_does_not_trip():
    """A decision already made is not an open fork."""
    body = """## Acceptance criteria

- [x] Decide scope: settled on high-value merges only.
- [ ] Implement the judge pass.
"""
    assert scan_not_pullable(body) is None


def test_clean_issue_returns_none():
    assert scan_not_pullable(BODY_CLEAN) is None


def test_fenced_code_block_is_ignored():
    """A marker inside a fenced block is sample text, not a gate."""
    body = """## Proposed

Match on these:

```
build only if
do not commit before
```

## Acceptance criteria

- [ ] Implement the matcher.
"""
    assert scan_not_pullable(body) is None


# --------------------------------------------------------------------------
# Degenerate input - the scan must never be the thing that breaks a board read
# --------------------------------------------------------------------------

@pytest.mark.parametrize("body", ["", None])
def test_empty_or_missing_body_is_pullable(body):
    assert scan_not_pullable(body) is None
