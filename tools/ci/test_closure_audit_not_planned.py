# tools/ci/test_closure_audit_not_planned.py
#
# Issue #1137 - the closure-audit bot applied the COMPLETED-close checklist to
# every close, so it flagged correctly-closed Issues:
#
#   #451  NOT_PLANNED, proper closing comment  -> flagged "8 unticked AC boxes"
#   #665  parent-epic rollup, no own deliverable -> flagged "missing lab-notebook entry"
#
# Both are false by construction. A not-planned close ships no work, so its AC
# boxes are SUPPOSED to stay unticked; a parent epic has no deliverable of its own,
# so it owes no lab-notebook entry. n=3 on this shape.
#
# The polarity here matters and is the whole point (shared rule: "fail safe, not
# fail open"). This is a LINTER, not a gate: its failure mode is the false ALARM.
# A gate that cries wolf on correct behavior trains the team to skim it, and the
# bot's signal-to-noise ratio IS its value. So every exemption below is paired with
# a control proving the check still FIRES on the shape it exists to catch -
# otherwise "no findings" would be indistinguishable from a bot that checks nothing.
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import closure_audit as ca  # noqa: E402

UNTICKED = "## Acceptance criteria\n- [ ] one\n- [ ] two\n\n**Priority rationale:** P2 because.\n"


def _issue(*, state_reason=None, body=UNTICKED, comments=(), sub_total=0, labels=("role:developer",)):
    return {
        "number": 1,
        "body": body,
        "labels": [{"name": n} for n in labels],
        "comments": [{"body": c} for c in comments],
        "closedAt": "2026-07-10T00:00:00Z",
        "stateReason": state_reason,
        "subIssuesSummary": {"total": sub_total, "completed": sub_total},
    }


class TestNotPlanned:
    def test_not_planned_close_is_not_flagged_for_unticked_acs(self):
        """The #451 shape. Unticked ACs on a not-planned close are CORRECT."""
        issue = _issue(state_reason="NOT_PLANNED", comments=["Closing: superseded by X. Outcome: no follow-up."])
        ac, pr, nb, cc = ca.collect_issue_gaps(issue, notebooks={})
        assert ac == [], "a not-planned close ships no work; its AC boxes are supposed to stay unticked"
        assert nb == [], "a not-planned close routes through a closing comment, not a notebook entry"
        assert cc == [], "it HAS a closing comment"

    def test_completed_close_with_unticked_acs_is_still_flagged(self):
        """CONTROL. Flip only stateReason: the AC check must still fire."""
        issue = _issue(state_reason="COMPLETED")
        ac, pr, nb, cc = ca.collect_issue_gaps(issue, notebooks={})
        assert ac, "the AC check must still catch a COMPLETED close with unticked boxes"

    def test_not_planned_without_any_closing_comment_is_flagged(self):
        """What a not-planned close DOES owe: a reason + outcome routing."""
        issue = _issue(state_reason="NOT_PLANNED", comments=[])
        ac, pr, nb, cc = ca.collect_issue_gaps(issue, notebooks={})
        assert cc, "a not-planned close with no closing comment at all has no recorded reason"


class TestParentRollup:
    def test_parent_rollup_owes_no_lab_notebook_entry(self):
        """The #665 shape: an epic closes as a rollup of children that each shipped."""
        issue = _issue(state_reason="COMPLETED", body="## Acceptance criteria\n- [x] done\n\n**Priority rationale:** P2 x.\n", sub_total=3)
        ac, pr, nb, cc = ca.collect_issue_gaps(issue, notebooks={"developer": ""})
        assert nb == [], "a parent epic has no deliverable of its own, so it owes no notebook entry"

    def test_leaf_completed_close_still_needs_a_notebook_entry(self):
        """CONTROL. Flip only sub_total: the notebook check must still fire."""
        issue = _issue(state_reason="COMPLETED", body="## Acceptance criteria\n- [x] done\n\n**Priority rationale:** P2 x.\n", sub_total=0)
        ac, pr, nb, cc = ca.collect_issue_gaps(issue, notebooks={"developer": ""})
        assert nb, "a leaf Issue that shipped work still owes its role's lab-notebook entry"


class TestPredicates:
    def test_is_not_planned_is_case_insensitive_and_none_safe(self):
        assert ca.is_not_planned({"stateReason": "not_planned"})
        assert ca.is_not_planned({"stateReason": "NOT_PLANNED"})
        assert not ca.is_not_planned({"stateReason": "COMPLETED"})
        assert not ca.is_not_planned({"stateReason": None})
        assert not ca.is_not_planned({})

    def test_is_parent_rollup_is_none_safe(self):
        assert ca.is_parent_rollup({"subIssuesSummary": {"total": 2}})
        assert not ca.is_parent_rollup({"subIssuesSummary": {"total": 0}})
        assert not ca.is_parent_rollup({"subIssuesSummary": None})
        assert not ca.is_parent_rollup({})


# --- Issue #1137 PR review: the subIssuesSummary fetch must never crash the audit ---


def test_fetch_sub_total_fails_open_on_gh_error(monkeypatch):
    """The bot's blocking scenario: if `gh issue view --json subIssuesSummary` is
    rejected (unknown field on some gh build), audit_issue must NOT crash on every
    closed issue. It degrades to total=0 (treat as leaf), never raises."""
    import subprocess as sp

    def boom(*a, **k):
        raise sp.CalledProcessError(1, "gh", stderr='Unknown JSON field: "subIssuesSummary"')

    monkeypatch.setattr(ca, "_gh", boom)
    assert ca._fetch_sub_total(1) == 0  # fail-open, not a raised exception


def test_fetch_issue_survives_a_subissues_field_rejection(monkeypatch):
    """End-to-end: the CORE fields still fetch even when the subIssues field errors.

    Matched-pair control against the fail-open above: the main fetch must succeed
    (returning a usable issue dict) rather than the whole call blowing up."""
    import subprocess as sp

    def fake_gh(*args, repo=None):
        if "subIssuesSummary" in args:  # the isolated second fetch
            raise sp.CalledProcessError(1, "gh", stderr='Unknown JSON field')
        return '{"number": 1, "body": "b", "labels": [], "comments": [], "closedAt": "2026-07-10T00:00:00Z", "stateReason": "COMPLETED"}'

    monkeypatch.setattr(ca, "_gh", fake_gh)
    issue = ca.fetch_issue(1)
    assert issue["number"] == 1
    assert issue["subIssuesSummary"]["total"] == 0  # degraded to leaf, no crash
    assert ca.is_parent_rollup(issue) is False
