"""Unit tests for scripts/pm/recheck_parent_status.py."""

import sys
from pathlib import Path

# Make the script importable as a module
SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts" / "pm"
sys.path.insert(0, str(SCRIPT_DIR))

import recheck_parent_status as rps
from closure_audit import scan_ac_boxes


class TestStatusLadder:
    def test_known_statuses_ranked(self):
        assert rps.rank("Backlog") == 0
        assert rps.rank("Ready") == 1
        assert rps.rank("In progress") == 2
        assert rps.rank("Ready for review") == 3
        assert rps.rank("In review") == 4
        assert rps.rank("Done") == 5

    def test_none_treated_as_backlog(self):
        assert rps.rank(None) == 0

    def test_unknown_status_treated_as_backlog(self):
        assert rps.rank("Mystery") == 0


class TestCollectiveState:
    def test_no_open_children_returns_done(self):
        children = []
        assert rps.collective_state(children) == "Done"

    def test_max_rank_across_children(self):
        children = [
            {"number": 100, "status": "Backlog"},
            {"number": 101, "status": "In progress"},
            {"number": 102, "status": "Ready"},
        ]
        assert rps.collective_state(children) == "In progress"

    def test_in_review_beats_in_progress(self):
        children = [
            {"number": 100, "status": "In progress"},
            {"number": 101, "status": "In review"},
        ]
        assert rps.collective_state(children) == "In review"

    def test_none_status_does_not_raise(self):
        children = [
            {"number": 100, "status": None},
            {"number": 101, "status": "Ready"},
        ]
        assert rps.collective_state(children) == "Ready"


class TestDriftClassification:
    def test_forward_drift_parent_ahead(self):
        # The case caught 2026-05-19: parent In progress, all children Backlog
        children = [{"number": 204, "status": "Backlog"}]
        result = rps.classify_drift(parent_status="In progress", open_children=children)
        assert result == "FORWARD DRIFT"

    def test_backward_drift_children_ahead(self):
        # Child In review but parent still Backlog
        children = [{"number": 100, "status": "In review"}]
        result = rps.classify_drift(parent_status="Backlog", open_children=children)
        assert result == "BACKWARD DRIFT"

    def test_completion_drift_all_children_closed_parent_not_done(self):
        # Empty open_children + parent not Done
        result = rps.classify_drift(parent_status="In progress", open_children=[])
        assert result == "COMPLETION DRIFT"

    def test_no_drift_when_aligned(self):
        children = [{"number": 100, "status": "In progress"}]
        result = rps.classify_drift(parent_status="In progress", open_children=children)
        assert result is None

    def test_no_drift_when_parent_done_and_all_children_closed(self):
        result = rps.classify_drift(parent_status="Done", open_children=[])
        assert result is None

    def test_forward_drift_requires_open_child(self):
        # Edge case: no open children but parent claims progress → COMPLETION not FORWARD
        result = rps.classify_drift(parent_status="In progress", open_children=[])
        assert result == "COMPLETION DRIFT"  # not FORWARD DRIFT

    def test_ready_parent_backlog_children_is_not_drift(self):
        # Normal state: epic groomed (Ready), sub-issues not yet groomed (Backlog).
        # Parent is not claiming active work, so not flagged as drift.
        children = [{"number": 100, "status": "Backlog"}]
        result = rps.classify_drift(parent_status="Ready", open_children=children)
        assert result is None

    # --- A2 epic-park (#776 / #794) -------------------------------------------
    # A parent parked in the off-ladder `Epic` Status no longer mirrors its
    # children's collective ladder rank. Because `Epic` is unranked, the old code
    # read it as rank-0 (Backlog) and spuriously flagged BACKWARD DRIFT against any
    # active child — fighting the park. Suppress the ladder mirror for Epic parents.
    def test_epic_parent_with_in_progress_child_is_not_drift(self):
        children = [{"number": 100, "status": "In progress"}]
        result = rps.classify_drift(parent_status="Epic", open_children=children)
        assert result is None

    def test_epic_parent_with_in_review_child_is_not_drift(self):
        children = [{"number": 100, "status": "In review"}]
        result = rps.classify_drift(parent_status="Epic", open_children=children)
        assert result is None

    def test_epic_parent_all_children_closed_still_flags_completion(self):
        # Preserved: A2 closes a completed parent to Done, so "all closed but not
        # Done" is still a close-the-parent signal — it is NOT a leaf-status mirror.
        result = rps.classify_drift(parent_status="Epic", open_children=[])
        assert result == "COMPLETION DRIFT"

    def test_epic_parent_all_closed_with_not_planned_child_still_reviews(self):
        # Preserved: the #632 scope-verification flag is not a status mirror, so it
        # must survive the A2 narrowing.
        result = rps.classify_drift(parent_status="Epic", open_children=[],
                                    has_not_planned=True)
        assert result == rps.NOT_PLANNED_REVIEW

    def test_epic_parent_with_backlog_only_children_is_not_drift(self):
        # Completeness guard: Epic + no-progress (Backlog) children. Never broke
        # (rank Epic == rank Backlog == 0), but documents the parked-no-progress case.
        children = [{"number": 100, "status": "Backlog"}]
        result = rps.classify_drift(parent_status="Epic", open_children=children)
        assert result is None

    def test_epic_parent_with_mixed_children_is_not_drift(self):
        # Belt-and-suspenders: collective_state() takes the max rank (In progress
        # here), which the pre-fix code flagged as BACKWARD DRIFT. The Epic guard
        # suppresses regardless of how the collective is computed.
        children = [{"number": 100, "status": "Backlog"},
                    {"number": 101, "status": "In progress"}]
        assert rps.collective_state(children) == "In progress"   # max-rank, would have drifted
        result = rps.classify_drift(parent_status="Epic", open_children=children)
        assert result is None


from unittest.mock import patch


class TestGhHelpers:
    @patch("recheck_parent_status.gh")
    def test_parent_issue_number_extracts_from_url(self, mock_gh):
        mock_gh.return_value = {
            "number": 204,
            "parent_issue_url": "https://api.github.com/repos/Jin-HoMLee/splice-neoepitope-pipeline/issues/86",
        }
        assert rps.parent_issue_number(204) == 86

    @patch("recheck_parent_status.gh")
    def test_parent_issue_number_returns_none_when_no_parent(self, mock_gh):
        mock_gh.return_value = {"number": 24, "parent_issue_url": None}
        assert rps.parent_issue_number(24) is None

    @patch("recheck_parent_status.gh")
    def test_open_sub_issues_filters_closed(self, mock_gh):
        mock_gh.return_value = [
            {"number": 204, "state": "open"},
            {"number": 205, "state": "closed"},
            {"number": 206, "state": "open"},
        ]
        result = rps.open_sub_issues(86)
        assert [c["number"] for c in result] == [204, 206]


class TestStatusLookup:
    @patch("recheck_parent_status.gh")
    def test_returns_status_name_when_present(self, mock_gh):
        mock_gh.return_value = {
            "data": {"repository": {"issue": {
                "projectItems": {"nodes": [{
                    "project": {"number": 9},
                    "fieldValues": {"nodes": [
                        {"field": {"name": "Status"}, "name": "In progress"},
                        {"field": {"name": "Size"}, "name": "M"},
                    ]},
                }]},
            }}}
        }
        assert rps.status_for_issue(86) == "In progress"

    @patch("recheck_parent_status.gh")
    def test_returns_none_when_not_on_project(self, mock_gh):
        mock_gh.return_value = {
            "data": {"repository": {"issue": {
                "projectItems": {"nodes": []},
            }}}
        }
        assert rps.status_for_issue(86) is None

    @patch("recheck_parent_status.gh")
    def test_skips_other_projects(self, mock_gh):
        mock_gh.return_value = {
            "data": {"repository": {"issue": {
                "projectItems": {"nodes": [{
                    "project": {"number": 99},
                    "fieldValues": {"nodes": [
                        {"field": {"name": "Status"}, "name": "In progress"},
                    ]},
                }]},
            }}}
        }
        assert rps.status_for_issue(86) is None


class TestAuditChain:
    @patch("recheck_parent_status.status_for_issue")
    @patch("recheck_parent_status.open_sub_issues")
    @patch("recheck_parent_status.parent_issue_number")
    def test_walks_two_level_chain(self, mock_parent, mock_subs, mock_status):
        # Simulate today's case: #204 is sub of #86 is sub of #24
        mock_parent.side_effect = lambda n: {204: 86, 86: 24, 24: None}[n]
        mock_subs.side_effect = lambda n: {
            86: [{"number": 204}, {"number": 205}, {"number": 206}],
            24: [{"number": 86}],
        }[n]
        # All sub-issues stale; both parents show In progress (drift)
        def status_side_effect(n):
            return {86: "In progress", 24: "In progress",
                    204: "Backlog", 205: "Backlog", 206: "Backlog"}.get(n)
        mock_status.side_effect = status_side_effect

        chain = rps.audit_parent_chain(204)
        # We expect 2 audit records: #86 and #24
        assert [r["issue"] for r in chain] == [86, 24]
        assert chain[0]["drift"] == "FORWARD DRIFT"
        # #24's immediate drift is None: its only open child #86 has status
        # "In progress" matching #24's "In progress". #86's drift is transitive
        # (caught at chain[0]) — classify_drift only flags per-level drift.
        assert chain[1]["drift"] is None

    @patch("recheck_parent_status.status_for_issue")
    @patch("recheck_parent_status.open_sub_issues")
    @patch("recheck_parent_status.parent_issue_number")
    def test_stops_at_root(self, mock_parent, mock_subs, mock_status):
        mock_parent.side_effect = lambda n: {100: None}[n]
        # Issue 100 has no parent → empty chain
        chain = rps.audit_parent_chain(100)
        assert chain == []


class TestFormatter:
    def test_no_drift_record_renders_no_change(self):
        record = {
            "issue": 24,
            "status": "Ready",
            "open_children": [{"number": 86, "status": "Ready"}],
            "collective": "Ready",
            "drift": None,
        }
        out = rps.format_record(record)
        assert "#24" in out
        assert "Ready" in out
        assert "[No change]" in out

    def test_forward_drift_record_renders_label(self):
        record = {
            "issue": 86,
            "status": "In progress",
            "open_children": [
                {"number": 204, "status": "Backlog"},
                {"number": 205, "status": "Backlog"},
            ],
            "collective": "Backlog",
            "drift": "FORWARD DRIFT",
        }
        out = rps.format_record(record)
        assert "#86" in out
        assert "[FORWARD DRIFT]" in out
        assert "#204 (Backlog)" in out


class TestCli:
    @patch("recheck_parent_status.audit_parent_chain")
    def test_issue_mode_exit_0_when_no_drift(self, mock_audit, capsys):
        mock_audit.return_value = [{
            "issue": 24, "status": "Ready",
            "open_children": [{"number": 86, "status": "Ready"}],
            "collective": "Ready", "drift": None,
        }]
        rc = rps.main(["--issue", "204"])
        assert rc == 0
        captured = capsys.readouterr()
        assert "#24" in captured.out
        assert "[No change]" in captured.out

    @patch("recheck_parent_status.audit_parent_chain")
    def test_issue_mode_exit_2_when_drift(self, mock_audit, capsys):
        mock_audit.return_value = [{
            "issue": 86, "status": "In progress",
            "open_children": [{"number": 204, "status": "Backlog"}],
            "collective": "Backlog", "drift": "FORWARD DRIFT",
        }]
        rc = rps.main(["--issue", "204"])
        assert rc == 2
        captured = capsys.readouterr()
        assert "[FORWARD DRIFT]" in captured.out

    @patch("recheck_parent_status.audit_parent_chain")
    def test_issue_with_no_parent_emits_message(self, mock_audit, capsys):
        mock_audit.return_value = []
        rc = rps.main(["--issue", "100"])
        assert rc == 0
        captured = capsys.readouterr()
        assert "no parent" in captured.out.lower()


class TestAllMode:
    @patch("recheck_parent_status.status_for_issue")
    @patch("recheck_parent_status.open_sub_issues")
    @patch("recheck_parent_status.all_parent_issues")
    def test_iterates_only_parents_skips_clean(self, mock_parents, mock_subs, mock_status):
        # Two parents: one drifted, one clean
        mock_parents.return_value = [86, 24]
        mock_subs.side_effect = lambda n: {
            86: [{"number": 204}],
            24: [{"number": 86}],
        }[n]
        status_map = {86: "In progress", 24: "Ready", 204: "Backlog"}
        mock_status.side_effect = lambda n: status_map.get(n)

        rc = rps.run_all_mode()
        assert rc == 2  # at least one drift

    @patch("recheck_parent_status.status_for_issue")
    @patch("recheck_parent_status.open_sub_issues")
    @patch("recheck_parent_status.all_parent_issues")
    def test_returns_0_when_no_drift(self, mock_parents, mock_subs, mock_status):
        mock_parents.return_value = [24]
        mock_subs.return_value = [{"number": 86}]
        mock_status.side_effect = lambda n: "Ready"  # parent + sub both Ready
        rc = rps.run_all_mode()
        assert rc == 0

    @patch("recheck_parent_status.status_for_issue")
    @patch("recheck_parent_status.open_sub_issues")
    @patch("recheck_parent_status.all_parent_issues")
    def test_terminal_gh_error_on_one_parent_does_not_abort_sweep(
        self, mock_parents, mock_subs, mock_status, capsys
    ):
        # Issue #1017 AC #2: a terminal gh failure (retries exhausted -> GhError)
        # on one parent must skip that parent, not abort the remaining audit.
        mock_parents.return_value = [86, 24]

        def status(n):
            if n == 86:  # this parent's probe fails terminally
                raise rps.GhError(1, ["gh"], stderr="HTTP 503: server error")
            return {24: "In progress", 204: "Backlog"}.get(n)

        mock_status.side_effect = status
        mock_subs.side_effect = lambda n: {24: [{"number": 204}]}.get(n, [])

        # If #86's GhError aborted the sweep this call would raise, not return.
        rc = rps.run_all_mode()

        # #24's In-progress-over-Backlog-child drift is still detected past the skip.
        assert rc == 2
        err = capsys.readouterr().err
        assert "skipped #86" in err  # the skip was surfaced, not silently dropped

    @patch("recheck_parent_status.status_for_issue")
    @patch("recheck_parent_status.open_sub_issues")
    @patch("recheck_parent_status.all_parent_issues")
    def test_all_parents_skipped_returns_error_not_zero(
        self, mock_parents, mock_subs, mock_status
    ):
        # Issue #1017 review: a broad outage that skips EVERY parent audits nothing;
        # exit 0 would be indistinguishable from a clean board. Must exit 1 instead.
        mock_parents.return_value = [86, 24]

        def boom(n):
            raise rps.GhError(1, ["gh"], stderr="HTTP 503: server error")

        mock_status.side_effect = boom
        mock_subs.side_effect = lambda n: []
        rc = rps.run_all_mode()
        assert rc == 1  # audited nothing -> error, not a false "clean"


class TestNotPlannedDrift:
    """Issue #632: a parent whose only remaining children are all closed must
    distinguish COMPLETED (delivered) from NOT_PLANNED (deferred) closes. A
    NOT_PLANNED child means the bare [COMPLETION DRIFT] rollup-close suggestion
    is unsafe — emit the softer verify prompt instead."""

    # The exact FINAL rendered flag the spec requires (single-bracketed by
    # format_record). Spelled out literally here so a typo / em-dash drift in
    # the source constant fails the test.
    RENDERED = (
        "Status: [REVIEW: parent has a not-planned child — "
        "verify scope was delivered, not deferred]"
    )

    def test_classify_not_planned_child_emits_review_not_completion(self):
        # The #24 mixed shape: parent claims progress, all children closed,
        # ≥1 closed NOT_PLANNED → softer REVIEW flag, not bare COMPLETION DRIFT.
        result = rps.classify_drift(
            parent_status="In progress", open_children=[], has_not_planned=True
        )
        assert result == rps.NOT_PLANNED_REVIEW
        assert "COMPLETION DRIFT" not in (result or "")

    def test_classify_all_completed_still_flags_completion_drift(self):
        # Regression guard: all children closed COMPLETED → COMPLETION DRIFT.
        result = rps.classify_drift(
            parent_status="In progress", open_children=[], has_not_planned=False
        )
        assert result == "COMPLETION DRIFT"

    def test_done_parent_with_not_planned_child_still_no_drift(self):
        # An already-Done parent never had completion drift; REVIEW only
        # substitutes for the flag that would otherwise fire.
        result = rps.classify_drift(
            parent_status="Done", open_children=[], has_not_planned=True
        )
        assert result is None

    @patch("recheck_parent_status.gh")
    def test_has_not_planned_child_true_for_mixed_set(self, mock_gh):
        # Real /sub_issues shape confirmed by live probe: state + state_reason.
        mock_gh.return_value = [
            {"number": 145, "state": "closed", "state_reason": "completed"},
            {"number": 86, "state": "closed", "state_reason": "completed"},
            {"number": 192, "state": "closed", "state_reason": "not_planned"},
        ]
        assert rps.has_not_planned_child(24) is True

    @patch("recheck_parent_status.gh")
    def test_has_not_planned_child_false_when_all_completed(self, mock_gh):
        mock_gh.return_value = [
            {"number": 145, "state": "closed", "state_reason": "completed"},
            {"number": 86, "state": "closed", "state_reason": "completed"},
        ]
        assert rps.has_not_planned_child(24) is False

    def test_format_record_single_brackets_review_flag(self):
        record = {
            "issue": 24,
            "status": "In progress",
            "open_children": [],
            "collective": "Done",
            "drift": rps.NOT_PLANNED_REVIEW,
        }
        out = rps.format_record(record)
        assert self.RENDERED in out
        assert "[[" not in out  # brackets owned solely by format_record

    @patch("recheck_parent_status.has_not_planned_child")
    @patch("recheck_parent_status.status_for_issue")
    @patch("recheck_parent_status.open_sub_issues")
    @patch("recheck_parent_status.parent_issue_number")
    def test_audit_chain_issue_path_emits_review(
        self, mock_parent, mock_subs, mock_status, mock_np
    ):
        # --issue path (the close-hook's production-critical route): closing
        # #86 walks to parent #24, whose children are all closed with one
        # NOT_PLANNED (#192) → REVIEW, not COMPLETION DRIFT.
        mock_parent.side_effect = lambda n: {86: 24, 24: None}[n]
        mock_subs.side_effect = lambda n: {24: []}[n]  # all #24 children closed
        mock_status.side_effect = lambda n: {24: "In progress"}.get(n)
        mock_np.return_value = True

        chain = rps.audit_parent_chain(86)
        assert [r["issue"] for r in chain] == [24]
        assert chain[0]["drift"] == rps.NOT_PLANNED_REVIEW
        assert self.RENDERED in rps.format_record(chain[0])

    @patch("recheck_parent_status.has_not_planned_child")
    @patch("recheck_parent_status.status_for_issue")
    @patch("recheck_parent_status.open_sub_issues")
    @patch("recheck_parent_status.parent_issue_number")
    def test_audit_chain_all_completed_still_completion_drift(
        self, mock_parent, mock_subs, mock_status, mock_np
    ):
        mock_parent.side_effect = lambda n: {86: 24, 24: None}[n]
        mock_subs.side_effect = lambda n: {24: []}[n]
        mock_status.side_effect = lambda n: {24: "In progress"}.get(n)
        mock_np.return_value = False

        chain = rps.audit_parent_chain(86)
        assert chain[0]["drift"] == "COMPLETION DRIFT"

    @patch("recheck_parent_status.has_not_planned_child")
    @patch("recheck_parent_status.status_for_issue")
    @patch("recheck_parent_status.open_sub_issues")
    @patch("recheck_parent_status.all_parent_issues")
    def test_all_mode_emits_review_for_not_planned(
        self, mock_parents, mock_subs, mock_status, mock_np, capsys
    ):
        # --all sweep (the second call site): same NOT_PLANNED parent must get
        # the REVIEW flag, proving run_all_mode is wired too.
        mock_parents.return_value = [24]
        mock_subs.side_effect = lambda n: {24: []}[n]
        mock_status.side_effect = lambda n: "In progress"
        mock_np.return_value = True

        rc = rps.run_all_mode()
        assert rc == 2  # REVIEW is a non-None drift → counted as a fire
        assert self.RENDERED in capsys.readouterr().out

    @patch("recheck_parent_status.has_not_planned_child")
    @patch("recheck_parent_status.status_for_issue")
    @patch("recheck_parent_status.open_sub_issues")
    @patch("recheck_parent_status.all_parent_issues")
    def test_all_mode_all_completed_still_completion_drift(
        self, mock_parents, mock_subs, mock_status, mock_np, capsys
    ):
        mock_parents.return_value = [24]
        mock_subs.side_effect = lambda n: {24: []}[n]
        mock_status.side_effect = lambda n: "In progress"
        mock_np.return_value = False

        rc = rps.run_all_mode()
        assert rc == 2
        assert "Status: [COMPLETION DRIFT]" in capsys.readouterr().out


# --- #1067: full native sub-bar, unticked body ACs --------------------------
#
# Under the A2 park, parent progress is read off the native sub-issue bar, so a
# full bar reads as "done" at a glance. Body scope routinely exceeds the filed
# sub-issues (design/docs children close first; implementation ACs remain), so a
# full bar with unticked body ACs is an INVISIBLE drift class - the bar says
# complete, the body says otherwise. Advisory, mirroring the #632 softening.


class TestBodyAcReview:
    def test_full_bar_with_unticked_body_boxes_flags_review(self):
        drift = rps.classify_drift(
            "Epic", [], body_unticked_count=4, body_unticked_where="4 under 'Sub-issues'"
        )
        assert drift is not None
        assert "4 body box(es) still unticked" in drift
        assert "Sub-issues" in drift          # names WHERE, so a human can judge
        assert drift != "COMPLETION DRIFT"

    def test_full_bar_with_everything_ticked_is_still_ready_to_close(self):
        # Matched-pair control: same inputs, box count flipped, opposite outcome.
        # If the new branch were dead, BOTH would return COMPLETION DRIFT.
        assert rps.classify_drift("Epic", [], body_unticked_count=0) == "COMPLETION DRIFT"

    def test_the_real_859_body_shape_fires(self):
        """The regression that the first design silently would NOT have caught.

        #859 is a motivating case named in Issue #1067, and its four unticked boxes
        live under `## Sub-issues`, not `## Acceptance criteria`. Keying the flag on
        AC boxes alone returned **zero on its own motivating example** - a no-op that
        looked green. Caught by running it against the live board before shipping.

        Falsifier: key `body_unticked` on `scan.ac_unticked` only and this goes red.
        """
        body = (
            "## Triage (28 sections, 4 verdicts)\n\nsome prose\n\n"
            "## Sub-issues (one PR each)\n\n"
            "- [ ] Pilot: authoring-research-decks + stub\n"
            "- [ ] running-snakemake, splice-pipeline-gotchas\n"
            "- [ ] Trim pass: Safety Wrappers to pointers\n"
            "- [ ] Final: prune AGENTS.md\n"
        )
        scan = scan_ac_boxes(body)
        assert scan.ac_unticked == 0, "the boxes are NOT under an AC heading"
        assert scan.stray_unticked == 4, "...but there are four of them"
        # So the flag must key on the total, not on ac_unticked.
        total = scan.ac_unticked + scan.stray_unticked
        drift = rps.classify_drift("Epic", [], body_unticked_count=total,
                                   body_unticked_where="4 under 'Sub-issues'")
        assert drift is not None and drift.startswith("REVIEW:")

    def test_not_planned_child_outranks_the_body_prompt(self):
        # Deferred scope is the stronger reason to distrust "done", so #632 wins.
        drift = rps.classify_drift(
            "Epic", [], has_not_planned=True, body_unticked_count=3
        )
        assert drift == rps.NOT_PLANNED_REVIEW

    def test_cross_repo_parent_with_invisible_children_is_not_flagged(self):
        # AC 3. "All children closed" and "no children visible" both arrive as an
        # empty list; without the guard we would flag a bar we never actually read.
        assert rps.classify_drift(
            "Epic", [], body_unticked_count=4, has_known_children=False
        ) == "COMPLETION DRIFT"

    def test_body_boxes_do_not_affect_a_parent_with_open_children(self):
        open_kids = [{"number": 1, "status": "In progress"}]
        assert rps.classify_drift("Epic", open_kids, body_unticked_count=9) is None

    def test_closed_parent_stays_clean(self):
        assert rps.classify_drift("Done", [], body_unticked_count=5) is None

    def test_advisory_never_blocks(self):
        drift = rps.classify_drift("Epic", [], body_unticked_count=1,
                                   body_unticked_where="1 under 'Tasks'")
        assert drift.startswith("REVIEW:")


class TestAllModeBodyScope:
    """The --all sweep must fire the #1067 flag too (PR #1132 review, finding 1).

    The flag was wired into `audit_parent_chain` (--issue) but NOT `run_all_mode`
    (--all), which inlined its own `classify_drift` call and so silently defaulted
    every new argument. --all is the **proactive board-wide discovery mode** - the
    mode that surfaces pre-existing full-bar/unticked-body parents, and the very one
    used to find this feature's motivating cases. Shipping it blind would have made a
    board sweep report clean while #859 and #527 sat right there.

    Mirrors `test_all_mode_emits_review_for_not_planned`, whose existence is exactly
    why the NOT_PLANNED flag did not slip the same way.
    """

    @patch("recheck_parent_status.body_unticked")
    @patch("recheck_parent_status.all_sub_issues")
    @patch("recheck_parent_status.has_not_planned_child")
    @patch("recheck_parent_status.status_for_issue")
    @patch("recheck_parent_status.open_sub_issues")
    @patch("recheck_parent_status.all_parent_issues")
    def test_all_mode_emits_body_scope_review(
        self, mock_parents, mock_open, mock_status, mock_np, mock_all, mock_body, capsys
    ):
        mock_parents.return_value = [859]
        mock_open.return_value = []                       # native bar full
        mock_all.return_value = [{"number": 1, "state": "closed"}]  # children ARE visible
        mock_status.return_value = "Epic"
        mock_np.return_value = False
        mock_body.return_value = (4, "4 under 'Sub-issues'")

        rc = rps.run_all_mode()
        out = capsys.readouterr().out
        assert rc == 2, "a REVIEW is a non-None drift and must be counted"
        assert "4 body box(es) still unticked" in out
        assert "Sub-issues" in out

    @patch("recheck_parent_status.body_unticked")
    @patch("recheck_parent_status.all_sub_issues")
    @patch("recheck_parent_status.has_not_planned_child")
    @patch("recheck_parent_status.status_for_issue")
    @patch("recheck_parent_status.open_sub_issues")
    @patch("recheck_parent_status.all_parent_issues")
    def test_all_mode_clean_body_still_completion_drift(
        self, mock_parents, mock_open, mock_status, mock_np, mock_all, mock_body, capsys
    ):
        # Matched pair: identical sweep, box count flipped -> the OTHER outcome.
        # If run_all_mode were still unwired, both cases would print COMPLETION DRIFT.
        mock_parents.return_value = [859]
        mock_open.return_value = []
        mock_all.return_value = [{"number": 1, "state": "closed"}]
        mock_status.return_value = "Epic"
        mock_np.return_value = False
        mock_body.return_value = (0, "")

        rps.run_all_mode()
        out = capsys.readouterr().out
        assert "COMPLETION DRIFT" in out
        assert "still unticked" not in out
