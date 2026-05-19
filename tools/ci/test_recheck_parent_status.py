"""Unit tests for scripts/pm/recheck_parent_status.py."""

import sys
from pathlib import Path

# Make the script importable as a module
SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts" / "pm"
sys.path.insert(0, str(SCRIPT_DIR))

import recheck_parent_status as rps


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


from unittest.mock import patch, MagicMock


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
