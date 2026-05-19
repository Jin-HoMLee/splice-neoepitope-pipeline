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
