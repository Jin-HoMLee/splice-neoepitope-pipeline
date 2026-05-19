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
