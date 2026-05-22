"""Unit tests for scripts/pm/recheck_milestone.py."""

import sys
from datetime import date, timedelta
from pathlib import Path

# Make the script importable as a module
SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts" / "pm"
sys.path.insert(0, str(SCRIPT_DIR))

import recheck_milestone as rm


class TestParseMilestoneTitle:
    def test_standard_s_stage(self):
        assert rm.parse_milestone_title("i3 - S5 - Modeling - HLA-Matched TCR Panel") == (3, 5)

    def test_s7_publication(self):
        assert rm.parse_milestone_title("i2 - S7 - Publication - HLA-Matched TCR Panel") == (2, 7)

    def test_pm_role_meta_returns_none(self):
        assert rm.parse_milestone_title("pm-i4 - PM Tooling, Memory & Methodology") is None

    def test_dev_role_meta_returns_none(self):
        assert rm.parse_milestone_title("dev-i3 - Pipeline Architecture & Dev Tooling") is None

    def test_legacy_m_milestone_returns_none(self):
        assert rm.parse_milestone_title("M1 - Problem Definition v1") is None

    def test_double_digit_iteration(self):
        assert rm.parse_milestone_title("i10 - S3 - Data Preparation - Arc") == (10, 3)
