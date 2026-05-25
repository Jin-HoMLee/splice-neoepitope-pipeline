"""Unit tests for scripts/pm/recheck_milestone.py."""

import sys
from datetime import date, timedelta
from pathlib import Path

import pytest

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


class TestFindPriorSameStage:
    """Fixture: 5 milestones across S3 spanning iterations 1-5, plus one S5."""

    FIXTURE = [
        {"number": 100, "title": "i1 - S3 - Data Prep - V1", "state": "closed", "due_on": "2026-03-01T00:00:00Z"},
        {"number": 101, "title": "i2 - S3 - Data Prep - GTEx", "state": "open", "due_on": "2026-06-03T00:00:00Z"},
        {"number": 102, "title": "i3 - S3 - Data Prep - Aligner", "state": "open", "due_on": "2026-06-05T00:00:00Z"},
        {"number": 103, "title": "i4 - S3 - Data Prep - nfcore", "state": "open", "due_on": "2026-06-09T00:00:00Z"},
        {"number": 104, "title": "i5 - S3 - Data Prep - STAR", "state": "open", "due_on": "2026-06-19T00:00:00Z"},
        {"number": 200, "title": "i2 - S5 - Modeling - TCR Panel", "state": "open", "due_on": "2026-05-31T00:00:00Z"},
        {"number": 300, "title": "pm-i4 - PM Tooling", "state": "open", "due_on": "2026-06-03T00:00:00Z"},
    ]

    def test_picks_highest_iteration_prior(self):
        prior = rm.find_prior_same_stage(5, 3, self.FIXTURE)
        assert prior["number"] == 103  # i4-S3, not i1/i2/i3

    def test_picks_immediate_prior_for_i3(self):
        prior = rm.find_prior_same_stage(3, 3, self.FIXTURE)
        assert prior["number"] == 101  # i2-S3

    def test_no_prior_for_i1(self):
        # i1-S3 exists but has no prior (i0 doesn't exist)
        prior = rm.find_prior_same_stage(1, 3, self.FIXTURE)
        assert prior is None

    def test_no_prior_for_first_in_stage(self):
        # No prior S5 before i2-S5
        prior = rm.find_prior_same_stage(2, 5, self.FIXTURE)
        assert prior is None

    def test_ignores_role_meta_in_chain(self):
        # pm-i4 doesn't pollute S3 chain
        prior = rm.find_prior_same_stage(5, 3, self.FIXTURE)
        assert prior["number"] != 300

    def test_includes_closed_milestones_as_prior(self):
        # i1-S3 is closed — should still be findable when no open prior exists
        fixture = [
            {"number": 100, "title": "i1 - S3 - V1", "state": "closed", "due_on": "2026-03-01T00:00:00Z"},
            {"number": 101, "title": "i2 - S3 - V2", "state": "open", "due_on": "2026-06-03T00:00:00Z"},
        ]
        prior = rm.find_prior_same_stage(2, 3, fixture)
        assert prior["number"] == 100


class TestFindOpenSameIterationS5:
    """Used for paired-S7 gating: an S7 milestone unblocks at its paired S5 close."""

    def test_finds_open_s5_in_same_iteration(self):
        fixture = [
            {"number": 200, "title": "i2 - S5 - Modeling - TCR Panel", "state": "open", "due_on": "2026-05-31T00:00:00Z"},
            {"number": 201, "title": "i2 - S7 - Publication - TCR Panel", "state": "open", "due_on": "2026-06-04T00:00:00Z"},
        ]
        paired = rm.find_open_same_iteration_S5(2, fixture)
        assert paired["number"] == 200

    def test_returns_none_when_paired_s5_closed(self):
        fixture = [
            {"number": 200, "title": "i2 - S5 - Modeling - TCR Panel", "state": "closed", "due_on": "2026-05-31T00:00:00Z"},
        ]
        assert rm.find_open_same_iteration_S5(2, fixture) is None

    def test_returns_none_when_no_s5_in_iteration(self):
        fixture = [
            {"number": 100, "title": "i3 - S3 - Data Prep", "state": "open", "due_on": "2026-06-05T00:00:00Z"},
        ]
        assert rm.find_open_same_iteration_S5(3, fixture) is None

    def test_loose_arc_match(self):
        # i4-S7 'TCR-pMHC Landscape' paired with i4-S5 'Google Batch' — arc mismatch
        # but iteration matches, so it pairs (data hygiene concern is separate).
        fixture = [
            {"number": 200, "title": "i4 - S5 - Modeling - Google Batch", "state": "open", "due_on": "2026-06-20T00:00:00Z"},
            {"number": 201, "title": "i4 - S7 - Publication - TCR-pMHC Landscape", "state": "open", "due_on": "2026-07-31T00:00:00Z"},
        ]
        paired = rm.find_open_same_iteration_S5(4, fixture)
        assert paired["number"] == 200


class TestComputeLayeredDueDate:
    """Covers all 7 logic branches of compute_layered_due_date."""

    def _today(self, monkeypatch, fake_today: date):
        class _FakeDate(date):
            @classmethod
            def today(cls):
                return fake_today
        monkeypatch.setattr(rm, "date", _FakeDate)

    def test_no_parse_pure_capacity(self, monkeypatch):
        # pm-i4 title gives (None, None) — falls to pure capacity
        self._today(monkeypatch, date(2026, 5, 22))
        proposed, note = rm.compute_layered_due_date(None, None, 5.0, [])
        # 5.0 / 5.0 * 7 = 7 days from today
        assert proposed == date(2026, 5, 29)
        assert note == ""

    def test_closed_prior_pure_capacity(self, monkeypatch):
        self._today(monkeypatch, date(2026, 5, 22))
        fixture = [
            {"number": 100, "title": "i1 - S3 - V1", "state": "closed", "due_on": "2026-03-01T00:00:00Z"},
        ]
        proposed, note = rm.compute_layered_due_date(2, 3, 5.0, fixture)
        assert proposed == date(2026, 5, 29)
        assert "M#100 closed" in note

    def test_undated_prior_skips_sequencing(self, monkeypatch):
        self._today(monkeypatch, date(2026, 5, 22))
        fixture = [
            {"number": 100, "title": "i2 - S4 - V1", "state": "open", "due_on": None},
        ]
        proposed, note = rm.compute_layered_due_date(3, 4, 2.5, fixture)
        # 2.5 / 5.0 * 7 = 3.5 → 4 days (rounded)
        assert proposed == date(2026, 5, 26)
        assert "M#100 undated" in note
        assert "sequencing skipped" in note

    def test_no_prior_pure_capacity(self, monkeypatch):
        self._today(monkeypatch, date(2026, 5, 22))
        # No S3 milestone exists in fixture
        fixture = [
            {"number": 200, "title": "i2 - S5 - Modeling", "state": "open", "due_on": "2026-05-31T00:00:00Z"},
        ]
        proposed, note = rm.compute_layered_due_date(2, 3, 8.5, fixture)
        # 8.5 / 5.0 * 7 = 11.9 → 12 days
        assert proposed == date(2026, 6, 3)
        assert "no prior" in note

    def test_normal_stack_after_open_prior(self, monkeypatch):
        self._today(monkeypatch, date(2026, 5, 22))
        fixture = [
            {"number": 100, "title": "i2 - S3 - V1", "state": "open", "due_on": "2026-06-03T00:00:00Z"},
        ]
        proposed, note = rm.compute_layered_due_date(3, 3, 1.0, fixture)
        # base = max(2026-06-03, 2026-05-22) = 2026-06-03
        # 1.0 / 5.0 * 7 = 1.4 → 1 day (rounded)
        assert proposed == date(2026, 6, 4)
        assert "stack after M#100" in note
        assert "2026-06-03" in note

    def test_overdue_prior_uses_today(self, monkeypatch):
        self._today(monkeypatch, date(2026, 5, 22))
        fixture = [
            {"number": 100, "title": "i2 - S3 - V1", "state": "open", "due_on": "2026-04-01T00:00:00Z"},
        ]
        proposed, note = rm.compute_layered_due_date(3, 3, 5.0, fixture)
        # Prior is overdue; max(today, overdue) = today
        # 5.0 / 5.0 * 7 = 7 days
        assert proposed == date(2026, 5, 29)
        assert "stack after M#100" in note

    def test_s7_paired_gates_on_s5(self, monkeypatch):
        self._today(monkeypatch, date(2026, 5, 22))
        fixture = [
            {"number": 200, "title": "i2 - S5 - Modeling - TCR Panel", "state": "open", "due_on": "2026-05-31T00:00:00Z"},
            {"number": 201, "title": "i2 - S7 - Publication - TCR Panel", "state": "open", "due_on": "2026-06-04T00:00:00Z"},
        ]
        proposed, note = rm.compute_layered_due_date(2, 7, 2.5, fixture)
        # base = max(2026-05-31, today) = 2026-05-31; +4 days (2.5/5*7 = 3.5 → 4)
        assert proposed == date(2026, 6, 4)
        assert "paired-S7" in note
        assert "M#200" in note

    def test_s7_standalone_no_paired_s5(self, monkeypatch):
        self._today(monkeypatch, date(2026, 5, 22))
        # i3-S7 Lit Review — no paired i3-S5
        fixture = [
            {"number": 300, "title": "i3 - S7 - Publication - Lit Review", "state": "open", "due_on": "2026-05-31T00:00:00Z"},
        ]
        proposed, note = rm.compute_layered_due_date(3, 7, 2.0, fixture)
        # Pure capacity: 2.0 / 5.0 * 7 = 2.8 → 3 days
        assert proposed == date(2026, 5, 25)
        assert "standalone S7" in note


class TestLiveIntegrationSmoke:
    """Live API smoke tests. Skipped by default; opt-in via pytest -m live."""

    @pytest.mark.live
    def test_seven_known_sequence_bound_milestones_show_no_change(self):
        """The 7 sequence-bound milestones from 2026-05-22 rate cascade should
        no longer flag UPDATE NEEDED after the sequencing-aware fix.
        """
        import subprocess
        expected_no_change = [10, 11, 13, 15, 18, 24, 30]
        failures = []
        for ms in expected_no_change:
            result = subprocess.run(
                ["python3", "scripts/pm/recheck_milestone.py", "--milestone", str(ms)],
                capture_output=True, text=True,
            )
            if result.returncode != 0:
                failures.append(f"M#{ms} exit={result.returncode}:\n{result.stdout}")
        assert not failures, "Sequence-bound milestones still flagging:\n" + "\n---\n".join(failures)

    @pytest.mark.live
    def test_eight_capacity_bound_milestones_still_no_change(self):
        """Regression check: the 8 capacity-bound milestones from the same
        cascade should still show [No change] after the fix.
        """
        import subprocess
        expected_no_change = [3, 5, 17, 18, 20, 21, 22, 26]
        # Note: #18 appears in both sets (capacity-bound AND prior in chain)
        failures = []
        for ms in expected_no_change:
            result = subprocess.run(
                ["python3", "scripts/pm/recheck_milestone.py", "--milestone", str(ms)],
                capture_output=True, text=True,
            )
            if result.returncode != 0:
                failures.append(f"M#{ms} exit={result.returncode}:\n{result.stdout}")
        assert not failures, "Capacity-bound milestones regressed:\n" + "\n---\n".join(failures)
