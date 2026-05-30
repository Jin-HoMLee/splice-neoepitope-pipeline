"""Unit tests for scripts/pm/recheck_milestone.py."""

import sys
from datetime import date, timedelta
from pathlib import Path

import pytest

# Make the script importable as a module
SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts" / "pm"
sys.path.insert(0, str(SCRIPT_DIR))

import recheck_milestone as rm


# --- Property-based live-smoke helpers (Issue #506) -------------------------
# The live smoke test validates that the recheck integration *works* on every
# open milestone, rather than freezing which specific milestones happen to show
# [No change] today. A milestone legitimately drifting to [UPDATE NEEDED] as the
# calendar advances is a correct recommendation, not a test failure — the old
# hardcoded baselines went red on exactly that drift.
#
# Exit codes (see recheck_milestone.py module docstring):
#   0 = [No change]      2 = [UPDATE NEEDED] / [UNSIZED]      1 = error/crash
WELL_FORMED_EXIT_CODES = (0, 2)
WELL_FORMED_STATUSES = ("[No change]", "[UPDATE NEEDED]", "[UNSIZED]")


def is_well_formed_recheck(returncode: int, stdout: str) -> bool:
    """True if a ``recheck_milestone.py --milestone N`` run completed its
    analysis and emitted a well-formed report.

    "Completed" means exit 0 (No change) or 2 (UPDATE NEEDED / UNSIZED) — both
    are legitimate recommendations the script is *supposed* to produce. Exit 1
    (error/crash) is NOT well-formed. "Well-formed" means the report carries its
    header (``Milestone:``) and a known ``Status:`` marker. This is the property
    the live smoke test asserts, in place of a frozen which-milestones-are-stable
    baseline list.
    """
    if returncode not in WELL_FORMED_EXIT_CODES:
        return False
    if "Milestone:" not in stdout:
        return False
    return any(status in stdout for status in WELL_FORMED_STATUSES)


def open_milestone_numbers() -> list[int]:
    """Live: all open milestone numbers, via the script's own ``gh()`` helper."""
    data = rm.gh(
        "api", "--paginate",
        f"repos/{rm.REPO}/milestones?state=open&per_page=100",
    )
    return [m["number"] for m in data]


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
    """Live API smoke test. Skipped by default; opt-in via ``pytest -m live``.

    Property-based (Issue #506). Asserts the recheck integration produces a
    well-formed recommendation for *every* open milestone, instead of pinning a
    hardcoded baseline of which milestones show [No change].

    Why property-based beats a baseline list: as the calendar advances,
    milestones legitimately flip between [No change] and [UPDATE NEEDED] (a
    re-date recommendation). A frozen baseline turns red on that normal drift —
    the old ``[10, 11, 13, 15, 18, 24, 30]`` and ``[3, 5, 17, 18, 20, 21, 22, 26]``
    lists did exactly this, producing recurring false-positive CI failures that
    gate nothing (the marker is ``live`` and ``ci-tools-pytest`` is not a
    required check). The due-date *logic* is fully covered by the unit-level
    ``TestComputeLayeredDueDate`` branch cases; this smoke test's only job is to
    validate the live *integration* — that the script runs end-to-end against
    the real board and emits a parseable report for whatever state it finds.
    """

    @pytest.mark.live
    def test_recheck_emits_well_formed_report_for_every_open_milestone(self):
        import subprocess

        milestones = open_milestone_numbers()
        assert milestones, "expected at least one open milestone on the board"

        malformed = []
        for ms in milestones:
            try:
                result = subprocess.run(
                    ["python3", "scripts/pm/recheck_milestone.py", "--milestone", str(ms)],
                    capture_output=True, text=True, timeout=60,
                )
            except subprocess.TimeoutExpired:
                # Bounds the blast radius now that the loop hits every open
                # milestone: a hung gh/network call fails one milestone cleanly
                # instead of hanging the whole suite.
                malformed.append(f"M#{ms} timed out after 60s (hung gh/network call)")
                continue
            if not is_well_formed_recheck(result.returncode, result.stdout):
                malformed.append(
                    f"M#{ms} exit={result.returncode}:\n{result.stdout}\n{result.stderr}"
                )
        assert not malformed, (
            "recheck produced a malformed report (exit 1 / missing header or "
            "Status line) — a crash or output-contract break, NOT a routine "
            "[UPDATE NEEDED] drift:\n" + "\n---\n".join(malformed)
        )


class TestRecheckOutputProperty:
    """Unit coverage proving the property-based smoke check above is not vacuous.

    It must ACCEPT both a stable [No change] run and a drifted [UPDATE NEEDED]
    run (two distinct milestone-state snapshots — Issue #506 AC), plus the
    [UNSIZED] case, and REJECT a crash or a report missing its Status line. No
    live API access required, so this runs in every ``ci-tools-pytest`` sweep.
    """

    NO_CHANGE_STDOUT = (
        "Milestone: i4 - S5 - Modeling - Google Batch\n"
        "Current due_on: 2026-06-20\n"
        "Open issues (1):\n  - #999 (M, ~2.5d)\n"
        "Remaining capacity: 2.5d\n"
        "Proposed due_on: 2026-06-22 (delta +2 days) (stack after M#100 close 2026-06-20)\n"
        "Status: [No change]\n"
    )
    UPDATE_NEEDED_STDOUT = (
        "Milestone: i3 - S5 - Modeling - HLA Panel\n"
        "Current due_on: 2026-05-10\n"
        "Open issues (3):\n  - #1 (L, ~3.5d)\n  - #2 (M, ~2.5d)\n  - #3 (M, ~2.5d)\n"
        "Remaining capacity: 8.5d\n"
        "Proposed due_on: 2026-06-15 (delta +36 days) (no prior same-S milestone)\n"
        "Status: [UPDATE NEEDED]\n"
    )
    UNSIZED_STDOUT = (
        "Milestone: pm-i4 - PM Tooling, Memory & Methodology\n"
        "Current due_on: 2026-06-03\n"
        "Open issues (2):\n  - #10 (no size, 0d)\n  - #11 (no size, 0d)\n"
        "Remaining capacity: 0.0d\n"
        "Proposed due_on: (cannot compute — 2 open issue(s) missing Size)\n"
        "Status: [UNSIZED] — assign Size on the project board, then re-run\n"
    )

    def test_stable_no_change_snapshot_is_well_formed(self):
        assert is_well_formed_recheck(0, self.NO_CHANGE_STDOUT)

    def test_drifted_update_needed_snapshot_is_well_formed(self):
        # exit 2 (UPDATE NEEDED) is a valid recommendation, not a failure.
        assert is_well_formed_recheck(2, self.UPDATE_NEEDED_STDOUT)

    def test_unsized_snapshot_is_well_formed(self):
        assert is_well_formed_recheck(2, self.UNSIZED_STDOUT)

    def test_crash_is_not_well_formed(self):
        assert not is_well_formed_recheck(
            1, "Traceback (most recent call last):\n  ...\nKeyError: 'milestone'\n"
        )

    def test_missing_status_line_is_not_well_formed(self):
        assert not is_well_formed_recheck(
            0, "Milestone: i4 - S5 - Modeling\nRemaining capacity: 1.0d\n"
        )

    def test_unexpected_exit_code_is_not_well_formed(self):
        # Even with a plausible-looking report, an exit code outside {0, 2}
        # signals the script did not complete its contract (e.g. SIGKILL=137).
        assert not is_well_formed_recheck(137, self.NO_CHANGE_STDOUT)
