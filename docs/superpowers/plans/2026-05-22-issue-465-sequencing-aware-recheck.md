# Sequencing-aware milestone recheck Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend `scripts/pm/recheck_milestone.py` so its proposed `due_on` dates respect same-S-stage sequencing and paired-S7 gating — eliminating false `[UPDATE NEEDED]` flags for the 7 sequence-bound milestones from the 2026-05-22 rate-change cascade.

**Architecture:** 3 new pure helper functions + 1 layered-compute function added to the existing single-file script. Integration point is the inline date calculation in `compute_recheck` (current line 120). Tests via pytest at `tools/ci/test_recheck_milestone.py` mirroring the sibling `test_recheck_parent_status.py` convention.

**Tech Stack:** Python 3 stdlib only (`re`, `datetime`, `subprocess`, `json`, `argparse`). Pytest for tests. `gh` CLI via subprocess for live API calls.

**Spec:** [docs/superpowers/specs/2026-05-22-issue-465-sequencing-aware-recheck-design.md](../specs/2026-05-22-issue-465-sequencing-aware-recheck-design.md)

---

## File Structure

**Modify:** `scripts/pm/recheck_milestone.py` — add `import re`, 3 helpers, 1 layered-compute function (~80 LOC added). Refactor lines 86-120 of `compute_recheck` to fetch `all_milestones` once and call the new layered function.

**Create:** `tools/ci/test_recheck_milestone.py` — pytest tests mirroring `tools/ci/test_recheck_parent_status.py`. Unit tests (no network, default-on) + 1 live integration smoke test (marked `live`, opt-in via `pytest -m live`).

**Modify:** `research/lab_notebook/pm.md` — append new time section under existing 2026-05-22 date block capturing the verification run.

**Out-of-repo (handled after PR merge, not part of this plan):**
- `~/.claude/projects/.../memory/feedback_milestones.md` — add cross-reference pointing to the script for sequencing math.

---

## Task 1: Add `re` import + `parse_milestone_title` helper (TDD)

**Files:**
- Create: `tools/ci/test_recheck_milestone.py`
- Modify: `scripts/pm/recheck_milestone.py:15-19` (imports block)

- [ ] **Step 1: Create the test file with parse-title tests**

Create `tools/ci/test_recheck_milestone.py`:

```python
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
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
workflow/tests/.venv/bin/python -m pytest tools/ci/test_recheck_milestone.py -v
```

Expected: FAIL with `AttributeError: module 'recheck_milestone' has no attribute 'parse_milestone_title'`.

- [ ] **Step 3: Add `re` import + `parse_milestone_title` to the script**

Modify `scripts/pm/recheck_milestone.py:15-19` (the imports block). Replace:

```python
import argparse
import json
import subprocess
import sys
from datetime import date, datetime, timedelta
```

with:

```python
import argparse
import json
import re
import subprocess
import sys
from datetime import date, datetime, timedelta
```

Then add the helper function immediately after line 27 (`SIZE_WEIGHTS = ...`):

```python


def parse_milestone_title(title: str) -> tuple[int, int] | None:
    """Parse 'i<N> - S<M> - ...' titles. Returns (iteration, stage) or None.

    Role-meta (pm-i*, dev-i*) and legacy (M1, M2) titles return None and fall
    through to pure-capacity behavior.
    """
    m = re.match(r"^i(\d+)\s*-\s*S(\d+)\s*-\s*", title)
    return (int(m.group(1)), int(m.group(2))) if m else None
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
workflow/tests/.venv/bin/python -m pytest tools/ci/test_recheck_milestone.py -v
```

Expected: 6 passed.

- [ ] **Step 5: Commit**

```bash
git add tools/ci/test_recheck_milestone.py scripts/pm/recheck_milestone.py
git commit -m "feat(pm): parse_milestone_title helper for sequencing-aware recheck

Pure regex parse of 'i<N> - S<M> - ...' shape. Returns None for
role-meta and legacy titles so they fall through to pure capacity.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 2: Add `find_prior_same_stage` helper (TDD)

**Files:**
- Modify: `tools/ci/test_recheck_milestone.py` (add test class)
- Modify: `scripts/pm/recheck_milestone.py` (add helper after `parse_milestone_title`)

- [ ] **Step 1: Add test class to test file**

Append to `tools/ci/test_recheck_milestone.py`:

```python


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
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
workflow/tests/.venv/bin/python -m pytest tools/ci/test_recheck_milestone.py::TestFindPriorSameStage -v
```

Expected: FAIL with `AttributeError: module 'recheck_milestone' has no attribute 'find_prior_same_stage'`.

- [ ] **Step 3: Implement `find_prior_same_stage`**

Add to `scripts/pm/recheck_milestone.py` immediately after `parse_milestone_title`:

```python


def find_prior_same_stage(
    iteration: int,
    stage: int,
    all_milestones: list[dict],
) -> dict | None:
    """Highest-iteration prior in the same S-stage chain. Includes closed."""
    candidates = []
    for ms in all_milestones:
        parsed = parse_milestone_title(ms["title"])
        if parsed is None:
            continue
        n, s = parsed
        if s == stage and n < iteration:
            candidates.append((n, ms))
    if not candidates:
        return None
    candidates.sort(key=lambda x: x[0])
    return candidates[-1][1]
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
workflow/tests/.venv/bin/python -m pytest tools/ci/test_recheck_milestone.py::TestFindPriorSameStage -v
```

Expected: 6 passed.

- [ ] **Step 5: Commit**

```bash
git add tools/ci/test_recheck_milestone.py scripts/pm/recheck_milestone.py
git commit -m "feat(pm): find_prior_same_stage helper

Same-S-stage chain lookup. Picks highest-iteration prior (i5-S3 returns
i4-S3, not i1). Includes closed milestones so first-iteration chains
can still anchor.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 3: Add `find_open_same_iteration_S5` helper (TDD)

**Files:**
- Modify: `tools/ci/test_recheck_milestone.py`
- Modify: `scripts/pm/recheck_milestone.py`

- [ ] **Step 1: Add test class**

Append to `tools/ci/test_recheck_milestone.py`:

```python


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
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
workflow/tests/.venv/bin/python -m pytest tools/ci/test_recheck_milestone.py::TestFindOpenSameIterationS5 -v
```

Expected: FAIL with `AttributeError`.

- [ ] **Step 3: Implement `find_open_same_iteration_S5`**

Add to `scripts/pm/recheck_milestone.py` immediately after `find_prior_same_stage`:

```python


def find_open_same_iteration_S5(
    iteration: int,
    all_milestones: list[dict],
) -> dict | None:
    """Find an OPEN i<N> - S5 - ... milestone. Used for paired-S7 gating.

    Loose match: same iteration number is enough. Arc-mismatch (e.g. i4-S7
    'TCR-pMHC Landscape' vs i4-S5 'Google Batch') is a separate data-hygiene
    concern not solvable here.
    """
    for ms in all_milestones:
        parsed = parse_milestone_title(ms["title"])
        if parsed is None:
            continue
        n, s = parsed
        if n == iteration and s == 5 and ms["state"] == "open":
            return ms
    return None
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
workflow/tests/.venv/bin/python -m pytest tools/ci/test_recheck_milestone.py::TestFindOpenSameIterationS5 -v
```

Expected: 4 passed.

- [ ] **Step 5: Commit**

```bash
git add tools/ci/test_recheck_milestone.py scripts/pm/recheck_milestone.py
git commit -m "feat(pm): find_open_same_iteration_S5 helper for paired-S7 gating

Loose match by iteration number only. Arc-mismatch tolerated (data
hygiene concern, separate).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 4: Add `compute_layered_due_date` core function (TDD)

**Files:**
- Modify: `tools/ci/test_recheck_milestone.py`
- Modify: `scripts/pm/recheck_milestone.py`

- [ ] **Step 1: Add comprehensive test class covering all 7 branches**

Append to `tools/ci/test_recheck_milestone.py`:

```python


class TestComputeLayeredDueDate:
    """Covers all 7 logic branches of compute_layered_due_date."""

    # Reference date for deterministic test expectations.
    # date.today() varies; we monkeypatch it in tests via freeze fixture.

    def _today(self, monkeypatch, fake_today: date):
        import recheck_milestone as rm_inner
        class _FakeDate(date):
            @classmethod
            def today(cls):
                return fake_today
        monkeypatch.setattr(rm_inner, "date", _FakeDate)

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
        # base = max(2026-05-31, today) = 2026-05-31; +4 days
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
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
workflow/tests/.venv/bin/python -m pytest tools/ci/test_recheck_milestone.py::TestComputeLayeredDueDate -v
```

Expected: FAIL with `AttributeError: module 'recheck_milestone' has no attribute 'compute_layered_due_date'`.

- [ ] **Step 3: Implement `compute_layered_due_date`**

Add to `scripts/pm/recheck_milestone.py` immediately after `find_open_same_iteration_S5`:

```python


def compute_layered_due_date(
    iteration: int | None,
    stage: int | None,
    capacity_days: float,
    all_milestones: list[dict],
) -> tuple[date, str]:
    """Return (proposed_due, reasoning_note).

    Handles 4 main branches:
      1. No title parse (role-meta etc.) -> pure capacity
      2. S7 paired with open S5 in same iteration -> stack on S5 close
      3. S7 standalone (no paired S5) -> pure capacity
      4. Same-S-stage stacking (closed/undated/normal/overdue prior cases)
    """
    today = date.today()

    if iteration is None or stage is None:
        # Non-S-stage milestone (pm-i*, dev-i*, M1, etc.) — pure capacity
        base = today
        note = ""
    elif stage == 7:
        paired = find_open_same_iteration_S5(iteration, all_milestones)
        if paired and paired.get("due_on"):
            paired_date = date.fromisoformat(paired["due_on"][:10])
            base = max(paired_date, today)
            note = f"(paired-S7: unblocks at M#{paired['number']} close {paired['due_on'][:10]})"
        else:
            # Standalone S7 (e.g. Lit Review i3-S7) — pure capacity
            base = today
            note = "(standalone S7 — no paired open S5)"
    else:
        prior = find_prior_same_stage(iteration, stage, all_milestones)
        if prior is None:
            base = today
            note = "(no prior same-S milestone)"
        elif prior["state"] == "closed":
            base = today
            note = f"(prior M#{prior['number']} closed)"
        elif prior.get("due_on") is None:
            base = today
            note = f"(prior M#{prior['number']} undated — sequencing skipped)"
        else:
            prior_date = date.fromisoformat(prior["due_on"][:10])
            base = max(prior_date, today)
            note = f"(stack after M#{prior['number']} close {prior['due_on'][:10]})"

    calendar_days = int(round(capacity_days / AVAILABILITY_RATE * 7))
    proposed = base + timedelta(days=calendar_days)
    return (proposed, note)
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
workflow/tests/.venv/bin/python -m pytest tools/ci/test_recheck_milestone.py::TestComputeLayeredDueDate -v
```

Expected: 9 passed.

- [ ] **Step 5: Run the full unit test file to verify no regressions**

```bash
workflow/tests/.venv/bin/python -m pytest tools/ci/test_recheck_milestone.py -v
```

Expected: 25 passed (6 parse + 6 prior + 4 paired-S5 + 9 layered-compute).

- [ ] **Step 6: Commit**

```bash
git add tools/ci/test_recheck_milestone.py scripts/pm/recheck_milestone.py
git commit -m "feat(pm): compute_layered_due_date — sequencing-aware proposed dates

Covers 7 branches: no-parse, closed-prior, undated-prior, no-prior,
normal-stack, overdue-prior (today guard), S7-paired, S7-standalone.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 5: Integrate `compute_layered_due_date` into `compute_recheck`

**Files:**
- Modify: `scripts/pm/recheck_milestone.py:86-129` (`compute_recheck` body)

- [ ] **Step 1: Read current `compute_recheck` to confirm line numbers**

```bash
sed -n '85,130p' scripts/pm/recheck_milestone.py
```

Expected: function `compute_recheck(milestone_number: int) -> int` body matching the spec.

- [ ] **Step 2: Replace the inline calculation with the layered call**

In `scripts/pm/recheck_milestone.py`, modify `compute_recheck` (currently lines 85-129). The full replacement function:

```python
def compute_recheck(milestone_number: int) -> int:
    meta = milestone_meta(milestone_number)
    title = meta["title"]
    current_due_raw = meta.get("due_on")
    current_due_date = (
        datetime.fromisoformat(current_due_raw.replace("Z", "+00:00")).date()
        if current_due_raw else None
    )

    issue_numbers = open_issues_in_milestone(title)
    sizes = sizes_for_issues(issue_numbers)

    print(f"Milestone: {title}")
    print(f"Current due_on: {current_due_date or '—'}")
    print(f"Open issues ({len(issue_numbers)}):")
    remaining = 0.0
    for n in sorted(issue_numbers):
        size = sizes.get(n)
        weight = SIZE_WEIGHTS.get(size or "", 0)
        remaining += weight
        size_disp = f"{size}, ~{weight}d" if size else "no size, 0d"
        print(f"  - #{n} ({size_disp})")
    print(f"Remaining capacity: {remaining}d")

    if remaining == 0:
        unsized_count = sum(1 for n in issue_numbers if sizes.get(n) is None)
        if unsized_count > 0:
            print(f"Proposed due_on: (cannot compute — {unsized_count} open issue(s) missing Size)")
            print("Status: [UNSIZED] — assign Size on the project board, then re-run")
            return 2
        print("Proposed due_on: (no open work)")
        print("Status: [No change] — milestone has no remaining capacity")
        return 0

    # Sequencing-aware: fetch all milestones once, derive (iteration, stage), apply layered logic
    all_milestones_raw = gh("api", f"repos/{REPO}/milestones?state=all&per_page=100")
    parsed = parse_milestone_title(title)
    iteration, stage = parsed if parsed else (None, None)
    proposed_due, note = compute_layered_due_date(iteration, stage, remaining, all_milestones_raw)

    delta = (proposed_due - current_due_date).days if current_due_date else None
    delta_str = f"{delta:+d}" if delta is not None else "n/a"
    note_suffix = f" {note}" if note else ""
    print(f"Proposed due_on: {proposed_due} (delta {delta_str} days){note_suffix}")

    if delta is None or abs(delta) <= THRESHOLD_DAYS:
        print("Status: [No change]")
        return 0
    print("Status: [UPDATE NEEDED]")
    return 2
```

- [ ] **Step 3: Smoke-test against one capacity-bound milestone (no regression)**

```bash
python3 scripts/pm/recheck_milestone.py --milestone 3
```

Expected output (M#3 i2-S3 GTEx is capacity-bound, no prior open S3):

```
Milestone: i2 - S3 - Data Preparation - GTEx Pan-Tissue Junction Filter
Current due_on: 2026-06-03
Open issues (3):
  - #126 (L, ~3.5d)
  - #211 (M, ~2.5d)
  - #212 (M, ~2.5d)
Remaining capacity: 8.5d
Proposed due_on: 2026-06-03 (delta +0 days) (no prior same-S milestone)
Status: [No change]
```

The new reasoning note `(no prior same-S milestone)` is the visible change. `[No change]` status unchanged.

- [ ] **Step 4: Smoke-test against one sequence-bound milestone (the fix)**

```bash
python3 scripts/pm/recheck_milestone.py --milestone 11
```

Expected:

```
Milestone: i3 - S3 - Data Preparation - Aligner & Input Format Improvements
Current due_on: 2026-06-05
Open issues (1):
  - #297 (S, ~1.0d)
Remaining capacity: 1.0d
Proposed due_on: 2026-06-04 (delta -1 days) (stack after M#3 close 2026-06-03)
Status: [No change]
```

The previously-flagged `[UPDATE NEEDED]` now reads `[No change]` because the layered date (2026-06-04) matches the stored due_on (2026-06-05) within the 7-day threshold.

- [ ] **Step 5: Commit**

```bash
git add scripts/pm/recheck_milestone.py
git commit -m "feat(pm): integrate layered date into compute_recheck

Fetches all milestones once, derives (iteration, stage), delegates
proposed_due computation to compute_layered_due_date. Reasoning note
appended to Proposed due_on output line.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 6: Add live integration smoke test

**Files:**
- Modify: `tools/ci/test_recheck_milestone.py` (add live class + pytest marker)
- Create: `tools/ci/pytest.ini` (if not already present — register `live` marker)

- [ ] **Step 1: Check if pytest.ini exists in tools/ci/**

```bash
ls tools/ci/pytest.ini 2>&1 || echo "not present"
```

- [ ] **Step 2: Register `live` marker (skip if pytest.ini already exists with markers)**

If `tools/ci/pytest.ini` does NOT exist, create it:

```ini
[pytest]
markers =
    live: tests requiring live gh API access (opt-in via -m live)
```

If it exists, ensure it has the `live` marker line under `markers =`.

- [ ] **Step 3: Append live integration smoke test**

Append to `tools/ci/test_recheck_milestone.py`:

```python


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
    def test_nine_capacity_bound_milestones_still_no_change(self):
        """Regression check: the 9 capacity-bound milestones from the same
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
```

Add the `import pytest` at the top of the file if not already present:

```python
import pytest
import sys
from datetime import date, timedelta
from pathlib import Path
```

- [ ] **Step 4: Run the live tests manually**

```bash
workflow/tests/.venv/bin/python -m pytest tools/ci/test_recheck_milestone.py -m live -v
```

Expected: 2 passed (both live smoke tests).

If any milestone flags `[UPDATE NEEDED]`, the failure output includes its full recheck output for debugging. Most likely cause: a prior milestone's `due_on` shifted since the spec was written.

- [ ] **Step 5: Run the full test suite (unit tests, live skipped)**

```bash
workflow/tests/.venv/bin/python -m pytest tools/ci/test_recheck_milestone.py -v
```

Expected: 25 unit tests pass, 2 live tests skipped (no `-m live`).

- [ ] **Step 6: Commit**

```bash
git add tools/ci/test_recheck_milestone.py
test -f tools/ci/pytest.ini && git add tools/ci/pytest.ini
git commit -m "test(pm): live integration smoke for sequencing-aware recheck

Verifies the 7 known sequence-bound milestones from 2026-05-22 cascade
now show [No change], and the 9 capacity-bound milestones still do.
Opt-in via pytest -m live.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 7: Add PM lab notebook entry

**Files:**
- Modify: `research/lab_notebook/pm.md` (insert new time section at TOP of 2026-05-22 date block)

- [ ] **Step 1: Get current UTC timestamp**

```bash
date -u +'%Y-%m-%d %H:%M UTC'
```

Example output: `2026-05-22 14:30 UTC`. Use this timestamp in the entry.

- [ ] **Step 2: Insert the new time section ABOVE the existing 13:04 UTC entry**

Locate the line `### 13:04 UTC — Editor: PM` in `research/lab_notebook/pm.md`. Insert immediately ABOVE it:

```markdown
### HH:MM UTC — Editor: PM

#### [Issue #465](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/465) — sequencing-aware milestone recheck ship

**Trigger.** Today's 13:04 UTC rate-change cascade created 7 false-flag `[UPDATE NEEDED]` milestones (sequence-bound, capacity-formula too early). [Issue #465](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/465) filed same session as the rung-3 mechanism fix per [[mechanism-over-memory]] — the sequencing rule lived in `feedback_milestones.md` prose only; moving it to the script silences the noise deterministically.

**Implementation.** Single-PR ship: spec → 3 helpers + 1 layered-compute function + integration → pytest unit tests + live integration smoke → this entry. ~80 LOC added to `scripts/pm/recheck_milestone.py`; ~150 LOC test file at `tools/ci/test_recheck_milestone.py`.

**Design choices** (from spec):
- **Single-level prior lookup** over recursive proposed-close — trusts GitHub's stored `due_on` as source of truth; hook's cascade-on-activity property converges naturally
- **Loose paired-S7 match** (same iteration, any arc) — arc-mismatch is a separate data-hygiene concern (e.g. [M#28](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/milestone/28) i4-S7 'TCR-pMHC Landscape' vs [M#13](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/milestone/13) i4-S5 'Google Batch')
- **Role-meta-axis explicitly skipped** — pm-i*/dev-i* run partly in parallel; strict stacking would create its own false flags
- **Report-only preserved** — operator runs PATCH manually per the existing script's design

**Verification.** Live integration smoke: all 7 sequence-bound milestones from morning's cascade now `[No change]`; 8 capacity-bound milestones unchanged (regression check). Unit tests cover 7 logic branches + edge cases (closed prior, undated prior, no prior, normal stack, overdue prior with today guard, S7-paired, S7-standalone, non-S-stage parse).

**Follow-ups.**
- **Memory update** (out-of-repo): point `feedback_milestones.md` at `scripts/pm/recheck_milestone.py` for the sequencing math instead of prose-only description
- **Arc-mismatch data hygiene**: i4-S7 ↔ i4-S5 arc mismatch is real; worth a future review to either rename milestones for arc-consistency or formalize the cross-iteration pairing exception ([M#28](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/milestone/28) probably should pair with i5-S5)
```

Replace `HH:MM` with the timestamp from Step 1.

- [ ] **Step 3: Commit**

```bash
git add research/lab_notebook/pm.md
git commit -m "docs(pm): lab notebook entry for Issue #465 ship

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 8: Tick Issue #465 acceptance criteria + open PR

**Files:**
- Modify: Issue #465 body (via `gh issue edit`)
- Open: PR via `gh pr create`

- [ ] **Step 1: Update Issue #465 ACs (tick the 6 acceptance items)**

```bash
gh issue view 465 --json body -q '.body' > /tmp/issue465_body.md
# Tick all 6 AC checkboxes
sed -i '' 's|- \[ \] All 7 sequence-bound|- [x] All 7 sequence-bound|' /tmp/issue465_body.md
sed -i '' 's|- \[ \] Capacity-bound|- [x] Capacity-bound|' /tmp/issue465_body.md
sed -i '' 's|- \[ \] Recompute paired-S7|- [x] Recompute paired-S7|' /tmp/issue465_body.md
sed -i '' 's|- \[ \] Lab notebook entry|- [x] Lab notebook entry|' /tmp/issue465_body.md
gh issue edit 465 --body-file /tmp/issue465_body.md
```

Also update the Tasks section (7 task checkboxes). Note: the task `- [ ] Update feedback_milestones.md` stays unticked since that's out-of-repo work post-merge. Tick the other 6:

```bash
sed -i '' 's|- \[ \] Parse `(iteration, stage)`|- [x] Parse `(iteration, stage)`|' /tmp/issue465_body.md
sed -i '' 's|- \[ \] Find same-S-stage prior|- [x] Find same-S-stage prior|' /tmp/issue465_body.md
sed -i '' 's|- \[ \] Compute layered|- [x] Compute layered|' /tmp/issue465_body.md
sed -i '' 's|- \[ \] Add S7-paired|- [x] Add S7-paired|' /tmp/issue465_body.md
sed -i '' 's|- \[ \] Update `\[UPDATE NEEDED\]`|- [x] Update `[UPDATE NEEDED]`|' /tmp/issue465_body.md
sed -i '' 's|- \[ \] Test against current state|- [x] Test against current state|' /tmp/issue465_body.md
gh issue edit 465 --body-file /tmp/issue465_body.md
```

- [ ] **Step 2: Push branch**

```bash
git push -u origin feat/pm/issue-465-sequencing-aware-recheck
```

Expected: branch pushed cleanly, new branch URL printed.

- [ ] **Step 3: Open PR**

```bash
gh pr create --title "feat(pm): sequencing-aware milestone recheck (Issue #465)" --body "$(cat <<'EOF'
**Created by:** PM

## Summary

Extends `scripts/pm/recheck_milestone.py` with same-S-stage sequencing and paired-S7 gating. Eliminates the 7 false-flag `[UPDATE NEEDED]` cases from today's 2026-05-22 rate-change cascade.

Closes [Issue #465](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/465).

## Design

[docs/superpowers/specs/2026-05-22-issue-465-sequencing-aware-recheck-design.md](docs/superpowers/specs/2026-05-22-issue-465-sequencing-aware-recheck-design.md)

## What's in this PR

- **`scripts/pm/recheck_milestone.py`** — 3 new helpers (`parse_milestone_title`, `find_prior_same_stage`, `find_open_same_iteration_S5`) + 1 layered-compute function (`compute_layered_due_date`) + integration into `compute_recheck`. ~80 LOC added.
- **`tools/ci/test_recheck_milestone.py`** — pytest unit tests (25 tests across 4 classes) + 2 live integration smoke tests (opt-in via `pytest -m live`). Mirrors `tools/ci/test_recheck_parent_status.py` convention.
- **`tools/ci/pytest.ini`** — register `live` marker (if not already present).
- **`research/lab_notebook/pm.md`** — 2026-05-22 entry capturing the implementation + verification.

## What's NOT in this PR

- **Memory update** (`feedback_milestones.md` cross-reference to the script) — out-of-repo; handled post-merge via direct memory edit.

## Test plan

- [ ] `workflow/tests/.venv/bin/python -m pytest tools/ci/test_recheck_milestone.py -v` — 25 unit tests pass
- [ ] `workflow/tests/.venv/bin/python -m pytest tools/ci/test_recheck_milestone.py -m live -v` — 2 live smoke tests pass
- [ ] Manual: `python3 scripts/pm/recheck_milestone.py --milestone 11` shows `[No change]` (was `[UPDATE NEEDED]` before)
- [ ] CI green (`pipeline-pytest`, `pipeline-snakemake-dry-run`, `ci-tools-pytest`)

## Verification of impact

After this PR, the 7 sequence-bound milestones (M#10, #11, #13, #15, #18, #24, #30) all report `[No change]` on hook fires, eliminating recurring false-flag noise. The 9 capacity-bound milestones (M#3, #5, #17, #20, #21, #22, #26 + #18 in both sets) still flag correctly on real capacity changes.
EOF
)"
```

- [ ] **Step 4: Confirm PR opened**

PR URL printed. Note the number for the merge step.

- [ ] **Step 5: Wait for CI + run audit_and_merge**

```bash
# Wait for CI (~2-3 min typically)
gh pr checks <PR_NUMBER> --watch
# Once green:
bash scripts/audit_and_merge.sh <PR_NUMBER>
```

Expected: `✓ PR #N merged (X test-plan boxes ticked, Y AC boxes across 1 linked issues ticked)`.

If audit fails on unticked Test plan boxes, tick them via `gh pr edit <N> --body-file` (mirror the PR #466 procedure from today's session).

- [ ] **Step 6: Sync main + delete local branch**

```bash
git checkout main && git pull --ff-only && git branch -D feat/pm/issue-465-sequencing-aware-recheck
```

---

## Task 9: Post-merge — out-of-repo memory update

**Files:**
- Modify: `~/.claude/projects/-Users-jin-holee-dev-GitHub-Jin-HoMLee-splice-neoepitope-pipeline-pm/memory/feedback_milestones.md`

- [ ] **Step 1: Update the sequencing prose section to reference the script**

In `feedback_milestones.md`, find the section starting `### Setting milestone due dates` (around line 116). After the "Formula" paragraph, add:

```markdown
**Sequencing implementation:** The sequencing logic above (same-S-stage stacking + paired-S7 gating) is implemented in [`scripts/pm/recheck_milestone.py`](scripts/pm/recheck_milestone.py) via the `compute_layered_due_date` function (added 2026-05-22 per [Issue #465](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/465)). When the recheck hook fires, proposed dates already account for sequencing — no need to mentally apply the rule on top of the script's output.
```

- [ ] **Step 2: No commit needed** (memory is out of repo)

Memory file changes propagate to the next session automatically.

---

## Self-Review Notes

**Spec coverage verified:**
- ✓ `parse_milestone_title` (Task 1)
- ✓ `find_prior_same_stage` (Task 2)
- ✓ `find_open_same_iteration_S5` (Task 3)
- ✓ `compute_layered_due_date` 7 branches (Task 4 + tests)
- ✓ 7 sequence-bound milestones → `[No change]` (Task 6 live smoke)
- ✓ Capacity-bound regression check (Task 6 live smoke)
- ✓ Lab notebook entry (Task 7)
- ✓ `feedback_milestones.md` reference (Task 9, post-merge)

**Type/name consistency check:**
- `compute_layered_due_date` signature matches between spec, Task 4 implementation, and Task 5 integration call
- `find_prior_same_stage` returns dict matching what `compute_layered_due_date` reads (`["state"]`, `["number"]`, `.get("due_on")`)
- All test fixture dicts include the keys read by helpers (`number`, `title`, `state`, `due_on`)

**Placeholder scan:** No TBDs, no "implement later", every code step has complete code.
