# Sequencing-aware milestone recheck — Issue #465

**Status:** Design approved 2026-05-22
**Issue:** [#465](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/465)
**Author:** PM (with Claude Opus 4.7)

## Problem

`scripts/pm/recheck_milestone.py` computes a milestone's "proposed" `due_on` using pure capacity:

```
proposed = today + (remaining_capacity / RATE * 7)
```

After today's rate change (1.5 → 5.0 d/wk) and resulting 14-milestone re-baseline, 7 of those milestones intentionally got dates *later* than pure capacity says because they stack behind upstream same-S-stage iterations. The hook fires on every subsequent activity (issue close, milestone edit, size change) and re-flags those 7 as `[UPDATE NEEDED]` — false alarm noise that tempts the operator to "fix" them by patching down to capacity-only dates, which would re-introduce sequencing violations.

The sequencing rule already lives in `feedback_milestones.md`: *"Same-stage iterations stack: `i3 - S5` starts when `i2 - S5` closes."* This spec moves that rule from memory-led to mechanism-led, per [[mechanism-over-memory]] and [[deterministic-before-semantic]].

## Goal

Extend `recheck_milestone.py` so its proposed dates respect **same-S-stage sequencing** and **paired-S7 gating**, eliminating the false-flag tax for the 7 sequence-bound milestones while preserving current behavior for the rest.

## Out of scope

- **Role-meta-axis stacking** (pm-i*, dev-i*): they run partly in parallel in practice; not worth the complexity.
- **Recursive proposed-close computation**: trust GitHub's stored `due_on` as the source of truth. The hook's cascade-on-activity property converges correctly over time.
- **Auto-PATCH**: keep the existing report-only design. Operator runs PATCH manually after reviewing the output.
- **Arc-mismatch enforcement** (e.g. i4-S7 "TCR-pMHC Scorer Landscape" not arc-matching i4-S5 "Google Batch"): a separate data-hygiene concern, not solvable here.

## Design

### New function: `compute_layered_due_date`

Signature:

```python
def compute_layered_due_date(
    iteration: int | None,
    stage: int | None,
    capacity_days: float,
    all_milestones: list[dict],
) -> tuple[date, str]:
    """Return (proposed_due, reasoning_note).

    Replaces the inline `proposed_due = today + capacity/RATE*7` calculation
    on the current line 120 of recheck_milestone.py.
    """
```

`all_milestones` is the result of `gh api repos/.../milestones?state=all&per_page=100` — fetched once at the top of `compute_recheck`, passed down so sub-helpers don't re-query.

### Logic flow

```
parse_milestone_title(meta["title"]) -> (iteration, stage) or None

if no parse:
    # Non-S-stage milestone (pm-i*, dev-i*, M1, etc.) — pure capacity
    base = today
    note = ""

elif stage == 7:
    # Paired-S7 sub-rule (takes precedence)
    paired = find_open_same_iteration_S5(iteration, all_milestones)
    if paired and paired.due_on:
        base = max(date.fromisoformat(paired.due_on), today)
        note = f"(paired-S7: unblocks at M#{paired.number} close {paired.due_on})"
    else:
        # Standalone S7 (e.g. Lit Review i3-S7) — pure capacity
        base = today
        note = "(standalone S7 — no paired open S5)"

else:
    # Same-S-stage stacking
    prior = find_prior_same_stage(iteration, stage, all_milestones)
    if prior is None:
        base = today
        note = "(no prior same-S milestone)"
    elif prior["state"] == "closed":
        base = today
        note = f"(prior M#{prior['number']} closed)"
    elif prior["due_on"] is None:
        base = today
        note = f"(prior M#{prior['number']} undated — sequencing skipped)"
    else:
        base = max(date.fromisoformat(prior["due_on"][:10]), today)
        note = f"(stack after M#{prior['number']} close {prior['due_on'][:10]})"

proposed = base + timedelta(days=round(capacity_days / AVAILABILITY_RATE * 7))
return (proposed, note)
```

### Helper functions

```python
def parse_milestone_title(title: str) -> tuple[int, int] | None:
    """Parse 'i<N> - S<M> - ...' titles. Returns (N, M) or None."""
    m = re.match(r"^i(\d+)\s*-\s*S(\d+)\s*-\s*", title)
    return (int(m.group(1)), int(m.group(2))) if m else None


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
    return candidates[-1][1]  # highest N


def find_open_same_iteration_S5(
    iteration: int,
    all_milestones: list[dict],
) -> dict | None:
    """Find an OPEN i<N> - S5 - ... milestone. Used for paired-S7 gating."""
    for ms in all_milestones:
        parsed = parse_milestone_title(ms["title"])
        if parsed is None:
            continue
        n, s = parsed
        if n == iteration and s == 5 and ms["state"] == "open":
            return ms
    return None
```

### Output format

Append the reasoning note to the existing `Proposed due_on:` line:

```
Proposed due_on: 2026-06-09 (delta -10 days) (stack after M#5 close 2026-05-31)
```

`Status: [No change]` continues to trigger when `abs(delta) <= THRESHOLD_DAYS` (7 days, unchanged).

### Integration point

`compute_recheck` (lines 85-129 of current script):

- Fetch `all_milestones` once near the top.
- Replace inline `proposed_due = date.today() + timedelta(days=calendar_days)` with `proposed_due, note = compute_layered_due_date(iteration, stage, remaining, all_milestones)`.
- Append `note` to the print of `Proposed due_on:` line.

Net diff: ~80 LOC added (3 helpers + 1 layered-compute function + integration), ~3 LOC modified.

## Tests

New file: `tools/ci/test_recheck_milestone.py` (pytest, runs via `workflow/tests/.venv/bin/python -m pytest tools/ci/test_recheck_milestone.py`). Mirrors the sibling `tools/ci/test_recheck_parent_status.py` convention for PM-script tests (distinct from `workflow/tests/` which is pipeline-only).

### Unit tests

**`parse_milestone_title`:**
- `"i3 - S5 - Modeling - ..."` → `(3, 5)`
- `"pm-i4 - PM Tooling, ..."` → `None`
- `"dev-i3 - ..."` → `None`
- `"i2 - S7 - Publication - ..."` → `(2, 7)`
- `"M1 - Problem Definition v1"` → `None`
- `"i10 - S3 - ..."` → `(10, 3)` (double-digit safety)

**`find_prior_same_stage`:**
- Fixture: 5 milestones spanning i1/i2/i3/i4/i5-S3
- For i5-S3 → returns i4-S3
- For i3-S3 → returns i2-S3
- For i1-S3 → returns `None`
- For i2-S5 with no S5 priors → returns `None`

**`find_open_same_iteration_S5`:**
- For iteration=2 with open i2-S5 in fixture → returns i2-S5
- For iteration=2 with i2-S5 closed → returns `None`
- For iteration=3 with no i3-S5 → returns `None`

**`compute_layered_due_date`:**
- Closed prior → pure capacity from today
- Undated prior → pure capacity + skip-note
- Normal stack with prior in future → stacks after prior.due_on
- Normal stack with prior overdue → uses today (max guard)
- S7 with paired open S5 → uses S5.due_on
- S7 with no paired S5 → standalone, pure capacity
- Non-S-stage title (pm-i4) → pure capacity, empty note

### Integration smoke test

Live API call against current state (not mocked; runs only when `--live` flag passed, skipped in CI):

```python
@pytest.mark.live
def test_seven_known_sequence_bound_milestones_now_no_change():
    """The 7 milestones from 2026-05-22's session should all report [No change]."""
    expected_no_change = [10, 11, 13, 15, 18, 24, 30]
    for ms in expected_no_change:
        result = subprocess.run(
            ["python3", "scripts/pm/recheck_milestone.py", "--milestone", str(ms)],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, (
            f"M#{ms} still flagging UPDATE NEEDED:\n{result.stdout}"
        )
```

This is the acceptance test from the Issue body, mechanized.

## Verification before completion

1. All unit tests pass.
2. Live smoke test confirms all 7 named milestones report `[No change]`.
3. The 9 capacity-bound milestones from today's session still report `[No change]` (no regression).
4. Manually verify role-meta milestones (pm-i2, pm-i4) unaffected.
5. Verify the hook still surfaces the recheck on test PATCH (don't re-test the hook itself; it's an unchanged integration point).

## Migration / rollout

- Single PR contains spec + implementation + tests + lab notebook entry.
- No `due_on` migrations needed — existing PATCH'd dates from today are correct under the new logic.
- No memory updates needed beyond pointing `feedback_milestones.md` at the script for the sequencing math (currently prose-only).

## Open questions

None at this point. All edge cases (closed prior, undated prior, no prior, S7-paired/standalone, role-meta) handled explicitly in the design.

## Acceptance criteria (mirrors Issue #465)

- [ ] `parse_milestone_title` handles all 6 title shapes correctly
- [ ] `find_prior_same_stage` picks highest-iteration prior, includes closed
- [ ] `find_open_same_iteration_S5` filters by open + iteration
- [ ] `compute_layered_due_date` covers all 7 branches (no-parse, closed-prior, undated-prior, no-prior, normal-stack, S7-paired, S7-standalone)
- [ ] All 7 sequence-bound milestones (M#10/#11/#13/#15/#18/#24/#30) report `[No change]` after the fix
- [ ] Capacity-bound milestones still flag correctly when capacity changes (regression check)
- [ ] Lab notebook entry capturing the verification run
- [ ] `feedback_milestones.md` updated to reference the script for sequencing math
