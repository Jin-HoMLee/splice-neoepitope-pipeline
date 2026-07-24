# Pullability Predicate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the scattered "can this Issue be worked now?" checks with one pullability predicate over natively-owned sources (Issue #1294, PR1).

**Architecture:** A pure `assess(item)` function in a new `scripts/pm/pullability.py` reads each cause from exactly one authoritative source - GitHub `blockedBy`, the `needs-design` label, a new `trigger-gated` label, and the `Start date` field - and returns one reason string from an enumerated taxonomy, or `None`. `board_open_items.py` computes its existing `not_pullable` field from the predicate instead of the prose scan; `not_pullable.py` is demoted to a low-recall proposer; a one-time backfill labels the currently-gated Issues for cutover parity.

**Tech Stack:** Python 3.13 (stdlib only, `re`/`typing`), pytest, bash + `jq`, `gh` CLI, GitHub ProjectV2 GraphQL.

## Global Constraints

- Pure functions only in `pullability.py`: no network, no stored state, re-derivable each run (mirrors `not_pullable.py`).
- Fail open everywhere: an unknown shape, missing field, or parse miss returns `None` (pullable), never a spurious gate - this feeds a board read and must never make a healthy queue look broken.
- Keep the output JSON key name `not_pullable` so `check_ready_queue.sh` needs no logic change.
- No em-dashes or en-dashes in any authored text; plain hyphen only.
- Pytest interpreter: `workflow/tests/.venv/bin/python -m pytest` (never a root `.venv`).
- New tests for the predicate and the proposer live in `workflow/tests/` (beside the existing `test_not_pullable.py` / `test_check_ready_queue.py`).
- `trigger-gated` is created with `gh label create` (a plain label), never via `updateProjectV2Field` single-select options (that path is destructive).

---

## File Structure

- Create: `scripts/pm/pullability.py` - the predicate + reason taxonomy (pure).
- Create: `workflow/tests/test_pullability.py` - predicate unit tests.
- Modify: `scripts/board_open_items.py` - add `blockedBy` to the query, capture `Start date`, compute `not_pullable` from the predicate.
- Modify: `scripts/pm/not_pullable.py` - demote to a low-recall proposer (docstring + a `propose_label` entry point).
- Modify: `workflow/tests/test_not_pullable.py` - proposer semantics + the four known-miss fixtures as accepted misses.
- Modify: `scripts/check_ready_queue.sh` - upgrade the cap comment to cite the WIP-limit convention.
- Operational (not committed code): `gh label create trigger-gated`; backfill labels/`Start date` on the standing gated Issue set.

---

### Task 1: The pullability predicate (pure function + taxonomy)

**Files:**
- Create: `scripts/pm/pullability.py`
- Test: `workflow/tests/test_pullability.py`

**Interfaces:**
- Produces: `assess(item: dict, *, today: Optional[str] = None) -> Optional[str]`. `item` keys read: `labels: list[str]`, `blocked_by: list[dict]` (each `{"number": int, "state": str}`), `start_date: Optional[str]` (ISO `YYYY-MM-DD`). Also module constants `NEEDS_DESIGN_LABEL = "needs-design"`, `TRIGGER_GATED_LABEL = "trigger-gated"`.

- [ ] **Step 1: Write the failing tests**

```python
# workflow/tests/test_pullability.py
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts", "pm"))

import pullability  # noqa: E402


def _item(**over):
    base = {"labels": [], "blocked_by": [], "start_date": None}
    base.update(over)
    return base


def test_pullable_when_nothing_gates():
    assert pullability.assess(_item(), today="2026-07-24") is None


def test_open_blocker_gates():
    it = _item(blocked_by=[{"number": 42, "state": "OPEN"}])
    assert pullability.assess(it, today="2026-07-24") == "blocked-by-issue: #42"


def test_closed_blocker_does_not_gate():
    it = _item(blocked_by=[{"number": 42, "state": "CLOSED"}])
    assert pullability.assess(it, today="2026-07-24") is None


def test_needs_design_label_gates():
    assert pullability.assess(_item(labels=["needs-design"]), today="2026-07-24") == "needs-design"


def test_trigger_gated_label_gates():
    assert pullability.assess(_item(labels=["trigger-gated"]), today="2026-07-24") == "trigger-gated"


def test_future_start_date_gates():
    it = _item(start_date="2026-12-01")
    assert pullability.assess(it, today="2026-07-24") == "date-gated: 2026-12-01"


def test_past_start_date_does_not_gate():
    it = _item(start_date="2026-01-01")
    assert pullability.assess(it, today="2026-07-24") is None


def test_ordering_blocker_wins_over_date():
    it = _item(blocked_by=[{"number": 9, "state": "OPEN"}], start_date="2026-12-01")
    assert pullability.assess(it, today="2026-07-24") == "blocked-by-issue: #9"


def test_fails_open_on_empty_item():
    assert pullability.assess({}, today="2026-07-24") is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_pullability.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'pullability'`.

- [ ] **Step 3: Write the predicate**

```python
# scripts/pm/pullability.py
"""One pullability predicate over natively-owned sources (Issue #1294).

Replaces four scattered "can this be worked now?" checks with a single
function. Each cause is read from exactly one authoritative source and never
copied, so no state can drift. Adding a fifth cause is one branch here.

The prose scanner (not_pullable.py) is NOT consulted here. It is a low-recall
proposer that suggests a label for a human to apply; the label - not the prose -
is what this predicate reads. That is the whole point of Issue #1294: the
decision has one home, and the sources stay where they are natively mastered.
"""
from typing import Optional

# The enumerated reason taxonomy: one authoritative source per cause. Adding a
# fifth cause is one constant plus one branch in assess(), in this one file.
NEEDS_DESIGN_LABEL = "needs-design"
TRIGGER_GATED_LABEL = "trigger-gated"


def _open_blockers(item):
    """Issue numbers of the still-open native blockedBy edges."""
    out = []
    for edge in item.get("blocked_by") or []:
        if (edge.get("state") or "").upper() != "CLOSED":
            num = edge.get("number")
            if num is not None:
                out.append(num)
    return out


def assess(item, today=None):
    # type: (dict, Optional[str]) -> Optional[str]
    """Return a short reason the Issue is not pullable now, or None if it is.

    `item` carries already-fetched structured fields:
      - "labels":     list[str] label names
      - "blocked_by": list[{"number", "state"}] native blockedBy edges
      - "start_date": Optional[str] ISO YYYY-MM-DD; a future value gates the Issue

    Sources are checked hardest-first: a hard dependency subsumes a softer gate,
    and a date gate ("not yet") is the weakest. Fails open on any unknown shape.
    """
    open_blockers = _open_blockers(item)
    if open_blockers:
        return "blocked-by-issue: #{}".format(open_blockers[0])

    labels = item.get("labels") or []
    if NEEDS_DESIGN_LABEL in labels:
        return "needs-design"
    if TRIGGER_GATED_LABEL in labels:
        return "trigger-gated"

    start = item.get("start_date")
    if start and today and start > today:
        return "date-gated: {}".format(start)

    return None
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_pullability.py -v`
Expected: PASS (9 passed).

- [ ] **Step 5: Commit**

```bash
git add scripts/pm/pullability.py workflow/tests/test_pullability.py
git commit -m "feat(pm): pullability predicate over natively-owned sources (#1294)"
```

---

### Task 2: Fetch `blockedBy` + `Start date` and wire the predicate into `board_open_items.py`

**Files:**
- Modify: `scripts/board_open_items.py` (QUERY Issue fragment ~line 96-101; `normalize()` fieldValues loop ~line 207-220; the return dict `not_pullable` ~line 284)
- Test: `scripts/tests/test_board_open_items_pullable.py` (new)

**Interfaces:**
- Consumes: `pullability.assess` from Task 1.
- Produces: each normalized item gains `"start_date"` and `"blocked_by"` keys; `"not_pullable"` is now `pullability.assess({...}, today=...)`.

- [ ] **Step 1: Write the failing cross-module test**

```python
# scripts/tests/test_board_open_items_pullable.py
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import board_open_items as boi  # noqa: E402


def _raw(number, labels=None, blocked=None, start=None):
    """A minimal GraphQL item node as normalize() expects it."""
    field_values = []
    field_values.append({
        "name": "Ready",
        "field": {"name": "Status"},
    })
    if start is not None:
        field_values.append({"date": start, "field": {"name": "Start date"}})
    return {
        "content": {
            "__typename": "Issue",
            "number": number,
            "title": f"#{number}",
            "state": "OPEN",
            "url": f"https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/{number}",
            "body": "",
            "labels": {"nodes": [{"name": n} for n in (labels or [])]},
            "blockedBy": {"nodes": blocked or []},
            "subIssuesSummary": {"total": 0},
        },
        "fieldValues": {"nodes": field_values},
    }


def test_needs_design_label_sets_not_pullable():
    it = boi.normalize(_raw(100, labels=["role:pm", "needs-design"]))
    assert it["not_pullable"] == "needs-design"


def test_clean_issue_is_pullable():
    it = boi.normalize(_raw(101, labels=["role:pm"]))
    assert it["not_pullable"] is None


def test_open_blocker_sets_not_pullable():
    it = boi.normalize(_raw(102, labels=["role:pm"], blocked=[{"number": 5, "state": "OPEN"}]))
    assert it["not_pullable"] == "blocked-by-issue: #5"
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest scripts/tests/test_board_open_items_pullable.py -v`
Expected: FAIL - `not_pullable` is computed from `scan_not_pullable(body)` (empty body -> None), so `test_needs_design_label_sets_not_pullable` and `test_open_blocker_sets_not_pullable` fail.

- [ ] **Step 3: Add `blockedBy` to the QUERY Issue fragment**

In `scripts/board_open_items.py`, the `... on Issue { ... }` fragment, add `blockedBy` after the `labels(first: 20) { nodes { name } }` line:

```graphql
              labels(first: 20) { nodes { name } }
              blockedBy(first: 50) { nodes { number state } }
```

(`first: 50` is GitHub's documented per-direction ceiling for blockers, matching `scan_unblocked.py`.)

- [ ] **Step 4: Capture `Start date` and `blocked_by` in `normalize()`**

Add a `start_date` accumulator beside `target_date`, capture it in the fieldValues loop, and extract `blocked_by` from content. Replace:

```python
    status = size = priority = target_date = None
    for fv in item["fieldValues"]["nodes"]:
        if not fv:
            continue
        field = fv.get("field", {}) or {}
        fname = field.get("name")
        if fname == "Status":
            status = fv.get("name")
        elif fname == "Size":
            size = fv.get("name")
        elif fname == "Priority":
            priority = fv.get("name")
        elif fname == "Target date":
            target_date = fv.get("date")
```

with:

```python
    status = size = priority = target_date = start_date = None
    for fv in item["fieldValues"]["nodes"]:
        if not fv:
            continue
        field = fv.get("field", {}) or {}
        fname = field.get("name")
        if fname == "Status":
            status = fv.get("name")
        elif fname == "Size":
            size = fv.get("name")
        elif fname == "Priority":
            priority = fv.get("name")
        elif fname == "Target date":
            target_date = fv.get("date")
        elif fname == "Start date":
            start_date = fv.get("date")
```

- [ ] **Step 5: Compute `not_pullable` from the predicate**

Add the import beside the existing one (`from not_pullable import scan_not_pullable`):

```python
from pullability import assess as assess_pullable  # noqa: E402
```

Add the blocked_by extraction near the `labels = [...]` line:

```python
    blocked_by = [b for b in (content.get("blockedBy") or {}).get("nodes", []) if b]
```

Replace the return dict's `not_pullable` line:

```python
        "not_pullable": scan_not_pullable(content.get("body")),
```

with (add `start_date` and `blocked_by` to the dict too, beside `target_date`):

```python
        "target_date": target_date,
        "start_date": start_date,
        "blocked_by": blocked_by,
        ...
        "not_pullable": assess_pullable(
            {"labels": labels, "blocked_by": blocked_by, "start_date": start_date},
            today=datetime.now(timezone.utc).date().isoformat(),
        ),
```

(Keep the `scan_not_pullable` import only if Task 3 still needs it here; Task 3 removes this caller. `datetime`/`timezone` are already imported.)

- [ ] **Step 6: Run the test to verify it passes**

Run: `workflow/tests/.venv/bin/python -m pytest scripts/tests/test_board_open_items_pullable.py -v`
Expected: PASS (3 passed).

- [ ] **Step 7: Run the full board_open_items + ready-queue suites for regressions**

Run: `workflow/tests/.venv/bin/python -m pytest scripts/tests/ workflow/tests/test_check_ready_queue.py -v`
Expected: PASS (no regressions).

- [ ] **Step 8: Commit**

```bash
git add scripts/board_open_items.py scripts/tests/test_board_open_items_pullable.py
git commit -m "feat(pm): wire pullability predicate into board_open_items (#1294)"
```

---

### Task 3: Demote `not_pullable.py` to a low-recall proposer

**Files:**
- Modify: `scripts/pm/not_pullable.py` (docstring + a `propose_label` entry point; the scan functions stay)
- Modify: `workflow/tests/test_not_pullable.py` (proposer semantics + four known-miss fixtures)
- Modify: `scripts/board_open_items.py` (drop the now-unused `scan_not_pullable` import)

**Interfaces:**
- Produces: `propose_label(body: Optional[str]) -> Optional[str]` returning `"needs-design"`, `"trigger-gated"`, or `None` - a suggestion, never an authority.

- [ ] **Step 1: Write the failing proposer tests**

```python
# add to workflow/tests/test_not_pullable.py
def test_proposer_suggests_needs_design_for_decision_ac():
    body = "## Acceptance criteria\n- [ ] Decide scope: A or B\n"
    assert not_pullable.propose_label(body) == "needs-design"


def test_proposer_suggests_trigger_gated_for_trigger_marker():
    body = "## Notes\nBuild only if it slips a 2nd time.\n"
    assert not_pullable.propose_label(body) == "trigger-gated"


def test_proposer_returns_none_for_clean_body():
    assert not_pullable.propose_label("## Acceptance criteria\n- [ ] Ship the thing\n") is None


# Known-miss fixtures: accepted low-recall behavior. A miss is asserted as an
# accepted non-event, because the label is the authority now (Issue #1294).
def test_known_miss_noun_form_decision_is_accepted():
    # #1111 / #353 noun-form: the proposer does NOT catch these, by design.
    body = "## Acceptance criteria\n- [ ] A decision (A or B) is recorded here.\n"
    assert not_pullable.propose_label(body) is None  # accepted miss


def test_known_miss_no_body_gate_is_accepted():
    # #876: no body gate at all; the proposer has nothing to see.
    assert not_pullable.propose_label("## Context\nArchive when needed.\n") is None
```

- [ ] **Step 2: Run to verify failure**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_not_pullable.py -k proposer -v`
Expected: FAIL with `AttributeError: module 'not_pullable' has no attribute 'propose_label'`.

- [ ] **Step 3: Add `propose_label` and rewrite the module docstring**

Add to `scripts/pm/not_pullable.py`, mapping the existing scan to a label suggestion:

```python
def propose_label(body):
    # type: (Optional[str]) -> Optional[str]
    """Suggest a gate label for a human to apply at Backlog -> Ready.

    LOW RECALL IS INTENDED, NOT A DEFECT. The label is the authority (read by
    scripts/pm/pullability.py); this only nudges a human to set it. A missed
    phrasing is a non-event, so never grow the pattern set to chase completeness
    - that is the denylist widen Issue #1294 exists to avoid.
    """
    reason = scan_not_pullable(body)
    if reason is None:
        return None
    if reason.startswith("trigger-gated"):
        return "trigger-gated"
    if reason.startswith("needs-design"):
        return "needs-design"
    return None
```

Replace the module docstring's framing so its first paragraph states it is a low-recall proposer, not the floor guard's authority (the authority is now `pullability.py`).

- [ ] **Step 4: Drop the unused authority wiring**

In `scripts/board_open_items.py`, remove the now-unused import `from not_pullable import scan_not_pullable` (Task 2 replaced its only caller with the predicate).

- [ ] **Step 5: Run to verify pass**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_not_pullable.py scripts/tests/test_board_open_items_pullable.py -v`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add scripts/pm/not_pullable.py workflow/tests/test_not_pullable.py scripts/board_open_items.py
git commit -m "refactor(pm): demote not_pullable prose scan to a low-recall proposer (#1294)"
```

---

### Task 4: Cap-comment WIP citation + create the `trigger-gated` label

**Files:**
- Modify: `scripts/check_ready_queue.sh` (the cap comment only)
- Test: `workflow/tests/test_check_ready_queue.py` (assert the WIP-convention wording is present, if the suite asserts comments; otherwise a shellcheck/no-op run)

**Interfaces:** none (comment + operational label).

- [ ] **Step 1: Upgrade the cap comment**

In `scripts/check_ready_queue.sh`, near the cap logic (around the `[CAP]` / pullable-count block ~line 159-179), replace the cap rationale comment with wording that cites the convention:

```bash
# Cap counts ALL Ready cards (gated included); the floor counts only pullable
# ones. This asymmetry is the documented Kanban WIP-limit convention - blocked
# items still consume a WIP slot, so they count against the cap, but they cannot
# satisfy the Ready floor because they are not workable. Not a local choice.
```

- [ ] **Step 2: Run the ready-queue suite to confirm no behavior change**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_check_ready_queue.py -v`
Expected: PASS (comment-only edit; behavior unchanged).

- [ ] **Step 3: Create the `trigger-gated` label (operational)**

Run:
```bash
gh label create trigger-gated \
  --repo Jin-HoMLee/splice-neoepitope-pipeline \
  --color B60205 \
  --description "Pullability gate: Issue waits on an external event/trigger (read by pullability.py)"
```
Expected: label created (or "already exists" - idempotent, safe to re-run).

- [ ] **Step 4: Commit the comment change**

```bash
git add scripts/check_ready_queue.sh
git commit -m "docs(pm): cap comment cites the WIP-limit convention (#1294)"
```

---

### Task 5: One-time backfill for cutover parity (operational, freshness-checked)

**Files:** none committed (board-state mutations). Verification recorded in the PR Test plan.

**Interfaces:** none.

- [ ] **Step 1: Freshness-check and label each standing gated Issue**

For each Issue below, `gh issue view <N>` first to confirm the gate is still live, then apply the marker. Do NOT bulk-apply.

```bash
# trigger-gated (external event / revival):
gh issue edit 841 --repo Jin-HoMLee/splice-neoepitope-pipeline --add-label trigger-gated
gh issue edit 876 --repo Jin-HoMLee/splice-neoepitope-pipeline --add-label trigger-gated
# needs-design (open decision fork):
gh issue edit 929 --repo Jin-HoMLee/splice-neoepitope-pipeline --add-label needs-design
# date-gated: set Start date (board field) via scripts/pm/set_status.sh sibling or a
# direct updateProjectV2ItemFieldValue on the Start date field:
#   #1247 -> 2026-08-03 ; #817 -> 2026-12-01
```

For the two decision-shaped Backlog items, apply `needs-design` only if the freshness read confirms the fork is still open:
```bash
gh issue edit 1111 --repo Jin-HoMLee/splice-neoepitope-pipeline --add-label needs-design  # if still open
gh issue edit 353  --repo Jin-HoMLee/splice-neoepitope-pipeline --add-label needs-design  # if still open
```

Re-read [Issue #585](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/585) to classify it (condition vs trigger vs date) before labeling.

- [ ] **Step 2: Set `Start date` on the date-gated Issues**

Get each card's project item id and set the `Start date` field. The Start date field id is discoverable via:
```bash
gh api graphql -f query='{ user(login:"Jin-HoMLee"){ projectV2(number:9){ f:field(name:"Start date"){ ... on ProjectV2FieldCommon { id } } } } }'
```
Then for each `(item_id, date)`:
```bash
gh api graphql -f query='mutation($p:ID!,$i:ID!,$f:ID!,$d:Date!){ updateProjectV2ItemFieldValue(input:{projectId:$p,itemId:$i,fieldId:$f,value:{date:$d}}){ projectV2Item{ id } } }' \
  -f p=PVT_kwHOB17eGc4BSomP -f i="<ITEM_ID>" -f f="<START_DATE_FIELD_ID>" -f d="2026-08-03"
```

- [ ] **Step 3: Live smoke - matched-pair control**

Run the queue check and confirm the backfilled Issues are now excluded from the pullable count and named on the `[NOT-PULLABLE ...]` line:
```bash
bash scripts/check_ready_queue.sh 2>&1 | grep -E "NOT-PULLABLE|Ready"
```
Expected: the labeled Issues appear in the `[NOT-PULLABLE ...]` line with their taxonomy reason; they do not count toward any role's pullable floor. (Matched-pair: before the backfill they leaked into the pullable count; after, they do not.)

- [ ] **Step 4: Record the smoke result in the PR Test plan**

No commit (board mutations only). Paste the `check_ready_queue.sh` before/after excerpt into the PR body Test plan.

---

## Self-Review

**Spec coverage:**
- Predicate + one taxonomy, sole decider -> Task 1 (function) + Task 2 (wiring; `check_ready_queue.sh` reads it via the kept `not_pullable` key).
- Each cause reads its native source, nothing copied -> Task 1 (`assess` reads labels/blockedBy/start_date directly).
- Body scan demoted to an intended-low-recall proposer -> Task 3.
- Four known misses as fixtures -> Task 3 (proposer tests) + Task 1 (the label-authority side).
- `trigger-gated` structured home -> Task 4 (label) + Task 5 (backfill).
- Cap comment cites WIP convention -> Task 4.
- `#715` constraint satisfied (User account, no Issue Types) -> spec-recorded; the plan uses labels + a date field, never a Type field.
- Deferred (`#876` structured state, MM DoR edit) -> out of this plan by design; `#876` still gets a `trigger-gated` label in Task 5 so the queue excludes it.

**Placeholder scan:** no TBD/TODO; every code step carries real code. The backfill Issue numbers are concrete; the per-Issue "if still open" is a required freshness gate, not a placeholder.

**Type consistency:** `assess(item, today=)` returns `Optional[str]`; `item` keys `labels`/`blocked_by`/`start_date` are used identically in Task 1, Task 2's dict construction, and the tests. `propose_label` returns the label strings `"needs-design"`/`"trigger-gated"` that Task 1's constants and Task 5's `--add-label` use verbatim.
</content>
