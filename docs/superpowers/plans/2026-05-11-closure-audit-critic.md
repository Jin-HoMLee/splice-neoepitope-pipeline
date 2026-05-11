# Closure-Audit Critic Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Auto-run the closure-ritual checks (AC ticked, lab notebook entry exists, Priority rationale present) on PR-merge and issue-close events, posting an idempotent comment listing any gaps.

**Architecture:** GitHub Actions workflow fires on `pull_request.closed` (gated by `merged == true`) and `issues.closed`. Workflow runs `workflow/scripts/closure_audit.py` which fetches state via `gh` CLI, calls pure-function checks in `workflow/scripts/closure_audit_lib.py`, and posts or edits a single marker-tagged comment on the closing PR/issue.

**Tech Stack:** Python 3.11 (CI), `gh` CLI (preinstalled on `ubuntu-latest`), pytest 8.x (existing test infra), GitHub Actions.

**Spec:** [docs/superpowers/specs/2026-05-11-post-merge-critic-design.md](../specs/2026-05-11-post-merge-critic-design.md)

---

## File Structure

**Created:**

- `workflow/scripts/closure_audit_lib.py` — pure functions: check_ac_boxes, check_priority_rationale, check_lab_notebook, is_lab_notebook_exempt, resolve_roles, format_comment
- `workflow/scripts/closure_audit.py` — orchestration: gh CLI I/O, comment posting/editing
- `workflow/tests/test_closure_audit.py` — pytest unit tests (≥1 per branch in lib)
- `workflow/tests/fixtures/closure_audit/` — fixture .md files (issue/PR bodies, lab notebooks)
- `.github/workflows/closure-audit.yml` — event triggers + dispatch

**Modified:**

- None. Test discovery is already wired (`pytest workflow/tests/ -v` in `tests.yml:59`).

---

## Task 1: Scaffold lib + test file + fixtures dir

**Files:**
- Create: `workflow/scripts/closure_audit_lib.py`
- Create: `workflow/tests/test_closure_audit.py`
- Create: `workflow/tests/fixtures/closure_audit/.gitkeep`

- [ ] **Step 1: Create empty lib file**

Write `workflow/scripts/closure_audit_lib.py` with this content:

```python
"""Pure-function library for the closure-audit critic.

Network-free. Consumes already-fetched issue/PR bodies, comments, labels,
and lab notebook text. Tested via workflow/tests/test_closure_audit.py.

See docs/superpowers/specs/2026-05-11-post-merge-critic-design.md for the
3-check specification.
"""

from __future__ import annotations
```

- [ ] **Step 2: Create empty test file**

Write `workflow/tests/test_closure_audit.py`:

```python
"""Tests for closure_audit_lib."""

import closure_audit_lib as lib
```

- [ ] **Step 3: Create fixtures directory**

```bash
mkdir -p workflow/tests/fixtures/closure_audit
touch workflow/tests/fixtures/closure_audit/.gitkeep
```

- [ ] **Step 4: Verify pytest discovers the new file (zero tests pass)**

Run: `pytest workflow/tests/test_closure_audit.py -v`
Expected: `no tests ran in <duration>` exit 5 (pytest exits 5 when no tests collected — that's OK at this stage).

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/closure_audit_lib.py workflow/tests/test_closure_audit.py workflow/tests/fixtures/closure_audit/.gitkeep
git commit -m "chore(closure-audit): scaffold lib + test file + fixtures dir (#325)"
```

---

## Task 2: Check 1 — AC checkboxes (TDD)

**Files:**
- Modify: `workflow/scripts/closure_audit_lib.py`
- Modify: `workflow/tests/test_closure_audit.py`
- Create: `workflow/tests/fixtures/closure_audit/issue_all_ticked.md`
- Create: `workflow/tests/fixtures/closure_audit/issue_some_unticked.md`
- Create: `workflow/tests/fixtures/closure_audit/issue_no_checkboxes.md`

- [ ] **Step 1: Create fixture — all ticked**

Write `workflow/tests/fixtures/closure_audit/issue_all_ticked.md`:

```markdown
**Created by:** Developer

## Acceptance criteria

- [x] Tests pass
- [x] Docs updated
- [x] Reviewer signed off

**Priority rationale:** P2 — example.
```

- [ ] **Step 2: Create fixture — some unticked**

Write `workflow/tests/fixtures/closure_audit/issue_some_unticked.md`:

```markdown
**Created by:** Developer

## Acceptance criteria

- [x] Tests pass
- [ ] Docs updated
- [ ] Reviewer signed off

**Priority rationale:** P2 — example.
```

- [ ] **Step 3: Create fixture — no checkboxes (legacy bullets)**

Write `workflow/tests/fixtures/closure_audit/issue_no_checkboxes.md`:

```markdown
**Created by:** Developer

## Acceptance criteria

- Tests pass
- Docs updated

**Priority rationale:** P2 — example.
```

- [ ] **Step 4: Write failing tests**

Append to `workflow/tests/test_closure_audit.py`:

```python
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures" / "closure_audit"


def _read(name: str) -> str:
    return (FIXTURES / name).read_text()


def test_check_ac_boxes_all_ticked_no_gap():
    body = _read("issue_all_ticked.md")
    has_gap, _ = lib.check_ac_boxes(body, comments=[])
    assert has_gap is False


def test_check_ac_boxes_some_unticked_no_deferral_gap():
    body = _read("issue_some_unticked.md")
    has_gap, detail = lib.check_ac_boxes(body, comments=[])
    assert has_gap is True
    assert "2/3 unticked" in detail or "2 unticked" in detail


def test_check_ac_boxes_some_unticked_with_deferral_no_gap():
    body = _read("issue_some_unticked.md")
    comments = ["❎ Criterion 'Docs updated' deferred to follow-up #99 — out of scope."]
    has_gap, _ = lib.check_ac_boxes(body, comments)
    assert has_gap is False


def test_check_ac_boxes_no_checkboxes_no_gap():
    """Legacy bodies with plain bullets have no checkboxes to fail on."""
    body = _read("issue_no_checkboxes.md")
    has_gap, _ = lib.check_ac_boxes(body, comments=[])
    assert has_gap is False
```

- [ ] **Step 5: Run tests to verify they fail**

Run: `pytest workflow/tests/test_closure_audit.py -v`
Expected: 4 failures, all `AttributeError: module 'closure_audit_lib' has no attribute 'check_ac_boxes'`.

- [ ] **Step 6: Implement**

Append to `workflow/scripts/closure_audit_lib.py`:

```python
import re

# Captures "- [ ]" or "- [x]" with flexible spacing. Case-insensitive for [x]/[X].
_UNTICKED_RE = re.compile(r"^\s*-\s*\[\s\]\s", re.MULTILINE)
_TICKED_RE = re.compile(r"^\s*-\s*\[[xX]\]\s", re.MULTILINE)


def check_ac_boxes(body: str, comments: list[str]) -> tuple[bool, str]:
    """Return (has_gap, detail) for AC checkbox state.

    v1 simplification: if any comment contains both '❎' and 'deferred' (case-
    insensitive), all unticked boxes are treated as deferred. See spec §
    "v1 simplification" for rationale.
    """
    unticked = len(_UNTICKED_RE.findall(body))
    ticked = len(_TICKED_RE.findall(body))
    total = unticked + ticked
    deferral_seen = any("❎" in c and "deferred" in c.lower() for c in comments)
    if unticked == 0 or deferral_seen:
        return False, ""
    return True, f"{unticked}/{total} unticked, no deferral comment found"
```

- [ ] **Step 7: Run tests to verify they pass**

Run: `pytest workflow/tests/test_closure_audit.py -v`
Expected: 4 passed.

- [ ] **Step 8: Commit**

```bash
git add workflow/scripts/closure_audit_lib.py workflow/tests/test_closure_audit.py workflow/tests/fixtures/closure_audit/
git commit -m "feat(closure-audit): Check 1 — AC checkbox audit (#325)"
```

---

## Task 3: Check 3 — Priority rationale (TDD)

**Files:**
- Modify: `workflow/scripts/closure_audit_lib.py`
- Modify: `workflow/tests/test_closure_audit.py`
- Create: `workflow/tests/fixtures/closure_audit/issue_no_rationale.md`
- Create: `workflow/tests/fixtures/closure_audit/issue_rationale_formatting_variant.md`

- [ ] **Step 1: Create fixture — no rationale**

Write `workflow/tests/fixtures/closure_audit/issue_no_rationale.md`:

```markdown
**Created by:** Developer

## Acceptance criteria

- [x] Tests pass
```

- [ ] **Step 2: Create fixture — formatting variant (no asterisks, header form)**

Write `workflow/tests/fixtures/closure_audit/issue_rationale_formatting_variant.md`:

```markdown
**Created by:** Developer

## Priority rationale

P1 — blocks the i3 milestone.

## Acceptance criteria

- [x] Tests pass
```

- [ ] **Step 3: Write failing tests**

Append to `workflow/tests/test_closure_audit.py`:

```python
def test_check_priority_rationale_present_no_gap():
    body = _read("issue_all_ticked.md")  # has "**Priority rationale:**" line
    has_gap, _ = lib.check_priority_rationale(body)
    assert has_gap is False


def test_check_priority_rationale_formatting_variant_no_gap():
    body = _read("issue_rationale_formatting_variant.md")  # header form
    has_gap, _ = lib.check_priority_rationale(body)
    assert has_gap is False


def test_check_priority_rationale_missing_gap():
    body = _read("issue_no_rationale.md")
    has_gap, detail = lib.check_priority_rationale(body)
    assert has_gap is True
    assert "priority rationale" in detail.lower()
```

- [ ] **Step 4: Run tests to verify they fail**

Run: `pytest workflow/tests/test_closure_audit.py -v -k priority_rationale`
Expected: 3 failures, `AttributeError: module 'closure_audit_lib' has no attribute 'check_priority_rationale'`.

- [ ] **Step 5: Implement**

Append to `workflow/scripts/closure_audit_lib.py`:

```python
def check_priority_rationale(body: str) -> tuple[bool, str]:
    """Loose case-insensitive substring match per PM audit method."""
    if "priority rationale" in body.lower():
        return False, ""
    return True, "no 'Priority rationale' line found in issue body"
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `pytest workflow/tests/test_closure_audit.py -v`
Expected: 7 passed.

- [ ] **Step 7: Commit**

```bash
git add workflow/scripts/closure_audit_lib.py workflow/tests/test_closure_audit.py workflow/tests/fixtures/closure_audit/
git commit -m "feat(closure-audit): Check 3 — Priority rationale audit (#325)"
```

---

## Task 4: Check 2 — Lab notebook entry (TDD)

**Files:**
- Modify: `workflow/scripts/closure_audit_lib.py`
- Modify: `workflow/tests/test_closure_audit.py`
- Create: `workflow/tests/fixtures/closure_audit/notebook_with_entry.md`
- Create: `workflow/tests/fixtures/closure_audit/notebook_date_only_no_ref.md`
- Create: `workflow/tests/fixtures/closure_audit/notebook_no_date.md`

- [ ] **Step 1: Create fixture — entry present (date header + sub-section + PR ref)**

Write `workflow/tests/fixtures/closure_audit/notebook_with_entry.md`:

```markdown
# Lab Notebook — Developer

---

## 2026-05-11

### 14:30 UTC — Editor: Developer

**Headline:** Shipped closure-audit critic ([PR #325](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/325)).

Some content here.

---

## 2026-05-10

### 17:04 UTC — Editor: Developer

Older entry.
```

- [ ] **Step 2: Create fixture — date present, no PR/issue reference**

Write `workflow/tests/fixtures/closure_audit/notebook_date_only_no_ref.md`:

```markdown
# Lab Notebook — Developer

---

## 2026-05-11

### 14:30 UTC — Editor: Developer

**Headline:** Some unrelated work today.

---
```

- [ ] **Step 3: Create fixture — no date header at all**

Write `workflow/tests/fixtures/closure_audit/notebook_no_date.md`:

```markdown
# Lab Notebook — Developer

---

## 2026-05-10

### 17:04 UTC — Editor: Developer

Older entry. Today's date 2026-05-11 not present.
```

- [ ] **Step 4: Write failing tests**

Append to `workflow/tests/test_closure_audit.py`:

```python
def test_check_lab_notebook_entry_present_no_gap():
    text = _read("notebook_with_entry.md")
    has_gap, _ = lib.check_lab_notebook(text, date="2026-05-11", number=325)
    assert has_gap is False


def test_check_lab_notebook_date_only_no_ref_gap():
    text = _read("notebook_date_only_no_ref.md")
    has_gap, detail = lib.check_lab_notebook(text, date="2026-05-11", number=325)
    assert has_gap is True
    assert "#325" in detail


def test_check_lab_notebook_no_date_gap():
    text = _read("notebook_no_date.md")
    has_gap, detail = lib.check_lab_notebook(text, date="2026-05-11", number=325)
    assert has_gap is True
    assert "2026-05-11" in detail


def test_check_lab_notebook_wrong_number_gap():
    """Date+sub-section exist but reference a different PR — should gap."""
    text = _read("notebook_with_entry.md")  # references #325
    has_gap, detail = lib.check_lab_notebook(text, date="2026-05-11", number=999)
    assert has_gap is True
    assert "#999" in detail
```

- [ ] **Step 5: Run tests to verify they fail**

Run: `pytest workflow/tests/test_closure_audit.py -v -k lab_notebook`
Expected: 4 failures, `AttributeError: module 'closure_audit_lib' has no attribute 'check_lab_notebook'`.

- [ ] **Step 6: Implement**

Append to `workflow/scripts/closure_audit_lib.py`:

```python
def check_lab_notebook(text: str, date: str, number: int) -> tuple[bool, str]:
    """Verify lab notebook has a date section that references the PR/issue.

    Args:
        text: full content of research/lab_notebook/<role>.md
        date: UTC date string 'YYYY-MM-DD' (merge or close date)
        number: PR number (for PR-merge) or issue number (for solo close)

    Returns (has_gap, detail).
    """
    date_header = f"## {date}"
    if date_header not in text:
        return True, f"no '## {date}' header in lab notebook"

    # Slice the date block: from this header to the next '## ' header or EOF.
    start = text.index(date_header)
    rest = text[start + len(date_header):]
    next_header_match = re.search(r"^## ", rest, re.MULTILINE)
    block = rest[:next_header_match.start()] if next_header_match else rest

    if not re.search(r"^### ", block, re.MULTILINE):
        return True, f"'## {date}' header exists but no '### HH:MM UTC — Editor: <Role>' sub-section under it"

    if f"#{number}" not in block:
        return True, f"'## {date}' block does not reference '#{number}'"

    return False, ""
```

- [ ] **Step 7: Run tests to verify they pass**

Run: `pytest workflow/tests/test_closure_audit.py -v`
Expected: 11 passed.

- [ ] **Step 8: Commit**

```bash
git add workflow/scripts/closure_audit_lib.py workflow/tests/test_closure_audit.py workflow/tests/fixtures/closure_audit/
git commit -m "feat(closure-audit): Check 2 — lab notebook entry audit (#325)"
```

---

## Task 5: Exemption-path filter (TDD)

**Files:**
- Modify: `workflow/scripts/closure_audit_lib.py`
- Modify: `workflow/tests/test_closure_audit.py`

- [ ] **Step 1: Write failing tests**

Append to `workflow/tests/test_closure_audit.py`:

```python
def test_is_lab_notebook_exempt_news_log_only():
    assert lib.is_lab_notebook_exempt(["research/news_log.md"]) is True


def test_is_lab_notebook_exempt_glossary_only():
    assert lib.is_lab_notebook_exempt(["research/glossary.md"]) is True


def test_is_lab_notebook_exempt_lab_notebook_files():
    """The lab notebook files themselves are exempt (avoids circular check)."""
    assert lib.is_lab_notebook_exempt(["research/lab_notebook/developer.md"]) is True
    assert lib.is_lab_notebook_exempt([
        "research/lab_notebook/developer.md",
        "research/lab_notebook/pm.md",
    ]) is True


def test_is_lab_notebook_exempt_combination_of_exempt_paths():
    assert lib.is_lab_notebook_exempt([
        "research/news_log.md",
        "research/glossary.md",
    ]) is True


def test_is_lab_notebook_exempt_mixed_not_exempt():
    """One non-exempt path disqualifies the whole PR."""
    assert lib.is_lab_notebook_exempt([
        "research/news_log.md",
        "workflow/scripts/run_mhcflurry.py",
    ]) is False


def test_is_lab_notebook_exempt_empty_list_not_exempt():
    """An empty changed-file list should NOT be treated as exempt."""
    assert lib.is_lab_notebook_exempt([]) is False
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest workflow/tests/test_closure_audit.py -v -k is_lab_notebook_exempt`
Expected: 6 failures, `AttributeError: module 'closure_audit_lib' has no attribute 'is_lab_notebook_exempt'`.

- [ ] **Step 3: Implement**

Append to `workflow/scripts/closure_audit_lib.py`:

```python
_EXEMPT_EXACT = {"research/news_log.md", "research/glossary.md"}
_EXEMPT_PREFIX = ("research/lab_notebook/",)


def is_lab_notebook_exempt(changed_files: list[str]) -> bool:
    """True iff every changed path is in the exempt set (and list is non-empty)."""
    if not changed_files:
        return False
    return all(
        p in _EXEMPT_EXACT or p.startswith(_EXEMPT_PREFIX)
        for p in changed_files
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest workflow/tests/test_closure_audit.py -v`
Expected: 17 passed.

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/closure_audit_lib.py workflow/tests/test_closure_audit.py
git commit -m "feat(closure-audit): exemption-path filter (#325)"
```

---

## Task 6: Role resolution (TDD)

**Files:**
- Modify: `workflow/scripts/closure_audit_lib.py`
- Modify: `workflow/tests/test_closure_audit.py`

- [ ] **Step 1: Write failing tests**

Append to `workflow/tests/test_closure_audit.py`:

```python
def test_resolve_roles_single_issue_single_label():
    issues = [{"labels": ["role:developer", "priority:p2"]}]
    assert lib.resolve_roles(issues) == ["developer"]


def test_resolve_roles_multiple_issues_distinct_roles():
    issues = [
        {"labels": ["role:developer"]},
        {"labels": ["role:scientist"]},
    ]
    assert lib.resolve_roles(issues) == ["developer", "scientist"]  # alphabetical


def test_resolve_roles_duplicates_deduped():
    issues = [
        {"labels": ["role:developer"]},
        {"labels": ["role:developer"]},
    ]
    assert lib.resolve_roles(issues) == ["developer"]


def test_resolve_roles_issue_with_multiple_role_labels_first_alphabetical():
    issues = [{"labels": ["role:scientist", "role:developer"]}]
    assert lib.resolve_roles(issues) == ["developer"]


def test_resolve_roles_no_role_label():
    issues = [{"labels": ["priority:p2"]}]
    assert lib.resolve_roles(issues) == []
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest workflow/tests/test_closure_audit.py -v -k resolve_roles`
Expected: 5 failures, `AttributeError: module 'closure_audit_lib' has no attribute 'resolve_roles'`.

- [ ] **Step 3: Implement**

Append to `workflow/scripts/closure_audit_lib.py`:

```python
def resolve_roles(issues: list[dict]) -> list[str]:
    """Pull deduped, alphabetically-sorted role names from issue labels.

    If an issue has multiple role:* labels, use the alphabetically first.
    """
    roles: set[str] = set()
    for issue in issues:
        role_labels = sorted(
            lbl[len("role:"):]
            for lbl in issue.get("labels", [])
            if lbl.startswith("role:")
        )
        if role_labels:
            roles.add(role_labels[0])
    return sorted(roles)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest workflow/tests/test_closure_audit.py -v`
Expected: 22 passed.

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/closure_audit_lib.py workflow/tests/test_closure_audit.py
git commit -m "feat(closure-audit): role resolution from labels (#325)"
```

---

## Task 7: Comment formatter (TDD)

**Files:**
- Modify: `workflow/scripts/closure_audit_lib.py`
- Modify: `workflow/tests/test_closure_audit.py`

- [ ] **Step 1: Write failing tests**

Append to `workflow/tests/test_closure_audit.py`:

```python
COMMENT_MARKER = "<!-- closure-audit -->"


def test_format_comment_no_gaps_returns_clean_state():
    gaps = {"ac": [], "priority_rationale": [], "lab_notebook": []}
    out = lib.format_comment(gaps, event_label="PR #325")
    assert COMMENT_MARKER in out
    assert "All clear" in out


def test_format_comment_with_gaps_lists_each():
    gaps = {
        "ac": [(999, "2/3 unticked, no deferral comment found")],
        "priority_rationale": [(999, "no 'Priority rationale' line found in issue body")],
        "lab_notebook": [("developer", "'## 2026-05-11' block does not reference '#325'")],
    }
    out = lib.format_comment(gaps, event_label="PR #325")
    assert COMMENT_MARKER in out
    assert "PR #325" in out
    assert "Issue #999" in out
    assert "2/3 unticked" in out
    assert "Priority rationale" in out
    assert "developer" in out
    assert "#325" in out


def test_format_comment_partial_gaps_only_lists_failing():
    gaps = {
        "ac": [(999, "2/3 unticked, no deferral comment found")],
        "priority_rationale": [],
        "lab_notebook": [],
    }
    out = lib.format_comment(gaps, event_label="PR #325")
    assert "AC checkboxes" in out
    assert "Priority rationale" not in out
    assert "Lab notebook" not in out
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest workflow/tests/test_closure_audit.py -v -k format_comment`
Expected: 3 failures, `AttributeError: module 'closure_audit_lib' has no attribute 'format_comment'`.

- [ ] **Step 3: Implement**

Append to `workflow/scripts/closure_audit_lib.py`:

```python
COMMENT_MARKER = "<!-- closure-audit -->"


def format_comment(gaps: dict, event_label: str) -> str:
    """Render the audit comment markdown.

    Args:
        gaps: dict with keys 'ac', 'priority_rationale', 'lab_notebook' →
              list of (subject, detail) tuples. Empty list = no gap.
              For 'ac' and 'priority_rationale': subject is issue number (int).
              For 'lab_notebook': subject is role name (str).
        event_label: human-readable closure context, e.g. "PR #325" or "Issue #200".
    """
    has_any = any(gaps.values())
    if not has_any:
        return (
            f"{COMMENT_MARKER}\n"
            f"✅ Closure audit — all clear\n\n"
            f"_— closure-audit bot_\n"
        )

    lines = [
        COMMENT_MARKER,
        f"## Closure audit — gaps found",
        "",
        f"{event_label} closed with the following items still open:",
        "",
    ]
    for issue_num, detail in gaps["ac"]:
        lines.append(f"- ⚠️ **AC checkboxes** on Issue #{issue_num} — {detail}")
    for issue_num, detail in gaps["priority_rationale"]:
        lines.append(f"- ⚠️ **Priority rationale** on Issue #{issue_num} — {detail}")
    for role, detail in gaps["lab_notebook"]:
        lines.append(f"- ⚠️ **Lab notebook entry** for `{role}` — {detail}")
    lines.append("")
    lines.append(
        "Per [closure ritual](shared/feedback_closure_ritual.md), either tick "
        "the boxes, add a comment-deferral, or add the missing entry/rationale. "
        "This comment auto-updates on next merge to a linked issue."
    )
    lines.append("")
    lines.append("_— closure-audit bot_")
    return "\n".join(lines) + "\n"
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest workflow/tests/test_closure_audit.py -v`
Expected: 25 passed.

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/closure_audit_lib.py workflow/tests/test_closure_audit.py
git commit -m "feat(closure-audit): comment formatter (#325)"
```

---

## Task 8: Main orchestration script

**Files:**
- Create: `workflow/scripts/closure_audit.py`
- Modify: `workflow/tests/test_closure_audit.py` (smoke test for arg parsing)

- [ ] **Step 1: Write the main script**

Write `workflow/scripts/closure_audit.py`:

```python
"""Closure-audit critic — orchestration.

Runs on PR-merge and issue-close GitHub Actions events. Fetches state via
`gh` CLI, calls pure-function checks in closure_audit_lib, posts (or edits)
a single marker-tagged comment on the closing PR/issue.

Exit code: always 0 (non-blocking). Failures are reported in the comment.

Usage:
    python closure_audit.py --event-type {pr|issue} --number <N>
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import closure_audit_lib as lib


REPO_ROOT = Path(__file__).resolve().parents[2]


def gh(*args: str) -> str:
    """Run `gh` and return stdout. Raises on non-zero exit."""
    result = subprocess.run(
        ["gh", *args],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout


def fetch_pr(number: int) -> dict:
    out = gh(
        "pr", "view", str(number),
        "--json",
        "number,mergedAt,closingIssuesReferences,files,labels",
    )
    return json.loads(out)


def fetch_issue(number: int) -> dict:
    out = gh(
        "issue", "view", str(number),
        "--json", "number,body,closedAt,labels,comments",
    )
    return json.loads(out)


def fetch_existing_audit_comment_id(target_type: str, number: int) -> int | None:
    """Find a prior closure-audit comment by marker. Returns its DB id or None."""
    out = gh(
        target_type, "view", str(number),
        "--json", "comments",
    )
    comments = json.loads(out).get("comments", [])
    for c in comments:
        if lib.COMMENT_MARKER in c.get("body", ""):
            return c["id"]
    return None


def post_or_update_comment(target_type: str, number: int, body: str) -> None:
    """Post a new comment, or edit the existing marker-tagged one."""
    existing_id = fetch_existing_audit_comment_id(target_type, number)
    if existing_id is None:
        # New comment
        subprocess.run(
            ["gh", target_type, "comment", str(number), "--body", body],
            check=True,
        )
        return
    # gh doesn't support comment-edit directly; use the REST API.
    subprocess.run(
        [
            "gh", "api", "--method", "PATCH",
            f"/repos/{{owner}}/{{repo}}/issues/comments/{existing_id}",
            "-f", f"body={body}",
        ],
        check=True,
    )


def run_pr_audit(pr_number: int) -> dict:
    pr = fetch_pr(pr_number)
    linked = pr.get("closingIssuesReferences", []) or []
    changed_files = [f["path"] for f in pr.get("files", [])]
    merge_date = pr["mergedAt"][:10]  # 'YYYY-MM-DD'

    # Fetch each linked issue's full body + comments + labels.
    linked_full = [fetch_issue(li["number"]) for li in linked]

    gaps = {"ac": [], "priority_rationale": [], "lab_notebook": []}

    for issue in linked_full:
        body = issue.get("body", "") or ""
        comments_bodies = [c.get("body", "") for c in issue.get("comments", [])]

        ac_gap, ac_detail = lib.check_ac_boxes(body, comments_bodies)
        if ac_gap:
            gaps["ac"].append((issue["number"], ac_detail))

        pr_gap, pr_detail = lib.check_priority_rationale(body)
        if pr_gap:
            gaps["priority_rationale"].append((issue["number"], pr_detail))

    # Lab notebook check (per distinct role)
    if not lib.is_lab_notebook_exempt(changed_files):
        normalized_issues = [
            {"labels": [lbl["name"] for lbl in i.get("labels", [])]}
            for i in linked_full
        ]
        roles = lib.resolve_roles(normalized_issues)
        for role in roles:
            notebook_path = REPO_ROOT / "research" / "lab_notebook" / f"{role}.md"
            if not notebook_path.exists():
                gaps["lab_notebook"].append(
                    (role, f"lab notebook file missing for role `{role}`")
                )
                continue
            text = notebook_path.read_text()
            nb_gap, nb_detail = lib.check_lab_notebook(text, merge_date, pr_number)
            if nb_gap:
                gaps["lab_notebook"].append((role, nb_detail))

    return gaps


def run_issue_audit(issue_number: int) -> dict:
    issue = fetch_issue(issue_number)
    body = issue.get("body", "") or ""
    close_date = issue["closedAt"][:10]
    comments_bodies = [c.get("body", "") for c in issue.get("comments", [])]

    gaps = {"ac": [], "priority_rationale": [], "lab_notebook": []}

    ac_gap, ac_detail = lib.check_ac_boxes(body, comments_bodies)
    if ac_gap:
        gaps["ac"].append((issue_number, ac_detail))

    pr_gap, pr_detail = lib.check_priority_rationale(body)
    if pr_gap:
        gaps["priority_rationale"].append((issue_number, pr_detail))

    # No PR → no path-based exemption applies; always run notebook check.
    normalized = [{"labels": [lbl["name"] for lbl in issue.get("labels", [])]}]
    roles = lib.resolve_roles(normalized)
    for role in roles:
        notebook_path = REPO_ROOT / "research" / "lab_notebook" / f"{role}.md"
        if not notebook_path.exists():
            gaps["lab_notebook"].append(
                (role, f"lab notebook file missing for role `{role}`")
            )
            continue
        text = notebook_path.read_text()
        nb_gap, nb_detail = lib.check_lab_notebook(text, close_date, issue_number)
        if nb_gap:
            gaps["lab_notebook"].append((role, nb_detail))

    return gaps


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--event-type", choices=["pr", "issue"], required=True)
    parser.add_argument("--number", type=int, required=True)
    args = parser.parse_args()

    if args.event_type == "pr":
        gaps = run_pr_audit(args.number)
        event_label = f"PR #{args.number}"
        target_type = "pr"
    else:
        gaps = run_issue_audit(args.number)
        event_label = f"Issue #{args.number}"
        target_type = "issue"

    body = lib.format_comment(gaps, event_label=event_label)
    post_or_update_comment(target_type, args.number, body)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 2: Write a smoke test for argparse**

Append to `workflow/tests/test_closure_audit.py`:

```python
import subprocess
import sys

SCRIPT = Path(__file__).parent.parent / "scripts" / "closure_audit.py"


def test_closure_audit_script_help_runs():
    """Smoke test — script's argparse help works (no import errors)."""
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--help"],
        capture_output=True,
        text=True,
        timeout=10,
    )
    assert result.returncode == 0
    assert "--event-type" in result.stdout
    assert "--number" in result.stdout
```

- [ ] **Step 3: Run all tests to verify they pass**

Run: `pytest workflow/tests/test_closure_audit.py -v`
Expected: 26 passed.

- [ ] **Step 4: Commit**

```bash
git add workflow/scripts/closure_audit.py workflow/tests/test_closure_audit.py
git commit -m "feat(closure-audit): main orchestration script (#325)"
```

---

## Task 9: GitHub Actions workflow

**Files:**
- Create: `.github/workflows/closure-audit.yml`

- [ ] **Step 1: Write the workflow**

Write `.github/workflows/closure-audit.yml`:

```yaml
name: Closure Audit

on:
  pull_request:
    types: [closed]
  issues:
    types: [closed]

permissions:
  contents: read
  issues: write
  pull-requests: write

jobs:
  audit-pr:
    if: github.event_name == 'pull_request' && github.event.pull_request.merged == true
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v6
      - uses: actions/setup-python@v6
        with:
          python-version: "3.11"
      - name: Run closure audit
        env:
          GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}
        run: |
          python workflow/scripts/closure_audit.py \
            --event-type pr \
            --number ${{ github.event.pull_request.number }}

  audit-issue:
    if: github.event_name == 'issues'
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v6
      - uses: actions/setup-python@v6
        with:
          python-version: "3.11"
      - name: Run closure audit
        env:
          GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}
        run: |
          python workflow/scripts/closure_audit.py \
            --event-type issue \
            --number ${{ github.event.issue.number }}
```

- [ ] **Step 2: Lint the YAML (basic syntax check)**

Run: `python -c "import yaml; yaml.safe_load(open('.github/workflows/closure-audit.yml'))"`
Expected: no output, exit 0.

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/closure-audit.yml
git commit -m "feat(closure-audit): GitHub Actions workflow (#325)"
```

---

## Task 10: Open PR with Test plan, request review

**Files:**
- None (PR ops only)

- [ ] **Step 1: Push branch**

```bash
git push -u origin chore/developer/issue-325-closure-audit-critic
```

- [ ] **Step 2: Create PR**

```bash
gh pr create \
  --title "feat(ci): closure-audit critic — auto-audit AC + lab notebook + Priority rationale (#325)" \
  --body "$(cat <<'EOF'
**Created by:** Developer

Closes #325.

## Summary

Adds a post-merge / post-close GitHub Actions workflow that automates the closure-ritual checks (AC ticked, lab notebook entry exists, Priority rationale present). Posts an idempotent marker-tagged comment on the closing PR or issue listing any gaps; silent on the clean path.

Spec: [docs/superpowers/specs/2026-05-11-post-merge-critic-design.md](docs/superpowers/specs/2026-05-11-post-merge-critic-design.md)
Plan: [docs/superpowers/plans/2026-05-11-closure-audit-critic.md](docs/superpowers/plans/2026-05-11-closure-audit-critic.md)

## Test plan

- [ ] All new tests pass (\`pytest workflow/tests/test_closure_audit.py -v\` — 26 tests)
- [ ] Full pytest suite still passes (\`pytest workflow/tests/ -v\`)
- [ ] YAML lint: \`python -c "import yaml; yaml.safe_load(open('.github/workflows/closure-audit.yml'))"\` returns clean
- [ ] (post-merge) Smoke test: this PR's own merge triggers the audit and posts a comment; verify expected behavior matches spec (silent if all 3 checks pass, comment if any gap)

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)" \
  --project "JH M Lee Lab" \
  --base main
```

- [ ] **Step 3: Move PR Status to "Ready for review"**

Look up the PR's project item ID:

```bash
PR_NUMBER=$(gh pr view --json number -q .number)
ITEM_ID=$(gh api graphql -f query="query {
  repository(owner: \"Jin-HoMLee\", name: \"splice-neoepitope-pipeline\") {
    pullRequest(number: $PR_NUMBER) {
      projectItems(first: 5) { nodes { id project { number } } }
    }
  }
}" | jq -r '.data.repository.pullRequest.projectItems.nodes[] | select(.project.number == 9) | .id')
```

Set status:

```bash
gh api graphql -f query="mutation {
  updateProjectV2ItemFieldValue(input: {
    projectId: \"PVT_kwHOB17eGc4BSomP\",
    itemId: \"$ITEM_ID\",
    fieldId: \"PVTSSF_lAHOB17eGc4BSomPzhAHFf8\",
    value: { singleSelectOptionId: \"8bf9192f\" }
  }) { projectV2Item { id } }
}"
```

- [ ] **Step 4: Wait for CI green**

```bash
gh pr checks --watch
```

Expected: `pipeline-pytest` and `pipeline-snakemake-dry-run` both green.

- [ ] **Step 5: Request bot review**

```bash
gh pr comment --body "@claude please review"
```

---

## Self-review checklist (executor: tick on completion)

- [ ] All 26 new pytest tests pass
- [ ] No new linter warnings
- [ ] Comment marker `<!-- closure-audit -->` is identical in lib.py and tests
- [ ] Spec sections all covered:
  - [ ] Architecture (Tasks 8 + 9)
  - [ ] AC check (Task 2)
  - [ ] Lab notebook check (Task 4)
  - [ ] Priority rationale check (Task 3)
  - [ ] Exemption logic (Task 5)
  - [ ] Role resolution (Task 6)
  - [ ] Comment format (Task 7)
  - [ ] Idempotency (Task 8: `post_or_update_comment`)
- [ ] Test plan boxes ticked before merge (per closure ritual)
