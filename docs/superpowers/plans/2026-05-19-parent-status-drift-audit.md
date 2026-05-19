# Parent-vs-Children Status Drift Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a parent-vs-children Status drift audit mechanism on project board #9, sibling to the existing milestone-capacity recheck hook.

**Architecture:** New standalone script `scripts/pm/recheck_parent_status.py` with `--issue N` / `--all` modes. Existing `.claude/hooks/recheck_milestone_dispatch.py` is renamed to `recheck_dispatch.py` and extended to dispatch parent-status rechecks on Status field mutations + `gh issue close`. Output emitted as PostToolUse `additionalContext`.

**Tech Stack:** Python 3.10+, `gh` CLI, GitHub REST + GraphQL APIs. Tests via pytest + `unittest.mock`.

**Spec reference:** [`docs/superpowers/specs/2026-05-19-parent-status-drift-audit-design.md`](../specs/2026-05-19-parent-status-drift-audit-design.md)

---

## File Structure

| Path | Action | Responsibility |
|---|---|---|
| `scripts/pm/recheck_parent_status.py` | CREATE | Recheck logic; CLI with `--issue N` / `--all`; exit codes 0/1/2 |
| `tools/ci/test_recheck_parent_status.py` | CREATE | Unit tests for ladder, collective state, drift detection, parent walk |
| `.claude/hooks/recheck_milestone_dispatch.py` | RENAME → `recheck_dispatch.py` | Multiplexer for both milestone + parent-status rechecks |
| `.claude/hooks/recheck_dispatch.py` | MODIFY (after rename) | Add `STATUS_FIELD_ID` watch + invoke parent-status recheck on `gh issue close` |
| `.claude/settings.local.json` | MODIFY | Update hook command path |
| `tools/ci/test_recheck_dispatch.py` | CREATE | Integration test: pipe PostToolUse JSON, assert dispatcher emits expected blocks |

**Constants used across files:**
- Status field ID: `PVTSSF_lAHOB17eGc4BSomPzhAHFf8`
- Status option IDs: Backlog=`f75ad846`, Ready=`61e4505c`, In progress=`47fc9ee4`, Ready for review=`8bf9192f`, In review=`df73e18b`, Done=`98236657`
- Project number: 9
- Repo: `Jin-HoMLee/splice-neoepitope-pipeline`

---

## Task 1: Status precedence ladder + rank function

**Files:**
- Create: `scripts/pm/recheck_parent_status.py`
- Create: `tools/ci/test_recheck_parent_status.py`

- [ ] **Step 1: Write the failing test**

Create `tools/ci/test_recheck_parent_status.py`:

```python
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
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tools/ci/test_recheck_parent_status.py::TestStatusLadder -v
```
Expected: FAIL with `ModuleNotFoundError: No module named 'recheck_parent_status'`

- [ ] **Step 3: Write minimal implementation**

Create `scripts/pm/recheck_parent_status.py`:

```python
#!/usr/bin/env python3
"""Recheck parent-vs-children Status drift on project board #9.

Implements the audit rule from docs/superpowers/specs/2026-05-19-parent-status-drift-audit-design.md.

Usage:
  scripts/pm/recheck_parent_status.py --issue N      # recheck this issue's parent chain
  scripts/pm/recheck_parent_status.py --all          # audit all parent issues on project #9

Exits 0 if no drift, 2 if drift detected, 1 on error.
"""
from __future__ import annotations

STATUS_LADDER = {
    "Backlog": 0,
    "Ready": 1,
    "In progress": 2,
    "Ready for review": 3,
    "In review": 4,
    "Done": 5,
}


def rank(status: str | None) -> int:
    """Map a Status name to its precedence rank. None/unknown → 0 (Backlog)."""
    return STATUS_LADDER.get(status or "", 0)
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest tools/ci/test_recheck_parent_status.py::TestStatusLadder -v
```
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add scripts/pm/recheck_parent_status.py tools/ci/test_recheck_parent_status.py
git commit -m "feat(pm): #407 recheck_parent_status — Status precedence ladder

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 2: Collective children state

**Files:**
- Modify: `scripts/pm/recheck_parent_status.py`
- Modify: `tools/ci/test_recheck_parent_status.py`

- [ ] **Step 1: Write the failing test**

Append to `tools/ci/test_recheck_parent_status.py`:

```python
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
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tools/ci/test_recheck_parent_status.py::TestCollectiveState -v
```
Expected: FAIL with `AttributeError: module 'recheck_parent_status' has no attribute 'collective_state'`

- [ ] **Step 3: Write minimal implementation**

Append to `scripts/pm/recheck_parent_status.py`:

```python
LADDER_INVERSE = {v: k for k, v in STATUS_LADDER.items()}


def collective_state(open_children: list[dict]) -> str:
    """Max-rank status across open children. Empty list → 'Done'."""
    if not open_children:
        return "Done"
    max_rank = max(rank(c.get("status")) for c in open_children)
    return LADDER_INVERSE[max_rank]
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest tools/ci/test_recheck_parent_status.py::TestCollectiveState -v
```
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add scripts/pm/recheck_parent_status.py tools/ci/test_recheck_parent_status.py
git commit -m "feat(pm): #407 recheck_parent_status — collective children state

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 3: Drift classification (3 classes)

**Files:**
- Modify: `scripts/pm/recheck_parent_status.py`
- Modify: `tools/ci/test_recheck_parent_status.py`

- [ ] **Step 1: Write the failing test**

Append to `tools/ci/test_recheck_parent_status.py`:

```python
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
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tools/ci/test_recheck_parent_status.py::TestDriftClassification -v
```
Expected: FAIL with `AttributeError: ... no attribute 'classify_drift'`

- [ ] **Step 3: Write minimal implementation**

Append to `scripts/pm/recheck_parent_status.py`:

```python
def classify_drift(parent_status: str | None, open_children: list[dict]) -> str | None:
    """Classify drift for a parent vs its open children.

    Returns one of: 'FORWARD DRIFT', 'BACKWARD DRIFT', 'COMPLETION DRIFT', or None.
    """
    p_rank = rank(parent_status)
    if not open_children:
        # All children closed; parent should be Done
        return None if p_rank == STATUS_LADDER["Done"] else "COMPLETION DRIFT"
    c_rank = rank(collective_state(open_children))
    if p_rank > c_rank:
        return "FORWARD DRIFT"
    if p_rank < c_rank:
        return "BACKWARD DRIFT"
    return None
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest tools/ci/test_recheck_parent_status.py::TestDriftClassification -v
```
Expected: PASS (6 tests)

- [ ] **Step 5: Commit**

```bash
git add scripts/pm/recheck_parent_status.py tools/ci/test_recheck_parent_status.py
git commit -m "feat(pm): #407 recheck_parent_status — drift classification

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 4: gh helpers — fetch issue, sub-issues, project Status

**Files:**
- Modify: `scripts/pm/recheck_parent_status.py`
- Modify: `tools/ci/test_recheck_parent_status.py`

- [ ] **Step 1: Write the failing test**

Append to `tools/ci/test_recheck_parent_status.py`:

```python
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
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tools/ci/test_recheck_parent_status.py::TestGhHelpers -v
```
Expected: FAIL with `AttributeError: ... no attribute 'parent_issue_number'`

- [ ] **Step 3: Write minimal implementation**

Append to `scripts/pm/recheck_parent_status.py`:

```python
import json
import subprocess

REPO = "Jin-HoMLee/splice-neoepitope-pipeline"
PROJECT_NUMBER = 9


def gh(*args: str, parse_json: bool = True):
    """Invoke `gh` and parse JSON output (or return raw text)."""
    result = subprocess.run(["gh", *args], capture_output=True, text=True, check=True)
    return json.loads(result.stdout) if parse_json else result.stdout


def parent_issue_number(issue_number: int) -> int | None:
    """Return parent issue number via REST parent_issue_url, or None."""
    data = gh("api", f"repos/{REPO}/issues/{issue_number}")
    url = data.get("parent_issue_url")
    if not url:
        return None
    # URL ends with /issues/{N}
    return int(url.rstrip("/").rsplit("/", 1)[-1])


def open_sub_issues(issue_number: int) -> list[dict]:
    """Return list of open sub-issues (each a dict with 'number' at minimum)."""
    data = gh("api", f"repos/{REPO}/issues/{issue_number}/sub_issues")
    return [c for c in data if c.get("state") == "open"]
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest tools/ci/test_recheck_parent_status.py::TestGhHelpers -v
```
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add scripts/pm/recheck_parent_status.py tools/ci/test_recheck_parent_status.py
git commit -m "feat(pm): #407 recheck_parent_status — gh helpers

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 5: Project Status lookup for an issue

**Files:**
- Modify: `scripts/pm/recheck_parent_status.py`
- Modify: `tools/ci/test_recheck_parent_status.py`

- [ ] **Step 1: Write the failing test**

Append to `tools/ci/test_recheck_parent_status.py`:

```python
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
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tools/ci/test_recheck_parent_status.py::TestStatusLookup -v
```
Expected: FAIL with `AttributeError: ... no attribute 'status_for_issue'`

- [ ] **Step 3: Write minimal implementation**

Append to `scripts/pm/recheck_parent_status.py`:

```python
def status_for_issue(issue_number: int) -> str | None:
    """Return the project #9 Status name for issue, or None if not on the project."""
    owner, name = REPO.split("/")
    query = (
        f'query {{ repository(owner: "{owner}", name: "{name}") {{ '
        f'issue(number: {issue_number}) {{ '
        f'projectItems(first: 5) {{ nodes {{ '
        f'project {{ number }} '
        f'fieldValues(first: 20) {{ nodes {{ '
        f'... on ProjectV2ItemFieldSingleSelectValue {{ '
        f'name field {{ ... on ProjectV2SingleSelectField {{ name }} }} '
        f'}} }} }} '
        f'}} }} }} }} }}'
    )
    data = gh("api", "graphql", "-f", f"query={query}")
    nodes = (
        data.get("data", {})
        .get("repository", {})
        .get("issue", {})
        .get("projectItems", {})
        .get("nodes", [])
    ) or []
    for pi in nodes:
        if (pi.get("project") or {}).get("number") != PROJECT_NUMBER:
            continue
        for fv in pi.get("fieldValues", {}).get("nodes", []):
            if (fv.get("field") or {}).get("name") == "Status":
                return fv.get("name")
    return None
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest tools/ci/test_recheck_parent_status.py::TestStatusLookup -v
```
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add scripts/pm/recheck_parent_status.py tools/ci/test_recheck_parent_status.py
git commit -m "feat(pm): #407 recheck_parent_status — project Status lookup

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 6: Per-issue audit + parent walk

**Files:**
- Modify: `scripts/pm/recheck_parent_status.py`
- Modify: `tools/ci/test_recheck_parent_status.py`

- [ ] **Step 1: Write the failing test**

Append to `tools/ci/test_recheck_parent_status.py`:

```python
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
        assert chain[1]["drift"] == "FORWARD DRIFT"

    @patch("recheck_parent_status.status_for_issue")
    @patch("recheck_parent_status.open_sub_issues")
    @patch("recheck_parent_status.parent_issue_number")
    def test_stops_at_root(self, mock_parent, mock_subs, mock_status):
        mock_parent.side_effect = lambda n: {100: None}[n]
        # Issue 100 has no parent → empty chain
        chain = rps.audit_parent_chain(100)
        assert chain == []
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tools/ci/test_recheck_parent_status.py::TestAuditChain -v
```
Expected: FAIL with `AttributeError: ... no attribute 'audit_parent_chain'`

- [ ] **Step 3: Write minimal implementation**

Append to `scripts/pm/recheck_parent_status.py`:

```python
def audit_parent_chain(issue_number: int) -> list[dict]:
    """Walk up the parent chain from issue_number; audit drift at each level.

    Returns list of records: {issue, status, open_children, collective, drift}.
    Empty list if issue has no parent.
    """
    chain: list[dict] = []
    cursor = parent_issue_number(issue_number)
    seen = {issue_number}
    while cursor is not None and cursor not in seen:
        seen.add(cursor)
        parent_status = status_for_issue(cursor)
        children = open_sub_issues(cursor)
        # Enrich children with their Status
        enriched = [{"number": c["number"], "status": status_for_issue(c["number"])}
                    for c in children]
        chain.append({
            "issue": cursor,
            "status": parent_status,
            "open_children": enriched,
            "collective": collective_state(enriched),
            "drift": classify_drift(parent_status, enriched),
        })
        cursor = parent_issue_number(cursor)
    return chain
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest tools/ci/test_recheck_parent_status.py::TestAuditChain -v
```
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add scripts/pm/recheck_parent_status.py tools/ci/test_recheck_parent_status.py
git commit -m "feat(pm): #407 recheck_parent_status — parent walk audit

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 7: Output formatter

**Files:**
- Modify: `scripts/pm/recheck_parent_status.py`
- Modify: `tools/ci/test_recheck_parent_status.py`

- [ ] **Step 1: Write the failing test**

Append to `tools/ci/test_recheck_parent_status.py`:

```python
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
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tools/ci/test_recheck_parent_status.py::TestFormatter -v
```
Expected: FAIL with `AttributeError: ... no attribute 'format_record'`

- [ ] **Step 3: Write minimal implementation**

Append to `scripts/pm/recheck_parent_status.py`:

```python
def format_record(record: dict) -> str:
    """Render one audit record as a multi-line block matching the spec format."""
    lines: list[str] = []
    lines.append(f"#{record['issue']} — Status: {record['status'] or 'NO STATUS'}")
    children = record["open_children"]
    lines.append(f"  Open sub-issues ({len(children)}):")
    for c in children:
        lines.append(f"    - #{c['number']} ({c['status'] or 'NO STATUS'})")
    lines.append(f"  Collective children state: {record['collective']}")
    drift = record["drift"]
    if drift is None:
        lines.append(f"  Status: [No change]")
    else:
        lines.append(f"  Status: [{drift}]")
    return "\n".join(lines)
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest tools/ci/test_recheck_parent_status.py::TestFormatter -v
```
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add scripts/pm/recheck_parent_status.py tools/ci/test_recheck_parent_status.py
git commit -m "feat(pm): #407 recheck_parent_status — output formatter

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 8: CLI — `--issue N` mode

**Files:**
- Modify: `scripts/pm/recheck_parent_status.py`

- [ ] **Step 1: Write the failing test**

Append to `tools/ci/test_recheck_parent_status.py`:

```python
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
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tools/ci/test_recheck_parent_status.py::TestCli -v
```
Expected: FAIL with `AttributeError: ... no attribute 'main'`

- [ ] **Step 3: Write minimal implementation**

Append to `scripts/pm/recheck_parent_status.py`:

```python
import argparse
import sys


def run_issue_mode(issue_number: int) -> int:
    chain = audit_parent_chain(issue_number)
    if not chain:
        print(f"Issue #{issue_number} has no parent — nothing to audit.")
        return 0
    print(f"Parent chain for #{issue_number} (walked {len(chain)} levels):\n")
    drifted = False
    for record in chain:
        print(format_record(record))
        print()
        if record["drift"] is not None:
            drifted = True
    return 2 if drifted else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--issue", type=int, help="Issue number; walk its parent chain")
    group.add_argument("--all", action="store_true", help="Audit all parent issues on project #9")
    args = parser.parse_args(argv)

    if args.issue:
        return run_issue_mode(args.issue)
    return run_all_mode()


def run_all_mode() -> int:
    # Implemented in Task 9
    print("--all not yet implemented", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest tools/ci/test_recheck_parent_status.py::TestCli -v
```
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add scripts/pm/recheck_parent_status.py tools/ci/test_recheck_parent_status.py
git commit -m "feat(pm): #407 recheck_parent_status — --issue CLI mode

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 9: CLI — `--all` mode

**Files:**
- Modify: `scripts/pm/recheck_parent_status.py`
- Modify: `tools/ci/test_recheck_parent_status.py`

- [ ] **Step 1: Write the failing test**

Append to `tools/ci/test_recheck_parent_status.py`:

```python
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
        mock_status.side_effect = lambda n: {
            86: "In progress", 24: "Ready",
            204: "Backlog", 86: "In progress",
        }.get(n) if n != 86 else "In progress"  # parent #86's own Status

        # Re-define the side_effect more cleanly:
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
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tools/ci/test_recheck_parent_status.py::TestAllMode -v
```
Expected: FAIL with `AttributeError: ... no attribute 'all_parent_issues'`

- [ ] **Step 3: Write minimal implementation**

Replace the stub `run_all_mode` in `scripts/pm/recheck_parent_status.py` and add `all_parent_issues`:

```python
def all_parent_issues() -> list[int]:
    """Return numbers of all open issues in the repo that have ≥1 sub-issue."""
    # GitHub REST search: filter open issues, then check sub_issues_summary.total > 0
    # via a per-issue follow-up. Use gh issue list + per-issue REST fetch.
    data = gh("issue", "list", "--repo", REPO, "--state", "open", "--limit", "200",
              "--json", "number")
    parents: list[int] = []
    for issue in data:
        n = issue["number"]
        meta = gh("api", f"repos/{REPO}/issues/{n}")
        if (meta.get("sub_issues_summary") or {}).get("total", 0) > 0:
            parents.append(n)
    return parents


def run_all_mode() -> int:
    parents = all_parent_issues()
    drifted_count = 0
    drift_blocks: list[str] = []
    for p in parents:
        parent_status = status_for_issue(p)
        children = open_sub_issues(p)
        enriched = [{"number": c["number"], "status": status_for_issue(c["number"])}
                    for c in children]
        record = {
            "issue": p,
            "status": parent_status,
            "open_children": enriched,
            "collective": collective_state(enriched),
            "drift": classify_drift(parent_status, enriched),
        }
        if record["drift"] is not None:
            drifted_count += 1
            drift_blocks.append(format_record(record))

    print(f"Audited {len(parents)} parent issues; {drifted_count} drifted.\n")
    for block in drift_blocks:
        print(block)
        print()
    return 2 if drifted_count > 0 else 0
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest tools/ci/test_recheck_parent_status.py::TestAllMode -v
```
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add scripts/pm/recheck_parent_status.py tools/ci/test_recheck_parent_status.py
git commit -m "feat(pm): #407 recheck_parent_status — --all CLI mode

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 10: Run full test suite + live `--all` smoke test

- [ ] **Step 1: Run full unit test suite**

```bash
pytest tools/ci/test_recheck_parent_status.py -v
```
Expected: ALL PASS

- [ ] **Step 2: Live `--all` smoke test against project #9**

```bash
python3 scripts/pm/recheck_parent_status.py --all
```
Expected output (after this morning's in-session flips of #24/#86/#126 → Ready):
- `Audited N parent issues; 0 drifted.` — no drift records emitted
- Exit code: 0

If drift IS surfaced, investigate immediately — either a real new drift or a bug in the audit logic.

- [ ] **Step 3: Live `--issue N` smoke test on a known-clean issue**

```bash
python3 scripts/pm/recheck_parent_status.py --issue 204
```
Expected output:
- Walks #204 → #86 → #24 (2-level chain)
- Both #86 and #24 show `[No change]` (Status: Ready, children's collective: Ready or Backlog at most)
- Exit code: 0

If `[FORWARD DRIFT]` appears, real drift slipped back in — investigate.

- [ ] **Step 4: No commit needed — these are verification steps only**

---

## Task 11: Rename dispatcher hook + add Status field trigger

**Files:**
- Move: `.claude/hooks/recheck_milestone_dispatch.py` → `.claude/hooks/recheck_dispatch.py`
- Modify: `.claude/hooks/recheck_dispatch.py` (add Status field watch + parent-status dispatch)

- [ ] **Step 1: Git-rename the file**

```bash
git mv .claude/hooks/recheck_milestone_dispatch.py .claude/hooks/recheck_dispatch.py
```

Verify:
```bash
git status
```
Expected: `renamed: .claude/hooks/recheck_milestone_dispatch.py -> .claude/hooks/recheck_dispatch.py`

- [ ] **Step 2: Update docstring + add new constants**

Edit `.claude/hooks/recheck_dispatch.py`. Replace the module docstring (lines 1-12) with:

```python
#!/usr/bin/env python3
"""PostToolUse hook: trigger milestone-capacity AND parent-status rechecks on capacity-change events.

Watches Bash commands for trigger shapes:
  1. gh issue close N                              -> recheck issue's milestone + parent-status chain
  2. gh issue edit N ... --milestone X             -> recheck old + new milestones (via /events history)
  3. project board Size mutation (Size field ID)   -> recheck affected issue's milestone
  4. gh api .../milestones/N PATCH (due_on edit)   -> recheck milestone N
  5. project board Status mutation (Status field ID) -> recheck affected issue's parent-status chain

For each match, invokes the relevant script and emits its output as
`additionalContext` so it appears in the next prompt.
"""
```

Add a new constant **immediately after** the existing `SIZE_FIELD_ID` line (currently line 22):

```python
STATUS_FIELD_ID = "PVTSSF_lAHOB17eGc4BSomPzhAHFf8"
PARENT_STATUS_SCRIPT = str(Path(__file__).resolve().parent.parent.parent / "scripts" / "pm" / "recheck_parent_status.py")
```

(For reference: `SIZE_FIELD_ID` is `PVTSSF_lAHOB17eGc4BSomPzhAHGiA`; `STATUS_FIELD_ID` is `PVTSSF_lAHOB17eGc4BSomPzhAHFf8`. Distinct values — last 5 chars differ.)

- [ ] **Step 3: Extend `dispatch()` to invoke parent-status on Status mutations**

In the `dispatch()` function, after the existing `if SIZE_FIELD_ID in cmd:` block, add:

```python
    if STATUS_FIELD_ID in cmd:
        item_match = PATTERN_ITEMID.search(cmd)
        if item_match:
            issue = lookup_issue_for_item(item_match.group(1))
            if issue:
                outputs.append(
                    f"[parent-status recheck — Status change on #{issue}]\n"
                    f"{run_parent_status_recheck('--issue', str(issue))}"
                )
```

- [ ] **Step 4: Extend `gh issue close` trigger to also recheck parent-status**

In the `PATTERN_CLOSE` block, after the milestone recheck append, add a parent-status append:

```python
    m = PATTERN_CLOSE.search(cmd)
    if m:
        outputs.append(
            f"[milestone recheck — close #{m.group(1)}]\n{run_recheck('--issue', m.group(1))}"
        )
        outputs.append(
            f"[parent-status recheck — close #{m.group(1)}]\n"
            f"{run_parent_status_recheck('--issue', m.group(1))}"
        )
```

- [ ] **Step 5: Add `run_parent_status_recheck` helper**

After the existing `run_recheck` function, add:

```python
def run_parent_status_recheck(*args: str) -> str:
    if not Path(PARENT_STATUS_SCRIPT).is_file():
        return f"(parent-status recheck error: script not found at {PARENT_STATUS_SCRIPT})"
    try:
        result = subprocess.run(
            ["python3", PARENT_STATUS_SCRIPT, *args],
            capture_output=True, text=True, timeout=30, check=False,
        )
    except (subprocess.TimeoutExpired, FileNotFoundError) as exc:
        return f"(parent-status recheck error: {exc})"
    out = result.stdout
    if result.stderr:
        out += result.stderr
    return out
```

- [ ] **Step 6: Smoke-check the file parses cleanly**

```bash
python3 -c "import ast; ast.parse(open('.claude/hooks/recheck_dispatch.py').read()); print('OK')"
```
Expected: `OK`

- [ ] **Step 7: Commit**

```bash
git add .claude/hooks/recheck_dispatch.py
git commit -m "feat(pm): #407 rename dispatcher + add Status field trigger

Rename .claude/hooks/recheck_milestone_dispatch.py → recheck_dispatch.py;
extend to multiplex milestone-capacity AND parent-status rechecks.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 12: Update `.claude/settings.local.json` hook path

**Files:**
- Modify: `.claude/settings.local.json`

- [ ] **Step 1: Read current hook config**

```bash
grep -A 3 "recheck_milestone_dispatch" .claude/settings.local.json
```
Expected: shows the existing PostToolUse hook entry referencing the old path.

- [ ] **Step 2: Update the path**

Use the Edit tool to change `recheck_milestone_dispatch.py` to `recheck_dispatch.py` in `.claude/settings.local.json`.

- [ ] **Step 3: Verify**

```bash
grep "recheck.*dispatch\.py" .claude/settings.local.json
```
Expected: shows `recheck_dispatch.py` (no more `recheck_milestone_dispatch.py`).

```bash
python3 -c "import json; json.load(open('.claude/settings.local.json'))" && echo "JSON valid"
```
Expected: `JSON valid`

- [ ] **Step 4: Commit**

```bash
git add .claude/settings.local.json
git commit -m "feat(pm): #407 settings — point hook to renamed recheck_dispatch.py

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 13: Hook integration smoke test

**Files:**
- Create: `tools/ci/test_recheck_dispatch.py`

- [ ] **Step 1: Write the test**

Create `tools/ci/test_recheck_dispatch.py`:

```python
"""Integration tests for `.claude/hooks/recheck_dispatch.py`.

Pipes a synthetic PostToolUse JSON payload to the hook and asserts that
the additionalContext includes the expected `[<recheck-name> — ...]` blocks.

These tests SHELL OUT TO `gh` via the scripts they dispatch — they will hit
the live GitHub API. Use only against benign queries that won't mutate state.
"""

import json
import subprocess
import sys
from pathlib import Path

HOOK = Path(__file__).parent.parent.parent / ".claude" / "hooks" / "recheck_dispatch.py"


def _run(cmd: str) -> tuple[int, str, str]:
    payload = json.dumps({"tool_input": {"command": cmd}})
    proc = subprocess.run(
        [sys.executable, str(HOOK)],
        input=payload,
        capture_output=True,
        text=True,
        timeout=60,
    )
    return proc.returncode, proc.stdout, proc.stderr


def _additional_context(stdout: str) -> str:
    if not stdout.strip():
        return ""
    return json.loads(stdout).get("hookSpecificOutput", {}).get("additionalContext", "")


class TestDispatcherSilent:
    def test_non_gh_command_silent(self):
        rc, out, _ = _run("ls -la")
        assert rc == 0
        assert out.strip() == ""

    def test_non_matching_gh_command_silent(self):
        rc, out, _ = _run("gh pr status")
        assert rc == 0
        assert out.strip() == ""


class TestStatusFieldTrigger:
    def test_status_mutation_emits_parent_status_block(self):
        # Use a synthetic command that contains both the Status field ID and an item ID.
        # Item ID is for issue #24 (known, audited 2026-05-19, currently Ready).
        cmd = (
            'gh api graphql -f query=\'mutation { updateProjectV2ItemFieldValue('
            'input: { projectId: "PVT_kwHOB17eGc4BSomP" '
            'itemId: "PVTI_lAHOB17eGc4BSomPzgpgr30" '
            'fieldId: "PVTSSF_lAHOB17eGc4BSomPzhAHFf8" '
            'value: { singleSelectOptionId: "61e4505c" } }) { projectV2Item { id } } }\''
        )
        rc, out, _ = _run(cmd)
        assert rc == 0
        ctx = _additional_context(out)
        assert "[parent-status recheck — Status change on #24]" in ctx
```

- [ ] **Step 2: Run the tests**

```bash
pytest tools/ci/test_recheck_dispatch.py -v
```
Expected: 3 PASS (uses live `gh` calls; ensure GitHub auth works).

If the third test takes >30s, the parent-walk for #24 may be slow — check that `gh` auth is cached (`gh auth status`).

- [ ] **Step 3: Commit**

```bash
git add tools/ci/test_recheck_dispatch.py
git commit -m "test(pm): #407 dispatcher integration smoke tests

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 14: Lab notebook entry

**Files:**
- Modify: `research/lab_notebook/pm.md`

- [ ] **Step 1: Check current pm.md structure**

```bash
head -30 research/lab_notebook/pm.md
```
Note the existing date sections and time-section format.

- [ ] **Step 2: Add a new time section under 2026-05-19**

Use the Edit tool to add a new time section at the TOP of the 2026-05-19 date block (lab notebook rule: new entries at TOP). Use today's actual UTC time:

```bash
date -u +"%Y-%m-%d %H:%M UTC"
```

Entry content (substitute the actual time):

```markdown
### HH:MM UTC — Parent-status drift audit mechanism (#407)

Built the sibling recheck mechanism to PR #397. Triggered by today's mid-day board sweep finding 3 epics drifted In progress with all open sub-issues in Backlog (#24, #86, #126 — flipped to Ready in-session).

**Shipped:**
- New `scripts/pm/recheck_parent_status.py` with `--issue N` / `--all` modes — Status precedence ladder, 3 drift classes (forward / backward / completion), parent-chain walk via REST `parent_issue_url`.
- Renamed `.claude/hooks/recheck_milestone_dispatch.py` → `recheck_dispatch.py`; now a multiplexer for milestone + parent-status rechecks.
- Hook fires parent-status recheck on (a) project Status field mutation + (b) `gh issue close`.
- Unit + integration tests in `tools/ci/`. `--all` smoke run against the live board: 0 drifted (clean post in-session flips).

**Mechanism dividend:** while filing #407 itself, the existing milestone-capacity recheck hook caught the pm-i2 capacity drift introduced by adding #407, prompted the due_on PATCH to 2026-07-09, and re-verified clean on the second fire. The new sibling hook will catch the parent-status equivalent automatically going forward.
```

- [ ] **Step 3: Commit**

```bash
git add research/lab_notebook/pm.md
git commit -m "docs(lab_notebook): #407 parent-status drift audit mechanism shipped

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 15: Push branch + open PR

- [ ] **Step 1: Push the branch**

```bash
git push -u origin feat/pm/issue-407-parent-status-drift-audit
```

- [ ] **Step 2: Verify Issue #407 ACs are still appropriate**

```bash
gh issue view 407 --repo Jin-HoMLee/splice-neoepitope-pipeline | head -80
```
Confirm the AC checklist is what was filed; no edits needed unless scope drifted.

- [ ] **Step 3: Open the PR**

```bash
gh pr create \
  --repo Jin-HoMLee/splice-neoepitope-pipeline \
  --base main \
  --head feat/pm/issue-407-parent-status-drift-audit \
  --title "feat(pm): parent-vs-children Status drift audit (#407)" \
  --label "role:pm" \
  --milestone "pm-i2 - PM Self-Improvement Tooling" \
  --project "JH M Lee Lab" \
  --body "$(cat <<'EOF'
## Summary
Sibling mechanism to [PR #397](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/397) milestone-capacity recheck hook. Catches parent-vs-children Status drift on project board #9 — the failure mode discovered during the 2026-05-19 mid-day board sweep where 3 epics ([Issue #24](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/24), [Issue #86](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/86), [Issue #126](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/126)) were In progress while all open sub-issues were Backlog.

- New `scripts/pm/recheck_parent_status.py` with `--issue N` and `--all` modes
- `.claude/hooks/recheck_milestone_dispatch.py` renamed to `recheck_dispatch.py`; multiplexes both rechecks
- Dispatcher triggers parent-status recheck on Status field mutation + `gh issue close`

Spec: [`docs/superpowers/specs/2026-05-19-parent-status-drift-audit-design.md`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/blob/feat/pm/issue-407-parent-status-drift-audit/docs/superpowers/specs/2026-05-19-parent-status-drift-audit-design.md)
Plan: [`docs/superpowers/plans/2026-05-19-parent-status-drift-audit.md`](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/blob/feat/pm/issue-407-parent-status-drift-audit/docs/superpowers/plans/2026-05-19-parent-status-drift-audit.md)

Closes [Issue #407](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/407).

## Test plan
- [ ] Unit tests pass: `pytest tools/ci/test_recheck_parent_status.py -v`
- [ ] Integration tests pass: `pytest tools/ci/test_recheck_dispatch.py -v`
- [ ] Live smoke: `python3 scripts/pm/recheck_parent_status.py --all` reports 0 drifted on current clean board
- [ ] Live smoke: `python3 scripts/pm/recheck_parent_status.py --issue 204` walks #204 → #86 → #24 with `[No change]` on both levels
- [ ] Hook re-arm verification: after merge, the dispatcher fires on a benign Status mutation and emits `[parent-status recheck — Status change on #N]` block
- [ ] CLAUDE.md / memory unchanged (mechanism is additive; no workflow change)

**Created by:** PM
EOF
)"
```

- [ ] **Step 4: Note the PR number returned**

The PR URL will be printed. Note the number for the final closure-ritual merge step.

---

## Task 16: Merge via closure-ritual gate

- [ ] **Step 1: Verify CI green + bot review (if any)**

```bash
gh pr checks <PR_NUMBER> --repo Jin-HoMLee/splice-neoepitope-pipeline
```
Wait for all required checks: `pipeline-pytest`, `pipeline-snakemake-dry-run`. Expected: PASS.

- [ ] **Step 2: Tick the Test plan checkboxes on the PR**

After verifying each item above:

```bash
gh pr edit <PR_NUMBER> --repo Jin-HoMLee/splice-neoepitope-pipeline --body "<UPDATED_BODY_WITH_TICKED_BOXES>"
```

Or use the GitHub web UI to tick the boxes. All 6 Test plan boxes must be `- [x]` before merge.

- [ ] **Step 3: Tick the AC boxes on Issue #407**

The Issue body has 7 AC checkboxes (see #407). Tick each that the implementation satisfied:

```bash
gh issue edit 407 --repo Jin-HoMLee/splice-neoepitope-pipeline --body "<UPDATED_BODY_WITH_TICKED_BOXES>"
```

- [ ] **Step 4: Run the closure-ritual gate**

```bash
bash scripts/audit_and_merge.sh <PR_NUMBER>
```
Expected: clean audit (no `- [ ]` remaining on PR Test plan or any linked Issue's AC), then auto-squash-merges + deletes branch.

If the gate refuses with unticked boxes, find and tick them before re-running.

---

## Self-Review

Done before handoff. Findings:

**1. Spec coverage** — All 7 ACs from Issue #407 map to tasks:
- AC "Audit script implemented" → Tasks 1-9
- AC "Runs against project #9 + report" → Tasks 8-9
- AC "Wired into PM morning routine / hook / CI" → Tasks 11-12 (hook)
- AC "Smoke-tested on known clean state" → Task 10
- AC "CLAUDE.md / memory updated if workflow changes" → Mechanism is additive; no workflow change needed (noted in Task 15 PR body)
- AC "Algorithm unit tests cover 3 drift classes + no-drift + no-open-children" → Task 3
- AC "Hook integration smoke confirms separate blocks" → Task 13

**2. Placeholder scan** — No TBD / TODO / vague phrases. Every step has either exact code or exact commands with expected output.

**3. Type consistency** — Function signatures used in later tasks match earlier definitions:
- `rank()` → returns `int` (Task 1) consumed by `collective_state()` (Task 2)
- `collective_state()` → returns `str` (Task 2) consumed by `classify_drift()` (Task 3)
- `audit_parent_chain()` → returns `list[dict]` with `{issue, status, open_children, collective, drift}` keys, consumed identically by `format_record()` (Task 7) and `run_issue_mode()` (Task 8)
- `parent_issue_number()`, `open_sub_issues()`, `status_for_issue()` all defined in Tasks 4-5 and consumed by `audit_parent_chain()` (Task 6) + `run_all_mode()` (Task 9)

**Field ID values explicitly listed** in Task 11 Step 2 to prevent the mix-up that almost slipped in during planning (Status = `...AHFf8`, Size = `...AHGiA`, differ only in the last 5 chars).

---
