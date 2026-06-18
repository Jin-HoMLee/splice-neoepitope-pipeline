# Prose-Dependency Reconciler Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `scripts/pm/scan_prose_deps.py` — a reconciler that finds prose-stated Issue dependencies ("depends on #M") and wires the corresponding native GitHub `blockedBy` edges — then run it to backfill the existing graph.

**Architecture:** Four thin layers — *fetch* (open issues via `gh issue list`), *parse* (pure regex: body → `(dependent, blocker)` pairs), *reconcile* (diff prose pairs vs the native graph, classify each), *act* (report / wire / exit-code). The parse layer is pure and carries the test coverage; the gh-backed layers are monkeypatched in tests. Mirrors `scripts/pm/recheck_milestone.py` (the `gh()` subprocess helper) and `scripts/board_open_items.py` (the inject-the-fetch test style).

**Tech Stack:** Python 3 (stdlib only: `argparse`, `re`, `subprocess`, `json`), `gh` CLI, pytest (`workflow/tests/.venv`).

**Spec:** `docs/superpowers/specs/2026-06-16-prose-dependency-reconciler-design.md`
**Issue:** [#722](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/722) — branch `design/pm/issue-722-prose-dependency-reconciler` (already checked out).

---

## File Structure

- **Create:** `scripts/pm/scan_prose_deps.py` — the reconciler CLI + library functions.
- **Create:** `workflow/tests/test_scan_prose_deps.py` — pytest (parse layer = bulk; reconcile/act = mocked gh).
- **Modify (personas memory, MM-committed):** `memory/shared/feedback_dependency_tracking.md` — cross-ref the scanner. Flagged for MM, not part of the code PR.

**Module API (locked here so later tasks match):**

```python
REPO = "Jin-HoMLee/splice-neoepitope-pipeline"

# fetch
def gh(*args, parse_json=True): ...                      # subprocess wrapper (mirror recheck_milestone.py)
def fetch_open_issues() -> list[dict]: ...               # [{number, title, body, state}, ...]

# parse (PURE — no I/O)
BLOCKER_RE: "re.Pattern"                                 # matches a blocker phrase + #<digits>
def parse_dependencies(number: int, body: str) -> list[tuple[int, int]]: ...   # [(dependent, blocker)]

# reconcile (gh-backed)
def issue_meta(number: int) -> dict: ...                 # {"state": "open"|"closed", "is_pr": bool}
def native_blockers(number: int) -> set[int]: ...        # existing blockedBy edge targets
def classify(dependent: int, blocker: int, *, blocker_meta: dict, existing: set[int]) -> str: ...
    # -> "needs-wiring" | "already-wired" | "closed-blocker" | "un-wireable-pr" | "self-ref"
def reconcile(pairs: list[tuple[int, int]]) -> list[dict]: ...   # [{dependent, blocker, state, action}]

# act
def render_report(records: list[dict]) -> str: ...
def wire(records: list[dict]) -> list[tuple[int, int]]: ...      # POSTs edges for needs-wiring; returns wired
def main() -> int: ...                                   # argparse: --report/--apply/--check/--issue/--only
```

---

## Task 1: Scaffold + fetch layer

**Files:**
- Create: `scripts/pm/scan_prose_deps.py`
- Test: `workflow/tests/test_scan_prose_deps.py`

- [ ] **Step 1: Write the failing test for the import + REPO constant**

```python
# workflow/tests/test_scan_prose_deps.py
"""Tests for scripts/pm/scan_prose_deps.py — the prose-dependency reconciler (Issue #722).

Parse layer is pure and gets the bulk of coverage; reconcile/act monkeypatch the
gh-backed helpers so nothing touches the network.
"""
import sys
from pathlib import Path

import pytest

_PM_DIR = Path(__file__).resolve().parents[2] / "scripts" / "pm"
sys.path.insert(0, str(_PM_DIR))

import scan_prose_deps as spd  # noqa: E402


def test_repo_constant():
    assert spd.REPO == "Jin-HoMLee/splice-neoepitope-pipeline"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_scan_prose_deps.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'scan_prose_deps'`

- [ ] **Step 3: Create the script skeleton with the fetch layer**

```python
#!/usr/bin/env python3
"""Reconcile prose-stated Issue dependencies into native GitHub blockedBy edges.

Finds bodies saying "depends on #M" / "blocked on #M" etc., diffs them against
the native blockedBy graph, and reports or wires the missing edges.

Usage:
  scripts/pm/scan_prose_deps.py                 # --report (default): drift table
  scripts/pm/scan_prose_deps.py --issue N       # single-issue scope
  scripts/pm/scan_prose_deps.py --check         # exit 2 if any drift
  scripts/pm/scan_prose_deps.py --apply         # wire all needs-wiring records
  scripts/pm/scan_prose_deps.py --apply --only 745 594   # wire only these dependents

Exits 0 clean / 2 drift-present (--check) / 1 on error.
"""
import argparse
import json
import re
import subprocess
import sys

REPO = "Jin-HoMLee/splice-neoepitope-pipeline"


def gh(*args, parse_json=True):
    """Run a gh command; return parsed JSON (default) or raw stdout text."""
    result = subprocess.run(["gh", *args], capture_output=True, text=True, check=True)
    return json.loads(result.stdout) if parse_json else result.stdout


def fetch_open_issues():
    """All open issues with bodies. Uses the issue list (NOT the project board),
    so the board's Done-first pagination trap does not apply."""
    return gh(
        "issue", "list", "--repo", REPO, "--state", "open", "--limit", "1000",
        "--json", "number,title,body,state",
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_scan_prose_deps.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add scripts/pm/scan_prose_deps.py workflow/tests/test_scan_prose_deps.py
git commit -m "feat(pm): scan_prose_deps scaffold + fetch layer (#722)"
```

---

## Task 2: Parse layer (the core — heavy TDD)

**Files:**
- Modify: `scripts/pm/scan_prose_deps.py`
- Test: `workflow/tests/test_scan_prose_deps.py`

- [ ] **Step 1: Write failing tests — allowlist hits, narrative excludes, multi-blocker, self-ref, dedupe**

```python
# append to workflow/tests/test_scan_prose_deps.py

@pytest.mark.parametrize("phrase", [
    "depends on #722", "blocked by #722", "blocked on #722",
    "blocked-by:722", "gated on #722", "requires #722",
])
def test_parse_allowlist_phrases_match(phrase):
    body = f"Some context. This {phrase} to proceed."
    assert spd.parse_dependencies(745, body) == [(745, 722)]


@pytest.mark.parametrize("phrase", [
    "informs #722", "consumes #722", "relates to #722", "related to #722",
    "follow-up to #722", "supersedes #722", "superseded by #722",
    "see #722", "cf. #722", "fixes #722",
])
def test_parse_narrative_phrases_excluded(phrase):
    # Narrative cross-refs are *why*-context, never blockers — must yield nothing.
    body = f"This PR {phrase}."
    assert spd.parse_dependencies(745, body) == []


def test_parse_multiple_blockers_deduped_and_sorted():
    body = "Depends on #708 and blocked by #707. Also depends on #708 again."
    assert spd.parse_dependencies(709, body) == [(709, 707), (709, 708)]


def test_parse_self_reference_dropped():
    body = "This depends on #745 (itself, nonsensically)."
    assert spd.parse_dependencies(745, body) == []


def test_parse_empty_or_none_body():
    assert spd.parse_dependencies(1, "") == []
    assert spd.parse_dependencies(1, None) == []


def test_parse_case_insensitive():
    assert spd.parse_dependencies(2, "DEPENDS ON #5") == [(2, 5)]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_scan_prose_deps.py -k parse -v`
Expected: FAIL — `AttributeError: module 'scan_prose_deps' has no attribute 'parse_dependencies'`

- [ ] **Step 3: Implement the parse layer**

```python
# add to scan_prose_deps.py (after fetch layer)

# Blocker-phrase allowlist — deliberately narrow: only verbs that mean "this is
# blocked until #M". Narrative verbs (informs/consumes/relates to/see/fixes/...)
# are NOT listed, so they never match — that IS the exclude mechanism (a narrow
# allowlist), proven by test_parse_narrative_phrases_excluded.
_BLOCKER_PHRASES = [
    r"depends on",
    r"blocked by",
    r"blocked on",
    r"blocked-by:",
    r"gated on",
    r"requires",
]
BLOCKER_RE = re.compile(
    r"(?:" + "|".join(_BLOCKER_PHRASES) + r")\s*#(\d+)",
    re.IGNORECASE,
)


def parse_dependencies(number, body):
    """Return sorted, deduped (dependent, blocker) pairs from one issue body.

    Pure: no I/O. Matches only the blocker-phrase allowlist; drops self-references.
    Blocker open/closed and PR-vs-issue filtering happen later, in reconcile().
    """
    if not body:
        return []
    blockers = {int(m.group(1)) for m in BLOCKER_RE.finditer(body)}
    blockers.discard(number)  # self-reference
    return [(number, b) for b in sorted(blockers)]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_scan_prose_deps.py -k parse -v`
Expected: PASS (all parametrized cases)

- [ ] **Step 5: Commit**

```bash
git add scripts/pm/scan_prose_deps.py workflow/tests/test_scan_prose_deps.py
git commit -m "feat(pm): scan_prose_deps parse layer — allowlist phrase extraction (#722)"
```

---

## Task 3: Reconcile layer (classify + gh-backed lookups)

**Files:**
- Modify: `scripts/pm/scan_prose_deps.py`
- Test: `workflow/tests/test_scan_prose_deps.py`

- [ ] **Step 1: Write failing tests for `classify()` (pure given inputs)**

```python
# append to workflow/tests/test_scan_prose_deps.py

def test_classify_needs_wiring():
    r = spd.classify(745, 722, blocker_meta={"state": "open", "is_pr": False}, existing=set())
    assert r == "needs-wiring"

def test_classify_already_wired():
    r = spd.classify(745, 722, blocker_meta={"state": "open", "is_pr": False}, existing={722})
    assert r == "already-wired"

def test_classify_closed_blocker():
    r = spd.classify(594, 211, blocker_meta={"state": "closed", "is_pr": False}, existing=set())
    assert r == "closed-blocker"

def test_classify_un_wireable_pr():
    r = spd.classify(680, 714, blocker_meta={"state": "open", "is_pr": True}, existing=set())
    assert r == "un-wireable-pr"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_scan_prose_deps.py -k classify -v`
Expected: FAIL — `AttributeError: ... 'classify'`

- [ ] **Step 3: Implement `classify` + the gh-backed lookups**

```python
# add to scan_prose_deps.py

def classify(dependent, blocker, *, blocker_meta, existing):
    """Classify one (dependent, blocker) pair. Pure given resolved inputs."""
    if dependent == blocker:
        return "self-ref"
    if blocker in existing:
        return "already-wired"
    if blocker_meta["state"] == "closed":
        return "closed-blocker"  # convention: only wire currently-open blockers
    if blocker_meta["is_pr"]:
        return "un-wireable-pr"  # native deps are issue<->issue
    return "needs-wiring"


def issue_meta(number):
    """{'state': 'open'|'closed', 'is_pr': bool} via the REST issues endpoint
    (which serves both issues and PRs; a PR carries a 'pull_request' key)."""
    obj = gh("api", f"repos/{REPO}/issues/{number}")
    return {"state": obj["state"], "is_pr": "pull_request" in obj}


def native_blockers(number):
    """Set of issue numbers this issue is already natively blockedBy."""
    q = (
        'query { repository(owner: "Jin-HoMLee", name: "splice-neoepitope-pipeline") {'
        f"  issue(number: {number}) {{ blockedBy(first: 50) {{ nodes {{ number }} }} }}"
        "}}"
    )
    data = gh("api", "graphql", "-f", f"query={q}")
    nodes = data["data"]["repository"]["issue"]["blockedBy"]["nodes"]
    return {n["number"] for n in nodes}
```

- [ ] **Step 4: Write failing test for `reconcile()` with monkeypatched gh lookups**

```python
# append to workflow/tests/test_scan_prose_deps.py

def test_reconcile_classifies_each_pair(monkeypatch):
    meta = {
        722: {"state": "open", "is_pr": False},   # needs-wiring
        211: {"state": "closed", "is_pr": False}, # closed-blocker
        714: {"state": "open", "is_pr": True},    # un-wireable-pr
        708: {"state": "open", "is_pr": False},   # already-wired (in existing)
    }
    monkeypatch.setattr(spd, "issue_meta", lambda n: meta[n])
    monkeypatch.setattr(spd, "native_blockers", lambda n: {708} if n == 709 else set())

    pairs = [(745, 722), (594, 211), (680, 714), (709, 708)]
    records = spd.reconcile(pairs)
    actions = {(r["dependent"], r["blocker"]): r["action"] for r in records}
    assert actions == {
        (745, 722): "needs-wiring",
        (594, 211): "closed-blocker",
        (680, 714): "un-wireable-pr",
        (709, 708): "already-wired",
    }
```

- [ ] **Step 5: Run to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_scan_prose_deps.py -k reconcile -v`
Expected: FAIL — `AttributeError: ... 'reconcile'`

- [ ] **Step 6: Implement `reconcile()`**

```python
# add to scan_prose_deps.py

def reconcile(pairs):
    """Resolve + classify each (dependent, blocker) pair into a record dict.
    Caches per-issue lookups so each blocker/dependent is fetched once."""
    meta_cache = {}
    edges_cache = {}

    def meta(n):
        if n not in meta_cache:
            meta_cache[n] = issue_meta(n)
        return meta_cache[n]

    def edges(n):
        if n not in edges_cache:
            edges_cache[n] = native_blockers(n)
        return edges_cache[n]

    records = []
    for dependent, blocker in pairs:
        bmeta = meta(blocker)
        action = classify(dependent, blocker, blocker_meta=bmeta, existing=edges(dependent))
        records.append({
            "dependent": dependent,
            "blocker": blocker,
            "state": bmeta["state"],
            "action": action,
        })
    return records
```

- [ ] **Step 7: Run all tests to verify pass**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_scan_prose_deps.py -v`
Expected: PASS (all)

- [ ] **Step 8: Commit**

```bash
git add scripts/pm/scan_prose_deps.py workflow/tests/test_scan_prose_deps.py
git commit -m "feat(pm): scan_prose_deps reconcile layer — classify + native-graph diff (#722)"
```

---

## Task 4: Act layer (report + wire) + `main()` CLI

**Files:**
- Modify: `scripts/pm/scan_prose_deps.py`
- Test: `workflow/tests/test_scan_prose_deps.py`

- [ ] **Step 1: Write failing tests for `render_report` and `main` wiring**

```python
# append to workflow/tests/test_scan_prose_deps.py

def _records():
    return [
        {"dependent": 745, "blocker": 722, "state": "open", "action": "needs-wiring"},
        {"dependent": 594, "blocker": 211, "state": "closed", "action": "closed-blocker"},
    ]

def test_render_report_lists_each_record():
    out = spd.render_report(_records())
    assert "745" in out and "722" in out and "needs-wiring" in out
    assert "594" in out and "closed-blocker" in out

def test_render_report_empty():
    assert "no prose-dependency drift" in spd.render_report([]).lower()

def _run_main(monkeypatch, capsys, argv, pairs, records):
    monkeypatch.setattr(spd, "fetch_open_issues",
                        lambda: [{"number": 745, "title": "x", "body": "depends on #722", "state": "OPEN"}])
    monkeypatch.setattr(spd, "parse_dependencies", lambda n, b: pairs)
    monkeypatch.setattr(spd, "reconcile", lambda p: records)
    monkeypatch.setattr(sys, "argv", ["scan_prose_deps.py", *argv])
    rc = spd.main()
    return rc, capsys.readouterr().out

def test_main_report_default_exit_zero(monkeypatch, capsys):
    rc, out = _run_main(monkeypatch, capsys, [], [(745, 722)], _records())
    assert rc == 0 and "needs-wiring" in out

def test_main_check_exits_2_on_drift(monkeypatch, capsys):
    rc, _ = _run_main(monkeypatch, capsys, ["--check"], [(745, 722)], _records())
    assert rc == 2  # a needs-wiring record is present

def test_main_check_exits_0_when_clean(monkeypatch, capsys):
    clean = [{"dependent": 594, "blocker": 211, "state": "closed", "action": "closed-blocker"}]
    rc, _ = _run_main(monkeypatch, capsys, ["--check"], [(594, 211)], clean)
    assert rc == 0  # no needs-wiring

def test_main_apply_wires_only_needs_wiring(monkeypatch, capsys):
    wired = []
    monkeypatch.setattr(spd, "wire", lambda recs: wired.extend((r["dependent"], r["blocker"]) for r in recs) or [])
    rc, _ = _run_main(monkeypatch, capsys, ["--apply"], [(745, 722)], _records())
    # wire() is handed only the needs-wiring subset
    assert rc == 0 and wired == [(745, 722)]
```

- [ ] **Step 2: Run to verify they fail**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_scan_prose_deps.py -k "report or main" -v`
Expected: FAIL — `AttributeError: ... 'render_report'` / `'main'`

- [ ] **Step 3: Implement act layer + `main()`**

```python
# add to scan_prose_deps.py

_ACTION_ORDER = ["needs-wiring", "un-wireable-pr", "already-wired", "closed-blocker", "self-ref"]


def render_report(records):
    if not records:
        return "scan_prose_deps: no prose-dependency drift found.\n"
    lines = [f"{'DEPENDENT':>9}  {'BLOCKER':>7}  {'STATE':<7}  ACTION",
             "-" * 44]
    key = lambda r: (_ACTION_ORDER.index(r["action"]), r["dependent"], r["blocker"])
    for r in sorted(records, key=key):
        lines.append(f"#{r['dependent']:<8} #{r['blocker']:<6} {r['state']:<7}  {r['action']}")
    return "\n".join(lines) + "\n"


def wire(records):
    """POST a native blockedBy edge for each record (must be needs-wiring).
    Uses the REST dependencies endpoint, which takes the blocker's numeric DB id."""
    wired = []
    for r in records:
        dependent, blocker = r["dependent"], r["blocker"]
        blocker_db_id = gh("api", f"repos/{REPO}/issues/{blocker}", parse_json=True)["id"]
        gh("api", "--method", "POST",
           f"repos/{REPO}/issues/{dependent}/dependencies/blocked_by",
           "-F", f"issue_id={blocker_db_id}", parse_json=False)
        print(f"  wired: #{dependent} blocked_by #{blocker}")
        wired.append((dependent, blocker))
    return wired


def _scan(issue_number=None):
    """Fetch -> parse -> reconcile. Single-issue if issue_number given."""
    if issue_number is not None:
        issues = [i for i in fetch_open_issues() if i["number"] == issue_number]
    else:
        issues = fetch_open_issues()
    pairs = []
    for i in issues:
        pairs.extend(parse_dependencies(i["number"], i.get("body")))
    return reconcile(pairs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--issue", type=int, help="restrict scan to a single dependent issue")
    parser.add_argument("--check", action="store_true", help="exit 2 if any needs-wiring drift")
    parser.add_argument("--apply", action="store_true", help="wire the needs-wiring edges")
    parser.add_argument("--only", type=int, nargs="*", default=None,
                        help="with --apply: wire only these dependent issue numbers")
    args = parser.parse_args()

    records = _scan(args.issue)
    needs = [r for r in records if r["action"] == "needs-wiring"]

    if args.apply:
        subset = needs if args.only is None else [r for r in needs if r["dependent"] in args.only]
        if not subset:
            print("scan_prose_deps: nothing to wire.")
            return 0
        wire(subset)
        return 0

    print(render_report(records), end="")
    if args.check:
        return 2 if needs else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run all tests to verify pass**

Run: `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_scan_prose_deps.py -v`
Expected: PASS (all)

- [ ] **Step 5: Make the script executable + commit**

```bash
chmod +x scripts/pm/scan_prose_deps.py
git add scripts/pm/scan_prose_deps.py workflow/tests/test_scan_prose_deps.py
git commit -m "feat(pm): scan_prose_deps act layer (report/wire) + CLI (#722)"
```

---

## Task 5: Live smoke test (read-only)

**Files:** none (operational verification).

- [ ] **Step 1: Run `--report` against the live board**

Run: `scripts/pm/scan_prose_deps.py --report`
Expected: a drift table. Sanity-check that the known #745→#722 pair appears as `needs-wiring`, and a known closed chain (e.g. #594→#211) shows `closed-blocker`.

- [ ] **Step 2: Run single-issue mode on #745**

Run: `scripts/pm/scan_prose_deps.py --issue 745`
Expected: exactly the #745→#722 `needs-wiring` row (this is the call-site the DoR/best-next backstop uses).

- [ ] **Step 3: Confirm `--check` reflects drift**

Run: `scripts/pm/scan_prose_deps.py --check; echo "exit=$?"`
Expected: `exit=2` (needs-wiring drift present at this point).

No commit (read-only).

---

## Task 6: Execute the backfill (AC items 1 + 2)

**Files:** none (board mutations + verification).

- [ ] **Step 1: Capture the review set**

Run: `scripts/pm/scan_prose_deps.py --report | tee /tmp/prose_dep_report.txt`
**PM review gate (non-negotiable):** read every `needs-wiring` row. Strike any false positive (a blocker phrase in negated/historical context) by noting its dependent number for exclusion. PR-blockers and closed-blockers are already excluded by classification.

- [ ] **Step 2: Wire the reviewed set**

If the full `needs-wiring` set is correct:
Run: `scripts/pm/scan_prose_deps.py --apply`
If excluding reviewed-out false positives, wire the rest explicitly:
Run: `scripts/pm/scan_prose_deps.py --apply --only <good dependent numbers…>`
Expected: one `wired: #N blocked_by #M` line per edge.

- [ ] **Step 3: Verify the native graph now answers `is:blocked` (AC item 2)**

Run: `gh search issues "is:open is:blocked" --repo Jin-HoMLee/splice-neoepitope-pipeline --json number,title`
Expected: the real blocked set (≥ the rows just wired), not ~0.

- [ ] **Step 4: Confirm `--check` is now clean**

Run: `scripts/pm/scan_prose_deps.py --check; echo "exit=$?"`
Expected: `exit=0` (no remaining needs-wiring drift).

No commit (board state only — there is no code change in this task).

---

## Task 7: Convention cross-ref (memory) + PR

**Files:**
- Modify (personas memory, MM-committed): `memory/shared/feedback_dependency_tracking.md`

- [ ] **Step 1: Add the scanner cross-ref to the dependency-tracking memory**

Under "How to apply", append a bullet:

```markdown
- **Reconcile prose → native in bulk / on demand:** `scripts/pm/scan_prose_deps.py` finds prose-stated dependencies (the blocker-phrase allowlist) and diffs them against the native graph. `--report` lists drift; `--apply` wires the missing edges (after a human review of the report); `--check` exits non-zero on drift (board-hygiene sweep); `--issue N` scopes to one candidate — the call the DoR commitment-act check (`feedback_board_hygiene.md`) and the best-next backstop (`feedback_best_next_issue.md` Step 1) use. Built for [Issue #722](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/722).
```

Flag for MM in the PR description (memory file → MM commits separately).

- [ ] **Step 2: Push the branch + open the PR**

```bash
git push -u origin design/pm/issue-722-prose-dependency-reconciler
gh pr create --repo Jin-HoMLee/splice-neoepitope-pipeline \
  --title "feat(pm): prose-dependency reconciler — backfill native blockedBy (#722)" \
  --body "$(cat <<'EOF'
Closes #722.

Builds `scripts/pm/scan_prose_deps.py` (fetch → parse → reconcile → act) and backfills the existing prose dependencies into native `blockedBy` edges.

## Test plan
- [ ] `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_scan_prose_deps.py -v` green
- [ ] `--report` lists drift; `--issue 745` returns the #745→#722 row
- [ ] backfill applied; `gh search issues "is:blocked"` returns the real set
- [ ] `--check` exits 0 after backfill

## Memory edits — for MM to commit
- `shared/feedback_dependency_tracking.md` — scanner cross-ref under "How to apply".

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

- [ ] **Step 3: Offer a bot review, write the lab-notebook entry post-review, then merge via the gate**

Per the closure ritual: offer `@claude review` on this non-trivial PR, write the `research/lab_notebook/pm.md` entry referencing the PR after review, then:
Run: `bash scripts/audit_and_merge.sh <PR_NUMBER> --squash --delete-branch`

---

## Self-Review (completed during authoring)

**Spec coverage:** §3 deliverable → Tasks 1-4. §4 four layers → Tasks 1 (fetch), 2 (parse), 3 (reconcile), 4 (act). §5 false-positive control → Task 2 (narrow allowlist + negative tests). §6 four modes → Task 4 (`main`). §7 backfill flow → Task 6. §8 integration → Task 7 cross-ref. §9 testing → Tasks 2-4. §10 doc update → Task 7. §11 AC → Tasks 4 (script+modes+tests), 6 (backfill+verify), 7 (convention doc). All covered.

**Placeholder scan:** none — every code step shows complete code; the `--only <numbers>` in Task 6 is an operational argument the PM fills from the live report, not a code placeholder.

**Type consistency:** `parse_dependencies(number, body)`, `classify(dependent, blocker, *, blocker_meta, existing)`, `reconcile(pairs) -> [{dependent, blocker, state, action}]`, `wire(records)`, `render_report(records)`, `main()` — names + signatures match across Tasks 2-4 and the API block.

**Deviation noted:** spec §5's "narrative-exclude list" is implemented as a *narrow allowlist + negative tests* (Task 2) rather than a second exclude pass — same guarantee, less code (DRY/YAGNI). Documented in the parse-layer docstring.
