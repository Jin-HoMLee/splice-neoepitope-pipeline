# GraphQL Budget Attribution Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give per-consumer visibility into GitHub GraphQL budget spend so the dominant consumer can be identified from data, not guessed ([Issue #1165](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1165), AC1 + AC3 + AC4).

**Architecture:** A single stdlib-only instrumentation helper (`graphql_meter.py`) exposes a rateLimit fragment + a fail-open `log_graphql_spend()` that appends one uniform line per request to `.agents/graphql_spend.jsonl`. Each GraphQL call site injects the fragment (query sites) or a cost-0 probe (mutation sites) and calls the meter, tagged by consumer name. Plus: document the separate REST budget as a fallback (AC4) and harden the poll-loop rule (AC3). AC2 (name top consumer) + AC5 (identities-vs-cap decision) are carved to a data-collection follow-up Issue because they need the meter running live first.

**Tech Stack:** Python 3 stdlib only (json, datetime, pathlib), GitHub GraphQL API via `gh api graphql`, pytest.

## Global Constraints

- **Fail-open everywhere:** telemetry must NEVER raise into the call site it measures. Every meter path swallows its own errors (`except Exception: pass`), matching the existing `_log_fire` fire-log idiom.
- **stdlib only** for `graphql_meter.py` (it is imported by dependency-light hooks that fire on every tool call).
- **No em-dash / en-dash** in any authored content (project rule; the pre-commit guard enforces it).
- **Log line schema** (exact keys, this order): `{"ts","consumer","cost","remaining","query_name"}`. `ts` is `datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")`. `cost` is an int for query sites, `null` for mutation-probe sites (the probe's own cost is 0 and is not the mutation's cost - a documented v1 limitation).
- **GraphQL only.** REST calls (`gh api repos/{o}/{r}/...`) draw a separate 5,000/hr budget and are OUT of scope for this instrumentation; do not meter them.
- **Keep `--jq` out of instrumented `gh` calls:** the meter reads `rateLimit` from the parsed response, so a `--jq`-filtered call must be converted to parse-in-Python first (this also aligns with the `gh_client.py` house rule, Issue #1011).
- Board #9 IDs (stable, user-level project): `PROJECT_ID = PVT_kwHOB17eGc4BSomP`.

---

### Task 1: `graphql_meter.py` module + tests + gitignore

**Files:**
- Create: `.agents/hooks/graphql_meter.py`
- Test: `tools/ci/test_graphql_meter.py`
- Modify: `.gitignore` (add `.agents/graphql_spend.jsonl` beside the existing `.agents/hook_fires.jsonl` entry at line 337)

**Interfaces:**
- Produces:
  - `RATE_LIMIT_FRAGMENT: str` = `"rateLimit { cost remaining }"` (inject as a sibling field at a query's top level)
  - `RATE_LIMIT_PROBE_QUERY: str` = `"query { rateLimit { cost remaining } }"` (cost-0 standalone probe for mutation sites)
  - `extract_rate_limit(response: dict) -> tuple[int | None, int | None]` returns `(cost, remaining)` from `response["data"]["rateLimit"]`, or `(None, None)` if absent/malformed
  - `log_graphql_spend(consumer: str, response: dict, *, query_name: str = "") -> None` (query sites: logs cost+remaining; fail-open)
  - `log_graphql_probe(consumer: str, probe_response: dict, *, query_name: str = "") -> None` (mutation sites: logs remaining, cost=null; fail-open)
  - `SPEND_LOG_PATH: Path` = `.agents/graphql_spend.jsonl`

- [ ] **Step 1: Write the failing tests**

```python
# tools/ci/test_graphql_meter.py
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / ".agents" / "hooks"))
import graphql_meter as gm  # noqa: E402


def test_fragment_and_probe_are_literal_strings():
    assert gm.RATE_LIMIT_FRAGMENT == "rateLimit { cost remaining }"
    assert gm.RATE_LIMIT_PROBE_QUERY == "query { rateLimit { cost remaining } }"


def test_extract_rate_limit_reads_data_ratelimit():
    resp = {"data": {"rateLimit": {"cost": 2, "remaining": 4998}}}
    assert gm.extract_rate_limit(resp) == (2, 4998)


def test_extract_rate_limit_missing_is_none_none():
    assert gm.extract_rate_limit({"data": {}}) == (None, None)
    assert gm.extract_rate_limit({}) == (None, None)
    assert gm.extract_rate_limit(None) == (None, None)
    assert gm.extract_rate_limit("not a dict") == (None, None)


def test_log_graphql_spend_appends_one_line(tmp_path, monkeypatch):
    log = tmp_path / "graphql_spend.jsonl"
    monkeypatch.setattr(gm, "SPEND_LOG_PATH", log)
    resp = {"data": {"rateLimit": {"cost": 26, "remaining": 4700}}}
    gm.log_graphql_spend("board_open_items", resp, query_name="board")
    lines = log.read_text().splitlines()
    assert len(lines) == 1
    rec = json.loads(lines[0])
    assert rec["consumer"] == "board_open_items"
    assert rec["cost"] == 26
    assert rec["remaining"] == 4700
    assert rec["query_name"] == "board"
    assert rec["ts"].endswith("Z")
    assert list(rec.keys()) == ["ts", "consumer", "cost", "remaining", "query_name"]


def test_log_graphql_probe_nulls_cost(tmp_path, monkeypatch):
    log = tmp_path / "graphql_spend.jsonl"
    monkeypatch.setattr(gm, "SPEND_LOG_PATH", log)
    probe = {"data": {"rateLimit": {"cost": 0, "remaining": 4321}}}
    gm.log_graphql_probe("post_gh_pr_create", probe, query_name="status_mutation")
    rec = json.loads(log.read_text().splitlines()[0])
    assert rec["cost"] is None
    assert rec["remaining"] == 4321


def test_log_is_fail_open_on_unwritable_path(monkeypatch):
    # A directory that cannot be created / written must not raise.
    monkeypatch.setattr(gm, "SPEND_LOG_PATH", Path("/proc/nonexistent/graphql_spend.jsonl"))
    gm.log_graphql_spend("x", {"data": {"rateLimit": {"cost": 1, "remaining": 1}}})  # no raise


def test_log_is_fail_open_on_malformed_response(tmp_path, monkeypatch):
    log = tmp_path / "graphql_spend.jsonl"
    monkeypatch.setattr(gm, "SPEND_LOG_PATH", log)
    gm.log_graphql_spend("x", {"garbage": True})  # no rateLimit; must not raise
    rec = json.loads(log.read_text().splitlines()[0])
    assert rec["cost"] is None and rec["remaining"] is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `workflow/tests/.venv/bin/python -m pytest tools/ci/test_graphql_meter.py -v`
Expected: FAIL - `ModuleNotFoundError: No module named 'graphql_meter'`

- [ ] **Step 3: Write the module**

```python
# .agents/hooks/graphql_meter.py
"""Single instrumentation point for GitHub GraphQL budget attribution (Issue #1165).

The GraphQL budget (5,000 pts/hr) is keyed per GitHub USER and shared across every
concurrent Claude session; none of them can see the others spending it, so
exhaustion keeps getting misattributed to the board. This module is the one place
per-consumer spend is recorded: inject RATE_LIMIT_FRAGMENT into a query (or run
RATE_LIMIT_PROBE_QUERY after a mutation), then call log_graphql_spend /
log_graphql_probe with the parsed response.

Design notes:
- stdlib only: imported by dependency-light hooks that fire on every tool call.
- FAIL-OPEN: telemetry must never raise into the call site it measures.
- Attribution, not enforcement: no shared live counter (that would need a daemon/
  Redis, unjustified for a handful of local agents on GitHub's own bucket).
"""
import json
from datetime import datetime, timezone
from pathlib import Path

RATE_LIMIT_FRAGMENT = "rateLimit { cost remaining }"
RATE_LIMIT_PROBE_QUERY = "query { rateLimit { cost remaining } }"

# Sibling of .agents/hook_fires.jsonl (both gitignored). __file__ is
# .agents/hooks/graphql_meter.py, so parents[1] is .agents/.
SPEND_LOG_PATH = Path(__file__).resolve().parents[1] / "graphql_spend.jsonl"


def extract_rate_limit(response):
    """Return (cost, remaining) from response['data']['rateLimit'], else (None, None)."""
    try:
        rl = response["data"]["rateLimit"]
        return rl.get("cost"), rl.get("remaining")
    except (TypeError, KeyError, AttributeError):
        return (None, None)


def _append(consumer, cost, remaining, query_name):
    payload = {
        "ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "consumer": consumer,
        "cost": cost,
        "remaining": remaining,
        "query_name": query_name,
    }
    line = json.dumps(payload, separators=(",", ":")) + "\n"
    try:
        SPEND_LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
        with open(SPEND_LOG_PATH, "a", encoding="utf-8") as f:
            f.write(line)
    except OSError:
        pass


def log_graphql_spend(consumer, response, *, query_name=""):
    """Query sites: log measured cost + remaining. Fail-open."""
    try:
        cost, remaining = extract_rate_limit(response)
        _append(consumer, cost, remaining, query_name)
    except Exception:
        pass


def log_graphql_probe(consumer, probe_response, *, query_name=""):
    """Mutation sites: log remaining only (cost=null). The probe's own cost (0) is
    not the mutation's cost - a documented v1 limitation. Fail-open."""
    try:
        _, remaining = extract_rate_limit(probe_response)
        _append(consumer, None, remaining, query_name)
    except Exception:
        pass
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `workflow/tests/.venv/bin/python -m pytest tools/ci/test_graphql_meter.py -v`
Expected: PASS (7 tests)

- [ ] **Step 5: Add the gitignore entry**

Add this line to `.gitignore` immediately after the `.agents/hook_fires.jsonl` line (currently line 337):

```
.agents/graphql_spend.jsonl
```

Verify: `git check-ignore .agents/graphql_spend.jsonl` prints the path (exit 0).

- [ ] **Step 6: Commit**

```bash
git add .agents/hooks/graphql_meter.py tools/ci/test_graphql_meter.py .gitignore
git commit -m "feat(infra): graphql_meter single instrumentation point for budget attribution (#1165)"
```

---

### Task 2: Instrument `board_open_items.py` (query site)

**Files:**
- Modify: `scripts/board_open_items.py` (QUERY block at lines 79-118; runner `fetch_all_items` around lines 149-176)

**Interfaces:**
- Consumes: `graphql_meter.RATE_LIMIT_FRAGMENT`, `graphql_meter.log_graphql_spend`

- [ ] **Step 1: Write the failing test**

```python
# append to tools/ci/test_graphql_meter.py
def test_board_query_carries_ratelimit_fragment():
    import importlib.util
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location("boi", root / "scripts" / "board_open_items.py")
    boi = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(boi)
    assert "rateLimit { cost remaining }" in boi.QUERY
```

- [ ] **Step 2: Run to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest tools/ci/test_graphql_meter.py::test_board_query_carries_ratelimit_fragment -v`
Expected: FAIL (`assert ... in QUERY` is False)

- [ ] **Step 3: Inject the fragment and meter each page**

In `scripts/board_open_items.py`, add the import near the other imports (after line 43 `import subprocess`):

```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / ".agents" / "hooks"))
import graphql_meter  # noqa: E402
```

Inject the fragment into `QUERY` as a sibling of the `user(...)` selection - change the opening of the query (line 80-81) from:

```
query($owner: String!, $number: Int!, $after: String) {
  user(login: $owner) {
```

to:

```
query($owner: String!, $number: Int!, $after: String) {
  rateLimit { cost remaining }
  user(login: $owner) {
```

In `fetch_all_items`, right after `data = json.loads(r.stdout)` (line 172) and before the `errors` check, add:

```python
        graphql_meter.log_graphql_spend("board_open_items", data, query_name="board_page")
```

- [ ] **Step 4: Run tests + a smoke import**

Run: `workflow/tests/.venv/bin/python -m pytest tools/ci/test_graphql_meter.py -v`
Expected: PASS. Also confirm no syntax/import regression: `workflow/tests/.venv/bin/python -c "import ast; ast.parse(open('scripts/board_open_items.py').read())"`

- [ ] **Step 5: Commit**

```bash
git add scripts/board_open_items.py tools/ci/test_graphql_meter.py
git commit -m "feat(infra): meter board_open_items GraphQL spend (#1165)"
```

---

### Task 3: Instrument `check_gh_issue_develop_parent.py` (query site)

**Files:**
- Modify: `.agents/hooks/check_gh_issue_develop_parent.py` (query at lines 195-210)

**Interfaces:**
- Consumes: `graphql_meter.RATE_LIMIT_FRAGMENT`, `graphql_meter.log_graphql_spend`

- [ ] **Step 1: Add import + inject fragment + meter**

The hook already `sys.path.insert(0, ... parent)` for its own helpers, so import the sibling module the same way (add near the existing imports):

```python
import graphql_meter  # noqa: E402   (same dir; path already inserted for _shell_parse siblings)
```

If no `sys.path.insert(0, str(Path(__file__).resolve().parent))` exists yet in this file, add it before the import.

Inject the fragment into the query string (lines 196-199), changing:

```python
    query = (
        "query($o:String!,$n:String!,$num:Int!){"
        "repository(owner:$o,name:$n){issue(number:$num){"
        "subIssuesSummary{total}}}}"
    )
```

to:

```python
    query = (
        "query($o:String!,$n:String!,$num:Int!){"
        "rateLimit{cost remaining}"
        "repository(owner:$o,name:$n){issue(number:$num){"
        "subIssuesSummary{total}}}}"
    )
```

After `data = json.loads(res.stdout)` (line 203) and before the return, add:

```python
        graphql_meter.log_graphql_spend("check_gh_issue_develop_parent", data, query_name="parent_check")
```

- [ ] **Step 2: Verify the existing hook test suite still passes**

Run: `workflow/tests/.venv/bin/python -m pytest tools/ci/ -k "develop_parent or hook_liveness" -v`
Expected: PASS (behavior unchanged; the added field + meter call are side-effect-free to the parent-total logic).

- [ ] **Step 3: Commit**

```bash
git add .agents/hooks/check_gh_issue_develop_parent.py
git commit -m "feat(infra): meter check_gh_issue_develop_parent GraphQL spend (#1165)"
```

---

### Task 4: Instrument the two mutation hooks via cost-0 probe

**Files:**
- Modify: `.agents/hooks/post_gh_pr_create.py` (`_set_status`, lines 273-282; reuses local `_gh`)
- Modify: `.agents/hooks/post_gh_pr_review_request.py` (status mutation near line 270; reuses local `_gh`)

**Interfaces:**
- Consumes: `graphql_meter.RATE_LIMIT_PROBE_QUERY`, `graphql_meter.log_graphql_probe`

- [ ] **Step 1: Add a probe-and-log helper call after each mutation**

Both hooks have a local `def _gh(*args) -> subprocess.CompletedProcess`. Add the meter import (both files already `sys.path.insert(0, ... parent)` for `_shell_parse` / siblings):

```python
import graphql_meter  # noqa: E402
```

In `post_gh_pr_create.py` `_set_status`, after the mutation `_gh(...)` call returns (end of the function body, line 282), add:

```python
    try:
        probe = json.loads(_gh("api", "graphql", "-f", f"query={graphql_meter.RATE_LIMIT_PROBE_QUERY}").stdout)
        graphql_meter.log_graphql_probe("post_gh_pr_create", probe, query_name="status_mutation")
    except Exception:
        pass
```

In `post_gh_pr_review_request.py`, immediately after the status-mutation `_gh("api", "graphql", ...)` call (near line 270), add the same block with `consumer="post_gh_pr_review_request"`.

- [ ] **Step 2: Verify the hooks' existing tests still pass**

Run: `workflow/tests/.venv/bin/python -m pytest tools/ci/ -k "pr_create or pr_review_request or hook_liveness" -v`
Expected: PASS (the probe is wrapped in its own `try/except`, so a probe failure cannot change mutation behavior).

- [ ] **Step 3: Commit**

```bash
git add .agents/hooks/post_gh_pr_create.py .agents/hooks/post_gh_pr_review_request.py
git commit -m "feat(infra): meter mutation-hook GraphQL spend via cost-0 probe (#1165)"
```

---

### Task 5: Instrument `recheck_dispatch.py` GraphQL calls (own task, higher blast radius)

**Files:**
- Modify: `.agents/hooks/recheck_dispatch.py` - GraphQL sites only: `lookup_issue_for_item` (lines 237-246, uses `--jq`), `get_issue_target_date` (lines 268-290, parses in Python), `_mutate_target_date` (lines 342-380, mutation)

**Interfaces:**
- Consumes: `graphql_meter.RATE_LIMIT_FRAGMENT`, `graphql_meter.RATE_LIMIT_PROBE_QUERY`, `graphql_meter.log_graphql_spend`, `graphql_meter.log_graphql_probe`

**Scope note:** This hook's `get_issue_milestone`, `prior_milestones_for_issue`, and `lookup_milestone_number_by_title` calls that hit `gh api repos/...` are REST (separate budget) - do NOT meter them. Only the three `gh api graphql` sites below.

- [ ] **Step 1: Add import**

Add `import graphql_meter` after the existing `sys.path.insert(0, str(Path(__file__).resolve().parent))` (the hook already inserts its own dir).

- [ ] **Step 2: `lookup_issue_for_item` - drop `--jq`, parse in Python, inject fragment + meter**

Replace lines 238-246:

```python
    query = f'query {{ node(id: "{item_id}") {{ ... on ProjectV2Item {{ content {{ ... on Issue {{ number }} }} }} }} }}'
    result = subprocess.run(
        ["gh", "api", "graphql", "-f", f"query={query}",
         "--jq", ".data.node.content.number"],
        capture_output=True, text=True, check=False,
    )
    out = result.stdout.strip()
    return int(out) if out.isdigit() else None
```

with:

```python
    query = f'query {{ rateLimit {{ cost remaining }} node(id: "{item_id}") {{ ... on ProjectV2Item {{ content {{ ... on Issue {{ number }} }} }} }} }}'
    result = subprocess.run(
        ["gh", "api", "graphql", "-f", f"query={query}"],
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        return None
    try:
        data = json.loads(result.stdout)
    except json.JSONDecodeError:
        return None
    graphql_meter.log_graphql_spend("recheck_dispatch", data, query_name="lookup_issue_for_item")
    num = data.get("data", {}).get("node", {}).get("content", {}).get("number")
    return int(num) if isinstance(num, int) else None
```

- [ ] **Step 3: `get_issue_target_date` - inject fragment + meter**

In the query string (lines 274-280) add `rateLimit { cost remaining } ` right after the opening `query { `. After `data = json.loads(result.stdout)` (line 288), add:

```python
    graphql_meter.log_graphql_spend("recheck_dispatch", data, query_name="get_issue_target_date")
```

- [ ] **Step 4: `_mutate_target_date` - probe after mutation**

After `data = json.loads(result.stdout)` (line 375), add:

```python
    try:
        probe = subprocess.run(
            ["gh", "api", "graphql", "-f", f"query={graphql_meter.RATE_LIMIT_PROBE_QUERY}"],
            capture_output=True, text=True, check=False,
        )
        graphql_meter.log_graphql_probe("recheck_dispatch", json.loads(probe.stdout), query_name="mutate_target_date")
    except Exception:
        pass
```

- [ ] **Step 5: Verify existing recheck_dispatch tests pass (the `--jq` removal is behavior-preserving)**

Run: `workflow/tests/.venv/bin/python -m pytest tools/ci/test_recheck_dispatch.py -v`
Expected: PASS. If any test asserted on the `--jq` arg string, update it to assert on the parsed-in-Python behavior (same returned value).

- [ ] **Step 6: Commit**

```bash
git add .agents/hooks/recheck_dispatch.py tools/ci/test_recheck_dispatch.py
git commit -m "feat(infra): meter recheck_dispatch GraphQL calls, drop --jq (#1165)"
```

---

### Task 6: REST-budget fallback - document + route (AC4)

**Files:**
- Modify: `.agents/memory/shared/feedback_board_queries.md` (add a REST-fallback section)
- Modify: `scripts/pm/scan_addressed_comments.py` (route its single-Issue reads that do NOT need board fields to REST)

- [ ] **Step 1: Document the fallback**

Append a section to `.agents/memory/shared/feedback_board_queries.md`:

```markdown
## REST carries a separate budget - use it as the GraphQL-exhaustion fallback

The GraphQL budget (5,000 pts/hr) is keyed per GitHub USER and shared across every
concurrent session (Issue #1165). REST carries a SEPARATE 5,000/hr budget - it sat
at 4,875/5,000 while GraphQL was at 0/5,000, and is what unblocked two stalls.

- `gh issue view N` / `gh pr view N` are GraphQL.
- `gh api repos/{owner}/{repo}/issues/{N}` (read) and `gh api ... -F body=@file` (edit) are REST.
- Board FIELD reads (Status / Priority / Size columns) genuinely need GraphQL; single-Issue/PR
  reads and edits do NOT - route them to REST when GraphQL is tight or exhausted.

Per-consumer GraphQL spend is logged to `.agents/graphql_spend.jsonl` (Issue #1165);
aggregate it to find the dominant consumer rather than assuming it is the board.
```

- [ ] **Step 2: Route the clear single-Issue reads in `scan_addressed_comments.py`**

Locate single-Issue reads that fetch only comment bodies (no board fields). Where the script calls `gh issue view N --json comments`, switch to REST `gh api repos/{owner}/{repo}/issues/{N}/comments` (which returns the same comment data on the separate budget). If a call also needs board Status/Priority/Size, leave it on GraphQL and add a one-line comment noting why.

Verify behavior is unchanged: `workflow/tests/.venv/bin/python -m pytest tools/ci/ scripts/tests/ -k "addressed or scan" -v` (run whichever suite covers this script; if none, do a live `--role developer` run and confirm the same pings surface as before).

- [ ] **Step 3: Commit**

```bash
git add .agents/memory/shared/feedback_board_queries.md scripts/pm/scan_addressed_comments.py
git commit -m "docs(infra): document REST-budget fallback and route single-issue reads (#1165)"
```

---

### Task 7: Harden the poll-loop rule (AC3)

**Files:**
- Modify: `.agents/memory/feedback_bot_review_poll_by_content.md`

**Audit result to record (verified during planning):** `poll_bot_review.sh` is already bounded (25-min default timeout, max-5-consecutive-failures, configurable interval); `run_cloud_gpu.sh`'s `while true` polls gcloud, not GitHub. No runaway exists in our tree. The runaway that triggered #1165 lived in an unrelated repo, outside our control.

- [ ] **Step 1: Add the bounded-poll mandate**

Append to `.agents/memory/feedback_bot_review_poll_by_content.md`:

```markdown
## Poll loops MUST be bounded, and prefer the harness's native re-invoke (Issue #1165)

A `while true` + `gh` poll on a SHARED GitHub budget is a standing drain every other
concurrent session pays for (the GraphQL budget is per-USER, #1165). So:

- **Bound every poll:** an explicit timeout AND an interval floor. `poll_bot_review.sh`
  is the reference shape (default 25-min timeout, max-5-consecutive-failures).
- **Prefer the harness's native task-completion re-invoke** over a hand-rolled `while true`:
  a background task re-invokes you on completion, so polling GitHub in a loop is usually
  redundant. Reach for a bounded poll only when there is genuinely no completion signal.
- **Residual (honest):** a poll loop in an UNRELATED repo is outside our control; this rule
  binds our own tree only.
```

- [ ] **Step 2: Commit**

```bash
git add .agents/memory/feedback_bot_review_poll_by_content.md
git commit -m "docs(infra): mandate bounded poll loops on the shared GitHub budget (#1165)"
```

---

### Task 8: File the AC2 + AC5 data-collection follow-up Issue

**Files:** none (GitHub Issue creation).

- [ ] **Step 1: Create the follow-up Issue**

```bash
gh issue create --repo Jin-HoMLee/splice-neoepitope-pipeline \
  --title "fix(infra): read GraphQL spend log - name top consumer + decide identities-vs-cap (#1165 AC2/AC5)" \
  --label "role:developer" --label "arc:board-governance" \
  --body "$(cat <<'EOF'
## Problem
#1165 shipped the GraphQL spend instrumentation (AC1/AC3/AC4). Its AC2 (name the top
consumer from measured data) and AC5 (decide per-role identities vs a concurrency cap)
need the meter to have run live across a real multi-session window first - the data did
not exist at #1165 merge time.

## Acceptance criteria
- [ ] After a real multi-session window, aggregate `.agents/graphql_spend.jsonl` and name
      the top GraphQL consumer with its measured share of the 5,000/hr budget (AC2).
- [ ] Record the decision on whether per-role bot identities (separate buckets) are worth
      the cost, or whether we simply cap concurrency (AC5). This is arc:board-governance and
      team-wide, so get PM concurrence rather than deciding unilaterally.
- [ ] If a dominant consumer is a hook, propose the concrete reduction (cache, REST reroute,
      debounce) with the measured before/after.

**Priority rationale:** P2 - depends on accumulated telemetry; not blocking once #1165 ships.

**Created by:** Developer
EOF
)"
```

- [ ] **Step 2: Confirm it was created and boarded**

Run: `gh issue list --repo Jin-HoMLee/splice-neoepitope-pipeline --label role:developer --search "GraphQL spend log" --limit 3`
Expected: the new Issue appears. Note its number for the #1165 PR body (comment-defer AC2/AC5 to it).

---

## Self-Review

**Spec coverage:**
- AC1 attribution -> Tasks 1-5 (meter + board + 4 hooks). Covered.
- AC2 name top consumer -> Task 8 follow-up. Covered (deferred by design).
- AC3 poll-loop class -> Task 7 (rule) + audit recorded. Covered.
- AC4 REST fallback -> Task 6 (doc + routing). Covered.
- AC5 identities-vs-cap decision -> Task 8 follow-up. Covered (deferred by design).
- gitignore entry for the new log -> Task 1 Step 5. Covered.
- Testing (meter unit tests, hook-liveness, fail-open) -> Task 1 tests + per-task suite runs. Covered.

**Type consistency:** `RATE_LIMIT_FRAGMENT`, `RATE_LIMIT_PROBE_QUERY`, `extract_rate_limit`, `log_graphql_spend`, `log_graphql_probe`, `SPEND_LOG_PATH` are named identically in Task 1 (definition) and Tasks 2-5 (consumption). Consistent.

**Placeholder scan:** every code step shows complete code; commands have expected output. No TBD/TODO.

## Merge-readiness (after all tasks)

1. Run the full local pytest sweep across the 5 CI dirs (`workflow/tests`, `tools/ci`, `tools/news`, `tools/project_map`, `scripts/tests`) - a meter/hook change can span dirs.
2. Push the branch; open the PR (the PR-open hook auto-requests the first bot review).
3. In the PR body: tick the Test plan, comment-defer AC2/AC5 to the Task 8 follow-up Issue, and reference #1165.
4. Address bot review; write the Developer lab-notebook entry (post-review, pre-merge).
5. STOP at the merge gate - get Jin-Ho's explicit merge go before `scripts/audit_and_merge.sh <PR>`.
