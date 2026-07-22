# GraphQL budget attribution ([Issue #1165](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1165))

## Problem

The GitHub **GraphQL** budget is 5,000 points/hr and is keyed **per GitHub user**, not per token, per session, or per repo.
We run many concurrent Claude sessions authenticated as one user; they all draw from one bucket that none of them can see the others spending.
The budget has hard-stopped two working sessions in a single day, and each time the cause was misattributed to our own board queries - because the board is the only consumer a session can *see*.

The measured reality (2026-07-14, from the Issue): a full 13-page board read costs **26 points**; a single-issue lookup costs **1**; an idle session drains **~600/hr**; concurrent sessions can drop **~1,100 points** between two of one session's own commands.
**Our board reads are cheap. Our concurrency is not.** The fix is *per-consumer visibility*, not a smaller board.

Recorded dead ends (do not re-walk): token rotation shares the same user bucket; a GitHub App buys ~nothing on a personal account; trimming query fields saves nothing (cost is per-request, not per-node).

## Scope of this PR

The Issue has five acceptance criteria. Two of them - naming the top consumer from measured data (AC2) and deciding per-role identities vs a concurrency cap (AC5) - **cannot close until the AC1 instrumentation has run live across several concurrent sessions**. That data does not exist yet.

This PR therefore delivers the shippable, no-wait subset:

- **AC1** - attribution instrumentation (the single instrumentation point + call-site adoption).
- **AC3** - the runaway-poll-loop class (rule hardening + audit of our own loops).
- **AC4** - the REST-budget fallback (documentation + routing the clearest calls).

**AC2 + AC5 are carved into a data-collection follow-up Issue** that reads the accumulated log after a real multi-session window.
This split is deliberate and honest about the temporal dependency; it is not descoping.

## Design

### 1. `graphql_meter.py` - the single instrumentation point (AC1)

Best practice for a shared-budget resource is a **single instrumentation point**, not per-site ad-hoc logging, so metrics are uniform instead of scattered (web cross-check: [Gravitee](https://www.gravitee.io/blog/rate-limiting-apis-scale-patterns-strategies), [SigNoz/OTel](https://signoz.io/blog/monitoring-graphql/)).
The heavyweight forms in the literature (Redis-backed shared counters, OpenTelemetry pipelines) are for **live enforcement/throttling** across instances - which AC1 does **not** ask for and our scale (4 local agents drawing on GitHub's own bucket) does not justify.
So the single point is a **thin, stdlib-only helper**, scoped to attribution, dependency-free so it works in both of our runtime contexts.

**Why a helper and not routing everything through `scripts/pm/gh_client.py`:** our GraphQL call sites live in two runtime contexts with different dependency budgets.
The hooks fire on every tool call and deliberately do **not** import the pm tooling; forcing them through `gh_client` would couple hot-path code to pm helpers.
A dependency-free helper is the single instrumentation point *both* contexts can adopt.

**Location:** `.agents/hooks/graphql_meter.py`, co-located with the existing shared hook helper `_shell_parse.py` (the established precedent for a shared, stdlib-only hook module), importable by pm/board scripts via the same `sys.path.insert` idiom they already use.

**API:**

```python
RATE_LIMIT_FRAGMENT = "rateLimit { cost remaining }"   # inject into a query document

def log_graphql_spend(consumer: str, response: dict, *, query_name: str = "") -> None:
    """Extract rateLimit from a GraphQL response and append one telemetry line.

    Appends {ts, consumer, cost, remaining, query_name} to .agents/graphql_spend.jsonl.
    FAIL-OPEN: never raises. Telemetry must never break the board read or hook it measures.
    """
```

**Log destination:** a dedicated `.agents/graphql_spend.jsonl` telemetry stream, sibling to the existing `.agents/hook_fires.jsonl` fire-log infra. It needs a **new `.gitignore` entry** (`.agents/hook_fires.jsonl` is ignored at `.gitignore:337`; `.agents/graphql_spend.jsonl` is not yet - add it alongside).
One line per request, uniform shape, so the AC2 follow-up aggregates one clean file.

This is a **deliberate, user-approved deviation** from the AC1 wording ("log rateLimit alongside each hook fire in `.agents/hook_fires.jsonl`").
Rationale: `hook_fires.jsonl` is a heterogeneous hook-event log; a dedicated stream is cleaner to aggregate and keeps pm-read telemetry out of a file named for hook fires.
Hooks keep their existing `hook_fires.jsonl` line unchanged.

**Fail-open discipline** (matches the whole hook design ethos): a malformed response, an absent `rateLimit` key, or an unwritable log path is swallowed - the instrumented call proceeds as if unmetered.

### 2. Instrumented call sites (AC1)

Inject `RATE_LIMIT_FRAGMENT` and call `log_graphql_spend()` at each GraphQL consumer, tagged by consumer name:

| Consumer | File | Kind |
|---|---|---|
| `board_open_items` | `scripts/board_open_items.py` | query (per page) |
| `recheck_dispatch` | `.agents/hooks/recheck_dispatch.py` | query |
| `check_gh_issue_develop_parent` | `.agents/hooks/check_gh_issue_develop_parent.py` | query |
| `post_gh_pr_create` | `.agents/hooks/post_gh_pr_create.py` | mutation |
| `post_gh_pr_review_request` | `.agents/hooks/post_gh_pr_review_request.py` | mutation |

**Query sites** select `rateLimit` inline and log `cost` + `remaining` with **zero extra calls**.

**Mutation-only sites:** `rateLimit` is a field on the GraphQL `Query` root, not `Mutation`, so it cannot be selected in a mutation response.
These sites log `remaining` via a single **cost-0** `rateLimit` probe (`query { rateLimit { cost remaining } }` does not deduct from the budget) after the mutation.
Exact per-mutation `cost` attribution is a **documented v1 limitation** - the log line carries `remaining` and a null `cost` for these consumers.

### 3. REST-budget fallback (AC4)

REST carries a **separate** 5,000/hr budget (observed at 4,875/5,000 while GraphQL sat at 0/5,000 - it is what unblocked both stalls).
`gh issue view` is GraphQL; `gh api repos/{owner}/{repo}/issues/{n}` is REST.

- **Document** the fallback in `.agents/memory/shared/feedback_board_queries.md`: REST has its own budget; single-Issue/PR reads + edits should prefer REST; only board-**field** reads (Status/Priority/Size columns) genuinely need GraphQL.
- **Route** the clearest single-issue GraphQL reads that do **not** need board fields to REST (starting with `scripts/pm/scan_addressed_comments.py`'s single-issue lookups), and list the remaining candidates for a follow-up sweep rather than converting everything in this PR.

### 4. Runaway-poll-loop class (AC3)

**Audit result:** our committed poll loops are already safe.
`.agents/skills/awaiting-bot-review/scripts/poll_bot_review.sh` is bounded (default 25-min timeout, max-5-consecutive-failures, configurable interval).
`scripts/run_cloud_gpu.sh`'s `while true` polls **gcloud**, not GitHub, so it is not a GraphQL/REST-budget drain.
There is **no runaway in our tree**.

**Rule hardening:** update the poll-loop guidance (`.agents/memory/feedback_bot_review_poll_by_content.md`) to **mandate bounded polls** - an explicit timeout and an interval floor - and to **prefer the harness's native task-completion re-invoke** (a background task re-invokes on completion) over a hand-rolled `while true` GitHub poll.

**Honest residual:** the runaway that triggered this investigation lived in an **unrelated repo**, outside our control.
A memory rule + our-own-tree audit cannot reach it; this is noted, not falsely claimed as fixed.

### 5. AC2 + AC5 follow-up Issue

File a `role:developer` follow-up Issue (titled around "read the GraphQL spend log; name the top consumer; decide identities-vs-cap") that, after the meter has run across a real multi-session window:

- aggregates `.agents/graphql_spend.jsonl` and **names the top consumer with its measured share** of the hourly budget (AC2);
- **records the decision** on whether per-role identities (separate buckets) are worth the cost, or whether we cap concurrency (AC5).

AC5 is `arc:board-governance` and affects the whole team, so the follow-up flags it for **PM concurrence** rather than a unilateral Developer call.

## Testing

- `tools/ci/test_graphql_meter.py` (pytest): fragment injection; `rateLimit` parsing (query and mutation-probe shapes); **fail-open** on malformed / absent-`rateLimit` / unwritable-path; log-line field shape.
- The meter changes hook query bodies, so verify/extend the hook-liveness contract (`tools/ci/test_hook_liveness_contract.py`, per the new-hook checklist) so a born-dead meter is caught.
- Run the affected hooks' existing test suites to confirm the injected fragment + meter call do not change their control-flow behavior.

## Acceptance criteria mapping

| AC | This PR | Follow-up |
|---|---|---|
| AC1 attribution instrumentation | Yes (§1, §2) | - |
| AC2 name top consumer from data | - | Yes (§5) |
| AC3 poll-loop class | Yes (§4) | - |
| AC4 REST fallback | Yes (§3) | - |
| AC5 identities-vs-cap decision | - | Yes (§5) |

## Non-goals

- No live shared counter / distributed throttler (Redis, daemon) - AC1 is attribution, not enforcement; unjustified at our scale.
- No OpenTelemetry pipeline - overkill for 4 local agents logging to jsonl.
- No wholesale conversion of every single-issue GraphQL read to REST - route the clear wins, sweep the rest later.
- No attempt to reach poll loops in unrelated repos - outside our control.
