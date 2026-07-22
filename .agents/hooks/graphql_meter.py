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
