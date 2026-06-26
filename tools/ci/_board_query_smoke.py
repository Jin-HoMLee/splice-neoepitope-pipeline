"""Helpers for the board GraphQL schema-drift smoke (Issue #771).

Not a ``test_*`` module, so pytest does not collect it.

This is the GitHub-API-boundary analogue of the Snakemake integration rule
(``feedback_integration_run_for_new_rules.md``): the board-query tooling
(``scripts/board_open_items.py``, ``scripts/pm/recheck_*.py``) builds GraphQL
queries that the unit tests validate only against hand-authored fixtures — and a
fixture *is* the assumed response, so it can never disagree with itself. A field
typo or upstream GitHub schema drift passes the hermetic suites green and only
surfaces at runtime.

``check_board_query_shape`` is a pure, structure-only checker (unit-tested
hermetically). ``LIVE_QUERIES`` registers the *real* queries — imported from
their source modules so a field rename flows in automatically — and the live
smoke (nightly, advisory) runs each against the real API and feeds the response
to the checker. It asserts schema validity only; data correctness (which items
are where, counts) stays the manual live-smoke convention
(``feedback_live_integration_smoke.md``).
"""
import importlib.util
import json
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

# tools/ci/ -> repo root
_REPO_ROOT = Path(__file__).resolve().parents[2]


def _load_module_from_path(name: str, relpath: str):
    """Import a non-packaged script (``scripts/`` has no ``__init__.py``) by path."""
    spec = importlib.util.spec_from_file_location(name, _REPO_ROOT / relpath)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# Single source of truth: the live smoke runs the SAME query string the board
# tooling ships, so a rename in board_open_items.py is reflected here with no
# hand-copied duplicate to drift.
_board_open_items = _load_module_from_path("board_open_items", "scripts/board_open_items.py")


@dataclass(frozen=True)
class LiveQuery:
    """A real GraphQL query to smoke-test against the live API."""

    name: str
    query: str
    variables: dict[str, Any] = field(default_factory=dict)
    # Path of keys into the parsed response down to the items connection, so the
    # checker can locate `nodes`/`pageInfo` for queries with a different envelope.
    items_path: tuple[str, ...] = ("data", "user", "projectV2", "items")
    # When True, the response must contain ≥1 Issue node — otherwise the
    # subIssuesSummary check never runs (a vacuous pass). The board always holds
    # issues, so this is a coverage guard, not a data assertion.
    require_issue_node: bool = True


LIVE_QUERIES: list[LiveQuery] = [
    LiveQuery(
        name="board_open_items",
        query=_board_open_items.QUERY,
        variables={
            "owner": _board_open_items.OWNER,
            "number": _board_open_items.PROJECT_NUMBER,
        },
    ),
]


def _dig(obj: Any, path: tuple[str, ...]) -> Any:
    """Walk a key path, returning None at the first missing/non-dict step."""
    cur = obj
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return None
        cur = cur[key]
    return cur


def check_board_query_shape(
    parsed: Any,
    *,
    require_issue_node: bool = True,
) -> list[str]:
    """Return a list of structural problems with a board GraphQL response.

    Empty list == schema-valid. Structure-only: asserts the response is
    schema-valid (no GraphQL ``errors``), parses to the expected shape, the
    pagination contract (``pageInfo``) survives, and ``subIssuesSummary { total }``
    is present and integer-typed on Issue nodes. Deliberately makes NO assertion
    on specific board data values (counts, titles, which items are present) —
    board state changes hourly.
    """
    problems: list[str] = []

    if not isinstance(parsed, dict):
        return [f"response is not a JSON object (got {type(parsed).__name__})"]

    # A field typo / removed field surfaces here — the core schema-drift signal.
    if parsed.get("errors"):
        problems.append(f"GraphQL errors present: {parsed['errors']}")

    if parsed.get("data") is None:
        problems.append("missing `data` (response has no top-level data object)")
        return problems

    items = _dig(parsed, ("data", "user", "projectV2", "items"))
    if items is None:
        problems.append(
            "missing items connection at data.user.projectV2.items "
            "(schema drift or wrong query envelope)"
        )
        return problems

    page_info = items.get("pageInfo")
    if not isinstance(page_info, dict):
        problems.append("missing or malformed pageInfo (board pagination contract)")
    else:
        for k in ("hasNextPage", "endCursor"):
            if k not in page_info:
                problems.append(f"pageInfo missing `{k}`")

    nodes = items.get("nodes")
    if not isinstance(nodes, list):
        problems.append(f"items.nodes is not a list (got {type(nodes).__name__})")
        return problems

    saw_issue = False
    for i, node in enumerate(nodes):
        content = node.get("content") if isinstance(node, dict) else None
        if not isinstance(content, dict):
            continue  # a null content node is legal (deleted/inaccessible item)
        if content.get("__typename") != "Issue":
            continue
        saw_issue = True
        if "subIssuesSummary" not in content:
            problems.append(f"node[{i}] (Issue) missing subIssuesSummary")
            continue
        sub = content["subIssuesSummary"]
        if not isinstance(sub, dict) or "total" not in sub:
            problems.append(f"node[{i}] subIssuesSummary missing `total`")
        elif not isinstance(sub["total"], int):
            problems.append(
                f"node[{i}] subIssuesSummary.total is not an int "
                f"(got {type(sub['total']).__name__})"
            )

    if require_issue_node and not saw_issue:
        problems.append(
            "no Issue nodes in response — subIssuesSummary went unverified "
            "(coverage guard; the board always holds issues)"
        )

    return problems


def run_graphql_with_retry(
    q: LiveQuery,
    *,
    max_attempts: int = 4,
    base_delay: float = 1.0,
    _runner: Callable[[list[str]], "subprocess.CompletedProcess"] | None = None,
    _sleep: Callable[[float], None] = time.sleep,
) -> Any:
    """Run a LiveQuery via ``gh api graphql``, retrying transient failures.

    Transient non-zero exits (TLS timeout, rate-limit, 5xx) are retried with
    capped exponential backoff so a GitHub blip doesn't red the nightly job.
    Returns the parsed JSON. ``_runner``/``_sleep`` are injectable for tests.
    """
    cmd = ["gh", "api", "graphql", "-f", f"query={q.query}"]
    for key, value in q.variables.items():
        cmd.extend(["-F", f"{key}={value}"])

    runner = _runner or (
        lambda c: subprocess.run(c, capture_output=True, text=True, timeout=60)
    )

    last_stderr = ""
    for attempt in range(max_attempts):
        try:
            result = runner(cmd)
        except subprocess.TimeoutExpired:
            last_stderr = "gh call timed out"
            result = None
        if result is not None and result.returncode == 0:
            return json.loads(result.stdout)
        if result is not None:
            last_stderr = result.stderr
        if attempt < max_attempts - 1:
            _sleep(base_delay * (2 ** attempt))

    raise RuntimeError(
        f"gh api graphql failed after {max_attempts} attempts for "
        f"query '{q.name}': {last_stderr}"
    )
