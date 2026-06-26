"""Tests for the board GraphQL schema-drift smoke (Issue #771).

Two layers, mirroring the recheck live-smoke split:

* **Hermetic unit tests** (default, ``-m "not live"``) exercise the pure
  ``check_board_query_shape`` structure-checker against curated well-formed and
  malformed responses. These run per-PR in ``ci-tools-pytest`` and are where the
  checker logic is proven.
* **One live test** (``@pytest.mark.live``, nightly via ``recheck-live-smoke.yml``)
  runs the *real* board query against the live GitHub Projects API and asserts
  the checker finds no structural problems. This is the only layer that can catch
  upstream schema drift or a field typo in the live query — a hand-authored
  fixture can never disagree with itself.
"""
import _board_query_smoke as bqs
import pytest
from _live_gh import REQUIRES_LIVE_GH


def _issue_node(with_sub_summary=True):
    node = {
        "content": {
            "__typename": "Issue",
            "number": 1,
            "title": "x",
            "state": "OPEN",
            "url": "u",
            "labels": {"nodes": []},
        },
        "fieldValues": {"nodes": []},
    }
    if with_sub_summary:
        node["content"]["subIssuesSummary"] = {"total": 0}
    return node


def _well_formed(nodes=None):
    if nodes is None:
        nodes = [_issue_node()]
    return {
        "data": {
            "user": {
                "projectV2": {
                    "items": {
                        "pageInfo": {"hasNextPage": False, "endCursor": None},
                        "nodes": nodes,
                    }
                }
            }
        }
    }


class TestCheckBoardQueryShape:
    def test_well_formed_response_has_no_problems(self):
        assert bqs.check_board_query_shape(_well_formed()) == []

    def test_graphql_errors_flagged(self):
        # A field typo / removed field surfaces as a top-level `errors` array —
        # the exact schema-drift signal fixtures structurally cannot catch.
        resp = _well_formed()
        resp["errors"] = [{"message": "Field 'subIssueSummary' doesn't exist"}]
        problems = bqs.check_board_query_shape(resp)
        assert any("error" in p.lower() for p in problems)

    def test_missing_data_flagged(self):
        assert bqs.check_board_query_shape({}) != []

    def test_missing_items_connection_flagged(self):
        resp = {"data": {"user": {"projectV2": {}}}}
        problems = bqs.check_board_query_shape(resp)
        assert any("items" in p.lower() for p in problems)

    def test_nodes_not_a_list_flagged(self):
        resp = _well_formed()
        resp["data"]["user"]["projectV2"]["items"]["nodes"] = None
        problems = bqs.check_board_query_shape(resp)
        assert any("nodes" in p.lower() for p in problems)

    def test_missing_pageinfo_flagged(self):
        # board_open_items.py paginates on pageInfo; losing it silently truncates.
        resp = _well_formed()
        del resp["data"]["user"]["projectV2"]["items"]["pageInfo"]
        problems = bqs.check_board_query_shape(resp)
        assert any("pageinfo" in p.lower() for p in problems)

    def test_missing_sub_issues_summary_on_issue_flagged(self):
        # The concrete trigger for #771: subIssuesSummary { total } must be present.
        resp = _well_formed(nodes=[_issue_node(with_sub_summary=False)])
        problems = bqs.check_board_query_shape(resp)
        assert any("subissuessummary" in p.lower() for p in problems)

    def test_sub_summary_total_wrong_type_flagged(self):
        resp = _well_formed()
        resp["data"]["user"]["projectV2"]["items"]["nodes"][0]["content"][
            "subIssuesSummary"
        ] = {"total": "not-an-int"}
        problems = bqs.check_board_query_shape(resp)
        assert any("total" in p.lower() for p in problems)

    def test_no_value_assertions(self):
        # Structure-only: the same well-formed shape passes regardless of the
        # specific board data (counts, titles, states). Two different data sets,
        # both structurally valid, must both pass.
        a = _well_formed(nodes=[_issue_node()])
        b = _well_formed(
            nodes=[
                {
                    "content": {
                        "__typename": "Issue",
                        "number": 999,
                        "title": "totally different",
                        "state": "CLOSED",
                        "url": "z",
                        "labels": {"nodes": [{"name": "role:developer"}]},
                        "subIssuesSummary": {"total": 42},
                    },
                    "fieldValues": {"nodes": []},
                }
            ]
        )
        assert bqs.check_board_query_shape(a) == []
        assert bqs.check_board_query_shape(b) == []

    def test_require_issue_node_coverage_guard(self):
        # A response with zero Issue nodes never exercises the subIssuesSummary
        # check — a vacuous pass. require_issue_node=True flags it so the live
        # test can't silently stop covering the field (annotate-canary-style
        # coverage guard).
        pr_only = _well_formed(
            nodes=[
                {
                    "content": {"__typename": "PullRequest", "number": 5, "isDraft": False},
                    "fieldValues": {"nodes": []},
                }
            ]
        )
        problems = bqs.check_board_query_shape(pr_only, require_issue_node=True)
        assert any("issue" in p.lower() for p in problems)
        # Without the guard, the same response is structurally fine.
        assert bqs.check_board_query_shape(pr_only, require_issue_node=False) == []


class TestLiveQueriesRegistry:
    def test_registry_includes_real_board_query(self):
        names = [q.name for q in bqs.LIVE_QUERIES]
        assert "board_open_items" in names

    def test_board_query_is_the_real_one(self):
        # Single source of truth: a field rename in board_open_items.py must
        # flow into the live smoke automatically, not via a hand-copied string.
        entry = next(q for q in bqs.LIVE_QUERIES if q.name == "board_open_items")
        assert "subIssuesSummary" in entry.query
        assert "projectV2" in entry.query


class _FakeProc:
    def __init__(self, returncode, stdout="", stderr=""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


class TestRunGraphqlWithRetry:
    _Q = bqs.LiveQuery(name="t", query="query{x}", variables={"owner": "o"})

    def test_returns_parsed_json_on_success(self):
        runner = lambda cmd: _FakeProc(0, stdout='{"data": {"ok": true}}')
        sleeps = []
        out = bqs.run_graphql_with_retry(
            self._Q, _runner=runner, _sleep=sleeps.append
        )
        assert out == {"data": {"ok": True}}
        assert sleeps == []  # no retry on first-try success

    def test_retries_transient_failure_then_succeeds(self):
        calls = {"n": 0}

        def runner(cmd):
            calls["n"] += 1
            if calls["n"] < 3:
                return _FakeProc(1, stderr="503 transient")
            return _FakeProc(0, stdout='{"data": {}}')

        sleeps = []
        out = bqs.run_graphql_with_retry(
            self._Q, _runner=runner, _sleep=sleeps.append
        )
        assert out == {"data": {}}
        assert sleeps == [1.0, 2.0]  # exponential backoff between the 3 attempts

    def test_raises_after_exhausting_attempts(self):
        runner = lambda cmd: _FakeProc(1, stderr="persistent failure")
        sleeps = []
        with pytest.raises(RuntimeError, match="failed after 4 attempts"):
            bqs.run_graphql_with_retry(
                self._Q, _runner=runner, _sleep=sleeps.append
            )
        assert sleeps == [1.0, 2.0, 4.0]  # backed off between each of the 4 tries

    def test_passes_variables_as_field_flags(self):
        seen = {}

        def runner(cmd):
            seen["cmd"] = cmd
            return _FakeProc(0, stdout="{}")

        bqs.run_graphql_with_retry(self._Q, _runner=runner, _sleep=lambda _: None)
        assert "-F" in seen["cmd"] and "owner=o" in seen["cmd"]


@REQUIRES_LIVE_GH
@pytest.mark.live
class TestBoardQueryLiveSmoke:
    """Runs every registered live query against the real GitHub API.

    Catches: schema drift, a renamed/removed field, a field typo in the live
    query, and auth/scope regressions. Does NOT check data correctness (counts,
    which items are where) — that stays the manual live-smoke convention
    (feedback_live_integration_smoke.md). Advisory/nightly only; promotion to a
    required check is a separate, deliberate decision (Issue #771 ACs).
    """

    def test_every_live_query_is_schema_valid(self):
        failures = []
        for q in bqs.LIVE_QUERIES:
            parsed = bqs.run_graphql_with_retry(q)
            problems = bqs.check_board_query_shape(
                parsed, require_issue_node=q.require_issue_node
            )
            if problems:
                failures.append(f"[{q.name}]\n  " + "\n  ".join(problems))
        assert not failures, (
            "live board GraphQL query is no longer schema-valid (drift / field "
            "typo / removed field) — NOT a data-state change:\n"
            + "\n".join(failures)
        )
