# scripts/tests/test_board_item.py
#
# Issue #1151: a board mutation that can only report success.
#
# Two roles independently derived the same silent no-op within 24 hours, because
# every failure mode of a board-item lookup returns EMPTY WITH EXIT 0:
#
#   1. the phantom key - `gh issue view N --json projectItems` exposes only
#      `status` and `title`, so `--jq '.projectItems[].id'` yields empty, exit 0;
#   2. the item genuinely is not on the board;
#   3. the read did not cover the whole board.
#
# All three are indistinguishable at the call site, so a resolver that returns
# None/empty hands a write a value it cannot vouch for. The helper must RAISE.
#
# WHY THE PER-ISSUE PATH (`Issue.projectItems`) RATHER THAN A BOARD SCAN:
# `ProjectV2.items` defaults to `archivedStates: [NOT_ARCHIVED]` and its
# `totalCount` is scoped to that same subset - so a board scan silently omits
# archived cards AND its completeness assertion passes anyway. Measured on board
# 9 (2026-07-23): 289 unarchived, 1,163 archived. The first draft of this helper
# scanned the board and answered "#569 is not on project (complete read of 289
# items)" for a card that exists. `Issue.projectItems` defaults
# `includeArchived: true` and is per-issue, so there is no board-wide read to
# truncate and no archived blind spot - both failure modes become
# unrepresentable rather than asserted-against.
#
# The tests below are the two-directional falsifier AC-3 demands: a resolver that
# raised on everything, or returned a value for everything, must go red.
import importlib.util
import json
import os
import subprocess

import pytest

# Load both modules by explicit path rather than mutating sys.path.
#
# The usual `sys.path.insert` convention in this directory is unsafe HERE because
# both dirs we need contain modules whose names collide with ones other suites
# import: `scripts/pm/recheck_parent_status.py` and `tools/ci/recheck_*.py`.
# Prepending either dir shadows them for the whole pytest session, which broke
# two unrelated tests in tools/ci when this file was first added - a failure that
# only appears in a FULL five-directory run, never when this file runs alone.
_HERE = os.path.dirname(__file__)


def _load(name, *parts):
    path = os.path.join(_HERE, *parts)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


bi = _load("board_item", "..", "pm", "board_item.py")
REQUIRES_LIVE_GH = _load("_live_gh", "..", "..", "tools", "ci", "_live_gh.py").REQUIRES_LIVE_GH


def _resp(nodes, *, has_next=False, issue_exists=True, kind="issue"):
    """A `repository.<kind>.projectItems` payload in the real response shape."""
    if not issue_exists:
        return {"data": {"repository": {kind: None}}}
    return {
        "data": {
            "repository": {
                kind: {
                    "projectItems": {
                        "pageInfo": {"hasNextPage": has_next, "endCursor": None},
                        "nodes": nodes,
                    }
                }
            }
        }
    }


def _item(item_id, project_number=9, archived=False, status=None):
    field_values = []
    if status is not None:
        field_values.append({"name": status, "field": {"name": "Status"}})
    return {
        "id": item_id,
        "isArchived": archived,
        "project": {"number": project_number},
        "fieldValues": {"nodes": field_values},
    }


def _fake_gh(response):
    calls = {"n": 0, "args": []}

    def _gh(*args, **kwargs):
        calls["n"] += 1
        calls["args"].append(args)
        return response

    _gh.calls = calls
    return _gh


def test_resolves_a_present_issue_to_its_board_item_id():
    gh = _fake_gh(_resp([_item("PVTI_present")]))

    assert bi.resolve_board_item_id(1151, _gh=gh) == "PVTI_present"


def test_resolves_an_archived_card():
    # THE case that killed the board-scan design. #569 is archived on board 9 and
    # a board scan reports it absent with a passing completeness check. A helper
    # that cannot resolve an archived card is silently wrong for 1,163 of the
    # board's 1,452 items.
    gh = _fake_gh(_resp([_item("PVTI_archived", archived=True)]))

    assert bi.resolve_board_item_id(569, _gh=gh) == "PVTI_archived"


def test_raises_when_the_issue_has_no_card_on_this_board():
    # The issue exists and carries project items, but none on board 9. This is
    # the only trustworthy "no" - and even here it raises rather than returning
    # None, so the value can never reach a mutation.
    gh = _fake_gh(_resp([_item("PVTI_other_board", project_number=42)]))

    with pytest.raises(bi.BoardItemNotFound):
        bi.resolve_board_item_id(1151, _gh=gh)


def test_raises_when_the_issue_has_no_project_items_at_all():
    gh = _fake_gh(_resp([]))

    with pytest.raises(bi.BoardItemNotFound):
        bi.resolve_board_item_id(1151, _gh=gh)


def test_truncated_project_items_raise_incomplete_not_not_found():
    # projectItems is tiny in practice, but if the connection ever reports more
    # pages we must not conclude absence from a partial list. Distinct from
    # BoardItemNotFound: "I did not see it" is not "it is not there".
    gh = _fake_gh(_resp([_item("PVTI_other_board", project_number=42)], has_next=True))

    with pytest.raises(bi.BoardReadIncomplete):
        bi.resolve_board_item_id(1151, _gh=gh)


def test_truncation_does_not_mask_a_hit():
    # A hit from a partial list is still trustworthy - the id matched is real.
    # Only the NEGATIVE conclusion needs completeness. Pins that we do not
    # over-raise and break every caller.
    gh = _fake_gh(_resp([_item("PVTI_present")], has_next=True))

    assert bi.resolve_board_item_id(1151, _gh=gh) == "PVTI_present"


def test_missing_issue_raises_lookup_error_not_not_found():
    # A wrong repo (or a PR number) returns issue: null. That is "I cannot
    # answer", not "it has no card" - conflating them is how a cross-repo typo
    # becomes a confident wrong answer.
    gh = _fake_gh(_resp([], issue_exists=False))

    with pytest.raises(bi.BoardLookupError) as exc:
        bi.resolve_board_item_id(1151, _gh=gh)
    assert not isinstance(exc.value, bi.BoardItemNotFound)


def test_matched_node_without_an_id_raises_rather_than_yielding_none():
    # The phantom key at node level: returning node.get("id") would hand None
    # into a mutation - the #1151 defect reintroduced inside its own fix.
    gh = _fake_gh(_resp([{"isArchived": False, "project": {"number": 9}}]))

    with pytest.raises(bi.BoardLookupError):
        bi.resolve_board_item_id(1151, _gh=gh)


def test_query_is_scoped_to_the_named_repo():
    # Board 9 spans two repos with INDEPENDENT numbering: a complete read (1,453
    # items, archived included) finds 56 numbers that are an Issue in BOTH repos,
    # #37 among them. The per-issue path makes this structurally safe - the repo
    # is in the query - but only if we actually thread it through, so pin that.
    #
    # Do not evidence the collision by scanning board content numbers: issues and
    # PRs share a number space, so #183 reads as a collision while being an Issue
    # in the pipeline repo and a PR in personas.
    gh = _fake_gh(_resp([_item("PVTI_x")]))

    bi.resolve_board_item_id(183, repo="Jin-HoMLee/claude-personas-splice-neoepitope-pipeline", _gh=gh)

    sent = " ".join(gh.calls["args"][0])
    assert "claude-personas-splice-neoepitope-pipeline" in sent
    assert "includeArchived" in sent, "archived cards must be in scope"


def test_resolve_board_item_returns_id_and_status():
    # The hooks need the current Status alongside the id (to stay idempotent and
    # to avoid overwriting a parked Epic / terminal Done card). Serving that here
    # is what lets them drop their own copy of this query.
    gh = _fake_gh(_resp([_item("PVTI_x", status="In progress")]))

    assert bi.resolve_board_item(1151, _gh=gh) == ("PVTI_x", "In progress")


def test_resolve_board_item_reports_none_status_when_unset():
    gh = _fake_gh(_resp([_item("PVTI_x")]))

    assert bi.resolve_board_item(1151, _gh=gh) == ("PVTI_x", None)


def test_pull_request_kind_queries_the_pull_request_field():
    # A PR has its own card. `repository.issue(N)` is Issue-only, so the field
    # name must switch - without this the hooks cannot delegate.
    gh = _fake_gh(_resp([_item("PVTI_pr")], kind="pullRequest"))

    assert bi.resolve_board_item(42, kind="pullRequest", _gh=gh) == ("PVTI_pr", None)
    assert "pullRequest(number:" in " ".join(gh.calls["args"][0])


def test_unknown_kind_is_rejected_rather_than_interpolated():
    # `kind` reaches a GraphQL field name by string building; an allow-list keeps
    # a caller typo (or anything worse) from being interpolated into the query.
    with pytest.raises(ValueError):
        bi.resolve_board_item(1, kind="issue) { x } #", _gh=_fake_gh(_resp([])))


def test_gh_failure_surfaces_as_a_board_lookup_error():
    # gh_client.GhError subclasses CalledProcessError. Without conversion the
    # caller gets an exception that is NOT catchable as BoardLookupError, which
    # breaks the module's single promise: every outcome is an id or a
    # BoardLookupError. (Hermetic mirror of the live PR-number test.)
    def _boom(*args, **kwargs):
        raise subprocess.CalledProcessError(
            1, "gh", stderr="Could not resolve to an Issue with the number of 183."
        )

    with pytest.raises(bi.BoardLookupError) as exc:
        bi.resolve_board_item_id(183, _gh=_boom)

    assert not isinstance(exc.value, bi.BoardItemNotFound)
    assert "PR" in str(exc.value), "should hint at the issue/PR number-space overlap"


def test_defaults_to_the_hardened_gh_client():
    # gh_client.gh carries the retry/backoff and the "no --jq inside gh" house
    # rule (#1017); re-rolling a bare subprocess call here would reintroduce the
    # raw-jq path that produces the phantom-key empty-with-exit-0.
    #
    # Identity comparison is impossible across two by-path loads, so pin what
    # matters: it resolves to gh_client's `gh`.
    resolved = bi._default_gh()

    assert resolved.__name__ == "gh"
    assert resolved.__module__ == "gh_client"
    assert os.path.samefile(
        resolved.__globals__["__file__"],
        os.path.join(_HERE, "..", "pm", "gh_client.py"),
    )


# ---------------------------------------------------------------------------
# AC-4: pin the live-API facts this helper depends on, so a schema change fails
# loudly instead of silently. Live-gated, matching the tools/ci convention.
# ---------------------------------------------------------------------------


def _gh_json(*args):
    out = subprocess.run(["gh", *args], capture_output=True, text=True, check=True)
    return json.loads(out.stdout)


@REQUIRES_LIVE_GH
@pytest.mark.live
def test_live_schema_archived_defaults_still_differ_between_the_two_paths():
    # The fact the whole design rests on. If `ProjectV2.items` ever defaults to
    # including archived, or `Issue.projectItems` stops doing so, revisit.
    payload = _gh_json(
        "api", "graphql", "-f",
        'query={ i: __type(name:"Issue"){ fields{ name args{ name defaultValue } } } '
        'p: __type(name:"ProjectV2"){ fields{ name args{ name defaultValue } } } }',
    )

    def _args(type_payload, field):
        return {
            a["name"]: a["defaultValue"]
            for f in type_payload["fields"] if f["name"] == field
            for a in f["args"]
        }

    issue_args = _args(payload["data"]["i"], "projectItems")
    project_args = _args(payload["data"]["p"], "items")

    assert issue_args.get("includeArchived") == "true", (
        "Issue.projectItems no longer includes archived cards by default; "
        "board_item.py would start missing archived items"
    )
    assert project_args.get("archivedStates") == "[NOT_ARCHIVED]", (
        "ProjectV2.items archived default changed; the board-scan blind spot "
        "documented in this module may no longer apply"
    )


@REQUIRES_LIVE_GH
@pytest.mark.live
def test_live_schema_gh_issue_view_still_has_no_board_item_id():
    payload = _gh_json(
        "issue", "view", "1151",
        "--repo", "Jin-HoMLee/splice-neoepitope-pipeline",
        "--json", "projectItems",
    )
    keys = set(payload["projectItems"][0])

    assert "id" not in keys, (
        "gh issue view --json projectItems now exposes 'id'; the workaround in "
        f"board_item.py may be simplifiable. Keys seen: {sorted(keys)}"
    )


@REQUIRES_LIVE_GH
@pytest.mark.live
def test_live_resolves_an_archived_card_the_board_scan_could_not():
    # End-to-end proof against real data: #569 is archived on board 9.
    assert bi.resolve_board_item_id(569) == "PVTI_lAHOB17eGc4BSomPzguPpXg"


@REQUIRES_LIVE_GH
@pytest.mark.live
def test_live_cross_repo_collision_resolves_distinctly():
    # #37 exists as an ISSUE in both repos, both carded on board 9. A complete
    # board read (1,453 items, archived included) found 56 such numbers - the
    # collision is substantial, and invisible in an unarchived-only view.
    #
    # Note the number matters: an earlier draft used #183, which is an Issue in
    # the pipeline repo but a PULL REQUEST in personas. Numbers are shared
    # between issues and PRs, so a "collision" found by scanning board content
    # numbers is not necessarily an Issue/Issue collision.
    pipeline = bi.resolve_board_item_id(37, repo="Jin-HoMLee/splice-neoepitope-pipeline")
    personas = bi.resolve_board_item_id(
        37, repo="Jin-HoMLee/claude-personas-splice-neoepitope-pipeline"
    )

    assert pipeline != personas


@REQUIRES_LIVE_GH
@pytest.mark.live
def test_live_pull_request_number_is_a_lookup_error_not_a_missing_card():
    # repository.issue(N) is Issue-only, so a PR number cannot resolve. That must
    # read as "I cannot answer" rather than "it has no card" - personas #183 is a
    # PR and is genuinely carded on board 9.
    with pytest.raises(bi.BoardLookupError) as exc:
        bi.resolve_board_item_id(
            183, repo="Jin-HoMLee/claude-personas-splice-neoepitope-pipeline"
        )

    assert not isinstance(exc.value, bi.BoardItemNotFound)
