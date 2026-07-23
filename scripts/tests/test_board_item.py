# scripts/tests/test_board_item.py
#
# Issue #1151: a board mutation that can only report success.
#
# Two roles independently derived the same silent no-op within 24 hours, because
# every failure mode of a board-item lookup returns EMPTY WITH EXIT 0:
#
#   1. the phantom key - `gh issue view N --json projectItems` has no `id` key at
#      all, so `--jq '.projectItems[].id'` yields empty, exit 0;
#   2. the item genuinely is not on the board;
#   3. the board read was TRUNCATED (`gh project item-list` defaults to 30 rows
#      against a board that held 286 items on 2026-07-23).
#
# All three are indistinguishable at the call site, so a resolver that returns
# None/empty hands a write a value it cannot vouch for. The helper must RAISE.
#
# The tests below are the two-directional falsifier the Issue's AC-3 demands: a
# resolver that raised on everything, or returned a value for everything, must go
# red. Fakes stand in for `gh` so the assertions are hermetic; the live-API schema
# facts are pinned separately in the `live`-marked tests at the bottom.
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


def _page(items, *, total, has_next=False, cursor=None):
    """One projectV2 items page in the shape the real GraphQL query returns."""
    return {
        "data": {
            "user": {
                "projectV2": {
                    "items": {
                        "totalCount": total,
                        "pageInfo": {"hasNextPage": has_next, "endCursor": cursor},
                        "nodes": [
                            {
                                "id": i["id"],
                                "content": {
                                    "number": i["number"],
                                    "repository": {
                                        "nameWithOwner": i.get(
                                            "repo", "Jin-HoMLee/splice-neoepitope-pipeline"
                                        )
                                    },
                                },
                            }
                            for i in items
                        ],
                    }
                }
            }
        }
    }


def _fake_gh(pages):
    """A gh() stand-in that returns the given pages in order."""
    calls = {"n": 0}

    def _gh(*args, **kwargs):
        page = pages[calls["n"]]
        calls["n"] += 1
        return page

    _gh.calls = calls
    return _gh


def test_resolves_a_present_issue_to_its_board_item_id():
    gh = _fake_gh([
        _page([{"id": "PVTI_present", "number": 1151}], total=1),
    ])

    assert bi.resolve_board_item_id(1151, _gh=gh) == "PVTI_present"


def test_raises_for_an_absent_issue_instead_of_returning_empty():
    # The complete board was read (1 of 1) and 9999 simply is not on it. This is
    # the ONLY case where "not on the board" is a trustworthy answer - and even
    # here the helper raises rather than returning None, so the value can never
    # reach a mutation.
    gh = _fake_gh([
        _page([{"id": "PVTI_other", "number": 1151}], total=1),
    ])

    with pytest.raises(bi.BoardItemNotFound):
        bi.resolve_board_item_id(9999, _gh=gh)


def test_follows_pagination_to_find_an_item_beyond_the_first_page():
    # THE case that motivated the Issue. #1135 sat at index 125 of 286 on the
    # live board, so a single-page read answers "not on the board" with exit 0.
    # A helper that only asserted non-emptiness would confidently return that
    # wrong answer.
    gh = _fake_gh([
        _page([{"id": "PVTI_a", "number": 1151}], total=2, has_next=True, cursor="CUR1"),
        _page([{"id": "PVTI_b", "number": 1135}], total=2),
    ])

    assert bi.resolve_board_item_id(1135, _gh=gh) == "PVTI_b"
    assert gh.calls["n"] == 2, "must have fetched the second page"


def test_incomplete_read_raises_incomplete_not_not_found():
    # The distinction the whole Issue turns on. The API claims 286 items but
    # delivered 1 and said there are no more pages - a degraded/truncated read.
    # Reporting that as "not on the board" is the false answer that cost a 2h
    # blackout, so it must be a DIFFERENT, louder error than a genuine absence.
    gh = _fake_gh([
        _page([{"id": "PVTI_a", "number": 1151}], total=286),
    ])

    with pytest.raises(bi.BoardReadIncomplete):
        bi.resolve_board_item_id(9999, _gh=gh)


def test_incomplete_read_is_not_masked_by_a_lucky_hit():
    # Subtle: if the target happens to be on the page we did get, returning it is
    # still correct - the id is real. Only the NEGATIVE conclusion needs a
    # complete read. This pins that we do not over-raise and break every caller.
    gh = _fake_gh([
        _page([{"id": "PVTI_a", "number": 1151}], total=286),
    ])

    assert bi.resolve_board_item_id(1151, _gh=gh) == "PVTI_a"


def test_matched_node_without_an_id_raises_rather_than_yielding_none():
    # The phantom key, at the node level: if a schema change ever drops `id` the
    # match must fail LOUDLY. Returning node.get("id") here would hand None into
    # a mutation - literally the #1151 defect, reintroduced inside its own fix.
    gh = _fake_gh([{
        "data": {"user": {"projectV2": {"items": {
            "totalCount": 1,
            "pageInfo": {"hasNextPage": False, "endCursor": None},
            "nodes": [{"content": {"number": 1151}}],  # no "id"
        }}}}
    }])

    with pytest.raises(bi.BoardLookupError):
        bi.resolve_board_item_id(1151, _gh=gh)


def test_disambiguates_the_same_issue_number_in_two_repos():
    # Board 9 spans TWO repos with independent numbering, and the collision is
    # live: #183 and #193 each exist in both splice-neoepitope-pipeline and
    # claude-personas-splice-neoepitope-pipeline, both carded on board 9
    # (verified 2026-07-23). Matching on content.number alone returns whichever
    # page-order happens to yield first - a silent wrong answer, which is the
    # very defect this module exists to prevent.
    gh = _fake_gh([
        _page(
            [
                {"id": "PVTI_personas", "number": 183,
                 "repo": "Jin-HoMLee/claude-personas-splice-neoepitope-pipeline"},
                {"id": "PVTI_pipeline", "number": 183,
                 "repo": "Jin-HoMLee/splice-neoepitope-pipeline"},
            ],
            total=2,
        ),
    ])

    got = bi.resolve_board_item_id(
        183, repo="Jin-HoMLee/splice-neoepitope-pipeline", _gh=gh
    )

    assert got == "PVTI_pipeline", "must not return the personas-repo card for #183"


def test_defaults_to_the_hardened_gh_client():
    # Without this the helper is unusable outside tests. gh_client.gh carries the
    # retry/backoff and the "no --jq inside gh" house rule (#1017) - re-rolling a
    # bare subprocess call here would reintroduce the raw-jq path that produces
    # the phantom-key empty-with-exit-0 in the first place.
    #
    # Identity comparison is impossible across two by-path loads (each exec makes
    # a fresh function object), so pin what actually matters: it resolves to
    # gh_client's `gh`, not some locally re-rolled subprocess call.
    resolved = bi._default_gh()

    assert resolved.__name__ == "gh"
    assert resolved.__module__ == "gh_client"
    assert os.path.samefile(
        resolved.__globals__["__file__"],
        os.path.join(_HERE, "..", "pm", "gh_client.py"),
    )


def test_unexpected_response_shape_raises_a_board_lookup_error():
    # A wrong owner/project number returns projectV2: null. Without this the
    # caller gets a bare TypeError from deep inside a dict walk.
    gh = _fake_gh([{"data": {"user": {"projectV2": None}}}])

    with pytest.raises(bi.BoardLookupError):
        bi.resolve_board_item_id(1151, _gh=gh)


# ---------------------------------------------------------------------------
# AC-4: pin the live-API schema facts this helper exists to work around.
#
# These are the facts that made the defect invisible. If GitHub ever changes
# them the workaround should be revisited - in EITHER direction. If `id` starts
# appearing under `gh issue view --json projectItems`, the one-call join gets
# simpler and these tests should fail loudly to tell us so, rather than leaving
# a now-unnecessary paginating helper in place forever.
#
# Live-gated: they skip without a `read:project`-scoped gh, matching the
# tools/ci convention rather than inventing a second one.
# ---------------------------------------------------------------------------


def _gh_json(*args):
    out = subprocess.run(["gh", *args], capture_output=True, text=True, check=True)
    return json.loads(out.stdout)


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
        "gh issue view --json projectItems now exposes 'id'. The phantom-key "
        f"workaround in board_item.py may be simplifiable. Keys seen: {sorted(keys)}"
    )


@REQUIRES_LIVE_GH
@pytest.mark.live
def test_live_schema_project_item_list_joins_item_id_to_issue_number():
    payload = _gh_json(
        "project", "item-list", "9",
        "--owner", "Jin-HoMLee", "--format", "json", "--limit", "5",
    )
    item = payload["items"][0]

    # The join that DOES work in one call (the Issue's correction to its parent).
    assert "id" in item
    assert "number" in item["content"]
    # ...and the key that is genuinely phantom.
    assert "id" not in item["content"], (
        "gh project item-list now returns content.id; revisit board_item.py"
    )


@REQUIRES_LIVE_GH
@pytest.mark.live
def test_live_resolver_agrees_with_an_independent_source():
    # Differential oracle: resolve via our paginating GraphQL helper, and via the
    # unrelated `gh project item-list` path. Two independent routes to the same
    # id. A silent no-op in either shows up as disagreement rather than as a
    # confident wrong answer.
    ours = bi.resolve_board_item_id(1151)

    payload = _gh_json(
        "project", "item-list", "9",
        "--owner", "Jin-HoMLee", "--format", "json", "--limit", "1000",
    )
    theirs = [i["id"] for i in payload["items"] if i["content"].get("number") == 1151]

    assert theirs, "independent source could not find #1151; the oracle itself is broken"
    assert ours == theirs[0]
