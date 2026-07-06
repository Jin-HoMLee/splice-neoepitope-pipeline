"""Tests for scripts/pm/request_review.sh - deterministic In-review flip (Issue #1063).

The wrapper posts `@claude review` on a PR and flips the PR card + every linked
Issue card on board #9 to "In review", independent of the version-fragile
`post_gh_pr_review_request.py` PostToolUse hook.

For deterministic offline tests we put a stub `gh` first on PATH. It is driven by
a per-test fixture directory ($RR_DIR) and emulates the real gh surfaces the
script touches, applying the script's OWN `--jq` filter to the read fixtures
(exactly as real gh would, so the @tsv parsing is exercised for real):

  - `gh pr view <pr> --json ...`        -> emits $RR_DIR/pr.json (or fails if
                                            $RR_PRVIEW_FAIL set).
  - `gh pr comment ...`                 -> records argv to comments.txt; fails if
                                            $RR_COMMENT_FAIL set.
  - `gh project item-add ...`           -> records argv to itemadds.txt, returns
                                            an id; fails if $RR_ITEMADD_FAIL set.
  - `gh api graphql` *read*             -> parses issueOrPullRequest(number: N),
                                            applies the script's --jq to
                                            $RR_DIR/card_<N>.json (missing => no
                                            card); fails for numbers in
                                            $RR_QUERY_FAIL_NUMS.
  - `gh api graphql` *mutation*         -> records full argv to mutations.txt so a
                                            test can assert which item ids were
                                            flipped; fails for item ids in
                                            $RR_MUTATION_FAIL_ITEMS.

No real board mutation or comment is ever performed.
"""
import json
import os
import subprocess
import tempfile
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "pm" / "request_review.sh"

# Canonical board-#9 In review option id (the assertion target - pinned here as
# the independent source against a silent fat-finger in the script).
IN_REVIEW_OPTION = "df73e18b"

_GH_STUB = r"""#!/usr/bin/env bash
sub="${1:-} ${2:-}"
case "$sub" in
  "pr view")
    if [ -n "${RR_PRVIEW_FAIL:-}" ]; then
      printf '%s\n' "${RR_PRVIEW_FAIL_MSG:-gh: no pull requests found}" >&2
      exit 1
    fi
    cat "$RR_DIR/pr.json"
    exit 0
    ;;
  "pr comment")
    printf '%s\n' "$@" >> "$RR_DIR/comments.txt"
    if [ -n "${RR_COMMENT_FAIL:-}" ]; then
      printf '%s\n' "${RR_COMMENT_FAIL_MSG:-gh: comment failed}" >&2
      exit 1
    fi
    echo "https://github.com/o/r/pull/${3}#issuecomment-1"
    exit 0
    ;;
  "project item-add")
    printf '%s\n' "$@" >> "$RR_DIR/itemadds.txt"
    if [ -n "${RR_ITEMADD_FAIL:-}" ]; then
      printf '%s\n' "${RR_ITEMADD_FAIL_MSG:-gh: item-add failed}" >&2
      exit 1
    fi
    echo '{"id":"PVTI_added"}'
    exit 0
    ;;
  "api graphql")
    qval=""; jqf=""; prev=""
    for a in "$@"; do
      case "$prev" in --jq) jqf="$a" ;; esac
      case "$a" in query=*) qval="${a#query=}" ;; esac
      prev="$a"
    done
    if printf '%s' "$qval" | grep -q "updateProjectV2ItemFieldValue"; then
      printf '%s\n' "$@" >> "$RR_DIR/mutations.txt"
      item=""
      for a in "$@"; do case "$a" in i=*) item="${a#i=}" ;; esac; done
      case ",${RR_MUTATION_FAIL_ITEMS:-}," in
        *",$item,"*) printf 'gh: mutation failed for %s\n' "$item" >&2; exit 1 ;;
      esac
      printf '{"data":{"updateProjectV2ItemFieldValue":{"projectV2Item":{"id":"%s"}}}}\n' "$item"
      exit 0
    fi
    num="$(printf '%s' "$qval" | sed -n 's/.*issueOrPullRequest(number: \([0-9]*\)).*/\1/p')"
    case ",${RR_QUERY_FAIL_NUMS:-}," in
      *",$num,"*) printf 'gh: board query failed for %s\n' "$num" >&2; exit 1 ;;
    esac
    fix="$RR_DIR/card_${num}.json"
    if [ -f "$fix" ]; then
      jq -r "$jqf" "$fix"
    else
      printf '%s' '{"data":{"repository":{"issueOrPullRequest":{"title":"","projectItems":{"nodes":[]}}}}}' | jq -r "$jqf"
    fi
    exit 0
    ;;
esac
printf 'stub gh: unexpected call: %s\n' "$*" >&2
exit 99
"""


def _card(item_id="PVTI_1", current="Ready for review", title="t", on_board=True):
    """A resolve_card read envelope (real-gh-shaped, issueOrPullRequest union)."""
    nodes = []
    if on_board:
        nodes = [{
            "id": item_id,
            "project": {"number": 9},
            "fieldValues": {"nodes": [
                {"name": current, "field": {"name": "Status"}},
            ]},
        }]
    return {"data": {"repository": {"issueOrPullRequest": {
        "title": title,
        "projectItems": {"nodes": nodes},
    }}}}


def _pr(number=500, linked=(401, 402),
        owner="Jin-HoMLee", repo="splice-neoepitope-pipeline", title="wrapper PR"):
    return {
        "url": f"https://github.com/{owner}/{repo}/pull/{number}",
        "number": number,
        "title": title,
        "closingIssuesReferences": [{"number": n} for n in linked],
    }


def _run(args, *, pr=None, cards=None, env_extra=None):
    """Run request_review.sh with a stub gh. `cards` maps number -> _card() dict
    (a missing number => no card on the board). Returns the CompletedProcess plus
    .comments / .itemadds / .mutations capture text."""
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        gh = d / "gh"
        gh.write_text(_GH_STUB)
        gh.chmod(0o755)
        rr = d / "fx"
        rr.mkdir()
        (rr / "pr.json").write_text(json.dumps(pr if pr is not None else _pr()))
        for num, card in (cards or {}).items():
            (rr / f"card_{num}.json").write_text(json.dumps(card))
        env = {
            **os.environ,
            "PATH": f"{d}:{os.environ['PATH']}",
            "RR_DIR": str(rr),
        }
        env.update(env_extra or {})
        proc = subprocess.run(
            ["bash", str(SCRIPT), *args],
            env=env, capture_output=True, text=True,
        )
        proc.comments = (rr / "comments.txt").read_text() if (rr / "comments.txt").exists() else ""
        proc.itemadds = (rr / "itemadds.txt").read_text() if (rr / "itemadds.txt").exists() else ""
        proc.mutations = (rr / "mutations.txt").read_text() if (rr / "mutations.txt").exists() else ""
        return proc


# --- Usage / arg contract (exit 2 before any gh call) ------------------------

def test_no_args_usage_exit_2():
    r = _run([])
    assert r.returncode == 2
    assert "Usage" in r.stderr


def test_non_numeric_pr_exit_2():
    r = _run(["abc"])
    assert r.returncode == 2
    assert "numeric" in r.stderr.lower()


def test_unknown_flag_exit_2():
    r = _run(["500", "--bogus"])
    assert r.returncode == 2
    assert "unknown flag" in r.stderr.lower()


def test_extra_positional_exit_2():
    r = _run(["500", "600"])
    assert r.returncode == 2


# --- Happy path: comment posted + all three cards flipped --------------------

def _happy_cards():
    return {
        500: _card(item_id="PVTI_pr", current="Ready for review"),
        401: _card(item_id="PVTI_401", current="In progress"),
        402: _card(item_id="PVTI_402", current="In progress"),
    }


def test_happy_path_posts_and_flips_all():
    r = _run(["500"], cards=_happy_cards())
    assert r.returncode == 0, r.stderr
    # comment posted with the exact canonical trigger
    assert "@claude review" in r.comments
    # all three cards flipped, each carrying the canonical In-review option id
    for item in ("PVTI_pr", "PVTI_401", "PVTI_402"):
        assert f"i={item}" in r.mutations, r.mutations
    assert f"o={IN_REVIEW_OPTION}" in r.mutations
    # exactly one mutation per card (3 total) - each records its own `i=<item>` arg
    assert r.mutations.count("i=PVTI_") == 3


def test_flip_only_skips_comment_still_flips():
    r = _run(["500", "--flip-only"], cards=_happy_cards())
    assert r.returncode == 0, r.stderr
    assert r.comments == ""  # no comment posted
    assert r.mutations.count("i=PVTI_") == 3  # cards still flipped


# --- Idempotency + skip rules ------------------------------------------------

def test_already_in_review_is_noop():
    cards = _happy_cards()
    cards[401] = _card(item_id="PVTI_401", current="In review")
    r = _run(["500", "--flip-only"], cards=cards)
    assert r.returncode == 0, r.stderr
    assert "i=PVTI_401" not in r.mutations  # already there -> no mutation
    assert "i=PVTI_pr" in r.mutations
    assert "i=PVTI_402" in r.mutations


def test_epic_and_done_cards_are_skipped():
    cards = _happy_cards()
    cards[401] = _card(item_id="PVTI_401", current="Epic")
    cards[402] = _card(item_id="PVTI_402", current="Done")
    r = _run(["500", "--flip-only"], cards=cards)
    assert r.returncode == 0, r.stderr  # a skip is not an error
    assert "i=PVTI_401" not in r.mutations
    assert "i=PVTI_402" not in r.mutations
    assert "i=PVTI_pr" in r.mutations   # the PR card still flips
    assert "Epic" in r.stdout
    assert "Done" in r.stdout


# --- PR add-if-missing (post_gh_pr_create hook may also have not fired) -------

def test_pr_not_on_board_gets_added_then_flipped():
    cards = {  # no card_500 -> PR not on board
        401: _card(item_id="PVTI_401", current="In progress"),
        402: _card(item_id="PVTI_402", current="In progress"),
    }
    r = _run(["500", "--flip-only"], cards=cards)
    assert r.returncode == 0, r.stderr
    assert "item-add" in r.itemadds          # PR was added to board #9
    assert "i=PVTI_added" in r.mutations     # then flipped (the added item id)
    assert "i=PVTI_401" in r.mutations
    assert "i=PVTI_402" in r.mutations


def test_untracked_repo_pr_not_added():
    # A PR in an unrelated repo must not be boarded onto #9 by the fallback.
    pr = _pr(owner="someone-else", repo="other-repo", linked=())
    r = _run(["500", "--flip-only"], pr=pr, cards={})  # no card_500
    assert r.returncode == 1  # PR card unresolved + not added -> surfaced
    assert r.itemadds == ""
    assert "no card" in r.stderr.lower()


# --- Per-card fail-safe: one failure does not abort the rest ------------------

def test_one_issue_query_failure_does_not_abort_others():
    r = _run(["500", "--flip-only"], cards=_happy_cards(),
             env_extra={"RR_QUERY_FAIL_NUMS": "401"})
    assert r.returncode == 1  # overall failure surfaced
    assert "board query failed" in r.stderr.lower()
    # the PR card and the OTHER issue still flipped despite #401 failing
    assert "i=PVTI_pr" in r.mutations
    assert "i=PVTI_402" in r.mutations
    assert "i=PVTI_401" not in r.mutations


def test_one_mutation_failure_does_not_abort_others():
    r = _run(["500", "--flip-only"], cards=_happy_cards(),
             env_extra={"RR_MUTATION_FAIL_ITEMS": "PVTI_401"})
    assert r.returncode == 1
    assert "flip to 'in review' failed" in r.stderr.lower()
    assert "i=PVTI_pr" in r.mutations
    assert "i=PVTI_402" in r.mutations


def test_comment_failure_still_flips_cards():
    r = _run(["500"], cards=_happy_cards(), env_extra={"RR_COMMENT_FAIL": "1"})
    assert r.returncode == 1  # comment failure surfaced
    assert "failed to post" in r.stderr.lower()
    # but every card still flipped
    assert r.mutations.count("i=PVTI_") == 3


# --- Linked-issue edge cases -------------------------------------------------

def test_no_linked_issues_flips_pr_only():
    pr = _pr(linked=())
    cards = {500: _card(item_id="PVTI_pr", current="Ready for review")}
    r = _run(["500", "--flip-only"], pr=pr, cards=cards)
    assert r.returncode == 0, r.stderr
    assert "i=PVTI_pr" in r.mutations
    assert r.mutations.count("i=PVTI_") == 1
    assert "no linked issues" in r.stdout.lower()


def test_linked_issue_not_on_board_is_warned_not_added():
    cards = _happy_cards()
    del cards[402]  # #402 has no card
    r = _run(["500", "--flip-only"], cards=cards)
    assert r.returncode == 1  # surfaced
    assert "no card on board" in r.stderr.lower()
    assert r.itemadds == ""  # issues are never auto-added
    assert "i=PVTI_pr" in r.mutations
    assert "i=PVTI_401" in r.mutations


# --- PR resolution failure ---------------------------------------------------

def test_pr_view_failure_exits_1():
    r = _run(["500"], env_extra={"RR_PRVIEW_FAIL": "1"})
    assert r.returncode == 1
    assert "could not resolve pr" in r.stderr.lower()
    assert r.mutations == ""  # nothing attempted
