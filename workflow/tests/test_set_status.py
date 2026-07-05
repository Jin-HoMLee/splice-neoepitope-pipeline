"""Tests for scripts/pm/set_status.sh - board #9 status-move wrapper (Issue #1024).

The script maps a human Status name to its board-#9 single-select option id and
performs the updateProjectV2ItemFieldValue mutation. For deterministic offline
tests we put a stub `gh` first on PATH that:

  - `gh repo view ...`         -> prints "TestOwner TestRepo" (owner/name).
  - `gh api graphql` *read*    -> applies the script's own `--jq` filter to the
                                   fixture envelope in $SET_STATUS_READ_FIXTURE,
                                   exactly as real gh would (so the @tsv parsing
                                   in the script is exercised for real).
  - `gh api graphql` *mutation* -> writes its full argv to $SET_STATUS_GH_CAPTURE
                                   so a test can assert which option id was sent.

The stub honours $SET_STATUS_GH_FAIL (non-empty => the read exits non-zero,
emitting $SET_STATUS_GH_FAIL_MSG) to exercise the query-failure / PR-number
branches. No real board mutation is ever performed.

Reviewer ask (PR #1039): lock the exit-code contract (usage / non-numeric /
unknown-status -> 2) and the Status-name -> option-id map against a silent
fat-fingered id. These are the canonical ids, mirrored from the board hooks.
"""
import os
import json
import subprocess
import tempfile
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "pm" / "set_status.sh"

# Canonical Status-name -> option-id map for board #9 (the assertion target -
# this list is the independent source the test pins the script's `case` against).
CANONICAL = {
    "Backlog": "f75ad846",
    "Ready": "61e4505c",
    "In progress": "47fc9ee4",
    "Ready for review": "8bf9192f",
    "In review": "df73e18b",
    "Done": "98236657",
    "Epic": "9f872564",
}

_GH_STUB = r"""#!/usr/bin/env bash
if [ "$1" = "repo" ] && [ "$2" = "view" ]; then
  echo "TestOwner TestRepo"
  exit 0
fi
if [ "$1" = "api" ] && [ "$2" = "graphql" ]; then
  # pull the query value and the --jq filter out of the -f/-flag argv
  qval=""; jqf=""; prev=""
  for a in "$@"; do
    case "$prev" in --jq) jqf="$a" ;; esac
    case "$a" in query=*) qval="${a#query=}" ;; esac
    prev="$a"
  done
  if printf '%s' "$qval" | grep -q "updateProjectV2ItemFieldValue"; then
    # mutation: record argv (so the option id sent is visible) and succeed
    printf '%s\n' "$@" > "$SET_STATUS_GH_CAPTURE"
    echo '{"data":{"updateProjectV2ItemFieldValue":{"projectV2Item":{"id":"PVTI_x"}}}}'
    exit 0
  fi
  # read query
  if [ -n "${SET_STATUS_GH_FAIL:-}" ]; then
    printf '%s\n' "${SET_STATUS_GH_FAIL_MSG:-gh: request failed}" >&2
    exit 1
  fi
  jq -r "$jqf" "$SET_STATUS_READ_FIXTURE"
  exit 0
fi
echo "stub gh: unexpected call: $*" >&2
exit 99
"""


def _envelope(item_id="PVTI_1", current="Backlog", title="Some issue", on_board=True):
    """The GraphQL read envelope the script consumes (real-gh-shaped)."""
    nodes = []
    if on_board:
        nodes = [{
            "id": item_id,
            "project": {"number": 9},
            "fieldValues": {"nodes": [
                {"name": current, "field": {"name": "Status"}},
            ]},
        }]
    return {"data": {"repository": {"issue": {
        "title": title,
        "projectItems": {"nodes": nodes},
    }}}}


def _run(args, *, envelope=None, gh_fail=False, gh_fail_msg=""):
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        gh = d / "gh"
        gh.write_text(_GH_STUB)
        gh.chmod(0o755)
        fix = d / "read.json"
        fix.write_text(json.dumps(envelope if envelope is not None else _envelope()))
        capture = d / "mutation_argv.txt"
        env = {
            **os.environ,
            "PATH": f"{d}:{os.environ['PATH']}",
            "SET_STATUS_READ_FIXTURE": str(fix),
            "SET_STATUS_GH_CAPTURE": str(capture),
        }
        if gh_fail:
            env["SET_STATUS_GH_FAIL"] = "1"
            env["SET_STATUS_GH_FAIL_MSG"] = gh_fail_msg
        proc = subprocess.run(
            ["bash", str(SCRIPT), *args],
            env=env, capture_output=True, text=True,
        )
        proc.mutation_argv = capture.read_text() if capture.exists() else ""
        return proc


# --- Exit-code contract (all three exit BEFORE any gh call) ------------------

def test_no_args_usage_exit_2():
    r = _run([])
    assert r.returncode == 2
    assert "Usage" in r.stderr


def test_one_arg_usage_exit_2():
    r = _run(["1024"])
    assert r.returncode == 2
    assert "Usage" in r.stderr


def test_non_numeric_issue_exit_2():
    r = _run(["abc", "Ready"])
    assert r.returncode == 2
    assert "numeric" in r.stderr.lower()


def test_unknown_status_exit_2():
    r = _run(["1024", "Bogus"])
    assert r.returncode == 2
    assert "unknown status" in r.stderr.lower()


# --- Status-name -> option-id map (the fat-finger guard) ---------------------

@pytest.mark.parametrize("status,option_id", sorted(CANONICAL.items()))
def test_option_id_map(status, option_id):
    # current status differs from target so a real move (mutation) is issued
    other = "Done" if status != "Done" else "Backlog"
    r = _run(["1024", status], envelope=_envelope(current=other))
    assert r.returncode == 0, r.stderr
    # the mutation must carry exactly the canonical option id for this name
    assert f"o={option_id}" in r.mutation_argv, r.mutation_argv


def test_map_covers_every_board_status():
    # the script must accept all seven names (none falls through to unknown)
    for status in CANONICAL:
        r = _run(["1024", status], envelope=_envelope(current="__none__"))
        assert r.returncode == 0, f"{status!r}: {r.stderr}"


# --- Idempotency: already-in-target is a no-op, no mutation ------------------

def test_already_in_target_is_noop_no_mutation():
    r = _run(["1024", "In review"], envelope=_envelope(current="In review"))
    assert r.returncode == 0, r.stderr
    assert "no-op" in r.stdout
    assert r.mutation_argv == ""  # no mutation issued


# --- Finding 1: query failure is distinct from an absent card ----------------

def test_no_card_on_board_exit_1():
    r = _run(["1024", "Ready"], envelope=_envelope(on_board=False))
    assert r.returncode == 1
    assert "no card" in r.stderr.lower()
    assert r.mutation_argv == ""


def test_query_failure_not_reported_as_no_card():
    r = _run(["1024", "Ready"], gh_fail=True, gh_fail_msg="gh: API rate limit exceeded")
    assert r.returncode == 1
    assert "board query failed" in r.stderr.lower()
    assert "no card" not in r.stderr.lower()  # must NOT misattribute
    assert "rate limit" in r.stderr           # surfaces the real error


def test_pr_number_gets_not_an_issue_hint():
    r = _run(["1039", "Ready"], gh_fail=True,
             gh_fail_msg="Could not resolve to an Issue with the number of 1039.")
    assert r.returncode == 1
    assert "not an issue" in r.stderr.lower()


# --- Resolved-target echo names repo + title (wrong-clone signal) ------------

def test_success_echo_names_repo_and_title():
    r = _run(["1024", "Ready"],
             envelope=_envelope(current="Backlog", title="my wrapper issue"))
    assert r.returncode == 0, r.stderr
    assert "TestOwner/TestRepo" in r.stdout
    assert "my wrapper issue" in r.stdout
