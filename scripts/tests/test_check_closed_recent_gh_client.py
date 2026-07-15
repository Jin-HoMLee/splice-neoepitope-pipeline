"""Issue #1062 - check_closed_recent's gh I/O now goes through the shared gh_client.

The local `gh_json()` was a bare `subprocess.run(..., check=True)` with **no retry**:
a transient 5xx / secondary-rate-limit / replication blip raised. In a *recap*, that
surfaces as "nothing closed recently" - which is precisely the silent-wrong-answer
this script exists to prevent (Issue #784). It was the last hand-rolled `gh` helper;
#1017 and #1055 converged the other four.

These tests prove the convergence is BEHAVIORAL, not cosmetic: a transient failure
must now be retried and survive. The matched-pair control is the deterministic 4xx,
which must still raise - otherwise "it retries" would be indistinguishable from
"it swallows every error".
"""
import functools
import json
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPTS_PM = Path(__file__).parent.parent / "pm"
sys.path.insert(0, str(SCRIPTS_PM))

import check_closed_recent as c  # noqa: E402
import gh_client  # noqa: E402

ROWS = [{"number": 1, "title": "t", "closedAt": "2026-07-14T00:00:00Z"}]


class _Result:
    def __init__(self, returncode, stdout="", stderr=""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


def _fake_runner(script):
    """A subprocess.run stand-in replaying `script` (a list of _Result), recording calls."""
    calls = []

    def run(cmd, **kwargs):
        calls.append(cmd)
        return script[min(len(calls) - 1, len(script) - 1)]

    run.calls = calls
    return run


def _inject(monkeypatch, runner):
    """Point the module's `gh` at the REAL gh_client.gh with the runner injected.

    Note gh()'s `_runner=subprocess.run` default is bound at import time, so
    monkeypatching `gh_client.subprocess.run` would patch nothing - the default has
    already captured the original callable. Injecting through the seam the wrapper
    exposes keeps the real retry/backoff logic under test rather than a stub of it.
    """
    monkeypatch.setattr(
        c, "gh", functools.partial(gh_client.gh, _runner=runner, _sleep=lambda _s: None)
    )


def test_transient_failure_is_retried_and_survives(monkeypatch):
    """The whole point of #1062: a blip no longer becomes a false 'nothing closed'."""
    runner = _fake_runner([
        _Result(1, "", "HTTP 502 Bad Gateway"),      # transient: must be retried
        _Result(0, json.dumps(ROWS), ""),            # ...and then succeed
    ])
    _inject(monkeypatch, runner)

    rows = c.list_closed_issues("closed:>=2026-07-13")

    assert rows == ROWS
    assert len(runner.calls) == 2, "a transient failure must be RETRIED, not raised"


def test_deterministic_4xx_still_raises(monkeypatch):
    """MATCHED-PAIR CONTROL. Without this, 'it retries' could mean 'it swallows
    everything' - and a recap that silently eats a real error is the original bug."""
    runner = _fake_runner([_Result(1, "", "HTTP 404 Not Found")])
    _inject(monkeypatch, runner)

    with pytest.raises(subprocess.CalledProcessError):
        c.list_closed_issues("closed:>=2026-07-13")

    assert len(runner.calls) == 1, "a deterministic 4xx must NOT be retried"


def test_merged_pr_listing_also_goes_through_the_wrapper(monkeypatch):
    """Both call sites converged, not just the first one."""
    prs = [
        {"number": 9, "title": "merged", "mergedAt": "2026-07-14T00:00:00Z"},
        {"number": 10, "title": "closed unmerged", "mergedAt": None},
    ]
    runner = _fake_runner([
        _Result(1, "", "HTTP 503 Service Unavailable"),
        _Result(0, json.dumps(prs), ""),
    ])
    _inject(monkeypatch, runner)

    rows = c.list_merged_prs("closed:>=2026-07-13")

    assert [r["number"] for r in rows] == [9], "unmerged PRs are still dropped"
    assert len(runner.calls) == 2, "list_merged_prs must retry too"


def test_no_hand_rolled_helper_survives():
    """`gh_json` is gone; the module now imports the shared wrapper."""
    assert not hasattr(c, "gh_json")
    assert c.gh is gh_client.gh
