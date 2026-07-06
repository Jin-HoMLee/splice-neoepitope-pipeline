"""Unit tests for scripts/pm/gh_client.py (the shared hardened gh() wrapper, Issue #1017)."""

import subprocess
import sys
from pathlib import Path

import pytest

# Make the shared module importable as a module
SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts" / "pm"
sys.path.insert(0, str(SCRIPT_DIR))

import gh_client as ghc


class TestGhRetry:
    """gh() retries transient non-zero exits with backoff, honoring Retry-After
    (ported from recheck_milestone.py's #711 logic, Issue #1017).

    GitHub secondary-rate-limit 403s, transient 5xx, and replication-lag errors
    surface as a non-zero ``gh`` exit unrelated to the request. Deterministic 4xx
    (404/422/…) must NOT be retried, and a terminal failure raises so the
    check=True contract callers depend on is preserved.
    """

    @staticmethod
    def _result(returncode, stdout="", stderr=""):
        return subprocess.CompletedProcess(["gh"], returncode, stdout, stderr)

    def _runner_seq(self, results):
        """A fake subprocess.run yielding the given results in order; records calls."""
        seq = list(results)
        calls = []

        def run(cmd, **kwargs):
            calls.append(cmd)
            return seq[len(calls) - 1]

        run.calls = calls
        return run

    def test_retries_transient_then_succeeds(self):
        runner = self._runner_seq([
            self._result(1, stderr="HTTP 403: You have exceeded a secondary rate limit"),
            self._result(0, stdout='{"ok": true}'),
        ])
        sleeps = []
        out = ghc.gh("api", "x", _runner=runner, _sleep=sleeps.append)
        assert out == {"ok": True}
        assert len(runner.calls) == 2   # retried once
        assert len(sleeps) == 1         # backed off once

    def test_no_retry_on_deterministic_404(self):
        runner = self._runner_seq([
            self._result(1, stderr="gh: Not Found (HTTP 404)"),
        ])
        sleeps = []
        with pytest.raises(ghc.GhError):
            ghc.gh("api", "missing", _runner=runner, _sleep=sleeps.append)
        assert len(runner.calls) == 1   # NOT retried
        assert sleeps == []

    def test_honors_retry_after_header(self):
        runner = self._runner_seq([
            self._result(1, stderr="HTTP 403: rate limited\nRetry-After: 7"),
            self._result(0, stdout="[]"),
        ])
        sleeps = []
        ghc.gh("api", "x", _runner=runner, _sleep=sleeps.append)
        assert sleeps == [7.0]          # used the Retry-After hint, not exp backoff

    def test_caps_excessive_retry_after(self):
        # GitHub can legally emit a large Retry-After during sustained degradation;
        # an uncapped hint would stall the nightly job. Cap it (Issue #711 review).
        runner = self._runner_seq([
            self._result(1, stderr="HTTP 403: rate limited\nRetry-After: 3600"),
            self._result(0, stdout="[]"),
        ])
        sleeps = []
        ghc.gh("api", "x", _runner=runner, _sleep=sleeps.append)
        assert sleeps == [ghc.GH_RETRY_AFTER_CAP_SECONDS]   # 3600 clamped to the ceiling

    def test_exhausts_attempts_then_raises(self):
        runner = self._runner_seq(
            [self._result(1, stderr="HTTP 503: server error")] * ghc.GH_MAX_ATTEMPTS
        )
        sleeps = []
        with pytest.raises(ghc.GhError):
            ghc.gh("api", "x", _runner=runner, _sleep=sleeps.append)
        assert len(runner.calls) == ghc.GH_MAX_ATTEMPTS
        assert len(sleeps) == ghc.GH_MAX_ATTEMPTS - 1  # no sleep after the last attempt

    def test_success_first_try_no_sleep(self):
        runner = self._runner_seq([self._result(0, stdout='{"a": 1}')])
        sleeps = []
        out = ghc.gh("api", "x", _runner=runner, _sleep=sleeps.append)
        assert out == {"a": 1}
        assert len(runner.calls) == 1
        assert sleeps == []

    def test_parse_json_false_returns_raw_stdout(self):
        runner = self._runner_seq([self._result(0, stdout="raw text")])
        out = ghc.gh("api", "x", parse_json=False, _runner=runner, _sleep=lambda _: None)
        assert out == "raw text"


class TestGhError:
    """The typed hard-failure exception (Issue #1017 AC #1).

    ``GhError`` subclasses ``subprocess.CalledProcessError`` so a caller can catch
    the specific type to isolate a per-item failure (the #989/#1012 shape) while
    every pre-existing ``except subprocess.CalledProcessError`` site keeps catching
    it unchanged (backward compatible).
    """

    @staticmethod
    def _result(returncode, stdout="", stderr=""):
        return subprocess.CompletedProcess(["gh"], returncode, stdout, stderr)

    def _runner_once(self, result):
        def run(cmd, **kwargs):
            return result
        return run

    def test_gherror_is_calledprocesserror_subclass(self):
        assert issubclass(ghc.GhError, subprocess.CalledProcessError)

    def test_terminal_failure_raises_gherror_caught_as_calledprocesserror(self):
        runner = self._runner_once(self._result(1, stderr="gh: Not Found (HTTP 404)"))
        # The specific type is GhError...
        with pytest.raises(ghc.GhError):
            ghc.gh("api", "missing", _runner=runner, _sleep=lambda _: None)
        # ...and it is still caught by the legacy contract.
        with pytest.raises(subprocess.CalledProcessError):
            ghc.gh("api", "missing", _runner=runner, _sleep=lambda _: None)

    def test_gherror_preserves_returncode_cmd_and_stderr(self):
        runner = self._runner_once(self._result(422, stdout="out", stderr="HTTP 422: bad"))
        with pytest.raises(ghc.GhError) as exc:
            ghc.gh("api", "bad", _runner=runner, _sleep=lambda _: None)
        assert exc.value.returncode == 422
        assert exc.value.cmd == ["gh", "api", "bad"]
        assert exc.value.stderr == "HTTP 422: bad"
