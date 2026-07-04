"""Tests for the awaiting-bot-review skill's poller (Issue #864).

Two layers, both hermetic (no live gh / no network):
  * the detection predicate `match_review.jq` - piped fixtures through real `jq`,
    which is the regression guard for the exact `gh --jq --arg` bug class that made
    an earlier inline poller silently time out (Issue #864 "why a tested script").
  * `poll_bot_review.sh` arg parsing + one-shot detection - driven with a fake `gh`
    on PATH so the real script runs end-to-end deterministically.
"""

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

SKILL = Path(__file__).parent.parent.parent / ".agents" / "skills" / "awaiting-bot-review"
MATCH_JQ = SKILL / "scripts" / "match_review.jq"
POLL_SH = SKILL / "scripts" / "poll_bot_review.sh"

requires_jq = pytest.mark.skipif(shutil.which("jq") is None, reason="jq not on PATH")

WATERMARK = "2026-07-04T01:49:00Z"


def _finished(created="2026-07-04T01:49:26Z", login="claude"):
    return {
        "author": {"login": login},
        "createdAt": created,
        "body": "**Claude finished @Jin-HoMLee's task in 5m 18s**\n\n### PR Review\nLGTM.",
    }


def _working(created="2026-07-04T01:49:20Z", login="claude"):
    return {
        "author": {"login": login},
        "createdAt": created,
        "body": "**Claude is working...** [View job]",
    }


def _comments(*nodes):
    return json.dumps({"comments": list(nodes)})


# --------------------------------------------------------------------------
# Layer 1: the jq detection predicate
# --------------------------------------------------------------------------
@requires_jq
class TestMatchReviewPredicate:
    def _run(self, comments_json, watermark=WATERMARK):
        """Pipe fixture JSON through the real predicate, exactly as the script does."""
        proc = subprocess.run(
            ["jq", "-r", "--arg", "w", watermark, "-f", str(MATCH_JQ)],
            input=comments_json, capture_output=True, text=True, timeout=15,
        )
        assert proc.returncode == 0, proc.stderr
        return proc.stdout.strip()

    def test_matches_finished_review_after_watermark(self):
        out = self._run(_comments(_finished()))
        assert "Claude finished" in out           # the verdict body is surfaced

    def test_ignores_interim_working_state(self):
        # A "working" comment is newer than the watermark and authored by the bot,
        # but lacks "Claude finished", so it must NOT be treated as landed.
        assert self._run(_comments(_working())) == ""

    def test_ignores_review_at_or_before_watermark(self):
        # A prior review on the same PR (createdAt <= watermark) must not re-match.
        old = _finished(created="2026-07-04T01:48:00Z")
        assert self._run(_comments(old)) == ""

    def test_ignores_non_bot_author(self):
        # A human quoting "Claude finished" must not trip detection.
        human = _finished(login="Jin-HoMLee")
        assert self._run(_comments(human)) == ""

    def test_picks_newest_when_multiple_finished(self):
        older = _finished(created="2026-07-04T01:49:26Z")
        older["body"] = "**Claude finished** older verdict"
        newer = _finished(created="2026-07-04T02:10:00Z")
        newer["body"] = "**Claude finished** newer verdict"
        out = self._run(_comments(older, newer))
        assert "newer verdict" in out and "older verdict" not in out

    def test_empty_comments_is_no_match(self):
        assert self._run(_comments()) == ""


# --------------------------------------------------------------------------
# Layer 2: the poll_bot_review.sh wrapper (arg parsing + one-shot detection)
# --------------------------------------------------------------------------
_FAKE_GH = (
    "#!/usr/bin/env python3\n"
    "import os, sys\n"
    "a = sys.argv[1:]\n"
    "if a[:2] == ['pr', 'view']:\n"
    "    sys.stdout.write(os.environ.get('FAKE_COMMENTS', '{\"comments\":[]}'))\n"
    "    sys.exit(0)\n"
    "sys.stderr.write('fake gh unmatched: %r\\n' % (a,)); sys.exit(1)\n"
)


# Fake gh that fails its first `pr view` call then serves FAKE_COMMENTS: models a
# transient blip (network / rate-limit / 5xx) that the poller must survive.
_FAKE_GH_FAIL_ONCE = (
    "#!/usr/bin/env python3\n"
    "import os, sys\n"
    "a = sys.argv[1:]\n"
    "if a[:2] != ['pr', 'view']:\n"
    "    sys.stderr.write('unmatched\\n'); sys.exit(1)\n"
    "c = os.environ['FAKE_COUNTER']\n"
    "try:\n"
    "    n = int(open(c).read())\n"
    "except Exception:\n"
    "    n = 0\n"
    "open(c, 'w').write(str(n + 1))\n"
    "if n == 0:\n"
    "    sys.stderr.write('gh: simulated transient failure\\n'); sys.exit(1)\n"
    "sys.stdout.write(os.environ.get('FAKE_COMMENTS', '{\"comments\":[]}')); sys.exit(0)\n"
)

# Fake gh that always fails: models a persistent error (bad auth / deleted PR).
_FAKE_GH_ALWAYS_FAIL = (
    "#!/usr/bin/env python3\n"
    "import sys\n"
    "sys.stderr.write('gh: simulated persistent failure\\n'); sys.exit(1)\n"
)


@requires_jq
class TestPollScript:
    def _run(self, args, comments_json="", tmp_path=None, gh_body=_FAKE_GH, extra_env=None):
        gh = tmp_path / "gh"
        gh.write_text(gh_body)
        gh.chmod(0o755)
        env = dict(os.environ)
        env["PATH"] = f"{tmp_path}:{env['PATH']}"
        if comments_json:
            env["FAKE_COMMENTS"] = comments_json
        if extra_env:
            env.update(extra_env)
        return subprocess.run(
            ["bash", str(POLL_SH), *args],
            capture_output=True, text=True, timeout=30, env=env,
        )

    def test_detects_landed_review_exit_zero(self, tmp_path):
        r = self._run(
            [str(985), "--since", WATERMARK, "--timeout-min", "0", "--interval-sec", "0"],
            comments_json=_comments(_finished()), tmp_path=tmp_path,
        )
        assert r.returncode == 0, r.stderr
        assert "Bot review landed on PR #985" in r.stdout
        assert "Claude finished" in r.stdout

    def test_times_out_on_working_only(self, tmp_path):
        # Only an interim "working" comment present, so no landed review: exit 3
        # after one poll (timeout 0), notify-and-stop.
        r = self._run(
            [str(985), "--since", WATERMARK, "--timeout-min", "0", "--interval-sec", "0"],
            comments_json=_comments(_working()), tmp_path=tmp_path,
        )
        assert r.returncode == 3
        assert "no bot review" in r.stderr

    def test_missing_pr_arg_is_usage_error(self, tmp_path):
        r = self._run(["--since", WATERMARK], tmp_path=tmp_path)
        assert r.returncode == 2
        assert "missing required" in r.stderr

    def test_non_numeric_pr_is_usage_error(self, tmp_path):
        r = self._run(["not-a-number"], tmp_path=tmp_path)
        assert r.returncode == 2
        assert "must be a number" in r.stderr

    def test_bad_timeout_is_usage_error(self, tmp_path):
        r = self._run([str(985), "--timeout-min", "soon"], tmp_path=tmp_path)
        assert r.returncode == 2
        assert "timeout-min" in r.stderr

    def test_unknown_option_is_usage_error(self, tmp_path):
        r = self._run([str(985), "--merge-when-done"], tmp_path=tmp_path)
        assert r.returncode == 2
        assert "unknown option" in r.stderr

    def test_help_renders_usage_and_exits_zero(self, tmp_path):
        # The help printer parses the script's own header; a non-portable sed form
        # broke it once, so pin exit 0 + real usage text (no shebang leakage).
        r = self._run(["--help"], tmp_path=tmp_path)
        assert r.returncode == 0, r.stderr
        assert "Usage:" in r.stdout
        assert not r.stdout.lstrip().startswith("!")   # shebang must not leak

    def test_survives_transient_gh_failure_then_detects(self, tmp_path):
        # PR #993 review, finding 1: a single transient `gh` failure must NOT abort
        # the poller (it would under `set -euo pipefail` without the guard). First
        # poll fails, second succeeds with a landed review -> exit 0.
        r = self._run(
            [str(985), "--since", WATERMARK, "--timeout-min", "5", "--interval-sec", "0"],
            comments_json=_comments(_finished()),
            gh_body=_FAKE_GH_FAIL_ONCE,
            extra_env={"FAKE_COUNTER": str(tmp_path / "count")},
            tmp_path=tmp_path,
        )
        assert r.returncode == 0, r.stderr
        assert "Claude finished" in r.stdout

    def test_gives_up_after_consecutive_failures(self, tmp_path):
        # A persistent failure surfaces loudly (exit 4) instead of masquerading as a
        # clean 25-min timeout. timeout-min large so the failure breaker, not the
        # deadline, is what stops it.
        r = self._run(
            [str(985), "--since", WATERMARK, "--timeout-min", "5", "--interval-sec", "0"],
            gh_body=_FAKE_GH_ALWAYS_FAIL, tmp_path=tmp_path,
        )
        assert r.returncode == 4
        assert "consecutive" in r.stderr

    def test_option_value_swallow_is_guarded(self, tmp_path):
        # PR #993 review, nit 1: `--since --timeout-min 5` must error, not take
        # "--timeout-min" as the watermark.
        r = self._run(["--since", "--timeout-min", "5", str(985)], tmp_path=tmp_path)
        assert r.returncode == 2
        assert "requires a value" in r.stderr
