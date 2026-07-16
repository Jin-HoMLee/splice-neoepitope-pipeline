"""Tests for `.agents/hooks/check_gh_issue_create_repo.py` (Issue #1083).

Pure-function unit tests import the module directly (no network). Orchestration
tests monkeypatch `resolve_target_repo` so `main()`'s ask/allow decision runs
without a live `gh` call. Subprocess tests pipe real PreToolUse JSON to the hook
(the "drive the real trigger, not a synthetic payload" discipline) and only
travel paths that return before any `gh` call.

Matched-pair controls throughout: the SAME create asks when target and shape
disagree AND allows when they agree - a check that could only fire is not a
check.
"""

import io
import json
import os
import subprocess
import sys
from pathlib import Path

HOOKS_DIR = Path(__file__).parent.parent.parent / ".agents" / "hooks"
HOOK = HOOKS_DIR / "check_gh_issue_create_repo.py"
sys.path.insert(0, str(HOOKS_DIR))

import check_gh_issue_create_repo as h  # noqa: E402

PROJECT = h.PROJECT_REPO
PERSONAS = h.PERSONAS_REPO


# --- the hook must actually be executable (feedback_hook_exec_bit) ---


def test_hook_is_executable():
    assert os.access(HOOK, os.X_OK), "hook must be committed 100755 or it never runs"


# --- create_args: command matching ---


class TestCreateArgs:
    def test_plain_invocation(self):
        assert h.create_args('gh issue create --title x') == ["--title", "x"]

    def test_after_separator(self):
        assert h.create_args('git add -A && gh issue create -t x') == ["-t", "x"]

    def test_trailing_separator_stops_collection(self):
        assert h.create_args('gh issue create -t x && echo done') == ["-t", "x"]

    def test_issue_list_not_matched(self):
        assert h.create_args("gh issue list") is None

    def test_issue_view_not_matched(self):
        assert h.create_args("gh issue view 549") is None

    def test_create_inside_quoted_body_not_matched(self):
        # a comment documenting the command must NOT match
        cmd = 'gh pr comment 1 --body "run gh issue create --repo foo/bar"'
        assert h.create_args(cmd) is None

    def test_unbalanced_quotes_fail_safe(self):
        assert h.create_args('gh issue create --title "oops') is None

    def test_heredoc_body_discussing_command_not_matched(self):
        cmd = (
            "cat > /tmp/c.md <<'EOF'\n"
            "example: gh issue create --repo x/y --title memory\n"
            "EOF\n"
            "gh pr comment 5 --body-file /tmp/c.md"
        )
        assert h.create_args(cmd) is None


# --- parse_fields: title / body / labels ---


class TestParseFields:
    def test_title_body_labels(self):
        t, b, ls = h.parse_fields(
            ["--title", "T", "--body", "B", "--label", "role:developer"]
        )
        assert t == "T" and b == "B" and ls == ["role:developer"]

    def test_equals_forms(self):
        t, b, ls = h.parse_fields(["--title=T", "--body=B", "--label=a,b"])
        assert t == "T" and b == "B" and ls == ["a", "b"]

    def test_repeated_and_comma_labels(self):
        _, _, ls = h.parse_fields(["--label", "a,b", "--label", "c"])
        assert ls == ["a", "b", "c"]

    def test_short_flags(self):
        t, b, ls = h.parse_fields(["-t", "T", "-b", "B", "-l", "x"])
        assert t == "T" and b == "B" and ls == ["x"]


# --- repo_signal: shape detection ---


class TestRepoSignal:
    def test_mm_label_is_deterministic_spine(self):
        side, why = h.repo_signal("anything", "", ["role:memory_manager"])
        assert side == "personas" and "memory_manager" in why

    def test_one_sided_personas_keywords(self):
        side, _ = h.repo_signal("episodic memory drain for the persona", "", [])
        assert side == "personas"

    def test_one_sided_project_keywords(self):
        side, _ = h.repo_signal("snakemake alignment junction rule", "", [])
        assert side == "project"

    def test_mixed_signal_is_none(self):
        # both sides present -> ambiguous -> allow
        side, _ = h.repo_signal("snakemake pipeline memory episode", "", [])
        assert side is None

    def test_single_hit_below_threshold_is_none(self):
        side, _ = h.repo_signal("a note about episode handling", "", [])
        assert side is None

    def test_empty_is_none(self):
        assert h.repo_signal("", "", []) == (None, "")

    def test_star_substring_false_positive_fixed(self):
        # #1210 review: "start ... workflow" must NOT read as project - `star` is
        # a substring of start/restart/startup. Word-boundary matching kills it:
        # only "workflow" hits (1 < 2 threshold) -> ambiguous -> allow.
        side, _ = h.repo_signal("start the board-hygiene workflow", "", [])
        assert side is None

    def test_star_whole_word_still_matches(self):
        # the fix must not lose the real signal: `star` as a whole word still hits
        side, _ = h.repo_signal("star aligner junction fix", "", [])
        assert side == "project"

    def test_fragment_keyword_still_substring(self):
        # FRAG keywords keep substring semantics (feedback_ / shared/)
        side, _ = h.repo_signal("update feedback_foo and shared/bar", "", [])
        assert side == "personas"


# --- ask_payload: shape ---


def test_ask_payload_shape():
    out = h.ask_payload(PROJECT, PERSONAS, "role:memory_manager label")["hookSpecificOutput"]
    assert out["hookEventName"] == "PreToolUse"
    assert out["permissionDecision"] == "ask"  # NOT deny
    assert PROJECT in out["permissionDecisionReason"]
    assert PERSONAS in out["permissionDecisionReason"]
    assert "CLAUDE_ALLOW_REPO_MISFILE" in out["permissionDecisionReason"]


# --- main(): orchestration (monkeypatched gh I/O, no network) ---


def _run_main(monkeypatch, capsys, cmd, *, target):
    monkeypatch.setattr(h, "resolve_target_repo", lambda args: target)
    payload = json.dumps({"tool_input": {"command": cmd}})
    monkeypatch.setattr(sys, "stdin", io.StringIO(payload))
    monkeypatch.delenv("CLAUDE_ALLOW_REPO_MISFILE", raising=False)
    rc = h.main()
    captured = capsys.readouterr().out.strip()
    return rc, (json.loads(captured) if captured else None)


# Matched pair 1: same MM-shaped create, target flips the outcome.
def test_asks_mm_shape_into_project(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(h, "LOG_PATH", tmp_path / "hook_fires.jsonl")
    cmd = 'gh issue create --title t --label role:memory_manager --repo Jin-HoMLee/' + PROJECT
    rc, out = _run_main(monkeypatch, capsys, cmd, target=("Jin-HoMLee", PROJECT))
    assert rc == 0
    assert out["hookSpecificOutput"]["permissionDecision"] == "ask"
    rec = json.loads((tmp_path / "hook_fires.jsonl").read_text().strip())
    assert rec["hook"] == "check_gh_issue_create_repo" and rec["target"] == PROJECT


def test_allows_mm_shape_into_personas(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(h, "LOG_PATH", tmp_path / "hook_fires.jsonl")
    cmd = 'gh issue create --title t --label role:memory_manager --repo Jin-HoMLee/' + PERSONAS
    rc, out = _run_main(monkeypatch, capsys, cmd, target=("Jin-HoMLee", PERSONAS))
    assert rc == 0 and out is None
    assert not (tmp_path / "hook_fires.jsonl").exists()  # no fire on allow


# Matched pair 2: pipeline shape into personas asks; into project allows.
def test_asks_project_shape_into_personas(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(h, "LOG_PATH", tmp_path / "hook_fires.jsonl")
    cmd = 'gh issue create -t "snakemake alignment junction fix" --repo Jin-HoMLee/' + PERSONAS
    rc, out = _run_main(monkeypatch, capsys, cmd, target=("Jin-HoMLee", PERSONAS))
    assert out["hookSpecificOutput"]["permissionDecision"] == "ask"


def test_allows_project_shape_into_project(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(h, "LOG_PATH", tmp_path / "hook_fires.jsonl")
    cmd = 'gh issue create -t "snakemake alignment junction fix" --repo Jin-HoMLee/' + PROJECT
    rc, out = _run_main(monkeypatch, capsys, cmd, target=("Jin-HoMLee", PROJECT))
    assert rc == 0 and out is None


def test_ambiguous_shape_fails_open(monkeypatch, capsys):
    cmd = 'gh issue create -t "small fix" --repo Jin-HoMLee/' + PROJECT
    rc, out = _run_main(monkeypatch, capsys, cmd, target=("Jin-HoMLee", PROJECT))
    assert rc == 0 and out is None


def test_start_workflow_into_personas_allows(monkeypatch, capsys):
    # #1210 review scenario end-to-end: an unlabelled governance-ish create must
    # NOT false-ask when it lands in personas ("start" no longer reads as "star").
    cmd = 'gh issue create -t "start the board-hygiene workflow" --repo Jin-HoMLee/' + PERSONAS
    rc, out = _run_main(monkeypatch, capsys, cmd, target=("Jin-HoMLee", PERSONAS))
    assert rc == 0 and out is None


def test_unresolvable_target_fails_open(monkeypatch, capsys):
    cmd = 'gh issue create --title t --label role:memory_manager'
    rc, out = _run_main(monkeypatch, capsys, cmd, target=None)
    assert rc == 0 and out is None


def test_third_repo_out_of_scope(monkeypatch, capsys):
    cmd = 'gh issue create --title t --label role:memory_manager --repo other/thing'
    rc, out = _run_main(monkeypatch, capsys, cmd, target=("other", "thing"))
    assert rc == 0 and out is None


def test_escape_hatch(monkeypatch, capsys):
    monkeypatch.setenv("CLAUDE_ALLOW_REPO_MISFILE", "1")
    monkeypatch.setattr(h, "resolve_target_repo", lambda args: ("Jin-HoMLee", PROJECT))
    cmd = 'gh issue create --title t --label role:memory_manager --repo Jin-HoMLee/' + PROJECT
    monkeypatch.setattr(sys, "stdin", io.StringIO(json.dumps({"tool_input": {"command": cmd}})))
    rc = h.main()
    assert rc == 0 and capsys.readouterr().out.strip() == ""


# --- subprocess fail-open / no-op (real trigger, no live gh reached) ---


def _run_subprocess(stdin_payload: str):
    proc = subprocess.run(
        [sys.executable, str(HOOK)],
        input=stdin_payload, capture_output=True, text=True, timeout=10,
    )
    return proc.returncode, proc.stdout, proc.stderr


def _payload(command: str) -> str:
    return json.dumps({"tool_input": {"command": command}})


class TestSubprocessNoOp:
    def test_empty_stdin_fail_open(self):
        rc, out, _ = _run_subprocess("")
        assert rc == 0 and out.strip() == ""

    def test_malformed_json_fail_open(self):
        rc, out, _ = _run_subprocess("{not valid json")
        assert rc == 0 and out.strip() == ""

    def test_non_create_command_noop(self):
        rc, out, _ = _run_subprocess(_payload("gh issue view 549"))
        assert rc == 0 and out.strip() == ""

    def test_ambiguous_shape_noop_before_gh(self):
        # no shape signal -> returns before any resolve_target_repo gh call
        rc, out, _ = _run_subprocess(_payload('gh issue create -t "tiny tweak"'))
        assert rc == 0 and out.strip() == ""
