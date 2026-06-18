"""Tests for the non-blocking stray-AC-box lint gate (Issue #730).

Covers the single-sourced orchestration (closure_audit.collect_stray_ac_warnings)
and the thin CLI wrapper's always-exit-0 (non-blocking) contract + fail-open.
"""

import json
import subprocess

import ac_section_lint as gate
import closure_audit as ca


def _fake_io(monkeypatch, pr, issues_by_num):
    monkeypatch.setattr(ca, "fetch_pr", lambda n, repo=None: pr)
    monkeypatch.setattr(ca, "fetch_issue", lambda n, repo=None: issues_by_num[n])


def _pr(refs=(42,)):
    return {"closingIssuesReferences": [{"number": n} for n in refs]}


def _issue(n, body):
    return {"number": n, "body": body}


# --- collect_stray_ac_warnings: the single-sourced check ---


def test_warns_when_issue_has_stray_boxes_without_ac(monkeypatch):
    body = "## Plan (phased)\n- [ ] P1\n- [ ] P2\n"
    _fake_io(monkeypatch, _pr(refs=(42,)), {42: _issue(42, body)})
    warnings = ca.collect_stray_ac_warnings(99)
    assert len(warnings) == 1
    num, msg = warnings[0]
    assert num == 42
    assert "Plan (phased)" in msg


def test_clean_when_issue_has_ac_section(monkeypatch):
    body = "## Plan\n- [ ] stray\n\n## Acceptance criteria\n- [ ] real\n"
    _fake_io(monkeypatch, _pr(refs=(42,)), {42: _issue(42, body)})
    assert ca.collect_stray_ac_warnings(99) == []


def test_clean_when_no_linked_issues(monkeypatch):
    _fake_io(monkeypatch, _pr(refs=()), {})
    assert ca.collect_stray_ac_warnings(99) == []


def test_collects_across_multiple_issues(monkeypatch):
    issues = {
        42: _issue(42, "## Tasks\n- [ ] a\n"),                       # stray → warn
        43: _issue(43, "## Acceptance criteria\n- [ ] ok\n"),        # clean
        44: _issue(44, "## Checklist\n- [ ] b\n- [ ] c\n"),          # stray → warn
    }
    _fake_io(monkeypatch, _pr(refs=(42, 43, 44)), issues)
    warnings = ca.collect_stray_ac_warnings(99)
    assert sorted(n for n, _ in warnings) == [42, 44]


# --- ac_section_lint.main: CLI is non-blocking (always exit 0) ---


def _run(monkeypatch, argv):
    monkeypatch.setattr(gate.sys, "argv", argv)
    return gate.main()


def test_cli_clean_returns_0(monkeypatch):
    monkeypatch.setattr(ca, "collect_stray_ac_warnings", lambda n, repo=None: [])
    assert _run(monkeypatch, ["ac_section_lint.py", "99"]) == 0


def test_cli_warning_is_non_blocking_returns_0(monkeypatch, capsys):
    monkeypatch.setattr(
        ca, "collect_stray_ac_warnings",
        lambda n, repo=None: [(42, "2 unticked checkbox(es) under 'Plan'")],
    )
    rc = _run(monkeypatch, ["ac_section_lint.py", "99"])
    assert rc == 0  # lint warns but never blocks the merge
    err = capsys.readouterr().err
    assert "Issue #42" in err
    assert "Plan" in err


def test_cli_gh_error_fails_open(monkeypatch):
    def boom(n, repo=None):
        raise subprocess.CalledProcessError(1, ["gh"])

    monkeypatch.setattr(ca, "collect_stray_ac_warnings", boom)
    assert _run(monkeypatch, ["ac_section_lint.py", "99"]) == 0


def test_cli_json_error_fails_open(monkeypatch):
    def boom(n, repo=None):
        raise json.JSONDecodeError("bad", "doc", 0)

    monkeypatch.setattr(ca, "collect_stray_ac_warnings", boom)
    assert _run(monkeypatch, ["ac_section_lint.py", "99"]) == 0


def test_cli_bad_usage_returns_2(monkeypatch):
    assert _run(monkeypatch, ["ac_section_lint.py"]) == 2
    assert _run(monkeypatch, ["ac_section_lint.py", "not-a-number"]) == 2


# --- #607 parity: REPO env read + forwarded ---


def test_cli_reads_repo_env_and_forwards(monkeypatch):
    captured = {}

    def fake(n, repo=None):
        captured["repo"] = repo
        return []

    monkeypatch.setenv("REPO", "fork/repo")
    monkeypatch.setattr(ca, "collect_stray_ac_warnings", fake)
    assert _run(monkeypatch, ["ac_section_lint.py", "99"]) == 0
    assert captured["repo"] == "fork/repo"


def test_cli_defaults_repo_when_env_unset(monkeypatch):
    captured = {}

    def fake(n, repo=None):
        captured["repo"] = repo
        return []

    monkeypatch.delenv("REPO", raising=False)
    monkeypatch.setattr(ca, "collect_stray_ac_warnings", fake)
    assert _run(monkeypatch, ["ac_section_lint.py", "99"]) == 0
    assert captured["repo"] == "Jin-HoMLee/splice-neoepitope-pipeline"
