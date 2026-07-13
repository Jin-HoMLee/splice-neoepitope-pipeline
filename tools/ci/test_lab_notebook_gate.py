"""Tests for the pre-merge lab-notebook gate (Issue #409).

Covers both the single-sourced orchestration (closure_audit.audit_pr_pre_merge)
and the thin CLI wrapper's exit codes (lab_notebook_gate.main).
"""

import json
import subprocess

import closure_audit as ca
import lab_notebook_gate as gate

DATE = "2026-06-01"


def _fake_io(monkeypatch, pr, issues_by_num, notebooks):
    monkeypatch.setattr(ca, "fetch_pr", lambda n, repo=None: pr)
    monkeypatch.setattr(ca, "fetch_issue", lambda n, repo=None: issues_by_num[n])
    monkeypatch.setattr(ca, "_load_notebook", lambda role: notebooks.get(role))


def _pr(refs=(42,), files=("workflow/x.py",), body=""):
    return {
        "closingIssuesReferences": [{"number": n} for n in refs],
        "files": [{"path": p} for p in files],
        "body": body,
    }


def _issue(n=42, roles=("scientist",)):
    return {"number": n, "labels": [{"name": f"role:{r}"} for r in roles]}


# --- audit_pr_pre_merge: the single-sourced check ---


def test_missing_entry_returns_gap(monkeypatch):
    # Entry exists for the date but references neither the PR (#99) nor Issue #42.
    nb = "## 2026-06-01\n### 10:00 UTC — Editor: X\nunrelated work"
    _fake_io(monkeypatch, _pr(), {42: _issue()}, {"scientist": nb})
    assert ca.audit_pr_pre_merge(99, DATE)


def test_present_entry_by_issue_number_returns_clean(monkeypatch):
    nb = "## 2026-06-01\n### 10:00 UTC — Editor: X\nclosing #42"
    _fake_io(monkeypatch, _pr(), {42: _issue()}, {"scientist": nb})
    assert ca.audit_pr_pre_merge(99, DATE) == []


def test_present_entry_by_pr_number_returns_clean(monkeypatch):
    nb = "## 2026-06-01\n### 10:00 UTC — Editor: X\nlanded in #99"
    _fake_io(monkeypatch, _pr(), {42: _issue(roles=("developer",))}, {"developer": nb})
    assert ca.audit_pr_pre_merge(99, DATE) == []


def test_entry_without_subsection_returns_gap(monkeypatch):
    nb = "## 2026-06-01\nclosing #99 but no '### HH:MM UTC' sub-header"
    _fake_io(monkeypatch, _pr(), {42: _issue(roles=("developer",))}, {"developer": nb})
    assert ca.audit_pr_pre_merge(99, DATE)


def test_entry_without_reference_returns_gap(monkeypatch):
    nb = "## 2026-06-01\n### 10:00 UTC — Editor: X\nsome work, no issue ref"
    _fake_io(monkeypatch, _pr(), {42: _issue(roles=("developer",))}, {"developer": nb})
    assert ca.audit_pr_pre_merge(99, DATE)


def test_missing_notebook_file_returns_gap(monkeypatch):
    _fake_io(monkeypatch, _pr(), {42: _issue()}, {"scientist": None})
    assert ca.audit_pr_pre_merge(99, DATE)


def test_pure_memory_manager_pr_returns_clean(monkeypatch):
    # #748: a pure role:memory_manager Issue needs no project-repo notebook entry
    # (MM's record is the personas-repo git log). End-to-end through the gate.
    _fake_io(monkeypatch, _pr(), {42: _issue(roles=("memory_manager",))}, {})
    assert ca.audit_pr_pre_merge(99, DATE) == []


def test_mixed_memory_manager_pr_still_requires_other_role(monkeypatch):
    # #748: dev+MM Issue with no developer entry still gaps (strip-role, not skip-PR).
    _fake_io(
        monkeypatch,
        _pr(),
        {42: _issue(roles=("developer", "memory_manager"))},
        {"developer": None},
    )
    assert ca.audit_pr_pre_merge(99, DATE)


def test_file_exempt_returns_clean(monkeypatch):
    # Diff is entirely the lab-notebook itself → no entry required for itself.
    pr = _pr(files=("research/lab_notebook/scientist.md",))
    _fake_io(monkeypatch, pr, {42: _issue()}, {"scientist": None})
    assert ca.audit_pr_pre_merge(99, DATE) == []


def test_skip_marker_returns_clean(monkeypatch):
    # The marker must be on its OWN line (Issue #1126). This fixture previously
    # trailed it after prose on the same line; that form no longer opts out, and
    # the change is deliberate - an unanchored match let a PR body that merely
    # *documents* the marker skip this gate silently. The failure direction is
    # safe: a mis-formatted marker now blocks the merge loudly rather than
    # bypassing the check quietly.
    pr = _pr(body="routine ship\n\n<!-- skip-lab-notebook: routine -->\n")
    _fake_io(monkeypatch, pr, {42: _issue()}, {"scientist": None})
    assert ca.audit_pr_pre_merge(99, DATE) == []


def test_skip_marker_trailing_after_prose_no_longer_opts_out(monkeypatch):
    """The behavior change, pinned (Issue #1126).

    A marker trailing on a prose line used to opt out; it no longer does. Pinning
    it makes the contract change explicit rather than an accident of a fixture
    edit - if someone re-loosens the regex, this test says so.
    """
    pr = _pr(body="routine ship <!-- skip-lab-notebook: routine -->")
    _fake_io(monkeypatch, pr, {42: _issue()}, {"scientist": None})
    assert ca.audit_pr_pre_merge(99, DATE) != []


def test_no_role_label_returns_clean(monkeypatch):
    issue = {"number": 42, "labels": [{"name": "enhancement"}]}
    _fake_io(monkeypatch, _pr(), {42: issue}, {})
    assert ca.audit_pr_pre_merge(99, DATE) == []


def test_no_linked_issues_returns_clean(monkeypatch):
    _fake_io(monkeypatch, _pr(refs=()), {}, {})
    assert ca.audit_pr_pre_merge(99, DATE) == []


def test_multi_role_satisfied_by_one_returns_clean(monkeypatch):
    # Issue tagged both roles; one role's notebook has the entry → clean
    # (any-role-satisfies, per closure_audit.check_lab_notebooks_for_issue).
    issue = {42: _issue(roles=("developer", "scientist"))}
    notebooks = {
        "developer": "## 2026-06-01\n### 10:00 UTC — Editor: X\nclosing #42",
        "scientist": "nothing here",
    }
    _fake_io(monkeypatch, _pr(), issue, notebooks)
    assert ca.audit_pr_pre_merge(99, DATE) == []


def test_multi_role_both_missing_returns_gap(monkeypatch):
    issue = {42: _issue(roles=("developer", "scientist"))}
    _fake_io(monkeypatch, _pr(), issue, {"developer": "nope", "scientist": "nope"})
    assert ca.audit_pr_pre_merge(99, DATE)


# --- lab_notebook_gate.main: CLI exit codes ---


def _run(monkeypatch, argv):
    monkeypatch.setattr(gate.sys, "argv", argv)
    return gate.main()


def test_cli_clean_returns_0(monkeypatch):
    monkeypatch.setattr(ca, "audit_pr_pre_merge", lambda n, today, repo=None: [])
    assert _run(monkeypatch, ["lab_notebook_gate.py", "99"]) == 0


def test_cli_gap_returns_1(monkeypatch, capsys):
    monkeypatch.setattr(
        ca, "audit_pr_pre_merge",
        lambda n, today, repo=None: [
            ("scientist", "no '## 2026-06-01' header in notebook")
        ],
    )
    rc = _run(monkeypatch, ["lab_notebook_gate.py", "99"])
    assert rc == 1
    assert "lab-notebook" in capsys.readouterr().err.lower()


def test_cli_gh_error_fails_open(monkeypatch):
    def boom(n, today, repo=None):
        raise subprocess.CalledProcessError(1, ["gh"])

    monkeypatch.setattr(ca, "audit_pr_pre_merge", boom)
    assert _run(monkeypatch, ["lab_notebook_gate.py", "99"]) == 0


def test_cli_json_error_fails_open(monkeypatch):
    def boom(n, today, repo=None):
        raise json.JSONDecodeError("bad", "doc", 0)

    monkeypatch.setattr(ca, "audit_pr_pre_merge", boom)
    assert _run(monkeypatch, ["lab_notebook_gate.py", "99"]) == 0


def test_cli_os_error_fails_open(monkeypatch):
    # _load_notebook can raise OSError if a notebook file is unreadable
    # (permissions). Fail open — a local FS hiccup must not block a merge.
    def boom(n, today, repo=None):
        raise OSError("permission denied")

    monkeypatch.setattr(ca, "audit_pr_pre_merge", boom)
    assert _run(monkeypatch, ["lab_notebook_gate.py", "99"]) == 0


def test_cli_bad_usage_returns_2(monkeypatch):
    assert _run(monkeypatch, ["lab_notebook_gate.py"]) == 2
    assert _run(monkeypatch, ["lab_notebook_gate.py", "not-a-number"]) == 2


# --- #607: REPO env var read + forwarded to audit_pr_pre_merge ---


def test_cli_reads_repo_env_and_forwards(monkeypatch):
    captured = {}

    def fake(n, today, repo):
        captured["repo"] = repo
        return []

    monkeypatch.setenv("REPO", "fork/repo")
    monkeypatch.setattr(ca, "audit_pr_pre_merge", fake)
    assert _run(monkeypatch, ["lab_notebook_gate.py", "99"]) == 0
    assert captured["repo"] == "fork/repo"


def test_cli_defaults_repo_when_env_unset(monkeypatch):
    """No REPO env → canonical repo default, matching stray_closers/bot_review_offer."""
    captured = {}

    def fake(n, today, repo):
        captured["repo"] = repo
        return []

    monkeypatch.delenv("REPO", raising=False)
    monkeypatch.setattr(ca, "audit_pr_pre_merge", fake)
    assert _run(monkeypatch, ["lab_notebook_gate.py", "99"]) == 0
    assert captured["repo"] == "Jin-HoMLee/splice-neoepitope-pipeline"
