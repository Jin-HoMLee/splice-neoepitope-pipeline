"""Tests for the pre-merge cross-repo AC gate (Issue #665, gap 2).

The same-repo AC check in audit_and_merge.sh keys off native
closingIssuesReferences, which GitHub only populates within one repo. A
cross-repo close (a personas-repo PR closing a project-repo Issue) is therefore
invisible to it — the project Issue's Acceptance criteria get no gate at all.
This gate parses the PR body for cross-repo closing forward-links and audits
each target Issue's ACs, blocking the merge (exit 1) on any unticked AC box.

Covers the thin CLI wrapper's exit codes; the parsing + collection logic is
single-sourced in closure_audit (tested in test_closure_audit.py).
"""

import json
import subprocess

import closure_audit as ca
import cross_repo_ac_gate as gate


def _run(monkeypatch, argv):
    monkeypatch.setattr(gate.sys, "argv", argv)
    return gate.main()


def test_cli_clean_returns_0(monkeypatch):
    monkeypatch.setattr(ca, "collect_cross_repo_ac_gaps", lambda n, repo=None: [])
    assert _run(monkeypatch, ["cross_repo_ac_gate.py", "99"]) == 0


def test_cli_gap_returns_1(monkeypatch, capsys):
    monkeypatch.setattr(
        ca, "collect_cross_repo_ac_gaps",
        lambda n, repo=None: [
            ("Jin-HoMLee/splice-neoepitope-pipeline#665", "2/3 unticked, no deferral comment found")
        ],
    )
    rc = _run(monkeypatch, ["cross_repo_ac_gate.py", "99"])
    assert rc == 1
    err = capsys.readouterr().err
    assert "#665" in err
    assert "cross-repo" in err.lower()


def test_cli_gh_error_fails_open(monkeypatch):
    def boom(n, repo=None):
        raise subprocess.CalledProcessError(1, ["gh"])

    monkeypatch.setattr(ca, "collect_cross_repo_ac_gaps", boom)
    assert _run(monkeypatch, ["cross_repo_ac_gate.py", "99"]) == 0


def test_cli_json_error_fails_open(monkeypatch):
    def boom(n, repo=None):
        raise json.JSONDecodeError("bad", "doc", 0)

    monkeypatch.setattr(ca, "collect_cross_repo_ac_gaps", boom)
    assert _run(monkeypatch, ["cross_repo_ac_gate.py", "99"]) == 0


def test_cli_bad_usage_returns_2(monkeypatch):
    assert _run(monkeypatch, ["cross_repo_ac_gate.py"]) == 2
    assert _run(monkeypatch, ["cross_repo_ac_gate.py", "not-a-number"]) == 2


# --- REPO env var read + forwarded (mirrors #607 threading) ---


def test_cli_reads_repo_env_and_forwards(monkeypatch):
    captured = {}

    def fake(n, repo=None):
        captured["repo"] = repo
        return []

    monkeypatch.setenv("REPO", "Jin-HoMLee/claude-personas-splice-neoepitope-pipeline")
    monkeypatch.setattr(ca, "collect_cross_repo_ac_gaps", fake)
    assert _run(monkeypatch, ["cross_repo_ac_gate.py", "99"]) == 0
    assert captured["repo"] == "Jin-HoMLee/claude-personas-splice-neoepitope-pipeline"


def test_cli_defaults_repo_when_env_unset(monkeypatch):
    captured = {}

    def fake(n, repo=None):
        captured["repo"] = repo
        return []

    monkeypatch.delenv("REPO", raising=False)
    monkeypatch.setattr(ca, "collect_cross_repo_ac_gaps", fake)
    assert _run(monkeypatch, ["cross_repo_ac_gate.py", "99"]) == 0
    assert captured["repo"] == "Jin-HoMLee/splice-neoepitope-pipeline"
