"""Tests for the pre-merge parent/child gate (Issue #1155).

Moves the parent guard from branch-time to merge-time: refuse a merge whose PR
closes a **parent** Issue (`subIssuesSummary.total > 0`) that still has >= 1
**open** child, naming the open children. A PR closing a parent auto-closes it on
merge and orphans the open children (the PR #543 -> Issue #538 harm), so the harm
lands at merge, which is where we now guard.

Fail posture is **B (fail-closed on the parent check)**, diverging from the
sibling gates: if a linked Issue's parent/child state cannot be determined (a
gh/network error), the gate BLOCKS with a "could not verify" message rather than
admitting the merge - the orphaning harm is irreversible, so a transient
false-block is cheaper than a false-pass. Fail-OPEN only where there is genuinely
nothing to check: a PR with no closing references, or a linked Issue that is a
leaf (not a parent).

The pure decision logic (`gate_decision`) is unit-tested with no I/O. The
two-directional falsifier pair (AC-3) drives `main()` with stubbed gh I/O: the
open-child case must go RED (exit 1), the all-children-closed case GREEN (exit 0)
- a gate that blocked everything, or blocked nothing, fails this pair.
"""

import json
import subprocess

import parent_child_gate as gate


# --- gate_decision: pure aggregation (no I/O) ---


class TestGateDecision:
    def test_open_children_blocks_and_names_them(self):
        blocked, messages = gate.gate_decision(
            [{"number": 1031, "status": "open_children", "open": [1155, 1160]}]
        )
        assert blocked is True
        assert len(messages) == 1
        assert "#1031" in messages[0]
        assert "#1155" in messages[0] and "#1160" in messages[0]

    def test_clean_parent_does_not_block(self):
        blocked, messages = gate.gate_decision(
            [{"number": 1026, "status": "clean"}]
        )
        assert blocked is False
        assert messages == []

    def test_leaf_does_not_block(self):
        blocked, messages = gate.gate_decision(
            [{"number": 549, "status": "leaf"}]
        )
        assert blocked is False
        assert messages == []

    def test_undetermined_fails_closed(self):
        blocked, messages = gate.gate_decision(
            [{"number": 538, "status": "undetermined"}]
        )
        assert blocked is True
        assert len(messages) == 1
        assert "#538" in messages[0]
        assert "verify" in messages[0].lower()

    def test_empty_does_not_block(self):
        blocked, messages = gate.gate_decision([])
        assert blocked is False
        assert messages == []

    def test_mixed_blocks_on_the_bad_one_only(self):
        blocked, messages = gate.gate_decision([
            {"number": 549, "status": "leaf"},
            {"number": 1026, "status": "clean"},
            {"number": 1031, "status": "open_children", "open": [1155]},
        ])
        assert blocked is True
        assert len(messages) == 1
        assert "#1031" in messages[0]


# --- resolve: classify one linked Issue from fetched state (I/O seam stubbed) ---


class TestResolve:
    def test_leaf_when_no_sub_issues(self, monkeypatch):
        monkeypatch.setattr(gate, "_fetch_parent_state", lambda n, repo: (0, []))
        assert gate.resolve(549, "o/r") == {"number": 549, "status": "leaf"}

    def test_clean_when_all_children_closed(self, monkeypatch):
        monkeypatch.setattr(
            gate, "_fetch_parent_state",
            lambda n, repo: (2, [{"number": 100, "state": "CLOSED"},
                                 {"number": 101, "state": "CLOSED"}]),
        )
        assert gate.resolve(1026, "o/r") == {"number": 1026, "status": "clean"}

    def test_open_children_lists_only_the_open_ones(self, monkeypatch):
        monkeypatch.setattr(
            gate, "_fetch_parent_state",
            lambda n, repo: (3, [{"number": 100, "state": "CLOSED"},
                                 {"number": 101, "state": "OPEN"},
                                 {"number": 102, "state": "OPEN"}]),
        )
        r = gate.resolve(1031, "o/r")
        assert r["number"] == 1031
        assert r["status"] == "open_children"
        assert r["open"] == [101, 102]

    def test_gh_error_maps_to_undetermined_fail_closed(self, monkeypatch):
        def boom(n, repo):
            raise subprocess.CalledProcessError(1, "gh")
        monkeypatch.setattr(gate, "_fetch_parent_state", boom)
        assert gate.resolve(538, "o/r") == {"number": 538, "status": "undetermined"}


# --- main(): the two-directional falsifier pair (AC-3) + exit codes ---


def _run_main(monkeypatch, *, linked, resolutions, tmp_path=None):
    """Drive main() with stubbed gh I/O; return exit code."""
    monkeypatch.setattr(gate.sys, "argv", ["parent_child_gate.py", "999"])
    monkeypatch.setattr(gate, "fetch_closing_issues", lambda pr, repo: linked)
    by_num = {r["number"]: r for r in resolutions}
    monkeypatch.setattr(gate, "resolve", lambda n, repo: by_num[n])
    if tmp_path is not None:
        monkeypatch.setattr(gate, "LOG_PATH", tmp_path / "hook_fires.jsonl")
    return gate.main()


def test_open_child_case_goes_red(monkeypatch, tmp_path):
    # Falsifier half A: a parent with an open child MUST block.
    rc = _run_main(
        monkeypatch,
        linked=[1031],
        resolutions=[{"number": 1031, "status": "open_children", "open": [1155]}],
        tmp_path=tmp_path,
    )
    assert rc == 1


def test_all_children_closed_case_goes_green(monkeypatch):
    # Falsifier half B: a parent whose children are all closed MUST pass.
    rc = _run_main(
        monkeypatch,
        linked=[1026],
        resolutions=[{"number": 1026, "status": "clean"}],
    )
    assert rc == 0


def test_no_closing_refs_fails_open(monkeypatch):
    rc = _run_main(monkeypatch, linked=[], resolutions=[])
    assert rc == 0


def test_leaf_linked_issue_fails_open(monkeypatch):
    rc = _run_main(
        monkeypatch,
        linked=[549],
        resolutions=[{"number": 549, "status": "leaf"}],
    )
    assert rc == 0


def test_undetermined_parent_state_fails_closed(monkeypatch, tmp_path):
    rc = _run_main(
        monkeypatch,
        linked=[538],
        resolutions=[{"number": 538, "status": "undetermined"}],
        tmp_path=tmp_path,
    )
    assert rc == 1


def test_pr_fetch_error_fails_closed(monkeypatch):
    # Cannot even enumerate the closing set -> cannot prove the merge is safe.
    monkeypatch.setattr(gate.sys, "argv", ["parent_child_gate.py", "999"])

    def boom(pr, repo):
        raise subprocess.CalledProcessError(1, "gh")
    monkeypatch.setattr(gate, "fetch_closing_issues", boom)
    assert gate.main() == 1


def test_block_writes_fire_log(monkeypatch, tmp_path):
    _run_main(
        monkeypatch,
        linked=[1031],
        resolutions=[{"number": 1031, "status": "open_children", "open": [1155]}],
        tmp_path=tmp_path,
    )
    logged = (tmp_path / "hook_fires.jsonl").read_text().strip()
    rec = json.loads(logged)
    assert rec["hook"] == "parent_child_gate"
    assert rec["issue"] == 1031


# --- usage ---


def test_usage_error_on_non_numeric(monkeypatch):
    monkeypatch.setattr(gate.sys, "argv", ["parent_child_gate.py", "nope"])
    assert gate.main() == 2


def test_usage_error_on_missing_arg(monkeypatch):
    monkeypatch.setattr(gate.sys, "argv", ["parent_child_gate.py"])
    assert gate.main() == 2
