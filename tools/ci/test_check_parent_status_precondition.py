"""Tests for `.claude/hooks/check_parent_status_precondition.py` (Issue #499).

Pure-function unit tests import the module directly (no network): command parsing
(`status_mutations`) and the precondition decision (`should_refuse`). Orchestration
tests monkeypatch the single `gh` I/O function (`resolve_parent_status`) so
`main()`'s deny/allow decision is exercised without a live `gh` call. The
subprocess fail-open tests mirror the harness invocation (pipe PreToolUse JSON to
stdin) and only travel paths that return before any `gh` call.

The full live sub-with-Backlog-parent refused path is an integration check,
verified by smoke-testing the hook against a real sub-issue whose parent sits in
Backlog (documented in the PR Test plan) — not in this suite.
"""

import io
import json
import subprocess
import sys
from pathlib import Path

HOOKS_DIR = Path(__file__).parent.parent.parent / ".claude" / "hooks"
HOOK = HOOKS_DIR / "check_parent_status_precondition.py"
sys.path.insert(0, str(HOOKS_DIR))

import check_parent_status_precondition as h  # noqa: E402

STATUS_FIELD = "PVTSSF_lAHOB17eGc4BSomPzhAHFf8"
IN_PROGRESS = "47fc9ee4"
IN_REVIEW = "df73e18b"
READY_FOR_REVIEW = "8bf9192f"
READY = "61e4505c"
BACKLOG = "f75ad846"
DONE = "98236657"


def _mut(item_id, option_id, field_id=STATUS_FIELD):
    """A single spaced updateProjectV2ItemFieldValue mutation command."""
    return (
        "gh api graphql -f query='mutation { updateProjectV2ItemFieldValue("
        f'input: {{projectId: "PVT_kwHOB17eGc4BSomP", itemId: "{item_id}", '
        f'fieldId: "{field_id}", value: {{singleSelectOptionId: "{option_id}"}}}})'
        " { projectV2Item { id } } }'"
    )


# --- status_mutations: command parsing + target filtering ---


class TestStatusMutations:
    def test_advanced_target_extracted(self):
        assert h.status_mutations(_mut("PVTI_abc", IN_PROGRESS)) == [
            ("PVTI_abc", IN_PROGRESS)
        ]

    def test_in_review_extracted(self):
        assert h.status_mutations(_mut("PVTI_xyz", IN_REVIEW)) == [
            ("PVTI_xyz", IN_REVIEW)
        ]

    def test_ready_for_review_extracted(self):
        assert h.status_mutations(_mut("PVTI_xyz", READY_FOR_REVIEW)) == [
            ("PVTI_xyz", READY_FOR_REVIEW)
        ]

    def test_backlog_target_ignored(self):
        # moving a sub → Backlog has no parent precondition
        assert h.status_mutations(_mut("PVTI_abc", BACKLOG)) == []

    def test_ready_target_ignored(self):
        assert h.status_mutations(_mut("PVTI_abc", READY)) == []

    def test_done_target_ignored(self):
        assert h.status_mutations(_mut("PVTI_abc", DONE)) == []

    def test_non_status_field_ignored(self):
        # a write to some other single-select field is not our concern
        other = "PVTSSF_someOtherFieldId00000000"
        assert h.status_mutations(_mut("PVTI_abc", IN_PROGRESS, field_id=other)) == []

    def test_compact_no_spaces(self):
        cmd = (
            'gh api graphql -f query="mutation{updateProjectV2ItemFieldValue('
            'input:{projectId:\\"PVT_kwHOB17eGc4BSomP\\",itemId:\\"PVTI_c\\",'
            f'fieldId:\\"{STATUS_FIELD}\\",value:{{singleSelectOptionId:'
            f'\\"{IN_PROGRESS}\\"}}}}){{projectV2Item{{id}}}}}}"'
        )
        assert h.status_mutations(cmd) == [("PVTI_c", IN_PROGRESS)]

    def test_aliased_dual_mutation_both_blocks(self):
        # the PR+Issue lifecycle flip: two updateProjectV2ItemFieldValue in one call
        cmd = (
            "gh api graphql -f query='mutation { "
            "setPR: updateProjectV2ItemFieldValue(input: {projectId: \"P\", "
            f'itemId: "PVTI_pr", fieldId: "{STATUS_FIELD}", '
            f'value: {{singleSelectOptionId: "{IN_REVIEW}"}}}}) {{ projectV2Item {{ id }} }} '
            "setIssue: updateProjectV2ItemFieldValue(input: {projectId: \"P\", "
            f'itemId: "PVTI_iss", fieldId: "{STATUS_FIELD}", '
            f'value: {{singleSelectOptionId: "{IN_REVIEW}"}}}}) {{ projectV2Item {{ id }} }} }}\''
        )
        assert h.status_mutations(cmd) == [
            ("PVTI_pr", IN_REVIEW),
            ("PVTI_iss", IN_REVIEW),
        ]

    def test_not_a_mutation(self):
        assert h.status_mutations("gh api graphql -f query='query { viewer { login } }'") == []

    def test_not_gh_api(self):
        assert h.status_mutations("gh issue view 499") == []

    def test_mutation_inside_quoted_body_not_matched(self):
        # a PR-comment body quoting the mutation must NOT match (command-start gate)
        cmd = (
            'gh pr comment 1 --body "we run updateProjectV2ItemFieldValue with '
            f'fieldId {STATUS_FIELD} and singleSelectOptionId {IN_PROGRESS}"'
        )
        assert h.status_mutations(cmd) == []

    def test_after_separator_matched(self):
        cmd = "echo hi && " + _mut("PVTI_abc", IN_PROGRESS)
        assert h.status_mutations(cmd) == [("PVTI_abc", IN_PROGRESS)]

    def test_unbalanced_quotes_fail_safe(self):
        assert h.status_mutations("gh api graphql -f query='mutation {") == []


# --- should_refuse: precondition decision ---


class TestShouldRefuse:
    def test_backlog_parent_refused(self):
        assert h.should_refuse(BACKLOG) is True

    def test_no_status_parent_refused(self):
        # parent on board but no Status value set (No Status) → blocking
        assert h.should_refuse(None) is True

    def test_in_progress_parent_allowed(self):
        assert h.should_refuse(IN_PROGRESS) is False

    def test_ready_parent_allowed(self):
        # parent committed (Ready) but not yet started — advancing a sub is allowed
        assert h.should_refuse(READY) is False

    def test_in_review_parent_allowed(self):
        assert h.should_refuse(IN_REVIEW) is False

    def test_done_parent_allowed(self):
        assert h.should_refuse(DONE) is False


# --- deny_payload: shape ---


def test_deny_payload_shape():
    payload = h.deny_payload(592, 547, "Backlog", "In progress")
    out = payload["hookSpecificOutput"]
    assert out["hookEventName"] == "PreToolUse"
    assert out["permissionDecision"] == "deny"
    reason = out["permissionDecisionReason"]
    assert "#592" in reason and "#547" in reason
    assert "Backlog" in reason
    assert "feedback_parent_sub_issues.md" in reason


# --- main(): orchestration (monkeypatched gh I/O, no network) ---


def _run_main(monkeypatch, capsys, cmd, resolver):
    """Drive main() with a stubbed resolve_parent_status; return (rc, stdout|None)."""
    monkeypatch.setattr(h, "resolve_parent_status", resolver)
    payload = json.dumps({"tool_input": {"command": cmd}})
    monkeypatch.setattr(sys, "stdin", io.StringIO(payload))
    rc = h.main()
    captured = capsys.readouterr().out.strip()
    return rc, (json.loads(captured) if captured else None)


def test_main_denies_sub_with_backlog_parent(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(h, "LOG_PATH", tmp_path / "hook_fires.jsonl")
    # resolver: this item is sub #592, parent #547 in Backlog
    resolver = lambda item_id: (592, 547, BACKLOG)
    rc, out = _run_main(monkeypatch, capsys, _mut("PVTI_592", IN_PROGRESS), resolver)
    assert rc == 0
    assert out["hookSpecificOutput"]["permissionDecision"] == "deny"
    reason = out["hookSpecificOutput"]["permissionDecisionReason"]
    assert "#592" in reason and "#547" in reason
    logged = json.loads((tmp_path / "hook_fires.jsonl").read_text().strip())
    assert logged["hook"] == "check_parent_status_precondition"
    assert logged["sub"] == 592 and logged["parent"] == 547


def test_main_allows_sub_with_in_progress_parent(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(h, "LOG_PATH", tmp_path / "hook_fires.jsonl")
    resolver = lambda item_id: (592, 547, IN_PROGRESS)
    rc, out = _run_main(monkeypatch, capsys, _mut("PVTI_592", IN_PROGRESS), resolver)
    assert rc == 0 and out is None
    assert not (tmp_path / "hook_fires.jsonl").exists()


def test_main_allows_when_no_parent(monkeypatch, capsys):
    # resolver returns None → standalone issue / not a sub / PR item → allow
    rc, out = _run_main(monkeypatch, capsys, _mut("PVTI_x", IN_PROGRESS),
                        lambda item_id: None)
    assert rc == 0 and out is None


def test_main_fails_open_on_resolver_error(monkeypatch, capsys):
    # resolver returns None on any gh/network error → allow
    rc, out = _run_main(monkeypatch, capsys, _mut("PVTI_x", IN_PROGRESS),
                        lambda item_id: None)
    assert rc == 0 and out is None


def test_main_noop_on_non_advanced_target(monkeypatch, capsys):
    # → Backlog: status_mutations yields nothing, resolver never consulted
    def _boom(item_id):
        raise AssertionError("resolver must not be called for a Backlog target")
    rc, out = _run_main(monkeypatch, capsys, _mut("PVTI_x", BACKLOG), _boom)
    assert rc == 0 and out is None


def test_main_dual_mutation_denies_on_offending_block(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(h, "LOG_PATH", tmp_path / "hook_fires.jsonl")
    # PR item (PVTI_pr) → no parent; Issue item (PVTI_iss) → sub with Backlog parent
    def resolver(item_id):
        if item_id == "PVTI_iss":
            return (592, 547, BACKLOG)
        return None  # PR item: content is a PullRequest, not an Issue
    cmd = (
        "gh api graphql -f query='mutation { "
        "setPR: updateProjectV2ItemFieldValue(input: {projectId: \"P\", "
        f'itemId: "PVTI_pr", fieldId: "{STATUS_FIELD}", '
        f'value: {{singleSelectOptionId: "{IN_REVIEW}"}}}}) {{ projectV2Item {{ id }} }} '
        "setIssue: updateProjectV2ItemFieldValue(input: {projectId: \"P\", "
        f'itemId: "PVTI_iss", fieldId: "{STATUS_FIELD}", '
        f'value: {{singleSelectOptionId: "{IN_REVIEW}"}}}}) {{ projectV2Item {{ id }} }} }}\''
    )
    rc, out = _run_main(monkeypatch, capsys, cmd, resolver)
    assert rc == 0
    assert out["hookSpecificOutput"]["permissionDecision"] == "deny"
    assert "#592" in out["hookSpecificOutput"]["permissionDecisionReason"]


# --- subprocess fail-open / no-op (no live gh reached) ---


def _run_subprocess(stdin_payload: str):
    proc = subprocess.run(
        [sys.executable, str(HOOK)],
        input=stdin_payload,
        capture_output=True,
        text=True,
        timeout=10,
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

    def test_non_mutation_command_noop(self):
        # returns before any gh call (status_mutations is empty)
        rc, out, _ = _run_subprocess(_payload("gh issue view 499"))
        assert rc == 0 and out.strip() == ""

    def test_backlog_target_noop(self):
        # advanced-target filter screens it out before any gh call
        rc, out, _ = _run_subprocess(_payload(_mut("PVTI_x", BACKLOG)))
        assert rc == 0 and out.strip() == ""
