"""Tests for recheck_dispatch.py fire-log helpers (Issue #453)."""
import json
import sys
from pathlib import Path

import pytest

# Make .agents/hooks importable
_HOOKS_DIR = Path(__file__).resolve().parents[2] / ".agents" / "hooks"
sys.path.insert(0, str(_HOOKS_DIR))


def test_log_fire_jsonl_append(tmp_path, monkeypatch):
    """_log_fire appends one valid JSONL line per call; metadata kwargs preserved."""
    import recheck_dispatch

    log_path = tmp_path / "hook_fires.jsonl"
    monkeypatch.setattr(recheck_dispatch, "LOG_PATH", log_path)

    recheck_dispatch._log_fire("hook_A", issue=1)
    recheck_dispatch._log_fire("hook_A", issue=2, extra="meta")
    recheck_dispatch._log_fire("hook_A", issue=3)

    lines = log_path.read_text().splitlines()
    assert len(lines) == 3
    for line in lines:
        rec = json.loads(line)
        assert rec["hook"] == "hook_A"
        assert "ts" in rec
        assert "issue" in rec
    rec2 = json.loads(lines[1])
    assert rec2.get("extra") == "meta"
    assert rec2.get("issue") == 2


def test_count_fires_filters_by_hook(tmp_path, monkeypatch):
    """_count_fires returns per-hook counts; unknown hook returns 0; missing log returns 0."""
    import recheck_dispatch

    log_path = tmp_path / "hook_fires.jsonl"
    monkeypatch.setattr(recheck_dispatch, "LOG_PATH", log_path)

    # Missing log → 0 for any hook
    assert recheck_dispatch._count_fires("hook_A") == 0

    recheck_dispatch._log_fire("hook_A", issue=1)
    recheck_dispatch._log_fire("hook_A", issue=2)
    recheck_dispatch._log_fire("hook_B", issue=3)

    assert recheck_dispatch._count_fires("hook_A") == 2
    assert recheck_dispatch._count_fires("hook_B") == 1
    assert recheck_dispatch._count_fires("hook_C") == 0


def test_threshold_prompt_equals_once(tmp_path, monkeypatch):
    """_threshold_prompt fires exactly at count == K, never before or after."""
    import recheck_dispatch

    log_path = tmp_path / "hook_fires.jsonl"
    monkeypatch.setattr(recheck_dispatch, "LOG_PATH", log_path)
    monkeypatch.setattr(recheck_dispatch, "HOOK_CONFIG", {
        "test_hook": {"threshold": 3, "dock": 999},
    })

    recheck_dispatch._log_fire("test_hook", issue=1)
    assert recheck_dispatch._threshold_prompt("test_hook") is None

    recheck_dispatch._log_fire("test_hook", issue=2)
    assert recheck_dispatch._threshold_prompt("test_hook") is None

    recheck_dispatch._log_fire("test_hook", issue=3)
    prompt = recheck_dispatch._threshold_prompt("test_hook")
    assert prompt is not None
    assert "test_hook" in prompt
    assert "3 times" in prompt
    assert "999" in prompt
    assert "check_hook_health.sh" in prompt

    recheck_dispatch._log_fire("test_hook", issue=4)
    # == K, not >=; 4 fires no longer matches.
    assert recheck_dispatch._threshold_prompt("test_hook") is None


def test_threshold_prompt_no_dock(tmp_path, monkeypatch):
    """_threshold_prompt emits '(no dock Issue filed yet)' when dock is None."""
    import recheck_dispatch

    log_path = tmp_path / "hook_fires.jsonl"
    monkeypatch.setattr(recheck_dispatch, "LOG_PATH", log_path)
    monkeypatch.setattr(recheck_dispatch, "HOOK_CONFIG", {
        "dockless_hook": {"threshold": 1, "dock": None},
    })

    recheck_dispatch._log_fire("dockless_hook", issue=1)
    prompt = recheck_dispatch._threshold_prompt("dockless_hook")
    assert prompt is not None
    assert "no dock Issue filed yet" in prompt


def test_is_fire_target_sync_check():
    """target_sync_check predicate: empty is no-fire; a manual-param warning is a fire;
    a successful auto-sync confirmation is NOT a fire (Route A #782 — the mechanism working
    as intended must not inflate the fire-log / trip the promotion prompt)."""
    import recheck_dispatch
    assert recheck_dispatch._is_fire("target_sync_check", "") is False
    assert recheck_dispatch._is_fire("target_sync_check", "[target re-sync needed — move on #11]") is True
    assert recheck_dispatch._is_fire("target_sync_check", '[target auto-synced — #11 → milestone "i5", Target 2026-07-15]') is False
    assert recheck_dispatch._is_fire("target_sync_check", "[target auto-synced — #11 demilestoned, Target cleared]") is False


def test_is_fire_recheck_milestone_no_change_paths():
    """recheck_milestone predicate: 'Status: [No change]' and 'has no milestone' are both no-fire."""
    import recheck_dispatch
    assert recheck_dispatch._is_fire("recheck_milestone", "Status: [No change]") is False
    assert recheck_dispatch._is_fire("recheck_milestone", "Status: [No change] — milestone has no remaining capacity") is False
    assert recheck_dispatch._is_fire("recheck_milestone", "error: issue #999 has no milestone") is False
    assert recheck_dispatch._is_fire("recheck_milestone", "Status: [UPDATE NEEDED]") is True
    assert recheck_dispatch._is_fire("recheck_milestone", "Status: [UNSIZED]") is True


def test_is_fire_recheck_parent_status_no_change_paths():
    """recheck_parent_status predicate: 'Status: [No change]' and 'has no parent' are both no-fire."""
    import recheck_dispatch
    assert recheck_dispatch._is_fire("recheck_parent_status", "Status: [No change]") is False
    assert recheck_dispatch._is_fire("recheck_parent_status", "Issue #999 has no parent — nothing to audit.") is False
    assert recheck_dispatch._is_fire("recheck_parent_status", "Status: [PARENT_AHEAD]") is True


def test_is_fire_unknown_hook_returns_false():
    """Unknown hook names default to no-fire (safe-by-default)."""
    import recheck_dispatch
    assert recheck_dispatch._is_fire("unknown_hook", "any text") is False


# ---------------------------------------------------------------------------
# Scope filter (Issue #454)
# ---------------------------------------------------------------------------

def test_parse_scope_mappings():
    """--scope {shared,pm} narrow the active set; absent/all/unknown fail open to everything."""
    import recheck_dispatch as rd
    assert rd._parse_scope([]) == {"shared", "pm"}                  # flagless → all (backward-compat)
    assert rd._parse_scope(["--scope", "all"]) == {"shared", "pm"}
    assert rd._parse_scope(["--scope", "shared"]) == {"shared"}
    assert rd._parse_scope(["--scope", "pm"]) == {"pm"}
    assert rd._parse_scope(["--scope=shared"]) == {"shared"}        # = form
    assert rd._parse_scope(["--scope", "bogus"]) == {"shared", "pm"}  # unknown → fail open


def test_in_scope_honors_active_scopes(monkeypatch):
    """_in_scope gates each check by its HOOK_CONFIG scope against ACTIVE_SCOPES."""
    import recheck_dispatch as rd

    monkeypatch.setattr(rd, "ACTIVE_SCOPES", {"shared"})
    assert rd._in_scope("target_sync_check") is True          # shared check runs in shared session
    assert rd._in_scope("recheck_milestone") is False         # pm check suppressed for Sci/Dev
    assert rd._in_scope("recheck_parent_status") is False

    monkeypatch.setattr(rd, "ACTIVE_SCOPES", {"pm"})
    assert rd._in_scope("target_sync_check") is False
    assert rd._in_scope("recheck_milestone") is True
    assert rd._in_scope("recheck_parent_status") is True

    monkeypatch.setattr(rd, "ACTIVE_SCOPES", {"shared", "pm"})  # PM session sees both
    assert rd._in_scope("target_sync_check") is True
    assert rd._in_scope("recheck_milestone") is True


def test_in_scope_unknown_hook_fails_open(monkeypatch):
    """A check with no 'scope' key (or unknown name) runs everywhere (fail-open)."""
    import recheck_dispatch as rd
    monkeypatch.setattr(rd, "ACTIVE_SCOPES", {"shared"})
    monkeypatch.setattr(rd, "HOOK_CONFIG", {"scopeless": {"threshold": 1, "dock": None}})
    assert rd._in_scope("scopeless") is True       # no scope key → always in scope
    assert rd._in_scope("never_configured") is True  # unknown name → always in scope


# ---------------------------------------------------------------------------
# Target-date auto-apply — Route A (Issue #782)
# ---------------------------------------------------------------------------

def _stub_getters(monkeypatch, *, ms_title, ms_due_on, target_date, item_id):
    """Stub get_issue_milestone / get_issue_target_date so apply_target_sync is offline."""
    import recheck_dispatch as rd
    monkeypatch.setattr(rd, "get_issue_milestone", lambda issue: (ms_title, ms_due_on))
    monkeypatch.setattr(rd, "get_issue_target_date", lambda issue: (target_date, item_id))


def test_apply_target_sync_updates_when_out_of_sync(monkeypatch):
    """Out-of-sync milestoned issue → executes the update mutation, returns a confirmation."""
    import recheck_dispatch as rd
    _stub_getters(monkeypatch, ms_title="i5 - S3", ms_due_on="2026-07-15",
                  target_date="2026-06-30", item_id="PVTI_x")
    calls = []
    monkeypatch.setattr(rd, "_mutate_target_date",
                        lambda item_id, due_on: calls.append((item_id, due_on)) or True)
    out = rd.apply_target_sync(782)
    assert calls == [("PVTI_x", "2026-07-15")]       # mutated with the new due_on
    assert out is not None
    assert "auto-synced" in out
    assert "2026-07-15" in out


def test_apply_target_sync_noop_when_already_synced(monkeypatch):
    """Target already == due_on → None, no mutation (idempotent)."""
    import recheck_dispatch as rd
    _stub_getters(monkeypatch, ms_title="i5", ms_due_on="2026-07-15",
                  target_date="2026-07-15", item_id="PVTI_x")
    calls = []
    monkeypatch.setattr(rd, "_mutate_target_date", lambda *a: calls.append(a) or True)
    assert rd.apply_target_sync(782) is None
    assert calls == []


def test_apply_target_sync_clears_on_demilestone(monkeypatch):
    """Un-milestoned but Target still set → clear mutation (due_on=None)."""
    import recheck_dispatch as rd
    _stub_getters(monkeypatch, ms_title=None, ms_due_on=None,
                  target_date="2026-06-30", item_id="PVTI_x")
    calls = []
    monkeypatch.setattr(rd, "_mutate_target_date",
                        lambda item_id, due_on: calls.append((item_id, due_on)) or True)
    out = rd.apply_target_sync(782)
    assert calls == [("PVTI_x", None)]               # cleared
    assert out is not None
    assert ("cleared" in out.lower()) or ("demilestone" in out.lower())


def test_apply_target_sync_not_on_board(monkeypatch):
    """Issue not on project board (item_id None) → None, no mutation."""
    import recheck_dispatch as rd
    _stub_getters(monkeypatch, ms_title="i5", ms_due_on="2026-07-15",
                  target_date=None, item_id=None)
    calls = []
    monkeypatch.setattr(rd, "_mutate_target_date", lambda *a: calls.append(a) or True)
    assert rd.apply_target_sync(782) is None
    assert calls == []


def test_apply_target_sync_fails_open_to_manual_surface(monkeypatch):
    """Mutation failure → falls back to target_sync_check's manual-param warning (nothing lost)."""
    import recheck_dispatch as rd
    _stub_getters(monkeypatch, ms_title="i5 - S3", ms_due_on="2026-07-15",
                  target_date="2026-06-30", item_id="PVTI_x")
    monkeypatch.setattr(rd, "_mutate_target_date", lambda *a: False)
    out = rd.apply_target_sync(782)
    assert out is not None
    assert "re-sync needed" in out                   # the manual surface, unchanged
    assert "updateProjectV2ItemFieldValue" in out


def test_mutate_target_date_success(monkeypatch):
    """_mutate_target_date → True on a clean gh api graphql response."""
    import recheck_dispatch as rd

    class R:
        returncode = 0
        stdout = '{"data":{"updateProjectV2ItemFieldValue":{"projectV2Item":{"id":"x"}}}}'
        stderr = ""
    monkeypatch.setattr(rd.subprocess, "run", lambda *a, **k: R())
    assert rd._mutate_target_date("PVTI_x", "2026-07-15") is True


def test_mutate_target_date_clear_success(monkeypatch):
    """_mutate_target_date(item_id, None) → exercises the clearProjectV2ItemFieldValue branch
    at function level (the demilestone path is otherwise only tested with the mutation patched out)."""
    import recheck_dispatch as rd

    seen = {}

    class R:
        returncode = 0
        stdout = '{"data":{"clearProjectV2ItemFieldValue":{"projectV2Item":{"id":"x"}}}}'
        stderr = ""

    def fake_run(args, **k):
        seen["query"] = args[args.index("-f") + 1]
        return R()
    monkeypatch.setattr(rd.subprocess, "run", fake_run)
    assert rd._mutate_target_date("PVTI_x", None) is True
    assert "clearProjectV2ItemFieldValue" in seen["query"]   # took the clear branch, not update


def test_mutate_target_date_rejects_malformed_item_id(monkeypatch):
    """A pathological item_id never reaches the GraphQL string — guarded to fail-open (False)."""
    import recheck_dispatch as rd
    called = []
    monkeypatch.setattr(rd.subprocess, "run", lambda *a, **k: called.append(a))
    assert rd._mutate_target_date('PVTI_x" injected', "2026-07-15") is False
    assert rd._mutate_target_date("not-an-id", "2026-07-15") is False
    assert called == []                                      # never built/fired a mutation


def test_apply_target_sync_milestoned_no_due_date(monkeypatch):
    """Milestoned but the milestone has no due_on → clear Target (intentional policy:
    Target tracks due_on, and an undated milestone has none). Assert the clear is issued."""
    import recheck_dispatch as rd
    _stub_getters(monkeypatch, ms_title="i5 - S3", ms_due_on=None,
                  target_date="2026-06-30", item_id="PVTI_x")
    calls = []
    monkeypatch.setattr(rd, "_mutate_target_date",
                        lambda item_id, due_on: calls.append((item_id, due_on)) or True)
    out = rd.apply_target_sync(782)
    assert calls == [("PVTI_x", None)]                       # cleared because due_on is None
    assert out is not None
    assert "auto-synced" in out


def test_mutate_target_date_graphql_error_is_failure(monkeypatch):
    """A GraphQL `errors` payload (e.g. insufficient scope) → False → triggers fail-open."""
    import recheck_dispatch as rd

    class R:
        returncode = 0
        stdout = '{"errors":[{"message":"insufficient scope"}]}'
        stderr = "insufficient scope"
    monkeypatch.setattr(rd.subprocess, "run", lambda *a, **k: R())
    assert rd._mutate_target_date("PVTI_x", "2026-07-15") is False


def test_mutate_target_date_nonzero_returncode_is_failure(monkeypatch):
    """Non-zero gh exit (HTTP error) → False."""
    import recheck_dispatch as rd

    class R:
        returncode = 1
        stdout = ""
        stderr = "HTTP 401"
    monkeypatch.setattr(rd.subprocess, "run", lambda *a, **k: R())
    assert rd._mutate_target_date("PVTI_x", "2026-07-15") is False


def test_dispatch_target_sync_covers_all_batched_moves(monkeypatch):
    """Batched `gh issue edit N --milestone && ...` → apply_target_sync runs for EVERY issue (AC3)."""
    import recheck_dispatch as rd
    monkeypatch.setattr(rd, "ACTIVE_SCOPES", {"shared"})    # isolate: PM capacity check off
    seen = []
    monkeypatch.setattr(rd, "apply_target_sync", lambda issue: seen.append(issue) or None)
    cmd = ('gh issue edit 11 --milestone "i5 - S3" && '
           'gh issue edit 22 --milestone "i5 - S3" && '
           'gh issue edit 33 --milestone "i5 - S3"')
    rd.dispatch(cmd)
    assert seen == [11, 22, 33]


def test_dispatch_target_sync_dedups_repeated_issue(monkeypatch):
    """Same issue edited twice in one command → apply_target_sync called once."""
    import recheck_dispatch as rd
    monkeypatch.setattr(rd, "ACTIVE_SCOPES", {"shared"})
    seen = []
    monkeypatch.setattr(rd, "apply_target_sync", lambda issue: seen.append(issue) or None)
    cmd = 'gh issue edit 11 --milestone "A" && gh issue edit 11 --milestone "B"'
    rd.dispatch(cmd)
    assert seen == [11]
