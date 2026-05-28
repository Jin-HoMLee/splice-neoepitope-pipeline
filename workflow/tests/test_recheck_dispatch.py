"""Tests for recheck_dispatch.py fire-log helpers (Issue #453)."""
import json
import sys
from pathlib import Path

import pytest

# Make .claude/hooks importable
_HOOKS_DIR = Path(__file__).resolve().parents[2] / ".claude" / "hooks"
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
