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
