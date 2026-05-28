"""Tests for recheck_dispatch.py fire-log helpers (Issue #453)."""
import json
import sys
from pathlib import Path

import pytest

# Make .claude/hooks importable
_HOOKS_DIR = Path(__file__).resolve().parents[2] / ".claude" / "hooks"
sys.path.insert(0, str(_HOOKS_DIR))


def test_placeholder():
    """Removed in Task 2; here to verify test discovery + import path setup."""
    import recheck_dispatch
    assert hasattr(recheck_dispatch, "dispatch")
