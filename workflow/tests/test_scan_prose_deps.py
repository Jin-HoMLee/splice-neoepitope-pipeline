"""Tests for scripts/pm/scan_prose_deps.py — the prose-dependency reconciler (Issue #722).

Parse layer is pure and gets the bulk of coverage; reconcile/act monkeypatch the
gh-backed helpers so nothing touches the network.
"""
import sys
from pathlib import Path

import pytest

_PM_DIR = Path(__file__).resolve().parents[2] / "scripts" / "pm"
sys.path.insert(0, str(_PM_DIR))

import scan_prose_deps as spd  # noqa: E402


def test_repo_constant():
    assert spd.REPO == "Jin-HoMLee/splice-neoepitope-pipeline"
