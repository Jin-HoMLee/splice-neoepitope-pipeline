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


@pytest.mark.parametrize("phrase", [
    "depends on #722", "blocked by #722", "blocked on #722",
    "blocked-by:722", "gated on #722", "requires #722",
])
def test_parse_allowlist_phrases_match(phrase):
    body = f"Some context. This {phrase} to proceed."
    assert spd.parse_dependencies(745, body) == [(745, 722)]


@pytest.mark.parametrize("phrase", [
    "informs #722", "consumes #722", "relates to #722", "related to #722",
    "follow-up to #722", "supersedes #722", "superseded by #722",
    "see #722", "cf. #722", "fixes #722",
])
def test_parse_narrative_phrases_excluded(phrase):
    # Narrative cross-refs are *why*-context, never blockers — must yield nothing.
    body = f"This PR {phrase}."
    assert spd.parse_dependencies(745, body) == []


def test_parse_multiple_blockers_deduped_and_sorted():
    body = "Depends on #708 and blocked by #707. Also depends on #708 again."
    assert spd.parse_dependencies(709, body) == [(709, 707), (709, 708)]


def test_parse_self_reference_dropped():
    body = "This depends on #745 (itself, nonsensically)."
    assert spd.parse_dependencies(745, body) == []


def test_parse_empty_or_none_body():
    assert spd.parse_dependencies(1, "") == []
    assert spd.parse_dependencies(1, None) == []


def test_parse_case_insensitive():
    assert spd.parse_dependencies(2, "DEPENDS ON #5") == [(2, 5)]


def test_classify_needs_wiring():
    r = spd.classify(745, 722, blocker_meta={"state": "open", "is_pr": False}, existing=set())
    assert r == "needs-wiring"

def test_classify_already_wired():
    r = spd.classify(745, 722, blocker_meta={"state": "open", "is_pr": False}, existing={722})
    assert r == "already-wired"

def test_classify_closed_blocker():
    r = spd.classify(594, 211, blocker_meta={"state": "closed", "is_pr": False}, existing=set())
    assert r == "closed-blocker"

def test_classify_un_wireable_pr():
    r = spd.classify(680, 714, blocker_meta={"state": "open", "is_pr": True}, existing=set())
    assert r == "un-wireable-pr"


def test_reconcile_classifies_each_pair(monkeypatch):
    meta = {
        722: {"state": "open", "is_pr": False},   # needs-wiring
        211: {"state": "closed", "is_pr": False}, # closed-blocker
        714: {"state": "open", "is_pr": True},    # un-wireable-pr
        708: {"state": "open", "is_pr": False},   # already-wired (in existing)
    }
    monkeypatch.setattr(spd, "issue_meta", lambda n: meta[n])
    monkeypatch.setattr(spd, "native_blockers", lambda n: {708} if n == 709 else set())

    pairs = [(745, 722), (594, 211), (680, 714), (709, 708)]
    records = spd.reconcile(pairs)
    actions = {(r["dependent"], r["blocker"]): r["action"] for r in records}
    assert actions == {
        (745, 722): "needs-wiring",
        (594, 211): "closed-blocker",
        (680, 714): "un-wireable-pr",
        (709, 708): "already-wired",
    }
