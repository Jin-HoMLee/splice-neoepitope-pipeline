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


def _records():
    return [
        {"dependent": 745, "blocker": 722, "state": "open", "action": "needs-wiring"},
        {"dependent": 594, "blocker": 211, "state": "closed", "action": "closed-blocker"},
    ]

def test_render_report_lists_each_record():
    out = spd.render_report(_records())
    assert "745" in out and "722" in out and "needs-wiring" in out
    assert "594" in out and "closed-blocker" in out

def test_render_report_empty():
    assert "no prose-dependency drift" in spd.render_report([]).lower()

def _run_main(monkeypatch, capsys, argv, pairs, records):
    monkeypatch.setattr(spd, "fetch_open_issues",
                        lambda: [{"number": 745, "title": "x", "body": "depends on #722", "state": "OPEN"}])
    monkeypatch.setattr(spd, "parse_dependencies", lambda n, b: pairs)
    monkeypatch.setattr(spd, "reconcile", lambda p: records)
    monkeypatch.setattr(sys, "argv", ["scan_prose_deps.py", *argv])
    rc = spd.main()
    return rc, capsys.readouterr().out

def test_main_report_default_exit_zero(monkeypatch, capsys):
    rc, out = _run_main(monkeypatch, capsys, [], [(745, 722)], _records())
    assert rc == 0 and "needs-wiring" in out

def test_main_check_exits_2_on_drift(monkeypatch, capsys):
    rc, _ = _run_main(monkeypatch, capsys, ["--check"], [(745, 722)], _records())
    assert rc == 2

def test_main_check_exits_0_when_clean(monkeypatch, capsys):
    clean = [{"dependent": 594, "blocker": 211, "state": "closed", "action": "closed-blocker"}]
    rc, _ = _run_main(monkeypatch, capsys, ["--check"], [(594, 211)], clean)
    assert rc == 0

def test_main_apply_wires_only_needs_wiring(monkeypatch, capsys):
    wired = []
    monkeypatch.setattr(spd, "wire", lambda recs: wired.extend((r["dependent"], r["blocker"]) for r in recs) or [])
    rc, _ = _run_main(monkeypatch, capsys, ["--apply"], [(745, 722)], _records())
    assert rc == 0 and wired == [(745, 722)]

def test_main_explicit_report_flag(monkeypatch, capsys):
    # The documented `--report` flag must be accepted and behave like the default.
    rc, out = _run_main(monkeypatch, capsys, ["--report"], [(745, 722)], _records())
    assert rc == 0 and "needs-wiring" in out


# --- markdown-aware parsing (real-world body formats; live-smoke regression) ---

def test_parse_bold_phrase():
    # "This Issue **depends on** #722" — bold markers between phrase and #N
    assert spd.parse_dependencies(745, "This Issue **depends on** #722 rather than folding.") == [(745, 722)]

def test_parse_markdown_link_ref():
    # "Blocked on [Issue #719](url)" — issue ref as a markdown link
    body = "Blocked on [Issue #719](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/719)."
    assert spd.parse_dependencies(725, body) == [(725, 719)]

def test_parse_bold_plus_link():
    # "Blocked by** [Issue #629](url)" — both bold and link
    body = "Blocked by** [Issue #629](https://github.com/x/y/issues/629) — htslib drift."
    assert spd.parse_dependencies(636, body) == [(636, 629)]

def test_parse_short_link_form():
    # "[#42](url)" link form
    assert spd.parse_dependencies(1, "depends on [#42](https://example/42)") == [(1, 42)]

def test_parse_tool_deps_still_ignored():
    # NON-issue deps must stay filtered (no #N after the phrase)
    assert spd.parse_dependencies(566, "requires NetMHCpan-4.0 + netMHCstabpan licenses") == []
    assert spd.parse_dependencies(708, "depends on paid-commercial tools (netMHCpan, netMHCstabpan)") == []
    assert spd.parse_dependencies(601, "blocked on P100: needs bf16 which Pascal lacks") == []

def test_parse_narrative_gap_not_caught_documented_limitation():
    # KNOWN LIMITATION (documented): a multi-word gap between phrase and #N is NOT caught.
    # Strict post-normalization adjacency is deliberate — loosening it would false-match
    # the tool/hardware lines above. These are surfaced by the human review gate instead.
    body = "Depends on the assembled registry ([Issue #732](https://example/732))."
    assert spd.parse_dependencies(737, body) == []
