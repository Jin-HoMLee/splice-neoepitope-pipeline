"""Tests for closure_audit. Lean — focuses on parsing & format edge cases."""

import closure_audit as ca


def test_lab_notebook_slices_block_correctly():
    """Block lookup must stop at next '## ' header — must not match later dates."""
    text = """# Lab Notebook

## 2026-05-11

### 14:30 UTC — Editor: Developer
Today's entry. Refs PR #100.

## 2026-05-10

### 17:00 UTC — Editor: Developer
Old entry. Refs PR #999.
"""
    # #999 is in the 05-10 block, not the 05-11 block → should report gap
    assert ca.check_lab_notebook(text, "2026-05-11", 999) is not None
    # #100 is in the 05-11 block → no gap
    assert ca.check_lab_notebook(text, "2026-05-11", 100) is None


def test_lab_notebook_no_date_header_is_gap():
    text = "# Lab Notebook\n\n## 2026-05-10\n\n### 17:00 UTC — Editor: Developer\nOld.\n"
    assert ca.check_lab_notebook(text, "2026-05-11", 100) is not None


def test_ac_deferral_comment_unblocks_unticked():
    body = "- [x] one\n- [ ] two\n"
    assert ca.check_ac(body, comments=[]) is not None
    assert ca.check_ac(
        body,
        comments=["❎ 'two' deferred to follow-up #99 — out of scope."],
    ) is None


def test_ac_all_ticked_no_gap():
    assert ca.check_ac("- [x] one\n- [x] two\n", comments=[]) is None


def test_ac_no_checkboxes_no_gap():
    """Legacy bodies with plain bullets have no boxes to fail on."""
    assert ca.check_ac("- one\n- two\n", comments=[]) is None


def test_priority_rationale_present_no_gap():
    body = "Some text.\n\n**Priority rationale:** P2 — example.\n"
    assert ca.check_priority_rationale(body) is None


def test_priority_rationale_header_form_no_gap():
    """Case-insensitive substring tolerates `## Priority rationale` header form."""
    body = "Some text.\n\n## Priority Rationale\n\nP1 — blocks i3.\n"
    assert ca.check_priority_rationale(body) is None


def test_priority_rationale_missing_is_gap():
    body = "Some text without the magic phrase.\n"
    assert ca.check_priority_rationale(body) is not None


def test_exempt_paths_all_or_nothing():
    assert ca.is_exempt(["research/news_log.md"]) is True
    assert ca.is_exempt([
        "research/glossary.md", "research/lab_notebook/pm.md",
    ]) is True
    assert ca.is_exempt([
        "research/news_log.md", "workflow/scripts/foo.py",
    ]) is False
    assert ca.is_exempt([]) is False


def test_resolve_roles_multi_dedupe_alphabetical():
    labels = [
        ["role:scientist", "priority:p2"],
        ["role:developer"],
        ["role:developer"],
        ["priority:p2"],
    ]
    assert ca.resolve_roles(labels) == ["developer", "scientist"]


def test_format_comment_clean_state():
    out = ca.format_comment("PR #1", [], [], [])
    assert "all clear" in out.lower()
    assert ca.COMMENT_MARKER in out


def test_format_comment_lists_only_failing_categories():
    out = ca.format_comment("PR #1", [(2, "2/3 unticked")], [], [])
    assert "AC checkboxes" in out
    assert "Issue #2" in out
    assert "Priority rationale" not in out
    assert "Lab notebook" not in out
