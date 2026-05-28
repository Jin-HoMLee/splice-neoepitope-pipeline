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
    assert ca.is_exempt(["research/glossary.md"]) is True
    assert ca.is_exempt([
        "research/glossary.md", "research/lab_notebook/pm.md",
    ]) is True
    assert ca.is_exempt([
        "research/glossary.md", "workflow/scripts/foo.py",
    ]) is False
    assert ca.is_exempt([]) is False


def test_resolve_roles_returns_per_issue_sets():
    """Per-Issue role sets — NOT flattened/deduped across Issues (#524 fix).

    Multi-role Issues must keep ALL their roles. Position is preserved so callers
    can zip with the original issue list.
    """
    labels = [
        ["role:scientist", "priority:p2"],
        ["role:developer", "role:scientist"],  # multi-role
        ["role:developer"],
        ["priority:p2"],  # no role
    ]
    assert ca.resolve_roles(labels) == [
        {"scientist"},
        {"developer", "scientist"},
        {"developer"},
        set(),
    ]


def test_lab_notebooks_multi_role_any_satisfies():
    """Multi-role Issue, at least one notebook references → no gap (PR #518 repro)."""
    dev_text = "## 2026-05-27\n\n### 14:00 UTC — Editor: Developer\nUnrelated work.\n"
    sci_text = "## 2026-05-27\n\n### 17:39 UTC — Editor: Scientist\nBoltz-2 eval. Refs PR #518.\n"
    notebooks = {"developer": dev_text, "scientist": sci_text}
    assert ca.check_lab_notebooks_for_issue(
        {"developer", "scientist"}, "2026-05-27", 518, notebooks
    ) is None


def test_lab_notebooks_multi_role_none_satisfies():
    """Multi-role Issue, neither notebook references → single combined gap entry."""
    dev_text = "## 2026-05-27\n\n### 14:00 UTC — Editor: Developer\nOther work.\n"
    sci_text = "## 2026-05-27\n\n### 17:39 UTC — Editor: Scientist\nOther work.\n"
    notebooks = {"developer": dev_text, "scientist": sci_text}
    gap = ca.check_lab_notebooks_for_issue(
        {"developer", "scientist"}, "2026-05-27", 518, notebooks
    )
    assert gap is not None
    role_label, desc = gap
    assert "developer" in role_label and "scientist" in role_label
    assert "#518" in desc


def test_lab_notebooks_single_role_no_ref_is_gap():
    """Single-role regression: legacy shape preserved."""
    text = "## 2026-05-27\n\n### 14:00 UTC — Editor: Developer\nOther.\n"
    gap = ca.check_lab_notebooks_for_issue(
        {"developer"}, "2026-05-27", 518, {"developer": text}
    )
    assert gap is not None
    assert gap[0] == "developer"


def test_lab_notebooks_single_role_with_ref_no_gap():
    """Single-role regression: ref present → no gap."""
    text = "## 2026-05-27\n\n### 14:00 UTC — Editor: Developer\nRefs PR #518.\n"
    assert ca.check_lab_notebooks_for_issue(
        {"developer"}, "2026-05-27", 518, {"developer": text}
    ) is None


def test_lab_notebooks_missing_file_is_gap():
    """Missing notebook file (None) reports a per-role gap."""
    gap = ca.check_lab_notebooks_for_issue(
        {"developer"}, "2026-05-27", 518, {"developer": None}
    )
    assert gap is not None
    assert "missing" in gap[1].lower()


def test_lab_notebooks_multi_role_one_missing_other_satisfies():
    """Missing file for one role, the other satisfies → no gap."""
    sci_text = "## 2026-05-27\n\n### 17:39 UTC — Editor: Scientist\nRefs PR #518.\n"
    notebooks = {"developer": None, "scientist": sci_text}
    assert ca.check_lab_notebooks_for_issue(
        {"developer", "scientist"}, "2026-05-27", 518, notebooks
    ) is None


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
