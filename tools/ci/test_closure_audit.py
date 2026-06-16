"""Tests for closure_audit. Lean — focuses on parsing & format edge cases."""

import json

import closure_audit as ca


# --- #607: REPO threading through the gh I/O layer ---


class _FakeCompleted:
    """Stand-in for subprocess.CompletedProcess — only .stdout is read."""

    def __init__(self, stdout):
        self.stdout = stdout


def _capture_gh(monkeypatch, stdout_for=None):
    """Patch subprocess.run to capture each argv.

    `stdout_for` is an optional callable (cmd argv) -> stdout str; when omitted
    every call returns empty-JSON (`"{}"`).
    """
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        return _FakeCompleted(stdout_for(cmd) if stdout_for else "{}")

    monkeypatch.setattr(ca.subprocess, "run", fake_run)
    return calls


def _repo_arg(cmd):
    """The value passed to `--repo` in a gh argv, or None if absent."""
    return cmd[cmd.index("--repo") + 1] if "--repo" in cmd else None


def test_fetch_pr_forwards_repo_when_set(monkeypatch):
    calls = _capture_gh(monkeypatch)
    ca.fetch_pr(99, repo="fork/repo")
    assert _repo_arg(calls[0]) == "fork/repo"


def test_fetch_issue_forwards_repo_when_set(monkeypatch):
    calls = _capture_gh(monkeypatch)
    ca.fetch_issue(42, repo="fork/repo")
    assert _repo_arg(calls[0]) == "fork/repo"


def test_fetch_pr_omits_repo_when_unset(monkeypatch):
    """Default (repo=None) → no --repo flag → gh resolves from git context.

    This is what keeps the post-hoc bot (audit_pr/audit_issue) on git-context
    resolution: they never pass a repo (AC4).
    """
    calls = _capture_gh(monkeypatch)
    ca.fetch_pr(99)
    assert "--repo" not in calls[0]


def test_fetch_issue_omits_repo_when_unset(monkeypatch):
    calls = _capture_gh(monkeypatch)
    ca.fetch_issue(42)
    assert "--repo" not in calls[0]


def test_audit_pr_pre_merge_forwards_repo_to_both_gh_calls(monkeypatch):
    """End-to-end: repo is forwarded to BOTH gh fetches (PR view + issue view).

    The PR fetch returns a closing-issue reference so the per-issue fetch
    actually fires — otherwise `refs == []` short-circuits and only `fetch_pr`
    is exercised (the prior test passed vacuously for the fetch_issue path; PR
    #615 review). The issue carries no role label, so the notebook check is a
    no-op and no _load_notebook FS access happens.
    """
    def stdout_for(cmd):
        if "pr" in cmd:  # `gh pr view ...`
            return json.dumps({"closingIssuesReferences": [{"number": 42}], "files": []})
        return json.dumps({"number": 42, "labels": []})  # `gh issue view ...`

    calls = _capture_gh(monkeypatch, stdout_for)
    ca.audit_pr_pre_merge(99, "2026-06-01", repo="fork/repo")

    # Both a `pr view` and an `issue view` actually happened...
    assert any("pr" in c for c in calls), "fetch_pr was not called"
    assert any("issue" in c for c in calls), "fetch_issue was not called"
    # ...and every gh call forwarded --repo.
    assert all(_repo_arg(c) == "fork/repo" for c in calls)


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


# --- #495: accept the PR number OR any closing-Issue number (`also_accept`) ---


def test_lab_notebook_accepts_pr_number_even_with_also_accept():
    """Only the PR # is referenced (closing Issue # absent) → pass."""
    text = "## 2026-05-26\n\n### 14:00 UTC — Editor: PM\nShipped. Refs PR #494.\n"
    assert ca.check_lab_notebook(text, "2026-05-26", 494, also_accept=[484]) is None


def test_lab_notebook_accepts_closing_issue_number_when_pr_absent():
    """Only a closing Issue # is referenced (PR # absent) → pass via also_accept.

    Repro of the PR #494 false positive: the entry journaled the work against
    Issue #484 but never named the PR number (it was written before the PR
    existed, per the PM/Sci before-PR flow).
    """
    text = "## 2026-05-26\n\n### 14:00 UTC — Editor: PM\nRetired news_log (Issue #484).\n"
    assert ca.check_lab_notebook(text, "2026-05-26", 494, also_accept=[484]) is None


def test_lab_notebook_accepts_either_pr_or_issue():
    """Both the PR # and a closing Issue # referenced → pass."""
    text = "## 2026-05-26\n\n### 14:00 UTC — Editor: PM\nRefs PR #494 / Issue #484.\n"
    assert ca.check_lab_notebook(text, "2026-05-26", 494, also_accept=[484]) is None


def test_lab_notebook_gap_when_neither_pr_nor_any_issue():
    """Neither the PR # nor any closing Issue # referenced → gap; message names all."""
    text = "## 2026-05-26\n\n### 14:00 UTC — Editor: PM\nUnrelated, refs #999.\n"
    gap = ca.check_lab_notebook(text, "2026-05-26", 494, also_accept=[484])
    assert gap is not None
    assert "#494" in gap and "#484" in gap


def test_collect_notebook_gaps_threads_also_accept_for_pr_path():
    """audit_pr path: block names the closing Issue # (not the PR #) → no gap."""
    text = "## 2026-05-26\n\n### 14:00 UTC — Editor: PM\nRetired news_log, Issue #484.\n"
    gaps = ca.collect_notebook_gaps(
        [{"pm"}], "2026-05-26", 494, {"pm": text}, also_accept=[484]
    )
    assert gaps == []


# --- #555: routine-ship opt-out marker skips the lab-notebook check ---


def test_skip_lab_notebook_marker_present_with_value():
    """Canonical marker with a rationale value → skip honored."""
    body = "Routine ship.\n\n<!-- skip-lab-notebook: routine -->\n"
    assert ca.skip_lab_notebook(body) is True


def test_skip_lab_notebook_marker_present_bare():
    """Marker with no `: value` suffix → still honored (presence is what matters)."""
    body = "Routine ship.\n\n<!-- skip-lab-notebook -->\n"
    assert ca.skip_lab_notebook(body) is True


def test_skip_lab_notebook_marker_tolerates_whitespace_and_case():
    """Comment spacing/case shouldn't matter — authors hand-type this."""
    assert ca.skip_lab_notebook("<!--   SKIP-LAB-NOTEBOOK : routine  -->") is True


def test_skip_lab_notebook_marker_absent():
    """No marker → not skipped (check runs normally)."""
    body = "Real feature work that needs a journal entry.\n"
    assert ca.skip_lab_notebook(body) is False


def test_skip_lab_notebook_empty_body():
    """Empty/None body → not skipped."""
    assert ca.skip_lab_notebook("") is False
    assert ca.skip_lab_notebook(None) is False


def test_skip_lab_notebook_does_not_match_unrelated_html_comment():
    """A different HTML comment must not trip the skip."""
    assert ca.skip_lab_notebook("<!-- closure-audit -->\nsome text") is False


def test_ac_deferral_comment_unblocks_unticked():
    # Boxes under a real AC heading (check_ac is AC-section-scoped, #726).
    body = "## Acceptance criteria\n- [x] one\n- [ ] two\n"
    assert ca.check_ac(body, comments=[]) is not None
    assert ca.check_ac(
        body,
        comments=["❎ 'two' deferred to follow-up #99 — out of scope."],
    ) is None


def test_ac_all_ticked_no_gap():
    # AC heading present so this exercises the all-ticked path, not the
    # no-AC-section path (check_ac is AC-section-scoped, #726).
    assert ca.check_ac("## Acceptance criteria\n- [x] one\n- [x] two\n", comments=[]) is None


def test_ac_no_checkboxes_no_gap():
    """Legacy bodies with plain bullets have no boxes to fail on."""
    assert ca.check_ac("- one\n- two\n", comments=[]) is None


def test_check_ac_ignores_unticked_boxes_outside_ac_section():
    """#726/#411: a fully-ticked AC section + an unticked NON-AC checklist
    (`## Flags to evaluate`) must NOT produce an AC gap. check_ac scopes to the
    AC section, matching the pre-merge gate's `unticked_under`."""
    body = (
        "## Acceptance criteria\n"
        "- [x] real ac one\n"
        "- [x] real ac two\n\n"
        "## Flags to evaluate\n"
        "- [x] flag a\n"
        "- [ ] flag b (out of scope)\n"
    )
    assert ca.check_ac(body, comments=[]) is None


def test_check_ac_count_scoped_to_ac_section():
    """A genuine unticked AC box is still flagged, and the count is scoped to
    the AC section (non-AC checklist boxes are excluded from the ratio)."""
    body = (
        "## Acceptance criteria\n- [x] done\n- [ ] not done\n\n"
        "## Tasks\n- [ ] noise one\n- [ ] noise two\n"
    )
    msg = ca.check_ac(body, comments=[])
    assert msg is not None
    assert "1/2" in msg  # 1 unticked of 2 AC boxes — Tasks excluded (not 1/4)


def test_check_ac_no_ac_section_with_stray_boxes_no_gap():
    """#726: with no `## Acceptance criteria` section, unticked boxes elsewhere
    are NOT an AC gap (the merge-time #730 lint surfaces those advisorily).
    Removes the whole-body false-positive that flagged non-AC checklists."""
    body = "## Plan (phased)\n- [ ] P1\n- [x] P2\n"
    assert ca.check_ac(body, comments=[]) is None


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


def test_collect_notebook_gaps_dedupes_identical_role_sets():
    """PR closes 2 Issues sharing the same role → one gap entry, not two."""
    dev_text = "## 2026-05-27\n\n### 14:00 UTC — Editor: Developer\nOther work.\n"
    gaps = ca.collect_notebook_gaps(
        [{"developer"}, {"developer"}], "2026-05-27", 518, {"developer": dev_text}
    )
    assert len(gaps) == 1
    assert gaps[0][0] == "developer"


def test_collect_notebook_gaps_keeps_distinct_role_sets():
    """Different role sets across Issues → distinct gap entries."""
    text = "## 2026-05-27\n\n### 14:00 UTC — Editor: X\nOther.\n"
    gaps = ca.collect_notebook_gaps(
        [{"developer"}, {"scientist"}],
        "2026-05-27", 518,
        {"developer": text, "scientist": text},
    )
    roles_in_gaps = {g[0] for g in gaps}
    assert roles_in_gaps == {"developer", "scientist"}


def test_collect_notebook_gaps_skips_empty_role_sets():
    """Issues with no role: label contribute no gap."""
    assert ca.collect_notebook_gaps(
        [set(), set()], "2026-05-27", 518, {}
    ) == []


# --- #743: not_planned (superseded) closes skip only the lab-notebook check ---


_ANNOTATED_BODY = (
    "- [superseded] one\n- [superseded] two\n\n**Priority rationale:** P2 — x.\n"
)


def _run_audit_issue(monkeypatch, *, state_reason, labels, notebook_text, body=None):
    """Drive audit_issue with a controlled issue + notebook; capture any post.

    Default body has its AC boxes annotated to a disposition (`- [superseded]`)
    so the AC check passes, and carries a Priority-rationale line — isolating the
    lab-notebook check as the only thing that could post a gap. The notebook
    text deliberately omits the Issue ref, so the notebook check WOULD fail if
    it ran. Pass `body` to exercise the AC check (e.g. with a real `- [ ]`).
    Returns the posted comment body, or None if nothing was posted.
    """
    issue = {
        "number": 743,
        "body": _ANNOTATED_BODY if body is None else body,
        "labels": [{"name": lbl} for lbl in labels],
        "comments": [],
        "closedAt": "2026-06-15T14:00:00Z",
        "stateReason": state_reason,
    }
    monkeypatch.setattr(ca, "fetch_issue", lambda n: issue)
    monkeypatch.setattr(ca, "_load_notebook", lambda role: notebook_text)
    posted = {}
    monkeypatch.setattr(
        ca, "post_comment", lambda target, n, body: posted.update(body=body)
    )
    ca.audit_issue(743)
    return posted.get("body")


def test_audit_issue_not_planned_skips_notebook_check(monkeypatch):
    """not_planned close + role label + notebook w/o ref → no notebook gap posted.

    AC boxes are annotated (pass) and a rationale is present, so with the
    notebook check skipped there are no gaps at all → nothing is posted.
    """
    body = _run_audit_issue(
        monkeypatch,
        state_reason="NOT_PLANNED",
        labels=["role:developer"],
        notebook_text="## 2026-06-15\n\n### 14:00 UTC — Editor: Developer\nUnrelated.\n",
    )
    assert body is None


def test_audit_issue_completed_still_flags_missing_notebook(monkeypatch):
    """Regression: COMPLETED close with no notebook ref still posts a notebook gap."""
    body = _run_audit_issue(
        monkeypatch,
        state_reason="COMPLETED",
        labels=["role:developer"],
        notebook_text="## 2026-06-15\n\n### 14:00 UTC — Editor: Developer\nUnrelated.\n",
    )
    assert body is not None and "Lab notebook" in body


def test_audit_issue_null_reason_still_flags_missing_notebook(monkeypatch):
    """Defensive: a null/absent stateReason behaves like COMPLETED (check runs)."""
    body = _run_audit_issue(
        monkeypatch,
        state_reason=None,
        labels=["role:developer"],
        notebook_text="## 2026-06-15\n\n### 14:00 UTC — Editor: Developer\nUnrelated.\n",
    )
    assert body is not None and "Lab notebook" in body


def test_unticked_regex_ignores_annotated_boxes():
    """The AC-annotation convention (#743) relies on _UNTICKED matching ONLY a
    single-space `- [ ]` — annotated dispositions must pass the AC check."""
    for form in ("- [superseded] x\n", "- [n/a] x\n", "- [deferred] x\n"):
        assert ca._UNTICKED.findall(form) == [], f"regex wrongly matched: {form!r}"
    assert ca._UNTICKED.findall("- [ ] x\n"), "regex must still match a real `- [ ]`"


def test_audit_issue_not_planned_still_flags_unticked_ac(monkeypatch):
    """not_planned skips ONLY the notebook check — a genuinely unticked AC
    (left as `- [ ]`, not annotated) is still flagged. Guards the skip's scope."""
    body = _run_audit_issue(
        monkeypatch,
        state_reason="NOT_PLANNED",
        labels=["role:developer"],
        notebook_text="## 2026-06-15\n\n### 14:00 UTC — Editor: Developer\nUnrelated.\n",
        body="## Acceptance criteria\n- [ ] one\n- [superseded] two\n\n**Priority rationale:** P2 — x.\n",
    )
    assert body is not None and "AC checkboxes" in body
    assert "Lab notebook" not in body  # notebook check still skipped


# --- #730: stray-AC-box lint (gating boxes outside an "Acceptance criteria" section) ---


def test_scan_ac_boxes_classifies_ac_section_boxes():
    """Boxes under `## Acceptance criteria` are counted as AC, not stray."""
    body = (
        "## Acceptance criteria\n"
        "- [ ] one\n"
        "- [x] two\n"
    )
    scan = ca.scan_ac_boxes(body)
    assert scan.has_ac_section is True
    assert scan.ac_unticked == 1
    assert scan.ac_total == 2
    assert scan.stray_unticked == 0
    assert scan.stray_headings == []


def test_scan_ac_boxes_flags_stray_boxes_under_non_ac_heading():
    """The #569 shape: gating boxes under `## Plan (phased)`, no AC section."""
    body = (
        "## Context\nsome prose\n\n"
        "## Plan (phased)\n"
        "- [ ] P1 do a thing\n"
        "- [x] P2 done\n"
        "- [ ] P3 another\n"
    )
    scan = ca.scan_ac_boxes(body)
    assert scan.has_ac_section is False
    assert scan.ac_unticked == 0
    assert scan.ac_total == 0
    assert scan.stray_unticked == 2
    assert scan.stray_headings == ["Plan (phased)"]


def test_scan_ac_boxes_stray_boxes_before_any_heading():
    """Boxes before any `## ` heading are stray under the `(top of body)`
    sentinel — guards the cur_heading initialisation (PR #761 review)."""
    body = "- [ ] orphan box\n- [x] done\n\n## Context\nprose\n"
    scan = ca.scan_ac_boxes(body)
    assert scan.has_ac_section is False
    assert scan.stray_unticked == 1
    assert scan.stray_headings == ["(top of body)"]


def test_scan_ac_boxes_no_checkboxes_anywhere():
    body = "## Context\njust prose, no boxes at all\n"
    scan = ca.scan_ac_boxes(body)
    assert scan.has_ac_section is False
    assert scan.ac_total == 0
    assert scan.stray_unticked == 0
    assert scan.stray_headings == []


def test_scan_ac_boxes_counts_both_ac_and_stray():
    """AC section present AND stray boxes elsewhere — both tracked separately."""
    body = (
        "## Plan\n- [ ] stray\n\n"
        "## Acceptance criteria\n- [ ] ac one\n- [x] ac two\n"
    )
    scan = ca.scan_ac_boxes(body)
    assert scan.has_ac_section is True
    assert scan.ac_unticked == 1
    assert scan.ac_total == 2
    assert scan.stray_unticked == 1
    assert scan.stray_headings == ["Plan"]


def test_check_stray_ac_boxes_warns_on_stray_boxes_without_ac_section():
    """No AC section + unticked boxes elsewhere → warning naming count + heading."""
    body = (
        "## Plan (phased)\n"
        "- [ ] P1\n"
        "- [ ] P2\n"
    )
    msg = ca.check_stray_ac_boxes(body)
    assert msg is not None
    assert "2" in msg                     # count of stray unticked boxes
    assert "Plan (phased)" in msg         # the non-AC heading
    assert "Acceptance criteria" in msg   # prompts the canonical-heading convention


def test_check_stray_ac_boxes_silent_when_ac_section_present():
    """An AC section present → the blocking AC gate owns it; lint stays silent
    even if there are unticked boxes under another heading."""
    body = (
        "## Plan\n- [ ] stray\n\n"
        "## Acceptance criteria\n- [ ] real ac\n"
    )
    assert ca.check_stray_ac_boxes(body) is None


def test_check_stray_ac_boxes_silent_when_no_checkboxes():
    body = "## Context\nprose only, nothing to gate\n"
    assert ca.check_stray_ac_boxes(body) is None


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
