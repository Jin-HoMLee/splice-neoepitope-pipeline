"""Tests for closure_audit. Lean — focuses on parsing & format edge cases."""

import json
import re

import pytest

import closure_audit as ca

REPO_ROOT = ca.REPO_ROOT


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


def test_lab_notebook_accepts_prior_day_block_within_window():
    """A block from a prior day, inside the window, satisfies the gate (#1092).

    This test previously asserted the OPPOSITE - that yesterday's block must not
    satisfy today's gate. That assertion encoded the #1092 bug, so it is inverted
    here rather than deleted, and renamed to describe the contract it now checks.
    """
    text = """# Lab Notebook

## 2026-05-11

### 14:30 UTC — Editor: Developer
Today's entry. Refs PR #100.

## 2026-05-10

### 17:00 UTC — Editor: Developer
Old entry. Refs PR #999.
"""
    # #999 sits in the 05-10 block, one day before the gate date. Under #1092 that
    # is ACCEPTED: the entry exists for that unit of work and lies in the lookback
    # window. (Pre-#1092 this asserted a gap - that assertion encoded the bug.)
    assert ca.check_lab_notebook(text, "2026-05-11", 999) is None
    # #100 is in the 05-11 block -> no gap
    assert ca.check_lab_notebook(text, "2026-05-11", 100) is None


def test_lab_notebook_no_date_header_is_gap():
    """No block references #100 in any scanned block -> gap."""
    text = "# Lab Notebook\n\n## 2026-05-10\n\n### 17:00 UTC - Editor: Developer\nOld.\n"
    assert ca.check_lab_notebook(text, "2026-05-11", 100) is not None


# --- #1092: the gate must not key on the *merge* date ------------------------
#
# Entry timing is commit -> push -> PR -> review -> entry -> merge, so the entry
# is written after review and before merge. Whenever those straddle a UTC
# midnight (the normal case, since merge is a human act and PRs sit at the gate)
# the entry is dated *yesterday* and a today-keyed gate false-blocks a PR whose
# author followed the rule exactly.


def test_lab_notebook_accepts_entry_written_before_the_merge_date():
    """The straddling-midnight case: entry dated earlier, merge lands today."""
    text = (
        "## 2026-07-09\n\n"
        "### 16:00 UTC - Editor: Scientist\n"
        "pVACtools spike shipped. Refs PR #1093.\n"
    )
    assert ca.check_lab_notebook(text, "2026-07-11", 1093) is None


def test_lab_notebook_accepts_prior_day_entry_via_also_accept():
    """A prior-day entry naming only the closing Issue still passes."""
    text = "## 2026-07-09\n\n### 16:00 UTC - Editor: Scientist\nSpike (Issue #1048).\n"
    assert ca.check_lab_notebook(text, "2026-07-11", 1093, also_accept=[1048]) is None


def test_lab_notebook_today_block_for_other_work_does_not_decide():
    """A today-dated block for *different* work must not decide the outcome.

    The exact #1093 shape: scientist.md had a `## 2026-07-11` block (for #1101)
    while #1093's real entry sat in the 07-09 block. Neither the presence of a
    today block nor its wrong reference may short-circuit the scan - only a block
    that actually references this unit of work passes.
    """
    text = (
        "## 2026-07-11\n\n### 09:00 UTC - Editor: Scientist\nRegistry schema. Refs PR #1101.\n\n"
        "## 2026-07-09\n\n### 16:00 UTC - Editor: Scientist\nSpike. Refs PR #1093.\n"
    )
    # #1093 is found in the older block -> pass
    assert ca.check_lab_notebook(text, "2026-07-11", 1093) is None
    # #1234 is referenced nowhere -> gap, despite a today block existing
    assert ca.check_lab_notebook(text, "2026-07-11", 1234) is not None


def test_lab_notebook_gap_when_wrong_reference_in_every_scanned_block():
    """An entry for a different unit of work still fails, in every block scanned."""
    text = (
        "## 2026-07-11\n\n### 09:00 UTC - Editor: PM\nRefs PR #111.\n\n"
        "## 2026-07-10\n\n### 09:00 UTC - Editor: PM\nRefs PR #222.\n\n"
        "## 2026-07-09\n\n### 09:00 UTC - Editor: PM\nRefs PR #333.\n"
    )
    gap = ca.check_lab_notebook(text, "2026-07-11", 999, also_accept=[888])
    assert gap is not None
    assert "#999" in gap and "#888" in gap


def test_lab_notebook_gap_when_no_entry_at_all():
    """A PR with no lab-notebook entry whatsoever still fails, unchanged."""
    assert ca.check_lab_notebook("# Lab Notebook\n", "2026-07-11", 1093) is not None


def test_lab_notebook_entry_older_than_the_window_is_a_gap():
    """The window is bounded: a months-old block must not satisfy the gate."""
    text = "## 2026-01-02\n\n### 16:00 UTC - Editor: Scientist\nAncient. Refs PR #1093.\n"
    assert ca.check_lab_notebook(text, "2026-07-11", 1093) is not None


def test_lab_notebook_matched_block_still_needs_a_time_subsection():
    """The '### HH:MM UTC - Editor: ...' requirement survives into the window scan."""
    text = "## 2026-07-10\n\nRefs PR #1093 but with no time sub-section.\n"
    gap = ca.check_lab_notebook(text, "2026-07-11", 1093)
    assert gap is not None
    assert "###" in gap


def test_lab_notebook_ignores_blocks_dated_after_the_gate_date():
    """A future-dated block is not a valid entry for this merge."""
    text = "## 2026-07-20\n\n### 16:00 UTC - Editor: PM\nRefs PR #1093.\n"
    assert ca.check_lab_notebook(text, "2026-07-11", 1093) is not None


# --- Both LIVE header conventions must parse ---------------------------------
#
# pm.md and scientist.md use a bare `## YYYY-MM-DD`. developer.md suffixes a
# description: `## 2026-07-09 - ship the cwd-drift guard pair ([PR #1088] ...)`.
# The old substring check (`f"## {date}" in text`) tolerated both. An anchored
# regex demanding a bare date silently skips every developer block, false-blocking
# every developer PR forever - strictly WORSE than the bug being fixed. Caught in
# review of PR #1121; every fixture above uses a bare date, which is exactly why
# it slipped. These lock both live shapes in.


def test_lab_notebook_accepts_description_suffixed_header():
    """The real developer.md shape: `## <date> - <description>`."""
    text = (
        "## 2026-07-09 - ship the cwd-drift guard pair "
        "([PR #1088](https://example.com) closes [Issue #1053](https://example.com))\n\n"
        "### 14:00 UTC - Editor: Developer\nShipped.\n"
    )
    assert ca.check_lab_notebook(text, "2026-07-11", 1088, also_accept=[1053]) is None


def test_lab_notebook_suffixed_header_still_gaps_on_wrong_ref():
    """A suffixed header must not become a blanket pass - the reference key still binds."""
    text = (
        "## 2026-07-09 - ship the cwd-drift guard pair ([PR #1088](https://example.com))\n\n"
        "### 14:00 UTC - Editor: Developer\nShipped.\n"
    )
    assert ca.check_lab_notebook(text, "2026-07-11", 9999, also_accept=[8888]) is not None


def test_dated_blocks_parses_both_live_conventions():
    """Both shapes are found, and a suffixed block does not swallow the next one."""
    text = (
        "## 2026-07-10 - suffixed header ([PR #1](https://example.com))\n\nbody A\n\n"
        "## 2026-07-09\n\nbody B\n"
    )
    blocks = ca._dated_blocks(text)
    assert [d.isoformat() for d, _ in blocks] == ["2026-07-10", "2026-07-09"]
    assert "body A" in blocks[0][1] and "body B" not in blocks[0][1]
    assert "body B" in blocks[1][1]


def test_dated_blocks_ignores_non_date_h2():
    """`## Acceptance criteria` and friends are not date blocks."""
    assert ca._dated_blocks("## Acceptance criteria\n\n- [ ] x\n") == []


def test_lab_notebook_reference_match_is_digit_bounded():
    """`#112` must not match `#1121` (decimal-prefix collision).

    Pre-existing in the old substring check, but amplified by the #1092 window:
    the scan went from one exact-date block to up to a week of them, so more
    blocks are exposed to a spurious prefix hit. A gate that passes because a
    *different* Issue's number happens to start with the gated one is not gating.
    """
    text = "## 2026-07-10\n\n### 09:00 UTC - Editor: PM\nRefs PR #1121.\n"
    assert ca.check_lab_notebook(text, "2026-07-11", 112) is not None   # prefix
    assert ca.check_lab_notebook(text, "2026-07-11", 1121) is None      # exact
    # also_accept must be boundary-bound too, not just `number`
    assert ca.check_lab_notebook(text, "2026-07-11", 999, also_accept=[112]) is not None


def test_lab_notebook_reference_match_allows_trailing_punctuation():
    """A boundary is a non-digit, so `#1121.` / `#1121)` / `#1121,` still match."""
    for tail in [".", ")", ",", "]", " "]:
        text = f"## 2026-07-10\n\n### 09:00 UTC - Editor: PM\nRefs PR #1121{tail}\n"
        assert ca.check_lab_notebook(text, "2026-07-11", 1121) is None, tail


def test_dated_blocks_skips_a_malformed_date_rather_than_raising():
    """A notebook typo must not crash the merge gate."""
    text = "## 2026-13-45\n\nbad month/day\n\n## 2026-07-10\n\n### 09:00 UTC - Editor: PM\nRefs #1.\n"
    parsed = [d.isoformat() for d, _ in ca._dated_blocks(text)]
    assert parsed == ["2026-07-10"]


@pytest.mark.parametrize("role", ["developer", "pm", "scientist"])
def test_every_dated_header_in_the_real_notebooks_parses(role):
    """Regression lock against the live notebooks, not curated fixtures.

    Asserts EVERY `## <date>...` header in each real notebook is parsed - not
    merely that *some* block was found. That weaker form passes vacuously on
    developer.md, whose older entries are bare-date while its recent ones carry a
    description suffix; only the suffixed ones regressed in PR #1121. A guard that
    can be satisfied by the un-regressed half of a file is not a guard.
    """
    nb = REPO_ROOT / "research" / "lab_notebook" / f"{role}.md"
    if not nb.exists():
        pytest.skip(f"{role}.md not present")
    text = nb.read_text()

    # Every line that a human would read as a dated entry header.
    written = re.findall(r"^## (\d{4}-\d{2}-\d{2})", text, re.MULTILINE)
    parsed = [d.isoformat() for d, _ in ca._dated_blocks(text)]
    assert written, f"no dated headers found in {role}.md at all"
    missed = [d for d in written if d not in parsed]
    assert not missed, f"{role}.md: _dated_blocks missed {len(missed)} header(s): {missed[:3]}"


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


# --- #748: Memory Manager role is lab-notebook-exempt (its record is the
# personas-repo git log, not research/lab_notebook/) ---


def test_collect_notebook_gaps_exempts_pure_memory_manager():
    """A pure role:memory_manager Issue needs no project-repo notebook entry."""
    assert ca.collect_notebook_gaps(
        [{"memory_manager"}], "2026-06-30", 748, {}
    ) == []


def test_collect_notebook_gaps_mixed_mm_still_requires_other_role():
    """dev+MM Issue with no dev entry → still a developer gap. The exemption
    strips only the MM role; it does not exempt the whole Issue."""
    gaps = ca.collect_notebook_gaps(
        [{"developer", "memory_manager"}], "2026-06-30", 748, {"developer": None}
    )
    assert len(gaps) == 1
    assert "developer" in gaps[0][0]
    assert "memory_manager" not in gaps[0][0]


def test_collect_notebook_gaps_mixed_mm_satisfied_by_other_role_entry():
    """dev+MM Issue with a developer entry → no gap (MM contributes nothing)."""
    dev_text = "## 2026-06-30\n\n### 12:00 UTC — Editor: Developer\nClosed #748.\n"
    assert ca.collect_notebook_gaps(
        [{"developer", "memory_manager"}], "2026-06-30", 748, {"developer": dev_text}
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


# --- #665: cross-repo closing forward-link AC coverage ---


THIS_REPO = "Jin-HoMLee/claude-personas-splice-neoepitope-pipeline"
PROJ = "Jin-HoMLee/splice-neoepitope-pipeline"


class TestParseCrossRepoAcTargets:
    def test_shorthand_owner_repo_hash(self):
        body = f"Closes {PROJ}#665 — cross-repo coverage fix."
        assert ca.parse_cross_repo_ac_targets(body, THIS_REPO) == [(PROJ, 665)]

    def test_full_issue_url(self):
        body = f"Fixes https://github.com/{PROJ}/issues/409 in this PR."
        assert ca.parse_cross_repo_ac_targets(body, THIS_REPO) == [(PROJ, 409)]

    def test_link_form_keyword_in_parens(self):
        # project convention: [Issue #N](url) (keyword) — keyword AFTER the link
        body = f"[Issue #665](https://github.com/{PROJ}/issues/665) (closes)"
        assert ca.parse_cross_repo_ac_targets(body, THIS_REPO) == [(PROJ, 665)]

    def test_excludes_same_repo_shorthand(self):
        # a forward-link to the PR's OWN repo is covered by native references
        body = f"Closes {THIS_REPO}#12 and {PROJ}#665"
        assert ca.parse_cross_repo_ac_targets(body, THIS_REPO) == [(PROJ, 665)]

    def test_same_repo_case_insensitive_exclusion(self):
        body = f"Closes {THIS_REPO.upper()}#12"
        assert ca.parse_cross_repo_ac_targets(body, THIS_REPO) == []

    def test_no_closing_keyword_on_line_ignored(self):
        # a bare cross-repo mention without a closing keyword is not a close intent
        body = f"See {PROJ}#665 for context (related, not closing)."
        assert ca.parse_cross_repo_ac_targets(body, THIS_REPO) == []

    def test_dedupes_repeated_target(self):
        body = (
            f"Closes {PROJ}#665.\n"
            f"Resolves https://github.com/{PROJ}/issues/665 too."
        )
        assert ca.parse_cross_repo_ac_targets(body, THIS_REPO) == [(PROJ, 665)]

    def test_multiple_distinct_targets_order_preserved(self):
        body = f"Closes {PROJ}#665\nFixes {PROJ}#409"
        assert ca.parse_cross_repo_ac_targets(body, THIS_REPO) == [
            (PROJ, 665),
            (PROJ, 409),
        ]

    def test_multiple_targets_one_closing_line(self):
        # two distinct cross-repo refs sharing one closing-keyword line are both found
        body = f"Closes {PROJ}#665 and {PROJ}#409 in a single line"
        assert ca.parse_cross_repo_ac_targets(body, THIS_REPO) == [
            (PROJ, 665),
            (PROJ, 409),
        ]

    def test_bare_same_repo_hash_not_matched(self):
        # `#665` with no owner/repo prefix is same-repo (native) — never matched
        body = "Closes #665"
        assert ca.parse_cross_repo_ac_targets(body, THIS_REPO) == []

    def test_empty_body(self):
        assert ca.parse_cross_repo_ac_targets("", THIS_REPO) == []

    def test_none_this_repo_includes_all_crossrepo_forms(self):
        # when the PR's own repo is unknown, don't drop anything (conservative)
        body = f"Closes {PROJ}#665"
        assert ca.parse_cross_repo_ac_targets(body, None) == [(PROJ, 665)]


class TestCollectCrossRepoAcGaps:
    def _io(self, monkeypatch, pr_body, issues_by_key):
        monkeypatch.setattr(
            ca, "fetch_pr", lambda n, repo=None: {"body": pr_body}
        )
        monkeypatch.setattr(
            ca, "fetch_issue",
            lambda n, repo=None: issues_by_key[(repo, n)],
        )

    def test_unticked_cross_repo_ac_is_a_gap(self, monkeypatch):
        pr_body = f"Closes {PROJ}#665"
        issue = {
            "number": 665,
            "body": "## Acceptance criteria\n- [ ] not done\n",
            "comments": [],
        }
        self._io(monkeypatch, pr_body, {(PROJ, 665): issue})
        gaps = ca.collect_cross_repo_ac_gaps(99, repo=THIS_REPO)
        assert len(gaps) == 1
        assert gaps[0][0] == f"{PROJ}#665"

    def test_ticked_cross_repo_ac_no_gap(self, monkeypatch):
        pr_body = f"Closes {PROJ}#665"
        issue = {
            "number": 665,
            "body": "## Acceptance criteria\n- [x] done\n",
            "comments": [],
        }
        self._io(monkeypatch, pr_body, {(PROJ, 665): issue})
        assert ca.collect_cross_repo_ac_gaps(99, repo=THIS_REPO) == []

    def test_fetches_target_from_its_own_repo(self, monkeypatch):
        # the cross-repo issue must be fetched with --repo <target>, not the PR repo
        seen = {}
        monkeypatch.setattr(ca, "fetch_pr", lambda n, repo=None: {"body": f"Closes {PROJ}#665"})

        def fake_issue(n, repo=None):
            seen["repo"] = repo
            return {"number": n, "body": "## Acceptance criteria\n- [x] ok\n", "comments": []}

        monkeypatch.setattr(ca, "fetch_issue", fake_issue)
        ca.collect_cross_repo_ac_gaps(99, repo=THIS_REPO)
        assert seen["repo"] == PROJ

    def test_cross_repo_issue_without_ac_section_no_gap(self, monkeypatch):
        # a target with no `## Acceptance criteria` section has 0 unticked AC boxes,
        # so it passes cleanly — consistent with same-repo check_ac behavior
        pr_body = f"Closes {PROJ}#665"
        issue = {"number": 665, "body": "Plain body, no AC heading.\n", "comments": []}
        self._io(monkeypatch, pr_body, {(PROJ, 665): issue})
        assert ca.collect_cross_repo_ac_gaps(99, repo=THIS_REPO) == []

    def test_no_targets_no_fetch_no_gap(self, monkeypatch):
        monkeypatch.setattr(ca, "fetch_pr", lambda n, repo=None: {"body": "no closers here"})

        def boom(*a, **k):
            raise AssertionError("fetch_issue must not be called with no targets")

        monkeypatch.setattr(ca, "fetch_issue", boom)
        assert ca.collect_cross_repo_ac_gaps(99, repo=THIS_REPO) == []
