"""Stubbed-gh tests for scripts/pm/apply_arc_labels.sh (Issue #973, model A').

Model A': membership is label-authoritative (the arc:<slug> label), the manifest
holds only <arc_slug> <phase>, and arc-phase labels are DERIVED. The script
discovers members live via ``gh issue list --label`` and never touches arc:<slug>
labels -- it only reconciles arc-phase.

A fake ``gh`` on PATH makes the reconcile deterministic with no live-board calls:
  - ``gh label list ... --jq ...``            -> $ARCLABELS (defined arc:* labels),
                                                 one per line (unknown-slug guard).
  - ``gh issue list --label <arc> ... --jq``  -> one ``<n>\t<space-joined labels>``
                                                 line per member of $MEMBERS_<arc>,
                                                 labels drawn from $LABELS_<n>
                                                 (labels are returned INLINE now, so
                                                 the script makes no per-member call).
  - ``gh issue edit N ...``                   -> the full argv, appended to $GH_EDIT_LOG.
The arc slug is sanitized to an env-safe name (non-alnum -> ``_``): e.g.
``arc:board-governance`` -> ``MEMBERS_arc_board_governance``.
"""

import os
import subprocess
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "pm" / "apply_arc_labels.sh"

_FAKE_GH = (
    "#!/usr/bin/env bash\n"
    'if [[ "$1" == "label" && "$2" == "list" ]]; then\n'
    "  printf '%s\\n' ${ARCLABELS:-}\n"
    "  exit 0\n"
    "fi\n"
    'if [[ "$1" == "issue" && "$2" == "list" ]]; then\n'
    '  arc=""\n'
    '  while [[ $# -gt 0 ]]; do [[ "$1" == "--label" ]] && arc="$2"; shift; done\n'
    "  san=$(printf '%s' \"$arc\" | sed 's/[^a-zA-Z0-9]/_/g')\n"
    '  mvar="MEMBERS_$san"\n'
    "  for n in ${!mvar:-}; do\n"
    '    lvar="LABELS_$n"\n'
    "    printf '%s\\t%s\\n' \"$n\" \"${!lvar:-}\"\n"
    "  done\n"
    "  exit 0\n"
    "fi\n"
    'if [[ "$1" == "issue" && "$2" == "edit" ]]; then\n'
    '  echo "$*" >> "$GH_EDIT_LOG"\n'
    '  failvar="FAIL_EDIT_$3"\n'
    '  [[ -n "${!failvar:-}" ]] && exit 1\n'
    "  exit 0\n"
    "fi\n"
    # The parent-set query (#1103). Without this branch the fake `gh` fell through
    # to the catch-all `exit 0`, PARENTS stayed empty, and the ENTIRE parent tier
    # was unreachable in every test - the leaf paths passed and proved nothing
    # about the half this feature actually changed. Caught in review of PR #1123.
    'if [[ "$1" == "api" && "$2" == "graphql" ]]; then\n'
    "  for n in ${PARENTS_STUB:-}; do printf '%s\\n' \"$n\"; done\n"
    "  exit 0\n"
    "fi\n"
    "exit 0\n"
)


def _run(tmp_path, manifest_text, *, arclabels="", members=None, labels=None, check=False, fail_edits=None, parents="", extra_env=None):
    gh = tmp_path / "gh"
    gh.write_text(_FAKE_GH)
    gh.chmod(0o755)
    manifest = tmp_path / "arc_taxonomy.tsv"
    manifest.write_text(manifest_text)
    log = tmp_path / "edit.log"
    env = dict(os.environ)
    env["PATH"] = f"{tmp_path}:{env['PATH']}"
    env["GH_EDIT_LOG"] = str(log)
    env["ARCLABELS"] = arclabels
    for arc, nums in (members or {}).items():
        san = "".join(c if c.isalnum() else "_" for c in arc)
        env[f"MEMBERS_{san}"] = nums
    for n, lbls in (labels or {}).items():
        env[f"LABELS_{n}"] = lbls
    for n in (fail_edits or []):
        env[f"FAIL_EDIT_{n}"] = "1"
    env["PARENTS_STUB"] = parents
    env.update(extra_env or {})
    cmd = ["bash", str(SCRIPT), str(manifest)]
    if check:
        cmd.insert(2, "--check")
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True)
    return proc, (log.read_text() if log.exists() else "")


_MANIFEST = "arc:board-governance\tactive\tBoard governance\n"


def test_phase_drift_fixed_arc_label_untouched(tmp_path):
    # #5 in board-governance carries a stale arc-phase:later -> reconcile to active.
    proc, log = _run(
        tmp_path, _MANIFEST,
        arclabels="arc:board-governance",
        members={"arc:board-governance": "5"},
        labels={5: "arc:board-governance arc-phase:later"},
    )
    assert proc.returncode == 0, proc.stderr
    assert "--add-label arc-phase:active" in log
    assert "--remove-label arc-phase:later" in log
    # membership is label-authoritative: the arc:<slug> label is never touched.
    assert "--remove-label arc:board-governance" not in log
    assert "--add-label arc:board-governance" not in log


def test_clean_member_no_edit(tmp_path):
    # Already-correct phase -> no board mutation at all (idempotent).
    proc, log = _run(
        tmp_path, _MANIFEST,
        arclabels="arc:board-governance",
        members={"arc:board-governance": "5"},
        labels={5: "arc:board-governance arc-phase:active"},
    )
    assert proc.returncode == 0, proc.stderr
    assert log == ""


def test_missing_phase_added(tmp_path):
    # Member carries the arc label but NO arc-phase yet -> phase added, nothing removed.
    proc, log = _run(
        tmp_path, _MANIFEST,
        arclabels="arc:board-governance",
        members={"arc:board-governance": "5"},
        labels={5: "arc:board-governance"},
    )
    assert proc.returncode == 0, proc.stderr
    assert "--add-label arc-phase:active" in log
    assert "--remove-label" not in log


def test_multi_arc_skipped_not_mutated(tmp_path):
    # >1 arc label -> ambiguous phase: skipped, never auto-set.
    proc, log = _run(
        tmp_path, _MANIFEST,
        arclabels="arc:board-governance arc:scoring-tcr-pmhc",
        members={"arc:board-governance": "5"},
        labels={5: "arc:board-governance arc:scoring-tcr-pmhc arc-phase:active"},
    )
    assert proc.returncode == 0, proc.stderr
    assert log == ""
    assert "multi-arc" in proc.stdout


def test_check_mode_exits_2_on_drift_and_mutates_nothing(tmp_path):
    proc, log = _run(
        tmp_path, _MANIFEST,
        arclabels="arc:board-governance",
        members={"arc:board-governance": "5"},
        labels={5: "arc:board-governance arc-phase:later"},
        check=True,
    )
    assert proc.returncode == 2
    assert "DRIFT phase" in proc.stdout
    assert log == ""  # read-only


def test_check_mode_clean_exits_0(tmp_path):
    proc, log = _run(
        tmp_path, _MANIFEST,
        arclabels="arc:board-governance",
        members={"arc:board-governance": "5"},
        labels={5: "arc:board-governance arc-phase:active"},
        check=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert "clean" in proc.stdout
    assert log == ""


def test_unknown_slug_flagged(tmp_path):
    # A defined arc:* label absent from the manifest roster is flagged as drift.
    proc, log = _run(
        tmp_path, _MANIFEST,
        arclabels="arc:board-governance arc:retired-arc",
        members={"arc:board-governance": "5"},
        labels={5: "arc:board-governance arc-phase:active"},
        check=True,
    )
    assert proc.returncode == 2
    assert "unknown-slug" in proc.stdout
    assert "arc:retired-arc" in proc.stdout


def test_transient_edit_failure_tolerated_batch_continues(tmp_path):
    # A 5xx on one issue must not abort the batch: #5 fails, #6 still gets fixed,
    # and the run exits 3 to signal a retry.
    proc, log = _run(
        tmp_path, _MANIFEST,
        arclabels="arc:board-governance",
        members={"arc:board-governance": "5 6"},
        labels={5: "arc:board-governance", 6: "arc:board-governance"},
        fail_edits=[5],
    )
    assert proc.returncode == 3
    assert "WARN: edit #5 failed" in proc.stdout
    # #6 was still attempted after #5 failed (batch did not abort).
    assert "issue edit 6" in log


def test_comment_and_blank_manifest_lines_skipped(tmp_path):
    # Comment/blank lines produce no membership queries / edits.
    proc, log = _run(
        tmp_path, "# comment\n\n" + _MANIFEST,
        arclabels="arc:board-governance",
        members={"arc:board-governance": "5"},
        labels={5: "arc:board-governance arc-phase:active"},
    )
    assert proc.returncode == 0, proc.stderr
    assert log == ""


def test_invalid_phase_in_manifest_fails_fast(tmp_path):
    # A phase typo would make every `gh issue edit --add-label arc-phase:<typo>`
    # fail permanently -> reject it up front, not as a spurious transient retry.
    proc, log = _run(
        tmp_path, "arc:board-governance\tactiv\tBoard governance\n",
        arclabels="arc:board-governance",
        members={"arc:board-governance": "5"},
        labels={5: "arc:board-governance"},
    )
    assert proc.returncode == 2
    assert "invalid phase" in (proc.stdout + proc.stderr)
    assert log == ""  # bailed before any edit


_FOUR_ACTIVE = (
    "arc:a\tactive\tA\n"
    "arc:b\tactive\tB\n"
    "arc:c\tactive\tC\n"
    "arc:d\tactive\tD\n"
)


def test_more_than_three_active_is_advisory_not_drift(tmp_path):
    # Issue #1102: the focus slate is a ~3 GUIDELINE, not a cap. Exceeding it prints
    # an advisory NOTE and must not affect the exit code -- `active` was measured not
    # to predict throughput, so a 4th active arc is a legitimate state, not drift.
    proc, log = _run(
        tmp_path, _FOUR_ACTIVE,
        arclabels="arc:a arc:b arc:c arc:d",
        members={},  # no members needed; the guideline is a manifest-level check
        check=True,
    )
    assert proc.returncode == 0
    assert "NOTE cap" in proc.stdout
    assert "4 arcs marked 'active'" in proc.stdout
    assert "DRIFT cap" not in proc.stdout


def test_active_slate_guideline_is_tunable(tmp_path):
    # Raising the guideline suppresses the NOTE entirely; exit code is unchanged.
    proc, log = _run(
        tmp_path, _FOUR_ACTIVE,
        arclabels="arc:a arc:b arc:c arc:d",
        members={},
        check=True,
        extra_env={"ACTIVE_SLATE_GUIDELINE": "99"},
    )
    assert proc.returncode == 0
    assert "NOTE cap" not in proc.stdout


def test_non_numeric_guideline_falls_back_to_default(tmp_path):
    # A garbage override must not make an *advisory* knob misbehave: fall back to 3,
    # warn, and still emit the NOTE for the 4-arc slate.
    proc, log = _run(
        tmp_path, _FOUR_ACTIVE,
        arclabels="arc:a arc:b arc:c arc:d",
        members={},
        check=True,
        extra_env={"ACTIVE_SLATE_GUIDELINE": "2x"},
    )
    assert proc.returncode == 0
    assert "not a non-negative integer" in proc.stderr
    assert "NOTE cap" in proc.stdout
    assert "value too great for base" not in proc.stderr  # no raw bash arith error


def test_over_guideline_does_not_mask_real_drift(tmp_path):
    # The advisory NOTE must coexist with genuine drift without swallowing its exit 2.
    # #5 is a member of arc:a but carries the wrong arc-phase, which is real drift.
    proc, log = _run(
        tmp_path, _FOUR_ACTIVE,
        arclabels="arc:a arc:b arc:c arc:d",
        members={"arc:a": "5"},
        labels={5: "arc:a arc-phase:later"},
        check=True,
    )
    assert proc.returncode == 2
    assert "NOTE cap" in proc.stdout
    assert "DRIFT phase" in proc.stdout


def test_multi_arc_reported_once_across_arcs(tmp_path):
    # #5 is a member of two arcs; it must be reported exactly once, not per-arc.
    manifest = "arc:board-governance\tactive\tBG\narc:scoring-tcr-pmhc\tactive\tScoring\n"
    proc, log = _run(
        tmp_path, manifest,
        arclabels="arc:board-governance arc:scoring-tcr-pmhc",
        members={"arc:board-governance": "5", "arc:scoring-tcr-pmhc": "5"},
        labels={5: "arc:board-governance arc:scoring-tcr-pmhc arc-phase:active"},
        check=True,
    )
    assert proc.returncode == 2
    assert proc.stdout.count("multi-arc: #5") == 1


def test_member_query_uses_limit_1000(tmp_path):
    # Guard against gh issue list's silent default --limit 30 truncation.
    assert "--limit 1000" in SCRIPT.read_text()


# --- #1103: the parent tier -------------------------------------------------
#
# Multi-arc is LEGAL on a parent (subIssuesSummary.total > 0) and a parent carries
# NO arc-phase (it is parked in `Epic` and never pulled, so a focus marker is
# meaningless). Both remain drift on a LEAF.
#
# Every test below was unreachable before PR #1123's review: the fake `gh` had no
# `api graphql` branch, so PARENTS was always empty and the parent code never ran.

_TWO_ARC_MANIFEST = (
    "arc:board-governance\tactive\tBoard governance\n"
    "arc:memory-methodology\tlater\tMemory methodology\n"
)


def test_multi_arc_parent_is_not_drift(tmp_path):
    """The #1036 case: a parent spanning two arcs is legal, not multi-arc drift."""
    proc, log = _run(
        tmp_path, _TWO_ARC_MANIFEST,
        arclabels="arc:board-governance arc:memory-methodology",
        members={"arc:board-governance": "1036", "arc:memory-methodology": "1036"},
        labels={"1036": "arc:board-governance arc:memory-methodology"},
        parents="1036",
        check=True,
    )
    assert "multi-arc" not in proc.stdout, proc.stdout
    assert proc.returncode == 0, f"guard should be CLEAN\n{proc.stdout}"
    assert log == "", "a clean parent must not be edited"


def test_multi_arc_LEAF_is_still_drift(tmp_path):
    """Control: the same two arcs on a LEAF are still drift. The relaxation is
    scoped to parents - if this ever passes, the rule has been over-relaxed."""
    proc, _ = _run(
        tmp_path, _TWO_ARC_MANIFEST,
        arclabels="arc:board-governance arc:memory-methodology",
        members={"arc:board-governance": "1036", "arc:memory-methodology": "1036"},
        labels={"1036": "arc:board-governance arc:memory-methodology"},
        parents="",  # NOT a parent
        check=True,
    )
    assert "multi-arc" in proc.stdout, proc.stdout
    assert proc.returncode == 2


def test_parent_carrying_a_phase_is_stripped(tmp_path):
    """A parent with an arc-phase is drift; apply-mode removes it."""
    proc, log = _run(
        tmp_path, _MANIFEST,
        arclabels="arc:board-governance",
        members={"arc:board-governance": "680"},
        labels={"680": "arc:board-governance arc-phase:active"},
        parents="680",
    )
    assert "strip arc-phase (parent" in proc.stdout, proc.stdout
    assert "--remove-label arc-phase:active" in log, log
    assert "--add-label" not in log, "a parent must never be GIVEN a phase"


def test_parent_carrying_a_phase_is_reported_in_check_mode(tmp_path):
    """--check reports it read-only, exits 2, and edits nothing."""
    proc, log = _run(
        tmp_path, _MANIFEST,
        arclabels="arc:board-governance",
        members={"arc:board-governance": "680"},
        labels={"680": "arc:board-governance arc-phase:active"},
        parents="680",
        check=True,
    )
    assert "DRIFT parent-phase: #680" in proc.stdout, proc.stdout
    assert proc.returncode == 2
    assert log == "", "--check must not mutate"


def test_multi_arc_parent_with_a_phase_is_reported_once(tmp_path):
    """seen_parent dedup: a parent found under 2 arcs is reported/stripped once."""
    proc, log = _run(
        tmp_path, _TWO_ARC_MANIFEST,
        arclabels="arc:board-governance arc:memory-methodology",
        members={"arc:board-governance": "1036", "arc:memory-methodology": "1036"},
        labels={"1036": "arc:board-governance arc:memory-methodology arc-phase:active"},
        parents="1036",
    )
    assert proc.stdout.count("strip arc-phase (parent") == 1, proc.stdout
    assert log.count("--remove-label") == 1, log


def test_clean_parent_is_a_no_op(tmp_path):
    """A parent with arcs and no phase: no drift, no edit, exit 0."""
    proc, log = _run(
        tmp_path, _MANIFEST,
        arclabels="arc:board-governance",
        members={"arc:board-governance": "680"},
        labels={"680": "arc:board-governance"},
        parents="680",
        check=True,
    )
    assert proc.returncode == 0, proc.stdout
    assert "DRIFT" not in proc.stdout, proc.stdout
    assert log == ""


def test_leaf_phase_reconcile_still_works_alongside_a_parent(tmp_path):
    """A parent and a leaf in the same arc: leaf gets its phase, parent gets none."""
    proc, log = _run(
        tmp_path, _MANIFEST,
        arclabels="arc:board-governance",
        members={"arc:board-governance": "680 999"},
        labels={"680": "arc:board-governance arc-phase:active", "999": "arc:board-governance"},
        parents="680",
    )
    assert "--add-label arc-phase:active" in log and "999" in log, log      # leaf phased
    assert "--remove-label arc-phase:active" in log and "680" in log, log   # parent stripped
