"""Stubbed-gh tests for scripts/pm/apply_arc_labels.sh (Issue #973, model A').

Model A': membership is label-authoritative (the arc:<slug> label), the manifest
holds only <arc_slug> <phase>, and arc-phase labels are DERIVED. The script
discovers members live via ``gh issue list --label`` and never touches arc:<slug>
labels -- it only reconciles arc-phase.

A fake ``gh`` on PATH makes the reconcile deterministic with no live-board calls:
  - ``gh label list ... --jq ...``            -> $ARCLABELS (defined arc:* labels),
                                                 one per line (unknown-slug guard).
  - ``gh issue list --label <arc> ... --jq``  -> $MEMBERS_<sanitized-arc> (numbers).
  - ``gh issue view N ... --json labels --jq``-> $LABELS_<N> (arc/arc-phase labels).
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
    '  varname="MEMBERS_$san"\n'
    "  printf '%s\\n' ${!varname:-}\n"
    "  exit 0\n"
    "fi\n"
    'if [[ "$1" == "issue" && "$2" == "view" ]]; then\n'
    '  varname="LABELS_$3"\n'
    "  printf '%s\\n' ${!varname:-}\n"
    "  exit 0\n"
    "fi\n"
    'if [[ "$1" == "issue" && "$2" == "edit" ]]; then\n'
    '  echo "$*" >> "$GH_EDIT_LOG"\n'
    '  failvar="FAIL_EDIT_$3"\n'
    '  [[ -n "${!failvar:-}" ]] && exit 1\n'
    "  exit 0\n"
    "fi\n"
    "exit 0\n"
)


def _run(tmp_path, manifest_text, *, arclabels="", members=None, labels=None, check=False, fail_edits=None):
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
