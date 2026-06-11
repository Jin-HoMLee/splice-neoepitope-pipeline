"""Stubbed-gh tests for scripts/pm/apply_arc_labels.sh re-sync (Issue #689).

The script is exercised with a fake ``gh`` on PATH so the re-sync logic is
verified deterministically, with no live-board mutation:
  - ``gh issue view N ... --json labels --jq ...`` → the canned arc/arc-phase
    labels for issue N (from a ``LABELS_<N>`` env var), one per line.
  - ``gh issue edit N ...``                        → the full argv, appended to
    ``$GH_EDIT_LOG``.
"""

import os
import subprocess
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "pm" / "apply_arc_labels.sh"

_FAKE_GH = (
    "#!/usr/bin/env bash\n"
    'if [[ "$1" == "issue" && "$2" == "view" ]]; then\n'
    '  varname="LABELS_$3"\n'
    "  printf '%s\\n' ${!varname:-}\n"
    "  exit 0\n"
    "fi\n"
    'if [[ "$1" == "issue" && "$2" == "edit" ]]; then\n'
    '  echo "$*" >> "$GH_EDIT_LOG"\n'
    "  exit 0\n"
    "fi\n"
    "exit 0\n"
)


def _run(tmp_path, manifest_text, labels):
    gh = tmp_path / "gh"
    gh.write_text(_FAKE_GH)
    gh.chmod(0o755)
    manifest = tmp_path / "arc_taxonomy.tsv"
    manifest.write_text(manifest_text)
    log = tmp_path / "edit.log"
    env = dict(os.environ)
    env["PATH"] = f"{tmp_path}:{env['PATH']}"
    env["GH_EDIT_LOG"] = str(log)
    for n, lbls in labels.items():
        env[f"LABELS_{n}"] = lbls
    subprocess.run(
        ["bash", str(SCRIPT), str(manifest)],
        env=env, check=True, capture_output=True, text=True,
    )
    return log.read_text() if log.exists() else ""


def test_removes_stale_arc_and_phase(tmp_path):
    # #5 mislabeled arc:old + arc-phase:later; manifest says board-governance/active.
    log = _run(tmp_path, "arc:board-governance\tactive\t5\n",
               {5: "arc:old arc-phase:later"})
    assert "--add-label arc:board-governance" in log
    assert "--add-label arc-phase:active" in log
    assert "--remove-label arc:old" in log
    assert "--remove-label arc-phase:later" in log


def test_noop_manifest_removes_nothing(tmp_path):
    # Already-correct issue → re-applies labels, removes nothing (idempotent).
    log = _run(tmp_path, "arc:board-governance\tactive\t5\n",
               {5: "arc:board-governance arc-phase:active"})
    assert "--add-label arc:board-governance" in log
    assert "--remove-label" not in log


def test_phase_flip_removes_only_old_phase(tmp_path):
    # Correct arc, stale phase (next → active): only the phase label churns.
    log = _run(tmp_path, "arc:board-governance\tactive\t5\n",
               {5: "arc:board-governance arc-phase:next"})
    assert "--remove-label arc-phase:next" in log
    assert "--remove-label arc:board-governance" not in log
    assert "--add-label arc-phase:active" in log


def test_retag_removes_only_old_arc(tmp_path):
    # Re-tagged to a different arc, same phase: only the arc label churns.
    log = _run(tmp_path, "arc:board-governance\tactive\t5\n",
               {5: "arc:memory-methodology arc-phase:active"})
    assert "--remove-label arc:memory-methodology" in log
    assert "--remove-label arc-phase:active" not in log
    assert "--add-label arc:board-governance" in log


def test_comment_and_blank_lines_skipped(tmp_path):
    # Only the real manifest row (#5) produces an edit; comment/blank lines don't.
    log = _run(tmp_path, "# comment\n\narc:board-governance\tactive\t5\n", {5: ""})
    edits = [ln for ln in log.splitlines() if ln.strip()]
    assert len(edits) == 1
    assert edits[0].split()[2] == "5"  # issue edit <N> ...
