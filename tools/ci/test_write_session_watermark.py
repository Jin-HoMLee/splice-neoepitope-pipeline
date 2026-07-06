"""Tests for `.agents/hooks/write_session_watermark.py` (Issue #820).

The hook is a Stop hook: it records a deterministic UTC high-water mark so the
morning routine can anchor its "activity since last check" window on the last
turn-end instead of inferring it from episode filenames or a reflexive
"yesterday".

Two layers:
- subprocess tests mirror how the harness invokes the hook (Stop JSON on stdin,
  always exit 0, never emit a stop decision);
- direct-import tests exercise the pure `write_watermark` writer with an
  injected timestamp for determinism.
"""

import importlib.util
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pytest

HOOK = Path(__file__).parent.parent.parent / ".agents" / "hooks" / "write_session_watermark.py"
MARKER_RELPATH = Path(".agents") / "last_session_marker.json"
ROUTINE_MARKER_RELPATH = Path(".agents") / "last_routine_marker.json"
ISO_Z = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")


def _load_module():
    spec = importlib.util.spec_from_file_location("write_session_watermark", HOOK)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run(stdin_payload: str, project_dir: Path):
    proc = subprocess.run(
        [sys.executable, str(HOOK)],
        input=stdin_payload,
        capture_output=True,
        text=True,
        timeout=5,
        env={"CLAUDE_PROJECT_DIR": str(project_dir), "PATH": ""},
    )
    return proc.returncode, proc.stdout, proc.stderr


def _stop_payload(**extra) -> str:
    base = {"hook_event_name": "Stop", "session_id": "abc", "stop_hook_active": False}
    base.update(extra)
    return json.dumps(base)


def _read_marker(project_dir: Path) -> dict:
    return json.loads((project_dir / MARKER_RELPATH).read_text())


class TestSubprocessContract:
    def test_exit_zero_and_no_stdout(self, tmp_path):
        (tmp_path / ".agents").mkdir()
        rc, out, _ = _run(_stop_payload(), tmp_path)
        # A Stop hook must never block stopping: exit 0, emit no decision.
        assert rc == 0
        assert out.strip() == ""

    def test_marker_written_with_valid_shape(self, tmp_path):
        (tmp_path / ".agents").mkdir()
        _run(_stop_payload(), tmp_path)
        marker = _read_marker(tmp_path)
        assert marker["schema"] == 1
        assert ISO_Z.match(marker["last_session_end_utc"])

    def test_timestamp_is_recent(self, tmp_path):
        (tmp_path / ".agents").mkdir()
        _run(_stop_payload(), tmp_path)
        ts = _read_marker(tmp_path)["last_session_end_utc"]
        written = datetime.strptime(ts, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
        assert abs((datetime.now(timezone.utc) - written).total_seconds()) < 120

    def test_empty_stdin_still_exits_zero_and_writes(self, tmp_path):
        # Fail-open: a malformed/empty payload must not crash the hook.
        (tmp_path / ".agents").mkdir()
        rc, out, _ = _run("", tmp_path)
        assert rc == 0
        assert (tmp_path / MARKER_RELPATH).exists()

    def test_overwrites_existing_marker(self, tmp_path):
        (tmp_path / ".agents").mkdir()
        marker_path = tmp_path / MARKER_RELPATH
        marker_path.write_text(json.dumps({"last_session_end_utc": "1999-01-01T00:00:00Z", "schema": 1}))
        _run(_stop_payload(), tmp_path)
        assert _read_marker(tmp_path)["last_session_end_utc"] != "1999-01-01T00:00:00Z"

    def test_missing_agents_dir_is_created(self, tmp_path):
        # No .agents/ pre-created: the hook should mkdir it rather than fail.
        rc, _, _ = _run(_stop_payload(), tmp_path)
        assert rc == 0
        assert (tmp_path / MARKER_RELPATH).exists()

    def test_direct_exec_via_shebang_and_exec_bit(self, tmp_path):
        # The harness invokes the hook as a bare executable (shebang + exec bit),
        # not as `python <path>`. Lock that contract in with a real PATH so
        # `/usr/bin/env python3` in the shebang can resolve.
        (tmp_path / ".agents").mkdir()
        import os

        proc = subprocess.run(
            [str(HOOK)],
            input=_stop_payload(),
            capture_output=True,
            text=True,
            timeout=5,
            env={"CLAUDE_PROJECT_DIR": str(tmp_path), "PATH": os.environ.get("PATH", "")},
        )
        assert proc.returncode == 0
        assert (tmp_path / MARKER_RELPATH).exists()


class TestAtomicWriteGuarantees:
    def test_no_stray_temp_after_successful_write(self, tmp_path):
        # The headline guarantee: the tempfile is renamed into place, never left.
        mod = _load_module()
        mod.write_watermark(tmp_path)
        leftover = list((tmp_path / ".agents").glob(".lsm-*"))
        assert leftover == []

    def test_cleanup_on_replace_failure_and_reraise(self, tmp_path, monkeypatch):
        # If os.replace fails mid-write, the temp must be unlinked and the error
        # re-raised (main() then swallows it for fail-open).
        mod = _load_module()

        def boom(_src, _dst):
            raise OSError("simulated replace failure")

        monkeypatch.setattr(mod.os, "replace", boom)
        with pytest.raises(OSError):
            mod.write_watermark(tmp_path)
        leftover = list((tmp_path / ".agents").glob(".lsm-*"))
        assert leftover == []


class TestWriteWatermark:
    def test_injected_timestamp_is_deterministic(self, tmp_path):
        mod = _load_module()
        path = mod.write_watermark(tmp_path, timestamp="2026-06-26T18:30:00Z")
        body = json.loads(path.read_text())
        assert body == {"last_session_end_utc": "2026-06-26T18:30:00Z", "schema": 1}

    def test_returns_marker_under_agents(self, tmp_path):
        mod = _load_module()
        path = mod.write_watermark(tmp_path, timestamp="2026-06-26T18:30:00Z")
        assert path == tmp_path / MARKER_RELPATH

    def test_default_timestamp_matches_iso_z(self, tmp_path):
        mod = _load_module()
        path = mod.write_watermark(tmp_path)
        ts = json.loads(path.read_text())["last_session_end_utc"]
        assert ISO_Z.match(ts)

    def test_cwd_from_payload_used_when_env_absent(self, tmp_path, monkeypatch):
        monkeypatch.delenv("CLAUDE_PROJECT_DIR", raising=False)
        mod = _load_module()
        root = mod._project_root({"cwd": str(tmp_path)})
        assert root == tmp_path

    def test_env_takes_precedence_over_payload_cwd(self, tmp_path, monkeypatch):
        monkeypatch.setenv("CLAUDE_PROJECT_DIR", str(tmp_path))
        mod = _load_module()
        root = mod._project_root({"cwd": "/some/other/path"})
        assert root == tmp_path


class TestRoutineMarker:
    """Issue #1002: a second, routine-only watermark distinct from the session
    marker. The Stop hook keeps writing the session marker every turn-end; the
    morning-routine recap beat writes the routine marker via `--marker routine`,
    and the recap/closure-audit lookback anchors on that one instead."""

    def test_default_marker_is_session(self, tmp_path):
        # Backward compat: no marker arg == the session marker the Stop hook writes.
        mod = _load_module()
        path = mod.write_watermark(tmp_path, timestamp="2026-07-06T12:00:00Z")
        assert path == tmp_path / MARKER_RELPATH
        assert json.loads(path.read_text()) == {
            "last_session_end_utc": "2026-07-06T12:00:00Z",
            "schema": 1,
        }

    def test_routine_marker_path_and_key(self, tmp_path):
        mod = _load_module()
        path = mod.write_watermark(tmp_path, timestamp="2026-07-06T12:00:00Z", marker="routine")
        assert path == tmp_path / ROUTINE_MARKER_RELPATH
        assert json.loads(path.read_text()) == {
            "last_routine_end_utc": "2026-07-06T12:00:00Z",
            "schema": 1,
        }

    def test_session_and_routine_markers_are_independent(self, tmp_path):
        # Writing one marker must never disturb the other (distinct files).
        mod = _load_module()
        mod.write_watermark(tmp_path, timestamp="2026-07-01T00:00:00Z", marker="session")
        mod.write_watermark(tmp_path, timestamp="2026-07-06T00:00:00Z", marker="routine")
        session = json.loads((tmp_path / MARKER_RELPATH).read_text())
        routine = json.loads((tmp_path / ROUTINE_MARKER_RELPATH).read_text())
        assert session["last_session_end_utc"] == "2026-07-01T00:00:00Z"
        assert routine["last_routine_end_utc"] == "2026-07-06T00:00:00Z"

    def test_unknown_marker_raises(self, tmp_path):
        mod = _load_module()
        with pytest.raises(KeyError):
            mod.write_watermark(tmp_path, marker="nope")

    def test_cli_marker_routine_writes_only_routine_file(self, tmp_path):
        (tmp_path / ".agents").mkdir()
        proc = subprocess.run(
            [sys.executable, str(HOOK), "--marker", "routine"],
            input=_stop_payload(),
            capture_output=True,
            text=True,
            timeout=5,
            env={"CLAUDE_PROJECT_DIR": str(tmp_path), "PATH": ""},
        )
        assert proc.returncode == 0
        assert proc.stdout.strip() == ""
        assert (tmp_path / ROUTINE_MARKER_RELPATH).exists()
        assert not (tmp_path / MARKER_RELPATH).exists()
        ts = json.loads((tmp_path / ROUTINE_MARKER_RELPATH).read_text())["last_routine_end_utc"]
        assert ISO_Z.match(ts)

    def test_bad_marker_arg_exits_zero_and_falls_back_to_session(self, tmp_path):
        # Stop-hook safety (review finding, PR #1057): argparse exits 2 on a bad
        # arg, and exit 2 is the "block the stop" signal a Stop hook must never
        # emit. main() must swallow it, write the default (session) marker, and
        # exit 0 - never propagate the non-zero.
        (tmp_path / ".agents").mkdir()
        proc = subprocess.run(
            [sys.executable, str(HOOK), "--marker", "bogus"],
            input=_stop_payload(),
            capture_output=True,
            text=True,
            timeout=5,
            env={"CLAUDE_PROJECT_DIR": str(tmp_path), "PATH": ""},
        )
        assert proc.returncode == 0
        assert (tmp_path / MARKER_RELPATH).exists()
        assert not (tmp_path / ROUTINE_MARKER_RELPATH).exists()

    def test_cli_default_writes_only_session_file(self, tmp_path):
        # No --marker == the Stop-hook path == session marker only.
        (tmp_path / ".agents").mkdir()
        proc = subprocess.run(
            [sys.executable, str(HOOK)],
            input=_stop_payload(),
            capture_output=True,
            text=True,
            timeout=5,
            env={"CLAUDE_PROJECT_DIR": str(tmp_path), "PATH": ""},
        )
        assert proc.returncode == 0
        assert (tmp_path / MARKER_RELPATH).exists()
        assert not (tmp_path / ROUTINE_MARKER_RELPATH).exists()
