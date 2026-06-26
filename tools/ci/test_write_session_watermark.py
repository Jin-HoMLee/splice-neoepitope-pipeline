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

HOOK = Path(__file__).parent.parent.parent / ".agents" / "hooks" / "write_session_watermark.py"
MARKER_RELPATH = Path(".agents") / "last_session_marker.json"
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
