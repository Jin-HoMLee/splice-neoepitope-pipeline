#!/usr/bin/env python3
"""Stop hook: write a deterministic last-session watermark (Issue #820).

Records the UTC timestamp at which an agent turn last ended, so the morning
routine can anchor its "activity since last check" window on a reliable
high-water mark instead of inferring it from `episodes/` filenames (which can be
a cross-midnight session tail) or a reflexive "yesterday" (which silently drops
a weekend/absence gap). The morning routine reads this marker as the window
floor with a ~1-day backward overlap, and falls back to a conservative 7-day
floor when the marker is absent (first run / lost clone).

Wired as a `Stop` hook: it fires at the end of every agent turn, so the marker
always reflects the most recent turn-end — strictly more robust than a
`SessionEnd` hook, which never fires if a session is killed/crashes/left open
(the exact stale-watermark failure this issue removes).

Contract: reads Stop hook JSON on stdin, writes the marker, always exits 0 with
no stdout (a Stop hook must never block the agent from stopping), and fails open
on any error (a watermark write must never break a session).

The marker (`.agents/last_session_marker.json`) is a gitignored per-clone local
artifact, mirroring `.agents/hook_fires.jsonl`.
"""
import json
import os
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

SCHEMA = 1
MARKER_RELPATH = (".agents", "last_session_marker.json")


def _project_root(payload):
    """Resolve the clone root: CLAUDE_PROJECT_DIR env > payload cwd > process cwd."""
    root = os.environ.get("CLAUDE_PROJECT_DIR")
    if root:
        return Path(root)
    cwd = payload.get("cwd")
    if cwd:
        return Path(cwd)
    return Path.cwd()


def _utc_now_iso():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def write_watermark(root, timestamp=None):
    """Atomically write the watermark JSON under `root`; return the marker path.

    Uses a temp file in the marker's own directory + os.replace so a concurrent
    reader never observes a partially written file.
    """
    marker = Path(root).joinpath(*MARKER_RELPATH)
    marker.parent.mkdir(parents=True, exist_ok=True)
    body = {
        "last_session_end_utc": timestamp or _utc_now_iso(),
        "schema": SCHEMA,
    }
    fd, tmp = tempfile.mkstemp(dir=str(marker.parent), prefix=".lsm-", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(body, fh)
            fh.write("\n")
        os.replace(tmp, marker)
    except Exception:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise
    return marker


def main():
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        payload = {}
    try:
        write_watermark(_project_root(payload))
    except Exception:
        # Fail open: a watermark failure must never block session stop.
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
