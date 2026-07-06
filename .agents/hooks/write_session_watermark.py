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

Two watermarks, one writer (Issue #1002). The Stop hook (no args) writes the
`session` marker every turn-end. The morning-routine recap beat invokes this
same script as a CLI with `--marker routine` to stamp a distinct
`.agents/last_routine_marker.json`, so the recap/closure-audit lookback can
anchor on "last routine" instead of "last session-end" (which advances after
*any* session, including a non-recap one). The two markers live in separate
files with separate timestamp keys and never disturb each other.
"""
import argparse
import json
import os
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

SCHEMA = 1

# Named watermarks. Each selects its own file + timestamp key so the session
# marker (Stop hook) and the routine marker (recap beat, Issue #1002) are
# independent: writing one never clobbers the other.
MARKERS = {
    "session": {
        "relpath": (".agents", "last_session_marker.json"),
        "key": "last_session_end_utc",
    },
    "routine": {
        "relpath": (".agents", "last_routine_marker.json"),
        "key": "last_routine_end_utc",
    },
}
DEFAULT_MARKER = "session"


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


def write_watermark(root, timestamp=None, marker=DEFAULT_MARKER):
    """Atomically write the named watermark JSON under `root`; return its path.

    `marker` selects which watermark (see MARKERS): "session" (default, the
    Stop-hook marker) or "routine" (the recap-beat marker, Issue #1002). Uses a
    temp file in the marker's own directory + os.replace so a concurrent reader
    never observes a partially written file. Raises KeyError on an unknown name.
    """
    spec = MARKERS[marker]
    marker = Path(root).joinpath(*spec["relpath"])
    marker.parent.mkdir(parents=True, exist_ok=True)
    body = {
        spec["key"]: timestamp or _utc_now_iso(),
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


def _parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Write a session/routine watermark marker (Issues #820, #1002).",
    )
    parser.add_argument(
        "--marker",
        choices=sorted(MARKERS),
        default=DEFAULT_MARKER,
        help="which watermark to write (default: session, the Stop-hook marker)",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = _parse_args(argv)
    payload = {}
    # As a Stop hook, stdin carries the hook JSON (used only for cwd fallback).
    # As a CLI (`--marker routine`), skip the read when stdin is a TTY so we
    # never block waiting for input. A closed/piped-empty stdin parses to {}.
    if not sys.stdin.isatty():
        try:
            payload = json.load(sys.stdin)
        except Exception:
            # Fail open on ANY read/parse error (malformed JSON, OSError on
            # stdin): an absent payload just falls back to env/cwd resolution.
            payload = {}
    try:
        write_watermark(_project_root(payload), marker=args.marker)
    except Exception:
        # Fail open: a watermark failure must never block session stop.
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
