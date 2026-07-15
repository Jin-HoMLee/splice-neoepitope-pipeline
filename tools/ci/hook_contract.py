#!/usr/bin/env python3
"""Shared harness for the hook-liveness contract suite (Issue #1140).

The contract suite proves each hook registered in `.agents/settings.json` still
*fires* when given its real trigger through its real entry path - the hook script
run as a subprocess reading the harness JSON envelope on stdin, exactly as the
harness invokes it. It exists because four governance hooks shipped completely
inert this month (each green on its unit tests, none ever invoked through its real
path); `post_gh_pr_create` was dead for months because `matches_pr_create()`
returned False for the heredoc-then-`gh pr create` shape that opens every real PR.

This module is the single-authored interface the per-hook contracts target:

- `registered_hooks()` enumerates the registry (the enumeration + drift AC).
- envelope builders produce the exact stdin payload per event/matcher.
- `drive()` runs a hook subprocess the way the harness does (stdin = envelope
  JSON), optionally with a stub `gh` first on PATH for hooks that reach a live
  `gh` call before their observable (the `test_set_status.py` pattern).
- observable detectors read the two artifact kinds a hook can produce: a
  `permissionDecision: deny` on stdout, or an appended `.agents/hook_fires.jsonl`
  line.

We deliberately do NOT re-implement the harness `if:` gate. That gate is
vendor-documented best-effort and every hook self-matches internally and fails
open, so it is a performance optimization, not a correctness boundary - and it is
the internal matcher, not the gate, where all four dead hooks actually failed. The
real entry path is the subprocess with the true command string.
"""
from __future__ import annotations

import json
import os
import stat
import subprocess
import sys
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterator, Optional

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
SETTINGS_PATH = REPO_ROOT / ".agents" / "settings.json"
HOOKS_DIR = REPO_ROOT / ".agents" / "hooks"
FIRE_LOG_PATH = REPO_ROOT / ".agents" / "hook_fires.jsonl"


# --- registry enumeration ----------------------------------------------------


@dataclass(frozen=True)
class Registration:
    """One `(event, matcher, if, command)` hook registration in settings.json."""

    event: str  # PreToolUse | PostToolUse | Stop
    matcher: str  # e.g. "Bash", "Edit|MultiEdit|Write", ""
    if_clause: str  # e.g. "Bash(gh *)", ""
    command: str  # the raw command string as registered
    basename: str  # the hook script filename, e.g. "check_no_force_push.py"

    @property
    def path(self) -> Path:
        return HOOKS_DIR / self.basename


def registered_hooks(settings_path: Path = SETTINGS_PATH) -> list[Registration]:
    """Parse settings.json and return one Registration per hook command.

    The suite is parametrized off this, so a newly registered hook with no
    contract entry fails the completeness test rather than being silently
    uncovered (the enumeration + drift AC of #1140).
    """
    data = json.loads(settings_path.read_text(encoding="utf-8"))
    out: list[Registration] = []
    for event, entries in (data.get("hooks") or {}).items():
        for entry in entries or []:
            matcher = entry.get("matcher", "") or ""
            for hk in entry.get("hooks") or []:
                cmd = hk.get("command", "") or ""
                out.append(
                    Registration(
                        event=event,
                        matcher=matcher,
                        if_clause=hk.get("if", "") or "",
                        command=cmd,
                        basename=_basename_of(cmd),
                    )
                )
    return out


def _basename_of(command: str) -> str:
    """The hook script filename from a registered command string.

    Registered commands look like `${CLAUDE_PROJECT_DIR}/.agents/hooks/foo.py`
    possibly with trailing args (`recheck_dispatch.py --scope shared`). We take
    the first token that ends in `.py`.
    """
    for tok in command.split():
        if tok.endswith(".py"):
            return Path(tok).name
    # Fall back to the last path-like token's name.
    return Path(command.split()[0]).name if command.split() else ""


def distinct_hook_basenames(settings_path: Path = SETTINGS_PATH) -> set[str]:
    """The set of distinct hook script filenames registered in settings.json.

    A single script may be registered under more than one matcher
    (`check_memory_path_cwd_drift.py` is on both Bash and Edit/Write); the
    contract is per-script, so we dedupe.
    """
    return {r.basename for r in registered_hooks(settings_path)}


# --- envelope builders -------------------------------------------------------
# Each returns the JSON dict the harness pipes on the hook's stdin. Keys match
# what the hooks actually read (`tool_name`, `tool_input.command`, etc.).


def pretooluse_bash(command: str, cwd: Optional[str] = None) -> dict:
    """PreToolUse Bash envelope. `cwd` populates the top-level session-cwd field
    the harness delivers, which `check_memory_path_cwd_drift.py` keys on to detect
    a drifted shell; omit it for hooks that only read `tool_input.command`."""
    env: dict = {
        "hook_event_name": "PreToolUse",
        "tool_name": "Bash",
        "tool_input": {"command": command},
    }
    if cwd is not None:
        env["cwd"] = cwd
    return env


def pretooluse_edit(file_path: str, old_string: str, new_string: str) -> dict:
    return {
        "hook_event_name": "PreToolUse",
        "tool_name": "Edit",
        "tool_input": {
            "file_path": file_path,
            "old_string": old_string,
            "new_string": new_string,
        },
    }


def pretooluse_write(file_path: str, content: str) -> dict:
    return {
        "hook_event_name": "PreToolUse",
        "tool_name": "Write",
        "tool_input": {"file_path": file_path, "content": content},
    }


def posttooluse_bash(command: str, tool_response: object) -> dict:
    return {
        "hook_event_name": "PostToolUse",
        "tool_name": "Bash",
        "tool_input": {"command": command},
        "tool_response": tool_response,
    }


def stop_envelope(cwd: Optional[str] = None) -> dict:
    env: dict = {"hook_event_name": "Stop"}
    if cwd is not None:
        env["cwd"] = cwd
    return env


# --- driving a hook the way the harness does ---------------------------------


@dataclass
class DriveResult:
    exit_code: int
    stdout: str
    stderr: str


def drive(
    hook_path: Path,
    envelope: dict,
    *,
    gh_stub: Optional[str] = None,
    env: Optional[dict] = None,
    extra_args: Optional[list[str]] = None,
    timeout: float = 30.0,
) -> DriveResult:
    """Run a hook script as a subprocess with `envelope` JSON on stdin.

    This is the real entry path: the harness invokes the registered command with
    the hook JSON on stdin. `gh_stub`, when given, is shell source written to a
    `gh` executable placed first on PATH, so hooks that reach a live `gh` call
    before their observable run end-to-end offline (the `test_set_status.py`
    pattern). `extra_args` mirrors a registration that carries args
    (`recheck_dispatch.py --scope shared`).
    """
    run_env = dict(os.environ if env is None else env)
    with _optional_gh_stub(gh_stub) as stub_dir:
        if stub_dir is not None:
            run_env["PATH"] = f"{stub_dir}{os.pathsep}{run_env.get('PATH', '')}"
        proc = subprocess.run(
            [sys.executable, str(hook_path), *(extra_args or [])],
            input=json.dumps(envelope),
            capture_output=True,
            text=True,
            env=run_env,
            timeout=timeout,
        )
    return DriveResult(proc.returncode, proc.stdout, proc.stderr)


@contextmanager
def _optional_gh_stub(source: Optional[str]) -> Iterator[Optional[Path]]:
    if source is None:
        yield None
        return
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        gh = Path(d) / "gh"
        gh.write_text(source, encoding="utf-8")
        gh.chmod(gh.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
        yield Path(d)


# --- observable detectors ----------------------------------------------------


def deny_decision(result: DriveResult) -> bool:
    """True iff the hook emitted a PreToolUse `permissionDecision: deny`."""
    if not result.stdout.strip():
        return False
    try:
        payload = json.loads(result.stdout)
    except json.JSONDecodeError:
        return False
    decision = payload.get("hookSpecificOutput") or payload
    return decision.get("permissionDecision") == "deny"


@contextmanager
def fire_log_guard(log_path: Path = FIRE_LOG_PATH) -> Iterator["FireLogWatcher"]:
    """Snapshot the real (gitignored) fire-log so the suite leaves it untouched.

    Yields a watcher whose `.fired()` reports whether new lines were appended
    across a `drive()` call. On exit the log is restored to its pre-test bytes
    (or removed if it did not exist), so contract tests never pollute the real
    `.agents/hook_fires.jsonl`.
    """
    existed = log_path.exists()
    original = log_path.read_bytes() if existed else b""
    watcher = FireLogWatcher(log_path, len(original.splitlines()))
    try:
        yield watcher
    finally:
        if existed:
            log_path.write_bytes(original)
        elif log_path.exists():
            log_path.unlink()


@dataclass
class FireLogWatcher:
    log_path: Path
    _baseline_lines: int

    def fired(self) -> bool:
        """True iff the fire-log has grown since the guard opened."""
        if not self.log_path.exists():
            return False
        current = len(self.log_path.read_bytes().splitlines())
        return current > self._baseline_lines

    def new_lines(self) -> list[dict]:
        if not self.log_path.exists():
            return []
        lines = self.log_path.read_text(encoding="utf-8").splitlines()[self._baseline_lines :]
        out = []
        for ln in lines:
            try:
                out.append(json.loads(ln))
            except json.JSONDecodeError:
                pass
        return out


# --- per-hook contract shape (populated by the suite's registry) -------------

DENY = "deny"
FIRE_LOG = "fire_log"
WATERMARK = "watermark"

_BUILDERS: dict[str, Callable[..., dict]] = {
    "pretooluse_bash": pretooluse_bash,
    "pretooluse_edit": pretooluse_edit,
    "pretooluse_write": pretooluse_write,
    "posttooluse_bash": posttooluse_bash,
    "stop_envelope": stop_envelope,
}


@dataclass(frozen=True)
class HookContract:
    """The liveness contract for one hook script.

    `fire_input` is the real trigger that MUST produce `observable`; `nofire_input`
    is the matched-pair counterexample that must produce nothing. Both are the
    kwargs dicts passed to the named `envelope_builder`. `heredoc_fire_input`, when
    set, is an additional fire trigger using the heredoc-then-command shape (the
    exact string class that killed `matches_pr_create`). `gh_stub_source` /
    `extra_args` cover the live-`gh` and arg-carrying registrations.
    """

    basename: str
    observable: str  # DENY | FIRE_LOG | WATERMARK
    envelope_builder: str  # a key of _BUILDERS
    fire_input: dict
    nofire_input: dict
    heredoc_fire_input: Optional[dict] = None
    gh_stub_source: Optional[str] = None
    extra_args: Optional[list[str]] = None
    notes: str = ""

    def envelope(self, spec_input: dict) -> dict:
        return _BUILDERS[self.envelope_builder](**spec_input)

    @property
    def hook_path(self) -> Path:
        return HOOKS_DIR / self.basename


def observe(contract: HookContract, spec_input: dict) -> bool:
    """Drive `contract`'s hook with `spec_input` and report whether its observable
    artifact appeared - the single place the three artifact kinds are detected.

    Runs inside a `fire_log_guard` so the real (gitignored) fire-log is restored
    afterward regardless of the observable kind.
    """
    if contract.observable == WATERMARK:
        return _observe_watermark(contract, spec_input)
    with fire_log_guard() as watcher:
        result = drive(
            contract.hook_path,
            contract.envelope(spec_input),
            gh_stub=contract.gh_stub_source,
            extra_args=contract.extra_args,
        )
        if contract.observable == DENY:
            return deny_decision(result)
        return watcher.fired()


WRITABLE_ROOT = "<TMP>"  # sentinel: resolve to a fresh writable tempdir at drive time


def _observe_watermark(contract: HookContract, spec_input: dict) -> bool:
    """Stop-hook watermark observable.

    `write_session_watermark` has no discriminating matcher - it writes the marker
    on every turn-end - so its only "no observable" state is a resolved root that
    cannot be written (it fails open on OSError). The matched pair is therefore a
    writable root (writes) vs an unwritable root (silent). The root is taken from
    `spec_input["root"]`: the `WRITABLE_ROOT` sentinel becomes a fresh writable
    tempdir, any other value (e.g. `/dev/null/sub`) is used literally so mkdir
    fails and no marker appears. The root is passed as `CLAUDE_PROJECT_DIR`, which
    the hook ranks first when resolving where to write.
    """
    import tempfile

    root_spec = spec_input.get("root", WRITABLE_ROOT)
    if root_spec == WRITABLE_ROOT:
        with tempfile.TemporaryDirectory() as d:
            return _drive_watermark(contract, d)
    return _drive_watermark(contract, root_spec)


def _drive_watermark(contract: HookContract, root: str) -> bool:
    env = dict(os.environ, CLAUDE_PROJECT_DIR=root)
    root_path = Path(root)
    before = set(root_path.rglob("*")) if root_path.exists() else set()
    drive(contract.hook_path, stop_envelope(), env=env, extra_args=contract.extra_args)
    if not root_path.exists():
        return False  # unwritable root: hook failed open, nothing written
    return bool(set(root_path.rglob("*")) - before)
