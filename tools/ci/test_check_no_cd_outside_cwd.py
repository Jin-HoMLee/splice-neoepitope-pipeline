"""Tests for `.agents/hooks/check_no_cd_outside_cwd.py` (Issue #626).

Pure-function unit tests pass an explicit project_root / env / session_cwd so
they don't depend on the real clone location. Subprocess tests pipe PreToolUse
JSON to stdin against the real hook (which derives PROJECT_ROOT from its own
path) and exercise deny / allow / escape-hatch / fail-open end-to-end.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

HOOKS_DIR = Path(__file__).parent.parent.parent / ".agents" / "hooks"
HOOK = HOOKS_DIR / "check_no_cd_outside_cwd.py"
REAL_ROOT = str(HOOKS_DIR.parent.parent)  # the actual project clone root
sys.path.insert(0, str(HOOKS_DIR))

import check_no_cd_outside_cwd as h  # noqa: E402

ENV = {"HOME": "/home/user"}
ROOT = "/proj"


# --- cd_target ---


class TestCdTarget:
    def test_not_a_cd(self):
        assert h.cd_target(["ls", "-la"]) is None

    def test_bare_cd_is_home(self):
        assert h.cd_target(["cd"]) == ""

    def test_plain_target(self):
        assert h.cd_target(["cd", "workflow"]) == "workflow"

    def test_dash_previous(self):
        assert h.cd_target(["cd", "-"]) == "-"

    def test_skips_options(self):
        assert h.cd_target(["cd", "-P", "workflow"]) == "workflow"

    def test_options_only_is_home(self):
        assert h.cd_target(["cd", "-L"]) == ""


# --- resolve_cd ---


class TestResolveCd:
    def test_absolute(self):
        assert h.resolve_cd("/a/b", "/proj", ENV) == "/a/b"

    def test_relative_against_base(self):
        assert h.resolve_cd("workflow", "/proj", ENV) == "/proj/workflow"

    def test_dotdot(self):
        assert h.resolve_cd("..", "/proj", ENV) == "/"

    def test_home_bare(self):
        assert h.resolve_cd("", "/proj", ENV) == "/home/user"

    def test_tilde(self):
        assert h.resolve_cd("~", "/proj", ENV) == "/home/user"
        assert h.resolve_cd("~/x", "/proj", ENV) == "/home/user/x"

    def test_dynamic_unresolvable(self):
        assert h.resolve_cd("$HOME/x", "/proj", ENV) is None
        assert h.resolve_cd("`pwd`/x", "/proj", ENV) is None

    def test_dash_unresolvable(self):
        assert h.resolve_cd("-", "/proj", ENV) is None

    def test_tilde_user_unresolvable(self):
        assert h.resolve_cd("~bob", "/proj", ENV) is None

    def test_relative_no_base(self):
        assert h.resolve_cd("workflow", None, ENV) is None


# --- is_inside_safe ---


class TestIsInsideSafe:
    def test_root_itself(self):
        assert h.is_inside_safe("/proj", ROOT, ENV)

    def test_subdir(self):
        assert h.is_inside_safe("/proj/workflow", ROOT, ENV)

    def test_outside(self):
        assert not h.is_inside_safe("/other/repo", ROOT, ENV)

    def test_sibling_prefix_not_inside(self):
        # "/projector" must not count as inside "/proj"
        assert not h.is_inside_safe("/projector", ROOT, ENV)

    def test_temp_root_allowed(self):
        assert h.is_inside_safe("/tmp/scratch", ROOT, ENV)
        assert h.is_inside_safe("/private/tmp/x", ROOT, ENV)


# --- first_escaping_cd ---


class TestFirstEscapingCd:
    def test_outside_absolute(self):
        assert h.first_escaping_cd("cd /other/repo", ROOT, ROOT, ENV) == "/other/repo"

    def test_subdir_allowed(self):
        assert h.first_escaping_cd("cd workflow", ROOT, ROOT, ENV) is None

    def test_abs_subdir_allowed(self):
        assert h.first_escaping_cd("cd /proj/workflow", ROOT, ROOT, ENV) is None

    def test_dotdot_escapes(self):
        assert h.first_escaping_cd("cd ..", ROOT, ROOT, ENV) == ".."

    def test_chained_net_inside_allowed(self):
        # cd into subdir, then back up to root -> still inside
        assert h.first_escaping_cd("cd workflow && cd ..", ROOT, ROOT, ENV) is None

    def test_temp_allowed(self):
        assert h.first_escaping_cd("cd /tmp/scratch", ROOT, ROOT, ENV) is None

    def test_dash_fails_open(self):
        assert h.first_escaping_cd("cd -", ROOT, ROOT, ENV) is None

    def test_dynamic_fails_open(self):
        assert h.first_escaping_cd("cd $HOME/other", ROOT, ROOT, ENV) is None

    def test_tilde_escapes(self):
        assert h.first_escaping_cd("cd ~", ROOT, ROOT, ENV) == "~"

    def test_bare_cd_home_escapes(self):
        assert h.first_escaping_cd("cd", ROOT, ROOT, ENV) == ""

    def test_no_cd(self):
        assert h.first_escaping_cd("ls -la && git status", ROOT, ROOT, ENV) is None

    def test_cd_within_after_unresolvable(self):
        # after an unresolvable cd, a later relative cd fails open (running unknown)
        assert h.first_escaping_cd("cd $X && cd workflow", ROOT, ROOT, ENV) is None


# --- end-to-end via stdin (uses the real project root) ---


def _run(payload: dict, env=None):
    return subprocess.run(
        [sys.executable, str(HOOK)],
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        env=env,
    )


def _bash(cmd: str, cwd: str = REAL_ROOT) -> dict:
    return {"tool_name": "Bash", "tool_input": {"command": cmd}, "cwd": cwd}


class TestEndToEnd:
    def test_outside_denied(self):
        proc = _run(_bash("cd /Users/someone/other-repo && git status"))
        assert proc.returncode == 0
        out = json.loads(proc.stdout)
        assert out["hookSpecificOutput"]["permissionDecision"] == "deny"

    def test_subdir_allowed(self):
        assert _run(_bash("cd workflow && ls")).stdout.strip() == ""

    def test_no_cd_allowed(self):
        assert _run(_bash("git status")).stdout.strip() == ""

    def test_escape_hatch(self):
        env = dict(os.environ, CLAUDE_ALLOW_CD="1")
        assert _run(_bash("cd /etc"), env=env).stdout.strip() == ""

    def test_bad_json_fails_open(self):
        proc = subprocess.run(
            [sys.executable, str(HOOK)], input="not json", capture_output=True, text=True
        )
        assert proc.returncode == 0
        assert proc.stdout.strip() == ""

    def test_non_bash_ignored(self):
        assert _run({"tool_name": "Edit", "tool_input": {"command": "cd /etc"}}).stdout.strip() == ""


def test_hook_is_executable():
    # The harness execs the hook by bare path; a non-executable file EACCESes
    # before the shebang (the blocker caught on PR #1029).
    assert os.access(HOOK, os.X_OK), f"{HOOK} must be committed executable (100755)"
