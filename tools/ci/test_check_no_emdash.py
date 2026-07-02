"""Tests for `.agents/hooks/check_no_emdash.py` (Issue #920).

Pure-function unit tests import the module directly (no I/O). The subprocess
tests mirror the harness invocation (pipe PreToolUse JSON to stdin, read the deny
decision from stdout; exit code is always 0).

The em-dash (U+2014) and en-dash (U+2013) are written as `\\u2014` / `\\u2013`
escapes throughout so this test file itself stays ASCII-clean and does not need
the guard's own escape hatch to be edited.
"""

import json
import subprocess
import sys
from pathlib import Path

HOOKS_DIR = Path(__file__).parent.parent.parent / ".agents" / "hooks"
HOOK = HOOKS_DIR / "check_no_emdash.py"
sys.path.insert(0, str(HOOKS_DIR))

import check_no_emdash as h  # noqa: E402

EM = "\u2014"  # em-dash
EN = "\u2013"  # en-dash


# --- pure: added_violations (net-added delta scan) ---


class TestAddedViolations:
    def test_plain_dash_only_clean(self):
        assert h.added_violations("a - b", "c - d - e") == []

    def test_em_added(self):
        assert "em-dash (U+2014)" in h.added_violations("plain", f"a{EM}b")

    def test_en_added(self):
        assert "en-dash (U+2013)" in h.added_violations("plain", f"S1{EN}S7")

    def test_both_added(self):
        got = h.added_violations("", f"{EM}{EN}")
        assert "em-dash (U+2014)" in got and "en-dash (U+2013)" in got

    def test_preserved_emdash_not_flagged(self):
        # old already has one em-dash; new keeps exactly one → net zero → clean.
        assert h.added_violations(f"keep {EM} me", f"keep {EM} me (edited)") == []

    def test_removed_emdash_not_flagged(self):
        assert h.added_violations(f"was {EM} here", "now plain -") == []

    def test_extra_emdash_beyond_preexisting_flagged(self):
        # old has one, new has two → one net-added → flagged.
        assert "em-dash (U+2014)" in h.added_violations(f"{EM}", f"{EM} and {EM}")


# --- pure: guarded_pairs (per-tool delta extraction) ---


class TestGuardedPairs:
    def test_edit_pair(self):
        ti = {"file_path": "x.py", "old_string": "a", "new_string": "b"}
        assert h.guarded_pairs("Edit", ti) == [("a", "b")]

    def test_multiedit_pairs(self):
        ti = {
            "file_path": "x.py",
            "edits": [
                {"old_string": "a", "new_string": "b"},
                {"old_string": "c", "new_string": "d"},
            ],
        }
        assert h.guarded_pairs("MultiEdit", ti) == [("a", "b"), ("c", "d")]

    def test_write_new_file_reads_empty(self):
        ti = {"file_path": "new.py", "content": f"x{EM}y"}
        pairs = h.guarded_pairs("Write", ti, read_file=lambda _p: "")
        assert pairs == [("", f"x{EM}y")]

    def test_write_existing_file_uses_prior_content(self):
        ti = {"file_path": "old.py", "content": f"keep {EM}"}
        pairs = h.guarded_pairs("Write", ti, read_file=lambda _p: f"keep {EM}")
        assert h.added_violations(*pairs[0]) == []

    def test_write_unreadable_fails_open(self):
        ti = {"file_path": "weird.bin", "content": f"x{EM}"}
        assert h.guarded_pairs("Write", ti, read_file=lambda _p: None) is None

    def test_non_guarded_tool_returns_none(self):
        assert h.guarded_pairs("Bash", {"command": "ls"}) is None


# --- pure: is_allowlisted (escape hatch) ---


class TestIsAllowlisted:
    def test_env_override_truthy(self):
        assert h.is_allowlisted("any/file.py", env={"CLAUDE_ALLOW_EMDASH": "1"}) is True

    def test_env_override_falsey(self):
        assert h.is_allowlisted("any/file.py", env={"CLAUDE_ALLOW_EMDASH": "0"}) is False

    def test_allowlisted_path_substring(self):
        assert h.is_allowlisted("tools/ci/test_check_no_emdash.py", env={}) is True

    def test_normal_path_not_allowlisted(self):
        assert h.is_allowlisted("workflow/scripts/foo.py", env={}) is False


# --- subprocess: end-to-end harness behavior ---


def _run(payload, env=None):
    proc = subprocess.run(
        [sys.executable, str(HOOK)],
        input=payload,
        capture_output=True,
        text=True,
        timeout=5,
        env=env,
    )
    return proc.returncode, proc.stdout, proc.stderr


def _edit_payload(old, new, file_path="workflow/scripts/foo.py"):
    return json.dumps(
        {
            "tool_name": "Edit",
            "tool_input": {"file_path": file_path, "old_string": old, "new_string": new},
        }
    )


def _is_deny(stdout):
    if not stdout.strip():
        return False
    return (
        json.loads(stdout).get("hookSpecificOutput", {}).get("permissionDecision")
        == "deny"
    )


class TestSubprocessDeny:
    def test_edit_adding_emdash_denied(self):
        rc, out, _ = _run(_edit_payload("plain", f"a{EM}b"))
        assert rc == 0 and _is_deny(out)

    def test_edit_adding_endash_denied(self):
        rc, out, _ = _run(_edit_payload("plain", f"S1{EN}S7"))
        assert rc == 0 and _is_deny(out)

    def test_write_new_file_with_emdash_denied(self, tmp_path):
        target = tmp_path / "brand_new.py"
        payload = json.dumps(
            {
                "tool_name": "Write",
                "tool_input": {"file_path": str(target), "content": f"x{EM}y"},
            }
        )
        rc, out, _ = _run(payload)
        assert rc == 0 and _is_deny(out)

    def test_deny_message_points_to_plain_dash_rule(self):
        rc, out, _ = _run(_edit_payload("plain", f"a{EM}b"))
        reason = json.loads(out)["hookSpecificOutput"]["permissionDecisionReason"]
        assert "plain" in reason.lower() and "CLAUDE_ALLOW_EMDASH" in reason


class TestSubprocessAllow:
    def test_edit_plain_dash_allowed(self):
        rc, out, _ = _run(_edit_payload("plain", "a - b - c"))
        assert rc == 0 and out.strip() == ""

    def test_edit_preserving_preexisting_emdash_allowed(self):
        rc, out, _ = _run(_edit_payload(f"keep {EM}", f"keep {EM} now"))
        assert rc == 0 and not _is_deny(out)

    def test_allowlisted_path_allowed(self):
        rc, out, _ = _run(
            _edit_payload("plain", f"x{EM}y", file_path="tools/ci/test_check_no_emdash.py")
        )
        assert rc == 0 and not _is_deny(out)

    def test_env_escape_hatch_allowed(self):
        import os

        env = dict(os.environ, CLAUDE_ALLOW_EMDASH="1")
        rc, out, _ = _run(_edit_payload("plain", f"x{EM}y"), env=env)
        assert rc == 0 and not _is_deny(out)


class TestFailOpen:
    def test_empty_stdin_allowed(self):
        rc, out, _ = _run("")
        assert rc == 0 and out.strip() == ""

    def test_malformed_json_allowed(self):
        rc, out, _ = _run("{not valid json")
        assert rc == 0 and out.strip() == ""

    def test_non_guarded_tool_allowed(self):
        payload = json.dumps({"tool_name": "Bash", "tool_input": {"command": "ls"}})
        rc, out, _ = _run(payload)
        assert rc == 0 and out.strip() == ""
