"""Integration tests for `.claude/hooks/recheck_dispatch.py`.

Pipes a synthetic PostToolUse JSON payload to the hook and asserts that
the additionalContext includes the expected `[<recheck-name> — ...]` blocks.

The "silent" tests are pure pattern-matching and run in any environment.
The Status-field-trigger test additionally SHELLS OUT TO `gh` for item-ID
resolution + parent-chain walk. It probes for `gh auth` with project read
scope at collection time and skips gracefully if unavailable (e.g., fork
PRs where the GH_PROJECT_TOKEN secret isn't exposed, or local environments
without `gh auth login`).
"""

import json
import subprocess
import sys
from pathlib import Path

from _live_gh import REQUIRES_LIVE_GH

HOOK = Path(__file__).parent.parent.parent / ".claude" / "hooks" / "recheck_dispatch.py"


def _run(cmd: str, *args: str) -> tuple[int, str, str]:
    payload = json.dumps({"tool_input": {"command": cmd}})
    proc = subprocess.run(
        [sys.executable, str(HOOK), *args],
        input=payload,
        capture_output=True,
        text=True,
        timeout=60,
    )
    return proc.returncode, proc.stdout, proc.stderr


def _additional_context(stdout: str) -> str:
    if not stdout.strip():
        return ""
    return json.loads(stdout).get("hookSpecificOutput", {}).get("additionalContext", "")


class TestDispatcherSilent:
    def test_non_gh_command_silent(self):
        rc, out, _ = _run("ls -la")
        assert rc == 0
        assert out.strip() == ""

    def test_non_matching_gh_command_silent(self):
        rc, out, _ = _run("gh pr status")
        assert rc == 0
        assert out.strip() == ""


class TestScopeFilter:
    """Issue #454 scope-aware split: --scope shared (Sci/Dev) runs only shared checks.

    Hermetic — the pm-scope checks (recheck_milestone, recheck_parent_status) are
    skipped *before* any `gh` subprocess on a `gh issue close`, so no live auth is
    needed to prove the suppression.
    """

    def test_shared_scope_suppresses_pm_checks_on_close(self):
        # gh issue close → only pm-scope checks; under --scope shared they are
        # gated out before run_recheck / run_parent_status_recheck run → silent.
        rc, out, _ = _run("gh issue close 123456", "--scope", "shared")
        assert rc == 0
        assert out.strip() == ""

    def test_non_gh_command_silent_with_scope_flag(self):
        # Scope flag must not change the non-matching-command silent contract.
        rc, out, _ = _run("ls -la", "--scope", "shared")
        assert rc == 0
        assert out.strip() == ""


class TestStatusFieldTrigger:
    @REQUIRES_LIVE_GH
    def test_status_mutation_emits_parent_status_block(self):
        # Use a synthetic command that contains both the Status field ID and an item ID.
        # Item ID is for issue #24 (known, audited 2026-05-19, currently Ready).
        cmd = (
            'gh api graphql -f query=\'mutation { updateProjectV2ItemFieldValue('
            'input: { projectId: "PVT_kwHOB17eGc4BSomP" '
            'itemId: "PVTI_lAHOB17eGc4BSomPzgpgr30" '
            'fieldId: "PVTSSF_lAHOB17eGc4BSomPzhAHFf8" '
            'value: { singleSelectOptionId: "61e4505c" } }) { projectV2Item { id } } }\''
        )
        rc, out, _ = _run(cmd)
        assert rc == 0
        ctx = _additional_context(out)
        assert "[parent-status recheck — Status change on #24]" in ctx
