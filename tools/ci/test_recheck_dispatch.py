"""Integration tests for `.agents/hooks/recheck_dispatch.py`.

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

HOOK = Path(__file__).parent.parent.parent / ".agents" / "hooks" / "recheck_dispatch.py"


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


class TestMovePathThreadsMovedIssue:
    """The `gh issue edit N --milestone X` move path threads `--moved-issue N`
    into every recheck it fires, so recheck_milestone.py can reconcile the laggy
    listing endpoint against a strongly-consistent read post-move (Issue #406).

    Imported in-process and driven through `dispatch()` with the gh-touching
    internals monkeypatched, so no live auth is needed to prove the wiring.
    """

    def _import_hook(self):
        import importlib.util

        spec = importlib.util.spec_from_file_location("recheck_dispatch_wiring", HOOK)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def _stub_common(self, monkeypatch, mod, calls):
        monkeypatch.setattr(mod, "run_recheck",
                            lambda *a: calls.append(a) or "Status: [No change]")
        monkeypatch.setattr(mod, "apply_target_sync", lambda n: None)

    def test_move_with_history_threads_moved_issue_to_each_milestone(self, monkeypatch):
        mod = self._import_hook()
        calls = []
        self._stub_common(monkeypatch, mod, calls)
        monkeypatch.setattr(mod, "prior_milestones_for_issue", lambda n: [17, 27])
        mod.dispatch('gh issue edit 381 --milestone "i3 - S4 - Dest"')
        assert ("--milestone", "17", "--moved-issue", "381") in calls
        assert ("--milestone", "27", "--moved-issue", "381") in calls

    def test_move_with_empty_history_threads_moved_issue(self, monkeypatch):
        mod = self._import_hook()
        calls = []
        self._stub_common(monkeypatch, mod, calls)
        monkeypatch.setattr(mod, "prior_milestones_for_issue", lambda n: [])
        mod.dispatch('gh issue edit 381 --milestone "i3 - S4 - Dest"')
        assert ("--issue", "381", "--moved-issue", "381") in calls


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
