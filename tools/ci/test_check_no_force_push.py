"""Tests for `.agents/hooks/check_no_force_push.py` (Issue #626).

Pure-function unit tests import the module directly (the guard does no I/O, only
string inspection). Subprocess tests mirror the harness invocation (pipe
PreToolUse JSON to stdin) and exercise the deny / allow / escape-hatch /
fail-open paths end-to-end.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

HOOKS_DIR = Path(__file__).parent.parent.parent / ".agents" / "hooks"
HOOK = HOOKS_DIR / "check_no_force_push.py"
sys.path.insert(0, str(HOOKS_DIR))

import check_no_force_push as h  # noqa: E402


# --- _is_force_flag / is_force_push ---


class TestForceFlagDetection:
    def test_long_force(self):
        assert h._is_force_flag("--force")

    def test_short_force(self):
        assert h._is_force_flag("-f")

    def test_short_bundle_with_f(self):
        assert h._is_force_flag("-fq")
        assert h._is_force_flag("-qf")

    def test_force_with_lease_bare(self):
        assert h._is_force_flag("--force-with-lease")

    def test_force_with_lease_ref(self):
        assert h._is_force_flag("--force-with-lease=origin/main")

    def test_force_if_includes(self):
        assert h._is_force_flag("--force-if-includes")

    def test_non_force_flags(self):
        for t in ("-u", "--verbose", "-v", "origin", "main", "--set-upstream", "-n", "-q"):
            assert not h._is_force_flag(t)

    def test_double_dash_foo_not_force(self):
        assert not h._is_force_flag("--foo")


class TestIsForcePush:
    def test_plain_force(self):
        assert h.is_force_push(["git", "push", "--force"])

    def test_force_flag_last(self):
        assert h.is_force_push(["git", "push", "origin", "main", "--force"])

    def test_short_f(self):
        assert h.is_force_push(["git", "push", "-f"])

    def test_lease(self):
        assert h.is_force_push(["git", "push", "--force-with-lease"])

    def test_push_no_force_allowed(self):
        assert not h.is_force_push(["git", "push", "origin", "main"])

    def test_non_push_with_force_token(self):
        # a force flag but no `push` subcommand -> not a force push
        assert not h.is_force_push(["git", "commit", "--force"])

    def test_requires_git_token(self):
        assert not h.is_force_push(["push", "--force"])

    # subcommand resolution: `push` must be the git subcommand, not any token
    def test_ref_named_push_not_flagged(self):
        assert not h.is_force_push(["git", "checkout", "-f", "push"])
        assert not h.is_force_push(["git", "branch", "-f", "push", "origin/main"])
        assert not h.is_force_push(["git", "switch", "-f", "push"])

    def test_force_push_through_global_opts(self):
        assert h.is_force_push(["git", "-C", "/other", "push", "--force"])
        assert h.is_force_push(["git", "-c", "k=v", "push", "-f"])

    # rewrite forms beyond force flags
    def test_force_refspec(self):
        assert h.is_force_push(["git", "push", "origin", "+main"])

    def test_mirror(self):
        assert h.is_force_push(["git", "push", "--mirror"])

    # --- Issue #1134: a branch DELETE is not a history rewrite -> allowed ---
    # These previously asserted is_force_push == True (the false positive). Inverted:
    # a delete removes a pointer, rewrites no history, and is reversible, so the
    # guard must let it through.
    def test_delete_refspec_allowed(self):
        assert not h.is_force_push(["git", "push", "origin", ":old"])

    def test_delete_flag_allowed(self):
        assert not h.is_force_push(["git", "push", "--delete", "origin", "x"])
        assert not h.is_force_push(["git", "push", "-d", "origin", "x"])

    def test_delete_and_force_together_is_still_denied(self):
        """Belt-and-suspenders: a delete does not LAUNDER a force flag in the same
        command. `--delete` is allowed on its own, but `+ref` alongside still fires."""
        assert h.is_force_push(["git", "push", "--delete", "origin", "+main"])

    def test_normal_refspec_allowed(self):
        assert not h.is_force_push(["git", "push", "origin", "main:main"])

    def test_compound_delete_after_branch_D_is_allowed(self):
        """The second #1134 false positive: a cleanup sweep does
        `git branch -D x` then `git push origin --delete x`. Neither is a rewrite;
        the whole command must pass (the -D and --delete must not be conflated)."""
        assert not h.command_force_pushes(
            "git branch -D feat/x && git push origin --delete feat/x"
        )

    def test_compound_force_after_delete_is_still_denied(self):
        """...but a real force push in a later segment is still caught."""
        assert h.command_force_pushes(
            "git push origin --delete old && git push --force origin main"
        )


class TestGitSubcommand:
    def test_plain(self):
        assert h.git_subcommand(["git", "push", "--force"]) == "push"

    def test_skips_value_opt(self):
        assert h.git_subcommand(["git", "-c", "user.name=x", "commit"]) == "commit"
        assert h.git_subcommand(["git", "-C", "/path", "push"]) == "push"

    def test_skips_valueless_flag(self):
        assert h.git_subcommand(["git", "--no-pager", "log"]) == "log"

    def test_ref_named_push_resolves_to_checkout(self):
        assert h.git_subcommand(["git", "checkout", "-f", "push"]) == "checkout"

    def test_non_git(self):
        assert h.git_subcommand(["ls", "-la"]) is None


# --- split_subcommands ---


class TestSplitSubcommands:
    def test_single(self):
        assert h.split_subcommands("git push --force") == [["git", "push", "--force"]]

    def test_chained(self):
        subs = h.split_subcommands("git commit -m x && git push --force")
        assert subs == [["git", "commit", "-m", "x"], ["git", "push", "--force"]]

    def test_quoted_collapses(self):
        # the quoted message is ONE token; no bare push/--force leaks out
        subs = h.split_subcommands('git commit -m "push --force please"')
        assert subs == [["git", "commit", "-m", "push --force please"]]

    def test_unbalanced_quotes_fail_open(self):
        assert h.split_subcommands('git push "--force') is None


# --- command_force_pushes ---


class TestCommandForcePushes:
    def test_standalone_force(self):
        assert h.command_force_pushes("git push --force")

    def test_chained_after_commit(self):
        assert h.command_force_pushes("git commit -m x && git push -f")

    def test_plain_push_allowed(self):
        assert not h.command_force_pushes("git push origin main")

    def test_quoted_message_not_flagged(self):
        assert not h.command_force_pushes('git commit -m "push --force"')

    def test_force_in_other_subcommand_not_attributed(self):
        # --force belongs to a different subcommand than the push
        assert not h.command_force_pushes("git push && ./deploy --force")

    def test_parse_miss_allows(self):
        assert not h.command_force_pushes('git push "--force')


# --- end-to-end via stdin (subprocess) ---


def _run(payload: dict, env=None):
    proc = subprocess.run(
        [sys.executable, str(HOOK)],
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        env=env,
    )
    return proc


def _bash(cmd: str) -> dict:
    return {"tool_name": "Bash", "tool_input": {"command": cmd}}


class TestEndToEnd:
    def test_force_push_denied(self):
        proc = _run(_bash("git push --force origin main"))
        assert proc.returncode == 0
        out = json.loads(proc.stdout)
        assert out["hookSpecificOutput"]["permissionDecision"] == "deny"

    def test_plain_push_allowed(self):
        proc = _run(_bash("git push origin main"))
        assert proc.returncode == 0
        assert proc.stdout.strip() == ""

    def test_escape_hatch_allows(self):
        env = dict(os.environ, CLAUDE_ALLOW_FORCE_PUSH="1")
        proc = _run(_bash("git push --force"), env=env)
        assert proc.stdout.strip() == ""

    def test_bad_json_fails_open(self):
        proc = subprocess.run(
            [sys.executable, str(HOOK)],
            input="not json",
            capture_output=True,
            text=True,
        )
        assert proc.returncode == 0
        assert proc.stdout.strip() == ""

    def test_non_bash_tool_ignored(self):
        proc = _run({"tool_name": "Edit", "tool_input": {"command": "git push --force"}})
        assert proc.stdout.strip() == ""


def test_hook_is_executable():
    # The harness execs the hook by bare path, so a non-executable file (100644)
    # EACCESes before the shebang runs - the blocker caught on PR #1029. The
    # subprocess tests above invoke via `sys.executable <hook>`, which bypasses
    # the exec bit, so this is the one check that exercises the real harness path.
    assert os.access(HOOK, os.X_OK), f"{HOOK} must be committed executable (100755)"


class TestNewlineForm:
    """Issue #1142: two commands on two lines MERGED into one token block.

    A newline is a command separator in shell but not in `shlex`, so before the
    `normalize_command()` fix `split_subcommands` returned a single merged block
    for a two-line command. That broke this guard in BOTH directions:

      - false negative: `git commit --amend` on line 1 made `git_subcommand`
        resolve the merged block to `commit`, so the `git push --force` on line 2
        was never examined - a force push walked straight through the guard;
      - false positive: `git push origin main` + `git branch -f tmp` merged into
        one block that read as a forced push, denying a legitimate command.

    Found by the PR #1145 bot review, which probed with a *second real command*
    where the author's probe had used a harmless `echo` prefix - the merge only
    bites when both lines are git commands.
    """

    def test_force_push_after_a_commit_line_is_caught(self):
        assert h.command_force_pushes(
            "git commit --amend\ngit push --force-with-lease origin main") is True

    def test_force_push_alone_still_caught(self):
        assert h.command_force_pushes("git push --force-with-lease origin main") is True

    def test_force_push_after_heredoc_is_caught(self):
        cmd = ("cat > /tmp/msg.txt <<'EOF'\ncommit message\nEOF\n"
               "git push --force origin main")
        assert h.command_force_pushes(cmd) is True

    def test_push_then_branch_force_is_not_a_force_push(self):
        # The false-positive half. `git branch -f` is not a force PUSH.
        assert h.command_force_pushes(
            "git push origin main\ngit branch -f tmp HEAD") is False

    def test_plain_push_still_allowed(self):
        assert h.command_force_pushes("git push origin main") is False
