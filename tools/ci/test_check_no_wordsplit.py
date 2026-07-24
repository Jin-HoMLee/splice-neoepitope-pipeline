"""Tests for `.agents/hooks/check_no_wordsplit.py` (Issue #1242).

Pure-function unit tests import the module directly (the guard does no I/O, only
string inspection). Subprocess tests mirror the harness invocation (pipe
PreToolUse JSON to stdin) and exercise the deny / allow / escape-hatch /
fail-open paths end-to-end.

The precision fixtures are the load-bearing half: an over-eager word-split guard
is worse than none, so the matched-pair PASS cases (`"${arr[@]}"`, `$(cmd)`,
`"$@"`, literal lists, `set -e`) matter as much as the FAIL cases.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

HOOKS_DIR = Path(__file__).parent.parent.parent / ".agents" / "hooks"
HOOK = HOOKS_DIR / "check_no_wordsplit.py"
sys.path.insert(0, str(HOOKS_DIR))

import check_no_wordsplit as h  # noqa: E402


# --- the anti-pattern FIRES ---


class TestFires:
    @pytest.mark.parametrize(
        "cmd",
        [
            "for x in $var; do echo $x; done",
            "for f in ${files}; do rm $f; done",
            "set -- $pair",
            "set -- ${blob}",
            # not at the very start of the command, but at a command position
            "cd /tmp && for x in $var; do echo $x; done",
            "echo start; set -- $args",
            # inside a loop body
            "for a in one two; do set -- $rest; done",
            # multi-line (normalize turns the newline into a separator)
            "cd /tmp\nfor x in $var; do echo $x; done",
            # a good operand first, a bad one second
            'for x in "${good[@]}" $bad; do :; done',
        ],
    )
    def test_wordsplit_fires(self, cmd):
        assert h.command_wordsplits(cmd) is True


# --- high precision: the correct forms PASS ---


class TestPasses:
    @pytest.mark.parametrize(
        "cmd",
        [
            # the correct forms named in the AC
            'for f in "${arr[@]}"; do echo "$f"; done',
            "for f in a b c; do echo $f; done",
            'set -- "$@"',
            # substitution and glob operands
            "for f in $(ls); do echo $f; done",
            "for n in $(seq 1 3); do echo $n; done",
            "for f in *.txt; do echo $f; done",
            # arithmetic substitution is masked, not a scalar
            "for i in $((1+2)); do echo $i; done",
            # unquoted ARRAY expansion is out of the named-scalar scope
            "for f in ${arr[@]}; do echo $f; done",
            # special params are not named scalars
            "set -- $@",
            "for i in $1 $2; do echo $i; done",
            # `set` without `--` is a shell-option set, never a positional reset
            "set -euo pipefail",
            "set -x",
            # a C-style for loop has no `in`
            "for ((i=0; i<3; i++)); do echo $i; done",
            # the trigger words appear only as literal arguments, not commands
            "echo for f in $x",
            'grep -r "for x in $pattern" .',
            # a var inside a quoted string is not a loop operand
            'echo "for x in $var"',
            # the fixed streamed form
            'printf "%s\\n" "$blob" | while read -r x; do echo $x; done',
        ],
    )
    def test_correct_form_passes(self, cmd):
        assert h.command_wordsplits(cmd) is False


# --- mask_opaque unit coverage ---


class TestMaskOpaque:
    def test_single_quote_masked(self):
        assert "$x" not in h.mask_opaque("echo 'a $x b'")

    def test_double_quote_masked(self):
        assert "$x" not in h.mask_opaque('echo "a $x b"')

    def test_command_sub_masked(self):
        assert "$x" not in h.mask_opaque("for f in $(echo $x); do :; done")

    def test_arith_sub_masked(self):
        assert "$" not in h.mask_opaque("echo $((1 + 2))")

    def test_backtick_masked(self):
        assert "$x" not in h.mask_opaque("echo `echo $x`")

    def test_unquoted_scalar_survives(self):
        assert "$x" in h.mask_opaque("for f in $x; do :; done")

    def test_length_preserved(self):
        cmd = 'a "b c" $(d) e'
        assert len(h.mask_opaque(cmd)) == len(cmd)


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
    def test_wordsplit_denied(self):
        proc = _run(_bash("for x in $var; do echo $x; done"))
        assert proc.returncode == 0
        out = json.loads(proc.stdout)
        assert out["hookSpecificOutput"]["permissionDecision"] == "deny"

    def test_correct_form_allowed(self):
        proc = _run(_bash('for f in "${arr[@]}"; do echo "$f"; done'))
        assert proc.returncode == 0
        assert proc.stdout.strip() == ""

    def test_escape_hatch_allows(self):
        env = dict(os.environ, CLAUDE_ALLOW_WORDSPLIT="1")
        proc = _run(_bash("set -- $pair"), env=env)
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
        proc = _run({"tool_name": "Edit", "tool_input": {"command": "for x in $v; do :; done"}})
        assert proc.stdout.strip() == ""


def test_hook_is_executable():
    # The harness execs the hook by bare path, so a non-executable file (100644)
    # EACCESes before the shebang runs - the blocker caught on PR #1029. The
    # subprocess tests above invoke via `sys.executable <hook>`, which bypasses
    # the exec bit, so this is the one check that exercises the real harness path.
    assert os.access(HOOK, os.X_OK), f"{HOOK} must be committed executable (100755)"
