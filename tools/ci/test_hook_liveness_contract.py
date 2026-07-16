"""Hook-liveness contract suite (Issue #1140, Sub D of #1135).

Every hook registered in `.agents/settings.json` gets a contract that drives its
REAL trigger through its REAL entry path (the hook script as a subprocess reading
the harness JSON envelope on stdin) and asserts its observable artifact appears -
a `permissionDecision: deny` on stdout, an appended `.agents/hook_fires.jsonl`
line, or a written watermark file. A matched-pair counterexample asserts the hook
stays silent when it should.

This exists because four governance hooks shipped completely inert this month,
each green on its unit tests, none ever invoked through its real path. The
unit tests called inner functions (`matches_pr_create`) or piped a synthetic
payload into `main()`; both bypass the layer where the hooks actually failed.
Mutation testing (Issue #1141) is structurally blind to this class - a hook that
is never invoked has no surviving mutant. This suite is the half that catches it.

The registry below is populated per hook by the analyze+verify fan-out of
Issue #1140; `test_every_registered_hook_has_a_contract` fails if the registry
and `.agents/settings.json` ever disagree, so a newly registered hook cannot be
silently uncovered.
"""
import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hook_contract as hc  # noqa: E402

# --- gh stubs for hooks that reach a live `gh` call before their observable ---
# Each returns just enough canned JSON for the hook to run its matcher-plus-handler
# path to the point of appending its fire-log line, with no network. Mirrors the
# `test_set_status.py` PATH-stub pattern; validated per hook by the #1140 fan-out.

_STUB_GH_ISSUE_DEVELOP_PARENT = """#!/bin/sh
case "$*" in
  *"repo view"*) echo "Jin-HoMLee/splice-neoepitope-pipeline" ;;
  *"api graphql"*) echo '{"data":{"repository":{"issue":{"subIssuesSummary":{"total":5}}}}}' ;;
esac
exit 0
"""

_STUB_PR_CREATE = """#!/bin/sh
case "$*" in
  *"pr view"*) echo '{"isDraft": true, "body": "draft body"}' ;;
  *"project item-add"*) echo '{"id": "ITEM_1"}' ;;
  *) echo '{}' ;;
esac
exit 0
"""

_STUB_PR_REVIEW_REQUEST = """#!/bin/bash
if [ "$1" = "pr" ] && [ "$2" = "view" ]; then
  echo '{"url":"https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/996","closingIssuesReferences":[{"number":123}]}'
  exit 0
fi
if [ "$1" = "api" ] && [ "$2" = "graphql" ]; then
  for a in "$@"; do
    case "$a" in
      *updateProjectV2ItemFieldValue*) echo '{"data":{"updateProjectV2ItemFieldValue":{"projectV2Item":{"id":"x"}}}}'; exit 0;;
    esac
  done
  echo '{"data":{"repository":{"pullRequest":{"projectItems":{"nodes":[{"id":"ITEM1","project":{"number":9},"fieldValues":{"nodes":[{"name":"Ready for review","field":{"name":"Status"}}]}}]}},"issue":{"projectItems":{"nodes":[{"id":"ITEM2","project":{"number":9},"fieldValues":{"nodes":[{"name":"Ready for review","field":{"name":"Status"}}]}}]}}}}}'
  exit 0
fi
exit 0
"""

_STUB_RECHECK_DISPATCH = """#!/bin/sh
args="$*"
case "$args" in
  *"repos/Jin-HoMLee/splice-neoepitope-pipeline/issues/"*)
    printf '%s' '{"title":"i5 - S3 - Data Preparation","due_on":"2026-08-01T00:00:00Z"}'
    ;;
  *updateProjectV2ItemFieldValue*)
    printf '%s' '{"errors":[{"message":"stub forced failure"}]}'
    ;;
  *projectItems*)
    printf '%s' '{"data":{"repository":{"issue":{"projectItems":{"nodes":[{"id":"PVTI_stubitem123","project":{"number":9},"fieldValues":{"nodes":[{"date":"2026-07-01","field":{"name":"Target date"}}]}}]}}}}}'
    ;;
  *)
    printf '%s' '{}'
    ;;
esac
exit 0
"""

_PR_URL = "https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/1140\n"
# cwd values for the drift hook, computed at runtime so the contract is portable:
# a cwd above the clone fires, the clone root itself stays silent.
_OUTSIDE_CLONE = str(Path(hc.REPO_ROOT).parent)
_INSIDE_CLONE = str(hc.REPO_ROOT)

# --- the registry ------------------------------------------------------------
# basename -> HookContract. `fire_input`/`nofire_input` are kwargs for the named
# envelope builder; `heredoc_fire_input` (when set) is the heredoc-then-command
# shape that killed matches_pr_create; `gh_stub_source` is set for board hooks
# that reach a live `gh` call before their fire-log observable.
#
# POPULATED FROM THE #1140 ANALYZE+VERIFY FAN-OUT. Until every registered hook is
# present, `test_every_registered_hook_has_a_contract` is red by design.

# The em-dash the check_no_emdash fixture must carry at runtime, built with
# chr() so this source file stays pure ASCII (and does not itself trip the guard).
_EMDASH = chr(0x2014)

CONTRACTS: dict[str, hc.HookContract] = {
    "check_gh_issue_develop_parent.py": hc.HookContract(
        basename="check_gh_issue_develop_parent.py",
        observable=hc.DENY,
        envelope_builder="pretooluse_bash",
        fire_input={"command": "gh issue develop 538 --name feat/x --checkout"},
        nofire_input={"command": "gh issue view 538"},
        heredoc_fire_input={
            "command": "cat > /tmp/plan.md <<'EOF'\nbranching notes for the epic\nEOF\ngh issue develop 538 --name feat/x --checkout"
        },
        gh_stub_source=_STUB_GH_ISSUE_DEVELOP_PARENT,
        notes="denies branching off a parent/epic (subIssuesSummary.total > 0)",
    ),
    "check_gh_issue_create_repo.py": hc.HookContract(
        basename="check_gh_issue_create_repo.py",
        observable=hc.FIRE_LOG,
        envelope_builder="pretooluse_bash",
        fire_input={
            "command": 'gh issue create --title "drain episodic memory" --label role:memory_manager --repo Jin-HoMLee/splice-neoepitope-pipeline'
        },
        nofire_input={
            "command": 'gh issue create --title "drain episodic memory" --label role:memory_manager --repo Jin-HoMLee/claude-personas-splice-neoepitope-pipeline'
        },
        heredoc_fire_input={
            "command": "cat > /tmp/note.md <<'EOF'\nfiling notes\nEOF\ngh issue create --title \"drain episodic memory\" --label role:memory_manager --repo Jin-HoMLee/splice-neoepitope-pipeline"
        },
        notes="asks (fire-log) on an MM-shaped create into the project repo; the same shape into personas passes (explicit --repo, no gh call)",
    ),
    "check_board_query_pagination.py": hc.HookContract(
        basename="check_board_query_pagination.py",
        observable=hc.DENY,
        envelope_builder="pretooluse_bash",
        fire_input={
            "command": "gh api graphql -f query='query { user(login:\"Jin-HoMLee\"){ projectV2(number:9){ items(first:100){ nodes { id } } } } }'"
        },
        nofire_input={
            "command": "gh api graphql -f query='query { user(login:\"Jin-HoMLee\"){ projectV2(number:9){ items(first:100){ pageInfo{ hasNextPage endCursor } nodes { id } } } } }'"
        },
        notes="denies an unpaginated projectV2 items(first:) board query (pure string inspection, no gh call)",
    ),
    "check_no_force_push.py": hc.HookContract(
        basename="check_no_force_push.py",
        observable=hc.DENY,
        envelope_builder="pretooluse_bash",
        fire_input={"command": "git push --force origin main"},
        nofire_input={"command": "git push origin main"},
        notes="denies a force push; a plain push passes",
    ),
    "check_commit_push_separation.py": hc.HookContract(
        basename="check_commit_push_separation.py",
        observable=hc.DENY,
        envelope_builder="pretooluse_bash",
        fire_input={"command": "git commit -m 'x' && git push"},
        nofire_input={"command": "git push origin main"},
        notes="denies a chained commit&&push; a lone push passes",
    ),
    "check_no_cd_outside_cwd.py": hc.HookContract(
        basename="check_no_cd_outside_cwd.py",
        observable=hc.DENY,
        envelope_builder="pretooluse_bash",
        fire_input={"command": "cd /etc/secret"},
        nofire_input={"command": "cd workflow"},
        notes="denies a cd outside the cwd subtree; a cd into a subdir passes",
    ),
    "check_memory_path_cwd_drift.py": hc.HookContract(
        basename="check_memory_path_cwd_drift.py",
        observable=hc.DENY,
        envelope_builder="pretooluse_bash",
        fire_input={
            "command": "grep foo .agents/memory/shared/MEMORY.md",
            "cwd": _OUTSIDE_CLONE,
        },
        nofire_input={
            "command": "grep foo .agents/memory/shared/MEMORY.md",
            "cwd": _INSIDE_CLONE,
        },
        notes="denies a relative memory-path op when the session cwd has drifted outside the clone",
    ),
    "check_no_emdash.py": hc.HookContract(
        basename="check_no_emdash.py",
        observable=hc.DENY,
        envelope_builder="pretooluse_edit",
        fire_input={
            "file_path": "research/notes.md",
            "old_string": "hello - world",
            "new_string": "hello " + _EMDASH + " world",
        },
        nofire_input={
            "file_path": "research/notes.md",
            "old_string": "hello - world",
            "new_string": "hello -- world",
        },
        notes="denies an edit that net-adds an em/en-dash; a plain-hyphen edit passes",
    ),
    "post_gh_pr_create.py": hc.HookContract(
        basename="post_gh_pr_create.py",
        observable=hc.FIRE_LOG,
        envelope_builder="posttooluse_bash",
        fire_input={
            "command": 'gh pr create --title "t" --body "closes #1140" --base main',
            "tool_response": {"stdout": _PR_URL},
        },
        nofire_input={
            "command": "gh pr view 1140",
            "tool_response": {"stdout": _PR_URL},
        },
        heredoc_fire_input={
            "command": "cat > /tmp/pr_1140.md <<'EOF'\nPR body with (parens) and prose.\nEOF\ngh pr create --title \"t\" --body-file /tmp/pr_1140.md --base main",
            "tool_response": {"stdout": _PR_URL},
        },
        gh_stub_source=_STUB_PR_CREATE,
        notes="board-adds a freshly created PR (the hook that sat dead for months on the heredoc shape)",
    ),
    "post_gh_pr_review_request.py": hc.HookContract(
        basename="post_gh_pr_review_request.py",
        observable=hc.FIRE_LOG,
        envelope_builder="posttooluse_bash",
        fire_input={
            "command": 'gh pr comment 996 --body "@claude review"',
            "tool_response": {"stdout": ""},
        },
        nofire_input={
            "command": 'gh pr comment 996 --body "@-claude review"',
            "tool_response": {"stdout": ""},
        },
        gh_stub_source=_STUB_PR_REVIEW_REQUEST,
        notes="advances a PR + its linked Issues to In review on the real review trigger; the hyphenated reference form must not fire",
    ),
    "recheck_dispatch.py": hc.HookContract(
        basename="recheck_dispatch.py",
        observable=hc.FIRE_LOG,
        envelope_builder="posttooluse_bash",
        fire_input={
            "command": 'gh issue edit 42 --milestone "i5 - S3 - Data Preparation"',
            "tool_response": {"stdout": ""},
        },
        nofire_input={
            "command": 'gh issue edit 42 --add-label "role:developer"',
            "tool_response": {"stdout": ""},
        },
        gh_stub_source=_STUB_RECHECK_DISPATCH,
        extra_args=["--scope", "shared"],
        notes="runs the milestone->Target-date sync on a milestone edit; a label edit does not trigger it",
    ),
    "write_session_watermark.py": hc.HookContract(
        basename="write_session_watermark.py",
        observable=hc.WATERMARK,
        envelope_builder="stop_envelope",
        # No discriminating matcher: writes on any writable root. The matched pair
        # is a writable root (writes) vs an unwritable root (fails open, silent).
        fire_input={"root": hc.WRITABLE_ROOT},
        nofire_input={"root": "/dev/null/sub"},
        notes="writes the session watermark on Stop; an unwritable root fails open with no marker",
    ),
}


# --- completeness / drift ----------------------------------------------------


def test_every_registered_hook_has_a_contract():
    """The registry and the settings.json registry must agree exactly.

    A registered-but-uncovered hook (drift toward an inert hook shipping unseen)
    and a contract for a hook no longer registered both fail here.
    """
    registered = hc.distinct_hook_basenames()
    covered = set(CONTRACTS)
    missing = registered - covered
    stale = covered - registered
    assert not missing, f"registered hooks with no liveness contract: {sorted(missing)}"
    assert not stale, f"contracts for unregistered hooks: {sorted(stale)}"


def test_every_contract_hook_script_exists():
    for name, c in CONTRACTS.items():
        assert c.hook_path.exists(), f"{name}: hook script missing at {c.hook_path}"


# --- the matched pair, per hook ----------------------------------------------


@pytest.mark.parametrize("name", sorted(CONTRACTS), ids=lambda n: n.replace(".py", ""))
def test_hook_fires_on_real_trigger(name):
    """The real trigger, driven through the real entry path, produces the artifact."""
    c = CONTRACTS[name]
    assert hc.observe(c, c.fire_input) is True, (
        f"{name}: real trigger produced no {c.observable} artifact - the hook is inert "
        f"on its own trigger (fire_input={c.fire_input})"
    )


@pytest.mark.parametrize("name", sorted(CONTRACTS), ids=lambda n: n.replace(".py", ""))
def test_hook_silent_on_counterexample(name):
    """A real, close counterexample produces no artifact (the matched-pair control)."""
    c = CONTRACTS[name]
    assert hc.observe(c, c.nofire_input) is False, (
        f"{name}: counterexample wrongly produced a {c.observable} artifact - the "
        f"contract cannot distinguish fire from no-fire (nofire_input={c.nofire_input})"
    )


# --- the shape that killed two hooks -----------------------------------------


@pytest.mark.parametrize(
    "name",
    sorted(n for n, c in CONTRACTS.items() if c.heredoc_fire_input is not None),
    ids=lambda n: n.replace(".py", ""),
)
def test_heredoc_shape_still_fires(name):
    """A heredoc-then-command trigger fires - the exact string that left
    `matches_pr_create` dead for months while every unit test passed."""
    c = CONTRACTS[name]
    assert hc.observe(c, c.heredoc_fire_input) is True, (
        f"{name}: the heredoc-then-command shape produced no {c.observable} artifact - "
        f"this is the shape that shipped two hooks inert (input={c.heredoc_fire_input})"
    )


# --- red-on-break: demonstrate the suite can fail ----------------------------


def test_breaking_a_matcher_turns_the_contract_red(tmp_path):
    """DEMONSTRATION (not assertion): neutralize `check_no_force_push`'s matcher in
    an isolated copy, drive its real force-push trigger through the copy, and confirm
    the deny observable DISAPPEARS. If breaking the matcher could not silence the
    hook, the whole suite would be a hollow check about hollow checks.

    The copy lives in a pytest tmp_path (never the tracked hooks dir), with
    `_shell_parse.py` copied alongside so the hook's relative import resolves - a
    lone copy would fail to import and vanish from an import CRASH, not the matcher,
    which is a *hollow* red (this exact trap surfaced during the #1140 fan-out)."""
    fire_env = hc.pretooluse_bash("git push --force origin main")

    # precondition: the real hook denies this force push
    real = hc.drive(hc.HOOKS_DIR / "check_no_force_push.py", fire_env)
    assert hc.deny_decision(real) is True, "precondition: real hook must deny a force push"

    # build an isolated broken copy: matcher neutralized to always return False
    src = (hc.HOOKS_DIR / "check_no_force_push.py").read_text(encoding="utf-8")
    marker = "def command_force_pushes(cmd: str) -> bool:"
    assert marker in src, "matcher signature changed; update the red-on-break demo"
    broken = src.replace(
        marker, marker + "\n    return False  # red-on-break: matcher neutralized", 1
    )
    assert broken != src

    broken_hook = tmp_path / "check_no_force_push.py"
    broken_hook.write_text(broken, encoding="utf-8")
    shutil.copy2(hc.HOOKS_DIR / "_shell_parse.py", tmp_path / "_shell_parse.py")

    # unbroken control: the same isolated copy (before breaking) must still fire,
    # so a silenced deny is attributable to the matcher break, not the relocation.
    control_hook = tmp_path / "control.py"
    control_hook.write_text(src, encoding="utf-8")
    control = hc.drive(control_hook, fire_env)
    assert hc.deny_decision(control) is True, (
        "the relocated but UNbroken copy must still deny - otherwise a silenced deny "
        "would be an artifact of the move (e.g. a failed import), a hollow red"
    )

    broken_result = hc.drive(broken_hook, fire_env)
    assert hc.deny_decision(broken_result) is False, (
        "breaking the matcher did NOT silence the hook - the contract cannot come "
        "back red, so it is a hollow check"
    )
