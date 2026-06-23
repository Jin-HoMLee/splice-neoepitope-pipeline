"""Tests for scripts/new_branch.sh — canonical branch-name helper (Issue #578).

The script derives <type>/<role>/issue-N-slug from an Issue and wraps
`gh issue develop`. For deterministic offline tests we put a stub `gh` first
on PATH that, for `gh api graphql ...`, prints canned GraphQL JSON read from
the file named by $NEW_BRANCH_GH_FIXTURE. Every case runs against --dry-run,
so the mutating path (`gh issue develop` / `git checkout -b`) is never invoked.

Spec: docs/superpowers/specs/2026-06-22-new-branch-helper-design.md
"""
import json
import os
import subprocess
import tempfile
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "new_branch.sh"

# A stub `gh`: for `api graphql` it cats $NEW_BRANCH_GH_FIXTURE; anything else errors
# (no non-graphql gh call should happen under --dry-run or on a clean error path).
_GH_STUB = """#!/usr/bin/env bash
if [ "$1" = "api" ] && [ "$2" = "graphql" ]; then
  cat "$NEW_BRANCH_GH_FIXTURE"
  exit 0
fi
echo "stub gh: unexpected call: $*" >&2
exit 99
"""

# A stub `git`: no test should reach a real git mutation. Any call is a failure
# signal (e.g. a swallowed --dry-run falling through to git fetch / checkout -b).
_GIT_STUB = """#!/usr/bin/env bash
echo "stub git: unexpected call: $*" >&2
exit 98
"""


def _graphql_json(title, roles, sub_total):
    """Build the GraphQL response envelope the script consumes."""
    return {
        "data": {
            "repository": {
                "issue": {
                    "title": title,
                    "labels": {"nodes": [{"name": n} for n in roles]},
                    "subIssuesSummary": {"total": sub_total},
                }
            }
        }
    }


def _run(args, *, title="feat(scripts): x", roles=("role:pm",), sub_total=0, raw=None):
    """Run new_branch.sh with stub gh + git on PATH returning a crafted issue.

    raw: if given, the gh-graphql fixture is this literal JSON string (used to
    inject edge envelopes like a null issue) instead of the built envelope.
    """
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        for name, body in (("gh", _GH_STUB), ("git", _GIT_STUB)):
            p = d / name
            p.write_text(body)
            p.chmod(0o755)
        fix = d / "issue.json"
        fix.write_text(raw if raw is not None
                       else json.dumps(_graphql_json(title, list(roles), sub_total)))
        env = {
            **os.environ,
            "PATH": f"{d}:{os.environ['PATH']}",
            "NEW_BRANCH_GH_FIXTURE": str(fix),
        }
        return subprocess.run(
            ["bash", str(SCRIPT), *args],
            env=env, capture_output=True, text=True,
        )


# --- Case 1: title-type parse → feat/pm/issue-578-branch-helper -------------

def test_title_type_parse_happy_path():
    r = _run(["578", "branch-helper", "--dry-run"],
             title="feat(scripts): branch-naming helper", roles=("role:pm",))
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "feat/pm/issue-578-branch-helper"


# --- Case 2: --type override -------------------------------------------------

def test_type_override():
    r = _run(["578", "branch-helper", "--type", "spike", "--dry-run"],
             title="feat(scripts): branch-naming helper", roles=("role:pm",))
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "spike/pm/issue-578-branch-helper"


# --- Case 3: multi-role issue requires --role --------------------------------

def test_multi_role_without_override_errors_naming_both():
    r = _run(["578", "branch-helper", "--dry-run"],
             title="feat(scripts): x", roles=("role:pm", "role:developer"))
    assert r.returncode == 1
    assert "pm" in r.stderr and "developer" in r.stderr


def test_multi_role_with_override_succeeds():
    r = _run(["578", "branch-helper", "--role", "pm", "--dry-run"],
             title="feat(scripts): x", roles=("role:pm", "role:developer"))
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "feat/pm/issue-578-branch-helper"


def test_no_role_label_without_override_errors():
    r = _run(["578", "branch-helper", "--dry-run"],
             title="feat(scripts): x", roles=())
    assert r.returncode == 1
    assert "role" in r.stderr.lower()


# --- Case 4: slug sanitization -----------------------------------------------

def test_slug_sanitization_spaces():
    r = _run(["578", "Branch Helper", "--dry-run"], roles=("role:pm",))
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "feat/pm/issue-578-branch-helper"


def test_slug_sanitization_unicode_arrow():
    r = _run(["578", "west1-b → west4-a", "--dry-run"], roles=("role:pm",))
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "feat/pm/issue-578-west1-b-west4-a"


def test_slug_sanitization_underscores_and_case():
    r = _run(["578", "P100__GPU", "--dry-run"], roles=("role:pm",))
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "feat/pm/issue-578-p100-gpu"


# --- Case 5: missing slug → error printing the title -------------------------

def test_missing_slug_errors_with_title():
    r = _run(["578", "--dry-run"],
             title="feat(scripts): branch-naming helper", roles=("role:pm",))
    assert r.returncode == 1
    assert "branch-naming helper" in r.stderr


# --- Case 6: no title type + no --type → error -------------------------------

def test_no_title_type_no_override_errors():
    r = _run(["578", "branch-helper", "--dry-run"],
             title="branch-naming helper with no CC prefix", roles=("role:pm",))
    assert r.returncode == 1
    assert "type" in r.stderr.lower()


# --- Case 7: parent issue → parent-refuse ------------------------------------

def test_parent_issue_refused():
    r = _run(["538", "branch-helper", "--dry-run"],
             title="feat(scripts): x", roles=("role:pm",), sub_total=3)
    assert r.returncode == 1
    assert "parent" in r.stderr.lower() or "epic" in r.stderr.lower()


# --- Case 8: --no-issue fallback ---------------------------------------------

def test_no_issue_fallback():
    r = _run(["--no-issue", "labnotebook", "pm", "memory-slim", "--dry-run"])
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "labnotebook/pm/memory-slim"


# --- Review 🟡#1: an option flag must not swallow the following flag ---------

def test_type_flag_without_value_rejected_no_mutation():
    # `--type --dry-run`: --dry-run must NOT be consumed as the type value
    # (that would clear the dry-run guard and attempt a real branch-create).
    r = _run(["578", "branch-helper", "--type", "--dry-run"], roles=("role:pm",))
    assert r.returncode == 1, r.stdout + r.stderr
    assert "--type" in r.stderr and "value" in r.stderr.lower()
    # the git stub exits 98 if a mutation was attempted — it must not have been
    assert "stub git" not in r.stderr


def test_role_flag_without_value_rejected():
    r = _run(["578", "branch-helper", "--role", "--dry-run"], roles=("role:pm",))
    assert r.returncode == 1, r.stdout + r.stderr
    assert "--role" in r.stderr and "value" in r.stderr.lower()


# --- Review 🟡#2: non-existent issue → clean error, not a raw jq trace -------

def test_nonexistent_issue_clean_error():
    r = _run(["999999", "branch-helper", "--dry-run"],
             raw='{"data":{"repository":{"issue":null}}}')
    assert r.returncode == 1
    assert "999999" in r.stderr and "not found" in r.stderr.lower()
    assert "Cannot iterate over null" not in r.stderr  # no raw jq leak


# --- Review nit: --no-issue sanitizes type + role ---------------------------

def test_no_issue_sanitizes_type_and_role():
    r = _run(["--no-issue", "Feat", "PM", "memory-slim", "--dry-run"])
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "feat/pm/memory-slim"


def test_no_issue_rejects_type_flag():
    # --type/--role are meaningless with --no-issue (type+role are positional);
    # silently consuming them then erroring on "missing slug" is confusing.
    r = _run(["--no-issue", "--type", "feat", "pm", "memory-slim", "--dry-run"])
    assert r.returncode == 1
    assert "--type" in r.stderr and "no-issue" in r.stderr.lower()


# --- Spec error-table gaps (lines 80-89) ------------------------------------

def test_non_integer_issue_errors():
    r = _run(["abc", "branch-helper", "--dry-run"])
    assert r.returncode == 1
    assert "integer" in r.stderr.lower()


def test_slug_sanitizes_to_empty_errors():
    r = _run(["578", "→→→", "--dry-run"], roles=("role:pm",))
    assert r.returncode == 1
    assert "empty" in r.stderr.lower()


def test_no_issue_missing_slug_errors():
    r = _run(["--no-issue", "labnotebook", "pm", "--dry-run"])
    assert r.returncode == 1
    assert "no-issue" in r.stderr.lower() or "slug" in r.stderr.lower()
