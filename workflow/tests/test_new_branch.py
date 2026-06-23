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
# (no non-graphql gh call should happen under --dry-run).
_GH_STUB = """#!/usr/bin/env bash
if [ "$1" = "api" ] && [ "$2" = "graphql" ]; then
  cat "$NEW_BRANCH_GH_FIXTURE"
  exit 0
fi
echo "stub gh: unexpected call: $*" >&2
exit 99
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


def _run(args, *, title="feat(scripts): x", roles=("role:pm",), sub_total=0):
    """Run new_branch.sh with a stub gh on PATH returning a crafted issue."""
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        # stub gh
        gh = d / "gh"
        gh.write_text(_GH_STUB)
        gh.chmod(0o755)
        # fixture JSON
        fix = d / "issue.json"
        fix.write_text(json.dumps(_graphql_json(title, list(roles), sub_total)))
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
