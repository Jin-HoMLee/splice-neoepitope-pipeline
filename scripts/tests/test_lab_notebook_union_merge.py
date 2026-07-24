"""Drive the real two-branch lab-notebook merge (Issue #1221, Option B).

The fix is one line in `.gitattributes`:

    research/lab_notebook/*.md merge=union

`merge=union` is git's built-in "keep both sides' lines" driver, so two
concurrent same-role PRs that each append a `### <time>` entry at the same
top-of-today anchor no longer produce a *blocking* conflict on merge.

These tests DRIVE that with real `git merge` invocations in throwaway repos - the
AC is explicit that verification must be the real case, not a synthetic claim.
The load-bearing pair is the matched-pair control: WITH the union attribute a
divergent-append merge auto-resolves and keeps both entries; WITHOUT it the exact
same merge conflicts. If the "without" half did not conflict, the "with" half
would prove nothing (the attribute would not be doing the work).
"""

import os
import re
import subprocess
from pathlib import Path

_PATH = os.environ.get("PATH", "/usr/bin:/bin")

REPO_ROOT = Path(__file__).resolve().parents[2]
GITATTRIBUTES = REPO_ROOT / ".gitattributes"

# The committed rule, so the test verifies OUR file, not union-in-the-abstract.
_UNION_RULE_RE = re.compile(
    r"^\s*research/lab_notebook/\*\.md\s+merge=union\s*$", re.MULTILINE
)

NOTEBOOK_REL = "research/lab_notebook/developer.md"
BASE = "# Lab Notebook - Developer\n\n## 2026-07-24\n"
ENTRY_A = "### 10:00 - PR A\n\nfirst concurrent entry.\n\n"
ENTRY_B = "### 11:00 - PR B\n\nsecond concurrent entry.\n\n"
ANCHOR = "## 2026-07-24\n"


def _git(repo: Path, *args, check=True):
    return subprocess.run(
        ["git", *args],
        cwd=repo,
        capture_output=True,
        text=True,
        check=check,
        env={"GIT_CONFIG_GLOBAL": "/dev/null", "GIT_CONFIG_SYSTEM": "/dev/null",
             "HOME": str(repo), "PATH": _PATH},
    )


def _init_repo(repo: Path, *, with_union: bool):
    _git(repo, "init", "-q", "-b", "main")
    _git(repo, "config", "user.email", "t@example.com")
    _git(repo, "config", "user.name", "T")
    _git(repo, "config", "commit.gpgsign", "false")
    if with_union:
        (repo / ".gitattributes").write_text(
            "research/lab_notebook/*.md merge=union\n", encoding="utf-8"
        )
        _git(repo, "add", ".gitattributes")
    nb = repo / NOTEBOOK_REL
    nb.parent.mkdir(parents=True, exist_ok=True)
    nb.write_text(BASE, encoding="utf-8")
    _git(repo, "add", NOTEBOOK_REL)
    _git(repo, "commit", "-qm", "base")


def _prepend_entry(repo: Path, entry: str):
    """Insert an entry right under the top-of-today anchor (the real hot spot)."""
    nb = repo / NOTEBOOK_REL
    text = nb.read_text(encoding="utf-8")
    text = text.replace(ANCHOR, ANCHOR + "\n" + entry, 1)
    nb.write_text(text, encoding="utf-8")


def _two_branch_append_merge(repo: Path):
    """Branch A and B each append at the same anchor; merge B into A.

    Returns the CompletedProcess of the `git merge` (check=False so the caller
    inspects success/conflict), with the notebook left in its post-merge state.
    """
    _git(repo, "checkout", "-qb", "branchA")
    _prepend_entry(repo, ENTRY_A)
    _git(repo, "commit", "-aqm", "A")

    _git(repo, "checkout", "-q", "main")
    _git(repo, "checkout", "-qb", "branchB")
    _prepend_entry(repo, ENTRY_B)
    _git(repo, "commit", "-aqm", "B")

    _git(repo, "checkout", "-q", "branchA")
    return _git(repo, "merge", "--no-edit", "branchB", check=False)


def test_committed_gitattributes_has_the_union_rule():
    """The repo's own .gitattributes carries the Option-B rule."""
    assert GITATTRIBUTES.exists(), ".gitattributes missing"
    content = GITATTRIBUTES.read_text(encoding="utf-8")
    assert _UNION_RULE_RE.search(content), (
        "research/lab_notebook/*.md merge=union not found in .gitattributes"
    )


def test_union_merge_auto_resolves_and_keeps_both(tmp_path):
    """WITH the union attribute: the divergent-append merge succeeds, no conflict
    markers, and BOTH entries survive."""
    repo = tmp_path / "with_union"
    repo.mkdir()
    _init_repo(repo, with_union=True)
    merge = _two_branch_append_merge(repo)

    assert merge.returncode == 0, (
        f"union merge should not conflict, but git merge failed:\n{merge.stdout}\n{merge.stderr}"
    )
    final = (repo / NOTEBOOK_REL).read_text(encoding="utf-8")
    assert "<<<<<<<" not in final and ">>>>>>>" not in final, "conflict markers present"
    assert "PR A" in final and "PR B" in final, "union merge dropped an entry"


def test_without_union_the_same_merge_conflicts(tmp_path):
    """Matched-pair CONTROL: WITHOUT the attribute, the identical merge conflicts.

    This is the falsifier - it proves the union attribute is what resolves the
    merge in the test above, not some artifact of the setup."""
    repo = tmp_path / "no_union"
    repo.mkdir()
    _init_repo(repo, with_union=False)
    merge = _two_branch_append_merge(repo)

    assert merge.returncode != 0, "expected a conflict without merge=union, got a clean merge"
    final = (repo / NOTEBOOK_REL).read_text(encoding="utf-8")
    assert "<<<<<<<" in final, "expected conflict markers in the un-attributed merge"
