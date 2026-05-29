"""Post-merge / post-close closure-audit critic.

Spec: docs/superpowers/specs/2026-05-11-post-merge-critic-design.md

Trying Sakana Fugu "critic-as-default" on PR-merge / issue-close. Posts a
single marker-tagged gap comment, or stays silent if all 3 checks pass.
No edit-in-place — comment is a post-only snapshot of close-time state.

Usage (called by .github/workflows/closure-audit.yml):
    python closure_audit.py --event-type {pr|issue} --number <N>
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from collections.abc import Collection
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
COMMENT_MARKER = "<!-- closure-audit -->"

# Matches only `- [ ]` / `- [x]` style — `*` and `+` bullet markers are not
# recognized. The repo's issue templates use `-` throughout; revisit if other
# styles show up.
_UNTICKED = re.compile(r"^\s*-\s*\[\s\]\s", re.MULTILINE)
_TICKED = re.compile(r"^\s*-\s*\[[xX]\]\s", re.MULTILINE)
_EXEMPT_FILES = {"research/glossary.md"}
_EXEMPT_PREFIX = ("research/lab_notebook/",)


# --- pure checks ---


def check_ac(body: str, comments: list[str]) -> str | None:
    unticked = len(_UNTICKED.findall(body))
    ticked = len(_TICKED.findall(body))
    if unticked == 0:
        return None
    if any("❎" in c and "deferred" in c.lower() for c in comments):
        return None
    return f"{unticked}/{unticked + ticked} unticked, no deferral comment found"


def check_priority_rationale(body: str) -> str | None:
    if "priority rationale" in body.lower():
        return None
    return "no 'Priority rationale' line in issue body"


def check_lab_notebook(
    text: str, date: str, number: int, also_accept: Collection[int] = ()
) -> str | None:
    """Gap unless the `## date` block references `#number` OR any `#also_accept`.

    `also_accept` carries the PR's closing-Issue numbers (#495): entries written
    before the PR exists reference the Issue, not the PR number — both are
    semantically the same unit of work.
    """
    header = f"## {date}"
    if header not in text:
        return f"no '## {date}' header in notebook"
    start = text.index(header)
    rest = text[start + len(header):]
    nxt = re.search(r"^## ", rest, re.MULTILINE)
    block = rest[:nxt.start()] if nxt else rest
    if not re.search(r"^### ", block, re.MULTILINE):
        return f"'## {date}' has no '### HH:MM UTC — Editor: …' sub-section"
    accepted = [number, *also_accept]
    if not any(f"#{n}" in block for n in accepted):
        refs = " or ".join(f"'#{n}'" for n in accepted)
        return f"'## {date}' block does not reference {refs}"
    return None


def is_exempt(changed_files: list[str]) -> bool:
    return bool(changed_files) and all(
        p in _EXEMPT_FILES or p.startswith(_EXEMPT_PREFIX)
        for p in changed_files
    )


def resolve_roles(labels_per_issue: list[list[str]]) -> list[set[str]]:
    """Per-Issue role sets, positionally aligned with the input list.

    Multi-role Issues keep ALL their roles (see #524). Issues with no `role:`
    label produce an empty set — callers should skip those entries.
    """
    return [
        {lbl[len("role:"):] for lbl in labels if lbl.startswith("role:")}
        for labels in labels_per_issue
    ]


def check_lab_notebooks_for_issue(
    roles: set[str],
    date: str,
    number: int,
    notebooks: dict[str, str | None],
    also_accept: Collection[int] = (),
) -> tuple[str, str] | None:
    """Any-role-satisfies notebook check for one Issue's role set.

    Returns None if at least one role's notebook references `#number` (or any
    `#also_accept`) in the `## date` block. Otherwise returns a single
    (role_label, description) tuple summarizing the gap across all roles.

    `notebooks` maps role → notebook text (None for missing file).
    """
    per_role_gaps: list[tuple[str, str]] = []
    for role in sorted(roles):
        text = notebooks.get(role)
        if text is None:
            per_role_gaps.append((role, "lab notebook file missing"))
            continue
        if (gap := check_lab_notebook(text, date, number, also_accept)) is None:
            return None
        per_role_gaps.append((role, gap))
    if len(per_role_gaps) == 1:
        return per_role_gaps[0]
    roles_str = " or ".join(r for r, _ in per_role_gaps)
    descs = "; ".join(f"{r}: {d}" for r, d in per_role_gaps)
    return (roles_str, descs)


def collect_notebook_gaps(
    role_sets_per_issue: list[set[str]],
    date: str,
    number: int,
    notebooks: dict[str, str | None],
    also_accept: Collection[int] = (),
) -> list[tuple[str, str]]:
    """Aggregate notebook gaps across closing Issues, deduped by role set.

    Two Issues with the same role set produce identical gap entries (the
    underlying check is deterministic on role-set/date/number/also_accept/
    notebooks — and `also_accept` is constant across one audit), so only the
    first is emitted.
    """
    gaps: list[tuple[str, str]] = []
    seen: set[frozenset[str]] = set()
    for roles in role_sets_per_issue:
        if not roles:
            continue
        key = frozenset(roles)
        if key in seen:
            continue
        seen.add(key)
        if gap := check_lab_notebooks_for_issue(
            roles, date, number, notebooks, also_accept
        ):
            gaps.append(gap)
    return gaps


def format_comment(
    event_label: str,
    ac_gaps: list[tuple[int, str]],
    pr_gaps: list[tuple[int, str]],
    nb_gaps: list[tuple[str, str]],
) -> str:
    if not (ac_gaps or pr_gaps or nb_gaps):
        return (
            f"{COMMENT_MARKER}\n"
            f"✅ Closure audit — all clear\n\n"
            f"_— closure-audit bot (experiment, see #325)_\n"
        )
    lines = [
        COMMENT_MARKER,
        "## Closure audit — gaps found",
        "",
        f"{event_label} closed with the following items still open:",
        "",
    ]
    for n, d in ac_gaps:
        lines.append(f"- ⚠️ **AC checkboxes** on Issue #{n} — {d}")
    for n, d in pr_gaps:
        lines.append(f"- ⚠️ **Priority rationale** on Issue #{n} — {d}")
    for r, d in nb_gaps:
        lines.append(f"- ⚠️ **Lab notebook entry** for `{r}` — {d}")
    lines += [
        "",
        "Per [closure ritual](shared/feedback_closure_ritual.md): tick the boxes, add a comment-deferral, or add the missing entry/rationale. This comment is a snapshot of close-time state — it won't auto-update after fixes.",
        "",
        "_— closure-audit bot (experiment, see #325)_",
    ]
    return "\n".join(lines) + "\n"


# --- gh I/O ---


def _load_notebook(role: str) -> str | None:
    path = REPO_ROOT / "research" / "lab_notebook" / f"{role}.md"
    return path.read_text(encoding="utf-8") if path.exists() else None


def _gh(*args: str) -> str:
    return subprocess.run(
        ["gh", *args], check=True, capture_output=True, text=True
    ).stdout


def fetch_pr(n: int) -> dict:
    return json.loads(_gh(
        "pr", "view", str(n),
        "--json", "mergedAt,closingIssuesReferences,files,number",
    ))


def fetch_issue(n: int) -> dict:
    return json.loads(_gh(
        "issue", "view", str(n),
        "--json", "number,body,labels,comments,closedAt",
    ))


def post_comment(target: str, n: int, body: str) -> None:
    subprocess.run(
        ["gh", target, "comment", str(n), "--body", body],
        check=True,
    )


# --- main orchestration ---


def audit_pr(n: int) -> None:
    pr = fetch_pr(n)
    date = pr["mergedAt"][:10]
    refs = pr.get("closingIssuesReferences") or []
    changed = [f["path"] for f in pr.get("files", [])]

    ac_gaps: list[tuple[int, str]] = []
    pr_gaps: list[tuple[int, str]] = []
    nb_gaps: list[tuple[str, str]] = []

    issues = [fetch_issue(r["number"]) for r in refs]
    for issue in issues:
        body = issue.get("body") or ""
        cmts = [c.get("body", "") for c in issue.get("comments", [])]
        if d := check_ac(body, cmts):
            ac_gaps.append((issue["number"], d))
        if d := check_priority_rationale(body):
            pr_gaps.append((issue["number"], d))

    if not is_exempt(changed):
        labels_per_issue = [
            [lbl["name"] for lbl in i.get("labels", [])] for i in issues
        ]
        role_sets = resolve_roles(labels_per_issue)
        all_roles = {r for rs in role_sets for r in rs}
        notebooks = {r: _load_notebook(r) for r in all_roles}
        # An entry referencing any closing Issue # is as valid as one naming the
        # PR # — entries are often written before the PR exists (#495).
        issue_numbers = [r["number"] for r in refs]
        nb_gaps.extend(
            collect_notebook_gaps(role_sets, date, n, notebooks, also_accept=issue_numbers)
        )

    if ac_gaps or pr_gaps or nb_gaps:
        post_comment("pr", n, format_comment(f"PR #{n}", ac_gaps, pr_gaps, nb_gaps))


def audit_issue(n: int) -> None:
    issue = fetch_issue(n)
    body = issue.get("body") or ""
    cmts = [c.get("body", "") for c in issue.get("comments", [])]
    date = issue["closedAt"][:10]

    ac_gaps: list[tuple[int, str]] = []
    pr_gaps: list[tuple[int, str]] = []
    nb_gaps: list[tuple[str, str]] = []

    if d := check_ac(body, cmts):
        ac_gaps.append((n, d))
    if d := check_priority_rationale(body):
        pr_gaps.append((n, d))

    labels = [lbl["name"] for lbl in issue.get("labels", [])]
    role_sets = resolve_roles([labels])
    all_roles = {r for rs in role_sets for r in rs}
    notebooks = {r: _load_notebook(r) for r in all_roles}
    nb_gaps.extend(collect_notebook_gaps(role_sets, date, n, notebooks))

    if ac_gaps or pr_gaps or nb_gaps:
        post_comment("issue", n, format_comment(f"Issue #{n}", ac_gaps, pr_gaps, nb_gaps))


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--event-type", choices=["pr", "issue"], required=True)
    p.add_argument("--number", type=int, required=True)
    args = p.parse_args()
    if args.event_type == "pr":
        audit_pr(args.number)
    else:
        audit_issue(args.number)
    return 0


if __name__ == "__main__":
    sys.exit(main())
