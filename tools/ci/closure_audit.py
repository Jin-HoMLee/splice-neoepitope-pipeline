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
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
COMMENT_MARKER = "<!-- closure-audit -->"

# Matches only `- [ ]` / `- [x]` style — `*` and `+` bullet markers are not
# recognized. The repo's issue templates use `-` throughout; revisit if other
# styles show up.
_UNTICKED = re.compile(r"^\s*-\s*\[\s\]\s", re.MULTILINE)
_TICKED = re.compile(r"^\s*-\s*\[[xX]\]\s", re.MULTILINE)
_EXEMPT_FILES = {"research/news_log.md", "research/glossary.md"}
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


def check_lab_notebook(text: str, date: str, number: int) -> str | None:
    header = f"## {date}"
    if header not in text:
        return f"no '## {date}' header in notebook"
    start = text.index(header)
    rest = text[start + len(header):]
    nxt = re.search(r"^## ", rest, re.MULTILINE)
    block = rest[:nxt.start()] if nxt else rest
    if not re.search(r"^### ", block, re.MULTILINE):
        return f"'## {date}' has no '### HH:MM UTC — Editor: …' sub-section"
    if f"#{number}" not in block:
        return f"'## {date}' block does not reference '#{number}'"
    return None


def is_exempt(changed_files: list[str]) -> bool:
    return bool(changed_files) and all(
        p in _EXEMPT_FILES or p.startswith(_EXEMPT_PREFIX)
        for p in changed_files
    )


def resolve_roles(labels_per_issue: list[list[str]]) -> list[str]:
    roles: set[str] = set()
    for labels in labels_per_issue:
        role_labels = sorted(
            lbl[len("role:"):]
            for lbl in labels
            if lbl.startswith("role:")
        )
        if role_labels:
            roles.add(role_labels[0])
    return sorted(roles)


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
        for role in resolve_roles(labels_per_issue):
            path = REPO_ROOT / "research" / "lab_notebook" / f"{role}.md"
            if not path.exists():
                nb_gaps.append((role, "lab notebook file missing"))
                continue
            if d := check_lab_notebook(path.read_text(), date, n):
                nb_gaps.append((role, d))

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
    for role in resolve_roles([labels]):
        path = REPO_ROOT / "research" / "lab_notebook" / f"{role}.md"
        if not path.exists():
            nb_gaps.append((role, "lab notebook file missing"))
            continue
        if d := check_lab_notebook(path.read_text(), date, n):
            nb_gaps.append((role, d))

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
