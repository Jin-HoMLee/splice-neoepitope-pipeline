"""Post-merge / post-close closure-audit critic.

Spec: docs/superpowers/specs/2026-05-11-post-merge-critic-design.md

Trying Sakana Fugu "critic-as-default" on PR-merge / issue-close. Posts a
single marker-tagged gap comment, or stays silent if all 3 checks pass.
No edit-in-place — comment is a post-only snapshot of close-time state.

The three checks: AC checkboxes, Priority-rationale line, lab-notebook entry.
The lab-notebook check is skipped when the PR is file-exempt (only glossary /
lab-notebook files touched) OR the PR body carries the routine-ship opt-out
marker `<!-- skip-lab-notebook: routine -->` (#555). The marker honors the
routine single-PR-closes-single-Issue skip that #483 declared optional, so the
bot stops coercing entries the lab-notebook rule says are unnecessary. The AC
and Priority-rationale checks are unaffected by the marker.

On the issue path, the lab-notebook check is also skipped when the Issue closed
as `not_planned` (descoped / superseded): such a close ships no work and routes
through a closing comment, not a notebook entry (closure ritual, #743). The AC
and Priority-rationale checks still run — superseded ACs are annotated to their
disposition (e.g. `- [superseded]`) rather than left as `- [ ]`.

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
from typing import NamedTuple

REPO_ROOT = Path(__file__).resolve().parents[2]
COMMENT_MARKER = "<!-- closure-audit -->"

# Matches only `- [ ]` / `- [x]` style — `*` and `+` bullet markers are not
# recognized. The repo's issue templates use `-` throughout; revisit if other
# styles show up.
_UNTICKED = re.compile(r"^\s*-\s*\[\s\]\s", re.MULTILINE)
_TICKED = re.compile(r"^\s*-\s*\[[xX]\]\s", re.MULTILINE)
_EXEMPT_FILES = {"research/glossary.md"}
_EXEMPT_PREFIX = ("research/lab_notebook/",)
# Roles that keep no project-repo lab notebook and so are skipped by the
# notebook check (see #748). The Memory Manager works in the personas (memory)
# repo, not a project-repo clone — its record is the personas-repo git log
# (`shared/feedback_lab_notebook.md` MM-exemption clause); there is no
# `research/lab_notebook/memory_manager.md`. Stripped per-role, so a mixed
# (e.g. developer + memory_manager) Issue still requires the other role's entry.
_NOTEBOOK_EXEMPT_ROLES = {"memory_manager"}
# Routine-ship lab-notebook opt-out (see #555): a PR author may skip the
# notebook check by placing this marker in the PR body. Honors the routine
# single-PR-closes-single-Issue skip that #483 declared optional — the bot must
# not coerce an entry the lab-notebook rule says is unnecessary. The value after
# the colon (e.g. `routine`) is free-text rationale with no `>` (the regex stops
# at the first `>`); presence of the marker is what matters, not the value.
_SKIP_LAB_NOTEBOOK = re.compile(r"<!--\s*skip-lab-notebook\b[^>]*-->", re.IGNORECASE)


# --- pure checks ---


def check_ac(body: str, comments: list[str]) -> str | None:
    """Gap when the `## Acceptance criteria` section has unticked boxes.

    Scoped to the AC section via scan_ac_boxes (Issue #726), mirroring the
    pre-merge gate's `unticked_under "Acceptance criteria"` so the two gates
    agree on what an AC checkbox is. Boxes under other headings (a non-AC
    `## Flags to evaluate` / `## Tasks` checklist) are NOT AC gaps — those are
    surfaced advisorily by the merge-time stray-box lint (check_stray_ac_boxes,
    Issue #730), not flagged here. Previously this scanned the whole body and
    false-flagged any non-AC checklist (Issue #411 via PR #720).
    """
    scan = scan_ac_boxes(body)
    if scan.ac_unticked == 0:
        return None
    if any("❎" in c and "deferred" in c.lower() for c in comments):
        return None
    return f"{scan.ac_unticked}/{scan.ac_total} unticked, no deferral comment found"


# A markdown `## ` heading (exactly two hashes — `### ` and deeper are sub-
# headings, not section boundaries, mirroring audit_and_merge.sh's `^## ` awk).
_H2 = re.compile(r"^##[ \t]+(.+?)\s*$")
# An `## Acceptance criteria` heading (case-insensitive; trailing content
# tolerated), mirroring the bash gate's `^## Acceptance criteria([[:space:]]|$)`.
_AC_HEADING = re.compile(r"acceptance criteria(\s|$)", re.IGNORECASE)


class AcBoxScan(NamedTuple):
    """Checkbox census of an Issue body, partitioned by AC-section membership.

    Shared scan used by the #730 stray-box lint (warns when gating-looking boxes
    live outside an Acceptance-criteria section) and reused by the #726 scoping
    of check_ac to the AC section. `stray_*` covers `- [ ]`/`- [x]` boxes NOT
    under an `## Acceptance criteria` heading; `ac_*` covers those that are.
    """

    has_ac_section: bool
    ac_unticked: int
    ac_total: int
    stray_unticked: int
    stray_headings: list[str]


def scan_ac_boxes(body: str) -> AcBoxScan:
    """Partition a body's `- [ ]`/`- [x]` boxes into AC-section vs stray.

    Walks lines tracking the current `## ` heading; a checkbox is an AC box when
    its enclosing heading matches `Acceptance criteria`, else stray. Boxes before
    any heading are stray under the label `(top of body)`. `stray_headings` lists
    the distinct headings carrying ≥1 *unticked* stray box, in first-seen order.
    """
    in_ac = False
    has_ac = False
    cur_heading = "(top of body)"
    ac_unticked = ac_total = stray_unticked = 0
    stray_headings: list[str] = []

    for line in body.splitlines():
        m = _H2.match(line)
        if m:
            cur_heading = m.group(1)
            in_ac = bool(_AC_HEADING.match(cur_heading))
            has_ac = has_ac or in_ac
            continue
        unticked = bool(_UNTICKED.match(line))
        ticked = bool(_TICKED.match(line))
        if not (unticked or ticked):
            continue
        if in_ac:
            ac_total += 1
            ac_unticked += int(unticked)
        elif unticked:
            stray_unticked += 1
            if cur_heading not in stray_headings:
                stray_headings.append(cur_heading)

    return AcBoxScan(has_ac, ac_unticked, ac_total, stray_unticked, stray_headings)


def check_stray_ac_boxes(body: str) -> str | None:
    """Non-blocking lint (#730): warn when unticked boxes live outside an AC
    section that doesn't exist.

    Returns None when an `## Acceptance criteria` section is present (the blocking
    AC gate owns those boxes) or when there are no stray unticked boxes. Otherwise
    returns a warning naming the count + non-AC heading(s), prompting the author
    to move any deliverable-gating boxes under a canonical heading.
    """
    scan = scan_ac_boxes(body)
    if scan.has_ac_section or scan.stray_unticked == 0:
        return None
    headings = ", ".join(f"'{h}'" for h in scan.stray_headings)
    return (
        f"{scan.stray_unticked} unticked checkbox(es) under {headings} but no "
        f"'## Acceptance criteria' section. If any are deliverable-gating, move "
        f"them under a canonical '## Acceptance criteria' heading so the closure "
        f"gate enforces them (non-AC checklists are intentionally exempt)."
    )


# --- #665: cross-repo closing forward-links ---

# GitHub's closing keywords (close/fix/resolve + their inflections), as whole
# words, case-insensitive. A line carrying one is read as a closing-intent line.
_CLOSING_KEYWORD = re.compile(r"\b(close[sd]?|fix(e[sd])?|resolve[sd]?)\b", re.IGNORECASE)
# A cross-repo Issue reference, `owner/repo#N`. The owner/repo segments mirror
# GitHub's allowed chars (alnum, `-`, `_`, `.`); `#` then the issue number.
_CROSS_REPO_SHORT = re.compile(
    r"\b([A-Za-z0-9][\w.-]*)/([A-Za-z0-9][\w.-]*)#(\d+)\b"
)
# The full-URL form, `https://github.com/owner/repo/issues/N`.
_CROSS_REPO_URL = re.compile(
    r"https?://github\.com/([\w.-]+)/([\w.-]+)/issues/(\d+)"
)


def parse_cross_repo_ac_targets(
    pr_body: str | None, this_repo: str | None
) -> list[tuple[str, int]]:
    """Cross-repo Issues a PR intends to close, as `(owner/repo, number)`.

    GitHub's `closingIssuesReferences` spans only the same repo, so a cross-repo
    close — e.g. a personas-repo PR closing a project-repo Issue — is expressed
    textually in the PR body. This finds, per line, a closing keyword
    (`close`/`fix`/`resolve` + inflections) co-occurring with a cross-repo Issue
    reference in `owner/repo#N` or full-URL form, and returns each distinct
    target. References to `this_repo` (the PR's own repo) are excluded — native
    `closingIssuesReferences` already covers same-repo, and the same-repo AC gate
    in audit_and_merge.sh audits those. A bare `#N` (no `owner/repo`) is
    inherently same-repo and never matches.

    The keyword↔reference association is line-scoped: it accepts both the GitHub
    form (`Closes owner/repo#N`) and the project's link+keyword form
    (`[Issue #N](url) (closes)`), at the cost of matching a non-closing line that
    merely co-mentions a closing word and a cross-repo ref (conservative — an
    extra target whose ACs are ticked passes anyway). To reference a cross-repo
    Issue without gating on it, keep the closing keyword off that line. Deduped,
    first-seen order. Returns [] for an empty body.
    """
    targets: list[tuple[str, int]] = []
    seen: set[tuple[str, int]] = set()
    this = this_repo.lower() if this_repo else None
    for line in (pr_body or "").splitlines():
        if not _CLOSING_KEYWORD.search(line):
            continue
        for rx in (_CROSS_REPO_SHORT, _CROSS_REPO_URL):
            for m in rx.finditer(line):
                slug = f"{m.group(1)}/{m.group(2)}"
                num = int(m.group(3))
                if this and slug.lower() == this:
                    continue  # same-repo → native references cover it
                key = (slug.lower(), num)
                if key in seen:
                    continue
                seen.add(key)
                targets.append((slug, num))
    return targets


def check_priority_rationale(body: str) -> str | None:
    if "priority rationale" in body.lower():
        return None
    return "no 'Priority rationale' line in issue body"


def check_lab_notebook(
    text: str, date: str, number: int, also_accept: Collection[int] = ()
) -> str | None:
    """Gap unless the `## date` block references `#number` OR any `#also_accept`.

    `also_accept` carries the PR's closing-Issue numbers (see #495): entries
    written before the PR exists reference the Issue, not the PR number — both
    are semantically the same unit of work.
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


def skip_lab_notebook(pr_body: str | None) -> bool:
    """True if the PR body carries the routine-ship lab-notebook opt-out marker.

    Marker form: `<!-- skip-lab-notebook: routine -->` (see #555). Matched
    case-insensitively and tolerant of comment whitespace. Skips only the
    lab-notebook check — AC-checkbox and priority-rationale checks still run.
    """
    return bool(_SKIP_LAB_NOTEBOOK.search(pr_body or ""))


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
        # Drop notebook-exempt roles (e.g. memory_manager, #748) before the
        # check; a pure-exempt Issue collapses to an empty set and is skipped,
        # while a mixed Issue still enforces its non-exempt roles' entries.
        roles = roles - _NOTEBOOK_EXEMPT_ROLES
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


def _gh(*args: str, repo: str | None = None) -> str:
    cmd = ["gh", *args]
    if repo:
        cmd += ["--repo", repo]
    return subprocess.run(
        cmd, check=True, capture_output=True, text=True
    ).stdout


def fetch_pr(n: int, repo: str | None = None) -> dict:
    return json.loads(_gh(
        "pr", "view", str(n),
        "--json", "mergedAt,closingIssuesReferences,files,number,body",
        repo=repo,
    ))


def fetch_issue(n: int, repo: str | None = None) -> dict:
    return json.loads(_gh(
        "issue", "view", str(n),
        "--json", "number,body,labels,comments,closedAt,stateReason",
        repo=repo,
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

    if not is_exempt(changed) and not skip_lab_notebook(pr.get("body")):
        labels_per_issue = [
            [lbl["name"] for lbl in i.get("labels", [])] for i in issues
        ]
        role_sets = resolve_roles(labels_per_issue)
        all_roles = {r for rs in role_sets for r in rs}
        notebooks = {r: _load_notebook(r) for r in all_roles}
        # An entry referencing any closing Issue # is as valid as one naming the
        # PR # — entries are often written before the PR exists (see #495).
        issue_numbers = [r["number"] for r in refs]
        nb_gaps.extend(
            collect_notebook_gaps(role_sets, date, n, notebooks, also_accept=issue_numbers)
        )

    if ac_gaps or pr_gaps or nb_gaps:
        post_comment("pr", n, format_comment(f"PR #{n}", ac_gaps, pr_gaps, nb_gaps))


def audit_pr_pre_merge(
    n: int, today: str, repo: str | None = None
) -> list[tuple[str, str]]:
    """Lab-notebook gaps for PR #n about to be merged on `today` (YYYY-MM-DD).

    Pre-merge sibling of audit_pr's notebook check (Issue #409), single-sourcing
    the same is_exempt / skip_lab_notebook / resolve_roles / collect_notebook_gaps
    logic. Two differences from audit_pr: the PR is not merged yet so `mergedAt`
    is null — the caller passes today's UTC date — and the gaps are *returned*
    for a blocking pre-merge gate rather than posted as a comment. Only the
    notebook check is mirrored; AC-checkbox + priority-rationale are already
    enforced by audit_and_merge.sh's bash checks.

    `repo` (Issue #607) is forwarded to the gh I/O layer so the gate composes
    with the `REPO` override that its sibling gates (stray_closers /
    bot_review_offer) honor in audit_and_merge.sh. When None (the post-hoc bot
    path, audit_pr / audit_issue), gh resolves the repo from git context — the
    production closure-audit bot is unaffected.

    Notebook text is read from the working tree (_load_notebook), so the gate
    must run from the PR branch where the entry was written — the documented
    closure-ritual flow (write the entry, then merge). Reading the entry from the
    PR head ref instead is a v2 candidate. Returns [] when clear, exempt, or
    skip-marked.
    """
    pr = fetch_pr(n, repo=repo)
    refs = pr.get("closingIssuesReferences") or []
    changed = [f["path"] for f in pr.get("files", [])]
    if is_exempt(changed) or skip_lab_notebook(pr.get("body")):
        return []
    issues = [fetch_issue(r["number"], repo=repo) for r in refs]
    labels_per_issue = [
        [lbl["name"] for lbl in i.get("labels", [])] for i in issues
    ]
    role_sets = resolve_roles(labels_per_issue)
    all_roles = {r for rs in role_sets for r in rs}
    notebooks = {r: _load_notebook(r) for r in all_roles}
    issue_numbers = [r["number"] for r in refs]
    return collect_notebook_gaps(
        role_sets, today, n, notebooks, also_accept=issue_numbers
    )


def collect_stray_ac_warnings(
    n: int, repo: str | None = None
) -> list[tuple[int, str]]:
    """Stray-AC-box warnings for PR #n's linked Issues (Issue #730).

    For each Issue in closingIssuesReferences, run the non-blocking
    check_stray_ac_boxes lint and collect `(issue_number, warning)` for any that
    have unticked boxes outside a (missing) `## Acceptance criteria` section.
    Single-sources the lint for both the merge-time gate (ac_section_lint.py) and
    its tests, mirroring audit_pr_pre_merge's shape. `repo` is forwarded to the gh
    I/O layer so the gate composes with audit_and_merge.sh's REPO override (#607).
    """
    pr = fetch_pr(n, repo=repo)
    refs = pr.get("closingIssuesReferences") or []
    warnings: list[tuple[int, str]] = []
    for r in refs:
        issue = fetch_issue(r["number"], repo=repo)
        if msg := check_stray_ac_boxes(issue.get("body") or ""):
            warnings.append((issue["number"], msg))
    return warnings


def collect_cross_repo_ac_gaps(
    n: int, repo: str | None = None
) -> list[tuple[str, str]]:
    """Blocking AC gaps for PR #n's cross-repo closing targets (Issue #665).

    The same-repo AC gate in audit_and_merge.sh keys off native
    `closingIssuesReferences`, which GitHub only populates within one repo — so a
    cross-repo close (a personas-repo PR closing a project-repo Issue) gets no AC
    gate at all. This parses the PR body for cross-repo closing forward-links
    (parse_cross_repo_ac_targets) and audits each target Issue's Acceptance
    criteria via the same check_ac the same-repo gate's logic mirrors, returning
    `(owner/repo#N, gap)` for any with unticked AC boxes.

    `repo` is the PR's own repo (the REPO override, #607); its own forward-links
    are excluded as same-repo. Each cross-repo target is fetched from *its* repo
    via fetch_issue(..., repo=owner/repo). Returns [] when there are no cross-repo
    closing targets or all their ACs are ticked.
    """
    pr = fetch_pr(n, repo=repo)
    targets = parse_cross_repo_ac_targets(pr.get("body") or "", repo)
    gaps: list[tuple[str, str]] = []
    for slug, num in targets:
        issue = fetch_issue(num, repo=slug)
        body = issue.get("body") or ""
        cmts = [c.get("body", "") for c in issue.get("comments", [])]
        if d := check_ac(body, cmts):
            gaps.append((f"{slug}#{num}", d))
    return gaps


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

    # A not_planned (descoped / superseded) close ships no work, so it routes
    # through a closing comment, not a lab-notebook entry (closure ritual, #743).
    # Skip ONLY the notebook check — AC + priority-rationale still run (the AC
    # boxes are annotated to their disposition, e.g. `- [superseded]`).
    if (issue.get("stateReason") or "").upper() != "NOT_PLANNED":
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
