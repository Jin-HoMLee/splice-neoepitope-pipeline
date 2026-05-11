---
title: Post-merge / post-close closure critic
issue: https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/325
status: draft
created: 2026-05-11
author: Developer
---

# Post-merge / post-close closure critic — design

## Background

The morning closure audit (PM, manual) verifies that closed issues and merged PRs satisfy the closure ritual:

1. All acceptance-criteria checkboxes ticked (or comment-deferred).
2. A per-role lab notebook entry exists under the close/merge date and references the PR/issue number.
3. The issue body contains a `**Priority rationale:**` line (case-insensitive substring).

Manual auditing produced false positives on 2026-05-10 (caught by `pm/feedback_closure_audit_method.md`). A standing post-merge / post-close critic automates the same checks, shortening the audit loop and removing human-grep error. Reference: [Sakana Fugu](https://sakana.ai/fugu-beta/) critic-as-default lesson; morning discussion 2026-05-11.

## Goals

- Run the 3 closure-ritual checks automatically on PR-merge and issue-close.
- Post a single non-blocking comment on the closing PR (or the closed issue) listing any gaps. Post-only (no edit-in-place); the comment is a snapshot of close-time state.
- Be silent on the clean path (no comment if all 3 checks pass).
- Keep the script small (~180 LOC) and unit-tested at the parsing edges. Easy to delete if the experiment doesn't pan out — see kill criterion on PR #332.

## Non-goals (v2 candidates)

- Per-checkbox deferral matching. v1 treats any `❎`-marker deferral comment as covering all unticked boxes on that issue.
- Semantic title-match for the lab notebook. v1 requires a literal `#<N>` reference under the date section.
- PR Test plan checkbox audit (closure ritual companion rule). Same shape as v1 AC check, deferred.
- Auto-fix actions (ticking boxes, generating notebook stubs). Critic only reports gaps.

## Architecture

```
┌─────────────────────────────┐         ┌─────────────────────────────┐
│ Event: pull_request closed  │         │ Event: issues closed        │
│ + merged == true            │         │                             │
└──────────────┬──────────────┘         └─────────────┬───────────────┘
               │                                       │
               ▼                                       ▼
   ┌─────────────────────────────────────────────────────────────┐
   │ .github/workflows/closure-audit.yml                          │
   │                                                              │
   │ - checkout (shallow)                                         │
   │ - setup-python                                               │
   │ - run: python tools/ci/closure_audit.py                      │
   │     --event-type {pr|issue}                                  │
   │     --number {N}                                             │
   └──────────────────────────────┬──────────────────────────────┘
                                  ▼
   ┌──────────────────────────────────────────────────────────────┐
   │ tools/ci/closure_audit.py (single file — no lib split)       │
   │                                                              │
   │ 1. Fetch payload via `gh` (issue/PR body, labels,            │
   │    closingIssuesReferences, changed files for PR, comments)  │
   │ 2. Exemption-path filter on PR (skip notebook check iff all  │
   │    touched paths match exemption list)                       │
   │ 3. Run the 3 checks → collect gaps                           │
   │ 4. If any gap: post comment via `gh {pr|issue} comment`.     │
   │    Silent if all 3 checks pass.                              │
   │ 5. Exit 0 (non-blocking)                                     │
   └──────────────────────────────────────────────────────────────┘
```

### Auth

Uses the default `GITHUB_TOKEN`. Required permissions on the workflow:
- `issues: write` (post comments on issues)
- `pull-requests: write` (post comments on PRs)
- `contents: read` (read lab notebook files at the merge SHA)

No PAT or OAuth required.

### Idempotency

v1 is **post-only** — comments are written but never edited or deleted. The audit fires exactly once per merge / close event in normal flow, so re-fire / drift handling is deferred. Posted comments are prefixed with the marker `<!-- closure-audit -->` so a future v2 could grep for them and update in place.

Silent on the clean path (no comment if all 3 checks pass). The all-clear case is a no-op — there is no `✅ All clear` notification.

## The 3 checks

### Check 1: AC checkboxes ticked

**Targets:** for PR-merge, each issue in `closingIssuesReferences`; for solo issue-close, the issue itself.

**Logic:**

```python
unticked = len(re.findall(r'(?mi)^\s*-\s*\[\s\]\s', body))
ticked   = len(re.findall(r'(?mi)^\s*-\s*\[x\]\s', body))
deferral_marker_seen = any('❎' in c and 'deferred' in c.lower() for c in comments)
gap = (unticked > 0) and not deferral_marker_seen
```

**v1 simplification:** any deferral comment ⇒ all unticked boxes considered deferred. False-negative risk acknowledged.

### Check 2: Lab notebook entry exists

**Gating:** runs only if (a) the closing PR is not exempted by file-path filter, and (b) role is determinable from a `role:*` label on the target issue.

**Role resolution for PR-merge with multiple linked issues:** if multiple linked issues carry different `role:*` labels, the critic runs Check 2 once per *distinct* role and reports per-role gaps. (Same role appearing on multiple linked issues: dedupe.) If any single linked issue has multiple `role:*` labels (rare), use the first in alphabetical order and record the choice in the comment for transparency.

**Exemption paths** (PR-merge event only; if PR's `changed_files` are all members of this set, skip the notebook check):

- `research/news_log.md`
- `research/glossary.md`
- `research/lab_notebook/*.md` (the notebook files themselves)

**Date:** UTC date of merge (`pull_request.merged_at`) or close (`issue.closed_at`), formatted `YYYY-MM-DD`.

**Logic:**

```python
path = f"research/lab_notebook/{role}.md"
text = read(path)
gap = (
    f"## {date}" not in text
    OR no "### " line in the block between "## {date}" and the next "## " header
    OR f"#{number}" not in the same block
)
```

`{number}` is the merged PR number for PR-merge events, the closed issue number for solo issue-close events.

### Check 3: Priority rationale present

**Targets:** same as Check 1.

**Logic:**

```python
gap = 'priority rationale' not in body.lower()
```

Matches PM's existing manual loose-substring method (`pm/feedback_closure_audit_method.md`).

## Comment format

Single Markdown comment, marker-prefixed:

```markdown
<!-- closure-audit -->
## Closure audit — gaps found

PR #325 merged with the following items still open:

- ⚠️ **AC checkboxes** on linked Issue #999 — 3/5 unticked, no deferral comment found
- ⚠️ **Lab notebook entry** missing for `developer` on 2026-05-11 (no `#325` reference under `## 2026-05-11` in `research/lab_notebook/developer.md`)
- ✅ **Priority rationale** present

Per [closure ritual](shared/feedback_closure_ritual.md): tick the boxes, add a comment-deferral, or add the missing entry/rationale. This comment is a snapshot of close-time state — it won't auto-update after fixes.

_— closure-audit bot (experiment, see #325)_
```

Clean state: no comment posted (silent no-op).

## File layout

```
.github/workflows/closure-audit.yml      # event triggers, dispatches script
tools/ci/closure_audit.py                # main: pure checks + gh I/O in one file
tools/ci/test_closure_audit.py           # pytest unit tests (parsing edge cases)
tools/ci/conftest.py                     # sys.path bootstrap for the non-package dir
```

v1 keeps everything in one ~180-LOC file (no lib/main split). The pure-function surface is still trivially testable since tests import `closure_audit` directly and call the check / formatter functions without going through `_gh()`.

## Test strategy

Lean pytest. Tests target the gnarly parsing and edge cases — block-slicing across multi-date notebooks, deferral-comment handling, exemption-path filtering, role resolution dedupe. Trivial one-liner checks (e.g. `"priority rationale" in body.lower()`) skipped intentionally.

CI integration: separate `ci-tools-pytest` job in `.github/workflows/tests.yml`, isolated from `pipeline-pytest` so a tooling-test failure doesn't block real pipeline PRs.

## Edge cases

- **PR merge with no `closingIssuesReferences`:** Check 1 + Check 3 (issue-body checks) skipped — no target issue to inspect. Check 2 still runs if the merged PR carries a `role:*` label directly; otherwise skipped. If all 3 checks end up skipped, no comment is posted (silent no-op).
- **Issue close where the issue was never opened-and-actioned by a role (e.g. closed as duplicate / wontfix):** the critic still runs the checks. Roles using `wontfix`/`duplicate` close paths should add a `❎ deferred — closed as <reason>` comment to suppress the AC gap. Spec change deferred to v2 if the noise level warrants.
- **Lab notebook file missing entirely for the role** (`research/lab_notebook/<role>.md` not found): treated as Check 2 gap with message "lab notebook file missing for role `<role>`". Should never happen in practice (all 3 role files exist as of 2026-05-11).

## Open questions

None at design time.

## References

- Issue [#325](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/325)
- `shared/feedback_closure_ritual.md` — the rule the critic enforces
- `pm/feedback_closure_audit_method.md` — manual audit prior art (substring grep, fresh GitHub state)
- `research/news_log.md` 2026-05-11 09:11 UTC entry — Sakana Fugu critic-as-default discussion
