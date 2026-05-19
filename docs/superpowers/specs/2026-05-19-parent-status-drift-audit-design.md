# Parent-vs-Children Status Drift Audit — Design

**Issue:** [#407](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/407) (feat(pm): parent-vs-children Status drift audit on project board)
**Date:** 2026-05-19
**Author:** PM

## Background

Parent-epic Issues on project board #9 can drift to `Status: In progress` during decomposition (when sub-issues are being added/triaged) and never flip back even after sub-issues settle to Backlog. Per the parent-vs-sub-issues rule (`memory/shared/feedback_parent_sub_issues.md`), parents should carry a **mirrored Status** reflecting their children's collective state — but no mechanism enforces this, so drift accumulates silently.

Today's PM mid-day board sweep (2026-05-19) found 3 epics drifted In progress with all open sub-issues in Backlog ([#24](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/24), [#86](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/86), [#126](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/126)). All 3 were manually corrected to `Ready`, but the drift was only caught by happening to look.

The existing milestone-capacity recheck hook ([PR #397](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/397) closing [#247](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/247)) is structurally identical to what's needed here. This design is a sibling mechanism following the same dispatcher/recheck-script split.

## Architecture & file layout

**New file:** `scripts/pm/recheck_parent_status.py` — recheck logic, callable standalone:
```
scripts/pm/recheck_parent_status.py --issue N      # recheck this issue's parent chain
scripts/pm/recheck_parent_status.py --all           # audit all parent issues on project #9
```
Mirrors `scripts/pm/recheck_milestone.py`'s structure (gh helpers, exit codes 0/1/2).

**Renamed file:** `.claude/hooks/recheck_milestone_dispatch.py` → `.claude/hooks/recheck_dispatch.py` — multiplexer that dispatches both milestone-capacity AND parent-status rechecks.

**Settings change:** `.claude/settings.local.json` hook command path updated from `recheck_milestone_dispatch.py` → `recheck_dispatch.py`. Same PostToolUse Bash matcher.

**Unchanged:** `scripts/pm/recheck_milestone.py` stays as-is — single-responsibility, dispatcher invokes both scripts independently.

## Algorithm — collective state & drift rule

### Status precedence ladder

Lowest → highest "progress":

| Rank | Status |
|------|--------|
| 0 | Backlog |
| 1 | Ready |
| 2 | In progress |
| 3 | Ready for review |
| 4 | In review |
| 5 | Done |

No Status / NO STATUS = treated as `Backlog` (rank 0). Conservative — if a sub-issue hasn't been triaged, it can't argue for parent progress.

### Collective children state

For issue X:
- Iterate native sub-issues via `GET /repos/{owner}/{repo}/issues/{N}/sub_issues`.
- Filter to **open** sub-issues only.
- `collective = max(rank(child.status))` across open children.
- Edge case: zero open sub-issues → `collective = Done` (rank 5), since all decomposed work is closed.

### Drift detection

Three drift classes, surfaced with distinct labels in output; none auto-fixed:

1. **Forward drift** (the case caught today): `rank(parent.status) > collective` AND ≥1 open child.
2. **Backward drift**: `rank(parent.status) < collective`. E.g., child In review but parent still Backlog.
3. **Completion drift**: all children closed (collective = Done) AND parent still open with Status < Done. Decomposition complete, parent unflipped.

### Parent walk

Starting from mutated issue X, fetch `parent_issue` via REST `GET .../issues/{N}` (response includes `parent_issue` if linked). Walk up until `parent_issue` is null. Audit drift at each level.

Caught today's 3-level chain (#204 → #86 → #24) where #86's drift transitively affected #24.

## Triggers — when the dispatcher fires

The renamed `recheck_dispatch.py` keeps all 4 existing milestone trigger shapes and adds these for parent-status:

**New trigger A — project board Status mutation:**
- Pattern: command contains `STATUS_FIELD_ID` constant `PVTSSF_lAHOB17eGc4BSomPzhAHFf8` AND an item ID matching `PVTI_[A-Za-z0-9_-]+`.
- Mirrors the existing `SIZE_FIELD_ID` pattern at `recheck_milestone_dispatch.py:123`. Reuses `lookup_issue_for_item()` to resolve item ID → issue number.
- Calls `recheck_parent_status.py --issue N` which then walks the parent chain.

**New trigger B — `gh issue close N`:**
- Already caught by the existing `PATTERN_CLOSE`. Adds a parent-status recheck alongside the milestone recheck (closing a sub-issue can shift the parent's collective state from "In progress" → "Done").

### Not added (out of scope for v1)

- Sub-issue link/unlink via REST (`POST .../sub_issues`, `DELETE .../sub_issues/{id}`). Rationale: link mutations are rare and already prompt manual board grooming; adding parsing for this surface is YAGNI for v1. Can be added later if drift on link-change becomes a pattern.
- Issue Status changes via `gh issue edit` flags — there isn't a `--status` flag; Status only lives on the project board.

## Output format

### Standalone `--issue N` invocation

```
Parent chain for #204 (walked 2 levels):

#86 (HLA-matched TCR panel) — Status: In progress
  Open sub-issues (3):
    - #204 (Backlog)
    - #205 (Backlog)
    - #206 (Backlog)
  Collective children state: Backlog
  Status: [FORWARD DRIFT] — parent In progress, children all Backlog

#24 (TRUST4 + ProTCR) — Status: In progress
  Open sub-issues (1):
    - #86 (In progress)  ← but #86 itself drifted (see above)
  Collective children state: In progress
  Status: [No change] — but transitive drift via #86
```

### Standalone `--all` invocation

Same per-issue block, repeated for every parent on project #9 (skipping those with no drift).
Header: `Audited N parent issues; M drifted.`

### Hook-fired (additionalContext)

Same body wrapped in `[parent-status recheck — <trigger description>]\n<body>` header to match the milestone hook's emit pattern.

The dispatcher emits **separate** `[milestone recheck ...]` and `[parent-status recheck ...]` blocks so they're independently scannable.

### Exit codes (mirror `recheck_milestone.py`)

- `0` — no drift found
- `1` — error (network / missing data / etc.)
- `2` — drift detected (UPDATE NEEDED equivalent)

## Testing strategy

Three test surfaces:

1. **Algorithm unit tests** (`workflow/tests/test_recheck_parent_status.py` or similar): mock the gh API responses and verify drift detection for each of the 3 classes (forward, backward, completion) plus the "no drift" and "no open children" cases.
2. **End-to-end smoke test:** run `recheck_parent_status.py --all` against the current live project board state. After today's in-session flips of #24/#86/#126 → Ready, the script should find **no drift** (clean baseline).
3. **Hook integration smoke:** issue a benign Status mutation on a known-clean issue and confirm the dispatcher emits the expected `[parent-status recheck — ...]` block with `Status: [No change]`.

## Out of scope (deferred)

- **Auto-mutating parent Status** — audit + surface only; manual confirm preserves judgment on edge cases like "parent intentionally In progress because of an unlinked in-flight PR".
- **Eventual-consistency retry logic** — separate concern in [#406](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/406). The recheck script reads board state via fresh gh calls; if a write hasn't propagated, the recheck will surface stale data. Acceptable for v1.
- **Sub-issue link/unlink triggers** — see "Not added" above.

## Acceptance criteria (from Issue #407)

- [ ] `scripts/pm/recheck_parent_status.py` implemented with `--issue N` and `--all` modes
- [ ] `.claude/hooks/recheck_milestone_dispatch.py` renamed to `recheck_dispatch.py`; extended to dispatch parent-status rechecks on Status mutations + `gh issue close`
- [ ] `.claude/settings.local.json` hook path updated
- [ ] Algorithm unit tests cover the 3 drift classes + no-drift + no-open-children cases
- [ ] `--all` smoke run on current board produces zero drift (post the 2026-05-19 in-session flips)
- [ ] Hook integration smoke confirms dispatcher emits separate `[parent-status recheck ...]` blocks
- [ ] CLAUDE.md / memory updated if the mechanism changes any existing PM workflow

## Related

- Parent-vs-sub-issue rule: `memory/shared/feedback_parent_sub_issues.md` ("mirrored Status")
- Sibling mechanism: [PR #397](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/397) (milestone-capacity recheck hook) closing [#247](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/247)
- Eventual-consistency companion: [#406](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/406)
- Today's standup post surfacing the pattern: `memory/shared/team_standup.md` (2026-05-19 11:56 UTC PM → Sci)
- Original drift-discovery board sweep that motivated this issue: same session 2026-05-19 11:30–12:00 UTC
