# `scripts/new_branch.sh` — branch-naming helper

**Issue:** [#578](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/578)
**Date:** 2026-06-22
**Status:** approved design

## Problem

Branch names have drifted into 5+ competing patterns. The dominant ~half uses the canonical `type/role/issue-N-slug` — which carries the `role` metadata the board + closure automation consume — but the other half (mostly from raw `gh issue develop`, which auto-slugs the full Issue title) drops it and produces long/mangled names (e.g. the literal `→` and mangled `run_cloud_gpush` in the abandoned [PR #523](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/523) branch).

The convention already exists; the tool that creates branches doesn't follow it. So the fix is a **mechanism, not a memory rule** (mechanism-over-memory ladder).

## Goal

A small `scripts/new_branch.sh` that derives `type` + `role` from an Issue, takes a **short, human-supplied** slug, and creates the canonical `type/role/issue-N-slug` branch by wrapping `gh issue develop` (preserving the Issue↔branch Development-panel link).

## Canonical pattern

```
<type>/<role>/issue-<N>-<slug>
```

- `<type>` — Conventional-Commit type (`feat`/`fix`/`refactor`/`docs`/… or a repo freeform type like `spike`/`design`/`tooling`/`process`/`labnotebook`/`research`/`eval`).
- `<role>` — the single implementer role (`pm` / `scientist` / `developer`), from the Issue's `role:<x>` label.
- `<N>` — the Issue number.
- `<slug>` — a **short, human-chosen** kebab descriptor (not an auto-slug of the title).

This matches the [Conventional Branch](https://conventional-branch.github.io/) best-practice format `<type>/<ticket>-<description>`, plus a `role` segment for the board automation.

## Interface

```bash
# Normal (issue-linked):
scripts/new_branch.sh <issue#> <short-slug> [--type T] [--role R] [--dry-run]
#   new_branch.sh 578 branch-helper              → feat/pm/issue-578-branch-helper
#   new_branch.sh 578 branch-helper --type spike → spike/pm/issue-578-branch-helper

# Issueless fallback (rare — e.g. a lab-notebook branch with no Issue):
scripts/new_branch.sh --no-issue <type> <role> <short-slug> [--dry-run]
#   new_branch.sh --no-issue labnotebook pm memory-slim → labnotebook/pm/memory-slim
```

## Derivation rules

| Segment | Source | Fallback / error |
|---|---|---|
| **type** | leading CC prefix of the Issue title (`feat(scripts):` → `feat`) | `--type` overrides; **no prefix AND no `--type` → error** (never guess) |
| **role** | the Issue's `role:<x>` label(s) | exactly one → use it; **multiple → require `--role`**; none → error |
| **slug** | the user's positional arg, sanitized | **omitted → error**: print the Issue title + a suggested invocation, create nothing |
| **N** | the issue number argument | must be a positive integer |

**Slug sanitization:** lowercase; spaces/underscores → `-`; strip every char not `[a-z0-9-]` (drops unicode like `→`); collapse repeated `-`; trim leading/trailing `-`. If the result is empty → error.

**type override note:** the Issue-title type is a *convenient default* (titles are conventionally prefixed for internal filings), never gospel — the eventual branch type can legitimately differ (a `research(...)` issue branched as `spike`). `--type` is the escape hatch. *Assumption:* internal Issue titles stay CC-prefixed; **revisit** the title-as-default if external/diverse Issue filings grow (then lean on `--type`).

## Creation mechanism — wrap `gh issue develop` (recorded decision)

The issue-linked branch is created via:

```bash
gh issue develop <N> --name <branch> --base main --checkout
```

**Decision: wrap, do not replace with `git checkout -b`.** Rationale: `gh issue develop` is the *only* way to link a branch to the Issue's Development panel — there is no retroactive link. The helper simply feeds it our canonical `--name`. `git fetch origin` runs first so the base is fresh.

The **`--no-issue`** form is the *single* sanctioned `git checkout -b` path (there is no Issue to link):

```bash
git checkout -b <type>/<role>/<slug>
```

## Safety: internal parent guard

The repo's `gh issue develop` parent-guard PreToolUse hook **will not fire** when the call is nested inside this script (the hook inspects the outer `bash new_branch.sh …` command, not the inner `gh issue develop`). So the script **replicates the guard**: the single issue fetch is a `gh api graphql` query returning `title`, the `role:` labels, **and** `subIssuesSummary.total` together; if `total > 0` the target is a parent/epic and the script **refuses** with the same guidance (branch off a leaf, or file a closure sub-issue), exit 1. Branching off a parent would create a PR↔parent `closingIssuesReferences` edge that auto-closes the epic and orphans its sub-issues.

Because the parent check rides the *same* fetch the script already needs for type/role, there is no separate fail-open/closed path: if that fetch fails the script errors out (it can derive nothing), so a parent can never be silently bypassed.

## Error handling (all exit 1, create nothing)

| Condition | Message |
|---|---|
| missing/invalid issue number | usage + "issue number must be a positive integer" |
| missing slug | Issue title + `suggested: new_branch.sh <N> <your-slug>` |
| no title type + no `--type` | "could not derive type from title; pass --type" |
| multiple `role:` labels, no `--role` | "multiple roles: pm, developer — pass --role <one>" |
| no `role:` label, no `--role` | "no role label on issue; pass --role <one>" |
| target is a parent/epic | parent-refuse guidance |
| slug sanitizes to empty | "slug is empty after sanitization" |
| branch already exists | git/gh's native error (not pre-checked) |

`--dry-run` prints the computed branch name and exits 0 **without** any git/gh mutation.

## Structure (for testability)

The script separates **pure logic** (no I/O) from **effects**:

- Pure functions: `parse_type_from_title`, `sanitize_slug`, `pick_role` (given a list of role labels + optional override), `assemble_branch`.
- Effectful: `fetch_issue` (one `gh api graphql` → `{title, roleLabels, subIssuesTotal}`), `create_branch` (`gh issue develop` / `git checkout -b`).

`--dry-run` exercises the full pure path + the `gh issue view` read, stopping before any mutation — which is what the tests drive.

## Testing

`workflow/tests/test_new_branch.py` (runs under the existing `pipeline-pytest` CI check) shells out to the script with a **stub `gh`** placed first on `PATH` that returns canned issue JSON per fixture. Cases:

1. title-type parse → `feat/pm/issue-578-branch-helper`
2. `--type spike` override → `spike/pm/issue-578-branch-helper`
3. multi-role issue, no `--role` → error mentioning both roles; with `--role pm` → succeeds
4. slug sanitization: `"Branch Helper"` → `branch-helper`; `"west1-b → west4-a"` → `west1-b-west4-a`; `"P100__GPU"` → `p100-gpu`
5. missing slug → error printing the title
6. no title type + no `--type` → error
7. parent issue (`subIssuesSummary.total > 0`) → parent-refuse
8. `--no-issue labnotebook pm memory-slim` → `labnotebook/pm/memory-slim`

All assertions run against `--dry-run` output (no branch is actually created in tests). `shellcheck` lints the script.

## Documentation

A short **"Branch naming"** subsection in `CLAUDE.md`: the canonical pattern, `new_branch.sh` as the way to produce it, the `--no-issue` fallback, and the wrap-not-replace decision (interim until the helper is habitual, per the AC).

## Out of scope (YAGNI)

- type allowlist/validation (repo uses freeform types) — accept any sanitized lowercase token.
- commit-message, PR-title, Issue-title conventions (already strong; explicitly out of scope per #578).
- any board mutation (Status moves stay manual / separate).

## Acceptance criteria (from #578)

- [x] derives `type` + `role` and creates `type/role/issue-N-slug` *(design)*
- [x] short slug is the user's, not an auto-slug *(design)*
- [x] handles issueless branches with a documented fallback (`--no-issue`) *(design)*
- [x] canonical pattern documented in `CLAUDE.md` *(design — implement in PR)*
- [x] decide + record: wrap `gh issue develop` (not replace) *(design)*
