# AGENTS.md slimming via skill extraction

**Epic:** [#859](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/859)
**Pilot:** [#860](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/860)
**Date:** 2026-06-24
**Status:** approved design

## Problem

`AGENTS.md` (the canonical file; `CLAUDE.md` is now a symlink to it) is loaded into context in full on every session, every turn.
It is ~468 lines / ~17.6k tokens across 28 top-level sections.
Most of that weight is *situational how-to*: the regtools BED12 anchor-outer gotcha only matters when editing `alignment.smk`; the four slide-deck-convention sections only matter when authoring a deck; the cloud-GPU infrastructure block only matters during a cloud run.
Paying ~17.6k always-on tokens to keep when-you-are-doing-X content resident is exactly the bloat that progressive disclosure (the design principle behind skills) exists to remove.

This is the sibling of the `MEMORY.md` slimming audit (`project_memory_md_slimming` in PM role memory).
Same technique (an always-on index plus on-demand detail), different scope and curator (project repo + PR here vs personas repo + Memory Manager there).
The two stay separate tracks, cross-referenced, not merged (see Boundary rule below).

## Goal

Extract situational, project-scoped how-to content out of `AGENTS.md` into on-demand **project skills** under `.claude/skills/` (committed, so they travel to all three clones exactly like the repo's hooks do).
For each extracted cluster, leave a one-line **pointer-stub** in `AGENTS.md` so the topic stays discoverable in always-on context even if a skill's own description-trigger does not fire.
Target end state: resident `AGENTS.md` drops from ~468 lines / ~17.6k tokens to roughly ~180 lines / ~7k tokens (about a 60% cut), with no content lost.

## Non-goals

- Rewriting or re-deciding any convention. This is pure relocation; the content moves verbatim.
- Touching `MEMORY.md` / shared memory content (different curator and repo; tracked separately).
- Extracting always-in-effect reflex rules (safety guards, branch naming, merge gate) into on-demand skills. Those stay resident by design (see Boundary rule).

## Boundary rule: where does a piece of knowledge live?

The deciding axes are **scope** (whose knowledge, who must see it) and **maturity** (how cemented / enforced it is), not format.
The pointer-stub and the `MEMORY.md` index look identical in shape, but they index different things.

| Home | Holds | Loaded by | Curator / repo |
|---|---|---|---|
| **`AGENTS.md`** (full, or stub -> skill) | Cemented **project** facts and conventions about the codebase | Every agent + human in the repo, **including non-Claude agents** (Codex/Cursor read AGENTS.md) | Project repo, via PR |
| **`.claude/skills/`** (stub points here) | Bulky, situational project how-to that does not need to be resident | Any agent, on-demand when the task matches the skill description | Project repo, via PR |
| **shared memory** (`shared/MEMORY.md` + files) | Cross-role conventions still evolving, or experiential "why we decided X" | This project's Claude personas only | Personas repo, via Memory Manager |
| **role memory** (`<role>/MEMORY.md` + files) | Persona-specific knowledge (e.g. how PM runs the morning routine) | That one role's Claude sessions only | Personas repo, via Memory Manager |

Decision test, in order:

1. Must a non-Claude agent or a human obey it, and is it cemented? -> `AGENTS.md` (full if it is a short reflex, or a stub -> skill if it is bulky how-to).
2. Is it project knowledge but situational and bulky? -> a skill, with a stub in `AGENTS.md`.
3. Is it a cross-role convention still evolving or mostly "why" context? -> shared memory.
4. Is it specific to one persona's way of working? -> that role's memory.

Why not migrate shared memory into `AGENTS.md` even though it is cross-role:

- **Maturity.** Memory is the staging ground (fast, experiential, "true when written, verify before relying"); `AGENTS.md` is the cemented, current-truth, PR-gated layer. They are two rungs of one ladder (memory -> inline Always-in-effect -> mechanism), not two bins for the same content.
- **Always-on cost.** `AGENTS.md` loads in full every turn; memory detail loads on-demand via recall. Migrating shared memory into `AGENTS.md` would balloon resident context, the opposite of this effort.
- **Curator boundary.** Personas repo + Memory Manager vs project repo + PR, on purpose.

The load-bearing, all-agents-must-obey subset has already correctly graduated into `AGENTS.md` (the no-`@claude` rule, branch naming, the merge gate). What remains in shared memory is the not-yet-graduated or experiential remainder.

## The pointer-stub pattern

When a section is extracted, `AGENTS.md` keeps a one-line stub instead of the full text:

> **Authoring an experiment notebook or Quarto slide deck?** Layout, cross-experiment sharing, size bands, and the 3 deck tiers (experiment / eval / research-decision) live in the `authoring-research-decks` skill.

Rationale: the per-section resident cost drops from tens of lines to ~1, while the topic stays visible in always-on context.
This neutralizes the main failure mode of on-demand skills: a description-trigger that does not fire would otherwise make the convention silently invisible.
A stub can point at a **skill** or, for a bloated graduated-summary section, back at a **shared-memory file** (the "trim to pointer" verdict below).

## Triage: 28 sections, 4 verdicts

Line counts are from the current file.

### Keep resident - all-agent reflexes + orientation

| Section | Lines | Why |
|---|---|---|
| Instructions for Claude | 3 | meta-rule about the file itself |
| Project Overview | 3 | orientation every agent needs |
| Branch naming | 20 | reflex; fires at branch-creation, no task-trigger to catch a miss |
| Merge workflow | 11 | reflex; the `audit_and_merge.sh` gate must stay top-of-mind |
| Three-axis work model | 20 | core board model all roles reason with |

### Extract to skill + stub - situational project how-to

| Cluster -> skill | Source sections | Lines |
|---|---|---|
| `authoring-research-decks` **(pilot)** | Experiment-notebook layout + Slide decks x3 | ~56 |
| `running-snakemake` | Snakemake 8 Gotchas + Conda Activation | ~59 |
| `splice-pipeline-gotchas` | regtools arg-order, regtools BED12, STAR strand, UCSC/ENSEMBL, sra-tools, HISAT2 cache | ~62 |
| `pipeline-conda-envs` | Known Dependency Issues + Python environments | ~47 |
| `running-the-pipeline` | Infrastructure (cloud GPU) + chr22 test dataset | ~37 |
| `pipeline-design-rationale` | Pipeline Design Decisions + Reference/Index layout | ~33 |
| `neoepitope-output-vocabulary` | MHC Presentation Vocabulary | ~13 |

### Trim to pointer - graduated summaries bloated into duplicates

| Section | Lines | Trim to |
|---|---|---|
| GitHub Safety Wrappers | 52 | ~6-line index; the hooks/CI self-enforce, full rationale already in shared memory + the hook files |
| Board status governance | 34 | ~10-line core rule + cite `feedback_board_hygiene.md` |
| WIP limits | 11 | ~3-line rule + cite `feedback_board_hygiene.md` |

### Leave / delete

| Section | Lines | Action |
|---|---|---|
| Config Migration Notes | 5 | stale historical; delete or fold into a CHANGELOG note |

**Projected result:** ~468 -> ~180 resident lines (~60% cut), every extracted cluster discoverable via its stub.

## Pilot design: `authoring-research-decks`

Lowest-risk cluster (zero safety-criticality: a missed trigger means a deck authored in a slightly-wrong folder, easy to course-correct), so it validates the pattern before higher-stakes clusters.

```
.claude/skills/authoring-research-decks/SKILL.md
```

Frontmatter (the `description` is the trigger; skill-creator can later optimize it):

```yaml
---
name: authoring-research-decks
description: >-
  Conventions for authoring research artifacts in this repo - experiment
  notebooks (research/experiments/issue_NNN_<short>/ layout, cross-experiment
  data sharing, size bands) and the three Quarto slide-deck tiers (experiment /
  eval / research-decision). Use when creating or editing a slide deck, an
  experiment notebook, or their outputs/ folder.
---
```

Body: the four current `AGENTS.md` sections moved verbatim:

1. Experiment notebooks live under `research/experiments/` (layout, cross-experiment data sharing, size bands).
2. Slide decks for experiment Issues.
3. Slide decks for eval Issues.
4. Slide decks for research-decision Issues.

(Render tooling - Quarto, the disabled-PDF / footnotehyper gotcha - travels with section 2 where it currently lives.)

The four sections in `AGENTS.md` are replaced by one pointer-stub (text above).

## Rollout

1. **Pilot PR (this branch, #860):** create the `authoring-research-decks` skill, replace the 4 sections with the stub, land this spec. Verify the skill loads after `/reload-plugins`.
2. **Per-cluster follow-up Issues** under epic #859: one PR each for `running-snakemake`, `splice-pipeline-gotchas`, `pipeline-conda-envs`, `running-the-pipeline`, `pipeline-design-rationale`, `neoepitope-output-vocabulary`. Filed once the pilot proves the pattern.
3. **Trim pass:** Safety Wrappers / Board governance / WIP to pointers.
4. **Final:** prune any leftover and confirm the resident-token target (~7k).

Cross-reference the `project_memory_md_slimming` audit so the two slimming efforts share the boundary rule above.

## Risks and mitigations

- **Trigger-miss (a skill does not surface when needed):** mitigated by the pointer-stub, which keeps the topic visible in resident context unconditionally.
- **Project-skill loading:** `.claude/skills/` is committed and active; verify in the pilot that the skill appears in the skills list after `/reload-plugins` (acceptance criterion), so the mechanism is proven before the other six clusters depend on it.
- **Content drift during the move:** every PR's diff must account for each moved line (no paraphrasing); the move is verbatim.
- **Docs back-references to `AGENTS.md` sections:** path references resolve via the symlink; section-anchor links (if any) are checked per PR.

## Verification

- Pilot: skill file exists with the content verbatim; the 4 sections are gone from `AGENTS.md` and replaced by the stub; `git diff` accounts for every moved line; the skill appears after `/reload-plugins`.
- CI green (no code paths touched).
- Resident-line count of `AGENTS.md` drops by ~56 lines after the pilot, toward the ~180-line end-state target as follow-ups land.
