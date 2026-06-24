# AGENTS.md slimming via skill extraction

**Epic:** [#859](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/859)
**Pilot:** [#860](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/860)
**Date:** 2026-06-24
**Status:** approved - pilot (#860) implemented, then reclassified (skill -> reference file) before merge after a `skill-creator` review - see "Pilot design" below

## Problem

`AGENTS.md` (the canonical file; `CLAUDE.md` is now a symlink to it) is loaded into context in full on every session, every turn.
It is ~468 lines / ~17.6k tokens across 28 top-level sections.
Most of that weight is *situational how-to*: the regtools BED12 anchor-outer gotcha only matters when editing `alignment.smk`; the four slide-deck-convention sections only matter when authoring a deck; the cloud-GPU infrastructure block only matters during a cloud run.
Paying ~17.6k always-on tokens to keep when-you-are-doing-X content resident is exactly the bloat that progressive disclosure (the design principle behind skills) exists to remove.

This is the sibling of the `MEMORY.md` slimming audit (`project_memory_md_slimming` in PM role memory).
Same technique (an always-on index plus on-demand detail), different scope and curator (project repo + PR here vs personas repo + Memory Manager there).
The two stay separate tracks, cross-referenced, not merged (see Boundary rule below).

## Goal

Move situational, project-scoped content out of always-on `AGENTS.md` into the *right* on-demand home, leaving a one-line **pointer-stub** behind so the topic stays discoverable in resident context.
The home is not always a skill.
Per Anthropic guidance (see "Skill vs reference file vs resident" below), the mechanism is a **mix**:

- **Skills** (`.claude/skills/`) for genuinely *procedural*, task-triggered clusters (how to do X).
- **Reference files** (`docs/`, linked from an `AGENTS.md` stub) for *reference* knowledge read while doing other work (facts, gotchas, rationale).
- **Resident** for short, always-relevant rules.

All homes are committed to the project repo, so they travel to all three clones like the repo's hooks do.
Target end state: resident `AGENTS.md` drops from ~468 lines / ~17.6k tokens toward roughly ~180 lines / ~7k tokens (about a 60% cut), with no content lost.

This corrects an earlier version of this spec that proposed extracting everything into seven skills.
That was a top-down bucket-sort of existing prose, which is the anti-pattern Anthropic warns against; the corrected approach matches content type to home and validates each non-pilot cluster against a real evaluation before building it.

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
| **`.claude/skills/`** (stub points here) | Bulky, situational project how-to that does not need to be resident | Any **Claude** agent, on-demand when the task matches the skill description (a non-Claude agent has no skill mechanism - it sees only the stub and can read the file path manually) | Project repo, via PR |
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

## Skill vs reference file vs resident (web-verified)

Anthropic's guidance on what belongs in a skill, confirmed against the engineering post and the skill-authoring best-practices doc (sources at the end):

- **Skills are for procedural knowledge** - "a way of doing something that activates based on context." The canonical examples are all task-shaped (processing PDFs, analyzing spreadsheets, filling forms, writing commit messages).
- **Pure reference facts are not automatically skills.** Reference content belongs in **reference files loaded on demand**. It rides *inside* a task-triggered skill (a `SKILL.md` overview pointing to reference files the agent reads while doing tasks, e.g. the BigQuery example) or lives as a plain linked file. It is not a skill on its own.
- **Build evaluations first.** Identify a real capability gap by running the task *without* the skill and observing failure, then build the minimal skill to close it. Do not pre-commit skills top-down from an existing document.
- **Conventions:** gerund-form names (`processing-pdfs`); third-person descriptions stating *what* it does and *when* to use it; `SKILL.md` body under ~500 lines with detail split into reference files; only add context Claude does not already have.

**Decision rule for each cluster:**

1. Is it an explicit task an agent performs, with procedural steps or commands? -> **skill** (validated by an eval).
2. Is it reference knowledge needed *while editing code or doing other work* (facts, gotchas, rationale)? -> **reference file** in `docs/`, linked from an `AGENTS.md` stub. A file behind an always-on stub triggers more reliably than a skill, because a file-edit is not a task invocation that would fire a skill description.
3. Is it a short rule relevant to almost everything (e.g. an output-vocabulary convention)? -> **keep resident**.

## The pointer-stub pattern

When a section is extracted, `AGENTS.md` keeps a one-line stub instead of the full text:

> **Authoring an experiment notebook or Quarto slide deck?** Layout, cross-experiment sharing, size bands, and the 3 deck tiers (experiment / eval / research-decision) live in [`docs/research_artifact_conventions.md`](docs/research_artifact_conventions.md).

Rationale: the per-section resident cost drops from tens of lines to ~1, while the topic stays visible in always-on context.
This neutralizes the main failure mode of on-demand skills: a description-trigger that does not fire would otherwise make the convention silently invisible.
A stub can point at a **skill** (procedural cluster), a **reference file** in `docs/` (reference knowledge), or, for a bloated graduated-summary section, back at a **shared-memory file** (the "trim to pointer" verdict below).

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
| Junction origin classification | ~12 | core tumor-specificity logic; foundational and short, not situational how-to |
| MHC Presentation Vocabulary | ~13 | naming rule that shapes all prediction prose/code; always-relevant (the textbook resident case), or a short reference file |

### Extract to skill - procedural, task-triggered (validate each via eval before building)

| Cluster -> skill (gerund name TBD at build) | Source sections | Lines | Status |
|---|---|---|---|
| `running-snakemake` | Snakemake 8 Gotchas + Conda Activation | ~59 | candidate - eval first; first skill-mechanism pilot once its eval passes |
| `running-the-pipeline` | Infrastructure (cloud GPU) + chr22 test dataset | ~37 | candidate - eval first |
| (conda-env setup recipe) | the procedural half of Python environments | part of ~47 | candidate - eval first |

The `authoring-research-decks` cluster started here as the pilot but was **reclassified to a reference file** before merge (see "Pilot design") - it has no procedural steps, so by the decision rule it is reference-shaped. It now appears in the reference-file table below.

### Extract to reference file - reference knowledge read while editing (stub -> `docs/` file)

| Cluster -> reference file | Source sections | Lines |
|---|---|---|
| research-artifact conventions (**pilot, shipped** -> `docs/research_artifact_conventions.md`) | Experiment-notebook layout + Slide decks x3 | ~56 |
| junction-extraction gotchas | regtools arg-order, regtools BED12, STAR strand, UCSC/ENSEMBL, sra-tools, HISAT2 cache | ~62 |
| conda dependency issues | Known Dependency Issues (the "known issues" half) | part of ~47 |
| tcrdock/structure rationale | TCRdock-via-Docker + PDB chain relabelling | ~14 |

### Trim to pointer - graduated summaries bloated into duplicates

| Section | Lines | Trim to |
|---|---|---|
| GitHub Safety Wrappers | 52 | ~6-line index; the hooks/CI self-enforce, full rationale already in shared memory + the hook files |
| Board status governance | 34 | ~10-line core rule + cite `feedback_board_hygiene.md` |
| WIP limits | 11 | ~3-line rule + cite `feedback_board_hygiene.md` |

The trimmed-resident core for each of these must stay **self-sufficient**: an agent (including a non-Claude one that cannot load shared memory) must be able to operate the board from the resident lines alone.
The shared-memory citation is *extra depth*, not load-bearing - the trim relocates rationale and edge cases, never the operative rule itself.

### Leave / delete

| Section | Lines | Action |
|---|---|---|
| Config Migration Notes | 5 | stale historical; delete or fold into a CHANGELOG note |

**Projected result:** ~468 -> ~180 resident lines (~60% cut), every extracted cluster discoverable via its stub.

## Pilot design: research-artifact conventions (reference file)

Lowest-risk cluster (zero safety-criticality: a missed reference means a deck authored in a slightly-wrong folder, easy to course-correct), so it validates the extraction pattern before higher-stakes clusters.

**Reclassified skill -> reference file before merge.**
The pilot was first built as a `.claude/skills/authoring-research-decks/` skill (verbatim move, proven byte-for-byte).
A pre-merge `skill-creator` review then reconsidered it against the decision rule above and found it **reference-shaped, not skill-shaped**: the body is pure layout facts, tier definitions, and rationale - no procedural steps or commands, which is rule-1's test for a skill.
Both `skill-creator` ("pure reference facts are not automatically skills") and this spec's own rule-2 ("reference knowledge read while doing other work -> a `docs/` reference file behind a stub") point the same way, and the cross-tool-asymmetry risk below favors the reference file (any agent can read it; a `.claude/skills/` skill is Claude-only).
So the cluster ships as a reference file:

```
docs/research_artifact_conventions.md
```

with the four current `AGENTS.md` sections moved verbatim under a short title + orienting intro:

1. Experiment notebooks live under `research/experiments/` (layout, cross-experiment data sharing, size bands).
2. Slide decks for experiment Issues.
3. Slide decks for eval Issues.
4. Slide decks for research-decision Issues.

(Render tooling - Quarto, the disabled-PDF / footnotehyper gotcha - travels with section 2 where it currently lives.)

The four sections in `AGENTS.md` are replaced by one pointer-stub (text above) pointing at the reference file.

**Consequence for the rollout:** this pilot now validates the **stub -> reference-file** extraction path (the most common home in the triage).
The **skill-loading mechanism** - whether a `.claude/skills/` skill surfaces and triggers - is no longer exercised here; it is first validated by the earliest *genuinely procedural* candidate (`running-snakemake`) once that clears its eval-first gate.

## Rollout

1. **Pilot PR (this branch, #860):** extract the 4 sections to `docs/research_artifact_conventions.md` (verbatim section bodies), replace them with the stub, land this spec. (Built first as a skill, then reclassified to a reference file before merge - see "Pilot design".)
2. **Per-cluster follow-up Issues** under epic #859, filed once the pilot proves the pattern, each routed by the decision rule:
   - *Skill candidates* (`running-snakemake`, `running-the-pipeline`, conda-env setup recipe): **build an eval first** - run the task without the skill, confirm a real gap and a reliable trigger - then build only if justified.
   - *Reference files* (junction-extraction gotchas, conda dependency issues, tcrdock/structure rationale): move to a `docs/` file + `AGENTS.md` stub. No eval needed; these are read-on-reference, not task-triggered.
3. **Trim pass:** Safety Wrappers / Board governance / WIP to pointers.
4. **Final:** prune any leftover and confirm the resident-token target (~7k).

Cross-reference the `project_memory_md_slimming` audit so the two slimming efforts share the boundary rule above.

## Risks and mitigations

- **Trigger-miss (a skill does not surface when needed):** mitigated by the pointer-stub, which keeps the topic visible in resident context unconditionally.
- **Project-skill loading:** the `.claude/skills/` loading mechanism is no longer exercised by this pilot (it ships a reference file, not a skill). It is first validated by the earliest *genuinely procedural* skill candidate (`running-snakemake`), which must verify the skill appears in the skills list after `/reload-plugins` before the remaining skill clusters depend on it.
- **Content drift during the move:** every PR's diff must account for each moved line (no paraphrasing); the move is verbatim.
- **Docs back-references to `AGENTS.md` sections:** path references resolve via the symlink; section-anchor links (if any) are checked per PR.
- **Cross-tool asymmetry of skills vs reference files:** `AGENTS.md` is cross-tool (Codex/Cursor read it), but a `.claude/skills/` skill is a Claude-only mechanism - a non-Claude agent only sees the stub. A `docs/` reference file behind a stub is readable by any agent. So when a cluster is a close call between "skill" and "reference file", **lean to the reference file** for cross-tool reach; reserve skills for genuinely procedural, task-triggered work where the trigger earns its keep. The directory-level fix for this asymmetry (a vendor-neutral `.agents/` config home) is tracked separately in #861, sequenced after this pilot.

## Verification

- Pilot: `docs/research_artifact_conventions.md` exists with the four section bodies extracted verbatim - byte-identical to the original `AGENTS.md` lines at extraction (sha256-verified) - with one intentional post-extraction delta: a `../` prefix on 3 relative-link targets (`research/slides/README.md`, `research/evals/issue_218_hermes/`, `research/decisions/issue_592/`), required because the file sits one directory deeper than the repo root, where the original root-relative links would 404. Link text is unchanged. The 4 sections are gone from `AGENTS.md` and replaced by the stub; the stub link resolves to the new file; `git diff` accounts for every moved line.
- CI green (no code paths touched).
- Resident-line count of `AGENTS.md` drops by ~56 lines after the pilot, toward the ~180-line end-state target as follow-ups land.

## Sources

The skill-vs-reference-vs-resident guidance was verified against:

- Anthropic, "Equipping agents for the real world with Agent Skills" - https://www.anthropic.com/engineering/equipping-agents-for-the-real-world-with-agent-skills
- Claude Docs, "Skill authoring best practices" - https://platform.claude.com/docs/en/agents-and-tools/agent-skills/best-practices
- MindStudio, "Claude Code Skills Architecture: skill.md process vs reference files" - https://www.mindstudio.ai/blog/claude-code-skills-architecture-skill-md-reference-files
