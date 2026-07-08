# Coordinator autonomy envelope - design

**Date:** 2026-07-08
**Author:** PM
**Status:** Draft (for review)
**Scope:** Shared (all roles), but PM-owned; it is the keystone gating [epic #1072](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1072) rungs 3-5.
**Issue:** [#1074](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1074) (Sub 2 of the epic).

## 1. Problem

Before any decide-and-dispatch automation (auto-routing pings, triggered Ready-pull, enveloped replenishment - rungs 3-5 of the epic), we need a formalized **autonomy envelope**: the pre-authorized policy stating what a coordinating agent may commit/dispatch **without** per-action human approval, and where it must stop.

Our governance already encodes this gate implicitly, in two rules:

- `autonomy-merge-gate-cadence` - run the reversible chain autonomously, stop at `gh pr merge`.
- `decision-ratification != action-authorization` - reversible work is autonomous; each irreversible outward act keeps its own gate.

The envelope's job is to make that policy **explicit and complete** - one named artifact every dispatch loop can obey - not to invent new policy. Without it, each dispatch rung would re-litigate "may I do this unattended?" ad hoc, which is exactly the action-time judgment our rules historically slip on.

## 2. Grounding and honest caveat

- **Inspired by, not copied from, the paper.** The SPM vision paper's **working modes** scale autonomy by *task complexity x risk* ([arXiv 2601.16392](https://arxiv.org/html/2601.16392v1)). Our per-act gate is **reversibility / blast-radius** instead. We keep the inspiration and state the difference plainly rather than importing a gate we do not actually use.
- **Altitude separation (load-bearing).** The paper's four working modes (assisted / supervised / delegated / autonomous) are **not** an act-level classifier - they describe a *posture over a class of work* and how it matures. That is the **epic's** altitude (the Supervised -> Delegated trajectory these tiers let us climb over time), not this spike's act-level gate. The Issue title's phrase "(working modes)" imported that posture vocabulary into an act-level spike; this spec separates the two explicitly (see section 6).

## 3. The 3-tier reversibility ladder

A **linear** ladder keyed on reversibility / blast-radius, chosen over a 2-axis (reversibility x novelty) matrix because **action-time classification speed** is where our rules slip - a linear read beats a grid lookup at the trigger moment.

| Tier | Name | Meaning | Where the stop is |
|------|------|---------|-------------------|
| **A** | Autonomous | Reversible, in-envelope. Do it. | No stop (but see the surface-if-novel clause, section 5). |
| **B** | Autonomous up to the gate | Irreversible-routine. Run the whole chain autonomously; the human performs the single irreversible act. | At the final act (e.g. the merge command). |
| **C** | Always stop | Novel / precedent / governance / cross-role-conflict / irreversible-outward beyond merge. Propose first. | Before committing effort down the path. |

The tiers differ by **where the stop sits**: Tier B does all the reversible preparation and hands over at the one irreversible act; Tier C stops at the *start*, before effort is spent down a novel path.

## 4. Act -> tier mapping

Each atomic coordination act has a **default tier**. The dispatch rungs (3-5) compose these atomic acts; they do not introduce new ones.

| Coordination act | Default tier | Rule |
|---|---|---|
| **Commit** (Backlog -> Ready) | **A** | Auto up to the per-role floor, on the active arc, DoR met. Reversible - moves back to Backlog cheaply. |
| **Dispatch** (route a `To:<role>` ping, trigger a Ready-pull) | **B now -> A on maturation** | The epic's target unlock. Reversible + low blast-radius; the dispatched work keeps its *own* downstream gates, so authorizing the dispatch does not authorize its merge. **Starts at Tier B** (propose the dispatch, human confirms) and **promotes to Tier A** once (1) the review-axis rung (epic rung 1) has shipped and (2) `hook_fires.jsonl` evidence shows dispatch behaving within a small blast radius. This is the one act whose tier is *not* fixed - it climbs as it earns trust (see section 7). |
| **Merge** (`gh pr merge`) | **B** | Irreversible outward act. Run branch -> build -> test -> review -> lab-notebook autonomously; stop at the merge command for the human's one touch. This *is* `autonomy-merge-gate-cadence`. **Flips to C when the PR itself trips a Tier-C trigger** (novel scope / precedent / governance / cross-role conflict) - then the merge is no longer a bounded ratification of verified work but carries novel judgment CI and review cannot check. |
| **Escalate** (flag / ask / surface a conflict) | **A** | Always permitted - an agent never needs approval to *raise* a concern. Escalation is also the mechanism a Tier-C trigger invokes. |

## 5. Two refinements (beyond the original handoff)

### 5.1 Tier C is an override condition, not a fourth act

Rather than enumerate "Tier-C acts" (which would duplicate the table), Tier C is a set of **triggers that bump *any* act to always-stop**:

- **Novel scope** - work with no established pattern to follow.
- **Precedent / governance-rule change** - anything that establishes or edits a rule, convention, or policy (including memory and CLAUDE.md edits - so ratifying *this* envelope is itself a Tier-C act).
- **Cross-role conflict** - two roles' interests or claims collide.
- **Irreversible-outward beyond merge** - publish, delete, external send, close-as-completed, or any outward act that is not the already-gated merge.

When a trigger fires on an act, its default tier is overridden to C: stop and escalate (which is itself Tier A) before proceeding.

### 5.2 Reversible-but-novel: surface vs stop

"Novel" was doing two jobs in the original handoff. Splitting it resolves the tension with the existing `ask-for-help` rule:

- Reversible + **novel instance** (a new *case* of a known pattern) -> **Tier A + surface**: do it, and explain the reasoning.
- Reversible + **precedent-setting** (establishes a new *rule or pattern*) -> **Tier C**: stop, because a precedent is governance.

Worked example on the commit act:
- Commit a routine flow item on the active arc -> silent **Tier A**.
- Commit the first item of a brand-new work-type -> **Tier A + surface** (novel instance).
- Commit in a way that sets a new commitment *policy* -> **Tier C** (precedent).

## 6. The human-gate boundary (stated unambiguously)

> **Reversible coordination is pre-authorized; the irreversible outward act and the novel/governance decision are the human's.**

Concretely, the human gate is the union of:
- every **Tier-B** irreversible act (today: the `gh pr merge` command), and
- every **Tier-C** trigger (novel scope, precedent/governance change, cross-role conflict, irreversible-outward beyond merge).

Everything below that line - Tier A commit, dispatch, escalate, and all reversible preparatory work - is inside the envelope and needs no per-action approval. This line is the same one `decision-ratification != action-authorization` draws; the envelope just makes it complete and act-indexed.

## 7. Act tiers vs working modes (the two altitudes, side by side)

| | **Act tiers** (this spike) | **Working modes** (the epic) |
|---|---|---|
| Question | May I do *this act* unattended? | What *posture* do we hold over this *class of work*, and is it maturing? |
| Axis | Reversibility / blast-radius | Complexity x risk (paper), tracked as a Supervised -> Delegated trajectory |
| Granularity | Per atomic act (commit/dispatch/merge/escalate) | Per workstream, over time |
| Changes | Fixed policy | Climbs as trust/evidence accrues |

The tiers are the *mechanism*; the working-mode trajectory is *how far up the tiers we let a given workstream operate* as it earns trust. Same vocabulary root, different altitude - kept separate on purpose.

**Dispatch is the crossover case.** Its act tier is Tier B today and Tier A once matured - i.e. the working-mode trajectory is exactly what moves dispatch up the ladder. Every other act's tier is fixed; dispatch is the one place the two altitudes touch, which is why it is worth calling out rather than hiding as an ordinary fixed assignment.

## 8. How rungs 3-5 reference it

- **Rung 3 (auto-route pings -> dispatch):** a routing dispatch is **Tier B until dispatch matures, then Tier A** - proposed-and-confirmed at first, firing within the envelope once promoted, logging to the fire-log, never crossing the human gate.
- **Rung 4 (event-triggered Ready-pull):** the pull composes a **Tier-A commit** with a **dispatch** (Tier B until matured, then A); the pulled work then proceeds under its own act tiers (its eventual merge is Tier B).
- **Rung 5 (enveloped replenishment):** autonomous `Backlog -> Ready` is the **Tier-A commit** act at scale, bounded by the floor/cap and the surface-if-novel / precedent-stop refinements.

Each rung's own sub-issue cites this envelope as its governing policy; none adds autonomy past the section-6 boundary.

## 9. Landing plan (NOT yet executed - awaits ratification of this full picture)

Because wiring this in is itself a Tier-C governance change, it stops for human ratification of the assembled spec first. On ratification:

1. Add a `shared/` memory file (`feedback_coordinator_autonomy_envelope.md`) with the ladder + mapping + boundary, indexed in `shared/MEMORY.md`.
2. Add a concise CLAUDE.md section (the boundary statement + the mapping table), pointing at this spec for the full rationale.
3. Tick #1074's ACs; the epic's rung-3-5 sub-issues (Dev/Sci-led) reference the landed policy.

## 10. Open questions for review

- **Dispatch tier - RESOLVED (2026-07-08):** starts at **Tier B** (propose-and-confirm), promotes to **Tier A** once the review-axis rung has shipped and fire-log evidence shows safe behavior. Chosen over an immediate conditional Tier A because it honors the epic's binding constraint (widen review before automating dispatch) and makes dispatch the first worked example of the Supervised -> Delegated trajectory rather than a fixed assignment.
- **Close-as-completed**: folded under Tier-C "irreversible-outward beyond merge" today. A PR-merge close is already Tier B; only a *direct* `gh issue close` is the Tier-C case. Confirm that split reads cleanly.

**Created by:** PM
