# Arc-Based Work Structuring — Design Spec

**Date:** 2026-06-05
**Status:** Approved — Plan 1 (foundation: labels + taxonomy + queryability) implemented in [PR #688](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/688); Plans 2–3 pending
**Issue:** tracked by [Issue #633](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/633) (`role:pm` — board governance review: parent-as-milestone-anchor, Milestones vs Iteration field, deferral tracking, WIP limits). This spec is the design output of that review's arc/milestone strand.
**Author:** PM session (Claude) + Jin-Ho Lee

> **Amended 2026-07-10 ([Issue #1102](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1102)):** the `active` slate is **no longer capped at 3**. It is an **advisory `~3` guideline** that prints a `NOTE` and never gates. `active` was measured not to predict throughput, and a WIP limit belongs on work in flight, not on categories. Where this spec says "cap `active` at 3" or "re-picks the ≤3 `active`" (§5.1, §6), read "aim for ~3, and never demote a productive arc merely to satisfy the number." The rest of the design stands.

## 1. Goal

Give the work a **narrative throughline that outlives a single milestone** and let the four asynchronous persona sessions (Dev, Sci, PM, MM) **pull from a shared current focus** without coupling them. Today work feels stand-alone day-to-day and the four sessions feel like separate workstreams — the symptom the user named.

The fix introduces one new dimension — the **arc** — and de-overloads the existing **milestone** so that three orthogonal axes each do exactly one job:

| Axis | Question it answers | Nature |
|---|---|---|
| **Stage** (`S#`) | "what *kind* of lifecycle work is this?" | general, fixed scaffold (DS lifecycle) |
| **Arc** | "which long-running concrete narrative does this serve?" | project-specific, fluid (NEW) |
| **Due date** | "when?" | the real chronological clock |

This spec is the *what/why*. The step-by-step build plan is produced next by the writing-plans skill.

## 2. Background — why now (the rationale, evidenced)

The current milestone scheme `i<iter> - S<N> - <Stage Name> - <Arc-descriptor>` encoded two best-practices: the **data-science lifecycle** as the stage axis (`S1 Domain Knowledge → S2 Data Collection → S3 Data Prep → S4 EDA → S5 Modeling → S6 Eval & Reporting → S7 Publication`, ≈ CRISP-DM/TDSP) and **iterating** through it as the `i<N>` axis. It was set up explicitly as a hypothesis to test and adjust. The verdict, from the live board:

- **The stage axis (S#) held up — keep it.** The arc partition test (§5) independently rediscovered it: the work genuinely clusters by lifecycle stage.
- **The iteration axis (i#) did not hold up as iterations.** An iteration is a *team-synchronization* device (one synchronized sprint, delivered at a boundary). The project runs **late-commitment Kanban** (flow-based, pull-when-ready) executed by **four async personas** — there is no synchronized cadence for an iteration to mean. Cut off from a global clock, `i<N>` mutated into a **per-track version counter** and fragmented into parallel counters (`i#`, `pm-i#`, `dev-i#`).
- **Evidence the number is not a clock:** sorting *open* milestones by `due_on`, the iteration numbers read `2, 5, 4, 5, 2, 6, 4` (and the `pm-i*` series reads `1, 3, 5, 4, 2, 6`). The field that *looks* chronological is not.
- **The arc is the iterative-DS-lifecycle idea expressed correctly.** "GTEx filter — pass 2, pass 3…" across iterations *is* iterating the data-prep stage on the GTEx theme. The intent was right; encoding it as a global synchronized sprint number (which fights async Kanban) was the part that broke. Expressed as a per-theme **arc lineage**, the same intent fits flow-based work.

What is being retired is **only the global `i<N>` as a fake clock** — not stages, not the DS-lifecycle intent, not Kanban flow.

## 3. How others handle "arc vs epic" (verified landscape)

A fact-checked sweep (Scrum/Cohn, SAFe, Jira, Linear, GitHub-native, Shape Up, Cagan/SVPG, Asana/ClickUp, ThoughtWorks) converged on:

- **The dividing line is closure, not size.** An **epic closes** when its children ship; a **theme/arc never closes** — successive epics keep rolling under it. *"If it ends when its children are done, it's an epic; if successive epics keep rolling under it, it's an arc."*
- **The long-lived tier is orthogonal, not a parent above epics.** Scrum (theme = label/grouping), SAFe (epics *tag* a Strategic Theme), Asana/ClickUp (rollup object), ThoughtWorks ("above nothing, beside everything") all model it as a cross-cutting lens. **Only Jira and Linear make it a literal parent tier** — Jira premium-gates it; Linear has no GitHub equivalent. Shape Up / Cagan say don't make it a tracked object at all.
- **A real 3rd structural rollup level is precisely what SaaS tools charge for** (automatic cross-level aggregation is expensive). The universal non-premium fallback is a label / field / single-select.

**Consequence for this repo (user project #9):** Issue Types are an organization-only feature and 404 on a *user* project, so a native "Arc" type is unavailable. The only way to fake a parent tier is a parent issue — which collides with the `recheck_parent_status` rollup (a never-closing arc-parent would emit perpetual forward-drift) and with the no-milestone-inheritance anchor model. **→ The arc must be an orthogonal label, not a tier.**

## 4. The model

### 4.1 Three things, three words already in the team's vocabulary

- **Epic** = a parent issue that **closes** when its sub-issues complete (unchanged — exactly what the repo already does).
- **Arc** = a never-closing `arc:<name>` **label** that any issue *or* epic carries, regardless of tree position. The narrative lens.
- **Milestone** = a dated, closeable bundle of "this stage-`S#` chunk of work." The time/stage slice.

Keep the word **"arc"** (it connotes a *narrative throughline* better than "theme"/"initiative", and the team already says it). Never call an arc an "epic that spans iterations" — that conflates a never-closing lens with a closeable container.

### 4.2 Stage × Arc is a grid, not a translation

Arcs are **not** the lifecycle stages relabeled. A real arc moves *through* stages; a stage *contains* many arcs. Two different arcs share a stage (e.g. Aligner and GTEx-filter are both `S3` but distinct threads) — that is the proof they are orthogonal axes. Drawing arcs at "upstream science / downstream science" granularity is the failure mode: at that coarseness the arc is just `S3`/`S5` with a new name and adds nothing. Arcs are drawn at **theme-lineage granularity** (§5).

### 4.3 Carriers — label-only MVP

1. **`arc:<name>` label** is the **load-bearing carrier.** The async cron/remote persona sessions read raw issue *text* (labels are in the REST issue payload; a Projects v2 board field is not). MVP = label only.
2. **`arc-phase:active | next | later` label** — the focus marker. **Aim for ~3 `active`** (advisory guideline since [#1102](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1102); originally a hard cap of 3).
3. *Optional later (only if missed):* a Projects v2 single-select `Arc` field (board grouping/views) and/or one never-closed `ARC: <name>` tracking issue per arc carrying the prose narrative. **Do not build all three carriers up front** — that is three sync-prone artifacts for one concept; start with the single text-readable label.

## 5. The arc taxonomy (validated against all 75 open issues, 2026-06-05)

Classification of every open issue into theme-lineage arcs. Largest arc is 23% — no catch-all (the coarse 4-arc version hit 42%; going *finer* fixed it).

| # | Arc (theme) | Threads | count |
|---|---|---|---|
| 1 | **Aligner & junction extraction** | STAR/HISAT2 verification & tuning, SJ/annotation correctness | 6 (8%) |
| 2 | **Junction filtering & tumor-specificity** | GTEx + matched-normal + AlphaGenome filter, integrity sweep | 6 (8%) |
| 3 | **Variant calling & cohort expansion** | somatic calling, SpliceAI/MMSplice prong, new patients | 6 (8%) |
| 4 | **Neoepitope scoring & TCR-pMHC modeling** | MHCflurry calibration, scorers, structures, immunogenicity | 7 (9%) |
| 5 | **Results, reporting & manuscript** | HLA-panel runs, report layer, decks, notebooks, writeup | 5 (7%) |
| 6 | **Cloud execution & reproducibility** | Batch, run registry, GPU fallback, `run_cloud_gpu.sh`, GCS/logs | 11 (15%) |
| 7 | **Memory & methodology** | `MEMORY.md` slimming/consolidation/audit, MM role, persona framing | 13 (17%) |
| 8 | **Board governance & enforcement** | Kanban governance, board hooks/guards, PM-tooling evals | 17 (23%) |
| — | **Unfiled** (not an arc) | dev-env hygiene + arc-agnostic refactors | 4 (5%) |

Membership (bare `#N` — code block):

```
1 Aligner & junction extraction:    297 375 377 378 411 636
2 Junction filtering:               126 212 304 381 594 663
3 Variant & cohort:                 413 416 436 437 438 440
4 Scoring & TCR-pMHC:               433 492 547 566 585 601 659
5 Results/reporting/manuscript:     193 198 233 435 455
6 Cloud execution & reproducibility:66 183 195 310 630 651 658 664 669 673 674
7 Memory & methodology:             248 265 324 326 346 353 527 538 539 540 541 542 672
8 Board governance & enforcement:   234 250 294 295 406 445 498 499 533 553 569 578 617 626 633 655 665
Unfiled:                            446 451 570 641
```

Two mirrors this surfaces, both worth a conscious decision rather than drift:

- **Science (Arcs 1–5) = 30 issues ≈ Persona-OS (Arcs 7–8) = 30 issues.** Half the open backlog is the personas building their own machinery. Investment, not waste — but now visible.
- **Arcs 7+8 are really one milestone lineage** (the `pm-i*` "PM Tooling, Memory & Methodology") that outgrew one thread; split at the seam memory/methodology vs board/enforcement. **First fission candidate** to watch.

### 5.1 Suggested opening active slate (PM call, not fixed)

Committed milestones currently touch **six** arcs at once — that scatter is the "too many parallel things" symptom, now measurable. A focused slate of three, one per persona-cluster:

- **Arc 1 · Aligner/STAR verification** — closest-due milestone ([Issue #24](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/24), 2026-06-19); gates downstream science; Sci+Dev+PM.
- **Arc 4 · Scoring & TCR-pMHC** — active modeling frontier ([Issue #29](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/29)); Sci+Dev.
- **Arc 8 · Board governance & enforcement** — carries [Issue #633](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/633) (this design) + the Kanban migration; PM+MM+Dev.

`next`: **Arc 6 (cloud)** — runs gate the science. Everything else `later`.

## 6. Arc lifecycle & governance

Arcs are the **most fluid** of the three axes by design. "Never closes" ≠ "fixed":

- **Closing** is mechanical (an epic closes when children ship). An arc has no such trigger.
- **Retiring** is editorial (a human decides the narrative has run its course). Arcs *can* be retired — Cagan's anti-immortal-bucket warning requires it.

**Stability gradient:** stage = fixed scaffold · arc = fluid narrative layer · due date = pointwise. The taxonomy in §5 is **v1, expected to evolve.**

**Operations are cheap because the carrier is a label** (this is a direct payoff of the label-only choice — a parent-issue or rigid field hierarchy would make these painful migrations):

| Operation | Mechanic |
|---|---|
| Rename | `gh label edit "arc:old" --name "arc:new"` — propagates to all issues |
| Split (fission) | create new label(s), re-tag the subset, delete old (as done for 7/8) |
| Merge (fusion) | re-tag issues `arc:B → arc:A`, delete `arc:B` |
| Retire | drop from the active/next/later rotation; optionally `arc:archived/<name>`; history preserved |
| Spawn | new thread emerges → new label, starts `later` |

**The one discipline — change on a cadence, not continuously.** Stable between reviews (that stability *is* the throughline); revised at a periodic **arc review** (monthly, or when the active slate feels stale) that deliberately decides splits/merges/retirements and re-picks the `active` slate (aim ~3, advisory since [#1102](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1102); promote on evidence, demote only a finished or stalled arc). This mirrors Shape Up's "re-bet each cycle." **PM-coordinated**, like the commitment act (portfolio-level call). The `arc-phase` markers churn fast (every review); the arc taxonomy churns slowly.

**The Unfiled lane is the nursery and the graveyard:** new arcs gestate there (when several Unfiled issues share a thread, an arc is born); a retired arc's stragglers drop back into it. **Never fold Unfiled items into the biggest arc** — that re-bloats it toward catch-all.

## 7. Milestone de-overloading

The milestone keeps its **correct, idiomatic GitHub use** (a dated, closeable deliverable bundle) and sheds its two passengers:

- **Drop the arc-descriptor from the name** (the arc is now its own label). Milestone name → just the stage chunk.
- **Drop / demote the `i<N>` global clock.** The real clock is `due_on` (already synced to the board `Target`). If a true time-box is ever wanted, use the **native Projects Iteration field**, not a name suffix — the open question in [Issue #633](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/633).
- **Per-arc "version" is expressed as the milestone's place in its arc lineage** ("the 5th milestone under `arc:aligner`"), not a global `i5`.
- **Cleanup:** delete dead legacy milestones (`M1 / M2 / M7` — empty, no due date); close/populate empty placeholder milestones.

Milestones remain the **commitment signal** (`Backlog → Ready`) and `Target` source — unchanged. One issue ↔ one milestone, sub-issues don't inherit — unchanged.

## 8. The real deliverable — the arc-aware daily pull

Labels are inert metadata until the sessions *use* them. The point of the whole exercise:

> Each persona session's pull rule becomes: **`status:Ready` AND `arc-phase:active` AND `role:<self>` → take the next.**

This is what turns four async sessions into a coordinated team pulling a shared current focus. It requires updating each persona's pull/cron prompt (PM, Sci, Dev, MM) + the morning routine to surface the active slate and each role's next pull. **This is the primary implementation work, not a downstream consequence.**

## 9. Rollout (phased, low-risk)

- **Phase 0** — create 8 `arc:*` labels + 3 `arc-phase:*` labels; batch-apply to the 75 open issues from §5; open an explicit Unfiled state for the ~4 arc-less.
- **Phase 1** — set the opening `active` slate (§5.1); surface it in the PM morning routine board-hygiene step (an un-arced `Ready` issue is a triage gap, like a missing role label).
- **Phase 2** — teach each persona pull/cron the arc-aware filter (§8). *The core deliverable.*
- **Phase 3** — milestone cleanup (§7): drop arc suffix, delete `M1/M2/M7`, resolve Milestones-vs-Iteration-field for [Issue #633](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/633).
- **Phase 4** — first **arc review** (§6) one month in: validate the taxonomy, split Arc 8 if it's grown, retire anything concluded.

## 10. Risks & mitigations

| Risk | Mitigation |
|---|---|
| **Label drift** — issues filed without an arc; arcs going stale | Morning board-hygiene sweep flags un-arced `Ready` issues; arc review refreshes the slate |
| **Arc 7/8 (persona-OS) re-bloating** | Unfiled lane + the 7/8 split are the guardrails; it's the named first fission candidate |
| **Cross-cutting issues** (~9, e.g. #665, #193, #455) | Resolve by *payload not label*; accept ≤10% need a judgment call |
| **Carrier sync** — a Projects field is not in issue text | Label is canonical for text-reading sessions; defer the field until needed |
| **Over-churn** erases the throughline | Change only at the cadenced arc review, not mid-session |
| **The 53%-meta mirror** | The structure makes it visible; that may itself prompt a deliberate rebalance |

## 11. Out of scope / explicitly NOT doing

- **No 3-level arc → epic → issue hierarchy** (premium-gated, blocked on a user project, collides with the rollup hook).
- **No per-arc tracking issue or Projects field in the MVP** (label-only first; add only if missed).
- **No change to epics, sub-issue non-inheritance, the commitment act, or `Target` sync.**
- **No retroactive renaming of closed milestones** (cosmetic; cleanup applies to dead/empty ones only).

## 12. Open questions (arc-rollout details resolved across Plans 2-3; the [Issue #633](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/633) governance remainder — parent-anchor / deferral / WIP — is tracked in [Issue #690](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/690))

1. **Milestones vs native Iteration field** — once `i<N>` leaves the name, do we adopt the Projects Iteration field for real time-boxes, or rely on `due_on` alone?
2. **`arc-phase` as label vs Projects single-select** — label is text-readable (MVP); a single-select is tidier on the board. Decide at Phase 1.
3. **Arc review cadence** — monthly vs per-"replenishment" — tune after Phase 4.

## 13. References

- Brainstorming dialogue 2026-06-05 (this session): framework landscape sweep + 75-issue partition test (theme granularity).
- [Issue #633](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/633) — board governance review (tracking issue for this strand; a `role:pm` leaf, not an epic).
- `CLAUDE.md` §"Board status governance — late-commitment Kanban" and `.claude/memory/feedback_milestones.md` — current milestone/commitment rules this spec amends.
