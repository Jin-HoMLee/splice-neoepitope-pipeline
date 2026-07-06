<!-- Author-owned narrative for i6 - S3 - Data Preparation. Sections 3/4/5 only.
     The script regenerates the HTML from this file + fresh board data;
     it never overwrites this sidecar once it exists. -->

## Deliverables (Review layer)

This milestone is the **data-preparation stage for the splice-immunogenicity ground-truth**: building, cleaning, and characterizing the curated registry of experimentally-validated splice-derived neoantigens that downstream calibration and benchmarking depend on.
All 10 issues were Scientist-led (one dual Sci/Dev). Two clusters are the milestone's core - the registry and its labeling/characterization; three further issues delivered supporting data-substrate and repo-hygiene work that closed under the milestone but sit adjacent to the ground-truth-preparation thesis.

### 1. The splice-immunogenicity registry (validated ground-truth substrate)

The core artifact of the stage: a provenance-tracked registry of real, T-cell-validated splice neoantigens, grown by systematic multi-database mining rather than opportunistic collection.

- **Registry growth from standing-watch + library sweep** - [Issue #733](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/733) ([PR #819](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/819)): expanded the seed registry to 44 rows.
- **4-DB splice-category gap analysis + CEDAR/IEDB free-text mine** - [Issue #734](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/734) ([PR #840](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/840)): a systematic recovery pass across four databases to find splice-category candidates the structured queries miss.
- **Verify + fold the recovered candidates** - [Issue #838](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/838) ([PR #905](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/905)): confirmed and folded the #734 IEDB-recovered set (long-read UM, Kwok GNAS/RPL22, POSTN, Kim extras) into the registry.
- **Data-integrity close-out** - [Issue #904](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/904) ([PR #913](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/913)): re-verified the two IRIS `candidate` rows against the located IRIS supplement, resolving the last open AC on #733.
- **Provenance + data-quality pass** - [Issue #823](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/823) ([PR #893](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/893)): added the `assay_context` column and upgraded Manoharan provenance, so each entry carries how its immunogenicity was established.

### 2. Labeling scheme + ground-truth sparsity characterization

Turning the registry into a usable train/eval substrate, and honestly quantifying how much truth it actually holds.

- **Junction mapping + documented pos/neg labeling scheme** - [Issue #735](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/735) ([PR #881](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/881)): mapped registry entries to junctions and documented the positive/negative labeling convention. The registry is a **two-class** substrate, not a positive-only list: alongside the T-cell-validated positives it carries curated negatives, including the scarce true-splice-junction-but-non-immunogenic hard-negatives (functional screens with no measured T-cell response) that a calibrator needs to learn the decision boundary.
- **Ground-truth sparsity quantification + writeup** - [Issue #737](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/737) ([PR #912](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/912)): quantified how sparse the validated ground truth is and shipped the writeup + experiment deck - establishing that the current set is a face-validity / stress-test probe rather than a statistically powered benchmark, and characterizing that limit head-on instead of papering over it.

### 3. Supporting / adjacent work

These three issues closed under i6-S3 but sit adjacent to the ground-truth-preparation thesis rather than delivering it. Two are reactive flow/hygiene work that our three-axis model commits milestone-free; the cohort refresh is genuine pipeline-data work, but on a different data axis than the registry.

- **STAR cohort re-run + RESULTS refresh** - [Issue #636](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/636) ([PR #943](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/943), dual Sci/Dev): re-ran patient_001 + patient_002 through the STAR path and refreshed `RESULTS.md` on corrected data. Pipeline-data work - the closest of the three to the data-preparation stage, though on the patient-cohort axis, not the ground-truth registry.
- **Widen Zotero capture to field-context + orphan notes** - [Issue #634](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/634) ([PR #670](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/670)): closed the retired news-log's lost "else" bucket by routing field-context updates into Zotero. Curation-process infrastructure.
- **Migrate per-Issue notebooks + slides into the `research/experiments/` convention** - [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455) ([PR #775](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/775)): consolidated scattered analysis artifacts under one convention. Repo hygiene.

## Carried-forward & routing

**Nothing carried forward** - 10 of 10 delivered, 0 descoped, 0 re-milestoned.

**Aggregate-outcome routing (milestone-close decision).**
The stage produced a standing artifact, not a one-off: the splice-immunogenicity registry (with its labeling scheme and characterized sparsity) is now the ground-truth substrate the downstream work consumes.
Routing is **(c) extend the workstream** - the data-prep outputs flow directly into the live `arc:immunogenicity-benchmark` program:

- The **simulated validation dataset** family ([Issue #1036](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1036) + sub-issues #1037/#1038, triaged 2026-07-06) draws its injected loci and two-tier truth labels from this registry.
- The **first open head-to-head caller benchmark** ([Issue #679](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/679)) uses the registry as its validated-positive set.
- The immunogenicity **calibrator** delivered in the parallel i5-S5 Modeling milestone is validated against this same registry.

No Publication (a) or Modeling (b) sibling is opened off this milestone directly - S3 is a foundation stage; its scientific payoff is realized in the benchmark + modeling work it feeds.

## Retrospective (process/health)

**Healthy close.**
10/10 delivered over a ~4-week window (2026-06-04 to 2026-07-02), 0 descoped, median cycle time 10 days (avg 12.4) - steady single-lane Scientist throughput with no stranded or aging items.

**What went well.**

- **Disciplined data curation.** The registry grew through systematic multi-DB mining (#733/#734/#838) with explicit provenance rigor (#823) and a data-integrity close-out that chased down the last unverified rows (#904) - not opportunistic collection.
- **Scientific honesty on the binding constraint.** #737 quantified the ground-truth *sparsity* head-on and shipped the characterization, rather than hiding a small validated set. That honesty is what makes the downstream benchmark trustworthy.
- **Clean handoff.** The stage leaves a labeled, provenance-tracked substrate that the benchmark and validation-dataset work can consume without rework.

**What to improve.**

- **Milestone composition swept in off-thesis work.** 3 of the 10 issues sit adjacent to the ground-truth-preparation thesis: [Issue #636](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/636) (cohort refresh, a different data axis) and especially [Issue #634](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/634) (curation infra) + [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455) (repo hygiene). By our three-axis model, #634 and #455 are flow/hygiene work that commits **milestone-free** - they should not have carried the S3 lifecycle milestone. Fix forward: at `Backlog -> Ready`, gate an `i<N> - S<N>` milestone to that stage's thematic deliverables only, and route flow/hygiene/tooling work milestone-free. All 10 are kept in this report (accurate to what was milestoned + merged) with the three reframed as "Supporting / adjacent." (Surfaced by the Scientist author pass.)
- **The report scaffolding itself was noisy.** `milestone_report.py`'s first-run seed over-collected - it pulled lab-notebook entries across the whole window and both roles instead of scoping to this milestone's 10 issues, so the Deliverables narrative had to be re-authored from the verified issue->PR set by hand. That is exactly the open [Issue #1005](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1005) (seed over-collects; scope Deliverables to the milestone's issues + lead role). This close is a second concrete datapoint for prioritizing #1005 before the next lifecycle-milestone report.
