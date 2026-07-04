<!-- Author-owned narrative for i5 - S5 - Modeling. Sections 3/4/5 only.
     The script regenerates the HTML from this file + fresh board data;
     it never overwrites this sidecar once it exists. -->

## Deliverables (Review layer)

<!-- Lead role: what shipped, grouped by deliverable, with PR + slide links.
     The 5 delivered issues are listed in the Inventory appendix
     below - narrate the highlights here, don't re-list them. -->

The milestone's throughline was the **immunogenicity calibrator**: turning MHCflurry's `presentation_score` (a presentation likelihood) into a **calibrated immunogenicity log-odds**, carried end to end from training-data selection through production pipeline wiring.
Four of the five deliverables are that chain; the fifth stands up our closest-peer cross-check.

### Immunogenicity calibrator (4 of 5 deliverables)

The calibrator moved from paper to pipeline in four gated steps:

- **Training-cohort selection** ([Issue #592](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/592), [PR #596](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/596)) landed the decision on *which* labelled immunogenicity cohorts to calibrate against, choosing four SNV / point-mutation neoantigen cohorts (NCI, TESLA, HiTIDE, IMPROVE).
  The choice was grounded in a fact-checked 109-agent deep-research sweep (26 sources, 80 claims extracted, 25 adversarially verified) rather than recall, and it satisfied AC 1 of the parent calibrator epic [Issue #547](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/547), unblocking the fit.
  Decision deck: `research/decisions/issue_592/slides.html`.
- **The calibrator itself** ([Issue #708](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/708), [PR #786](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/786)) reimplemented NeoGuider's rank-calibration ML *from the Wei et al. (Genome Medicine 2026) paper* rather than vendoring its code (NeoGuider is AGPL-3.0 with an explicit non-commercial restriction and depends on paid tools we do not use).
  The module chains adaptive KDE of positive-vs-negative `presentation_score` densities into a density ratio, then isotonic and finally **centered isotonic regression** (Oron & Flournoy 2017), emitting a calibrated immunogenicity log-odds and validated leave-one-cohort-out (LOCO).
  Ships the fitted `calibrator_v1.joblib`; experiment at `research/experiments/issue_547_immunogenicity_calibration`.
- **Applicability gate** ([Issue #826](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/826), [PR #873](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/873)) is the honest guardrail on the SNV-to-splice transfer.
  No labelled splice-neoepitope immunogenicity data exists yet (that gap is [Issue #680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680)), so a full validation cannot run; this is the *label-free* applicability check that can.
  It quantifies the fraction of real patient_001 / patient_002 splice `presentation_score`s falling outside the calibrator's fitted support - where the curve flat-line extrapolates and the log-odds is unreliable - and surfaces those out-of-support presenters instead of silently trusting them.
  This gate had to pass before the production wiring below.
- **Pipeline wiring** ([Issue #709](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/709), [PR #907](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/907)) is the single Developer deliverable: a new Snakemake rule between MHCflurry presentation and TCRdock that applies the serialized calibrator and emits the `calibrated_immunogenicity_log_odds` column, making the calibrated signal a first-class pipeline output.

### Peer cross-check setup

- **ASNEO cross-check - setup + decision slice** ([Issue #566](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/566), [PR #849](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/849)).
  ASNEO (Zhang et al., *Aging* 2020) is our closest published peer (RNA-seq to splice junctions to neoepitopes).
  This deliverable resolved the open MHC-path question with an **option-B decision** (standardize both pipelines on MHCflurry so call-concordance isolates the junction-detection signal and holds the MHC predictor constant as a nuisance variable) and shipped a **turnkey chr22 runner** so the actual cross-check is one command.
  The VM-bound execution (run, concordance, writeup) is carved to [Issue #848](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/848); dir `research/experiments/issue_566_asneo_crosscheck`.

## Carried-forward & routing

<!-- PM: issues that didn't close + where they went (carve / arc) + the
     closure-routing decision (a/b/c/d). -->

None carried forward - all 5 issues delivered COMPLETED and i5-S5 closes clean. The sibling Data-Preparation milestone i6-S3 is a **separate** carve/decommit decision, tracked via [Discussion #964](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/discussions/964) (not part of this milestone's scope).

## Retrospective (process/health)

<!-- PM: was it healthy? what to improve? WIP/aging observations. -->

i5-S5 closed **healthy**: 5/5 delivered (4 Scientist + 1 Developer), **0 carried forward**, cycle time avg 11.2d / median 8.3d. The Scientist calibrator workstream (KDE + centered isotonic + LOCO validation, then within-cohort CV, Jeffreys/boundary CIs, monotone extrapolation, and report surfacing) drove the milestone; the Developer contribution wired the calibrated immunogenicity log-odds through the pipeline into the report. No aging-WIP or blocker incidents inside the milestone window.

One process note routed to a follow-up: the closure-report **seed over-collects** (dumps the whole date-window across all roles rather than the milestone's issues), so this Deliverables section needed heavier hand-pruning than it should - scoped as a fix in [Issue #1005](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1005).
