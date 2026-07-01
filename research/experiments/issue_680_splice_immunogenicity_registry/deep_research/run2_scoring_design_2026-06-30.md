# Deep-research run #2 - #736 scoring-harness DESIGN brief

**Date:** 2026-06-30 (~23:55 UTC) · **Run ID:** wf_0e9ba4be-d18 · **Task:** wh236wtj2
**Raw output (durable):** `tasks/wh236wtj2.output` (session f1a362d9…)
**Scale:** 27 agents, ~1.6M tokens, ~14 min. Completed locally (no stall/resume needed).
**Pipeline:** frame (live counts) -> 5 candidate scoring methodologies -> adversarial 4-lens judge -> synthesize.

## Judged ranking
1. **uncertainty-first (7.5)** - winner (best on statistical-validity, reviewer-proofness, honesty)
2. positive-only-enrichment (7.5)
3. ranking-only (7.0)
4. tiered-stratified (6.75)
5. negative-sensitivity (6.0)

## Recommended design (winner + grafts)
**#736 ships as a ranking/enrichment benchmark on HLA-A\*02:01 positives, with NO scalar output.** The deliverable is a **power ledger + a facet grid**, not "the predictor scores 0.NN."

1. **HEADLINE = power ledger, computed first, on EFFECTIVE (cluster-deflated ~3.84) sample sizes.** Red light: 9 Tier-1 negatives vs 19-31 needed -> specificity **UNMEASURABLE**; negative-set expansion named as the blocking prerequisite.
2. **Non-pooling is the SHAPE of the output** - long-form rows keyed (metric x neg_tier_set x allele_stratum x evidence_stratum); a `pool_guard` REFUSES any unregistered cross-tier union in code (LABELING_SCHEME sec-7 enforced, not remembered). `power_flag = POWERED | UNPOWERED | N/A-structural`.
3. **Primary positive readout = positive-rank-percentile ECDF (A2-only)** - threshold/ratio/k-free; displaces the winner's confounded Tier-3 enrichment (grafted from positive-only-enrichment).
4. **Tier-3 demoted to a labeled presentation floor**, two sub-schemes: T3-shuffle (composition floor, "measures presentation NOT immunogenicity") + T3-dude (DUD-E-style presentation-matched decoys = the honest "floats above presentation" test). Rank-concordance C never relabeled AUROC.
5. **Uncertainty: LOSO first** - leave-one-study-out fold-range as the PRIMARY interval (drop-Bigot swing = the instability); cluster bootstrap relabeled "coverage-unknown"; i.i.d. peptide bootstrap only as a labeled foil; cluster-restricted permutation; NO permutation p on T3-shuffle (exchangeability violated).
6. **Discrete-negative cells (T1/T2) report effect size, not AUC** below the effective-power floor - Cliff's delta / rank-biserial + AUPRC-with-prevalence-baseline. Tier-1 (VELEDHVML A\*11:01 + 8 soft) is a permanent UNPOWERED directional anecdote; VELEDHVML quarantined as off-A2-manifold.
7. **Strong-vs-weak** = internal assay-modality consistency check (NOT calibration); check source-confounding first.
8. **Allele: A2-only stamp everywhere.** Non-A2 (<=4 each) appendix-ranks only; A2 x Tier-1-hard specificity cell printed empty (n=0), never borrowing the A\*11:01 negative.

## Claims FORBIDDEN (first-class manifest artifact)
No specificity/ROC-AUC/FPR/PPV on T1 or T1+2 (underpowered -> returns UNPOWERED); no pooled cross-tier negative metric; no point estimate without interval+verdict; no nominal 95% CI from a ~3.8-cluster bootstrap; no i.i.d. peptide bootstrap as a result; no per-allele non-A2 metric; no "T3 enrichment = specificity"; no scalar headline; no calibrated immunogenicity probability; no max-min cross-scheme swing headline; never count the 6 non-scorable rows into N.

## 7 open questions for Jin-Ho
1. Confirm THE headline = red-light power ledger, with the A2 rank-percentile ECDF as the primary positive readout beneath it. (Recommended yes.)
2. Tier-3 v1 scope: composition-shuffle only, or invest in the presentation-matched T3-dude path now? (T3-dude needs a sanctioned non-immunogenic presented-background pool - which one?)
3. **Which predictor score does the harness score with: raw MHCflurry `presentation_score`, or `calibrated_immunogenicity_log_odds` (wired via #907)?** Changes what "the predictor" means in every claim.
4. Power-ledger assumptions: target AUC (default 0.75), effective-power floor, confirm the 19-31 negatives-needed band.
5. Cross-scheme sensitivity in v1, or defer until >1 predictor to compare?
6. File the negative-set-expansion data ask as a board Issue now (Tier-2 functional testing + more hard true-negatives, non-A2-balanced)? Overlaps **#911**.
7. Confirm the canonical-source collapse map (SNAF x2, IRIS x2) as a shared #737 artifact - effective-source 3.84 is load-bearing for the ledger + every interval.

## 11 implementation next steps (for #736)
Canonical-source map (assert 10 sources/3.84 eff/Bigot 43%); loader + assert tripwires (scorable_pos==81, A2==72, evidence 40/41 full & 31/41 A2-only as SEPARATE guards, exclude 6 non-scorable, T1==9, T2==13 blank-HLA); SCHEME_REGISTRY + pool_guard (unit-tested to raise); predictor-IN-LOOP scoring bridge; power-ledger module FIRST on deflated eff-n; metric layer; uncertainty stack (LOSO primary); Tier-3 generator (v1 shuffle, stub T3-dude, seeded RNG in provenance); emit 3 artifacts (facet-grid TSV + results.json + human Markdown grid) + the claims_forbidden manifest; unit tests; file the negative-set-expansion Issue + back all deck decisions on the board before the #736 deck PR merges.
