# Issue #708 — single-feature KDE + centered-isotonic immunogenicity calibrator

**Date:** 2026-06-18
**Issue:** [#708](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/708) — `feat(scoring): single-feature KDE + centered isotonic calibrator + leave-one-cohort-out validation`
**Parent epic:** [#547](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/547) · **data dep:** [#707](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/707) (landed) · **downstream:** [#709](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/709) (Snakemake wiring, Developer)
**Author:** Scientist

## Goal

Reimplement NeoGuider's [Wei et al., *Genome Medicine* 2026] rank-calibration ML **from the paper** as a standalone Python module: adaptive kernel density estimation (KDE) of positive-vs-negative `genotype_presentation_score` densities → log density-ratio + prior → isotonic regression → **centered isotonic regression** (Oron & Flournoy 2017) → a `calibrated_immunogenicity_log_odds` column. Validate the centered-isotonic port against the R [`cir`](https://cran.r-project.org/web/packages/cir/cir.pdf) package, and assess generalization with **leave-one-cohort-out** validation on the harmonized training cohorts acquired in #707.

**Reimplement from the paper — do not vendor NeoGuider's code.** The [NeoGuider repo](https://github.com/XuegongLab/neoguider) is AGPL-3.0 + non-commercial and depends on paid-commercial tools (netMHCpan, netMHCstabpan, MixCR, PRIME). We use MHCflurry; none of those deps are pulled in.

## Non-goals (out of scope for #708)

- **Snakemake rule wiring** (`calibrated_immunogenicity_log_odds` into the pipeline) — Developer sub-issue #709.
- **Multi-feature extension** (foreignness, agretopicity, abundance, motif) — follow-up after the single-feature MVP.
- **patient_001 evaluation** — deferred until a validated cohort-level calibrator exists (single-patient calibration is meaningless; same shape as the AlphaGenome chr22 PoC #393).
- **Production-code change to `run_mhcflurry.py`** — the genotype-scoring formula is replicated (and test-pinned) in the precompute, not refactored out, to keep #708 off the production path.

## Approach (chosen: A — module in the epic folder + thin orchestrating notebook)

All work lands under the epic-keyed folder `research/experiments/issue_547_immunogenicity_calibration/`. The importable calibrator module + serialized artifact are written so #709 can consume them; promotion to `workflow/scripts/` is #709's decision. Rejected alternatives: (B) writing the production module straight into `workflow/scripts/` now — puts splice-unvalidated research code on the production path before any rule calls it and blurs the #708/#709 seam; (C) notebook-centric with no importable module — fails the "standalone module" + "serialized artifact for the wiring sub-issue" ACs.

## Components

### `calibrator.py` (importable, pure — no cohort I/O)

- **`centered_isotonic(x, y_iso, w) -> (cx, cy)`** — the CIR port. Collapse each flat level-set of the PAVA isotonic fit to its sample-size-weighted centroid (`cx = Σ x·w / Σ w` over the level-set), keep the isotonic `y`, and linearly interpolate between centred knots (Oron & Flournoy 2017). A from-scratch reference already exists in `research/evals/issue_258_neoguider/figures/_regenerate_figures.py`; this is the production-quality, validated version.
- **`class PresentationCalibrator`**
  - `fit(scores, labels)` →
    1. split `scores` by binary `labels` into positive / negative arrays;
    2. fit class-conditional KDEs — **adaptive variable-bandwidth (sample-point) KDE** (the AC-named method) **and** a fixed-bandwidth Silverman baseline, both retained for the diagnostics comparison;
    3. on a score grid, compute pointwise log density-ratio `log p(s|pos) − log p(s|neg)` and add an **explicit prior log-odds** → raw log-odds curve. **The prior uses the TRUE pre-subsample assayed counts** (`log(n_pos_true / n_neg_true)`), passed in by the caller — *not* the artificially-balanced subsample counts. The class KDEs estimate each density within-class (unaffected by how many negatives were scored), so combining within-class shapes with a true-base-rate prior makes the subsample unbiased for the **absolute** calibrated log-odds, not merely its shape;
    4. enforce monotonicity with `sklearn.isotonic.IsotonicRegression` (weighted by local sample count);
    5. smooth flat plateaus with `centered_isotonic` → store knots `(cx, cy)`.
  - `transform(scores)` → monotone interpolation of `calibrated_immunogenicity_log_odds` at the given scores; scores outside the fitted range clip to the boundary knots (documented behavior).
  - `save(path)` / `load(path)` → joblib. The artifact persists the **knots + prior + score-range + fit-cohort metadata + KDE-mode flag + code version** — not pickled KDE/estimator objects (lightweight, version-robust).

### `score_cohort.py` (the precompute — MHCflurry, one-time, cached)

1. Load NeoRanking `Neopep_data_org.txt` + IMPROVE `In_house_neoepitope_for_CV.tsv` (paths from `data_manifest.yaml`).
2. Filter to the **assayed** subset and map labels to binary: NeoRanking `response_type` `CD8`→1 / `negative`→0 / drop `not_tested`; IMPROVE `response` already 1/0.
3. Normalize HLA to MHCflurry canonical `HLA-A*01:01` (NeoRanking `A*01:01`, IMPROVE `HLA-A01:01` → canonical).
4. **Subsample (per the scoring-venue decision):** keep **all** positives (178 NeoRanking CD8 + 467 IMPROVE = 645); take a **cohort-stratified random sample of ~50K negatives** with a fixed seed. Rationale: `genotype_presentation_score` is 1-D on [0,1]; 50K negatives estimate that density essentially exactly. The subsampling risk only bites if positives are thinned, which we never do. **Record the true pre-subsample assayed positive/negative counts per cohort** in the parquet metadata so the calibrator's prior uses the true base rate (see `fit` step 3) — this makes the subsample unbiased for the absolute log-odds, not just its shape. The negative subsample size is a parameter; full-scale rerun is a config change if the curve proves sample-sensitive.
5. Compute `genotype_presentation_score` per peptide via `Class1PresentationPredictor`, scoring each peptide against its patient's ≤6 alleles (`HLA_allotypes.txt`) and combining with the **same `1 − ∏(1 − wᵢ·pᵢ)` HLA-C-weighted formula as `workflow/scripts/run_mhcflurry.py`** so the calibrator's input distribution matches production.
6. Write `outputs/scored_cohort_subsample.parquet` (columns: `peptide`, `cohort`, `label`, `genotype_presentation_score`, `best_presentation_percentile`); **skip-if-cached** keyed on the manifest checksums + subsample params.

### `notebook.ipynb` (orchestration)

Load the scored parquet, then:
- **Primary: leave-one-cohort-out (LOCO)** — 4 folds over `{NCI, TESLA, HiTIDE, IMPROVE}`; train `PresentationCalibrator` on three, evaluate on the held-out cohort. Per-cohort diagnostics: reliability curve (predicted log-odds bin → observed positive rate, with CIs), log-odds monotonicity check, and a ranking metric (AUPRC / positives-in-top-k).
- **Contrast: within-cohort k-fold** — stratified k-fold *within* each cohort, the optimistic in-distribution bound. Report the **LOCO-vs-within-cohort gap** per cohort = the cohort-shift penalty (the most decision-relevant number for #547's "trust on a new context?" question).
- **Adaptive-vs-fixed KDE comparison** — both bandwidth modes on the same folds, to see whether "adaptive" buys anything on this data.
- **Final fit** on all 4 cohorts → serialize `outputs/calibrator_v1.joblib`.

### Tests

- **`test_centered_isotonic.py`** — assert the port matches **frozen R `cir` fixtures** within tolerance. Fixtures generated once by running R `cir` (locally / throwaway script) on: the Oron & Flournoy 2017 worked example + edge cases (multiple level-sets, boundary plateaus, strictly-monotone no-collapse, single point). A committed generator script documents the R + `cir` version and the inputs; no standing R dependency in CI.
- **`test_calibrator.py`** — monotonicity invariant (`transform` non-decreasing in score), fit/transform output shapes, `save`/`load` round-trip equality, prior-odds handling, out-of-range clipping.
- **`test_score_formula.py`** — a few peptides scored by `score_cohort.py` match the documented `1 − ∏(1 − wᵢ·pᵢ)` HLA-C-weighted formula (faithfulness to `run_mhcflurry.py`).

## Data flow

```
raw cohort tables (gitignored; on this clone + GCS mirror)
   └─ score_cohort.py   [MHCflurry, one-time, subsample, cached]
        └─ outputs/scored_cohort_subsample.parquet
             └─ notebook.ipynb  → LOCO (primary) + within-cohort k-fold (contrast) + KDE compare
                  ├─ outputs/*.png             (reliability, monotonicity, shift-gap, KDE compare)
                  └─ outputs/calibrator_v1.joblib   → consumed by #709
```

## Validation design rationale (LOCO + the two folds)

**Why leave-one-cohort-out is primary:**

1. **Estimand matches deployment.** The calibrator runs on a distribution it never trained on (our splice neoepitopes, our HLA panel), so the target is *out-of-cohort* generalization. Pooled random k-fold leaks cohort-level confounds (assay, selection bias, HLA distribution, tumor type) across the train/test boundary and is optimistically biased.
2. **Non-IID group structure → grouped CV is the textbook remedy.** Cohort is a hard grouping variable (different assays: mini-gene IFN-γ ELISpot vs HLA-multimer vs ELISpot; different selection criteria, tumor types, label encodings). The standard correction for grouped data is blocked / grouped CV (`LeaveOneGroupOut`, `group = cohort`); same principle as blocked CV for spatially/hierarchically structured data (Roberts et al. 2017, *Ecography*).
3. **Field precedent — including the source paper.** Müller 2023 (the paper we port from) trains on NCI and tests on held-out TESLA + HiTIDE, explicitly to avoid learning cohort-specific selection bias ("ML classifiers could easily capture inherent biases related to the selection of neo-peptides for screening assays … This justifies our approach of training our classifiers on the NCI dataset"). The TESLA consortium and neoantigen-benchmark literature evaluate cross-dataset for the same reason.
4. **Conservative direction.** LOCO is the pessimistic bound relative to in-distribution k-fold, which matches #592's caution that splice transfer is unproven. Over-stating calibration quality is the dangerous error.

**Fold 1 — proxy-not-direct-test caveat (must be stated wherever results are cited):** none of the four cohorts is a splice cohort — all are SNV point-mutation neoantigens (splice-validation data is essentially nonexistent; the #680 sparsity problem). LOCO answers *"does calibration survive a real cross-lab / cross-assay distribution shift?"* — the best available stand-in — but is **not** a direct test of point-mutation → splice transfer. LOCO success must not be read as splice validation.

**Fold 2 — within-cohort k-fold contrast:** report stratified within-cohort k-fold alongside LOCO. The LOCO-vs-within-cohort gap quantifies the cohort-shift penalty and is nearly free to compute.

## Key decisions & risks

- **Class imbalance + negative subsampling:** handled by fitting each class KDE separately (within-class densities are unaffected by the subsample ratio) and folding an explicit prior log-odds computed from the **true pre-subsample** counts into the artifact metadata, so the subsample is unbiased for the absolute log-odds and #709 can re-anchor the prior to a different base rate. Surfaced, not hidden.
- **Small held-out folds:** TESLA (~34 pos) and HiTIDE (~41 pos) give wide CIs on per-fold reliability. LOCO trades per-fold precision for distributional honesty; CIs are reported, not hidden.
- **Output is log-odds, not probability:** the column name `calibrated_immunogenicity_log_odds` enforces the #592 caveat (presentation-anchored prior, conservative under-call on splice).
- **Environment (REVISED during implementation — two envs, not one):** the original plan was to add `mhcflurry` to `research/.venv`, but Python 3.14 has no TF/torch wheels for mhcflurry. Resolved by a **two-env split**: scoring (`score_cohort.py`) runs under a dedicated `mhcflurry-scoring` conda env (Python 3.11); the calibrator, notebook, and tests run under `research/.venv` (3.14). `score_cohort.py` imports `mhcflurry` lazily so the pure formula/HLA unit tests run under `research/.venv` with no mhcflurry. The **scored parquet is the clean handoff** between the two envs — arguably cleaner than one env. Documented in the #547 README "How to reproduce".
- **Label/HLA reconciliation across cohorts:** documented in the #547 README; encoded in `score_cohort.py`.

## Acceptance-criteria coverage

| AC (#708) | Where satisfied |
|---|---|
| Standalone adaptive-KDE + centered-isotonic module on `genotype_presentation_score` | `calibrator.py` |
| Centered isotonic in Python, validated vs R `cir` oracle | `test_centered_isotonic.py` + frozen fixtures + generator script |
| Leave-one-cohort-out validation (not within-cohort k-fold) | notebook LOCO (primary); within-cohort k-fold added as labeled contrast |
| Calibration-curve diagnostics per cohort (reliability + log-odds monotonicity) | notebook diagnostics |
| Output named as presentation-anchored prior, not validated probability | `calibrated_immunogenicity_log_odds` |
| Serialized fitted-calibrator artifact for #709 | `outputs/calibrator_v1.joblib` |
| Lab-notebook entry (per-cohort diagnostics) | `research/lab_notebook/scientist.md`, written post-review / pre-merge |

## References

- Wei et al., *Genome Medicine* 2026, [10.1186/s13073-025-01592-9](https://doi.org/10.1186/s13073-025-01592-9) — NeoGuider (Zotero `Z8FJSDVT`); method = custom KDE → centered isotonic regression.
- Oron & Flournoy, *Statistics in Biopharmaceutical Research* 2017, [10.1080/19466315.2017.1286256](https://doi.org/10.1080/19466315.2017.1286256) — centered isotonic regression.
- Müller et al., *Immunity* 2023, [10.1016/j.immuni.2023.09.002](https://doi.org/10.1016/j.immuni.2023.09.002) — NeoRanking harmonized cohorts (Zotero `TJ3RK87P`); cross-cohort hold-out precedent.
- Roberts et al., *Ecography* 2017 — cross-validation for non-independent (grouped/structured) data; blocked-CV rationale.
- [#592](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/592) — cohort + hold-out design decision (the splice-transfer caveat).
- [#119](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/119) — `genotype_presentation_score` (input feature).
- `research/experiments/issue_547_immunogenicity_calibration/README.md` + `data_manifest.yaml` — cohort schema, label encodings, HLA normalization, provenance.
