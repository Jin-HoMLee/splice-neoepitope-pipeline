# issue_547 — immunogenicity calibration

**Goal:** build the post-MHCflurry **immunogenicity calibration** layer (KDE + centered isotonic regression) that maps `presentation_score` → `calibrated_immunogenicity_log_odds`. This folder is the **epic-level home** for that work; the first sub-issue ([Issue #707](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/707)) acquires the labelled training cohorts chosen in [Issue #592](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/592) and confirms their schema is usable. Nothing downstream (#708 calibrator, #709 Snakemake wiring) starts until the tables are in hand and parseable.

**Parent epic:** [Issue #547](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/547) · **acquisition sub-issue:** [Issue #707](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/707)

Folder is **epic-keyed** (`issue_547`, not `issue_707`) because the cohort data + calibrator are shared across the sub-issues; it mirrors the GCS prefix `gs://splice-neoepitope-project/experiments/issue_547/`.

## Layout

```
issue_547_immunogenicity_calibration/
├── README.md                  # this file (committed)
├── data_manifest.yaml         # pinned schema v1 — paths + checksums + fetch cmds + license (committed)
├── calibrator.py              # PresentationCalibrator (adaptive+fixed KDE → density-ratio + true-base-rate
│                              #   prior → isotonic → centered isotonic) + centered_isotonic();
│                              #   .fit / .transform / .save / .load (committed)
├── score_cohort.py            # MHCflurry precompute (run with `mhcflurry-scoring` conda env):
│                              #   genotype_score, normalize_hla, build_scored_cohort; has __main__ (committed)
├── notebook.ipynb             # LOCO + within-cohort k-fold validation + diagnostics + final artifact fit (committed)
├── applicability_notebook.ipynb   # Issue #826 label-free SNV→splice applicability gate (score-support + covariate-shift) (committed)
├── tests/                     # unit tests (committed)
│   ├── test_centered_isotonic.py  # vs frozen R cir package fixtures
│   ├── test_calibrator.py         # fit/transform/save/load round-trip + prior injection
│   ├── test_score_formula.py      # genotype_score + normalize_hla
│   ├── generate_cir_fixtures.R    # R oracle script that produced fixtures/cir_fixtures.json
│   ├── fixtures/
│   │   └── cir_fixtures.json      # frozen R oracle outputs (committed)
│   └── conftest.py
├── outputs/                   # calibration artifacts + diagnostics (see Outputs index below)
│   ├── calibrator_v1.joblib           # deliverable artifact (committed)
│   ├── pr_reliability.png             # logit-scale calibration panels (45° + Cox fit) (committed)
│   ├── shift_gap.png                  # LOCO vs within-cohort gap (committed)
│   ├── kde_compare.png                # adaptive vs fixed KDE (logit-scale reliability) (committed)
│   ├── kde_fit_compare.png            # adaptive vs fixed KDE — densities + resulting calibration map on pooled data (committed)
│   ├── method_fitted.png              # the actual fitted pipeline on our pooled cohorts (KDE→ratio+prior→isotonic→CIR→transform) (committed)
│   ├── auprc_explained.png            # teaching: how AUPRC is built (ranking→precision/recall→area), real PR curves TESLA vs IMPROVE (committed)
│   ├── concept_presentation_funnel.png         # teaching schematic: presentation ≠ immunogenicity (committed)
│   ├── concept_discrimination_vs_calibration.png  # teaching schematic: the two eval axes (committed)
│   ├── panel_anatomy.png              # teaching schematic: how to read a calibration panel (committed)
│   ├── conclusion_scorecard.png       # conclusion: per-cohort discrimination × calibration verdict (committed)
│   ├── conclusion_ladder.png          # conclusion: evidence ladder (within → LOCO → splice = untested) (committed)
│   ├── scored_cohort_subsample.parquet        # GITIGNORED — LICR-derived subsample
│   └── scored_cohort_subsample.parquet.true_counts.csv  # prior counts per cohort (committed)
└── data/                      # raw cohort tables — GITIGNORED; not on fresh clones
    ├── neoranking/            #   primary: Neopep + Mutation + HLA_allotypes (+ source .zip)
    └── improve_borch/         #   augment: In_house_neoepitope_for_CV + CEDAR benchmark
```

Raw tables are **never committed**: NeoRanking is LICR copyright (no redistribution) and the tables are large. `data/` is gitignored repo-wide; recreate it with `mkdir -p data/{neoranking,improve_borch}` after clone and fetch per `data_manifest.yaml`. The GCS mirror at `gs://…/experiments/issue_547/` uses the **same two subdirs** (`neoranking/` + `improve_borch/`).

## Status (2026-06-18)

| Step | State |
|------|-------|
| Sources located + recorded | ✅ figshare share links + filenames + LICR license (`data_manifest.yaml`) |
| Byte-download (browser route) | ✅ 3 NeoRanking tables downloaded + unzipped into `data/` |
| Checksums + sizes | ✅ sha256 + byte sizes in `data_manifest.yaml` |
| Schema + join-feasibility doc | ✅ see "Schema & join feasibility" below |
| GCS mirror (>100 MB) | ✅ `gs://splice-neoepitope-project/experiments/issue_547/` |
| Cohort-composition reconcile vs #592 | ✅ paper+data confirmed: mutation 131 pt / neo-pep 99 pt, 178 CD8 neo-peptides, **no Bjerregaard** — #592 correcting note posted |
| IMPROVE/Borch augment | ✅ downloaded — CV training table (17,520 / 467) + CEDAR benchmark from `SRHgroup/IMPROVE_paper` (the repo named in the paper's Data Availability Statement) |
| MHCflurry scoring (`score_cohort.py`) | ✅ all 4 cohorts scored; `scored_cohort_subsample.parquet` + `.true_counts.csv` in `outputs/` |
| Calibrator build (`calibrator.py` + `notebook.ipynb`) | ✅ `calibrator_v1.joblib` in `outputs/`; LOCO + within-cohort validation done — see "Validation result" below |
| Unit tests | ✅ `tests/` — centered-isotonic vs R oracle, calibrator round-trip, score formula |

## Schema & join feasibility (NeoRanking, confirmed from the downloaded tables 2026-06-17)

Two parallel TSVs, same patient cohort at two granularities; both keyed on `patient` + the
mutation/peptide identity columns. HLA per patient in `HLA_allotypes.txt`.

| table | rows | cols | grain |
|---|---|---|---|
| `Neopep_data_org.txt` | 1,787,710 | 57 | one row per candidate **8–12mer neo-peptide** |
| `Mutation_data_org.txt` | 48,306 | 59 | one row per **mutation** (long/mut-seq level) |
| `HLA_allotypes.txt` | 150 patients | 2 | `Patient` ⇥ comma-separated 4-digit alleles |

**Columns that matter for calibration** (neo-peptide table):

| role | column | notes |
|---|---|---|
| peptide | `mutant_seq` | the 8–12mer; 0 missing |
| immunogenicity label | `response_type` | `CD8` = positive, `negative` = assayed-no-response, `not_tested` = exclude |
| HLA allele | `mutant_best_alleles` | MixMHCpred best allele for the peptide; 0 missing |
| cohort | `dataset` | `NCI` / `TESLA` / `HiTIDE` |
| split | `train_test` | author train/test partition |

**Label encoding** (neo-peptide table): `not_tested` 1,364,625 · `negative` 422,907 · **`CD8` 178**.
The calibration target is the assayed subset (`CD8` ∪ `negative` = 423,085 peptides); `not_tested`
candidates are dropped (neither positive nor negative). Mutation-level immunogenic count = 213 CD8.

**Peptide-length distribution** (`mutant_seq`): 8mer 297,537 · 9mer 337,905 · 10mer 377,640 ·
11mer 417,204 · 12mer 357,424 — **all within MHCflurry's 8–15 supported range**.

**Join to `genotype_presentation_score` — FEASIBLE.** Every neo-peptide row carries both a peptide
(`mutant_seq`) and an HLA allele (`mutant_best_alleles`), 0 missing on either. So each assayed row is
directly scoreable by `Class1PresentationPredictor` (peptide + allele → `presentation_score` /
`presentation_percentile`), which is exactly the input the calibrator (#708) maps to an immunogenicity
log-odds. Per-patient genotype scoring can use the ≤6 alleles from `HLA_allotypes.txt`.

**HLA coverage (verified 2026-06-17).** `HLA_allotypes.txt` lists **150** patients vs **131** in the
data tables — but the 131 (mutation ∪ neo-pep) are a **strict subset** of the 150: every patient with
mutations/peptides has an HLA genotype (0 missing), so full-genotype scoring is safe for all of them.
The 19 surplus are HLA-only patients with no rows in either data table → drop silently on join.

### IMPROVE / Borch augment (`data/improve_borch/`, confirmed 2026-06-17)

Source: `github.com/SRHgroup/IMPROVE_paper` @ `c670942` (the repo named in the paper's **Data
Availability Statement** — Hadrup group, DTU; *not* `mnielLab`). Two TSVs copied; the larger
`data.zip`/`results.zip` bundles stay in the repo (reproducible via `git clone`).

| file | rows | cols | grain |
|---|---|---|---|
| `In_house_neoepitope_for_CV.tsv` | 17,520 | 3 | the training set — one row per candidate peptide |
| `Neoepitopes_CEDAR_benchmark_data.tsv` | 2,436 | 3 | external CEDAR held-out benchmark (548 positive) |

Schema (both): `Mut_peptide` · `HLA_allele` · `response`. **`response` is BINARY: 1 = immunogenic,
0 = not** (17,053 neg / 467 pos in the CV table — matches the paper's "17,520 / 467"). No `not_tested`
class, unlike NeoRanking. Lengths 8–11mers (8:515 · 9:10,561 · 10:4,716 · 11:1,728), all in range;
0 missing peptide/HLA → fully scoreable.

**Join nuance:** IMPROVE HLA is `HLA-A01:01` (no asterisk); NeoRanking is `A*01:01`. Both must be
normalized to MHCflurry's canonical `HLA-A*01:01` before scoring. Pooling the two cohorts also means
reconciling the 3-way (`CD8`/`negative`/`not_tested`) vs binary (`1`/`0`) label encodings — the
calibrator (#708) consumes assayed-positive vs assayed-negative, so map NeoRanking `CD8`→1,
`negative`→0, drop `not_tested`; IMPROVE is already in that form.

## Outputs index

- `data_manifest.yaml` — **pinned schema v1** (Issue #707): paths + checksums + fetch commands + license. Raw tables are **not** committed (LICR copyright, no redistribution).
- `outputs/calibrator_v1.joblib` — **deliverable artifact** (Issue #708): `PresentationCalibrator` fit with adaptive KDE mode, prior log-odds −6.525 (derived from `true_counts.csv` pooled across all 4 cohorts), trained on NCI + TESLA + HiTIDE + IMPROVE. Consumed by [Issue #709](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/709) (Snakemake wiring). Load: `PresentationCalibrator.load("outputs/calibrator_v1.joblib")`.
- `outputs/scored_cohort_subsample.parquet.true_counts.csv` — per-cohort true positive and true negative counts used to compute the base-rate prior (committed; the `.parquet` itself is gitignored as LICR-derived).
- `outputs/pr_reliability.png` — per-cohort **logit-scale calibration plot**: predicted log-odds (x) vs empirical logit of the observed rate (y, +0.5-corrected on the Kish-effective-n scale), with Wilson CIs, the 45° perfect-calibration diagonal, and the fitted **Cox recalibration** line. Headline numbers are the **Cox slope** (≈1 = well-calibrated spread; <1 = over-confident) and **intercept** (calibration-in-the-large / base-rate offset). Calibration is a held-out *diagnostic* on the SNV cohorts, not a deployment probability claim (#592); discrimination lives in the AUPRC-lift + top-k tables. The empirical curve is not expected to be monotone (sampling noise, not miscalibration); the monotonicity guarantee is structural (isotonic step, grid-checked at final fit).
- `outputs/shift_gap.png` — LOCO vs within-cohort validation gap: visualises the cross-lab/assay transfer drop (the "proxy caveat" gap).
- `outputs/kde_compare.png` — adaptive vs fixed-bandwidth KDE, logit-scale reliability for NCI + IMPROVE; confirms the two modes are indistinguishable on the large cohorts.
- `outputs/kde_fit_compare.png` — the **same pooled data fitted both ways**: (left) class-conditional densities under fixed (Scott) vs adaptive (Abramson) KDE; (right) the resulting calibration maps. Shows what "adaptive" changes (wider kernels in the data-sparse tail, smoother map) and why the ranking/AUPRC barely moves (both maps monotone, track closely).
- `outputs/method_fitted.png` — the **actual fitted pipeline** on our four pooled SNV cohorts (not a schematic): (1) class-conditional adaptive-KDE densities, (2) log density-ratio → base-rate prior shift, (3) isotonic → centered isotonic, (4) the fitted lookup curve `transform()` applies at inference. Replaces the earlier synthetic NeoGuider schematics; regenerated by the notebook's "Method, on our data" cell.
- `outputs/auprc_explained.png` — teaching figure: **how AUPRC is built**, on real LOCO predictions. Left: ranking → threshold sweep → precision/recall (TESLA candidates). Right: the precision–recall curves whose **area = AUPRC**, with the prevalence baseline, contrasting TESLA (signal, ~4× lift) vs IMPROVE (near-random, ~1× lift). Numbers reproduce the LOCO/lift tables exactly.
- `outputs/concept_presentation_funnel.png` — teaching schematic (illustrative, not data): presentation ⊃ immunogenic — why the presentation score's power is capped.
- `outputs/concept_discrimination_vs_calibration.png` — teaching schematic: the two evaluation axes — discrimination (separate the classes) vs calibration (are the numbers right).
- `outputs/panel_anatomy.png` — teaching schematic: how to read one calibration panel (45° = perfect, Cox-fit tilt = slope, offset = intercept) + four archetypes.
- `outputs/conclusion_scorecard.png` — conclusion visual: per-cohort verdict, discrimination (AUPRC lift) × calibration (Cox slope), bubble ∝ n positives; numbers pulled live from the result tables.
- `outputs/conclusion_ladder.png` — conclusion visual: the evidence ladder — within-cohort → LOCO (cross-SNV-cohort) → splice-junction (untested extrapolation, #680) — showing where the evidence runs out.
- `outputs/applicability/splice_applicability_support_map.png` — Issue #826 applicability gate (committed): left = splice vs SNV score distributions (overlap 0.92, KS 0.084); right = support map placing splice mass on the calibration curve (21% floor-clipped → flat constant log-odds). Reproduced by `applicability_notebook.ipynb`.

## Cohort-composition reconcile (primary-source confirmed)

[Issue #592](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/592) described NeoRanking as *"NCI-train + TESLA + HiTIDE + Bjerregaard, ~165 positives / 74 patients."* The **paper itself** (Müller 2023 Methods + Data S1, read from the Zotero PDF) states the cohort is:

| cohort | patients (mutation) | patients (neo-pep) | CD8 (mut / neo-pep) | source accession |
|---|---|---|---|---|
| NCI (Rosenberg/Gartner) | 112 | 80 | 147 / 103 | dbGaP `phs001003.v1.p1` |
| TESLA | 8 | 8 | 36 / 34 | Synapse `syn21048999` |
| HiTIDE (in-house) | 11 | 11 | 30 / 41 | EGA `EGAS00001007101` |
| **total** | **131** | **99** | **213 / 178** | — |

So #592 **diverges**: it's **131 patients, not 74**, and the cohorts are NCI + TESLA + HiTIDE — **no Bjerregaard** (that + the "74 pt / ~165 pos" framing is the **NeoGuider** 7-cohort benchmark bleeding in). Now **confirmed from the downloaded tables** (2026-06-17): mutation-level 131 patients / 213 immunogenic mutations; neo-peptide-level 99 patients (32 NCI patients have no testable 8–12mer) / **178 immunogenic neo-peptides** — matching the paper's "~178". A single consolidated correcting note was posted on #592 (AC #5).

## Download

### NeoRanking primary → `data/neoranking/` (browser route)

**⚠️ The NeoRanking README mislabels its figshare links** — it points the data-table names at the *classifier* records. Use the links from the **paper's** Data and Code Availability statement (Müller 2023, attachment `96V52NKW`) instead:

| data matrix (feature values + immunogenicity annotations) | real figshare link |
|---|---|
| neo-peptide → `Neopep_data_org.txt` | [147e67dde683fb769908](https://figshare.com/s/147e67dde683fb769908) |
| mutation → `Mutation_data_org.txt` | [2462b62bb6630fe2d257](https://figshare.com/s/2462b62bb6630fe2d257) |
| `HLA_allotypes.txt` | [35361871fdad4d1754d7](https://figshare.com/s/35361871fdad4d1754d7) |

### IMPROVE augment → `data/improve_borch/` (git clone)

Source from the **paper's** Data Availability Statement (Borch 2024, Frontiers full text): `github.com/SRHgroup/IMPROVE_paper` (Hadrup group, DTU) — **not** `mnielLab`.

```bash
git clone --depth 1 https://github.com/SRHgroup/IMPROVE_paper.git
cp IMPROVE_paper/neoepitope_tabels/{In_house_neoepitope_for_CV,Neoepitopes_CEDAR_benchmark_data}.tsv data/improve_borch/
```

### Verify + mirror (run from this folder)

```bash
sha256sum data/neoranking/*.txt data/improve_borch/*.tsv          # → match the manifest
gsutil -m cp data/neoranking/*.txt    gs://splice-neoepitope-project/experiments/issue_547/neoranking/
gsutil    cp data/improve_borch/*.tsv gs://splice-neoepitope-project/experiments/issue_547/improve_borch/
```

The `Classifier_neopeptide.zip` / `Classifiers_mutation.zip` from the NeoRanking README's `a000…`/`3c27…` links are the trained NCI-train LR/XGBoost **classifiers**, not data tables (kept out of `data/` here; fetch only if reproducing #708 figures).

## How to reproduce

Two separate Python environments are required because MHCflurry (TensorFlow/PyTorch) has no Python 3.14 wheels.

**Step 1 — Score cohorts (one-time; `mhcflurry-scoring` conda env)**

```bash
# From this folder. Use `conda activate` (NOT `conda run`, which buffers stdout and
# hides the per-allele progress logs during the 30+ min run — see CLAUDE.md).
conda activate mhcflurry-scoring && python score_cohort.py
# Writes: outputs/scored_cohort_subsample.parquet  (GITIGNORED)
#         outputs/scored_cohort_subsample.parquet.true_counts.csv
```

The `mhcflurry-scoring` env hosts MHCflurry 2.x + its model weights. The scored parquet is the clean handoff to the calibrator — no MHCflurry dependency past this point.

**Step 2 — Run tests (`research/.venv`, Python 3.14)**

```bash
cd research/experiments/issue_547_immunogenicity_calibration
../../../research/.venv/bin/python -m pytest tests/ -v
```

**Step 3 — Execute the notebook (`research/.venv`, Python 3.14)**

```bash
../../../research/.venv/bin/jupyter nbconvert --to notebook --execute --inplace notebook.ipynb
# Writes: outputs/calibrator_v1.joblib  outputs/*.png
```

**Step 4 — Splice applicability gate ([Issue #826](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/826); `research/.venv`, Python 3.14)**

```bash
# Fetch the splice presentation scores from the production patient runs (GITIGNORED dir).
mkdir -p data/splice_scores
for p in patient_001 patient_002; do
  gsutil cp "gs://splice-neoepitope-project/results/$p/predictions/mhc_presentation.tsv" \
            "data/splice_scores/${p}_mhc_presentation.tsv"
done
../../../research/.venv/bin/jupyter nbconvert --to notebook --execute --inplace applicability_notebook.ipynb
# Writes: outputs/applicability/splice_applicability_support_map.png
```

## Validation result

**LOCO (leave-one-cohort-out) prevalence-relative lift** — primary cross-lab/assay transfer metric:

| cohort left out | LOCO lift |
|---|---|
| NCI | ~111× |
| TESLA | ~4.4× |
| HiTIDE | ~2.7× |
| IMPROVE | ~1.2× (near baseline) |

NCI is the strongest discriminator; IMPROVE shows essentially no lift above the prevalence baseline (1×).

**Proxy caveat:** all four cohorts are **SNV point-mutation neoantigens**. LOCO tests cross-lab/assay generalization as a stand-in for cross-antigen-type generalization; LOCO success does **not** validate performance on splice-junction-derived neoantigens. That transfer gap is a known open question — see [Issue #547](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/547) and the design spec at [`docs/superpowers/specs/2026-06-18-issue-708-calibrator-design.md`](../../../docs/superpowers/specs/2026-06-18-issue-708-calibrator-design.md). The label-free **applicability** assessment of that transfer is below.

## Splice applicability — SNV→splice transfer ([Issue #826](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/826))

The evidence ladder's top rung — *splice = untested* — could not be turned into an *accuracy* validation, because no labelled splice-immunogenicity data exists ([Issue #680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680); TESLA excluded splice isoforms).
[Issue #826](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/826) instead runs the **label-free applicability gate** that *can* run now: does splice presentation-score mass fall where the calibrator actually learned a curve, or out in its flat extrapolation tails?
Notebook: [`applicability_notebook.ipynb`](applicability_notebook.ipynb); figure: [`outputs/applicability/splice_applicability_support_map.png`](outputs/applicability/splice_applicability_support_map.png).
Splice scores are the `genotype_presentation_score`s from the production patient runs (patient_001 n=395 + patient_002 n=1761 = **2156**), vs the **50,645** pooled SNV fit-cohort scores.

| Check | Metric | Value | Reading |
|---|---|---|---|
| Covariate-shift | overlap coefficient | **0.92** | splice scores overlap the SNV fit distribution almost entirely |
| Covariate-shift | KS statistic | **0.084** | small max-CDF gap; the two distributions track each other |
| Score-support | in-support mass `[0.0163, 0.9839]` | **78.4%** | interpolated — real fitted calibration |
| Score-support | floor-clipped mass (`< cx_[0] ≈ 0.0163`) | **21.4%** | flat constant −9.71 log-odds; near-non-presenters |
| Score-support | ceil-clipped mass (`> cx_[-1] ≈ 0.9839`) | **0.19%** | negligible |

**Verdict — provisional GO.** The SNV→splice transfer is *in-support*: there is no covariate-shift obstacle (overlap 0.92, KS 0.084 — splice scores live almost exactly where the calibrator was fit) and 78.4% of splice mass is interpolated.
The 21.4% floor-clip is **not** a covariate-shift artifact but an extrapolation-flatness property of the calibrator itself: the SNV fit-cohort *also* floor-clips **13.5%** of its own mass at the same `cx_[0] ≈ 0.0163`, so the flat low log-odds beyond the knot is intrinsic to the fitted curve, not a splice-specific gap. (Splice floor-clips modestly more — 21.4% vs 13.5% — reflecting somewhat more near-non-presenters, but well within the populated SNV range.)
The floor-clip lands on peptides with `genotype_presentation_score < 0.0163` — near-non-presenters that a flat low log-odds scores *correctly* (lowest presenters → lowest immunogenicity), so it is benign for ranking rather than misleading.
The two checks are in fact the same difference seen twice: the KS statistic (**0.084**) is essentially the floor-region CDF gap (splice CDF at the floor knot ≈ 0.214 vs SNV ≈ 0.135, gap ≈ **0.079** ≈ KS), so the *only* material covariate difference is the ~8pp extra floor mass — and it sits **exactly in the flat-extrapolation region** that `out_of_calibration_support` already flags, where the calibrated log-odds is constant regardless. The shift therefore lands in the most benign place it could, making the transfer benign for ranking *by construction*, not merely empirically.
This is a **label-free applicability gate, not an accuracy validation** — it certifies the transfer is in-support, *not* that the calibrated values are correct for splice neoepitopes; accuracy on measured splice labels remains [Issue #870](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/870), under the open benchmark [Issue #680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680).

> **Reference note.** The SNV comparand is `outputs/scored_cohort_subsample.parquet` (n=50,645) — the exact data `calibrator_v1` was fit on, so it is the correct "where did the calibrator see data" reference by construction. It is label-stratified (all 645 positives + 50,000 cohort-stratified negatives), but positives are only **1.27%** of rows, so the score marginal is ~98.7% set by the negative subsample and the enrichment does not meaningfully skew the overlap/KS comparison. The KS p-value (4.86e-13) is a large-`n` artifact — at n=50,645 the test rejects on any trivial difference — so the KS *statistic* (0.084) and overlap (0.92) are the load-bearing readings, not the p-value.

**Conditions routed to [Issue #709](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/709)** (production wiring):
1. The calibrated column is a **provisional SECONDARY ranking signal** — `genotype_presentation_score` stays the **primary** ranker; the calibrated immunogenicity log-odds augments, never replaces it. (Load-bearing guardrail: splice immunogenicity is unvalidated, so it must not drive primary ranking.)
2. Rows whose `genotype_presentation_score` falls **outside** the interpolation support `[cx_[0], cx_[-1]] ≈ [0.0163, 0.9839]` carry an **`out_of_calibration_support`** flag — covering both the 21.4% floor-clip **and** the 0.19% ceil-clip, where `transform()` returns a flat constant — so the flat-extrapolation region stays transparent downstream and is never mistaken for a fitted value.
3. The **provisional** label is discharged by [Issue #870](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/870) (validate the SNV→splice transfer on measured splice labels), not by a general [Issue #680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680) pointer.

## Cross-experiment deps

- Consumes the cohort selection decision from [Issue #592](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/592) (the cohort-calibration gate deck).
- Produces the cohort tables consumed by [Issue #708](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/708) (calibrator) and [Issue #709](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/709) (Snakemake wiring).
