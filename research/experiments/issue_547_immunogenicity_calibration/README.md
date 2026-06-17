# issue_547 — immunogenicity calibration

**Goal:** build the post-MHCflurry **immunogenicity calibration** layer (KDE + centered isotonic regression) that maps `presentation_score` → `calibrated_immunogenicity_log_odds`. This folder is the **epic-level home** for that work; the first sub-issue ([Issue #707](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/707)) acquires the labelled training cohorts chosen in [Issue #592](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/592) and confirms their schema is usable. Nothing downstream (#708 calibrator, #709 Snakemake wiring) starts until the tables are in hand and parseable.

**Parent epic:** [Issue #547](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/547) · **acquisition sub-issue:** [Issue #707](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/707)

Folder is **epic-keyed** (`issue_547`, not `issue_707`) because the cohort data + calibrator are shared across the sub-issues; it mirrors the GCS prefix `gs://splice-neoepitope-project/experiments/issue_547/`.

## Layout

```
issue_547_immunogenicity_calibration/
├── README.md            # this file (committed)
├── data_manifest.yaml   # pinned schema v1 — paths + checksums + fetch cmds + license (committed)
├── data/                # raw cohort tables — GITIGNORED (root .gitignore `data/`); not on fresh clones
└── (later) notebook.ipynb + outputs/   # the #708 calibrator analysis
```

Raw tables are **never committed**: NeoRanking is LICR copyright (no redistribution) and the tables are large. `data/` is gitignored repo-wide; recreate it with `mkdir -p data/` after clone and fetch per `data_manifest.yaml`.

## Status (2026-06-17)

| Step | State |
|------|-------|
| Sources located + recorded | ✅ figshare share links + filenames + LICR license (`data_manifest.yaml`) |
| Byte-download (browser route) | ✅ 3 NeoRanking tables downloaded + unzipped into `data/` |
| Checksums + sizes | ✅ sha256 + byte sizes in `data_manifest.yaml` |
| Schema + join-feasibility doc | ✅ see "Schema & join feasibility" below |
| GCS mirror (>100 MB) | ✅ `gs://splice-neoepitope-project/experiments/issue_547/` |
| Cohort-composition reconcile vs #592 | ✅ paper+data confirmed: mutation 131 pt / neo-pep 99 pt, 178 CD8 neo-peptides, **no Bjerregaard** — #592 correcting note drafted (pending post) |
| IMPROVE/Borch augment | 📄 access path documented (not blocking; byte-download deferred to VM) |

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

## Outputs index

- `data_manifest.yaml` — **pinned schema v1** (Issue #707): paths + checksums + fetch commands + license. Raw tables are **not** committed (LICR copyright, no redistribution).

## Cohort-composition reconcile (primary-source confirmed)

[Issue #592](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/592) described NeoRanking as *"NCI-train + TESLA + HiTIDE + Bjerregaard, ~165 positives / 74 patients."* The **paper itself** (Müller 2023 Methods + Data S1, read from the Zotero PDF) states the cohort is:

| cohort | patients (mutation) | patients (neo-pep) | CD8 (mut / neo-pep) | source accession |
|---|---|---|---|---|
| NCI (Rosenberg/Gartner) | 112 | 80 | 147 / 103 | dbGaP `phs001003.v1.p1` |
| TESLA | 8 | 8 | 36 / 34 | Synapse `syn21048999` |
| HiTIDE (in-house) | 11 | 11 | 30 / 41 | EGA `EGAS00001007101` |
| **total** | **131** | **99** | **213 / 178** | — |

So #592 **diverges**: it's **131 patients, not 74**, and the cohorts are NCI + TESLA + HiTIDE — **no Bjerregaard** (that + the "74 pt / ~165 pos" framing is the **NeoGuider** 7-cohort benchmark bleeding in). Now **confirmed from the downloaded tables** (2026-06-17): mutation-level 131 patients / 213 immunogenic mutations; neo-peptide-level 99 patients (32 NCI patients have no testable 8–12mer) / **178 immunogenic neo-peptides** — matching the paper's "~178". A single consolidated correcting note for #592 is drafted (AC #5; pending post).

## Download (browser route — into `data/`)

**⚠️ The NeoRanking README mislabels its figshare links** — it points the data-table names at the *classifier* records. Use the links from the **paper's** Data and Code Availability statement (Müller 2023, attachment `96V52NKW`) instead:

| data matrix (feature values + immunogenicity annotations) | real figshare link |
|---|---|
| neo-peptide → `Neopep_data_org.txt` | [147e67dde683fb769908](https://figshare.com/s/147e67dde683fb769908) |
| mutation → `Mutation_data_org.txt` | [2462b62bb6630fe2d257](https://figshare.com/s/2462b62bb6630fe2d257) |
| `HLA_allotypes.txt` | [35361871fdad4d1754d7](https://figshare.com/s/35361871fdad4d1754d7) |

Then (run from this folder):

```bash
sha256sum data/Neopep_data_org.txt data/Mutation_data_org.txt data/HLA_allotypes.txt   # → fill PENDING in manifest
gsutil cp data/*.txt gs://splice-neoepitope-project/experiments/issue_547/
```

The `Classifier_neopeptide.zip` / `Classifiers_mutation.zip` already in `data/` came from the README's `a000…`/`3c27…` links — those are the trained NCI-train LR/XGBoost **classifiers**, not data tables (kept for potential #708 figure-repro).

## Cross-experiment deps

- Consumes the cohort selection decision from [Issue #592](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/592) (the cohort-calibration gate deck).
- Produces the cohort tables consumed by [Issue #708](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/708) (calibrator) and [Issue #709](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/709) (Snakemake wiring).
