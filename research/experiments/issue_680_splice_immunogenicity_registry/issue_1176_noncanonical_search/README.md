# Issue #1176 - non-canonical public-spectra re-search

Analysis set for [Issue #1176](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1176), the leg of best-bet 3 that re-searches public tumor immunopeptidome spectra against a junction-derived non-canonical database.
The design/scoping front-half lives one level up in [`../deep_research/run4_issue1176_noncanonical_search_scoping_2026-07-15.md`](../deep_research/run4_issue1176_noncanonical_search_scoping_2026-07-15.md); this directory holds the back-half's decision records and their supporting measurements.

## Contents

- **[`decision_normal_filter_2026-07-23.md`](decision_normal_filter_2026-07-23.md)** - the normal-filter decision.
  The Courcelles cohort has no matched normal, so the GTEx pan-tissue filter becomes the sole specificity control.
  The decision is to proceed with that substitution and report it as a limitation; the brief carries the measured cost, the limitations, and the re-open triggers.
- **`normal_vs_gtex_control.py`** - the matched-pair control behind that decision.
- **`outputs/`** - its results (`normal_vs_gtex.json` plus the per-run classified and stats TSVs).

## Reproducing the control

The control runs the pipeline's own `classify_junctions` twice over the chr22 gastric tumor/normal pair, flipping exactly one variable (the matched-normal arm) and holding everything else identical.
It asserts that both funnels reconcile against `junctions_raw`, so a silent accounting error fails the run rather than producing a plausible number.

It needs the GTEx chr22 blacklist, which is gitignored and regenerable.
**Use `--region`, not `--restrict-chrom` alone** - the latter filters rows after download and leaves the query genome-wide:

```bash
conda activate snakemake
python workflow/scripts/build_gtex_pan_tissue_ref.py \
  --output-bed resources/test/gtex_gtexv2_pan_tissue_junctions.chr22.bed \
  --output-qc  resources/test/gtex_gtexv2_pan_tissue_junctions.chr22.qc.json \
  --min-samples 1 --region chr22 --restrict-chrom chr22

python research/experiments/issue_680_splice_immunogenicity_registry/issue_1176_noncanonical_search/normal_vs_gtex_control.py
```

The rebuild is verifiable against `resources/test/gtex_gtexv2_data_manifest.yaml`: 880769 junctions, 47561526 bytes, sha256 `d727ce5c8489f5940f69df229c3c9acfcb8d8ff071a6ac61131510b760da7f58`.
Note that the manifest's `gcs_path` is stale (that bucket was deleted in the June GCP decommission); the artifact is recoverable by rebuild, not by fetch.

## Result in one line

Of the 4 junctions a matched normal removes, GTEx independently blacklists 3; the candidate set inflates by 1 junction (0.8%).
Read the decision brief for why that is reassuring rather than conclusive, and for the limitations that keep it from being a calibrated estimate.
