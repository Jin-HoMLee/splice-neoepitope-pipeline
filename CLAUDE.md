# Splice Neoepitope Pipeline — Project Notes

## Instructions for Claude
Keep this file up to date. Whenever changes are made during a session — bug fixes, dependency updates, new pipeline steps, config changes, decisions made — append a summary to the relevant section or add a new section. Update this file before ending the conversation.

## Project Overview
Modernised reimplementation of a 2015 cancer neoepitope prediction pipeline (Jin-Ho Lee, Seoul National University). Identifies novel splice junctions from RNA-Seq data and predicts MHC-binding neoepitopes.

Active branch: `main` (feature work done on issue branches, e.g. `3-add-chr22-test-dataset-for-local-macos-runs`)

## Infrastructure
- Running on a GCP Compute Engine VM (`splice-pipeline`, `us-central1-a`)
- See `docs/google_cloud_guide.md` for full setup instructions
- Pipeline is run with `snakemake --cores $(nproc) --use-conda` inside a `tmux` session

## History
- Completed initial GCP setup and debugging
- Completed a first pipeline run using sample `SRR37781424` (Luminal A breast cancer tumor)
- Several runs failed due to conda dependency issues (see below)
- Full end-to-end cloud run succeeded (2026-04-07)
- Set up chr22 test dataset for local macOS development (see below)

## Known Dependency Issues (Fixed)

### 1. `hisat2.yaml` — samtools/libdeflate conflict (resolved)
`regtools >= 1.0.0` and `hisat2` pull in `libdeflate >= 1.26`; older `samtools`/`htslib` builds required `libdeflate < 1.26`.
**Fix:** pin `samtools >= 1.20` in `hisat2.yaml`. samtools 1.20+ updated its libdeflate dependency and is compatible with the newer libdeflate that regtools requires. Earlier workaround (removing samtools and using the base env PATH) was dropped as it broke on macOS where the base env is not on PATH inside activated conda envs. `samtools` and `pysam` remain removed from `biotools.yaml` — neither is imported by any script in that env.

### 2. `biotools.yaml` — missing pandas
`workflow/scripts/assemble_contigs.py` imports `pandas` but it was not listed in `workflow/envs/biotools.yaml`.
**Fix:** added `pandas >= 2.0` to `workflow/envs/biotools.yaml`.

### 3. `run_mhcflurry.py` — mhcflurry 2.2.0 API change
mhcflurry 2.2.0 changed `predict()` to return a raw numpy array of affinities instead of a DataFrame. Column names also changed: `mhcflurry_affinity` → `prediction`, `mhcflurry_affinity_percentile` → `prediction_percentile`.
**Fix:** switched to `predict_to_dataframe()` (returns a proper DataFrame with percentile ranks) and updated column name references in `workflow/scripts/run_mhcflurry.py`.

### 4. `statistical_analysis.py` — sample_type NaN in local mode
In the `else` branch (local mode / no normal samples), the script merged `source_header` against `manifest["file_id"]` — these never match because `source_header` is a full `junc_id|coords|sample_type|frame` string. `sample_type` stayed NaN (float), crashing `generate_report.py` on `.str.contains()`.
**Fix:** parse `sample_type` directly from the 3rd pipe-separated field of `source_header` in `workflow/scripts/statistical_analysis.py`.

## auto_stop.sh
Shuts down the VM after the pipeline finishes (success or error) to save costs.
Uses `sudo shutdown -h now` — **not** `gcloud compute instances stop`, which fails with
`ACCESS_TOKEN_SCOPE_INSUFFICIENT` because the VM service account lacks the compute API scope.
Guest-initiated OS shutdown stops the instance without triggering automatic restart.

Run as:
```bash
snakemake --cores $(nproc) --use-conda --rerun-triggers mtime 2>&1 | tee pipeline.log ; bash auto_stop.sh
```

## sra-tools Note
Use version `3.1.1` on GCP VMs — newer versions (3.4.x) have a segfault bug.
On macOS arm64, sra-tools conda installation is unreliable due to libcurl/openssl
conflicts. Use ENA HTTPS download instead (see `scripts/prepare_test_data.sh`).

## regtools Argument Order
`regtools junctions extract` requires all options before the positional BAM
argument. Placing `-o` after the BAM causes `Error parsing inputs!(2)`.
Correct order:
```bash
regtools junctions extract -s XS -a 8 -m 50 -M 500000 -o out.bed input.bam
```

## Local Test Dataset (chr22)
For development and testing on macOS (M1, 8 GB RAM), use the chr22 subset:
```bash
bash scripts/prepare_test_data.sh   # one-time: downloads reference + FASTQs
snakemake --cores 4 --use-conda --configfile config/test_config.yaml
```
- Reference: chr22 FASTA (UCSC hg38) + GENCODE v47 GTF filtered to chr22
- FASTQ: 500K read pairs from ERR188273 (GEUVADIS LCL, ENA HTTPS — no sra-tools)
- Runtime: ~2 min; produces ~400 novel junctions, ~80 strong MHC binders
- HISAT2 index stored in `resources/test/hisat2_index/` (separate from production)
- All test outputs go to `results/test/` and `logs/test/`
