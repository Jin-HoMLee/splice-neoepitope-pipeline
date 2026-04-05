# Splice Neoepitope Pipeline — Project Notes

## Instructions for Claude
Keep this file up to date. Whenever changes are made during a session — bug fixes, dependency updates, new pipeline steps, config changes, decisions made — append a summary to the relevant section or add a new section. Update this file before ending the conversation.

## Project Overview
Modernised reimplementation of a 2015 cancer neoepitope prediction pipeline (Jin-Ho Lee, Seoul National University). Identifies novel splice junctions from RNA-Seq data and predicts MHC-binding neoepitopes.

Active branch: `copilot/modernize-cancer-neoepitope-pipeline`

## Infrastructure
- Running on a GCP Compute Engine VM (`splice-pipeline`, `us-central1-a`)
- See `docs/google_cloud_guide.md` for full setup instructions
- Pipeline is run with `snakemake --cores $(nproc) --use-conda` inside a `tmux` session

## History
- Completed initial GCP setup and debugging
- Completed a first pipeline run using sample `SRR37781424` (Luminal A breast cancer tumor)
- Several runs failed due to conda dependency issues (see below)

## Known Dependency Issues (Fixed)

### 1. `hisat2.yaml` — samtools/libdeflate conflict
`regtools >= 1.0.0` and `hisat2` pull in `libdeflate >= 1.26`; all `samtools`/`htslib` versions require `libdeflate < 1.26` — irreconcilable in the same env.
**Fix:** removed `samtools` from `hisat2.yaml` entirely. `samtools` from the base conda env (`~/miniforge3/bin/samtools 1.23.1`) is used via PATH instead. Same issue would affect `biotools.yaml` — `samtools` and `pysam` were also removed there (neither is imported by any script in that env).

### 2. `biotools.yaml` — missing pandas
`workflow/scripts/assemble_contigs.py` imports `pandas` but it was not listed in `workflow/envs/biotools.yaml`.
**Fix:** added `pandas >= 2.0` to `workflow/envs/biotools.yaml`.

### 3. `run_mhcflurry.py` — mhcflurry 2.2.0 API change
mhcflurry 2.2.0 changed `predict()` to return a raw numpy array of affinities instead of a DataFrame. Column names also changed: `mhcflurry_affinity` → `prediction`, `mhcflurry_affinity_percentile` → `prediction_percentile`.
**Fix:** switched to `predict_to_dataframe()` (returns a proper DataFrame with percentile ranks) and updated column name references in `workflow/scripts/run_mhcflurry.py`.

## sra-tools Note
Use version `3.1.1` — newer versions (3.4.x) have a segfault bug on GCP VMs.
```
conda install -c bioconda sra-tools=3.1.1 -y
```
