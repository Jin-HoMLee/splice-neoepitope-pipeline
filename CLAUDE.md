# Splice Neoepitope Pipeline — Project Notes

## Instructions for Claude
Update this file only when something is **not derivable from the code or git history** — non-obvious decisions, known gotchas, workarounds, and infrastructure facts. Do not use it as a change log; that's what `git log` is for.

## Project Overview
Modernised reimplementation of a 2015 cancer neoepitope prediction pipeline (Jin-Ho Lee, Seoul National University). Identifies tumor-specific splice junctions from RNA-Seq data and predicts MHC-binding neoepitopes.

## Infrastructure
- Running on a GCP Compute Engine VM (`splice-pipeline`, `us-central1-a`)
- See `docs/google_cloud_guide.md` for full setup instructions
- Pipeline is run with `snakemake --cores $(nproc) --use-conda` inside a `tmux` session

## Pipeline Design Decisions

### Junction origin classification
Normal samples are used to filter tumor junctions at the junction level, not at the prediction level. The hierarchy:
```
all junctions
  └─ annotated        (in GENCODE)            → discard
  └─ unannotated      (not in GENCODE)
       ├─ patient_specific  (also in normal)  → discard (kept in TSV for reference)
       └─ tumor_specific    (absent in normal) → neoepitope prediction
```
This is the clinically correct approach: a junction present in matched normal tissue is not tumor-specific and should not be a neoepitope target. The Fisher's exact test (end-of-pipeline statistical comparison) was removed in favour of this upstream filtering step.

When no normal sample is present, all unannotated junctions are labeled `tumor_specific` with a warning — the pipeline still runs.

## Known Dependency Issues (Fixed)

### `hisat2.yaml` — samtools/libdeflate conflict
`regtools >= 1.0.0` and `hisat2` pull in `libdeflate >= 1.26`; older `samtools`/`htslib` builds required `libdeflate < 1.26`.
**Fix:** pin `samtools >= 1.20` in `hisat2.yaml`. Earlier workaround (removing samtools and using base env PATH) broke on macOS where the base env is not on PATH inside activated conda envs.

### `run_mhcflurry.py` — mhcflurry 2.2.0 API change
mhcflurry 2.2.0 changed `predict()` to return a raw numpy array instead of a DataFrame.
**Fix:** use `predict_to_dataframe()` instead.

## auto_stop.sh
Shuts down the VM after the pipeline finishes to save costs.
Uses `sudo shutdown -h now` — **not** `gcloud compute instances stop`, which fails with
`ACCESS_TOKEN_SCOPE_INSUFFICIENT` because the VM service account lacks the compute API scope.

Run as:
```bash
snakemake --cores $(nproc) --use-conda --rerun-triggers mtime 2>&1 | tee pipeline.log ; bash auto_stop.sh
```

## sra-tools Note
Use version `3.1.1` on GCP VMs — newer versions (3.4.x) have a segfault bug.
On macOS arm64, sra-tools conda installation is unreliable due to libcurl/openssl conflicts. Use ENA HTTPS download instead (see `scripts/prepare_test_data.sh`).

## regtools Argument Order
`regtools junctions extract` requires all options before the positional BAM argument. Placing `-o` after the BAM causes `Error parsing inputs!(2)`.
```bash
regtools junctions extract -s XS -a 8 -m 50 -M 500000 -o out.bed input.bam
```

## Local Test Dataset (chr22)
For development and testing on macOS (M1, 8 GB RAM):
```bash
bash scripts/prepare_test_data.sh   # one-time: downloads reference + FASTQs
snakemake --cores 4 --use-conda --configfile config/test_config.yaml
```
- Reference: chr22 FASTA (UCSC hg38) + GENCODE v47 GTF filtered to chr22
- FASTQs: 500K reads each from a matched gastric cancer pair via ENA HTTPS (no sra-tools):
  - **SRR9143066** — Primary Tumor (gastric cancer surgical section)
  - **SRR9143065** — Solid Tissue Normal (adjacent stomach tissue)
- Both samples are single-end Illumina HiSeq 3000; HISAT2 handles this via `-U` mode
- HISAT2 index stored in `resources/test/hisat2_index/` (separate from production)
- All test outputs go to `results/test/` and `logs/test/`
