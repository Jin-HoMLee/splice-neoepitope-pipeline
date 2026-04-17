# Lab Notebook

---

## 2026-04-17

### Patient_001 (gastric cancer) — first full production run with HLA typing + TCRdock

**Patient:** SRR9143065 (Solid Tissue Normal) / SRR9143066 (Primary Tumor) — gastric cancer surgical section, single-end Illumina HiSeq 3000.

**Junction filtering results:**
- Unannotated junctions: 30,029
- Normal-shared (filtered out): 2,682 (8.9%)
- Tumor-exclusive candidates: 27,347

**HLA typing (OptiType):** HLA-A\*26:01, HLA-A\*31:01, HLA-B\*15:05, HLA-B\*18:20, HLA-C\*03:21, HLA-C\*07:01 — no ground-truth alleles available for this patient, so typing cannot be validated.

**MHCflurry predictions:** 8,226 peptide × allele pairs (1,371 9-mers across 6 alleles), 54 strong binders (IC50 ≤ 50 nM), 317 weak binders (IC50 ≤ 500 nM).

**Top neoepitope candidate:** EVAEYNASF / HLA-A\*26:01, IC50 16.5 nM (strong binder). TCRdock structure predicted on P100 GPU VM; Mol\* viewer renders ternary complex with correct chain labels (A=MHC, B=peptide, C=TCR-α, D=TCR-β). Results archived to `gs://splice-neoepitope-project-tcrdock-handoff/results/`.

**Infrastructure issues resolved during this run:**
- P100 GPU VM required proprietary nvidia kernel modules (`linux-modules-nvidia-570-server-6.8.0-1053-gcp`); open-source modules (`linux-modules-nvidia-570-server-open-*`) do not support Pascal (P100) GPUs.
- GCS download of pipeline outputs onto the GPU VM fresh-stamped all files, causing Snakemake's `--rerun-triggers mtime` to re-run MHCflurry unnecessarily. Fixed by switching the GPU phase to `--rerun-triggers code params`.
- Mol\* viewer silently broken by unpinned CDN URL (`unpkg.com/molstar` → v5.8.0, breaking API). Pinned to `molstar@4.9.0`.

---

### Issue #50 — parallel MHCflurry allele predictions (ProcessPoolExecutor)

**Problem:** Serial per-allele loop was the bottleneck for patients with 6+ HLA alleles.

**First attempt (ThreadPoolExecutor) — failed:** `predict_to_dataframe()` is not thread-safe due to shared TensorFlow state. Testing revealed identical IC50 values (267.12 nM) across different alleles for the same peptide — confirmed state corruption. Threads release the GIL during I/O but TensorFlow's internal state is mutated during inference.

**Fix:** Switched to `ProcessPoolExecutor` with `initializer=_worker_init`. Each worker process loads its own predictor copy, sets `TF_NUM_INTRAOP_THREADS=1` / `OMP_NUM_THREADS=1` before TensorFlow imports (prevents CPU oversubscription when running multiple workers), and returns a lean per-allele DataFrame (no `peptides_df` pickling across process boundaries).

**Test result:** 8,226 rows, 54 strong, 317 weak — max IC50 diff vs. serial baseline = 0.00e+00. Local run: 6 alleles, 4 workers, ~12 s total.

**Also renamed:** `predict.smk` → `mhcflurry.smk`, `predictions.tsv` → `mhc_affinity.tsv` (tool-agnostic naming).

---

### Patient_002 planning — osteosarcoma IPISRC044

**Dataset:** Publicly available osteosarcoma dataset (https://osteosarc.com/data/). Patient IPISRC044, multi-institutional (UCLA / UCSF / Boston Gene / Tempus). GCS bucket: `gs://osteosarc-genomics`.

**Plan:** Start with T0 tumor (Boston Gene, Nov 2022) as the baseline timepoint.
- FASTQs: `rna-seq/fastq/bostongene_2022/202211_bostongene_tumor_rna_BG003082_R1.fastq.gz` + `_R2.fastq.gz` (paired-end, ~10 GB)
- Run OptiType for HLA typing — ground-truth Class I alleles are known from Red Cross serology (A\*01:01/01:11N, B\*08:01/27:05, C\*01:02/07:01), giving us a validation opportunity we didn't have for patient_001.

**No matched RNA-seq normal available.** Blood WGS DNA cannot be used as a substitute — `regtools junctions extract` requires spliced RNA-seq reads (`N` CIGAR operations); DNA-seq reads map continuously without splicing and produce no junctions. Pipeline runs in warning mode, labelling all unannotated junctions `tumor_exclusive`. Based on patient_001 statistics, approximately 9% of unannotated junctions may be normal-shared and would be misclassified. The downstream MHCflurry + TCRdock funnel is expected to absorb most of this noise.

---

## 2026-04-16

**Goal:** Implement HLA typing with OptiType and stabilise the cloud pipeline for a production run.

**Done:**
- Implemented HLA typing step (issue #42): OptiType runs on each sample's FASTQ, typing results aggregated per patient into `alleles.tsv`, passed to MHCflurry in place of fallback alleles when `hla.enabled: true`.
- Added configurable CBC ILP solver for OptiType (issue #49) — reduces OptiType runtime significantly on VMs with multiple cores.
- Fixed `PYTHONUNBUFFERED=1` for OptiType log flushing (issue #52): without it, log output was buffered and appeared to stall.
- Fixed HLA concordance display in report (issue #53): loci with no normal/tumor discrepancy now show ✓ concordant instead of blank.
- Upgraded CPU VM to n1-standard-8 and made thread counts config-driven (issue #40).
- Fixed Snakemake unlock on interrupted restarts (issue #47).
- Suppressed interactive "Next steps" prompts when `run_cloud_gpu.sh` calls sub-scripts (issue #48).
- Dropped all GDC/TCGA code (issue #44) — pipeline is now fully open-access, no registration-gated data sources.
- Renamed junction origin labels to `tumor_exclusive` / `normal_shared` (issue #38).

---

## 2026-04-14

**Goal:** Refactor pipeline to be patient-centric and unify alignment rules.

**Done:**
- Replaced `{cancer_type}` wildcard with `{patient_id}` throughout (issue #26). Config key `cancer_types: [local]` → `patient_id` string. Routing in `filter.smk` now switches on `config["data_source"]` rather than hardcoded `"local"` path segments.
- Unified STAR and HISAT2 alignment into a single `alignment.smk` module (issue #35). Both aligners produce junction TSVs in the same format; downstream rules are aligner-agnostic.
- `PATIENT_IDS` now derived from `samples_tsv` at DAG construction time; sample IDs renamed to SRR accessions throughout.

---

## 2026-04-13

**Goal:** Merge TCRdock structural validation and get a clean production baseline.

**Done:**
- Merged PR #27 (issue #25): TCRdock step, Mol\* 3D viewer, `run_cloud_gpu.sh` three-phase lifecycle (CPU → GCS → GPU Spot VM). Full end-to-end test passed.

---

## 2026-04-10

**Goal:** Get TCRdock structural validation running end-to-end via Docker on the GCP GPU Spot VM.

**Done:**
- Fixed Docker build failure: the official TCRdock Dockerfile fails at `conda install openmm=7.7.0` (no longer available). Created `docker/Dockerfile.pipeline` which uses pip only and omits openmm/pdbfixer — these are only needed for optional Amber relaxation, not for `run_prediction.py`.
- Successfully built `tcrdock:latest` on GPU VM and ran TCRdock end-to-end. Pipeline completed 3/3 steps, report generated.
- Identified visual issue in report: AlphaFold outputs all residues as a single chain A, so Mol* rendered the ternary complex without chain-colour distinction. Added `relabel_pdb_chains()` to `run_tcrdock.py`, which reassigns chain IDs (A=MHC, B=peptide, C=TCR-alpha, D=TCR-beta) using the chain lengths from TCRdock's `alphafold_setup/targets.tsv`.
- Updated `setup_tcrdock_vm.sh` to build from `docker/Dockerfile.pipeline` instead of cloning TCRdock separately for the (broken) official Dockerfile.

### 2026-04-10 15:30 — Mol* COMPND fix

Fixed Mol* sequence panel showing generic "Polymer 1/2/3/4" instead of meaningful chain names. Root cause: PDB COMPND records had off-by-one column positions (continuation number in cols 9–11 instead of PDB-standard 8–10) and the first line incorrectly included a continuation number. Also padded lines to 80 chars and transliterated Unicode α/β to ASCII for PDB compliance. Confirmed fix locally — Mol* now shows "MHC heavy chain", "Peptide", "TCR alpha", "TCR beta".

### 2026-04-10 17:00 — Pre-PR refactoring

Code cleanup before creating PR for #25:
- `generate_report.py`: fixed `html` module shadowing (local variable named `html` overrode the stdlib import → renamed to `report_html`, import aliased to `html_mod`). Moved `import json` from function body to top-level. Extracted COMPND record building into `_build_compnd_records()` helper.
- `run_tcrdock.py`: replaced `assert` with `raise ValueError` for input validation.
- All 57 tests passing.

### 2026-04-10 18:00 — Documentation updates

Updated all project documentation for the branch:
- README.md: pipeline diagram (now 7 steps), TOC, config table with `tcrdock.*` parameters, output tree with `tcrdock/` directory, project structure with `docker/` and new files, citations for TCRdock and Mol*.
- `docs/google_cloud_guide.md`: new "Automated GPU Pipeline (TCRdock)" section with quick start, retrieval, how-it-works, detached mode, cost estimate. Reordered TOC so Troubleshooting is last.
- CLAUDE.md: added "TCRdock via Docker" and "PDB chain relabelling" decision notes.

### 2026-04-10 19:00 — Final cloud test

Kicked off `run_cloud_gpu.sh` for end-to-end validation on GCP. Pending result before creating PR for #25.

---

## 2026-04-09

**Goal:** Implement TCRdock structural validation (issue #25) and automate the full CPU→GPU cloud pipeline.

**Done:**
- Merged issues #20 (filter junction-spanning 9-mers) and #22 (move junction-spanning filter to translation step; output peptides as TSV). Closed #21.
- Implemented TCRdock step (`workflow/rules/tcrdock.smk`, `workflow/scripts/run_tcrdock.py`) and Mol* 3D viewer in the HTML report. Verified TCRdock input column format and two-step workflow (`setup_for_alphafold.py` → `run_prediction.py`) against the real API before writing code.
- Wrote `scripts/run_cloud_gpu.sh`: three-phase lifecycle (CPU VM steps 1–5 → GCS handoff → GPU Spot VM TCRdock). Used GCS bucket (`tcrdock-handoff`) for VM-to-VM transfer.
- Hit a series of dependency issues with the conda-based TCRdock env: wrong script name (`predict.py` vs `run_prediction.py`), missing tensorflow, Python 3.8 incompatibility with dm-haiku, cuDNN 8.9 vs 9.10 mismatch on the Deep Learning VM image.
- Decided to switch to Docker to sidestep the dependency issues entirely. Rewrote `run_tcrdock.py` to call TCRdock via `docker run --gpus all` instead of a conda env; updated `tcrdock_gpu.yaml` and `tcrdock.smk` accordingly.
- Official TCRdock Dockerfile fails at `conda install openmm=7.7.0`. Confirmed openmm/pdbfixer are only used in `alphafold/relax/` (Amber relaxation), not in `run_prediction.py`. Session ended with this as the open problem.

**Key decisions:**
- Use Docker for TCRdock rather than a conda env — eliminates host-side cuDNN/JAX version management.
- Use `--new_docking` flag (1 AlphaFold run per target instead of 3) to reduce GPU time.
- GCS bucket for VM-to-VM result transfer rather than direct SCP, so the two VMs don't need to be up simultaneously.

---

## 2026-04-07

**Goal:** Run full pipeline on real cancer data; add local test dataset for macOS development.

**Done:**
- Full cloud pipeline run succeeded on SRR37781424 (Luminal A breast cancer, GEO GSE119889). All steps 1–5 completed on `splice-prod-test` VM.
- Added `scripts/prepare_test_data.sh` and chr22 test config for local macOS runs (M1, 8 GB RAM). Downloads chr22 reference + 500K-read subsets of a matched gastric cancer pair (SRR9143066 tumor / SRR9143065 normal) via ENA HTTPS to avoid sra-tools arm64 issues.
- Merged PR #1 (Copilot-assisted modernisation baseline).

---

## 2026-04-05

**Goal:** Fix pipeline bugs found during first cloud run.

**Done:**
- Fixed `auto_stop.sh`: `gcloud compute instances stop` fails with `ACCESS_TOKEN_SCOPE_INSUFFICIENT` from inside the VM. Switched to `sudo shutdown -h now`.
- Fixed conda env dependency conflicts (samtools/libdeflate, regtools version pinning).
- Fixed mhcflurry 2.2.0 API change: `predict()` now returns a numpy array, not a DataFrame. Switched to `predict_to_dataframe()`.
- Fixed `statistical_analysis.py`: extract `sample_type` from `source_header` correctly.
- Auto-download MHCflurry models as a pipeline step rather than requiring manual setup.

---

## 2026-04-04

**Goal:** Document GCP deployment so the pipeline can be handed off or reproduced.

**Done:**
- Wrote `docs/google_cloud_guide.md`: full step-by-step from project creation to pipeline run, including VM setup, conda, sra-tools version pinning (3.1.1 — newer versions segfault), and regtools argument order gotcha.

---

## 2026-04-03

**Goal:** Make alignment work on the GCP VM (8 GB RAM) without STAR.

**Done:**
- Added HISAT2 as a low-memory alternative to STAR. STAR requires ~30 GB RAM for hg38; HISAT2 indexes fit in ~8 GB.
- Added DAG diagram (PDF) for workflow visualisation.

---

## 2026-04-02

**Goal:** Replace NetMHCPan with an open-source MHC binding predictor.

**Done:**
- Replaced NetMHCPan (registration-gated, no programmatic access) with MHCflurry 2.x. MHCflurry is fully open-source, pip-installable, and produces IC50 predictions compatible with the 500 nM strong-binder threshold used in the original 2015 paper.

---

## 2026-03-25

**Goal:** Implement the modernised pipeline from scratch.

**Done:**
- Initial Snakemake pipeline implementing all steps: STAR/HISAT2 alignment → regtools junction extraction → GENCODE annotation filtering → normal-sample filtering → peptide translation → MHCflurry binding prediction → HTML report.
- Junction origin classification hierarchy: annotated → discard; unannotated + in normal → patient-specific (discard); unannotated + absent in normal → tumor-specific (predict). This replaces the original Fisher's exact test with an upstream biological filter.

---
