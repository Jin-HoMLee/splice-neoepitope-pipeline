# Lab Notebook

---

## 2026-04-10

**Goal:** Get TCRdock structural validation running end-to-end via Docker on the GCP GPU Spot VM.

**Done:**
- Fixed Docker build failure: the official TCRdock Dockerfile fails at `conda install openmm=7.7.0` (no longer available). Created `docker/Dockerfile.pipeline` which uses pip only and omits openmm/pdbfixer — these are only needed for optional Amber relaxation, not for `run_prediction.py`.
- Successfully built `tcrdock:latest` on GPU VM and ran TCRdock end-to-end. Pipeline completed 3/3 steps, report generated.
- Identified visual issue in report: AlphaFold outputs all residues as a single chain A, so Mol* rendered the ternary complex without chain-colour distinction. Added `relabel_pdb_chains()` to `run_tcrdock.py`, which reassigns chain IDs (A=MHC, B=peptide, C=TCR-alpha, D=TCR-beta) using the chain lengths from TCRdock's `alphafold_setup/targets.tsv`.
- Updated `setup_tcrdock_vm.sh` to build from `docker/Dockerfile.pipeline` instead of cloning TCRdock separately for the (broken) official Dockerfile.

**Next:** Re-run pipeline via `run_cloud_gpu.sh` to validate the chain-relabelling fix in the full report.

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
