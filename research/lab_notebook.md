# Lab Notebook

---

## 2026-04-28

### 13:30 UTC — Editor: Developer

#### Issue #181 — add research/glossary.md

Added `research/glossary.md` as a living dictionary of project-relevant abbreviations. Seeded with three entries that came up in conversation today (AFDB, GKE, SLURM). Format: alphabetical sections, `**ABBREV** — Expansion. One-line context. *Domain: bio | ml | cloud | pipeline | stats.*` Going forward, new abbreviations are added when they come up in conversation (rule saved in memory). Merged via PR #182.

---

## 2026-04-27

### 17:17 UTC — Editor: Developer

#### Issue #158 — rename config/gpu.yaml → config/gpu_config.yaml

Renamed `config/gpu.yaml` to `config/gpu_config.yaml` for naming consistency with other `*_config.yaml` files and to avoid confusion with conda env files in `workflow/envs/*.yaml`. Updated all live references across 9 files: `config/config.yaml`, `scripts/run_cloud_gpu.sh`, `scripts/setup_vm.sh`, `.github/workflows/tests.yml`, `README.md`, `CLAUDE.md`, `docs/configuration.md`, `docs/google_cloud_guide.md`. Historical entries in this notebook (lines 173, 572, 576, 708) intentionally left unchanged — they are point-in-time records.

Evaluated Snakemake profile migration (Issue #158 step 2) and deferred: `gpu_config.yaml` holds pipeline-level config values (TCRdock settings, fallback HLA/TCR alleles), not executor/resource parameters — the `--configfile` overlay is the correct abstraction.

Merged via PR #174 (replaced PR #170 which closed during a branch rename operation).

---

### 16:58 UTC — Editor: Scientist

#### PR #173 — patient_002 Junction Filtering correction merged

Reviewed by Claude (automated, LGTM, math verified). Both CI checks passed. Squash merged, branch deleted. Closes #172.

---

### 15:43 UTC — Editor: Scientist

#### Issue #172 — RESULTS.md patient_002 Junction Filtering corrected to valid post-#148 run

The Junction Filtering table previously reflected the invalid pre-#148 run (total 347,046, annotated 291,131, unannotated 55,915, normal-shared 3, tumor-exclusive 55,912). Corrected to the valid post-#148 run numbers derived from `alignment/BG003082_T0/junctions.tsv` and `reports/report.tsv` on GCS:

| Stage | Old (invalid) | New (valid) |
|-------|--------------|-------------|
| Total extracted | 347,046 | 364,168 |
| Annotated | 291,131 (83.9%) | 305,254 (83.8%) |
| Unannotated | 55,915 (16.1%) | 58,914 (16.2%) |
| Normal-shared | 3 (~0.0%) | 0 (0.0%) |
| Tumor-exclusive | 55,912 (16.1%) | 58,914 (16.2%) |

`normal_shared = 0` is correct and expected: BG003082_N0_WES was used as the normal input and yielded valid HLA typing results, but WES data contains no RNA splice junctions so junction-level normal subtraction was not effective. Also updated Dataset section to reflect this accurately.

Developer standup (2026-04-26 21:31 UTC, Q4) confirmed the valid run used no RNA-seq normal; the WES-yielding-0-junctions behaviour was separately clarified by the user.

---

## 2026-04-26

### 22:40 UTC — Editor: Developer

#### Issue #123 — M1 production cloud run: acceptance criteria verified (both patients)

Verified all GCS outputs against the Issue #123 acceptance criteria. All criteria passed.

| Criterion | patient_001 | patient_002 |
|---|---|---|
| Pipeline completes end-to-end without errors | ✅ | ✅ |
| `mhc_presentation.tsv` has breadth columns (`best_allele`, `genotype_presentation_score`, `n_strong_alleles`, `best_presentation_percentile`) | ✅ | ✅ |
| `genotype_presentation_score` in [0, 1] and scientifically sensible | ✅ [0.0147, 0.9999] — 1,286,492 rows | ✅ [0.0122, 0.9999] — 2,330,687 rows |
| `report.tsv` and `report.html` generated | ✅ | ✅ |
| Results uploaded to GCS | ✅ | ✅ |

Note: two column names in the acceptance criteria are stale — `strong_alleles` (actual: `n_strong_alleles`) and `presentation_percentile_strong` (actual: `best_presentation_percentile`). The equivalent columns are present; these are pre-rename leftovers and do not block closure.

Issue #123 ready to close.

---

### 22:18 UTC — Editor: Scientist

#### PR #168 — review fixes (docs/scientist/issue-165-patient-002-documentation)

Three items addressed before merge:

1. **Weak threshold label clarified** — `Weak (percentile ≤ 2.0%)` → `Weak (0.5% < percentile ≤ 2.0%)` to make the exclusive band unambiguous.
2. **GPS phrasing fixed** — `exceeded GPS > 0.9` → `had GPS > 0.9` (redundant double-negative removed).
3. **9-mer discrepancy** — added WES proxy normal usage as a second possible cause alongside the #148 fix.

---

### 21:50 UTC — Editor: Scientist

#### PR #167 — review fixes (feat/scientist/issue-164-patient-002-notebook)

Three issues addressed in code review before merge:

1. **GPS subsection headings moved to markdown cells** — headings were embedded as comments inside code cells; split into proper markdown cells and renumbered 5.1–5.4 → 6.1–6.4 to match the parent Section 6.
2. **Inflation check comment corrected** — comment previously said candidates *are* caught by the quality gate; corrected to say they are *not* (best_percentile ~0.5–0.55%, below the 2% threshold — GPS inflates from allele breadth, not per-allele strength).
3. **Allele list derived programmatically** — removed hardcoded patient-specific allele list; now auto-detected from `_presentation_score` column names, making the notebook reusable across patients.

---

### 21:28 UTC — Editor: Scientist

#### PR #167 — patient_002 results analysis notebook (feat/scientist/issue-164-patient-002-notebook)

Added `research/notebooks/patient_002_results.ipynb` with full scientific analysis of the post-#148 valid run:

- **Junction analysis:** 55,912 tumor-exclusive junctions (no matched RNA-seq normal; WES proxy filtered 3 junctions)
- **Peptide translation:** 703,106 (8-mer) / 781,159 (9-mer) / 846,422 (10-mer) = 2,330,687 total, 2,313,700 unique
- **MHC presentation:** 67,935 strong presenters (2.9%); HLA-C\*01:02 + HLA-C\*07:01 account for ~69% of strong-presenting candidates
- **GPS validation:** mean 0.101, median 0.026 — discriminating; 174 inflation edge cases (GPS > 0.9 but n_strong_alleles = 0, best_percentile ~0.5–0.55%)
- **Top candidate:** FADLRPLLL / HLA-C\*01:02 (IC50 = 33.2 nM, percentile = 0.0045%, GPS = 0.9999, n_strong = 4)

Open questions posted to Developer in team standup: min read support filter, GPS n_strong_alleles gate, HLA-C calibration, WES normal confirmation.

---

### 21:24 UTC — Editor: Scientist

#### PR #166 — Jupyter notebook environment setup (feat/scientist/issue-163-notebook-env)

Merged notebook environment scaffold for `research/notebooks/`:

- `research/requirements.txt` — lower-bound pins: jupyterlab ≥ 4, pandas ≥ 2, matplotlib ≥ 3.7
- `.python-version` added to `.gitignore` (personal env config, not committed); `pyenv local 3.14.4` kept in README as the explicit setup step
- `research/README.md` updated with setup instructions: pyenv → venv → pip install → VSCode kernel selection; GCS auth note

---

### 18:15 UTC — Editor: Scientist

#### Issues #162–165 — patient_002 results analysis (valid run, post Issue #148)

First scientific analysis session for patient_002 (BG003082, osteosarcoma). Full analysis in `research/notebooks/patient_002_results.ipynb`.

**Key findings:**

- Top candidate: **FADLRPLLL / HLA-C\*01:02** — IC50 33.2 nM, presentation_percentile 0.0045%, GPS 0.9999, strong presenter across 4/5 alleles.
- **HLA-C dominance:** HLA-C\*01:02 + HLA-C\*07:01 account for ~69% of all 67,935 strong-presenting peptides despite 0.5× locus weight in GPS. Three hypotheses: (a) broader HLA-C groove/promiscuity, (b) osteosarcoma splice motifs enriched for HLA-C, (c) MHCflurry percentile calibration less precise for HLA-C (sparser training data). Sent to Developer for investigation.
- **HLA-A\*01:01 nearly silent:** median presentation percentile 8.5% among strong presenters — contributes almost nothing to GPS for this patient.
- **GPS is discriminating:** mean 0.101, median 0.026; only 1.1% of 2.3M predictions exceed GPS > 0.9. Formula works as designed.
- **GPS inflation edge case (minor):** 174 candidates (0.7% of GPS > 0.9) have `n_strong_alleles = 0` with best_percentile ~0.5–0.55%. Quality gate (2%) does not catch these. Flagged to Developer to consider adding `n_strong_alleles ≥ 1` as a hard filter.
- **Minimum junction read support is 174** — suspiciously high. Awaiting Developer confirmation of whether a min-reads filter is applied.
- **Peptide translation:** 2,330,687 total (703,106 × 8-mer, 781,159 × 9-mer, 846,422 × 10-mer), 2,313,700 unique. 9-mer count differs by 265 from a prior run (781,424) — untracked discrepancy, motivates run registry (Issue TBD).

**Actions taken:**

- Created `research/notebooks/patient_002_results.ipynb` (Issue #164).
- Set up `research/.venv` (Python 3.14.4 via pyenv) + `research/requirements.txt` + `research/README.md` notebook section (Issue #163).
- Updated `research/manuscript/RESULTS.md` patient_002 section: MHC presentation + GPS, peptide translation, top candidate FADLRPLLL, TCRdock. Junction Filtering table on hold pending Developer confirmation of WES normal usage (Issue #165).
- Opened Issue #161 (local GCS cache for notebooks); drafted run registry issue for Developer review in standup.
- Terminology convention established: "presenters/presenting" replaces "binders/binding" throughout.

**Open / pending:**
- Developer standup: min read filter, GPS gate, HLA-C calibration, WES normal confirmation, run registry engineering input.
- RESULTS.md Junction Filtering: on hold.

---

### 17:03 UTC — Editor: Developer

#### Issue #148 — PR #152 review fixes (assemble_contigs.py + alignment.smk)

Two issues raised in code review addressed before merge:

1. **Silent failure mode fixed (`assemble_contigs.py`)** — `bedtools getfasta` stderr is now always logged as WARNING (was only logged on non-zero exit, so chromosome-not-found warnings were silently discarded). Added a skip-rate check at the summary: if >10% of junctions are skipped due to short sequences, the log is upgraded to WARNING with an explicit message pointing to UCSC/ENSEMBL chr naming as the likely cause.

2. **Docstring clarified (`alignment.smk`)** — `hisat2_download_index` rule docstring and inline comment updated from "GRCh38" (ambiguous) to "hg38 (UCSC naming)" to remove the ambiguity that caused Issue #148 in the first place.

---

### 15:25 UTC — Editor: Developer

#### Issue #148 — patient_002 rerun results (after hg38_tran HISAT2 index fix)

Full pipeline rerun on patient_002 after merging `fix/issue-148-hisat2-chr-naming` (PR #152 — https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/152), which switched the HISAT2 prebuilt index from `grch38_tran` (ENSEMBL naming, no `chr` prefix) to `hg38_tran` (UCSC naming, `chr` prefix) to match the genome FASTA. Previous results were invalid due to the chromosome naming mismatch causing `bedtools getfasta` to return empty sequences for 56,735 of 56,774 junctions.

**Junction filtering**

| Metric | Count |
|---|---|
| Unannotated junctions | 58,914 |
| Tumor-exclusive | 58,914 |
| Normal-shared | 0 |

No matched normal sample for patient_002 → all unannotated junctions labeled tumor_exclusive (expected behaviour with a warning in the pipeline log).

**MHC predictions**

| Class | Count |
|---|---|
| Total predictions | 2,330,687 |
| Strong binders (≤ 0.5%) | 67,935 |
| Weak binders (≤ 2.0%) | 222,823 |
| Non-binders (> 2.0%) | 2,039,929 |

**HLA typing (OptiType, tumor)**

| Locus | Allele(s) |
|---|---|
| HLA-A | A\*01:01 |
| HLA-B | B\*08:01 / B\*27:05 |
| HLA-C | C\*01:02 / C\*07:01 |

**Top neoepitope candidate**

| Field | Value |
|---|---|
| Peptide | FADLRPLLL |
| Best allele | HLA-C\*01:02 |
| Presentation percentile | 0.0045% |
| Genotype presentation score | 0.9999 |
| IC50 | 33.2 nM |
| Presentation class | strong |

TCRdock PDB available: yes.

---

### ~14:54 UTC — Editor: Developer

#### Issue #118 — CI dry-run with GPU config overlay (PR #159)

Added a `pipeline-snakemake-dry-run` CI job to `.github/workflows/tests.yml` that runs `snakemake -n --configfile config/config.yaml config/gpu.yaml` on every push/PR to main. This catches config key errors in `structure.smk` TCRdock rule blocks — previously invisible because `gpu.yaml` was never loaded in CI. The root motivation was the `KeyError: 'ic50_strong'` on the first cloud prod run (fixed in #117), which this job would have caught. Placeholder input files (`touch`) are created before the dry-run so Snakemake can build the DAG without real data. Also renamed the existing `test` job to `pipeline-pytest` for clarity, and updated the branch protection required status checks accordingly.

---

### ~10:01 UTC — Editor: Scientist

#### Issue #153 — Add --update-note flag to zotero_add.py

Added `update_note()` helper and `--update-note ITEM_KEY` flag to `research/scripts/zotero_add.py`. Allows updating the note on an existing Zotero entry without re-adding the paper — needed for the morning reading routine when the note format evolves after initial add (e.g. adding a Limitations section). Uses optimistic locking via `If-Unmodified-Since-Version`. Falls back to creating a new note if none exists.

---

### ~09:30 UTC — Editor: Scientist

#### Issue #154 — Move zotero_add.py to research/scripts/ and add research/README.md

Moved `scripts/zotero_add.py` to `research/scripts/zotero_add.py` to clarify that the script is a research support tool, not part of the bioinformatics pipeline. Created `research/README.md` with a folder overview and full usage docs for `zotero_add.py` including the three-section note format convention (Results / Methods / Limitations).

---

### ~08:30 UTC — Editor: Developer

#### Issue #148 — patient_002 full rerun completed (first-pass success)

Full rerun of patient_002 completed successfully after merging the `hg38_tran` HISAT2 index fix (Issue #148, branch `fix/issue-148-hisat2-chr-naming`). All Snakemake steps finished without errors. At first view the results look good — junction counts and peptide counts appear in the expected range, unlike the previous invalid run (which produced only 783 peptides due to the ENSEMBL/UCSC chr naming mismatch). Detailed result review and top candidate recording deferred to next session.

---

## 2026-04-25

### ~16:56 UTC — Editor: Developer

#### Issue #136 — VM data cleanup + NVIDIA driver fixes for common-cu129 image

**Context:** patient_002 run restarted after morning disk resize. Left running while user went to lunch (~12:00 UTC).

**NVIDIA driver failures — three iterations (~13:00–15:00 UTC)**

The `ubuntu-accelerator-2204-amd64-with-nvidia-570` Deep Learning VM image was retired by Google; the replacement is `common-cu129-ubuntu-2204-nvidia-580`. The new image ships driver 580, which lacks GSP firmware support for P100 (Pascal SM 6.0). Three driver downgrade attempts were needed before a working configuration was found:

| Attempt | Package | Failure reason |
|---------|---------|----------------|
| 1 | `nvidia-driver-570-server` | Pulls display dependencies (`libgbm1`, `libxcb-*`, `libnvidia-egl-wayland1`) absent on headless VMs |
| 2 | `nvidia-headless-no-dkms-570-server` | Pre-compiled kernel modules don't match GCP kernel (`6.8.0-1053-gcp`) on `common-cu129` image — `nvidia-smi` fails with "driver not loaded" |
| 3 | `nvidia-headless-570-server` (DKMS) | **Works** — DKMS recompiles modules against the running kernel at install time |

`IMAGE_FAMILY` in `run_cloud_gpu.sh` updated to `common-cu129-ubuntu-2204-nvidia-580`. Driver downgrade block updated to install the DKMS variant. Three commits on `feat/issue-136-vm-data-cleanup`; PR #151 opened.

**patient_002 pipeline completed (~15:00–15:46 UTC)**

TCRdock rerun manually with `snakemake --forcerun run_tcrdock` after DKMS driver confirmed active. All 22 Snakemake steps finished. Results and logs uploaded to GCS at ~15:46 UTC.

Top candidate: **NSISRPSSL / HLA-C\*01:02, ic50\_nM=39.77, pLDDT=91.99**

**Issue #148 discovered — results invalid: chromosome naming mismatch (~16:00 UTC)**

Report showed only ~800 total neoepitopes (expected thousands). Root cause: `hisat2_prebuilt_url` in `config/config.yaml` used `grch38_tran.tar.gz` (ENSEMBL naming, no `chr` prefix) while the genome FASTA uses UCSC naming (`chr` prefix). `bedtools getfasta` returned empty sequences → 56,735 of 56,774 junctions skipped in `assemble_contigs.py` → only 18 contigs → 783 peptides. patient_001 was unaffected as it predates the config key and built the index locally from the UCSC FASTA. Fix: `hg38_tran.tar.gz` on `fix/issue-148-hisat2-chr-naming`. **Today's patient_002 results are invalid; full rerun required after Issue #148 is merged.**

---

### ~15:00 UTC — Editor: Scientist

#### TCR-pMHC binding prediction — field overview and structural improvement plan (Issue #86)

**Field overview:**

Two main paradigms:
- **Sequence/ML-based** (NetTCR-2.x, pMTnet, ERGO, MixTCRpred, TULIP): treat TCR+peptide as a sequence classification problem; generalisation to unseen epitopes remains a key limitation across all methods — to revisit in a future Scientist session
- **Structure-based** (TCRdock, AlphaFold3): model the 3D TCR-pMHC complex; TCRdock (Alam et al., *Science* 2023) is already integrated; AF3 (Abramson et al., *Nature* 2024) is a newer competitor

Key databases: VDJdb (curated TCR-pMHC specificity), IEDB, McPAS-TCR.

**Three axes for improving the structural approach — agreed execution order:**

1. **Better inputs (Issue #86):** patient-specific HLA-matched TCR panel from VDJdb — prerequisite for all downstream improvements
2. **Model upgrade:** benchmark AlphaFold3 vs. TCRdock on a known-binding panel once Axis 1 is in place
3. **Rescoring:** post-process PDB outputs with Rosetta `InterfaceAnalyzer` or FoldX `AnalyseComplex` for interface ΔΔG as a secondary ranking signal

**Issue #86 — VDJdb TCR panel design (first-pass, conservative):**

| Parameter | Value |
|---|---|
| HLA matching | Exact 4-digit |
| MHC class | Class I only |
| Paired α/β chains | Required |
| VDJdb confidence score | ≥ 2 |
| Panel size | Top 10 per allele |
| Redundancy reduction | None |
| Antigen category filter | None |

Key non-obvious dependency: **`stitchr`** (Peacock et al. 2022, *Bioinformatics*) — VDJdb stores CDR3 + V/J gene assignments only; TCRdock needs full α/β sequences; stitchr reconstructs them from IMGT germline references.

Pipeline integration: new Snakemake rule `fetch_vdjdb_panel` (input: patient HLA alleles; output: `results/{sample}/tcrdock/vdjdb_panel.tsv`) feeds into existing `run_tcrdock`. VDJdb data from `antigenomics/vdjdb-db` flat TSV (pinned release). Issue #86 body updated to reflect this design.

---

### ~12:00 UTC — Editor: Scientist

#### Zotero integration for morning science reading routine (Issue #137, PR #138)

Established daily science reading habit for Scientist sessions. Each morning warm-up ("good morning") now produces a Zotero entry rather than a markdown log.

**Setup:**
- Created Zotero collection "Splice Neoepitope Pipeline" (key `Z38GTJNW`, library `lee.jin-ho`, user ID 10082130)
- `scripts/zotero_add.py`: CLI tool to add a paper by DOI — fetches metadata from CrossRef (authors, journal, ISSN) with PubMed fallback for full date, structured abstract, and PMID; supports `--note` (child note), `--tags`, `--dry-run`
- Credentials in `.env` (gitignored)

**First entry:**
Weber et al. (2024) KEYNOTE-942, *The Lancet* 403:632–644. DOI: 10.1016/S0140-6736(23)02268-7. Phase 2b trial: personalised mRNA neoepitope vaccine (mRNA-4157/V940) + pembrolizumab → ~44% reduction in recurrence/death vs. pembrolizumab alone in resected melanoma. Targets SNV/indel neoantigens only; splice-junction neoepitopes (our focus) are absent — directly motivates this project.

**PR #138 merged.**

---

### ~10:30 UTC — Editor: Developer

#### Issue #123 — M1 production cloud run: sub-issue retrospective, VM disk fix, branch rebase

**patient_002 run failure — root cause investigation (~09:30 UTC)**

patient_002 run launched the previous evening (~22:12 UTC, 2026-04-24) got stuck overnight.
SSH inspection of the VM showed `samtools sort` hanging with no output: the HISAT2 index was
missing (`resources/hisat2_index/genome_tran` not found), causing HISAT2 to exit immediately.
`samtools sort` was reading from HISAT2's stdout pipe and hung instead of propagating the error.
`set -euo pipefail` was insufficient because the right-hand side of the pipe does not exit when
the left-hand side fails — it just blocks waiting for more input.

Immediate fix: killed the tmux session, deleted the stale `resources/hisat2_index/` directory
(contained `genome.*.ht2` from an old build, but `index.done` sentinel was present so Snakemake
skipped re-download). Pipeline restarted by user after disk issue resolved (see below).

Code fix (Issue #131, commit `e2b1a65`): added a pre-check at the top of `hisat2_align` shell
block — exits with a clear error and writes to `{log}` before `samtools` is ever launched:
```bash
if [[ ! -f "{params.index_prefix}.1.ht2" ]]; then
    echo "ERROR: HISAT2 index not found at {params.index_prefix}.*.ht2" | tee -a {log}
    exit 1
fi
```
This ensures the error propagates to `pipeline.log` and the orchestrator can detect it.

**VM disk full**

VM (100 GB SSD) was full — patient_001 data had not been cleaned up post-GCS upload. 300 GB
resize rejected (europe-west1 SSD quota is 250 GB; `neoepitope-pipeline` 200 GB +
`orchestrator` 10 GB = 210 GB in use). Resized to 200 GB (`gcloud compute disks resize`),
then `sudo resize2fs /dev/sda1` to extend the filesystem. patient_001 data deleted from VM
(already safe in GCS).

**Sub-issue retrospective — PRs #132–#135 merged**

10 commits accumulated on `feat/issue-123-prod-cloud-run` since patient_001 re-run were
retrospectively organised into 4 focused sub-issues, each cherry-picked to its own branch
and merged into main independently:

| Issue | PR | Commit(s) | Description |
|-------|----|-----------|-------------|
| #128 | #132 | `6dc0fdc` | fix(cloud): upload pipeline.log + orchestrator.log to GCS |
| #129 | #133 | `8d3ff15` | fix(tcrdock): match `_pdb_file` suffix to avoid pdbid false-match |
| #130 | #134 | `57c1f3f` | fix(cloud): remove `--conda-cleanup-envs` (deletes all envs on empty DAG) |
| #131 | #135 | `e2b1a65` | fix(alignment): fail fast on missing HISAT2 index |

`feat/issue-123-prod-cloud-run` now contains only 3 docs commits
(`bbd88c5`, `f2eb363`, `86e7e15`) on top of the latest main. Branch rebased and
force-pushed. PR for this branch will be opened after patient_002 run completes.

**patient_002 run status:** restarted after disk resize; in progress.

---

## 2026-04-24

### ~22:12 UTC — Editor: Developer

#### Issue #123 — M1 production cloud run: patient_002 started (WES normal, interim)

Launched production pipeline run for patient_002 (osteosarcoma IPISRC044, BG003082) via
`run_cloud_gpu.sh --detach` on `neoepitope-pipeline` VM. Branch: `feat/issue-123-prod-cloud-run`.

Normal sample: Blood Derived WES (`BG003082_N0_WES`) — interim proxy; yields near-zero junction
overlap (~3 junctions). All unannotated tumour junctions will be labeled `tumor_exclusive` with
a pipeline warning. Proper normal (Jan 2025 CD3+ PBMC Cell Ranger BAM) pending Issue #127
implementation.

Also includes two bug fixes committed since patient_001 re-run:
- `2625ec7` fix(tcrdock): match `_pdb_file` suffix in column search (avoids `pdbid` false-match)
- `117d174` fix(cloud): remove `--conda-cleanup-envs` (was deleting all 4 envs when DAG empty)

Issue #127 opened: support pre-aligned BAM as normal input for junction filtering.

---

### ~22:00 UTC

#### Issue #123 / #126 — patient_002 normal sample decision + GTEx filter design (Researcher session)

**Context:** patient_002 (osteosarcoma IPISRC044) has no matched RNA-seq normal. Current
production run uses Blood WES as proxy, yielding effectively no junction filtering (3 junctions
overlapping tumour set — consistent with WES spliced reads being alignment artefacts).

**Normal sample decision for current production run:**

Longitudinal single-cell RNA-seq (PBMC, sorted CD3+ T cells) is available from the Hudson Lab:

- Protocol: 10x Chromium (confirmed by Cell Ranger output + Seurat objects)
- Time points: Jan 2025 – Dec 2025 (growing dataset)
- Location: `hudson_lab/PBMC_scRNAseq/FASTQ` and `hudson_lab/PBMC_scRNAseq/cellranger`

Decision: use **Jan 2025 Cell Ranger BAM** as the normal input for the current patient_002 run.
Rationale:
- CD3+ T cells are the primary TIL population contaminating solid tumour RNA-seq; filtering their
  junctions directly targets the most relevant contamination source
- Jan 2025 is the earliest (most treatment-naive) time point — least treatment-modified T cell state
- 10x 3' capture gives sparse junction coverage, but still produces genuine biological splice junctions
  unlike WES (which gave near-zero overlap)
- This is an interim measure pending GTEx pan-tissue filter implementation (issue #126)

Developer action: run regtools on the Jan 2025 Cell Ranger BAM (already genome-aligned; no
re-alignment needed) and use resulting junction set as `normal_junctions` for patient_002.

**GTEx pan-tissue filter — scientific rationale documented (issue #126):**

A pan-tissue GTEx filter (not tissue-matched) is the correct long-term approach for a vaccination
application. Key argument: vaccine-induced CTLs are systemic — they patrol all tissues, not only
the tumour site. A junction present in any normal tissue could be presented on those cells, creating
off-tumour autoimmune toxicity risk. Pan-tissue filtering serves dual purpose:

1. Safety — excludes junctions expressed in any of ~54 GTEx tissue types
2. Candidate quality — pan-tissue absence is a stronger tumour-exclusivity claim; precision over
   recall is appropriate given 10–20 vaccine slot constraint

Full rationale added to `research/manuscript/DISCUSSIONS.md` (section "Normal sample filtering →
GTEx pan-tissue filter"). Issue #126 opened for Developer implementation.

---

### ~18:30 UTC

#### Issue #123 — M1 production cloud run: patient_001 started

Launched production pipeline run for patient_001 (gastric cancer, SRR9143066 T / SRR9143065 N) via `run_cloud_gpu.sh --detach` on `neoepitope-pipeline` VM (n1-highmem-8 + P100, europe-west1-b). Branch: `main`. patient_002 run to follow after patient_001 completes.

---

### ~17:30 UTC

#### Issue #124 — CLI flag polish (Developer session, PR #125, branch `feat/issue-124-cli-flag`)

**Goal:** Expose `--presentation-percentile-weak` via CLI parsers in both `generate_report.py` and `run_tcrdock.py`; fix docstring omission.

**Done:**

- Added `--presentation-percentile-weak` (type=float, default=2.0) to `_cli_main()` in both scripts and threaded through to the respective `generate_report()` / `run_structural_validation()` calls.
- Fixed `generate_report()` docstring: `presentation_percentile_weak` was present in the function signature but absent from the `Args:` section.
- 4 new `TestCLIParser` tests (default + custom value, one class per script). 200 tests passing.

**Commits (2, branch `feat/issue-124-cli-flag`):**
- `6cea0dd` feat(cli): add --presentation-percentile-weak flag to generate_report and run_tcrdock CLI parsers
- `7c6dd35` test(cli): add CLI parser tests for --presentation-percentile-weak flag

**PR #125 opened.** Awaiting merge.

---

### ~16:42 UTC

#### Issue #119 — PR #122 re-review fix + merge prep

**Re-review verdict:** Ready to merge. Three minor observations: (1) `presentation_notes` strings in `_build_report_tsv` still hardcoded `"2%"` for `weak`/`non` entries despite `presentation_percentile_weak` parameter being available; (2) CLI parsers missing `--presentation-percentile-weak` flag; (3) `generate_report()` docstring missing the parameter.

Fixed (1) now — `weak`/`non` entries converted to f-strings using `presentation_percentile_weak` (`18e1fac`). Items (2) and (3) deferred to follow-up issue.

**29 tests passing** (generate_report suite).

---

### ~16:15 UTC

#### Issue #119 — PR #122 review cycle (Developer session)

**Reviewer comments (6):** stale docstring, GPS out-of-range clamping, hardcoded quality-gate threshold, GPS test gap in `generate_report.py`, hardcoded 9-mer `end_nt_incl`, misleading log message on quality-gate rejection.

**Key design decisions made during review:**
- `hla_c_weight` outside [0,1] → hard `ValueError` (own config, must be correct).
- MHCflurry `presentation_score` outside [0,1] → warning only, no clamp (third-party output; out-of-range values cannot be reliably interpreted and should be fixed upstream).
- "binder" → "presenter" terminology throughout `run_tcrdock.py` and tests.

**Fix commits (8):**
- `4366af0` fix(mhc): fix stale docstring and add hla_c_weight validation + score warning
- `9f0b6c6` test(mhc): add tests for hla_c_weight validation and out-of-range score warning
- `b73b27b` fix(rules): thread presentation_percentile_weak into report and tcrdock params
- `94a480a` fix(report): replace hardcoded quality-gate 2.0 with configurable presentation_percentile_weak
- `34f7aa5` fix(tcrdock): replace hardcoded quality-gate 2.0 with configurable presentation_percentile_weak
- `1b08ebe` test(report): add GPS quality-gate and ranking tests for generate_report.py
- `fe6e188` fix(report): compute end_nt_incl from peptide length instead of hardcoded 9-mer
- `d0e7456` fix(tcrdock): distinguish no-presenters from quality-gate-empty in log messages

**196 tests passing.** Re-review requested.

---

### ~14:38 UTC

#### Issue #119 — Two-level MHC presentation prediction: implementation (Developer session, PR #122, branch `feat/issue-119-allele-breadth`)

**Goal:** Implement the allele breadth model designed in the Researcher session (~11:30 UTC) and open PR #122.

**Done:**

- Replaced `_compute_strong_alleles()` in `run_mhcflurry.py` with `_compute_per_allele_features()`: one `Class1PresentationPredictor.predict()` call per allele (single-element genotype list) to recover per-allele `presentation_score` and `presentation_percentile` discarded by the genotype API.
- New `genotype_presentation_score` formula: `1 − ∏(1 − wᵢ·pᵢ)`, with `w(HLA-A) = w(HLA-B) = 1.0`, `w(HLA-C) = hla_c_weight` (default 0.5). New supporting columns: `n_strong_alleles`, `best_presentation_percentile`.
- `allele` column renamed to `best_allele` throughout (`run_mhcflurry.py`, `generate_report.py`, `run_tcrdock.py`, tests).
- Updated ranking in `generate_report.py` and `run_tcrdock.py`: GPS ↓ → `n_strong_alleles` ↓ → `best_presentation_percentile` ↑. Quality gate: `best_presentation_percentile > 2%` excluded from top-candidates list.
- Added `mhcflurry.hla_c_weight: 0.5` to `config/config.yaml` and propagated through `mhc_affinity.smk` params.
- `METHODS.md` Section 6 rewritten: two-level architecture, GPS formula in LaTeX, quality gate, ranking rule, full output column table.
- `docs/configuration.md`: `hla_c_weight` documented with HLA-C surface density rationale.
- Test suite: 186 passed, 1 skipped (stale integration pipeline output — schema-version skip guard added).

**Key implementation detail:** The genotype API (`predict()` with all alleles at once) reports only the best-allele scores and silently discards all others. Per-allele scores must be recovered via separate single-allele calls. This dual-call design is why `_compute_per_allele_features` iterates over `resolved_alleles` individually.

**Column name change mid-session:** `breadth_score` → `genotype_presentation_score` after Researcher consultation. All code, tests, and docs use the new name.

**Local test run (chr22, patient_001_test):** Full pipeline completed. `mhc_presentation.tsv` has 24 columns (6 allele pairs + 3 breadth columns). Top hit: DVFGTPFSR / HLA-A\*33:03, GPS = 0.9996, best-allele percentile = 0.001%.

**Commits (9, all on `feat/issue-119-allele-breadth`):**
- `530ddfc` docs(manuscript): rename breadth_score → genotype_presentation_score
- `a13fd30` feat(config): add mhcflurry.hla_c_weight config key and Snakemake param
- `38b6ac8` feat(mhc): add per-allele genotype_presentation_score model (Issue #119)
- `49b8e96` feat(report): update ranking to genotype_presentation_score, add quality gate
- `b90c75d` feat(tcrdock): update candidate selection to genotype_presentation_score ranking
- `291fd23` test(mhc): update tests for new output schema and add breadth feature tests
- `3d8b903` test(report,tcrdock): update test data to best_allele; add breadth-aware tcrdock tests
- `260da3a` test(integration): update required columns for new schema with skip guards
- `f792ac1` docs(methods): update Section 6 for two-level MHC presentation prediction
- `ac89bd2` docs(config): document hla_c_weight and update mhcflurry section

---

### ~11:30 UTC

#### Issue #119 — Allele breadth model: scientific design (Researcher session, branch `feat/issue-119-allele-breadth`)

**Goal:** Develop and document the scientific rationale for a multi-allele breadth scoring
model for `mhc_presentation.tsv`. No code changes — manuscript only.

**Key scientific decisions:**

- **Two-level architecture:** MHCflurry is a molecular-level predictor (one peptide × one
  MHC allele → `presentation_score`). The breadth model is a separate genotype-level
  combiner: `breadth_score = 1 − ∏(1 − wᵢ·pᵢ)`. These are distinct modelling layers.
- **HLA-C weight:** enters at the genotype-combination level (surface density ~50% of
  HLA-A/B), not at the molecular prediction level. Configurable (`w_C`, default 0.5).
- **`presentation_score` in formula (not `presentation_percentile`):** `presentation_score`
  is a calibrated absolute probability — the correct input for a complementary probability
  formula. `presentation_percentile` is a rank statistic; using it would require an
  arbitrary mapping. The two metrics can disagree (high breadth_score with all alleles in
  non-binder percentile territory is possible for promiscuous alleles) — this is the
  intended behaviour of the two-tier system.
- **Committed ranking for cancer vaccine application:**
  1. `breadth_score` — primary (multi-allele coverage, LOH robustness, vaccine slot efficiency)
  2. `n_strong_alleles` — secondary (count of alleles at ≤ 0.5% percentile)
  3. `best_presentation_percentile` — minimum quality gate only (not a ranking dimension);
     filters candidates where no allele reaches weak-binder territory
- **Immunodominance acknowledged but not dominant in vaccine context:** in natural
  anti-tumour immunity, one very strong allele can dominate via intramolecular MHC groove
  competition and immunodomination (Yewdell & Bennink, *Annu Rev Immunol* 1999; Chen &
  McCluskey, *Adv Cancer Res* 2006). In therapeutic vaccination, immunodomination is largely
  bypassed; HLA LOH and vaccine slot efficiency make breadth the primary criterion.

**Manuscript changes (committed `c9d5eb3`, branch `feat/issue-119-allele-breadth`):**
- `INTRODUCTION.md`: new section "HLA Genotype, Surface Expression, and Allele Breadth"
- `DISCUSSIONS.md`: new section "Allele breadth and immunodominance: two complementary
  ranking signals" (subsections: two-level architecture, immunodominance mechanisms,
  vaccine application with committed table, calibration note)

**METHODS.md:** not updated — Developer session responsibility after implementation.

**Design spec for Developer session:** `memory/project_allele_breadth_design.md`

---

### ~10:02 UTC

#### Issue #115 — Rename `mhc_affinity.tsv` → `mhc_presentation.tsv` (PR #121, branch `feat/issue-115-rename-mhc-affinity-tsv`)

**Motivation:** The output file was named after `Class1AffinityPredictor`, but the pipeline now uses `Class1PresentationPredictor` (since Issue #85). The name `mhc_affinity.tsv` was misleading — the primary ranking metric is `presentation_percentile`, not IC50.

**Changes:** Pure mechanical rename across 9 files — output path in `mhc_affinity.smk`, Snakemake output key `mhc_affinity_tsv` → `mhc_presentation_tsv` in `mhc_affinity.smk`, `run_mhcflurry.py`, `structure.smk`; path references in `analysis.smk`, `generate_report.py`, `run_tcrdock.py`, `test_integration.py`; docs in `METHODS.md`, `configuration.md`, `README.md`. The `.smk` rule module filename (`mhc_affinity.smk`) is unchanged — it is the rule module for the MHC prediction step, not named after the output.

**Testing:** 25/25 integration tests pass. Local test result file renamed from `mhc_affinity.tsv` → `mhc_presentation.tsv`.

---

### ~08:49 UTC

#### Hotfix — rename `ic50_strong` → `presentation_percentile_strong` in `structure.smk` (PR #117, branch `fix/structure-smk-ic50-key`)

**Problem:** First prod cloud run after Issue #85 failed immediately with `KeyError: 'ic50_strong'` in `workflow/rules/structure.smk` line 77. The `generate_report_with_structure` rule was still referencing the removed config key.

**Root cause:** The TCRdock rule block in `structure.smk` is guarded by `if _TCRDOCK_ENABLED:`, which is only true when `gpu.yaml` is merged in. Local tests and CI never enable TCRdock, so the stale key was never evaluated during development.

**Fix:** One-line rename of `ic50_strong` → `presentation_percentile_strong` in the `params:` block of `generate_report_with_structure`, matching the already-corrected `analysis.smk`.

**Follow-up filed:** Issue #118 — add a Snakemake dry-run with `gpu.yaml` overlay to CI to catch this class of bug automatically.

**Prod cloud run** (patient_001, full genome, GPU): completed successfully after fix. UTC: 2026-04-24T00:12:39Z. Key results vs. baseline (pre-Issue #85):
- Junction counts unchanged: 30,029 unannotated → 27,348 tumor-exclusive
- Total predictions: 7,718,952 → 1,286,492 (6× reduction — one row per peptide, not per peptide×allele)
- Strong presenters: 15,880 (IC50 < 50 nM) → 44,916 (presentation_percentile ≤ 0.5%)
- Top candidate changed: EVAETLSLF/HLA-A\*26:01 (IC50 16.6 nM) → FAFPFAQTL/HLA-C\*03:03 (percentile 0.0007%)

---

### ~00:30 UTC

#### Issue #85 — Switch to Class1PresentationPredictor (PR #114, branch `feat/issue-85-prediction-mode`)

**Goal:** Replace `Class1AffinityPredictor` with `Class1PresentationPredictor`, which integrates binding affinity with antigen processing (proteasomal cleavage, TAP transport) trained on mass-spec MHC ligand data.

**API surprises encountered:**
1. `predict_to_dataframe()` does not exist on `Class1PresentationPredictor` (affinity predictor only) — use `predict()` instead.
2. `Class1PresentationPredictor.predict()` is a *genotype-level* API: takes all patient HLA alleles at once (≤6) and returns **one best-allele prediction per peptide**. Passing the same allele repeated N times (the `AffinityPredictor` convention) raises `ValueError: alleles list must have at most 6 elements`.

**Design decisions:**
- Dropped `affinity_percentile` / `binder_class`: `predict()` does not return `affinity_percentile`, and a second inference pass would double runtime with no clear biological gain.
- Single `presentation_class` label from `presentation_percentile`: strong ≤ 0.5%, weak ≤ 2%, non > 2% — thresholds from Jiang et al. 2024 (*Communications Biology*).
- Merged `_load_predictor_for_gpu()` / `_load_predictor_for_cpu()` into `_load_predictor()`: both were identical after parallel-worker removal.

**Final output schema (`mhc_affinity.tsv`):** `contig_key, start_nt, peptide, allele, ic50_nM, processing_score, presentation_score, presentation_percentile, presentation_class`. One row per peptide; `allele` = best presenter in patient genotype.

**Local test run** (patient_001_test, chr22, 6-allele genotype): 4045 unique peptides → 147 strong / 401 weak / 3571 non. UTC: 2026-04-23T21:43:51Z.

182 tests pass (unit + integration).

---

## 2026-04-23

### ~20:00 UTC

#### Issue #110 — Fix stale integration test assertions (PR #113, branch `fix/issue-110-contig-length-test`)

**Goal:** Update two integration test assertions that became stale after PR #99 changed flank size to `3 × (max_length − 1)` and introduced multi-length peptides.

**Changes:**
- `test_all_contigs_are_50nt` → `test_all_contigs_are_54nt`: contig length is now 54 nt (= 27 + 27 flanks with `peptide_lengths=[8, 9, 10]`, `max_length=10`)
- `test_all_peptides_are_9mers` → `test_all_peptides_are_8_to_10mers`: pipeline produces 8-, 9-, and 10-mer peptides per config

No functional code changes; tests updated to match current pipeline behaviour.

**Note on branch base:** this branch was originally created off `feat/issue-97-cds-reading-frame` before #112 merged, making the PR appear to contain ~300 lines of #97 reading frame code. After #112 merged to main, the branch was rebased onto main — the diff is now exactly the 4-line test fix.

---

### ~19:30 UTC

#### PR #112 review fixes — issue #97

- Removed redundant `gene_cds` intermediate `defaultdict` from `_build_cds_donor_lookup`; donor frames now computed in a single GTF parsing pass (halves peak memory on large GTFs)
- Fixed misleading comment on GTF end coordinate: `# GTF end is 1-based inclusive = 0-based exclusive (no adjustment needed)`
- Added `test_minus_strand_reading_frame_annotated` to `TestClassifyJunctionsReadingFrame` — covers the `donor_coord = junction.end` path through the full `classify_junctions` pipeline for a − strand junction
- 157 tests pass

### ~18:30 UTC

#### Issue #97 — CDS reading frame annotation for novel junctions (branch `feat/issue-97-cds-reading-frame`)

**Goal:** Annotate each tumor-exclusive junction in `novel_junctions.tsv` with its canonical CDS reading frame, derived from the GENCODE GTF. This supports downstream biological interpretation without restricting translation.

**Key design decisions:**

1. **Annotation only, not restriction.** The `reading_frame` column is informational metadata. All three sense-strand frames are still translated downstream. Restricting translation to the CDS-derived frame would introduce false negatives in hypermutated tumors (MSI-high, POLE-mutant) where upstream frameshift indels, SVs, or chained novel junctions shift the active frame — precisely the tumors most likely to harbour actionable junction neoepitopes. These frame-altering events cannot be inferred from RNA-seq junction data alone.

2. **Protein-coding transcripts only.** CDS records are filtered to `transcript_type "protein_coding"` in the GENCODE GTF. NMD and other non-translated isoforms are excluded to avoid spurious frame assignments from transcripts that are not translated by the ribosome.

3. **Union of frames across transcripts.** When a donor site is annotated in multiple protein-coding transcripts with different reading frames, all attested frame offsets are retained (e.g. `"0,2"`). This is strictly more informative than falling back to all three frames when ambiguous.

4. **No chained frame propagation.** If a second upstream novel junction shifts the reading frame before reaching the junction of interest, the propagated frame could in principle be computed as `(phase_B + L) % 3`. However, phasing two junctions to the same transcript is impossible from short-read RNA-seq without transcript assembly — left for a future improvement.

5. **Sense-strand translation only.** For strand-specific RNA-seq (dUTP / first-strand), `bedtools getfasta -s` already reverse-complements minus-strand contigs so they represent the correct transcript sequence. Antisense translation of a strand-corrected contig has no established biological basis. Confirmed strand-specific for patient_001 (KAPA RNA HyperPrep Kit with RiboErase HMR, dUTP second-strand marking). **patient_002 strandedness not yet verified.**

**Frame offset formula:**
```
phase_at_donor = (exon_length − gtf_frame) % 3
frame_offset   = (−phase_at_donor) % 3          → 0, 1, or 2
```
`upstream_nt` is always divisible by 3 by pipeline design, so `frame_offset` equals the start position within the upstream codon context (0 = in-frame, 1 = +1 shift, 2 = +2 shift). Donor coordinate: `junction.start` for + strand, `junction.end` for − strand.

**Changes:**
- `workflow/scripts/filter_junctions.py`: new `_build_cds_donor_lookup(gtf_path)` + `reading_frame` column in output TSV; `classify_junctions()` accepts optional `gencode_gtf` arg; `--gencode-gtf` CLI flag added
- `workflow/rules/filter_junctions.smk`: `gencode_gtf` added as explicit input to `filter_junctions` rule
- `research/manuscript/METHODS.md` §5: updated to describe reading frame annotation subsection
- `research/manuscript/DISCUSSIONS.md`: new section "Reading frame annotation: why translation is not restricted to the canonical frame"
- `workflow/tests/test_filter_junctions.py`: 13 new tests (9 × `_build_cds_donor_lookup`, 4 × `reading_frame` column integration); all 156 tests pass

### ~13:30 UTC

#### Issue #16 — Pre-built HISAT2 index (PR #111)

**Goal:** Replace the 60–90 min `hisat2-build` step with a download of the official pre-built `grch38_tran` index (~10–15 min).

**Changes:**
- Added `alignment.hisat2_prebuilt_url` config key; production config points to `grch38_tran.tar.gz` on the HISAT2 S3 mirror; test config sets it empty to keep building from the chr22 FASTA
- Replaced `hisat2_index` rule with conditional `hisat2_download_index` (URL set) / `hisat2_index` (URL empty) branching; introduced `_HISAT2_INDEX_PREFIX` (`genome_tran` vs `genome`) used by `hisat2_align`
- `grch38_tran` includes GENCODE splice sites baked in — strictly better than the plain genome index previously built from scratch

**Validation:**
- Tarball URL confirmed reachable; structure verified locally (`tar -tz`): top-level `grch38_tran/` with `genome_tran.{1..8}.ht2`
- Snakemake dry runs confirmed correct rule selection for both empty and non-empty URL
- Cloud extraction test on `hisat2-index-test` (e2-standard-4, europe-west1-b): all 8 `.ht2` files extracted correctly via `--strip-components=1`
- 143 unit tests pass

**Post-review fixes:** wrapped `curl | tar` in subshell for correct log capture; added `--retry 3`; added `resources: mem_mb=1000`; added re-download comment in config.

### ~12:52 UTC

#### Known limitation — contig assembly uses reference genome, not patient reads

`assemble_contigs.smk` uses `bedtools getfasta` on GRCh38 to extract flanking sequences around novel junction coordinates. It does not assemble from patient RNA-seq reads.

**Implication:** somatic SNVs or indels in the exon flanks are ignored. Predicted peptide sequences assume wild-type exonic context around the junction.

**Why acceptable for now:** the neoepitope signal is primarily from the novel exon–exon junction itself. Flank mutations are a second-order effect; reference extraction is standard in comparable published pipelines.

**Future improvement:** in hypermutated tumors (MSI-high, POLE-mutant), flank mutations could meaningfully alter predicted peptides. A future direction would be to extract junction-spanning reads from the BAM and assemble the local haplotype directly.

### ~12:41 UTC

#### Issue #107 — Single-VM consolidation merged (PR #109)

**Goal:** Eliminate the two-VM / three-phase architecture (CPU alignment VM + GPU TCRdock VM with GCS handoff) in favour of a single consolidated VM.

**Changes:**
- Deleted `scripts/setup_cloud.sh` and `scripts/setup_tcrdock_vm.sh`; replaced with unified `scripts/setup_vm.sh` (8-step idempotent setup: deps → GPU check → Docker + NVIDIA Container Toolkit → Miniforge3 → snakemake env → repo clone → TCRdock image build → reference data)
- Rewrote `scripts/run_cloud_gpu.sh` to manage a single `neoepitope-pipeline` VM (n1-highmem-8 + P100, 200 GB); Snakemake runs alignment → MHCflurry (GPU) → TCRdock (Docker/GPU) → report in one session with no GCS handoff
- Renamed `config/tcrdock_gpu.yaml` → `config/gpu.yaml` (overlay now enables MHCflurry GPU acceleration in addition to TCRdock); introduced `GPU_CONFIG_FILE` variable in `run_cloud_gpu.sh`; factored out mode-independent settings (`RESULTS_DIR`, `GCS_PATH`) from `case` block
- Fixed `--conda-cleanup-envs` call missing the GPU overlay; added `nvidia-container-toolkit` apt install fallback; forwarded `--keep-vm` through detach re-invocation

**Validated:** end-to-end cloud run on `neoepitope-pipeline` VM (2026-04-23, alignment → MHCflurry 7.72M predictions → TCRdock → HTML report).

### 09:40 UTC

#### Issue #105 — GPU MHCflurry validation (patient_001)

**Goal:** Validate that MHCflurry runs correctly on the P100 GPU and measure the speedup vs CPU baseline.

**Run:** patient_001 (SRR37781424, Luminal A breast cancer) on `neoepitope-pipeline` VM (n1-highmem-8 + P100, europe-west1-b).

| Timestamp (UTC) | Event |
|---|---|
| 09:16:25 | Snakemake started |
| 09:18:42 | Conda env activated, alleles loaded |
| 09:18:46 | GPU detected — sequential execution path selected |
| 09:18:54 | HLA-A\*31:01 started |
| 09:21:40 | HLA-A\*26:01 started (2m46s) |
| 09:24:24 | HLA-B\*18:01 started (2m44s) |
| 09:27:09 | HLA-B\*15:63 started (2m45s) |
| 09:29:53 | HLA-C\*07:01 started (2m44s) |
| 09:32:37 | HLA-C\*03:03 started (2m44s) |
| 09:36:15 | All 6 alleles complete (3m38s for last) |
| 09:36:33 | Report generated, pipeline done |

**Results:**

- Total predictions: 7,718,952 (1,260,074 unique peptides × 6 alleles)
- Strong binders (IC50 ≤ 50 nM): 15,880
- Weak binders (IC50 ≤ 500 nM): 149,662
- **Total MHCflurry time: ~17m21s vs CPU baseline ~1.5 hours → ~5.2× speedup on P100**

**Why not 10–50× as estimated?** Sequential allele execution (one CUDA context, one model in VRAM) — can't parallelise alleles because each would require a separate model copy in the 16 GB P100 VRAM. GPU parallelism applies within each allele's `predict_to_dataframe()` call (all 1.26M peptides at once), not across alleles.

**Infrastructure issues resolved during development (branch `feat/issue-105-gpu-mhcflurry`):**

1. **`hisat2.yaml` samtools/libdeflate conflict:** `regtools >=1.0.0` requires `libdeflate >=1.26`; no bioconda linux-64 samtools/htslib build is compiled against it. Removed samtools from conda env; pipeline uses system apt samtools 1.13 instead.
2. **`ProcessPoolExecutor` OOM on CPU:** 6 workers × full model (~8 GB RAM each) = ~48 GB on a 52 GB VM → BrokenProcessPool. Fixed by removing the process pool entirely — sequential execution on both CPU and GPU paths.
3. **NVIDIA driver 580 incompatible with P100:** Driver 580 requires GSP firmware; P100 (Pascal, SM 6.0) lacks it. Fixed by in-place downgrade to driver 570 (`apt-get install nvidia-driver-570-server && purge *580* && reboot`) — no VM deletion needed.
4. **PyTorch SM 6.0 mismatch:** PyTorch 2.5+ dropped SM 6.0 (P100/Pascal) support. Pinned `torch>=2.0,<2.5` in `python.yaml` (installs 2.4.1). Also rewrote `_has_gpu()` to use a PyTorch smoke-test kernel instead of TensorFlow — TF reported GPU available even when PyTorch kernels would fail on SM mismatch.
5. **Orphan sentinel file:** `.snakemake/conda/080a2daa..._.env_setup_done` existed without the actual env directory → Snakemake skipped rebuild → `ModuleNotFoundError: No module named 'mhcflurry'`. Fixed by deleting the orphan sentinel.

---

## 2026-04-22

### Issue #105 — GPU MHCflurry: CUDA architecture notes

**Host driver vs CUDA toolkit — how they relate:**

The **host NVIDIA driver** lives on the VM OS and is the only software that directly controls the GPU hardware. Its version determines the highest CUDA toolkit version it can support (e.g. driver 570 → CUDA 12.8).

The **CUDA toolkit** (cuDNN, cuBLAS, etc.) is what a specific application was compiled against. It lives wherever the application lives — inside a Docker container, inside a pip package — completely independent of the host.

**NVIDIA's backward-compatibility rule:** a host driver supports all CUDA toolkit versions ≤ its own. Driver 570 supports CUDA 12.8, 12.0, 11.8, 10.x, etc.

In our pipeline, two applications run on the same VM with different CUDA requirements:

```
Host (Deep Learning VM image)
  └─ NVIDIA driver 570 (CUDA 12.8 capable)
       ├─ conda env: python.yaml
       │    └─ tensorflow[and-cuda] pip package → bundles CUDA 12.x libs → MHCflurry
       └─ Docker container: tcrdock:latest → bundles CUDA 11.8 libs → TCRdock
```

Neither application depends on the host CUDA toolkit — each carries its own. The host driver just needs to be ≥ the highest toolkit version in use. 12.8 covers both.

**Why `ProcessPoolExecutor` crashes on GPU:**

MHCflurry's CPU path spawns multiple worker processes (one per allele), each importing TensorFlow and initialising its own CUDA context. Multiple CUDA contexts on the same GPU conflict → `BrokenProcessPool` crash. Fix: detect GPU at startup, load the predictor once in the main process, and iterate over alleles sequentially. The GPU handles massive internal parallelism within each `predict_to_dataframe()` call — sequential allele iteration adds negligible overhead.

---

### Issue #98 — Proteome k-mer filter (PR #106, merged)

Rewrote `blastp_filter.py` as `proteome_filter.py`. Instead of running blastp as a subprocess, the script parses the Swiss-Prot FASTA once to build a `set[str]` of all k-mers from every canonical human protein, then checks each query peptide via O(1) set lookup. Drops the BLAST conda env entirely — runtime goes from ~2 hours (424K peptides × 4 threads) to seconds.

The `matched_accessions` column in `peptides_excluded.tsv` carries all source proteins for each excluded peptide (semicolon-separated), replacing the separate `blastp_hits.tsv` audit file. Rule, script, test, and config key all renamed `blastp_filter` → `proteome_filter`. 19 unit tests passing.

**Local test results (chr22, patient_001_test, 8/9/10-mers):**

| | Count | % |
|---|---|---|
| Total peptides into filter | 4,163 | 100% |
| Novel (passed) | 4,119 | 98.9% |
| Excluded (self-peptides) | 44 | 1.1% |

---

### PR #101 — Move research artefacts to `research/` (Issue #100)

Separated operational documentation from research outputs. `docs/lab_notebook.md` → `research/lab_notebook.md`; `docs/manuscript/` → `research/manuscript/`. `docs/` now contains only software/infrastructure documentation (installation, configuration, cloud guide, etc.).

---

### PR #102 — Remove automatic Claude PR review CI workflow

Removed `.github/workflows/claude-code-review.yml`, which triggered a full Claude review on every push. In a solo project this was noisy and consumed usage unnecessarily. On-demand review is still available via `@claude` mentions through `claude.yml`.

---

### Issue #82 — Multi-length peptide extraction (8-mer, 9-mer, 10-mer)

**Motivation:** MHC-I presents 8–10-mer peptides. Extracting only 9-mers misses a significant fraction of true binders.

**Design:**
- Replaced hardcoded 9-mer extraction in `translate_peptides.py` with a loop over a configurable `translation.peptide_lengths: [8, 9, 10]` list.
- Contig flank is now auto-derived: `flank_nt = 3 * (max_length - 1)`. For max length 10 this gives symmetric 27/27 nt flanks (54 nt contigs), up from 26/24 (50 nt). The old 50 nt contigs were 3 nt too short to cover all 10-mer junction-spanning windows on the downstream side.
- MHCflurry natively accepts 8–15-mers; blastp operates on whatever peptides it receives. No changes needed in either step.
- 16 unit tests updated/added covering 8/9/10-mer start ranges and known-start sequences.

**Local test results (chr22, patient_001_test):**

| | 9-mer only | 8/9/10-mer |
|---|---|---|
| Unique peptides into MHCflurry | 1,339 | 4,045 (~3×) |
| Total predictions (6 alleles) | 8,172 | 24,714 (~3×) |
| Strong binders (≤50 nM) | 39 | 37 |
| Weak binders (≤500 nM) | 290 | 502 (+73%) |

**Follow-on issues filed:**
- **Issue #97** — Use GENCODE CDS reading frame to restrict translation to the biologically relevant frame. Currently all 3 reading frames are used per junction; the correct frame can be determined from the GENCODE GTF for junctions where (a) the upstream donor maps unambiguously to a single CDS exon end, and (b) no other tumor-exclusive junction is detected upstream in the same gene. Both conditions must hold; overlap between genes and upstream frame-shifting junctions are edge cases requiring fallback to all 3 frames.
- **Issue #98** — Replace blastp exact-match filter with a Python set-based k-mer lookup. On the production patient (424K unique peptides), blastp ran for >2 hours. A set of all k-mers from the Swiss-Prot FASTA built once at startup reduces this to seconds via O(1) lookup per peptide. Also supports near-self filtering in future (generate all 1-mismatch variants of query, check against exact set).

---

### Issue #89 — blastp filter against human reference proteome (PR #96, merged)

**Overview:** Added `blastp_filter_peptides` rule between `translate_peptides` and `run_mhcflurry`. Excludes peptides with a 100% identity / 100% query coverage match in UniProt Swiss-Prot (self-peptides the immune system is tolerized to).

**Cloud run debugging:**
- *Run 1 — blastp silently skipped:* `mhc_affinity.tsv` existed from a pre-blastp run. With `--rerun-triggers mtime`, Snakemake never propagated upward to check for the missing `peptides_novel.tsv`. Fixed by deleting downstream outputs before re-running.
- *Run 2 — `Error: Unknown argument: "perc_identity"`:* `-perc_identity` is a `blastn`-only flag; blastp 2.16.0 rejects it. Fixed: removed flag, now uses `-outfmt "6 qseqid sseqid pident"` + `-qcov_hsp_perc 100` (full coverage at search time) + Python `int(float(pident)) == 100` filtering. Added 15 unit tests.

**Local test (9-mers):** 1,371 peptides → 9 self-peptides excluded (0.7%) → 1,362 novel passed to MHCflurry.

---

### PR #93 — Remove `gcs:` block and BAM upload fix

Removed the `gcs:` block from `config/config.yaml` and `config/test_config.yaml` (BAM upload had been moved to `scripts/run_cloud_cpu.sh` and no longer read config). Also removed `temp()` + `upload_bam` from `alignment.smk` which was preventing BAM/BAI/BED files from reaching GCS.

---

### PR #94 — Snakefile refactor (Issue #92)

Reduced the Snakefile to a minimal entry point: patient ID derivation and shared config moved into `common.smk`. Added a zero-rows guard in `_read_samples_tsv` that raises a clear `ValueError` on empty or missing TSV rather than a bare `IndexError`.

---

### PR #95 — Claude GitHub Actions workflows

Added Claude Code Review and Claude PR Assistant GitHub Actions workflows (`.github/workflows/`).

---

## 2026-04-21

### Patient_001 (gastric cancer) — GPU run completed

Patient_001 TCRdock run completed successfully overnight. Top candidate EVAEYNASF / HLA-A\*26:01 (IC50 = 16.5 nM) was run through TCRdock on the P100 GPU VM. Outputs archived to `gs://splice-neoepitope-project/results/patient_001/`: `docking_scores.tsv`, `top_candidate.pdb`, `report.html`. Both CPU and GPU VMs stopped cleanly (TERMINATED).

This is the second confirmed end-to-end run for patient_001 (first was on the old `tcrdock-handoff` bucket; this run used the new `splice-neoepitope-project` bucket with the updated `run_cloud_gpu.sh`).

### Infrastructure fix — `.snakemake/metadata` sync via GCS (commit `3ddfe96`)

**Problem:** The GPU VM was re-running all CPU pipeline steps (alignment, filtering, MHCflurry) on every invocation, even when results were already present. Root cause: without `.snakemake/metadata/` on the GPU VM, Snakemake has no baseline for `--rerun-triggers code params` and triggers re-runs of every rule.

**Fix:** Extended `run_cloud_gpu.sh` to sync `.snakemake/metadata/` from the CPU VM to GCS (Phase 2 upload) and then from GCS to the GPU VM (Phase 3 download), alongside `results/` and `logs/`. Only `.snakemake/metadata/` is synced — the `conda/` subdir is gigabytes and is not transferred.

**Related investigation — "Incomplete files" warning:**
The GPU VM also showed an `IncompleteFilesException` on startup despite all result files being present. Root cause: stale `.snakemake/incomplete/` markers from a previous Spot VM preemption. These are separate from `.snakemake/metadata/` — incomplete markers are only cleared by `--rerun-incomplete`, `--cleanup-metadata`, or manual removal. Resolved by manually running `rm -rf .snakemake/incomplete/` on the GPU VM as a one-time fix. Decision not to add this to the script: (a) the VM is now STANDARD (not Spot), so future preemptions are unlikely; (b) clearing incomplete markers blindly would mask legitimately partial TCRdock outputs from future interrupted runs, since `gcloud storage cp` does not delete local files absent from GCS.

### Issue #73 — patient_002 manuscript results (PR #74, merged)

Added patient_002 osteosarcoma results to `research/manuscript/`:
- **RESULTS.md:** full patient_002 section — dataset, HLA typing + serology validation, junction funnel (347,046 raw → 55,912 tumor-exclusive), peptide translation (781,424 9-mers), MHC predictions (12,430 strong binders), top candidate TTDPVQALY / HLA-A\*01:01 (IC50 = 23.9 nM), TCRdock caveat (fallback allele bug).
- **CONCLUSIONS.md:** added patient_002 Key Findings (points 5–8); updated Limitations HLA section (confirmed match, no longer future tense); updated Future Directions (T0 complete, focus on T1/T2).
- **DISCUSSIONS.md:** updated "Impact of missing matched normal: patient_002" section with real WES proxy counts (106,474 apparent junctions, only 3 overlap tumor).

PR #74 (`docs/issue-73-manuscript-patient002` → `main`) opened and merged.

### Issue #59 — Known HLA serology input via samples TSV (PR #69, merged)

**Motivation:** Clinical HLA serotyping (Red Cross / WGS) is more reliable than OptiType from RNA-seq for germline alleles. Patient_002 has confirmed Red Cross serology (A\*01:01/01:11N, B\*08:01/27:05, C\*01:02/07:01) which should take priority over the blood WES OptiType call when no tumor call is available.

**Design decision:** Following nf-core convention (samplesheet as single source of truth), serology alleles are stored as six inline columns in each patient's samples TSV (`serology_A1/A2`, `serology_B1/B2`, `serology_C1/C2`). The tumor row is left empty; only the normal/blood row carries the germline values.

**Priority order** implemented in `aggregate_hla_alleles.py`:
1. **Tumor OptiType** — preferred because HLA LOH can occur in the tumor, altering what is actually presented
2. **Serology** — germline ground truth from clinical typing
3. **Normal/blood OptiType**
4. **Config fallback**

**Null allele handling:** A\*01:11N is a null allele (not expressed at the cell surface). It is excluded from the prediction allele list but retained in `hla_qc.tsv` for transparency. The serology validation result is `match (null allele excluded)` when OptiType finds only the expressed allele.

**QC output extended:** `hla_qc.tsv` now includes `serology_allele1`, `serology_allele2`, and `serology_validation` columns. The HTML report's HLA section shows a "Known alleles / validation" column when serology data is present.

---

### Issue #78 / PR #77 — `gcloud storage rsync` + `mtime` rerun trigger (merged)

**Problem:** `gcloud storage cp` fresh-timestamps all downloaded files. With `--rerun-triggers code params` on the GPU VM, this is fine as long as only code/param changes matter. But if `mhc_affinity.tsv` changes content (e.g. HLA alleles change), TCRdock would not re-run under `code params` alone because TCRdock's script and params are unchanged.

**Fix:** Switched Phase 3 GPU download from `gcloud storage cp` to `gcloud storage rsync --recursive` (CRC32C checksum comparison, only downloads changed files, preserves mtime of unchanged files). Also added `mtime` to GPU Snakemake's `--rerun-triggers` so file content changes cascade correctly. Metadata download remains `cp` (always overwrite with CPU VM's authoritative `.snakemake/metadata/`).

---

### Issue #57 Phase 1 — `report.tsv` structured summary artifact (PR #80 merged, PR #81 open)

Added `_build_report_tsv()` to `generate_report.py`. Emitted alongside `report.html` whenever `output_tsv` is provided. Schema: `patient_id | stage | metric | value | notes`. Stages:

| Stage | Metrics written |
|---|---|
| `junction_filtering` | `unannotated`, `tumor_exclusive`, `normal_shared` per sample |
| `mhc_prediction` | `total_predictions`, `strong`, `weak`, `non` binder counts |
| `top_candidate` | `peptide`, `allele`, `ic50_nM`, `binder_class` |
| `hla_typing` | `HLA-A/B/C` alleles with source and read count |
| `tcrdock` | `pdb_available` true/false |

`report_tsv` added as a declared Snakemake output in both `analysis.smk` (non-GPU path) and `structure.smk` (GPU path with TCRdock). The `structure.smk` output declaration was missing from PR #80, causing an `AttributeError` on the GPU VM. Fixed in PR #81.

Phase 2 (refactor `report.html` to read from `report.tsv` instead of recomputing) deferred to Issue #79.

---

### Patient_001 (gastric cancer) — re-run with tumor-first HLA (successful)

Re-ran patient_001 with all merged changes (tumor-first HLA priority, serology columns, rsync, report.tsv). TCRdock re-ran despite existing results because `aggregate_hla_alleles.py` code changed → HLA re-typed → `mhc_affinity.tsv` regenerated with fresh mtime → `mtime` trigger correctly fired TCRdock (~4 min GPU time).

**HLA-B alleles updated** (tumor-first policy now active): B\*15:01/B\*18:02 → **B\*15:63/B\*18:01**. The A alleles were unchanged so the top candidate is the same.

**`report.tsv` read directly from GCS** (`gcloud storage cat gs://splice-neoepitope-project/results/patient_001/reports/report.tsv`):

| Stage | Metric | Value |
|---|---|---|
| junction_filtering | tumor_exclusive | 27,348 |
| junction_filtering | normal_shared | 2,681 |
| mhc_prediction | total_predictions | 2,598,882 |
| mhc_prediction | strong | 13,139 |
| top_candidate | peptide | EVAEYNASF |
| top_candidate | allele | HLA-A\*26:01 |
| top_candidate | ic50_nM | 16.51 |
| hla_typing | HLA-B | HLA-B\*18:01 / HLA-B\*15:63 |

patient_002 re-run planned next (with all latest changes including PR #81).

---

### Documentation update (PR #87)

Updated three docs to reflect recent pipeline changes:
- `docs/data_preparation.md` — serology columns added to sample manifest format and column table; HLA typing roles table updated with tumor-first priority order
- `docs/configuration.md` — HLA section: added allele priority order note and serology column cross-reference
- `docs/google_cloud_guide.md` — fixed stale bucket name (`<PROJECT_ID>-tcrdock-handoff` → `splice-neoepitope-project`); added `{patient_id}` to result paths; added `report.tsv` to retrieval example; noted `gcloud storage rsync` in "How it works"

---

## 2026-04-20

### Patient_002 (osteosarcoma BG003082) — first full production run

**Patient:** BG003082 T0 tumor (Boston Gene, Nov 2022, paired-end RNA-seq ~10 GB) + BG003082 N0 WES normal (blood-derived, used for HLA typing only).

**HLA typing:** A\*01:01/A\*01:01, B\*08:01/B\*27:05, C\*07:01/C\*01:02 — confirmed match to Red Cross serology (A\*01:01/01:11N, B\*08:01/27:05, C\*01:02/07:01). First patient with ground-truth HLA validation.

**Results:** Run completed end-to-end: alignment → HLA typing → MHCflurry → TCRdock → HTML report with Mol\* 3D viewer. Final outputs archived to `gs://splice-neoepitope-project/results/patient_002/`.

**Infrastructure bugs discovered and fixed (PR #69, branch `feat/issue-65-patient002-cloud-run`):**

- **mtime re-run cascade:** Re-downloaded temp FASTQs had newer mtime than existing `junctions.tsv`, causing unnecessary re-alignment and OptiType re-runs. Fixed with `ancient()` on FASTQ inputs in both `hisat2_align` and `run_optitype`.
- **OptiType OOM:** razers3 peaks at ~36 GB on full RNA-seq FASTQs. CPU VM upgraded from n1-standard-8 → n1-highmem-8 (52 GB) → n2-highmem-8 (64 GB; n1 unavailable in zone). OptiType threads capped at 5 to force sequential sample execution and prevent concurrent OOM.
- **samtools sort OOM:** `-m 3G` caused OOM on the WES normal sample (8 threads × 3 GB). Reduced to `-m 1G`.
- **Boot disk:** 50 GB → 100 GB pd-ssd to handle reference index + paired-end FASTQ staging.
- **`--rerun-incomplete`:** Added to orchestrator snakemake invocation so killed runs resume cleanly instead of raising `IncompleteFilesException`.
- **MHCflurry re-run on GPU VM:** `resources/mhcflurry_models.done` sentinel absent on GPU VM → cascade triggered `run_mhcflurry` re-run. Fixed by running `snakemake resources/mhcflurry_models.done --use-conda` before the main TCRdock run (uses the correct `python.yaml` env; direct `mhcflurry-downloads fetch` in the `snakemake` bootstrap env would fail on a fresh VM).
- **GPU VM GCS upload permission:** Existing TERMINATED GPU VM lacked `--scopes=cloud-platform`. `gcloud storage cp` failed with permission denied; TCRdock results were not uploaded. Fixed by adding `gcloud compute instances set-service-account ... --scopes=cloud-platform` before `instances start` in orchestrator.
- **tmux not installed on GPU VM:** Deep Learning VM image does not include tmux by default. Added idempotent `apt-get install -y -q tmux` to GPU provisioning block.
- **VM auto-stop:** CPU and GPU VMs now unconditionally stop on pipeline exit (success or failure).

**TCRdock result:** pLDDT 92.25 for top candidate `FMSGFLYFV` on `HLA-A*02:01` (fallback allele — note: for patient_002 the actual alleles are A\*01:01, not A\*02:01; fallback was used due to a sentinel issue, now fixed for future runs).

---

### Patient_001 (gastric cancer) — cloud run started

Updated `config/samples/patient_001.tsv` to use ENA HTTPS URLs instead of local `data/` paths, enabling cloud runs without pre-staging FASTQs. Run started; `n2-highmem-8` required after `n1-highmem-8` hit `ZONE_RESOURCE_POOL_EXHAUSTED` in `europe-west1-b`.

---

### Documentation update

README slimmed from ~600 to 337 lines. Detailed content moved to new dedicated docs:
- `docs/installation.md` — full conda/Snakemake setup
- `docs/data_preparation.md` — aligner selection, FASTQ sources, manifest format
- `docs/configuration.md` — full `config.yaml` parameter reference
- `docs/google_cloud_guide.md` — added manual TCRdock run section

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
