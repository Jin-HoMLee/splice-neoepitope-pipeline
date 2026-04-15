# Splice Neoepitope Pipeline

A modernised, reproducible reimplementation of the 2015 research pipeline for
predicting cancer neoepitopes arising from alternative splicing detected by
RNA-Seq.

> **Original work**: Jin-Ho Lee, Seoul National University, 2015.
> *"Identification of Cancer-Specific Neoepitopes Arising from Alternative
> Splicing Detected by RNA-Seq."*

---

## Table of Contents

1. [Scientific Background](#scientific-background)
2. [Pipeline Overview](#pipeline-overview)
3. [Installation and Setup](#installation-and-setup)
4. [Data Source Options](#data-source-options)
   - [Option A: Local Alignment (Recommended)](#option-a-local-alignment-recommended---no-institutional-access-required)
   - [Option B: Cloud-Based Execution (Google Cloud)](#option-b-cloud-based-alignment-for-users-without-local-resources)
   - [Option C: GDC Download](#option-c-gdc-download-requires-institutional-access)
5. [MHCflurry](#mhcflurry-epitope-predictor)
6. [TCRdock Structural Validation (Optional)](#tcrdock-structural-validation-optional)
7. [Reference Data](#reference-data)
8. [Local Testing (chr22 subset)](#local-testing-chr22-subset)
9. [Running the Pipeline](#running-the-pipeline)
10. [Configuration](#configuration)
11. [Output Description](#output-description)
12. [Modernisation Changelog](#modernisation-changelog)
13. [Project Structure](#project-structure)
14. [Citation](#citation)
15. [Further Reading](#further-reading)

---

## Scientific Background

Tumors frequently exhibit splice junctions absent from matched normal tissue.
These tumor-specific junctions can produce novel peptide sequences — **neoepitopes** —
that are recognisable as foreign by the immune system and are candidate targets for
cancer immunotherapy.

This pipeline identifies those junctions from RNA-seq data, filters them against the
matched normal sample, and predicts which resulting peptides bind MHC class I molecules.

> For full biological background and study design, see [`docs/INTRODUCTION.md`](docs/INTRODUCTION.md).

---

## Pipeline Overview

```
RNA-Seq data (local FASTQ files or GDC API)
        │
        ▼ Step 1: Align/Download
  Splice junction quantification files (.tsv)
        │
        ▼ Step 2: Classify by origin
  - Remove annotated junctions (GENCODE reference)
  - Compare tumor vs. matched normal:
      patient_specific (in normal) → excluded
      tumor_specific   (not in normal) → keep
        │
        ▼ Step 3: Assemble
  50 nt contigs (26 nt upstream + 24 nt downstream, bedtools getfasta)
        │
        ▼ Step 4: Translate
  16-mer peptides (3 reading frames, truncated at stop codons)
        │
        ▼ Step 5: Predict
  Junction-spanning 9-mer filter → MHCflurry 2.x → IC50 binding affinities (HLA-A*02:01)
        │
        ▼ Step 6: TCRdock structural validation (optional, GPU)
  TCR-peptide-MHC ternary complex 3D structure (AlphaFold v2 backend)
        │
        ▼ Step 7: Report
  Junction origin summary + top binders + Mol* 3D viewer HTML report
```

---

## Installation and Setup

### 1. System Requirements

| Requirement | Minimum | Notes |
|-------------|---------|-------|
| OS | Linux (x86-64) or macOS | Windows is not supported |
| CPU | 4 cores | More cores speed up parallel steps |
| RAM | **8 GB** | Using HISAT2 aligner (or 32 GB for STAR) |
| Disk | 50 GB free | Reference genome + data files |
| Python | 3.11+ | Managed automatically via conda |
| Git | any | For cloning this repository |

> **Note**: The default aligner is now HISAT2, which requires only ~8 GB RAM.
> If you have a high-memory system (32+ GB), you can optionally use STAR for
> maximum accuracy.

---

### 2. Install Conda (Miniforge — recommended)

Miniforge is a minimal Conda installer that defaults to the `conda-forge`
channel and ships with the faster `mamba` solver.

```bash
# Download the installer for your OS
# Linux:
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
  -o Miniforge3.sh
# macOS (Apple Silicon):
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh \
  -o Miniforge3.sh
# macOS (Intel):
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh \
  -o Miniforge3.sh

# Run the installer (accept the licence; let it initialise your shell)
bash Miniforge3.sh -b -p "$HOME/miniforge3"
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda init bash   # or: conda init zsh

# Reload your shell, then verify
conda --version   # should print: conda 24.x.x or similar
```

> **Already have Anaconda or Miniconda?**  That works too.
> Make sure `conda-forge` and `bioconda` channels are added:
> ```bash
> conda config --add channels bioconda
> conda config --add channels conda-forge
> conda config --set channel_priority strict
> ```

---

### 3. Install Snakemake

Create a dedicated conda environment that contains Snakemake 7.x and the
`snakemake-executor-plugin-cluster-generic` plugin (needed for cluster runs):

```bash
conda create -n snakemake -c conda-forge -c bioconda \
  "snakemake>=7.0,<9" \
  snakemake-executor-plugin-cluster-generic \
  python=3.11 \
  -y

conda activate snakemake

# Verify
snakemake --version   # expected output: 7.x.x or 8.x.x
```

> The pipeline conda environments (`workflow/envs/python.yaml`,
> `workflow/envs/biotools.yaml`) are created automatically the first time
> Snakemake runs each rule — you do **not** need to install biopython,
> bedtools, etc. manually.

---

### 4. Clone the Repository

```bash
git clone https://github.com/Jin-HoMLee/splice-neoepitope-pipeline.git
cd splice-neoepitope-pipeline
```

Verify the expected layout is in place:

```bash
ls -1
# Expected output:
# LICENSE  README.md  Snakefile  config/  docs/  resources/  workflow/
```

---

### 5. Quick Sanity Check

With the `snakemake` environment active, confirm Snakemake can parse the
workflow without errors (make sure to download the `gencode.v47.annotation.gtf.gz` and `GRCh38.primary_assembly.genome.fa` files according to [/resources/README.md](/resources/README.md) before running):

```bash
conda activate snakemake
snakemake --cores 1 --use-conda -n 2>&1 | head -20
```

You should see a list of jobs (or a note that all outputs are up to date).
If Snakemake prints a Python traceback, check that the `snakemake` conda
environment is active and that all files were cloned correctly.

---

## Data Source Options

This pipeline supports two data source modes:

| Mode | Institutional Access | Description |
|------|---------------------|-------------|
| **Local Alignment** | ❌ Not required | Align your own FASTQ files using HISAT2 or STAR |
| **GDC Download** | ✅ Required | Download pre-computed files from GDC |

---

### Option A: Local Alignment (Recommended) — No Institutional Access Required

This mode generates splice junction quantification from raw FASTQ files using
either **STAR** or **HISAT2**. You can use:

- Your own RNA-Seq data
- Publicly available datasets from SRA, GEO, or ENCODE
- Any GRCh38/hg38-aligned RNA-Seq data

#### Choose Your Aligner

| Aligner | RAM Required | Index Size | Best For |
|---------|--------------|------------|----------|
| **HISAT2** | ~8 GB | ~8 GB | Laptops, small servers, limited resources |
| **STAR** | ~32 GB | ~30 GB | Full accuracy, high-memory systems |

> **Recommendation**: Start with **HISAT2** if you have limited RAM. HISAT2
> produces compatible splice junction output with ~8 GB RAM vs STAR's ~32 GB.

#### Step 1: Set Data Source Mode and Aligner

Edit `config/config.yaml`:

```yaml
data_source: "fastq"   # Align your own FASTQ files

samples_tsv: "config/samples.tsv"

alignment:
  aligner: "hisat2"    # "hisat2" (8 GB RAM) or "star" (32 GB RAM)
```

#### Step 2: Obtain RNA-Seq FASTQ Data

**Public RNA-Seq datasets (no access restrictions):**

| Source | URL | Description |
|--------|-----|-------------|
| **SRA/ENA** | https://www.ncbi.nlm.nih.gov/sra | Millions of publicly available RNA-Seq runs |
| **GEO** | https://www.ncbi.nlm.nih.gov/geo | Gene Expression Omnibus |
| **ENCODE** | https://www.encodeproject.org | High-quality RNA-Seq from cell lines |
| **GTEx** | https://gtexportal.org | Normal tissue RNA-Seq (open access) |

**Downloading public RNA-Seq data:**

ENA (European Nucleotide Archive) provides direct HTTPS access to FASTQ files
and is the most reliable option, especially on macOS:

```bash
# Stream from ENA — no extra tools needed, works on all platforms
# Replace ERR188273 with your accession; adjust the FTP path accordingly
curl -L "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188273/ERR188273_1.fastq.gz" \
    -o data/sample_R1.fastq.gz
curl -L "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188273/ERR188273_2.fastq.gz" \
    -o data/sample_R2.fastq.gz
```

Alternatively, use `fasterq-dump` from sra-tools if you already have it
installed (note: conda installation of sra-tools can be unreliable on macOS
arm64 — ENA download is preferred in that case):

```bash
fasterq-dump SRR12345678 --split-files --outdir data/
```

> **Tip**: Search GEO for `"RNA-Seq" AND "cancer" AND "Homo sapiens"` to find
> publicly available cancer RNA-Seq datasets. Check ENA for direct FASTQ URLs.

#### Step 3: Create Sample Manifest

Create a tab-separated value (TSV) file at `config/samples.tsv` with **no comments or headers**:

```tsv
patient_id	sample_id	sample_type	fastq1	fastq2
patient_001	tumor_01	Primary Tumor	data/tumor_01_R1.fastq.gz	data/tumor_01_R2.fastq.gz
patient_001	normal_01	Solid Tissue Normal	data/normal_01_R1.fastq.gz	data/normal_01_R2.fastq.gz
patient_002	tumor_02	Primary Tumor	data/tumor_02_R1.fastq.gz	data/tumor_02_R2.fastq.gz
patient_002	normal_02	Solid Tissue Normal	data/normal_02_R1.fastq.gz
```

**Column descriptions:**

| Column | Required | Description |
|--------|----------|-------------|
| `patient_id` | Yes | Unique patient/experiment identifier. All rows with the same `patient_id` is processed together as a matched set (e.g. tumor + normal). Used as the output directory name. |
| `sample_id` | Yes | Unique sample identifier (e.g., SRR9143066, tumor_sample_1). Used in output filenames. |
| `sample_type` | Yes | Either `"Primary Tumor"` or `"Solid Tissue Normal"`— used for junction classification (tumor samples provide junctions; normal samples filter out patient-specific junctions). |
| `fastq1` | Yes | Path to FASTQ file (gzip-compressed or uncompressed). For paired-end, this is read 1. |
| `fastq2` | No | Path to read 2 FASTQ for paired-end sequencing. Leave empty for single-end data. |

**Important notes:**
- The file **must have exactly one header row** with column names in the order shown above.
- Do **not** include comment rows (lines starting with `#`) — comments belong in this README, not in the TSV.
- Each `patient_id` should match all related samples (tumor + normal, or multiple tumors from the same patient).
- Paths in `fastq1` and `fastq2` are interpreted relative to the pipeline working directory (usually where you run `snakemake`).

#### Step 4: Run the Pipeline

```bash
snakemake --cores 8 --use-conda
```

The pipeline will:
1. Build genome index (first run only: ~10 min for HISAT2, ~30 min for STAR)
2. Align each sample (~10-30 min per sample)
3. Continue with filtering, translation, and prediction

#### System Requirements by Aligner

**HISAT2 (Recommended for most users):**

| Resource | Minimum | Recommended | Notes |
|----------|---------|-------------|-------|
| RAM | **8 GB** | 16 GB | HISAT2 indexing needs only ~8 GB |
| Disk | 50 GB | 100 GB | Smaller index (~8 GB) |
| CPU | 4 cores | 8 cores | Alignment is parallelised |

**STAR (Full accuracy mode):**

| Resource | Minimum | Recommended | Notes |
|----------|---------|-------------|-------|
| RAM | 32 GB | 64 GB | STAR indexing needs ~32 GB |
| Disk | 100 GB | 200 GB | Larger index (~30 GB) |
| CPU | 4 cores | 16 cores | Alignment is parallelised |

---

### Option B: Cloud-Based Alignment (For users without local resources)

If you don't have a machine with sufficient RAM, you can run the **entire pipeline**
on cloud platforms. We provide a comprehensive guide for Google Cloud Platform.

| Platform | Cost | Description |
|----------|------|-------------|
| **Google Cloud** ⭐ | Pay-per-use | **[Full guide](docs/google_cloud_guide.md)** — Recommended |
| **Galaxy** | Free | [usegalaxy.org](https://usegalaxy.org) — Run HISAT2/STAR in browser |
| **AWS** | Pay-per-use | EC2 instances |

#### Google Cloud Platform (Recommended)

We provide a **[complete step-by-step guide](docs/google_cloud_guide.md)** for running
the pipeline on Google Cloud, including:

- **VM setup** — One-click instance creation with correct specs
- **Cost estimates** — ~$1–5 for small test runs, ~$10–30 for full analysis
- **Data transfer** — Upload/download FASTQ files and results
- **Spot VMs** — Save 60–80% with preemptible instances
- **Troubleshooting** — Memory, disk, and SSH issues

**Quick start (Google Cloud):**

```bash
# 1. Create a VM (16 GB RAM, sufficient for HISAT2)
gcloud compute instances create splice-pipeline \
    --zone=us-central1-a \
    --machine-type=n2-standard-4 \
    --boot-disk-size=100GB \
    --image-family=ubuntu-2204-lts \
    --image-project=ubuntu-os-cloud

# 2. Connect
gcloud compute ssh splice-pipeline --zone=us-central1-a

# 3. On the VM, install conda, clone repo, configure samples, then run:
snakemake --cores $(nproc) --use-conda 2>&1 | tee pipeline.log

# 4. Download results when done
gcloud compute scp --recurse splice-pipeline:~/splice-neoepitope-pipeline/results/ ./results/ --zone=us-central1-a

# 5. Delete VM to stop billing
gcloud compute instances delete splice-pipeline --zone=us-central1-a
```

📖 **See [docs/google_cloud_guide.md](docs/google_cloud_guide.md) for the full guide.**

#### Galaxy Workflow (Free, browser-based)

1. Upload your FASTQ files to [usegalaxy.org](https://usegalaxy.org)
2. Search for "HISAT2" or "STAR" in the tool panel
3. Run alignment with GRCh38/hg38 genome
4. Download the junction output file (`.tsv` or `.bed`)
5. Place in `results/raw_data/local/files/` and continue with the pipeline locally

---

### Option C: GDC Download (Requires Institutional Access)

This mode downloads pre-computed splice junction quantification files from the
GDC Data Portal. **TCGA data is controlled-access** and requires:

1. An eRA Commons account
2. dbGaP access approved for TCGA
3. A GDC authentication token

#### Step 1: Set Data Source Mode

Edit `config/config.yaml`:

```yaml
data_source: "gdc"   # Download from GDC (requires authentication)
```

#### Step 2: Obtain dbGaP Access

1. Apply for controlled-access TCGA data via [dbGaP](https://dbgap.ncbi.nlm.nih.gov/)
2. Your institution's signing official must approve the Data Access Request (DAR)
3. This process typically takes 1–4 weeks

> **Note**: If you already have TCGA access through another project (e.g.,
> via your institution's blanket approval), you can skip to Step 3.

#### Step 3: Download Your GDC Token

1. Go to **https://portal.gdc.cancer.gov/**
2. Click **Login** (top right) and authenticate via eRA Commons / NIH
3. After login, click your username (top right) → **Download Token**
4. Save the downloaded file (e.g., `gdc-user-token.txt`)

The token is valid for **30 days**; you'll need to download a new one after
it expires.

#### Step 4: Configure the Pipeline

Edit `config/config.yaml` and set the path to your token file:

```yaml
gdc:
  # ... other settings ...
  token_file: "~/.gdc-user-token.txt"   # path to your GDC token
```

Alternatively, place the token in your home directory:

```bash
mv ~/Downloads/gdc-user-token*.txt ~/.gdc-user-token.txt
chmod 600 ~/.gdc-user-token.txt   # restrict permissions
```

#### Troubleshooting GDC Downloads

| Error | Cause | Solution |
|-------|-------|----------|
| `403 Forbidden` | No token or expired | Download a fresh token from GDC |
| `403 Forbidden` | No dbGaP access | Apply for TCGA access via dbGaP |
| Token file not found | Wrong path | Check `gdc.token_file` in config.yaml |

---

## MHCflurry (Epitope Predictor)

This pipeline uses **MHCflurry 2.x**, an open-source, state-of-the-art MHC-I
binding predictor that does not require academic registration.  MHCflurry is
installed automatically via the conda environment — no manual setup needed.


### Reference

> O'Donnell TJ et al. (2020). MHCflurry 2.0: Improved Pan-Allele Prediction
> of MHC Class I-Presented Peptides by Incorporating Antigen Processing.
> *Cell Systems*, 11(1), 42-48.e7.

---

## TCRdock Structural Validation (Optional)

When enabled, the pipeline predicts the 3D structure of the TCR-peptide-MHC
ternary complex for the top neoepitope candidate using **TCRdock** (Bradley
et al.), a modified AlphaFold v2 multimer backend adapted for TCR:pMHC
complexes.

- **Requires**: Linux x86-64, NVIDIA GPU (T4 or better), Docker with NVIDIA Container Toolkit
- **Not compatible** with macOS or CPU-only machines
- **Disabled by default** — local / CPU-only runs are unaffected

The predicted structure is embedded in the HTML report as an interactive
[Mol*](https://molstar.org) 3D viewer with labeled chains (MHC heavy chain,
peptide, TCR α-chain, TCR β-chain).

### Running TCRdock

**Automated (recommended):** Use the cloud GPU script which handles the full
CPU → GPU lifecycle:

```bash
bash scripts/run_cloud_gpu.sh --mode test    # chr22 test run
bash scripts/run_cloud_gpu.sh --mode prod    # full production run
bash scripts/run_cloud_gpu.sh --mode test --detach  # detached (close laptop)
```

See [`docs/google_cloud_guide.md`](docs/google_cloud_guide.md) for details.

**Manual (on a GPU VM):**

```bash
# 1. Set up Docker + TCRdock image
bash scripts/setup_tcrdock_vm.sh config/config.yaml

# 2. Run pipeline with TCRdock overlay
snakemake --cores $(nproc) --use-conda --rerun-triggers mtime \
    --configfile config/config.yaml config/tcrdock_gpu.yaml
```

### Reference

> Bradley P (2023). Structure-based prediction of T cell receptor:peptide-MHC
> interactions. *eLife*, 12, e82813.

---

## Reference Data

Place the following files in the `resources/` directory before running
(see `resources/README.md` for download commands):

| File | Description |
|------|-------------|
| `GRCh38.primary_assembly.genome.fa` | GRCh38 reference genome FASTA + `.fai` index |
| `gencode.v47.annotation.gtf.gz` | GENCODE v47 GTF annotation |

The reference junction BED file (`resources/reference_junctions.bed`) is
generated automatically by the pipeline.

---

## Local Testing (chr22 subset)

Before running a full analysis, you can smoke-test the pipeline end-to-end
using a tiny chr22-only reference and a small FASTQ subset.  The test runs in
~2 minutes on a MacBook Air M1 (8 GB RAM) and exercises every pipeline step.

### Step 1: Download test data (one-time setup, ~15–30 min)

```bash
bash scripts/prepare_test_data.sh
```

This downloads:
- chr22 FASTA from UCSC hg38 (~52 MB)
- GENCODE v47 GTF filtered to chr22 (streamed, no full file stored)
- 500K reads each from a matched gastric cancer pair via ENA HTTPS (no sra-tools needed):
  - **SRR9143066** — Primary Tumor (gastric cancer surgical section)
  - **SRR9143065** — Solid Tissue Normal (adjacent stomach tissue)

### Step 2: Run the test pipeline

```bash
conda activate snakemake
snakemake --cores 4 --use-conda --configfile config/test_config.yaml
```

Expected output: ~234 unannotated junctions (231 tumor_specific, 3 patient_specific) → ~79 contigs → ~11 strong MHC binders.

> **Note**: Only reads mapping to chr22 are used, so junction counts are much
> lower than a full-genome run — this is expected.

---

## Running the Pipeline

### Dry run (check the DAG without executing)

```bash
snakemake --cores 1 --use-conda -n
```

### Full run

```bash
snakemake --cores <N> --use-conda
```

Replace `<N>` with the number of CPU cores to use.

### Run a specific step

```bash
# Only download data for TCGA-BRCA
snakemake --cores 4 --use-conda \
  results/raw_data/TCGA-BRCA/download.done

# Only build the reference junction list
snakemake --cores 2 --use-conda \
  resources/reference_junctions.bed
```

### Cluster execution

Snakemake supports many cluster backends.  Example for SLURM:

```bash
snakemake --cores 100 --use-conda \
  --executor cluster-generic \
  --cluster-generic-submit-cmd "sbatch --mem={resources.mem_mb}M -c {threads}"
```

---

## Configuration

All parameters are in `config/config.yaml`.  Key options:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cancer_types` | BRCA, LUAD, LAML | TCGA project IDs to analyse |
| `reference.genome_fasta` | `resources/…` | Path to GRCh38 FASTA |
| `reference.gencode_gtf` | `resources/…` | Path to GENCODE GTF |
| `mhcflurry.hla_allele` | `HLA-A*02:01` | HLA allele for prediction |
| `mhcflurry.ic50_strong` | `50` | Strong binder threshold (nM) |
| `mhcflurry.ic50_weak` | `500` | Weak binder threshold (nM) |
| `assembly.upstream_nt` | `26` | Nucleotides upstream of junction |
| `assembly.downstream_nt` | `24` | Nucleotides downstream of junction |
| `filtering.min_normal_reads` | `2` | Min reads in normal to trust a junction as patient-specific |
| `tcrdock.enabled` | `false` | Enable TCRdock structural validation (GPU required) |
| `tcrdock.docker_image` | `tcrdock:latest` | Docker image built by `setup_tcrdock_vm.sh` |
| `tcrdock.n_candidates` | `1` | Number of top candidates to model |
| `tcrdock.fallback_hla` | `{A: HLA-A*02:01, B: HLA-B*07:02, C: HLA-C*07:02}` | Fallback HLA alleles as A/B/C mapping (until HLA typing #23) |
| `tcrdock.fallback_tcr` | DMF5 TCR | Fallback TCR sequences (until TRUST4 #24) |

---

## Output Description

All outputs are written to the `results/` directory:

```
results/
├── raw_data/
│   └── {cancer_type}/
│       ├── manifest.tsv          # Sample manifest (file_id → sample_type)
│       └── files/                # Junction quantification TSVs per sample
├── junctions/
│   └── {cancer_type}/
│       └── novel_junctions.tsv   # Classified junctions (junction_origin column:
│                                 #   tumor_specific | patient_specific)
├── contigs/
│   └── {cancer_type}/
│       └── contigs.fa            # 50 nt FASTA contigs (tumor_specific only)
├── peptides/
│   └── {cancer_type}/
│       └── peptides.fa           # 16-mer peptide FASTA
├── predictions/
│   └── {cancer_type}/
│       ├── predictions.tsv       # MHCflurry results
│       └── tcrdock/              # (when TCRdock is enabled)
│           ├── top_candidate.pdb # Predicted TCR-pMHC ternary complex
│           └── docking_scores.tsv # pLDDT/PAE quality metrics
└── reports/
    └── {cancer_type}/
        └── report.html           # Junction origin summary + top binders
                                  # (+ Mol* 3D viewer when TCRdock is enabled)
```

---

## Modernisation Changelog

The following table documents every substantive change from the original 2015
pipeline to this modernised implementation.  See
[`docs/modernization_notes.md`](docs/modernization_notes.md) for full details.

| Component | Original (2015) | Modernised | Reason |
|-----------|----------------|------------|--------|
| Reference genome | hg19 | **GRCh38/hg38** | Current standard assembly |
| Reference annotation | UCSC RefSeq hg19 | **GENCODE v47 GRCh38** | Comprehensive, programmatically reproducible |
| Data source | TCGA HTTP directory (retired) | **GDC Data Portal REST API** | TCGA HTTP retired in 2016 |
| RNA-Seq aligner | TopHat2 | **HISAT2** (default) / **STAR** (optional) | HISAT2 for low-memory runs; STAR planned for production (issue #17) |
| Epitope predictor | NetMHCPan **2.8** | **MHCflurry 2.x** | Open source; no registration; SOTA accuracy |
| Junction-spanning filter | None (all 9-mers included) | **Complete-codon rule** — only 9-mers with ≥1 full codon from each exon retained (issue #18); applied at translation step (issue #20) | Eliminates purely exonic false positives (e.g. `YLADLYHFV` = SH3BP1) |
| Biopython API | `Bio.Alphabet` | **`Bio.Seq` only** | `Bio.Alphabet` removed in Biopython ≥1.78 |
| Workflow management | Manual shell scripts | **Snakemake** | Reproducibility, parallelism, DAG tracking |
| Environment management | None | **Conda** (per rule) | Reproducible software environments |
| Reference junction list | Manual hg19 list | **Derived from GENCODE GTF** | Reproducible, versioned, easily updated |

---

## Project Structure

```
splice-neoepitope-pipeline/
├── README.md
├── LICENSE
├── .gitignore
├── Snakefile                         # Main Snakemake workflow
├── config/
│   ├── config.yaml                   # Production configuration
│   ├── tcrdock_gpu.yaml              # TCRdock GPU overlay (enables tcrdock.enabled)
│   ├── samples.tsv                   # Sample manifest (production)
│   ├── test_config.yaml              # chr22 test configuration (local testing)
│   └── test_samples.tsv              # Sample manifest for test run
├── scripts/
│   ├── prepare_test_data.sh          # One-time setup: download chr22 test data
│   ├── run_cloud_gpu.sh              # Automated CPU→GPU cloud lifecycle
│   ├── setup_cloud.sh                # CPU VM setup (conda, snakemake)
│   └── setup_tcrdock_vm.sh           # GPU VM setup (Docker, TCRdock image)
├── docker/
│   └── Dockerfile.pipeline           # TCRdock Docker image (CUDA 11.8 + AlphaFold)
├── workflow/
│   ├── rules/
│   │   ├── download.smk              # Step 1a: GDC data download
│   │   ├── local_alignment.smk       # Step 1b: STAR local alignment
│   │   ├── hisat2_alignment.smk      # Step 1b: HISAT2 local alignment
│   │   ├── filter.smk                # Step 2: Novel junction filtering
│   │   ├── assemble.smk              # Step 3: Contig assembly
│   │   ├── translate.smk             # Step 4: Peptide translation
│   │   ├── predict.smk               # Step 5: MHCflurry prediction
│   │   ├── tcrdock.smk               # Step 6: TCRdock structural validation (optional, GPU)
│   │   └── analysis.smk              # Step 7: Report generation
│   ├── envs/
│   │   ├── hisat2.yaml               # hisat2, samtools, regtools
│   │   ├── biotools.yaml             # bedtools, biopython, pandas
│   │   ├── star.yaml                 # STAR aligner
│   │   └── python.yaml               # mhcflurry, pandas, scipy, ...
│   └── scripts/
│       ├── download_gdc_data.py      # GDC API download
│       ├── filter_junctions.py       # Novel junction filtering
│       ├── build_reference_junctions.py  # GENCODE reference junction list
│       ├── assemble_contigs.py       # 50 nt contig assembly
│       ├── translate_peptides.py     # In-silico translation
│       ├── run_mhcflurry.py          # MHCflurry 2.x wrapper + parser
│       ├── run_tcrdock.py            # TCRdock Docker wrapper + PDB chain relabelling
│       └── generate_report.py        # HTML report + Mol* 3D viewer
├── docs/
│   ├── INTRODUCTION.md               # Biological background and study design
│   ├── METHODS.md                    # Technical pipeline description
│   ├── DISCUSSIONS.md                # Design tradeoffs and future directions
│   ├── google_cloud_guide.md         # GCP setup and cost guide
│   └── modernization_notes.md        # Detailed change log from 2015 pipeline
└── resources/
    └── README.md                     # Instructions for reference data
```

---

## Citation

If you use this pipeline, please cite the original work:

> Jin-Ho Lee. *"Identification of Cancer-Specific Neoepitopes Arising from
> Alternative Splicing Detected by RNA-Seq."* Seoul National University, 2015.

And the key tools used:

- **Snakemake**: Mölder et al. (2021). Sustainable data analysis with
  Snakemake. *F1000Research*, 10, 33.
- **bedtools**: Quinlan & Hall (2010). BEDTools: a flexible suite of utilities
  for comparing genomic features. *Bioinformatics*, 26(6), 841–842.
- **MHCflurry 2.0**: O'Donnell TJ et al. (2020). MHCflurry 2.0: Improved
  Pan-Allele Prediction of MHC Class I-Presented Peptides by Incorporating
  Antigen Processing. *Cell Systems*, 11(1), 42-48.e7.
- **Biopython**: Cock et al. (2009). Biopython: freely available Python tools
  for computational molecular biology and bioinformatics. *Bioinformatics*,
  25(11), 1422–1423.
- **TCRdock**: Bradley P (2023). Structure-based prediction of T cell
  receptor:peptide-MHC interactions. *eLife*, 12, e82813.
- **Mol\***: Sehnal D et al. (2021). Mol\* Viewer: modern web app for 3D
  visualization and analysis of large biomolecular structures.
  *Nucleic Acids Research*, 49(W1), W431–W437.
- **GENCODE**: Frankish et al. (2023). GENCODE: reference annotation for the
  human and mouse genomes in 2023. *Nucleic Acids Research*, 51(D1),
  D942–D949.

---

## Further Reading

Detailed documentation for the pipeline is available in the `docs/` directory:

| Document | Contents |
|----------|----------|
| [`docs/INTRODUCTION.md`](docs/INTRODUCTION.md) | Biological background, motivation, and study design |
| [`docs/METHODS.md`](docs/METHODS.md) | Technical pipeline description with diagrams and formulas |
| [`docs/DISCUSSIONS.md`](docs/DISCUSSIONS.md) | Design tradeoffs, known limitations, and future directions |
| [`docs/google_cloud_guide.md`](docs/google_cloud_guide.md) | Step-by-step GCP setup and run guide |
