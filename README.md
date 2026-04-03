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
   - [Option B: GDC Download](#option-b-gdc-download-requires-institutional-access)
5. [MHCflurry](#mhcflurry-epitope-predictor)
6. [Reference Data](#reference-data)
7. [Running the Pipeline](#running-the-pipeline)
8. [Configuration](#configuration)
9. [Output Description](#output-description)
10. [Modernisation Changelog](#modernisation-changelog)
11. [Project Structure](#project-structure)
12. [Citation](#citation)

---

## Scientific Background

Cancer cells frequently harbour somatic mutations and aberrant splicing events.
Novel splice junctions — those not present in the reference transcriptome —
can produce peptide sequences that are absent from normal tissues and therefore
recognisable as foreign by the immune system.  These **neoepitopes** are
candidate targets for cancer immunotherapy.

This pipeline:
1. Generates splice junction quantification from RNA-Seq data (local alignment or GDC download).
2. Identifies novel (non-reference) splice junctions enriched in tumour samples.
3. Constructs short nucleotide contigs spanning each junction.
4. Translates the contigs into peptides in all three reading frames.
5. Predicts MHC-I binding affinity using MHCflurry 2.x.
6. Statistically evaluates enrichment of predicted epitopes in tumour vs. normal.

Cancer types analysed: **BRCA** (breast adenocarcinoma), **LUAD** (lung
adenocarcinoma), **LAML** (acute myeloid leukemia).

---

## Pipeline Overview

```
RNA-Seq data (local FASTQ files or GDC API)
        │
        ▼ Step 1: Align/Download
  Splice junction quantification files (.tsv)
        │
        ▼ Step 2: Filter
  Novel junctions (read count > mean; not in GENCODE reference)
        │
        ▼ Step 3: Assemble
  50 nt contigs (26 nt upstream + 24 nt downstream, bedtools getfasta)
        │
        ▼ Step 4: Translate
  16-mer peptides (3 reading frames, truncated at stop codons)
        │
        ▼ Step 5: Predict
  MHCflurry 2.x → IC50 binding affinities (HLA-A*02:01, 9-mers)
        │
        ▼ Step 6: Analyse
  Fisher's exact test, summary statistics, HTML report
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
| **Local Alignment** | ❌ Not required | Align your own FASTQ files using STAR |
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
data_source: "local"   # Use local FASTQ alignment

local_samples:
  samples_tsv: "config/samples.tsv"
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

**Example: Download SRA data with `fasterq-dump`:**

```bash
# Install SRA toolkit
conda install -c conda-forge -c bioconda sra-tools

# Download example breast cancer RNA-Seq (SRR12345678 is a placeholder)
fasterq-dump SRR12345678 --split-files --outdir data/

# This creates:
#   data/SRR12345678_1.fastq
#   data/SRR12345678_2.fastq
```

**Example datasets for testing:**

| Accession | Description | Sample Type |
|-----------|-------------|-------------|
| SRP064305 | TCGA-BRCA-adjacent (public) | Breast cancer |
| SRP066790 | Lung adenocarcinoma cell lines | Lung cancer |
| SRP055401 | AML cell lines | Leukemia |

> **Tip**: Search GEO for `"RNA-Seq" AND "cancer" AND "Homo sapiens"` to find
> publicly available cancer RNA-Seq datasets.

#### Step 3: Create Sample Manifest

Create a TSV file listing your samples at `config/samples.tsv`:

```tsv
sample_id	sample_type	fastq1	fastq2
tumor_01	Primary Tumor	data/tumor_01_R1.fastq.gz	data/tumor_01_R2.fastq.gz
tumor_02	Primary Tumor	data/tumor_02_R1.fastq.gz	data/tumor_02_R2.fastq.gz
normal_01	Solid Tissue Normal	data/normal_01_R1.fastq.gz	data/normal_01_R2.fastq.gz
```

**Column descriptions:**

| Column | Required | Description |
|--------|----------|-------------|
| `sample_id` | Yes | Unique identifier for the sample |
| `sample_type` | Yes | `"Primary Tumor"` or `"Solid Tissue Normal"` for Fisher's test |
| `fastq1` | Yes | Path to read 1 FASTQ (can be gzipped) |
| `fastq2` | No | Path to read 2 FASTQ for paired-end data (leave empty for single-end) |

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

If you don't have a machine with sufficient RAM, you can run the alignment
step on cloud platforms:

| Platform | Cost | Description |
|----------|------|-------------|
| **Galaxy** | Free | [usegalaxy.org](https://usegalaxy.org) — Run HISAT2/STAR in browser |
| **Google Colab** | Free tier | Use free GPU instances with ~12 GB RAM |
| **AWS** | Pay-per-use | r5.xlarge instance (~$0.25/hr, 32 GB RAM) |
| **Google Cloud** | Pay-per-use | n2-highmem-4 (~$0.20/hr, 32 GB RAM) |

**Galaxy Workflow (Free):**
1. Upload your FASTQ files to [usegalaxy.org](https://usegalaxy.org)
2. Search for "HISAT2" or "STAR" in the tool panel
3. Run alignment with GRCh38/hg38 genome
4. Download the junction output file (`.tsv` or `.bed`)
5. Place in `results/raw_data/local/files/` and continue with the pipeline

---

### Option B: GDC Download (Requires Institutional Access)

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

### First-time setup: Download models

After the first Snakemake run creates the conda environment, download the
MHCflurry trained models:

```bash
# Activate the pipeline's Python environment
conda activate .snakemake/conda/<hash>   # or run within a Snakemake job
mhcflurry-downloads fetch
```

Alternatively, you can run this once before the pipeline:

```bash
pip install mhcflurry
mhcflurry-downloads fetch
```

The models (~1 GB) are cached in `~/.local/share/mhcflurry/` and reused.

### Reference

> O'Donnell TJ et al. (2020). MHCflurry 2.0: Improved Pan-Allele Prediction
> of MHC Class I-Presented Peptides by Incorporating Antigen Processing.
> *Cell Systems*, 11(1), 42-48.e7.

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
| `filtering.strategy` | `mean` | Read-count filter strategy |

---

## Output Description

All outputs are written to the `results/` directory:

```
results/
├── raw_data/
│   └── {cancer_type}/
│       ├── manifest.tsv          # GDC file manifest
│       └── files/                # Downloaded junction quantification TSVs
├── junctions/
│   └── {cancer_type}/
│       └── novel_junctions.tsv   # Filtered novel junctions
├── contigs/
│   └── {cancer_type}/
│       └── contigs.fa            # 50 nt FASTA contigs
├── peptides/
│   └── {cancer_type}/
│       └── peptides.fa           # 16-mer peptide FASTA
├── predictions/
│   └── {cancer_type}/
│       └── predictions.tsv       # MHCflurry results
├── analysis/
│   └── {cancer_type}/
│       └── statistics.tsv        # Fisher's test + epitope counts
└── reports/
    ├── summary_table.tsv         # Cross-cancer summary
    └── {cancer_type}/
        └── report.html           # Per-cancer HTML report
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
| RNA-Seq aligner | TopHat2 | **STAR** (GDC harmonised) | GDC re-aligned all TCGA data with STAR |
| Epitope predictor | NetMHCPan **2.8** | **MHCflurry 2.x** | Open source; no registration; SOTA accuracy |
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
│   └── config.yaml                   # All configurable parameters
├── workflow/
│   ├── rules/
│   │   ├── download.smk              # Step 1: GDC data download
│   │   ├── filter.smk                # Step 2: Novel junction filtering
│   │   ├── assemble.smk              # Step 3: Contig assembly
│   │   ├── translate.smk             # Step 4: Peptide translation
│   │   ├── predict.smk               # Step 5: MHCflurry prediction
│   │   └── analysis.smk              # Step 6: Statistics + reports
│   ├── envs/
│   │   ├── biotools.yaml             # bedtools, samtools, pysam
│   │   └── python.yaml               # Python, biopython, mhcflurry, pandas, scipy, ...
│   └── scripts/
│       ├── download_gdc_data.py      # GDC API download
│       ├── filter_junctions.py       # Novel junction filtering
│       ├── build_reference_junctions.py  # GENCODE reference junction list
│       ├── assemble_contigs.py       # 50 nt contig assembly
│       ├── translate_peptides.py     # In-silico translation (modern Biopython)
│       ├── run_mhcflurry.py          # MHCflurry 2.x wrapper + parser
│       ├── statistical_analysis.py   # Fisher's exact test + summary stats
│       └── generate_report.py        # HTML report generation
├── docs/
│   └── modernization_notes.md        # Detailed change log
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
- **GENCODE**: Frankish et al. (2023). GENCODE: reference annotation for the
  human and mouse genomes in 2023. *Nucleic Acids Research*, 51(D1),
  D942–D949.
