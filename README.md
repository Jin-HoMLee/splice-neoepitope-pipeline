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
4. [GDC Authentication](#gdc-authentication-required-for-tcga-data)
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
1. Downloads TCGA RNA-Seq splice-junction quantification data.
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
TCGA RNA-Seq data (GDC API)
        │
        ▼ Step 1: Download
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
| RAM | 16 GB | 32 GB recommended for LUAD/BRCA |
| Disk | 50 GB free | Reference genome + TCGA downloads |
| Python | 3.11+ | Managed automatically via conda |
| Git | any | For cloning this repository |

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

## GDC Authentication (Required for TCGA Data)

**TCGA data on the GDC Data Portal is controlled-access** and requires
authentication. You will need an eRA Commons account linked to dbGaP access
for TCGA. Without a valid token, downloads will fail with `403 Forbidden`.

### Step 1: Obtain dbGaP Access

1. Apply for controlled-access TCGA data via [dbGaP](https://dbgap.ncbi.nlm.nih.gov/)
2. Your institution's signing official must approve the Data Access Request (DAR)
3. This process typically takes 1–4 weeks

> **Note**: If you already have TCGA access through another project (e.g.,
> via your institution's blanket approval), you can skip to Step 2.

### Step 2: Download Your GDC Token

1. Go to **https://portal.gdc.cancer.gov/**
2. Click **Login** (top right) and authenticate via eRA Commons / NIH
3. After login, click your username (top right) → **Download Token**
4. Save the downloaded file (e.g., `gdc-user-token.txt`)

The token is valid for **30 days**; you'll need to download a new one after
it expires.

### Step 3: Configure the Pipeline

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

### Troubleshooting

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
