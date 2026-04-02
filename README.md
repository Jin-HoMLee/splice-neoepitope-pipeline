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
4. [NetMHCPan 4.1 Installation](#netmhcpan-41-installation)
5. [Reference Data](#reference-data)
6. [Running the Pipeline](#running-the-pipeline)
7. [Configuration](#configuration)
8. [Output Description](#output-description)
9. [Modernisation Changelog](#modernisation-changelog)
10. [Project Structure](#project-structure)
11. [Citation](#citation)

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
5. Predicts MHC-I binding affinity using NetMHCPan 4.1.
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
  NetMHCPan 4.1 → IC50 binding affinities (HLA-A*02:01, 9-mers)
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

## NetMHCPan 4.1 Installation

NetMHCPan 4.1 is **free for academic use** but requires registration with DTU
Bioinformatics.

1. **Register and download** — fill in the form at
   **<https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/>**.
   You will receive a download link by e-mail within minutes.

2. **Unpack the archive**:
   ```bash
   tar -xzf netMHCpan-4.1b.Linux.tar.gz   # filename may differ slightly
   cd netMHCpan-4.1
   ```

3. **Set the install directory** — edit the first line of the `netMHCpan`
   script to point to its location:
   ```bash
   # Open the script in a text editor and set:
   #   setenv  NMHCPAN  /full/path/to/netMHCpan-4.1
   # e.g.
   sed -i "s|/usr/cbs/packages/netMHCpan/4.1|$(pwd)|" netMHCpan
   ```

4. **Add to PATH** (or set `netmhcpan.executable` in `config/config.yaml`):
   ```bash
   export PATH="$(pwd):$PATH"
   # Add this line to ~/.bashrc or ~/.bash_profile to make it permanent
   ```

5. **Verify**:
   ```bash
   netMHCpan -h 2>&1 | head -5
   # Should print: # NetMHCpan version 4.1b ...
   ```



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
| `netmhcpan.executable` | `netMHCpan` | Path to NetMHCPan binary |
| `netmhcpan.hla_allele` | `HLA-A02:01` | HLA allele for prediction |
| `netmhcpan.ic50_strong` | `50` | Strong binder threshold (nM) |
| `netmhcpan.ic50_weak` | `500` | Weak binder threshold (nM) |
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
│       └── predictions.tsv       # NetMHCPan results
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
| Epitope predictor | NetMHCPan **2.8** | **NetMHCPan 4.1** | Improved accuracy; eluted-ligand training data |
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
│   │   ├── predict.smk               # Step 5: NetMHCPan prediction
│   │   └── analysis.smk              # Step 6: Statistics + reports
│   ├── envs/
│   │   ├── biotools.yaml             # bedtools, samtools, pysam
│   │   └── python.yaml               # Python, biopython, pandas, scipy, ...
│   └── scripts/
│       ├── download_gdc_data.py      # GDC API download
│       ├── filter_junctions.py       # Novel junction filtering
│       ├── build_reference_junctions.py  # GENCODE reference junction list
│       ├── assemble_contigs.py       # 50 nt contig assembly
│       ├── translate_peptides.py     # In-silico translation (modern Biopython)
│       ├── run_netmhcpan.py          # NetMHCPan 4.1 wrapper + parser
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
- **NetMHCPan 4.1**: Reynisson et al. (2020). NetMHCpan-4.1 and
  NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by
  concurrent motif deconvolution and integration of MS MHC eluted ligand data.
  *Nucleic Acids Research*, 48(W1), W449–W454.
- **Biopython**: Cock et al. (2009). Biopython: freely available Python tools
  for computational molecular biology and bioinformatics. *Bioinformatics*,
  25(11), 1422–1423.
- **GENCODE**: Frankish et al. (2023). GENCODE: reference annotation for the
  human and mouse genomes in 2023. *Nucleic Acids Research*, 51(D1),
  D942–D949.
