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
3. [Quick Start](#quick-start)
4. [Installation](#installation)
5. [Data Preparation](#data-preparation)
6. [TCRdock Structural Validation (Optional)](#tcrdock-structural-validation-optional)
7. [Running the Pipeline](#running-the-pipeline)
8. [Configuration](#configuration)
9. [Output](#output)
10. [Project Structure](#project-structure)
11. [Citation](#citation)
12. [Further Reading](#further-reading)

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
RNA-Seq FASTQ files
        │
        ▼ Step 1: Align (HISAT2 or STAR)
  Splice junction quantification files (.tsv)
        │
        ├─────────────────────────────────────────────────────┐
        ▼                                                     ▼
  Step 2: Classify by origin                        Step 2b: HLA typing (optional)
  - Remove annotated junctions (GENCODE)            OptiType on all sample FASTQs
  - Compare tumor vs. matched normal:               → patient-specific A/B/C alleles
      normal_shared (in normal) → excluded
      tumor_exclusive (not in normal) → keep
        │                                                     │
        ▼─────────────────────────────────────────────────────┘
        │
        ▼ Step 3: Assemble
  50 nt contigs (26 nt upstream + 24 nt downstream, bedtools getfasta)
        │
        ▼ Step 4: Translate
  Junction-spanning 9-mers (3 reading frames, complete-codon filter)
        │
        ▼ Step 5: Predict
  MHCflurry 2.x Class1PresentationPredictor → presentation_score + IC50 per peptide (best allele)
        │
        ▼ Step 6: TCRdock structural validation (optional, GPU)
  TCR-peptide-MHC ternary complex 3D structure (AlphaFold v2 backend)
        │
        ▼ Step 7: Report
  Junction origin summary + HLA typing QC + top binders + Mol* 3D viewer HTML report
```

---

## Quick Start

### Local test run (chr22, ~2 min on MacBook Air M1)

```bash
bash scripts/prepare_test_data.sh   # one-time: downloads chr22 reference + FASTQs
conda activate snakemake
snakemake --cores 4 --use-conda --configfile config/test_config.yaml
```

Expected output: ~372 unannotated junctions → ~75 contigs → ~147 strong presenters
(best allele across 6 patient-specific HLA-A/B/C alleles).

### Cloud run (full genome + TCRdock, ~4–6 hours)

```bash
bash scripts/run_cloud_gpu.sh \
    --samples config/samples/patient_001.tsv \
    --mode prod \
    --detach
```

> See [`docs/google_cloud_guide.md`](docs/google_cloud_guide.md) for prerequisites and cost estimates.

---

## Installation

Requires Conda (Miniforge recommended) and Snakemake 8.x. Rule-specific
environments are created automatically on first use — no manual dependency
installation needed.

```bash
git clone https://github.com/Jin-HoMLee/splice-neoepitope-pipeline.git
cd splice-neoepitope-pipeline
conda create -n snakemake -c conda-forge -c bioconda "snakemake>=8.0,<9" python=3.11 -y
conda activate snakemake
```

See [`docs/installation.md`](docs/installation.md) for the full setup guide
(multi-platform conda install, reference data download, sanity check).

---

## Data Preparation

Per-patient sample manifests live in `config/samples/`. FASTQ paths can be
local files, `gs://` URIs, or `https://` URLs — all downloaded automatically.

```tsv
patient_id  sample_id   sample_type          fastq1                    fastq2
patient_001 SRR9143066  Primary Tumor        https://ftp.sra.ebi.ac.uk/...
patient_001 SRR9143065  Solid Tissue Normal  https://ftp.sra.ebi.ac.uk/...
```

See [`docs/data_preparation.md`](docs/data_preparation.md) for aligner selection
(HISAT2 vs STAR), FASTQ sources, and full manifest format.

---

## TCRdock Structural Validation (Optional)

When enabled, predicts the 3D TCR-pMHC ternary complex for the top neoepitope
candidate using TCRdock (AlphaFold v2 backend). Requires a Linux x86-64 NVIDIA
GPU. Disabled by default — local / CPU-only runs are unaffected.

The predicted structure is embedded in the HTML report as an interactive
[Mol*](https://molstar.org) 3D viewer.

**Automated cloud run (recommended):**

```bash
bash scripts/run_cloud_gpu.sh --samples config/samples/patient_001.tsv --mode prod --detach
```

See [`docs/google_cloud_guide.md`](docs/google_cloud_guide.md) for prerequisites,
cost estimates, and manual GPU VM setup.

---

## Running the Pipeline

### Dry run

```bash
conda activate snakemake
snakemake --cores 1 --use-conda -n \
    --configfile config/config.yaml \
    --config samples_tsv=config/samples/patient_001.tsv
```

### Full run

```bash
conda activate snakemake
snakemake --cores <N> --use-conda \
    --rerun-triggers mtime \
    --rerun-incomplete \
    --configfile config/config.yaml \
    --config samples_tsv=config/samples/patient_001.tsv
```

### Run a specific target

```bash
snakemake --cores 8 --use-conda \
    --configfile config/config.yaml \
    --config samples_tsv=config/samples/patient_001.tsv \
    results/patient_001/alignment/download.done
```

### Cluster execution (SLURM example)

```bash
snakemake --cores 100 --use-conda \
    --rerun-triggers mtime --rerun-incomplete \
    --configfile config/config.yaml \
    --config samples_tsv=config/samples/patient_001.tsv \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "sbatch --mem={resources.mem_mb}M -c {threads}"
```

---

## Configuration

All parameters live in `config/config.yaml`. For TCRdock, merge
`config/gpu.yaml` via a single `--configfile` invocation.

See [`docs/configuration.md`](docs/configuration.md) for the full parameter reference.

---

## Output

All outputs are written to `results/{patient_id}/`:

```
results/
└── {patient_id}/
    ├── alignment/
    │   ├── manifest.tsv              # Sample manifest (file_id → sample_type)
    │   └── {sample}/
    │       └── junctions.tsv         # Junction read counts from aligner
    ├── hla_typing/                   # (when hla.enabled: true)
    │   ├── {sample}/
    │   │   └── {sample}_result.tsv   # Per-sample OptiType output
    │   ├── alleles.tsv               # Aggregated patient HLA-A/B/C alleles
    │   └── hla_qc.tsv                # Per-locus source, read counts, discrepancies
    ├── junctions/
    │   └── novel_junctions.tsv       # Classified junctions (tumor_exclusive | normal_shared)
    ├── contigs/
    │   └── contigs.fa                # 50 nt FASTA contigs (tumor_exclusive only)
    ├── peptides/
    │   └── peptides.tsv              # Junction-spanning 9-mers
    ├── predictions/
    │   ├── mhc_presentation.tsv      # MHCflurry results (one row per peptide, best allele)
    │   └── tcrdock/                  # (when TCRdock is enabled)
    │       ├── top_candidate.pdb     # Predicted TCR-pMHC ternary complex
    │       └── docking_scores.tsv    # pLDDT/PAE quality metrics
    └── reports/
        └── report.html               # Summary report (+ Mol* 3D viewer when TCRdock enabled)
```

---

## Project Structure

```
splice-neoepitope-pipeline/
├── README.md
├── LICENSE
├── Snakefile                         # Main Snakemake workflow
├── config/
│   ├── config.yaml                   # Production configuration
│   ├── gpu.yaml                      # GPU overlay (MHCflurry + TCRdock)
│   ├── test_config.yaml              # chr22 test configuration
│   └── samples/
│       ├── patient_001.tsv           # Gastric cancer (SRR9143066/SRR9143065)
│       ├── patient_001_test.tsv      # chr22 subset for local testing
│       └── patient_002.tsv           # Osteosarcoma (BostonGene BG003082)
├── scripts/
│   ├── prepare_test_data.sh          # One-time: download chr22 test data
│   ├── run_cloud_gpu.sh              # Automated CPU→GPU cloud lifecycle
│   ├── setup_cloud.sh                # CPU VM setup (conda, snakemake)
│   └── setup_tcrdock_vm.sh           # GPU VM setup (Docker, TCRdock image)
├── docker/
│   └── Dockerfile.pipeline           # TCRdock image (CUDA 11.8 + AlphaFold)
├── workflow/
│   ├── rules/
│   │   ├── alignment.smk             # Step 1: HISAT2 or STAR alignment
│   │   ├── filter_junctions.smk      # Step 2: Novel junction filtering
│   │   ├── hla_typing.smk            # Step 2b: OptiType HLA typing (optional)
│   │   ├── assemble_contigs.smk      # Step 3: Contig assembly
│   │   ├── translate_peptides.smk    # Step 4: Peptide translation
│   │   ├── mhc_affinity.smk          # Step 5: MHCflurry prediction
│   │   ├── structure.smk             # Step 6: TCRdock (optional, GPU)
│   │   └── analysis.smk              # Step 7: Report generation
│   ├── envs/
│   │   ├── hisat2.yaml               # hisat2, samtools, regtools
│   │   ├── biotools.yaml             # bedtools, biopython, pandas
│   │   ├── star.yaml                 # STAR aligner
│   │   ├── optitype.yaml             # OptiType, razers3, glpk
│   │   └── python.yaml               # mhcflurry, pandas, scipy, ...
│   └── scripts/
│       ├── filter_junctions.py
│       ├── build_reference_junctions.py
│       ├── assemble_contigs.py
│       ├── translate_peptides.py
│       ├── aggregate_hla_alleles.py
│       ├── run_mhcflurry.py
│       ├── run_tcrdock.py
│       └── generate_report.py
├── docs/
│   ├── INTRODUCTION.md               # Biological background and study design
│   ├── DISCUSSIONS.md                # Design tradeoffs and future directions
│   ├── google_cloud_guide.md         # GCP setup, cost guide, GPU run instructions
│   ├── installation.md               # Full installation guide
│   ├── data_preparation.md           # Aligner selection, FASTQ sources, manifest format
│   ├── configuration.md              # Full config.yaml parameter reference
│   └── modernization_notes.md        # Changes from the 2015 pipeline
└── resources/
    └── README.md                     # Reference data download instructions
```

---

## Citation

If you use this pipeline, please cite the original work:

> Jin-Ho Lee. *"Identification of Cancer-Specific Neoepitopes Arising from
> Alternative Splicing Detected by RNA-Seq."* Seoul National University, 2015.

And the key tools used:

- **Snakemake**: Mölder et al. (2021). Sustainable data analysis with Snakemake. *F1000Research*, 10, 33.
- **HISAT2**: Kim et al. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nature Biotechnology*, 37, 907–915.
- **bedtools**: Quinlan & Hall (2010). BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*, 26(6), 841–842.
- **OptiType**: Szolek A et al. (2014). OptiType: precision HLA typing from next-generation sequencing data. *Bioinformatics*, 30(23), 3310–3316.
- **MHCflurry 2.0**: O'Donnell TJ et al. (2020). MHCflurry 2.0: Improved Pan-Allele Prediction of MHC Class I-Presented Peptides by Incorporating Antigen Processing. *Cell Systems*, 11(1), 42-48.e7.
- **TCRdock**: Bradley P (2023). Structure-based prediction of T cell receptor:peptide-MHC interactions. *eLife*, 12, e82813.
- **Mol\***: Sehnal D et al. (2021). Mol\* Viewer: modern web app for 3D visualization and analysis of large biomolecular structures. *Nucleic Acids Research*, 49(W1), W431–W437.
- **GENCODE**: Frankish et al. (2023). GENCODE: reference annotation for the human and mouse genomes in 2023. *Nucleic Acids Research*, 51(D1), D942–D949.

---

## Further Reading

| Document | Contents |
|----------|----------|
| [`docs/installation.md`](docs/installation.md) | Full installation guide (conda, Snakemake, reference data) |
| [`docs/data_preparation.md`](docs/data_preparation.md) | Aligner selection, FASTQ sources, sample manifest format |
| [`docs/configuration.md`](docs/configuration.md) | Full `config/config.yaml` parameter reference |
| [`docs/google_cloud_guide.md`](docs/google_cloud_guide.md) | GCP setup, cost estimates, GPU VM lifecycle |
| [`docs/INTRODUCTION.md`](docs/INTRODUCTION.md) | Biological background, motivation, and study design |
| [`docs/DISCUSSIONS.md`](docs/DISCUSSIONS.md) | Design tradeoffs, known limitations, and future directions |
| [`docs/modernization_notes.md`](docs/modernization_notes.md) | Detailed changes from the 2015 pipeline |
