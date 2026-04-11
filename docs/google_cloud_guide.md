# Running the Pipeline on Google Cloud Platform

This guide provides step-by-step instructions for running the splice neoepitope
pipeline on Google Cloud Platform (GCP), ideal for users who:
- Don't have sufficient local resources (32+ GB RAM for STAR, or 8 GB for HISAT2)
- Need to process large-scale datasets
- Want reproducible cloud-based analyses

---

## Table of Contents

1. [Overview](#overview)
2. [Cost Estimation](#cost-estimation)
3. [Prerequisites](#prerequisites)
4. [Option A: Compute Engine VM (Recommended)](#option-a-compute-engine-vm-recommended)
5. [Option B: Google Cloud Life Sciences API](#option-b-google-cloud-life-sciences-api)
6. [Option C: Google Batch](#option-c-google-batch)
7. [Data Transfer](#data-transfer)
8. [Automated GPU Pipeline (TCRdock)](#automated-gpu-pipeline-tcrdock)
9. [Troubleshooting](#troubleshooting)

---

## Overview

| Approach | Best For | Complexity | Cost |
|----------|----------|------------|------|
| **Compute Engine VM** | Interactive development, small datasets | Low | ~$0.20–$1.00/hr |
| **Life Sciences API** | Automated batch workflows | Medium | ~$0.10–$0.50/hr |
| **Google Batch** | Large-scale parallelisation | Medium-High | Variable |

For most users, we recommend **Option A (Compute Engine VM)** as it's the
simplest to set up and provides a familiar Linux environment.

---

## Cost Estimation

| Machine Type | vCPUs | RAM | Hourly Cost* | Best For |
|--------------|-------|-----|--------------|----------|
| `e2-standard-2` | 2 | 8 GB | ~$0.07 | HISAT2 alignment (small samples) |
| `n2-standard-4` | 4 | 16 GB | ~$0.19 | HISAT2 alignment (typical) |
| `n2-highmem-4` | 4 | 32 GB | ~$0.26 | STAR alignment |
| `n2-highmem-8` | 8 | 64 GB | ~$0.52 | STAR alignment (large samples) |

*Prices as of 2024; actual costs vary by region. Spot/preemptible VMs are ~60–80% cheaper.

**Estimated total costs per analysis:**
- Small test run (2 samples, HISAT2): ~$1–5
- Full analysis (10 samples, HISAT2): ~$10–30
- Full analysis (10 samples, STAR): ~$30–100

---

## Prerequisites

### 1. Google Cloud Account

1. Go to [console.cloud.google.com](https://console.cloud.google.com/)
2. Create a Google Cloud account (if you don't have one)
3. New users get **$300 free credits** for 90 days

### 2. Create a Project

```bash
# Install Google Cloud CLI (if not already installed)
# https://cloud.google.com/sdk/docs/install

# Authenticate
gcloud auth login

# Create a new project (or use existing)
gcloud projects create splice-neoepitope-project --name="Splice Neoepitope"

# Set as active project
gcloud config set project splice-neoepitope-project

# Enable billing (required for Compute Engine)
# Do this in the Cloud Console: https://console.cloud.google.com/billing
```

### 3. Enable Required APIs

```bash
gcloud services enable compute.googleapis.com
gcloud services enable storage.googleapis.com
```

---

## Option A: Compute Engine VM (Recommended)

This is the simplest approach — create a Linux VM and run the pipeline as you
would on any local machine.

### Step 1: Create a VM Instance

**Using Cloud Console (GUI):**

1. Go to **Compute Engine** → **VM instances** → **Create Instance**
2. Configure:
   - **Name**: `splice-pipeline`
   - **Region**: Choose closest to you (or `us-central1` for lowest cost)
   - **Machine type**: 
     - For HISAT2: `n2-standard-4` (4 vCPUs, 16 GB RAM) — ~$0.19/hr
     - For STAR: `n2-highmem-4` (4 vCPUs, 32 GB RAM) — ~$0.26/hr
   - **Boot disk**: 
     - OS: Ubuntu 22.04 LTS
     - Size: **100 GB** (or more for large datasets)
     - Type: SSD persistent disk (for faster I/O)
   - **Firewall**: Allow HTTP/HTTPS (optional, for downloading data)
3. Click **Create**

**Using gcloud CLI:**

> **Zone availability:** `us-central1-a` is used as the example zone below, but
> capacity for specific machine types varies. If you get a quota or availability
> error, try another zone (e.g. `us-central1-c`, `europe-west1-b`, `us-east1-b`).
> Run `gcloud compute zones list` to see all available zones.

```bash
# For HISAT2 (8 GB RAM sufficient, using 16 GB for comfort)
gcloud compute instances create splice-pipeline \
    --zone=us-central1-a \
    --machine-type=n2-standard-4 \
    --boot-disk-size=100GB \
    --boot-disk-type=pd-ssd \
    --image-family=ubuntu-2204-lts \
    --image-project=ubuntu-os-cloud

# For STAR (32 GB RAM required)
gcloud compute instances create splice-pipeline \
    --zone=us-central1-a \
    --machine-type=n2-highmem-4 \
    --boot-disk-size=150GB \
    --boot-disk-type=pd-ssd \
    --image-family=ubuntu-2204-lts \
    --image-project=ubuntu-os-cloud
```

**💡 Cost-Saving Tip: Use Spot VMs**

Spot (preemptible) VMs are ~60–80% cheaper but can be terminated with 30s notice:

```bash
gcloud compute instances create splice-pipeline \
    --zone=us-central1-a \
    --machine-type=n2-highmem-4 \
    --boot-disk-size=150GB \
    --boot-disk-type=pd-ssd \
    --image-family=ubuntu-2204-lts \
    --image-project=ubuntu-os-cloud \
    --provisioning-model=SPOT \
    --instance-termination-action=STOP
```

### Step 2: Connect to the VM

```bash
gcloud compute ssh splice-pipeline --zone=us-central1-a
```

Or use the **SSH** button in the Cloud Console.

### Step 3: Install Dependencies

Once connected to the VM:

```bash
# Update system packages
sudo apt update && sudo apt upgrade -y

# Install required system packages
sudo apt install -y git curl wget

# Install Miniforge (conda)
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b -p $HOME/miniforge3
rm Miniforge3-*.sh

# Initialize conda
~/miniforge3/bin/conda init bash
source ~/.bashrc

# Install Snakemake and samtools
conda install -n base -c conda-forge -c bioconda "snakemake>=7.0,<9" python=3.11 samtools -y
```

### Step 4: Clone and Configure the Pipeline

```bash
# Clone the repository
git clone https://github.com/Jin-HoMLee/splice-neoepitope-pipeline.git
cd splice-neoepitope-pipeline

# Verify structure
ls -la

# Verify Snakemake
snakemake --version
```

### Step 5: Download Reference Files

Download the GRCh38 genome and GENCODE annotation into the `resources/` directory.
These files are required before running the pipeline (~3.1 GB uncompressed + 57 MB).

```bash
cd ~/splice-neoepitope-pipeline/resources

# Download GENCODE v47 GTF annotation (~57 MB)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz

# Download GRCh38 primary assembly genome (~840 MB compressed, ~3.1 GB uncompressed)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

# Index the genome (required by bedtools getfasta in the assembly step)
samtools faidx GRCh38.primary_assembly.genome.fa

cd ..
```

### Step 6: Upload Your Data

**Option A: Upload from local machine using gcloud**

```bash
# From your LOCAL machine (not the VM):
gcloud compute scp /path/to/your/sample_R1.fq.gz splice-pipeline:~/splice-neoepitope-pipeline/data/ --zone=us-central1-a
gcloud compute scp /path/to/your/sample_R2.fq.gz splice-pipeline:~/splice-neoepitope-pipeline/data/ --zone=us-central1-a
```

**Option B: Download directly from SRA (recommended)**

```bash
# On the VM — install SRA Toolkit
# Note: use version 3.1.1; newer versions (3.4.x) have a segfault bug on GCP VMs
conda install -c bioconda sra-tools=3.1.1 -y

# Download example breast cancer RNA-Seq data
mkdir -p data

# Download a sample (replace SRR10971381 with your accession)
# Example: SRR37781424 (Luminal A breast cancer tumor, ~3.3 GB SRA)
fasterq-dump SRR10971381 --split-files --outdir data/

# Rename to match expected format
mv data/SRR10971381_1.fastq data/tumor_01_R1.fastq
mv data/SRR10971381_2.fastq data/tumor_01_R2.fastq
```

**Option C: Use Google Cloud Storage**

```bash
# Create a bucket (from your local machine or Cloud Shell)
gsutil mb gs://splice-pipeline-data

# Upload data to bucket
gsutil cp /path/to/your/*.fq.gz gs://splice-pipeline-data/

# On the VM, download from bucket
gsutil cp gs://splice-pipeline-data/*.fq.gz ~/splice-neoepitope-pipeline/data/
```

### Step 7: Configure samples.tsv

```bash
cd ~/splice-neoepitope-pipeline

# Edit samples.tsv
nano config/samples.tsv
```

Add your samples:
```tsv
sample_id	sample_type	fastq1	fastq2
tumor_01	Primary Tumor	data/tumor_01_R1.fastq	data/tumor_01_R2.fastq
normal_01	Solid Tissue Normal	data/normal_01_R1.fastq	data/normal_01_R2.fastq
```

### Step 8: Run the Pipeline

The pipeline takes several hours. Use `tmux` so it keeps running if your SSH connection drops:

```bash
# Start a tmux session
tmux new -s pipeline

# Inside tmux: dry run first
snakemake --cores 4 --use-conda -n

# Inside tmux: full run — pipeline logs to pipeline.log
snakemake --cores $(nproc) --use-conda --rerun-triggers mtime 2>&1 | tee pipeline.log
```

To detach from tmux without stopping the pipeline: press `Ctrl+B`, then `D`.
To reattach after reconnecting via SSH: `tmux attach -t pipeline`

> **Important:** The VM will NOT shut down automatically — remember to stop it
> when the run finishes (see [Step 10](#step-10-stopdelete-the-vm)).

> **Why tmux instead of nohup?** With `nohup`, you lose visibility into the running pipeline after disconnecting. With `tmux`, you can reattach at any time and see live output.

> **Prefer a fully automated workflow?** `run_cloud_gpu.sh` handles the full
> lifecycle — CPU steps, GCS handoff, GPU TCRdock, and automatic VM shutdown:
> ```bash
> # From your local machine (VMs are managed automatically):
> bash scripts/run_cloud_gpu.sh --mode prod
>
> # Detached mode — runs on an orchestrator VM so you can close your laptop:
> bash scripts/run_cloud_gpu.sh --mode prod --detach
> ```
> See [Automated GPU Pipeline](#automated-gpu-pipeline-tcrdock) for details.

### Step 9: Download Results

```bash
# From your LOCAL machine:
gcloud compute scp --recurse splice-pipeline:~/splice-neoepitope-pipeline/results/ ./results/ --zone=us-central1-a
```

### Step 10: Stop/Delete the VM

**Important: Stop your VM when not in use to avoid charges!**

```bash
# Stop (preserves disk, ~$5/month for 100 GB SSD)
gcloud compute instances stop splice-pipeline --zone=us-central1-a

# Start again later
gcloud compute instances start splice-pipeline --zone=us-central1-a

# Delete completely (when finished)
gcloud compute instances delete splice-pipeline --zone=us-central1-a
```

---

## Option B: Google Cloud Life Sciences API

For automated batch processing without maintaining a VM. This approach submits
the pipeline as a job that runs to completion.

### Step 1: Enable Life Sciences API

```bash
gcloud services enable lifesciences.googleapis.com
```

### Step 2: Create a Snakemake Profile for GCP

Create `~/.config/snakemake/gcp/config.yaml`:

```yaml
# Google Cloud Life Sciences profile
executor: google-lifesciences
google-lifesciences-region: us-central1
default-resources:
  - runtime=120  # minutes
  - mem_mb=16000
jobs: 50
use-conda: true
```

### Step 3: Run with Life Sciences

```bash
# Set up Google credentials
gcloud auth application-default login

# Run pipeline
snakemake --profile gcp --google-lifesciences-region us-central1
```

---

## Option C: Google Batch

Google Batch is the newer, more scalable option for large workloads.

### Step 1: Enable Batch API

```bash
gcloud services enable batch.googleapis.com
```

### Step 2: Create Snakemake Profile

Create `~/.config/snakemake/googlebatch/config.yaml`:

```yaml
executor: googlebatch
googlebatch-region: us-central1
googlebatch-project: splice-neoepitope-project
default-resources:
  - googlebatch_boot_disk_size=50
  - googlebatch_machine_type=n2-standard-4
use-conda: true
```

### Step 3: Run with Google Batch

```bash
snakemake --profile googlebatch
```

---

## Data Transfer

### Uploading Large Datasets

For datasets > 10 GB, use `gsutil` with parallel uploads:

```bash
# Enable parallel composite uploads
gsutil -o GSUtil:parallel_composite_upload_threshold=150M \
    cp -r ./large_data_folder gs://your-bucket/
```

### Downloading from SRA on GCP

The SRA maintains a mirror on Google Cloud for faster downloads:

```bash
# Use the GCP SRA mirror
prefetch --location GCP SRR10971381
```

---

## Automated GPU Pipeline (TCRdock)

The `scripts/run_cloud_gpu.sh` script automates the full three-phase lifecycle
for running the pipeline with TCRdock structural validation:

| Phase | VM | What happens |
|-------|-----|-------------|
| **1 — CPU** | `neoepitope-predict-cpu` | Steps 1–5 (alignment → MHCflurry) |
| **2 — Copy** | GCS bucket | Results uploaded; CPU VM stopped |
| **3 — GPU** | `mhc-p-tcr-structure-spot-gpu` (Spot T4) | TCRdock + report; results uploaded; VM auto-stops |

### Quick start

```bash
# Test run (chr22, ~1 hour)
bash scripts/run_cloud_gpu.sh --mode test

# Production run
bash scripts/run_cloud_gpu.sh --mode prod

# Detached mode — closes laptop safely, orchestrator VM manages everything
bash scripts/run_cloud_gpu.sh --mode test --detach

# Specify branch and zone
bash scripts/run_cloud_gpu.sh --mode prod --branch main --zone europe-west1-b
```

### Retrieving results

The bucket name is derived from your GCP project ID (`<PROJECT_ID>-tcrdock-handoff`):

```bash
# Test
gcloud storage cp -r gs://<PROJECT_ID>-tcrdock-handoff/results/test/reports ./tcrdock_report
open tcrdock_report/local/report.html

# Production
gcloud storage cp -r gs://<PROJECT_ID>-tcrdock-handoff/results/reports ./tcrdock_report
open tcrdock_report/local/report.html
```

### How it works

1. **CPU VM** is created (or started) automatically. `setup_cloud.sh` installs
   conda + Snakemake. The pipeline runs steps 1–5 in a tmux session while the
   script polls `pipeline.log` for completion.
2. Results are uploaded to `gs://<PROJECT_ID>-tcrdock-handoff/` and the CPU VM is stopped.
3. A **GPU Spot VM** (n1-standard-4 + NVIDIA T4) is created.
   `setup_tcrdock_vm.sh` builds the TCRdock Docker image (~25 GB: CUDA 11.8,
   AlphaFold params, BLAST). Results are downloaded from GCS, TCRdock runs,
   and the HTML report is regenerated with the embedded Mol* 3D viewer. Final
   results are uploaded to GCS and the VM auto-stops.

### Detached mode

With `--detach`, the script launches a tiny `e2-micro` orchestrator VM that
runs the three phases via internal networking. You can close your laptop — the
orchestrator manages everything and shuts itself down when done.

### Cost estimate

| Component | Duration | Hourly rate | Total |
|-----------|----------|-------------|-------|
| CPU VM (n1-standard-4) | ~1–3 hr | ~$0.19 | ~$0.20–0.60 |
| GPU Spot VM (n1-standard-4 + T4) | ~30 min | ~$0.35 (Spot) | ~$0.20 |
| GCS bucket storage | — | negligible | — |
| Orchestrator (e2-micro, detach only) | ~2–4 hr | ~$0.01 | ~$0.02 |
| **Total (test run)** | | | **~$0.40–0.80** |

---

## Troubleshooting

### VM runs out of memory

```bash
# Check memory usage
free -h

# Solutions:
# 1. Switch to HISAT2 (requires only 8 GB)
# Edit config/config.yaml:
#   local_samples:
#     aligner: "hisat2"

# 2. Use a larger machine type
gcloud compute instances set-machine-type splice-pipeline \
    --machine-type=n2-highmem-8 \
    --zone=us-central1-a
```

### Pipeline crashes with "disk full"

```bash
# Check disk usage
df -h

# Resize disk (VM must be stopped first)
gcloud compute instances stop splice-pipeline --zone=us-central1-a
gcloud compute disks resize splice-pipeline --size=200GB --zone=us-central1-a
gcloud compute instances start splice-pipeline --zone=us-central1-a
```

### SSH connection drops during long runs

Use `tmux` or `screen` to keep sessions alive:

```bash
# Install tmux
sudo apt install tmux

# Start a tmux session
tmux new -s pipeline

# Run your pipeline
snakemake --cores $(nproc) --use-conda

# Detach: Ctrl+B, then D
# Reattach later:
tmux attach -t pipeline
```

### Spot VM was preempted

Spot VMs can be terminated at any time. To handle this:

1. Use Snakemake's built-in checkpointing (jobs restart from last checkpoint)
2. Use regular VMs for the final/critical stages
3. Save intermediate results to Cloud Storage periodically

---

## Quick Reference

```bash
# Create VM (HISAT2, 16 GB RAM)
gcloud compute instances create splice-pipeline \
    --zone=us-central1-a \
    --machine-type=n2-standard-4 \
    --boot-disk-size=100GB \
    --image-family=ubuntu-2204-lts \
    --image-project=ubuntu-os-cloud

# Connect
gcloud compute ssh splice-pipeline --zone=us-central1-a

# Stop (saves money, keeps disk)
gcloud compute instances stop splice-pipeline --zone=us-central1-a

# Start
gcloud compute instances start splice-pipeline --zone=us-central1-a

# Delete (removes everything)
gcloud compute instances delete splice-pipeline --zone=us-central1-a

# Transfer files TO VM
gcloud compute scp localfile.txt splice-pipeline:~/path/ --zone=us-central1-a

# Transfer files FROM VM
gcloud compute scp splice-pipeline:~/path/file.txt ./ --zone=us-central1-a
```

---

## Example: Full Workflow

Here's a complete example from start to finish:

```bash
# 1. Create VM
gcloud compute instances create splice-pipeline \
    --zone=us-central1-a \
    --machine-type=n2-standard-4 \
    --boot-disk-size=100GB \
    --image-family=ubuntu-2204-lts \
    --image-project=ubuntu-os-cloud

# 2. Connect
gcloud compute ssh splice-pipeline --zone=us-central1-a

# 3. Setup (run on VM)
sudo apt update && sudo apt install -y git curl wget
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
rm Miniforge3-Linux-x86_64.sh
~/miniforge3/bin/conda init bash
source ~/.bashrc
conda install -n base -c conda-forge -c bioconda "snakemake>=7.0,<9" python=3.11 samtools -y
conda install -c bioconda sra-tools=3.1.1 -y

# 4. Clone pipeline
git clone https://github.com/Jin-HoMLee/splice-neoepitope-pipeline.git
cd splice-neoepitope-pipeline

# 5. Download reference files
cd resources
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
samtools faidx GRCh38.primary_assembly.genome.fa
cd ..

# 6. Download RNA-Seq data from SRA
mkdir -p data
fasterq-dump SRR10971381 --split-files --outdir data/
mv data/SRR10971381_1.fastq data/tumor_01_R1.fastq
mv data/SRR10971381_2.fastq data/tumor_01_R2.fastq

# 7. Configure samples
cat > config/samples.tsv << 'EOF'
sample_id	sample_type	fastq1	fastq2
tumor_01	Primary Tumor	data/tumor_01_R1.fastq	data/tumor_01_R2.fastq
EOF

# 8. Run pipeline (VM will NOT auto-stop — remember to delete it when done)
snakemake --cores $(nproc) --use-conda --rerun-triggers mtime 2>&1 | tee pipeline.log

# 9. Exit VM (Ctrl+D or type 'exit')

# 10. Download results (from LOCAL machine)
gcloud compute scp --recurse splice-pipeline:~/splice-neoepitope-pipeline/results/ ./gcp_results/ --zone=us-central1-a

# 11. Delete VM when done to avoid charges
gcloud compute instances delete splice-pipeline --zone=us-central1-a

# Or: use the automated workflow instead (handles VM lifecycle + TCRdock):
# bash scripts/run_cloud_gpu.sh --mode prod --detach
```

---

## See Also

- [Google Cloud Compute Engine Documentation](https://cloud.google.com/compute/docs)
- [Snakemake Cloud Execution](https://snakemake.readthedocs.io/en/stable/executing/cloud.html)
- [GCP Free Tier](https://cloud.google.com/free)
