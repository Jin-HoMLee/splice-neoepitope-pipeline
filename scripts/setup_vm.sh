#!/usr/bin/env bash
# =============================================================================
# setup_vm.sh — One-time setup for the consolidated pipeline + TCRdock VM
# =============================================================================
#
# What this does:
#   1. Installs system dependencies
#   2. Verifies GPU (nvidia-smi required)
#   3. Installs Docker + NVIDIA Container Runtime
#   4. Installs Miniforge3 (conda)
#   5. Creates the 'snakemake' conda environment
#   6. Clones / updates the pipeline repository
#   7. Builds the TCRdock Docker image (~20-30 min)
#   8. Downloads GRCh38 reference FASTA + GENCODE v47 GTF (~20-60 min)
#
# Usage (fresh VM — run from home directory):
#   curl -fsSL https://raw.githubusercontent.com/Jin-HoMLee/splice-neoepitope-pipeline/main/scripts/setup_vm.sh \
#       -o setup_vm.sh && bash setup_vm.sh [--repo-branch <branch>]
#
# After the repo is already cloned:
#   bash scripts/setup_vm.sh [--repo-branch <branch>]
#
# Expected runtime: ~60-90 min (Docker image build + reference download)
# Expected disk:    ~30 GB (image + params + reference data)
# =============================================================================

set -euo pipefail

REPO_URL="https://github.com/Jin-HoMLee/splice-neoepitope-pipeline.git"
REPO_DIR="$HOME/splice-neoepitope-pipeline"
REPO_BRANCH="main"
STANDALONE=true

while [[ $# -gt 0 ]]; do
    case $1 in
        --repo-branch) REPO_BRANCH="${2:?--repo-branch requires a value}"; shift 2 ;;
        --no-next-steps) STANDALONE=false; shift ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

log "=== Splice Neoepitope Pipeline — VM Setup ==="
log ""

# ---------------------------------------------------------------------------
# 1. System dependencies
# ---------------------------------------------------------------------------
log "[1/8] Installing system dependencies..."
sudo apt-get update -qq
sudo apt-get install -y git curl wget tmux samtools
log "    Done."

# ---------------------------------------------------------------------------
# 2. GPU check
# ---------------------------------------------------------------------------
log ""
log "[2/8] Checking GPU..."
if command -v nvidia-smi &>/dev/null; then
    log "  GPU:    $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
    log "  Driver: $(nvidia-smi --query-gpu=driver_version --format=csv,noheader 2>/dev/null | head -1)"
else
    echo "ERROR: nvidia-smi not found — this VM requires NVIDIA drivers (use the Deep Learning image)." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# 3. Docker + NVIDIA Container Runtime
# ---------------------------------------------------------------------------
log ""
log "[3/8] Checking Docker..."
if ! command -v docker &>/dev/null; then
    log "  Installing Docker..."
    curl -fsSL https://get.docker.com | sh
    sudo usermod -aG docker "${USER}"
    log "  Docker installed."
else
    log "  Docker already installed: $(docker --version)"
fi

if ! sudo docker info 2>/dev/null | grep -q "nvidia"; then
    log "  Configuring NVIDIA Container Runtime..."
    sudo nvidia-ctk runtime configure --runtime=docker
    sudo systemctl restart docker
fi

log "  GPU passthrough: $(sudo docker run --rm --gpus all nvidia/cuda:11.8.0-base-ubuntu22.04 \
    nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"

# ---------------------------------------------------------------------------
# 4. Miniforge3
# ---------------------------------------------------------------------------
log ""
if [[ -f "$HOME/miniforge3/bin/conda" ]]; then
    log "[4/8] conda already installed — skipping."
else
    log "[4/8] Installing Miniforge3..."
    MINIFORGE_INSTALLER="/tmp/Miniforge3.sh"
    curl --fail --location --retry 3 --retry-delay 2 --progress-bar \
        "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" \
        -o "$MINIFORGE_INSTALLER"
    if [[ ! -s "$MINIFORGE_INSTALLER" ]]; then
        echo "Error: Miniforge installer download failed." >&2
        rm -f "$MINIFORGE_INSTALLER"; exit 1
    fi
    bash "$MINIFORGE_INSTALLER" -b -p "$HOME/miniforge3"
    rm "$MINIFORGE_INSTALLER"
    # shellcheck source=/dev/null
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
    "$HOME/miniforge3/bin/conda" init bash
    log "  Miniforge3 installed at: $HOME/miniforge3"
fi

if [[ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]]; then
    # shellcheck source=/dev/null
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
fi

# Strict channel priority prevents conda-forge/defaults mixing that causes
# subtle env conflicts. Required for robust Snakemake conda environments.
conda config --set channel_priority strict

# ---------------------------------------------------------------------------
# 5. snakemake conda env
# ---------------------------------------------------------------------------
log ""
if [[ -d "$HOME/miniforge3/envs/snakemake" ]]; then
    log "[5/8] 'snakemake' environment already exists — skipping."
else
    log "[5/8] Creating 'snakemake' conda environment..."
    conda create -n snakemake -c conda-forge -c bioconda \
        "snakemake>=8.0,<9" python=3.11 -y
    log "  Created 'snakemake' environment."
fi

# ---------------------------------------------------------------------------
# 6. Clone / update repo
# ---------------------------------------------------------------------------
log ""
if [[ -d "$REPO_DIR/.git" ]]; then
    log "[6/8] Repository already cloned — pulling branch: ${REPO_BRANCH}..."
    git -C "$REPO_DIR" fetch origin
    git -C "$REPO_DIR" checkout "$REPO_BRANCH"
    git -C "$REPO_DIR" pull origin "$REPO_BRANCH"
else
    log "[6/8] Cloning repository (branch: ${REPO_BRANCH})..."
    git clone --branch "$REPO_BRANCH" "$REPO_URL" "$REPO_DIR"
fi

# ---------------------------------------------------------------------------
# 7. Build TCRdock Docker image
# ---------------------------------------------------------------------------
log ""
log "[7/8] Building TCRdock Docker image (~20-30 min)..."
log "  Image bundles CUDA 11.8 + cuDNN 8, Python 3.10, JAX 0.3.25, AlphaFold params, BLAST."
log "  Running CUDA 11.8 inside the container on this host's newer driver is supported"
log "  by NVIDIA's forward-compatibility guarantee."

if sudo docker image inspect "tcrdock:latest" &>/dev/null; then
    log "  Image 'tcrdock:latest' already exists — skipping."
else
    sudo docker build -t "tcrdock:latest" \
        -f "${REPO_DIR}/docker/Dockerfile.pipeline" \
        "${REPO_DIR}"
    log "  Image built: tcrdock:latest"
fi

# ---------------------------------------------------------------------------
# 8. Reference data
# ---------------------------------------------------------------------------
log ""
log "[8/8] Downloading reference data (20–60 min)..."

RESOURCES="$REPO_DIR/resources"
mkdir -p "$RESOURCES"

FASTA_GZ="$RESOURCES/GRCh38.primary_assembly.genome.fa.gz"
FASTA="$RESOURCES/GRCh38.primary_assembly.genome.fa"

if [[ -f "$FASTA" ]]; then
    log "  GRCh38 FASTA already exists — skipping."
else
    log "  Downloading GRCh38 primary assembly FASTA (~830 MB compressed)..."
    curl --fail -L --progress-bar \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz" \
        -o "$FASTA_GZ"
    log "  Decompressing FASTA (~3.1 GB uncompressed)..."
    gunzip "$FASTA_GZ"
    [[ -s "$FASTA" ]] || { echo "ERROR: FASTA decompression failed." >&2; exit 1; }
    log "  Saved: $FASTA"
fi

FAI="$FASTA.fai"
if [[ -f "$FAI" ]]; then
    log "  FASTA index (.fai) already exists — skipping."
else
    log "  Indexing FASTA with samtools faidx..."
    samtools faidx "$FASTA"
    log "  Saved: $FAI"
fi

GTF="$RESOURCES/gencode.v47.annotation.gtf.gz"
if [[ -f "$GTF" ]]; then
    log "  GENCODE GTF already exists — skipping."
else
    log "  Downloading GENCODE v47 annotation GTF (~50 MB)..."
    curl --fail -L --progress-bar \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz" \
        -o "$GTF"
    gzip -t "$GTF" || { echo "ERROR: GTF download appears corrupt." >&2; exit 1; }
    log "  Saved: $GTF"
fi

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
log ""
log "=== Setup complete! ==="

if [[ "${STANDALONE}" == true ]]; then
    log ""
    log "Next steps:"
    log "  0. Open a new shell (or: source ~/.bashrc) so 'conda' is on PATH"
    log "  1. Edit config/samples.tsv with your sample FASTQ paths"
    log "  2. conda activate snakemake"
    log "  3. From a tmux session, run the full pipeline:"
    log "       cd $REPO_DIR"
    log "       snakemake --cores \$(nproc) --use-conda --rerun-triggers mtime \\"
    log "           --configfile config/config.yaml config/gpu.yaml \\"
    log "           2>&1 | tee pipeline.log"
fi
