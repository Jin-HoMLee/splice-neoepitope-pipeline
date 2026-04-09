#!/usr/bin/env bash
# =============================================================================
# setup_tcrdock_vm.sh — Set up TCRdock + AlphaFold v2 on a GCP GPU VM
# =============================================================================
#
# Designed for use with GCP Deep Learning VM images (common-cu128-ubuntu-2204)
# which have CUDA 12.8 + NVIDIA drivers pre-installed. Skips driver setup.
#
# Run this script once on a fresh GPU VM (Linux x64, NVIDIA T4) to install:
#   1. Miniconda (if not already present)
#   2. TCRdock and its dependencies
#   3. AlphaFold v2 parameters (~3.5 GB, CC BY 4.0)
#
# After this script completes, update config/config.yaml:
#   tcrdock:
#     enabled: true
#     tcrdock_dir: "$TCRDOCK_DIR"
#     alphafold_params_dir: "$AF_PARAMS_DIR"
#
# Then run the pipeline:
#   snakemake --cores $(nproc) --use-conda --rerun-triggers mtime
#
# Usage:
#   bash scripts/setup_tcrdock_vm.sh
#
# Expected runtime: ~20-30 min (AlphaFold weights download dominates)
# Expected disk: ~6 GB (TCRdock + weights)
# =============================================================================

set -euo pipefail

TCRDOCK_DIR="${HOME}/tcrdock"
AF_PARAMS_DIR="${HOME}/af_params"
CONDA_ENV="tcrdock"
PYTHON_VERSION="3.8"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ---------------------------------------------------------------------------
# Step 1 — Verify GPU is visible (CUDA pre-installed on Deep Learning VM image)
# ---------------------------------------------------------------------------
log "Step 1: Checking GPU..."
if command -v nvidia-smi &>/dev/null; then
    log "  GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
    log "  CUDA: $(nvcc --version 2>/dev/null | grep release | awk '{print $6}' | tr -d ',')"
else
    echo "ERROR: nvidia-smi not found. Use a Deep Learning VM image with CUDA pre-installed." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Step 2 — Conda (prefer Miniforge already installed by run_cloud_gpu.sh)
# ---------------------------------------------------------------------------
log "Step 2: Checking conda..."

# Prefer Miniforge (installed by run_cloud_gpu.sh); fall back to installing it
if [[ -f "${HOME}/miniforge3/bin/conda" ]]; then
    log "  Using Miniforge at ${HOME}/miniforge3"
    source "${HOME}/miniforge3/etc/profile.d/conda.sh"
elif command -v conda &>/dev/null; then
    log "  Using conda: $(conda --version)"
    eval "$(conda shell.bash hook)"
else
    log "  Installing Miniforge3..."
    curl -fsSL https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
        -o /tmp/miniforge.sh
    bash /tmp/miniforge.sh -b -p "${HOME}/miniforge3"
    rm /tmp/miniforge.sh
    source "${HOME}/miniforge3/etc/profile.d/conda.sh"
    log "  Miniforge installed at ${HOME}/miniforge3"
fi

# ---------------------------------------------------------------------------
# Step 3 — TCRdock
# ---------------------------------------------------------------------------
log "Step 3: Installing TCRdock..."

if [[ ! -d "${TCRDOCK_DIR}" ]]; then
    git clone https://github.com/phbradley/TCRdock.git "${TCRDOCK_DIR}"
    log "TCRdock cloned to ${TCRDOCK_DIR}"
else
    log "TCRdock already present at ${TCRDOCK_DIR}, pulling latest..."
    git -C "${TCRDOCK_DIR}" pull
fi

# Create dedicated conda environment for TCRdock (Python 3.8 required)
if ! conda env list | grep -q "^${CONDA_ENV} "; then
    log "Creating conda environment '${CONDA_ENV}' (Python ${PYTHON_VERSION})..."
    conda create -y -n "${CONDA_ENV}" python="${PYTHON_VERSION}"
fi

log "Installing TCRdock Python dependencies..."
conda run -n "${CONDA_ENV}" pip install -r "${TCRDOCK_DIR}/requirements.txt" -q

log "Installing BLAST (required by TCRdock)..."
conda run -n "${CONDA_ENV}" python "${TCRDOCK_DIR}/download_blast.py"

# Install AlphaFold Python dependencies into the TCRdock environment
# JAX with CUDA 12 support (compatible with Deep Learning VM CUDA 12.8 image)
log "Installing AlphaFold dependencies..."
conda run -n "${CONDA_ENV}" pip install \
    "jax[cuda12_pip]" \
    -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html \
    -q
conda run -n "${CONDA_ENV}" pip install \
    dm-haiku \
    dm-tree \
    ml-collections \
    biopython==1.79 \
    -q

log "TCRdock environment ready."

# ---------------------------------------------------------------------------
# Step 4 — AlphaFold v2 parameters
# ---------------------------------------------------------------------------
log "Step 4: Downloading AlphaFold v2 parameters (~3.5 GB)..."

mkdir -p "${AF_PARAMS_DIR}"

AF_PARAMS_URL="https://storage.googleapis.com/alphafold/alphafold_params_colab_2022-12-06.tar"

if [[ ! -f "${AF_PARAMS_DIR}/params/params_model_1_multimer_v3.npz" ]]; then
    log "Downloading AlphaFold params (this may take 10-20 min)..."
    wget -q --show-progress \
        "${AF_PARAMS_URL}" \
        -O /tmp/af_params.tar
    log "Extracting AlphaFold params..."
    # AlphaFold expects data_dir/params/*.npz — extract into the params/ subdirectory.
    mkdir -p "${AF_PARAMS_DIR}/params"
    tar -xf /tmp/af_params.tar -C "${AF_PARAMS_DIR}/params"
    rm /tmp/af_params.tar
    log "AlphaFold params extracted to ${AF_PARAMS_DIR}/params"
else
    log "AlphaFold params already present at ${AF_PARAMS_DIR}/params"
fi

# ---------------------------------------------------------------------------
# Done — print config snippet
# ---------------------------------------------------------------------------
log "================================================================"
log "Setup complete!"
log ""
log "Update the paths in config/tcrdock_gpu.yaml if they differ:"
log "  tcrdock_dir:          ${TCRDOCK_DIR}"
log "  alphafold_params_dir: ${AF_PARAMS_DIR}"
log ""
log "Then run the pipeline with TCRdock enabled:"
log ""
log "  snakemake --cores \$(nproc) --use-conda --rerun-triggers mtime \\"
log "      --configfile config/config.yaml \\"
log "      --configfile config/tcrdock_gpu.yaml"
log ""
log "To run without TCRdock (CPU-only, local):"
log ""
log "  snakemake --cores \$(nproc) --use-conda --rerun-triggers mtime"
log "================================================================"
