#!/usr/bin/env bash
# =============================================================================
# setup_tcrdock_vm.sh — Set up TCRdock + AlphaFold v2 on a GCP GPU VM
# =============================================================================
#
# Run this script once on a fresh GPU VM (Linux x64, NVIDIA T4) to install:
#   1. NVIDIA drivers + CUDA toolkit
#   2. Miniconda (if not already present)
#   3. TCRdock and its dependencies
#   4. AlphaFold v2 parameters (~3.5 GB, CC BY 4.0)
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
# Step 1 — NVIDIA driver + CUDA
# ---------------------------------------------------------------------------
log "Step 1: Checking NVIDIA driver..."

if ! command -v nvidia-smi &>/dev/null; then
    log "Installing NVIDIA driver and CUDA toolkit..."
    # Install CUDA 11.8 (compatible with AlphaFold v2 + TCRdock)
    wget -q https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-keyring_1.0-1_all.deb
    sudo dpkg -i cuda-keyring_1.0-1_all.deb
    sudo apt-get update -q
    sudo apt-get install -y cuda-11-8 2>&1 | tail -5
    rm cuda-keyring_1.0-1_all.deb
    log "NVIDIA driver installed. A reboot may be required if this is a first install."
else
    log "NVIDIA driver already present: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
fi

# ---------------------------------------------------------------------------
# Step 2 — Miniconda
# ---------------------------------------------------------------------------
log "Step 2: Checking Miniconda..."

if ! command -v conda &>/dev/null; then
    log "Installing Miniconda..."
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p "${HOME}/miniconda3"
    rm /tmp/miniconda.sh
    eval "$("${HOME}/miniconda3/bin/conda" shell.bash hook)"
    conda init bash
    log "Miniconda installed at ${HOME}/miniconda3"
else
    log "Miniconda already present: $(conda --version)"
    eval "$(conda shell.bash hook)"
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
log "Installing AlphaFold dependencies..."
conda run -n "${CONDA_ENV}" pip install \
    "jax[cuda11_pip]==0.3.25" \
    -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html \
    -q
conda run -n "${CONDA_ENV}" pip install \
    dm-haiku==0.0.9 \
    dm-tree==0.1.8 \
    ml-collections==0.1.1 \
    biopython==1.79 \
    -q

log "TCRdock environment ready."

# ---------------------------------------------------------------------------
# Step 4 — AlphaFold v2 parameters
# ---------------------------------------------------------------------------
log "Step 4: Downloading AlphaFold v2 parameters (~3.5 GB)..."

mkdir -p "${AF_PARAMS_DIR}"

AF_PARAMS_URL="https://storage.googleapis.com/alphafold/alphafold_params_colab_2022-12-06.tar"

if [[ ! -f "${AF_PARAMS_DIR}/params_model_1_multimer_v3.npz" ]]; then
    log "Downloading AlphaFold params (this may take 10-20 min)..."
    wget -q --show-progress \
        "${AF_PARAMS_URL}" \
        -O /tmp/af_params.tar
    log "Extracting AlphaFold params..."
    tar -xf /tmp/af_params.tar -C "${AF_PARAMS_DIR}"
    rm /tmp/af_params.tar
    log "AlphaFold params extracted to ${AF_PARAMS_DIR}"
else
    log "AlphaFold params already present at ${AF_PARAMS_DIR}"
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
