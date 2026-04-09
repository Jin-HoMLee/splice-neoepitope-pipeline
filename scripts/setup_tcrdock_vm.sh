#!/usr/bin/env bash
# =============================================================================
# setup_tcrdock_vm.sh — Set up TCRdock via Docker on a GCP GPU VM
# =============================================================================
#
# Installs Docker + NVIDIA Container Toolkit, then builds the official
# TCRdock Docker image (CUDA 11.8 + cuDNN 8 + Python 3.10 + JAX 0.3.25).
# Running a CUDA 11.8 container on a host with a newer driver (e.g. 12.8)
# is supported by NVIDIA's forward-compatibility guarantee.
#
# The image bundles TCRdock, AlphaFold params, and BLAST — no host-side
# dependency management required.
#
# Usage:
#   bash scripts/setup_tcrdock_vm.sh
#
# Expected runtime: ~30-45 min (Docker image build + param downloads)
# Expected disk:    ~20 GB (base image + pip deps + AlphaFold params)
# =============================================================================

set -euo pipefail

DOCKER_IMAGE="tcrdock:latest"
TCRDOCK_DIR="${HOME}/tcrdock"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ---------------------------------------------------------------------------
# Step 1 — Verify GPU
# ---------------------------------------------------------------------------
log "Step 1: Checking GPU..."
if command -v nvidia-smi &>/dev/null; then
    log "  GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
    log "  Driver: $(nvidia-smi --query-gpu=driver_version --format=csv,noheader 2>/dev/null | head -1)"
else
    echo "ERROR: nvidia-smi not found." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Step 2 — Docker
# ---------------------------------------------------------------------------
log "Step 2: Checking Docker..."

if ! command -v docker &>/dev/null; then
    log "  Installing Docker..."
    curl -fsSL https://get.docker.com | sh
    sudo usermod -aG docker "${USER}"
    log "  Docker installed."
else
    log "  Docker already installed: $(docker --version)"
fi

# Configure NVIDIA Container Runtime if not already done
if ! sudo docker info 2>/dev/null | grep -q "nvidia"; then
    log "  Configuring NVIDIA Container Runtime..."
    sudo nvidia-ctk runtime configure --runtime=docker
    sudo systemctl restart docker
fi

log "  GPU passthrough: $(sudo docker run --rm --gpus all nvidia/cuda:11.8.0-base-ubuntu22.04 \
    nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"

# ---------------------------------------------------------------------------
# Step 3 — TCRdock repo (for the Dockerfile)
# ---------------------------------------------------------------------------
log "Step 3: Checking TCRdock repo..."

if [[ ! -d "${TCRDOCK_DIR}" ]]; then
    git clone https://github.com/phbradley/TCRdock.git "${TCRDOCK_DIR}"
    log "  Cloned to ${TCRDOCK_DIR}"
else
    git -C "${TCRDOCK_DIR}" pull
    log "  Updated at ${TCRDOCK_DIR}"
fi

# ---------------------------------------------------------------------------
# Step 4 — Build Docker image
# ---------------------------------------------------------------------------
log "Step 4: Building TCRdock Docker image (this takes ~30-45 min)..."
log "  The image includes CUDA 11.8 + cuDNN 8, Python 3.10, JAX 0.3.25,"
log "  AlphaFold params, and BLAST — no host-side deps needed."

if sudo docker image inspect "${DOCKER_IMAGE}" &>/dev/null; then
    log "  Image '${DOCKER_IMAGE}' already exists — skipping build."
else
    sudo docker build -t "${DOCKER_IMAGE}" -f "${TCRDOCK_DIR}/docker/Dockerfile" "${TCRDOCK_DIR}"
    log "  Image built: ${DOCKER_IMAGE}"
fi

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
log "================================================================"
log "Setup complete!"
log ""
log "TCRdock Docker image: ${DOCKER_IMAGE}"
log ""
log "Run the pipeline with TCRdock enabled:"
log ""
log "  snakemake --cores \$(nproc) --use-conda --rerun-triggers mtime \\"
log "      --configfile config/config.yaml \\"
log "      --configfile config/tcrdock_gpu.yaml"
log "================================================================"
