#!/usr/bin/env bash
# =============================================================================
# setup_tcrdock_vm.sh — Set up TCRdock via Docker on a GCP GPU VM
# =============================================================================
#
# Installs Docker + NVIDIA Container Toolkit, then builds the pipeline
# TCRdock Docker image (CUDA 11.8 + cuDNN 8 + Python 3.10 + JAX 0.3.25).
# Running a CUDA 11.8 container on a host with a newer driver (e.g. 12.8)
# is supported by NVIDIA's forward-compatibility guarantee.
#
# Uses docker/Dockerfile.pipeline (omits openmm/pdbfixer — not needed for
# structure prediction, only for optional Amber relaxation).
#
# The image bundles TCRdock, AlphaFold params, and BLAST — no host-side
# dependency management required.
#
# Usage:
#   bash scripts/setup_tcrdock_vm.sh [config/config.yaml]
#
# Expected runtime: ~20-30 min (Docker image build + param downloads)
# Expected disk:    ~25 GB (base image + pip deps + AlphaFold params)
# =============================================================================

set -euo pipefail

DOCKER_IMAGE="tcrdock:latest"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BASE_CONFIG="${1:-config/config.yaml}"   # caller can override, e.g. config/test_config.yaml

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ---------------------------------------------------------------------------
# Step 0.5 — Install NumPy in gcloud SDK Python to silence IAP tunnel warnings
# ---------------------------------------------------------------------------
log "Step 0.5: Installing NumPy in gcloud SDK Python..."
GCLOUD_ROOT="$(gcloud info --format='value(installation.sdk_root)' 2>/dev/null || true)"
GCLOUD_PYTHON=""

# Try standard location first (e.g., /opt/google-cloud-sdk)
if [[ -n "${GCLOUD_ROOT}" && -x "${GCLOUD_ROOT}/bin/python3" ]]; then
    GCLOUD_PYTHON="${GCLOUD_ROOT}/bin/python3"
# Try snap location (e.g., /snap/google-cloud-cli/XXX)
elif [[ -n "${GCLOUD_ROOT}" && -x "${GCLOUD_ROOT}/platform/bundledpythonunix/bin/python3" ]]; then
    GCLOUD_PYTHON="${GCLOUD_ROOT}/platform/bundledpythonunix/bin/python3"
fi

if [[ -n "${GCLOUD_PYTHON}" ]]; then
    "${GCLOUD_PYTHON}" -m pip install --quiet numpy 2>/dev/null || \
        log "  Note: NumPy installation failed (non-critical; tunnel will still work)."
else
    log "  gcloud SDK Python not found (non-critical; will proceed anyway)."
fi

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
# Step 3 — Build Docker image
# ---------------------------------------------------------------------------
log "Step 3: Building TCRdock Docker image (this takes ~20-30 min)..."
log "  Using docker/Dockerfile.pipeline (no openmm/pdbfixer — not needed"
log "  for prediction, only for optional Amber relaxation)."
log "  The image includes CUDA 11.8 + cuDNN 8, Python 3.10, JAX 0.3.25,"
log "  AlphaFold params, and BLAST — no host-side deps needed."

if sudo docker image inspect "${DOCKER_IMAGE}" &>/dev/null; then
    log "  Image '${DOCKER_IMAGE}' already exists — skipping build."
else
    sudo docker build -t "${DOCKER_IMAGE}" \
        -f "${REPO_DIR}/docker/Dockerfile.pipeline" \
        "${REPO_DIR}"
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
log "      --configfile ${BASE_CONFIG} config/tcrdock_gpu.yaml"
log "================================================================"
