#!/usr/bin/env bash
# =============================================================================
# test_run_cloud_gpu.sh — Test pipeline run: chr22 + 500K reads + TCRdock
# =============================================================================
#
# Same three-phase lifecycle as run_cloud_gpu.sh, but uses the test config
# (chr22, 500K reads, test_samples.tsv) so the run completes in ~1 hour
# instead of several hours. Use this to validate TCRdock and report changes
# before committing to a full production run.
#
# Uses:
#   config/test_config.yaml  — chr22 reference, results/test/ output paths
#   config/tcrdock_gpu.yaml  — TCRdock overlay (GPU VM only)
#   scripts/prepare_test_data.sh — downloads chr22 reference + 500K FASTQs
#
# Usage:
#   bash scripts/test_run_cloud_gpu.sh [--branch <branch>] [--zone <zone>]
#
# After the run, copy the report:
#   gcloud compute scp --tunnel-through-iap --recurse \
#       mhc-p-tcr-structure-spot-gpu:~/splice-neoepitope-pipeline/results/test/reports \
#       ./tcrdock_report --zone=europe-west1-b
#   gcloud compute instances delete mhc-p-tcr-structure-spot-gpu --zone=europe-west1-b --quiet
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
CPU_VM="neoepitope-predict-cpu"
GPU_VM="mhc-p-tcr-structure-spot-gpu"
ZONE="europe-west1-b"
MACHINE_TYPE="n1-standard-4"
ACCELERATOR="type=nvidia-tesla-t4,count=1"
DISK_SIZE="100GB"
IMAGE_FAMILY="common-cu128-ubuntu-2204-nvidia-570"  # CUDA 12.8 pre-installed
IMAGE_PROJECT="deeplearning-platform-release"
REPO_URL="https://github.com/Jin-HoMLee/splice-neoepitope-pipeline.git"

# GCS bucket for VM-to-VM result transfer (created automatically if absent)
GCS_BUCKET="tcrdock-handoff"
GCS_PATH="gs://${GCS_BUCKET}/results/test"

BRANCH="$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "main")"

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --branch) BRANCH="${2:?--branch requires a value}"; shift 2 ;;
        --zone)   ZONE="${2:?--zone requires a value}";   shift 2 ;;
        *)
            echo "Unknown argument: $1" >&2
            echo "Usage: $0 [--branch <branch>] [--zone <zone>]" >&2
            exit 1 ;;
    esac
done

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

ssh_cmd() {
    local vm="$1"; shift
    gcloud compute ssh "${vm}" --zone="${ZONE}" --tunnel-through-iap "$@"
}

wait_for_ssh() {
    local vm="$1"
    local max_wait=300
    local elapsed=0
    log "  Waiting for SSH on ${vm} (up to ${max_wait}s)..."
    until ssh_cmd "${vm}" --command="echo ssh-ok" --quiet 2>/dev/null; do
        if [[ ${elapsed} -ge ${max_wait} ]]; then
            echo "ERROR: SSH on ${vm} not available after ${max_wait}s." >&2; exit 1
        fi
        sleep 15; elapsed=$((elapsed + 15))
        log "  Still waiting... (${elapsed}s)"
    done
    log "  SSH ready on ${vm}."
}

vm_status() {
    gcloud compute instances describe "$1" --zone="${ZONE}" \
        --format="value(status)" 2>/dev/null || echo "NOT_FOUND"
}

log "================================================================"
log "Cloud GPU test pipeline run (chr22 + 500K reads)"
log "  CPU VM:  ${CPU_VM}"
log "  GPU VM:  ${GPU_VM}"
log "  Zone:    ${ZONE}"
log "  Branch:  ${BRANCH}"
log "  Bucket:  gs://${GCS_BUCKET}"
log "================================================================"

# ---------------------------------------------------------------------------
# Ensure GCS handoff bucket exists
# ---------------------------------------------------------------------------
if gsutil ls "gs://${GCS_BUCKET}" &>/dev/null; then
    log "GCS bucket gs://${GCS_BUCKET} already exists."
else
    log "Creating GCS bucket gs://${GCS_BUCKET}..."
    gsutil mb -l "${ZONE%-*}" "gs://${GCS_BUCKET}"
    log "  Bucket created."
fi

PROJECT_NUMBER="$(gcloud projects describe "$(gcloud config get-value project)" \
    --format='value(projectNumber)' 2>/dev/null)"
COMPUTE_SA="${PROJECT_NUMBER}-compute@developer.gserviceaccount.com"
log "Granting ${COMPUTE_SA} objectAdmin on gs://${GCS_BUCKET}..."
gsutil iam ch "serviceAccount:${COMPUTE_SA}:objectAdmin" "gs://${GCS_BUCKET}" 2>/dev/null || \
    log "  Warning: could not set bucket IAM (may already be set)."

# ===========================================================================
# Phase 1 — CPU pipeline (steps 1-5, test config)
# ===========================================================================
log ""
log "=== Phase 1: CPU pipeline — test config (chr22, 500K reads) ==="

CPU_STATUS="$(vm_status "${CPU_VM}")"
if [[ "${CPU_STATUS}" == "TERMINATED" ]]; then
    log "Starting ${CPU_VM}..."
    gcloud compute instances start "${CPU_VM}" --zone="${ZONE}"
elif [[ "${CPU_STATUS}" == "NOT_FOUND" ]]; then
    log "CPU VM ${CPU_VM} not found — creating it..."
    gcloud compute instances create "${CPU_VM}" \
        --zone="${ZONE}" \
        --machine-type="${MACHINE_TYPE}" \
        --image-family="ubuntu-2204-lts" \
        --image-project="ubuntu-os-cloud" \
        --boot-disk-size="50GB" \
        --boot-disk-type=pd-standard \
        --scopes=cloud-platform \
        --metadata=enable-oslogin=FALSE
fi

wait_for_ssh "${CPU_VM}"

log "Running setup_cloud.sh on ${CPU_VM}..."
ssh_cmd "${CPU_VM}" -- bash -s -- --repo-branch "${BRANCH}" < scripts/setup_cloud.sh
log "  Setup complete."

log "Downloading test data (chr22 + 500K reads) on ${CPU_VM}..."
ssh_cmd "${CPU_VM}" -- bash -c "cd \$HOME/splice-neoepitope-pipeline && bash -s" < scripts/prepare_test_data.sh
log "  Test data ready."

log "Pulling branch '${BRANCH}' and starting pipeline on ${CPU_VM}..."
ssh_cmd "${CPU_VM}" -- bash -s <<EOF
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"

git fetch --all --prune
git checkout "${BRANCH}"
git pull origin "${BRANCH}"

tmux kill-session -t pipeline 2>/dev/null || true

tmux new-session -d -s pipeline "
    source \"\$HOME/miniforge3/etc/profile.d/conda.sh\" 2>/dev/null || true
    conda activate snakemake
    snakemake \\
        --cores \$(nproc) \\
        --use-conda \\
        --rerun-triggers mtime \\
        --configfile config/test_config.yaml \\
        2>&1 | tee pipeline.log
    echo 'Pipeline finished.'
"

echo "Pipeline started in tmux session 'pipeline'."
EOF

log "Polling pipeline.log for completion (this may take ~1 hour)..."
log "  Monitor with:"
log "    gcloud compute ssh ${CPU_VM} --zone=${ZONE} --tunnel-through-iap -- tail -f splice-neoepitope-pipeline/pipeline.log"
while true; do
    DONE="$(ssh_cmd "${CPU_VM}" --command \
        "grep -cE 'steps \(100%\) done|Nothing to be done' \$HOME/splice-neoepitope-pipeline/pipeline.log 2>/dev/null || echo 0" \
        2>/dev/null | tail -1)"
    if [[ "${DONE}" -ge 1 ]]; then
        log "  Pipeline finished (100% done detected in log)."
        break
    fi
    FAILED="$(ssh_cmd "${CPU_VM}" --command \
        "grep -c 'Error\|Exiting because' \$HOME/splice-neoepitope-pipeline/pipeline.log 2>/dev/null || echo 0" \
        2>/dev/null | tail -1)"
    if [[ "${FAILED}" -ge 1 ]]; then
        echo "ERROR: Pipeline appears to have failed. Check pipeline.log on ${CPU_VM}." >&2
        exit 1
    fi
    PROGRESS="$(ssh_cmd "${CPU_VM}" --command \
        "tail -3 \$HOME/splice-neoepitope-pipeline/pipeline.log 2>/dev/null" \
        2>/dev/null | grep -v '^$' | sed 's/^/    /')"
    log "  Still running — checking again in 60s..."
    [[ -n "${PROGRESS}" ]] && echo "${PROGRESS}"
    sleep 60
done

# ===========================================================================
# Phase 2 — Copy results to GCS, stop CPU VM
# ===========================================================================
log ""
log "=== Phase 2: Copying results from ${CPU_VM} ==="

log "Uploading results/test/ from ${CPU_VM} to ${GCS_PATH}..."
ssh_cmd "${CPU_VM}" -- bash -s <<EOF
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"
if [[ -d "results/test" ]]; then
    gcloud storage cp -r "results/test" "${GCS_PATH%/test}"
    echo "Upload complete."
else
    echo "ERROR: results/test/ not found." >&2; exit 1
fi
EOF

log "Stopping ${CPU_VM} in background to save cost..."
gcloud compute instances stop "${CPU_VM}" --zone="${ZONE}" &

# ===========================================================================
# Phase 3 — GPU Spot VM: TCRdock structural validation
# ===========================================================================
log ""
log "=== Phase 3: GPU Spot VM — TCRdock ==="

GPU_STATUS="$(vm_status "${GPU_VM}")"
if [[ "${GPU_STATUS}" == "NOT_FOUND" ]]; then
    log "Creating GPU Spot VM ${GPU_VM}..."
    gcloud compute instances create "${GPU_VM}" \
        --zone="${ZONE}" \
        --machine-type="${MACHINE_TYPE}" \
        --accelerator="${ACCELERATOR}" \
        --maintenance-policy=TERMINATE \
        --provisioning-model=SPOT \
        --instance-termination-action=STOP \
        --image-family="${IMAGE_FAMILY}" \
        --image-project="${IMAGE_PROJECT}" \
        --boot-disk-size="${DISK_SIZE}" \
        --boot-disk-type=pd-ssd \
        --scopes=cloud-platform \
        --metadata=enable-oslogin=FALSE
    log "  GPU VM created."
elif [[ "${GPU_STATUS}" == "TERMINATED" ]]; then
    log "GPU VM ${GPU_VM} exists (stopped) — starting it..."
    gcloud compute instances start "${GPU_VM}" --zone="${ZONE}"
else
    log "GPU VM ${GPU_VM} already running (status: ${GPU_STATUS})."
fi

wait_for_ssh "${GPU_VM}"

log "Cloning / updating repo on GPU VM (branch: ${BRANCH})..."
ssh_cmd "${GPU_VM}" -- bash -s <<EOF
set -euo pipefail
REPO_DIR="\$HOME/splice-neoepitope-pipeline"
if [[ -d "\${REPO_DIR}/.git" ]]; then
    git -C "\${REPO_DIR}" fetch --all --prune
    git -C "\${REPO_DIR}" checkout "${BRANCH}"
    git -C "\${REPO_DIR}" pull origin "${BRANCH}"
else
    git clone --branch "${BRANCH}" "${REPO_URL}" "\${REPO_DIR}"
fi
echo "  Branch: \$(git -C \${REPO_DIR} rev-parse --abbrev-ref HEAD)"
EOF

log "Installing conda and Snakemake on GPU VM..."
ssh_cmd "${GPU_VM}" -- bash -s <<'EOF'
set -euo pipefail

if [[ ! -f "$HOME/miniforge3/bin/conda" ]]; then
    echo "  Installing Miniforge3..."
    curl -fsSL https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
        -o /tmp/miniforge.sh
    bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"
    rm /tmp/miniforge.sh
fi

source "$HOME/miniforge3/etc/profile.d/conda.sh"

if [[ ! -d "$HOME/miniforge3/envs/snakemake" ]]; then
    echo "  Creating snakemake conda environment..."
    conda create -n snakemake -c conda-forge -c bioconda \
        "snakemake>=8.0,<9" python=3.11 -y
fi
echo "  Snakemake environment ready."
EOF

log "Running setup_tcrdock_vm.sh on GPU VM..."
log "  This downloads ~3.5 GB of AlphaFold params — may take 15-20 min."
ssh_cmd "${GPU_VM}" -- bash -s <<'EOF'
set -euo pipefail
cd "$HOME/splice-neoepitope-pipeline"
bash scripts/setup_tcrdock_vm.sh
EOF

log "Downloading results/test/ from ${GCS_PATH} onto GPU VM..."
ssh_cmd "${GPU_VM}" -- bash -s <<EOF
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"
mkdir -p results
gcloud storage cp -r "${GCS_PATH}" results/
echo "Download complete."
EOF

log "Running TCRdock step on GPU VM..."
log "  The VM will shut down automatically when finished."
log "  Monitor with:"
log "    gcloud compute ssh ${GPU_VM} --zone=${ZONE} --tunnel-through-iap -- tail -f splice-neoepitope-pipeline/tcrdock.log"

ssh_cmd "${GPU_VM}" -- bash -s <<'EOF'
set -euo pipefail
cd "$HOME/splice-neoepitope-pipeline"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate snakemake

snakemake \
    --cores "$(nproc)" \
    --use-conda \
    --rerun-triggers mtime \
    --configfile config/test_config.yaml \
    --configfile config/tcrdock_gpu.yaml \
    2>&1 | tee tcrdock.log

echo "TCRdock finished — shutting down GPU VM..."
sudo shutdown -h now
EOF

# ===========================================================================
# Done
# ===========================================================================
log "================================================================"
log "All phases complete. GPU VM will stop shortly."
log ""
log "To retrieve the report:"
log "  Start the VM first if needed:"
log "    gcloud compute instances start ${GPU_VM} --zone=${ZONE}"
log "  Then copy:"
log "    gcloud compute scp --tunnel-through-iap --recurse \\"
log "        ${GPU_VM}:~/splice-neoepitope-pipeline/results/test/reports \\"
log "        ./tcrdock_report --zone=${ZONE}"
log ""
log "To delete the GPU VM and stop disk charges:"
log "    gcloud compute instances delete ${GPU_VM} --zone=${ZONE} --quiet"
log "================================================================"
