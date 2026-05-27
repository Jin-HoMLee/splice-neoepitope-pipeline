#!/usr/bin/env bash
# =============================================================================
# run_cloud_gpu.sh — Single-VM pipeline: alignment → MHCflurry (GPU) → TCRdock → report
# =============================================================================
#
# All pipeline steps run on one VM (neoepitope-pipeline):
#   alignment → HLA typing → junction filtering → translation →
#   proteome filter → MHCflurry (GPU) → TCRdock (Docker/GPU) → report
#
# Results are uploaded to GCS when the pipeline finishes, then the VM stops.
#
# Modes:
#   --mode prod   Full genome, production config (default)
#   --mode test   chr22, 500K reads — fast validation
#
# Detached mode (--detach):
#   Launches a tiny e2-micro "orchestrator" VM that runs this script in a
#   tmux session. You can close your laptop — the orchestrator manages the
#   pipeline VM and auto-stops when done.
#
# Requirements:
#   - gcloud CLI authenticated: gcloud auth login
#   - Active project: gcloud config set project <PROJECT_ID>
#   - P100 GPU quota in the target zone (NVIDIA_P100_GPUS >= 1)
#
# Usage:
#   bash scripts/run_cloud_gpu.sh --samples config/samples/patient_002.tsv \
#       [--mode test|prod] [--branch <branch>] [--zone <zone>] [--detach] [--keep-vm]
#
# After the run, download the report:
#   gcloud storage cp -r gs://splice-neoepitope-project/results/<patient_id>/reports ./report
#   open report/report.html
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PIPELINE_VM="neoepitope-pipeline"
ZONE="europe-west4-a"
MACHINE_TYPE="n1-highmem-8"   # n1 required for P100; 52 GB RAM for OptiType (~36 GB peak)
ACCELERATOR="type=nvidia-tesla-p100,count=1"  # T4 quota is 0/0 in europe-west1; P100 standard quota (NVIDIA_P100_GPUS >= 1) required
DISK_SIZE="200GB"              # pipeline data + Docker image (~25 GB) + reference data
IMAGE_FAMILY="common-cu129-ubuntu-2204-nvidia-580"  # only available Ubuntu 22.04 DL image; driver downgraded to 570-server below for P100 (Pascal) compatibility
IMAGE_PROJECT="deeplearning-platform-release"
REPO_URL="https://github.com/Jin-HoMLee/splice-neoepitope-pipeline.git"
ORCH_VM="neoepitope-orchestrator"

GCS_BUCKET="splice-neoepitope-project"

BRANCH="$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "main")"
MODE="prod"
SAMPLES=""
DETACH=false
KEEP_VM=false
CLOUD_INTERNAL=false

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --mode)    MODE="${2:?--mode requires test or prod}";    shift 2 ;;
        --branch)  BRANCH="${2:?--branch requires a value}";    shift 2 ;;
        --zone)    ZONE="${2:?--zone requires a value}";        shift 2 ;;
        --samples) SAMPLES="${2:?--samples requires a path}";   shift 2 ;;
        --detach)  DETACH=true;                                  shift ;;
        --keep-vm) KEEP_VM=true;                                 shift ;;
        --_cloud-internal) CLOUD_INTERNAL=true;                  shift ;;
        *)
            echo "Unknown argument: $1" >&2
            echo "Usage: $0 [--mode test|prod] [--branch <branch>] [--zone <zone>] [--samples <path>] [--detach] [--keep-vm]" >&2
            exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
RESULTS_DIR="results"
GCS_PATH="gs://${GCS_BUCKET}/results"
GPU_CONFIG_FILE="config/gpu_config.yaml"

case "${MODE}" in
    test)
        CONFIG_FILE="config/test_config.yaml"
        PREPARE_DATA_SCRIPT="scripts/prepare_test_data.sh"
        [[ -z "${SAMPLES}" ]] && SAMPLES="config/samples/patient_001_test.tsv"
        ;;
    prod|production)
        CONFIG_FILE="config/config.yaml"
        PREPARE_DATA_SCRIPT=""
        [[ -z "${SAMPLES}" ]] && { echo "ERROR: --samples required for prod mode (e.g. --samples config/samples/patient_002.tsv)" >&2; exit 1; }
        ;;
    *)
        echo "ERROR: --mode must be 'test' or 'prod', got '${MODE}'" >&2
        exit 1 ;;
esac
PATIENT_ID=$(basename "${SAMPLES}" .tsv)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

ssh_cmd() {
    local vm="$1"; shift
    if [[ "${CLOUD_INTERNAL}" == true ]]; then
        gcloud compute ssh "${vm}" --zone="${ZONE}" --internal-ip \
            --ssh-key-expire-after=24h --quiet "$@"
    else
        gcloud compute ssh "${vm}" --zone="${ZONE}" --tunnel-through-iap "$@"
    fi
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
log "Cloud pipeline run (mode: ${MODE})"
log "  VM:      ${PIPELINE_VM}"
log "  Zone:    ${ZONE}"
log "  Branch:  ${BRANCH}"
log "  Config:  ${CONFIG_FILE}"
log "  Bucket:  gs://${GCS_BUCKET}"
[[ "${DETACH}" == true ]] && log "  Detach:  yes (orchestrator: ${ORCH_VM})"
[[ "${KEEP_VM}" == true ]] && log "  Keep VM: yes (VM will NOT be stopped on exit — remember to stop manually)"
log "================================================================"

# ===========================================================================
# Detached mode — offload to an orchestrator VM so you can close your laptop
# ===========================================================================
if [[ "${DETACH}" == true ]]; then
    log ""
    log "=== Detached mode: launching orchestrator VM ==="

    ORCH_STATUS="$(vm_status "${ORCH_VM}")"
    if [[ "${ORCH_STATUS}" == "NOT_FOUND" ]]; then
        log "Creating orchestrator VM ${ORCH_VM} (e2-micro)..."
        gcloud compute instances create "${ORCH_VM}" \
            --zone="${ZONE}" \
            --machine-type="e2-micro" \
            --image-family="ubuntu-2204-lts" \
            --image-project="ubuntu-os-cloud" \
            --boot-disk-size="10GB" \
            --boot-disk-type=pd-standard \
            --scopes=cloud-platform \
            --metadata=enable-oslogin=FALSE
        log "  Orchestrator VM created."
    elif [[ "${ORCH_STATUS}" == "TERMINATED" ]]; then
        log "Starting orchestrator VM ${ORCH_VM}..."
        gcloud compute instances start "${ORCH_VM}" --zone="${ZONE}"
    else
        log "Orchestrator VM ${ORCH_VM} already running (status: ${ORCH_STATUS})."
    fi

    wait_for_ssh "${ORCH_VM}"

    log "Cloning / updating repo on orchestrator VM (branch: ${BRANCH})..."
    ssh_cmd "${ORCH_VM}" -- bash -s <<ORCH_SETUP
set -euo pipefail
sudo apt-get update -qq && sudo apt-get install -y -qq tmux git >/dev/null 2>&1 || true
if [[ ! -f "\$HOME/.ssh/google_compute_engine" ]]; then
    ssh-keygen -t rsa -f "\$HOME/.ssh/google_compute_engine" -N "" -q
fi
REPO_DIR="\$HOME/splice-neoepitope-pipeline"
if [[ -d "\${REPO_DIR}/.git" ]]; then
    git -C "\${REPO_DIR}" fetch --all --prune
    git -C "\${REPO_DIR}" checkout "${BRANCH}"
    git -C "\${REPO_DIR}" pull origin "${BRANCH}"
else
    git clone --branch "${BRANCH}" "${REPO_URL}" "\${REPO_DIR}"
fi
ORCH_SETUP

    log "Starting pipeline in detached tmux session on ${ORCH_VM}..."
    ssh_cmd "${ORCH_VM}" -- bash -s <<ORCH_RUN
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"
tmux kill-session -t orchestrator 2>/dev/null || true
tmux new-session -d -s orchestrator "
    bash scripts/run_cloud_gpu.sh \\
        --mode ${MODE} \\
        --samples ${SAMPLES} \\
        --branch ${BRANCH} \\
        --zone ${ZONE} \\
        --_cloud-internal \\
        $([[ "${KEEP_VM}" == true ]] && echo "--keep-vm") \\
        2>&1 | tee orchestrator.log
    gcloud storage cp orchestrator.log "gs://${GCS_BUCKET}/logs/orchestrator.log" 2>/dev/null || true
    echo 'Orchestrator finished — shutting down.'
    sudo shutdown -h now
"
ORCH_RUN

    log "================================================================"
    log "Pipeline is running on ${ORCH_VM}. You can close your laptop."
    log ""
    log "Monitor progress:"
    log "  gcloud compute ssh ${ORCH_VM} --zone=${ZONE} --tunnel-through-iap -- tail -f splice-neoepitope-pipeline/orchestrator.log"
    log ""
    log "Attach to tmux session:"
    log "  gcloud compute ssh ${ORCH_VM} --zone=${ZONE} --tunnel-through-iap -- tmux attach -t orchestrator"
    log ""
    log "Results (after completion):"
    log "  gcloud storage cp -r ${GCS_PATH}/<patient_id>/reports ./report"
    log "  open report/report.html"
    log ""
    log "Delete orchestrator when done:"
    log "  gcloud compute instances delete ${ORCH_VM} --zone=${ZONE} --quiet"
    log "================================================================"
    exit 0
fi

# Stop VM on exit (success or failure) to avoid charges.
# Pass --keep-vm to skip this for debugging sessions.
trap '
    if [[ "${KEEP_VM}" == true ]]; then
        log "Cleanup: --keep-vm set — VM left running. Stop manually when done:"
        log "  gcloud compute instances stop ${PIPELINE_VM} --zone=${ZONE}"
    else
        log "Cleanup: stopping ${PIPELINE_VM}..."
        gcloud compute instances stop "${PIPELINE_VM}" --zone="${ZONE}" 2>/dev/null || true
    fi
' EXIT

# ---------------------------------------------------------------------------
# Ensure GCS bucket exists
# ---------------------------------------------------------------------------
if gcloud storage ls "gs://${GCS_BUCKET}" &>/dev/null; then
    log "GCS bucket gs://${GCS_BUCKET} already exists."
else
    log "Creating GCS bucket gs://${GCS_BUCKET}..."
    gcloud storage buckets create "gs://${GCS_BUCKET}" --location="${ZONE%-*}"
    log "  Bucket created."
fi

PROJECT_NUMBER="$(curl -sH 'Metadata-Flavor: Google' \
    http://metadata.google.internal/computeMetadata/v1/project/numeric-project-id 2>/dev/null \
    || gcloud projects describe "$(gcloud config get-value project 2>/dev/null)" \
       --format='value(projectNumber)' 2>/dev/null \
    || true)"
if [[ -n "${PROJECT_NUMBER}" ]]; then
    COMPUTE_SA="${PROJECT_NUMBER}-compute@developer.gserviceaccount.com"
    log "Granting ${COMPUTE_SA} objectAdmin on gs://${GCS_BUCKET}..."
    gcloud storage buckets add-iam-policy-binding "gs://${GCS_BUCKET}" \
        --member="serviceAccount:${COMPUTE_SA}" --role="roles/storage.objectAdmin" \
        2>/dev/null || log "  Warning: could not set bucket IAM (may already be set)."
else
    log "  Warning: could not determine project number — skipping bucket IAM grant."
fi

# ---------------------------------------------------------------------------
# Create or start pipeline VM
# ---------------------------------------------------------------------------
VM_STATUS="$(vm_status "${PIPELINE_VM}")"
if [[ "${VM_STATUS}" == "TERMINATED" ]]; then
    log "Starting ${PIPELINE_VM}..."
    gcloud compute instances start "${PIPELINE_VM}" --zone="${ZONE}"
elif [[ "${VM_STATUS}" == "NOT_FOUND" ]]; then
    log "Pipeline VM ${PIPELINE_VM} not found — creating it..."
    gcloud compute instances create "${PIPELINE_VM}" \
        --zone="${ZONE}" \
        --machine-type="${MACHINE_TYPE}" \
        --accelerator="${ACCELERATOR}" \
        --maintenance-policy=TERMINATE \
        --image-family="${IMAGE_FAMILY}" \
        --image-project="${IMAGE_PROJECT}" \
        --boot-disk-size="${DISK_SIZE}" \
        --boot-disk-type=pd-ssd \
        --scopes=cloud-platform \
        --metadata=enable-oslogin=FALSE
    log "  VM created."
else
    log "${PIPELINE_VM} already running (status: ${VM_STATUS})."
fi

wait_for_ssh "${PIPELINE_VM}"

# ---------------------------------------------------------------------------
# Verify the image's pre-installed NVIDIA driver works. The image family
# `common-cu129-ubuntu-2204-nvidia-580` ships 535.288.01 (matched
# userspace + kernel module) — Pascal-compatible (535 < 575). Do NOT
# attempt to install nvidia-headless-570-server: as of 2026-05-27,
# Ubuntu has reduced 570-server to a hard-Depends wrapper for 580-server,
# which drops SM 6.0 support and breaks P100. Do NOT run `apt-get update
# && apt-get install nvidia-utils-535-server` either — the archive churns
# the userspace 535 past the image's kernel-module version (535.288 →
# 535.309 observed), and the source package doesn't ship a usable
# dkms.conf so DKMS can't rebuild to align. Leave the image alone.
# See CLAUDE.md "Infrastructure" + Issue #522.
# ---------------------------------------------------------------------------
log "Verifying NVIDIA driver on ${PIPELINE_VM}..."
if ! ssh_cmd "${PIPELINE_VM}" -- nvidia-smi --query-gpu=name,driver_version --format=csv,noheader &>/dev/null; then
    log "ERROR: nvidia-smi failed on ${PIPELINE_VM} — image's driver not loading."
    log "  Likely a new image-family regression. Investigate before retrying:"
    log "    gcloud compute ssh ${PIPELINE_VM} --zone=${ZONE} --tunnel-through-iap -- nvidia-smi"
    exit 1
fi
DRIVER_INFO=$(ssh_cmd "${PIPELINE_VM}" -- nvidia-smi --query-gpu=name,driver_version --format=csv,noheader 2>&1 | head -1)
log "  Driver OK: ${DRIVER_INFO}"

# ---------------------------------------------------------------------------
# Setup VM (idempotent — skips steps already done)
# ---------------------------------------------------------------------------
log "Running setup_vm.sh on ${PIPELINE_VM}..."
ssh_cmd "${PIPELINE_VM}" -- bash -s -- --repo-branch "${BRANCH}" --no-next-steps < scripts/setup_vm.sh

if [[ -n "${PREPARE_DATA_SCRIPT}" ]]; then
    log "Preparing test data on ${PIPELINE_VM}..."
    ssh_cmd "${PIPELINE_VM}" -- bash -s <<REMOTE
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"
bash ${PREPARE_DATA_SCRIPT} --no-next-steps
REMOTE
fi

# ---------------------------------------------------------------------------
# Pull branch and start pipeline
# ---------------------------------------------------------------------------
log "Pulling branch '${BRANCH}' and starting pipeline on ${PIPELINE_VM}..."
ssh_cmd "${PIPELINE_VM}" -- bash -s <<EOF
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"

git fetch --all --prune
git checkout "${BRANCH}"
git pull origin "${BRANCH}"

tmux kill-session -t pipeline 2>/dev/null || true

tmux new-session -d -s pipeline "
    source \"\$HOME/miniforge3/etc/profile.d/conda.sh\" 2>/dev/null || true
    conda activate snakemake
    snakemake --unlock \\
        --configfile ${CONFIG_FILE} ${GPU_CONFIG_FILE} \\
        --config samples_tsv=${SAMPLES} 2>/dev/null || true
    snakemake \\
        --cores \$(nproc) \\
        --use-conda \\
        --rerun-triggers mtime \\
        --rerun-incomplete \\
        --configfile ${CONFIG_FILE} ${GPU_CONFIG_FILE} \\
        --config samples_tsv=${SAMPLES} \\
        2>&1 | tee pipeline.log
    echo 'Pipeline finished.'
"

echo "Pipeline started in tmux session 'pipeline'."
EOF

# ---------------------------------------------------------------------------
# Poll for completion
# ---------------------------------------------------------------------------
log "Polling pipeline.log for completion..."
log "  Monitor with:"
log "    gcloud compute ssh ${PIPELINE_VM} --zone=${ZONE} --tunnel-through-iap -- tail -f splice-neoepitope-pipeline/pipeline.log"
while true; do
    DONE="$(ssh_cmd "${PIPELINE_VM}" --command \
        "grep -cE 'steps \(100%\) done|Nothing to be done' \$HOME/splice-neoepitope-pipeline/pipeline.log 2>/dev/null || echo 0" \
        2>/dev/null | tail -1)"
    if [[ "${DONE}" -ge 1 ]]; then
        log "  Pipeline finished (100% done detected in log)."
        break
    fi
    FAILED="$(ssh_cmd "${PIPELINE_VM}" --command \
        "grep -c 'Error\|Exiting because' \$HOME/splice-neoepitope-pipeline/pipeline.log 2>/dev/null || echo 0" \
        2>/dev/null | tail -1)"
    if [[ "${FAILED}" -ge 1 ]]; then
        echo "ERROR: Pipeline appears to have failed. Check pipeline.log on ${PIPELINE_VM}." >&2
        exit 1
    fi
    PROGRESS="$(ssh_cmd "${PIPELINE_VM}" --command \
        "tail -3 \$HOME/splice-neoepitope-pipeline/pipeline.log 2>/dev/null" \
        2>/dev/null | grep -v '^$' | sed 's/^/    /')"
    log "  Still running — checking again in 60s..."
    [[ -n "${PROGRESS}" ]] && echo "${PROGRESS}"
    sleep 60
done

# ---------------------------------------------------------------------------
# Upload results to GCS
# ---------------------------------------------------------------------------
log ""
log "=== Uploading results to GCS ==="

ssh_cmd "${PIPELINE_VM}" -- bash -s <<EOF
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"
[[ -d "${RESULTS_DIR}" ]] || { echo "ERROR: ${RESULTS_DIR}/ not found." >&2; exit 1; }
gcloud storage cp -r "${RESULTS_DIR}" "${GCS_PATH%/*}/"
echo "Results upload complete."
[[ -d "logs" ]] && gcloud storage cp -r "logs" "gs://${GCS_BUCKET}/" && echo "Logs upload complete."
[[ -f "pipeline.log" ]] && gcloud storage cp "pipeline.log" "gs://${GCS_BUCKET}/logs/pipeline.log" && echo "Pipeline log uploaded."
EOF

# ---------------------------------------------------------------------------
# Clean up per-patient data from VM (runs only after confirmed GCS upload)
# ---------------------------------------------------------------------------
log ""
log "=== Cleaning up per-patient data from VM ==="

ssh_cmd "${PIPELINE_VM}" -- bash -s <<EOF
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"
if [[ -d "data/${PATIENT_ID}" ]]; then
    rm -rf "data/${PATIENT_ID}"
    echo "Deleted data/${PATIENT_ID}/"
fi
if [[ -d "results/${PATIENT_ID}/alignment" ]]; then
    find "results/${PATIENT_ID}/alignment" -type f \
        \( -name "*.bam" -o -name "*.bam.bai" -o -name "*.bed" \) -delete
    echo "Deleted alignment intermediates for ${PATIENT_ID}"
fi
echo "VM cleanup complete."
EOF

# ===========================================================================
# Done
# ===========================================================================
log "================================================================"
log "Pipeline complete. Results at gs://${GCS_BUCKET}/results/"
log ""
log "Download the report:"
log "  gcloud storage cp -r ${GCS_PATH}/<patient_id>/reports ./report"
log "  open report/report.html"
log ""
log "Delete the VM when no longer needed:"
log "  gcloud compute instances delete ${PIPELINE_VM} --zone=${ZONE} --quiet"
log "================================================================"
