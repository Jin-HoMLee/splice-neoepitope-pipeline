#!/usr/bin/env bash
# =============================================================================
# run_cloud_gpu.sh — Full pipeline lifecycle: CPU steps → TCRdock GPU → auto-stop
# =============================================================================
#
# Three phases:
#
#   Phase 1 — CPU VM (splice-prod-test)
#     Start VM, pull branch, run pipeline steps 1-5 (alignment → MHCflurry)
#     in a tmux session, poll pipeline.log until 100% done (VM stays running).
#
#   Phase 2 — Copy results
#     Copy the results directories needed by TCRdock (junctions, contigs,
#     predictions — not raw BAMs) while VM is still up, then stop CPU VM.
#
#   Phase 3 — GPU Spot VM (splice-tcrdock-spot)
#     Create VM, install conda + Snakemake + CUDA + TCRdock + AlphaFold params,
#     copy results in, run only the TCRdock + report rules (earlier steps are
#     already satisfied by the copied results), then VM auto-stops.
#
# Requirements:
#   - gcloud CLI authenticated: gcloud auth login
#   - Active project: gcloud config set project <PROJECT_ID>
#   - T4 GPU quota in europe-west1 (PREEMPTIBLE_NVIDIA_T4_GPUS >= 1)
#
# Usage:
#   bash scripts/run_cloud_gpu.sh [--branch <branch>] [--zone <zone>]
#
# Options:
#   --branch <branch>   Git branch to use (default: current local branch)
#   --zone   <zone>     GCP zone (default: europe-west1-b)
#
# After the run, copy results from the GPU VM before deleting it:
#   gcloud compute scp --tunnel-through-iap --recurse \
#       splice-tcrdock-spot:~/splice-neoepitope-pipeline/results/local/reports \
#       ./tcrdock_report --zone=europe-west1-b
#   gcloud compute instances delete splice-tcrdock-spot --zone=europe-west1-b --quiet
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
CPU_VM="splice-prod-test"
GPU_VM="splice-tcrdock-spot"
ZONE="europe-west1-b"
MACHINE_TYPE="n1-standard-4"
ACCELERATOR="type=nvidia-tesla-t4,count=1"
DISK_SIZE="50GB"
IMAGE_FAMILY="ubuntu-2204-lts"
IMAGE_PROJECT="ubuntu-os-cloud"
REPO_URL="https://github.com/Jin-HoMLee/splice-neoepitope-pipeline.git"

# GCS bucket for VM-to-VM result transfer (created automatically if absent)
GCS_BUCKET="tcrdock-handoff"
GCS_PATH="gs://${GCS_BUCKET}/results"

# Results subdirectories to upload (excludes raw_data/ which contains large BAMs)
RESULT_DIRS=("junctions" "contigs" "peptides" "predictions" "reports")

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
    # Wrapper: always use IAP tunnel so the script works from any network
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
log "Cloud GPU pipeline run"
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

# Grant the default compute service account write access so VMs can upload/download
PROJECT_NUMBER="$(gcloud projects describe "$(gcloud config get-value project)" \
    --format='value(projectNumber)' 2>/dev/null)"
COMPUTE_SA="${PROJECT_NUMBER}-compute@developer.gserviceaccount.com"
log "Granting ${COMPUTE_SA} objectAdmin on gs://${GCS_BUCKET}..."
gsutil iam ch "serviceAccount:${COMPUTE_SA}:objectAdmin" "gs://${GCS_BUCKET}" 2>/dev/null || \
    log "  Warning: could not set bucket IAM (may already be set)."

# ===========================================================================
# Phase 1 — CPU pipeline (steps 1-5)
# ===========================================================================
log ""
log "=== Phase 1: CPU pipeline (steps 1-5) ==="

# Start CPU VM if it is stopped
CPU_STATUS="$(vm_status "${CPU_VM}")"
if [[ "${CPU_STATUS}" == "TERMINATED" ]]; then
    log "Starting ${CPU_VM}..."
    gcloud compute instances start "${CPU_VM}" --zone="${ZONE}"
elif [[ "${CPU_STATUS}" == "NOT_FOUND" ]]; then
    echo "ERROR: CPU VM '${CPU_VM}' not found in zone ${ZONE}." >&2
    echo "       Create it with scripts/setup_cloud.sh first." >&2
    exit 1
fi

wait_for_ssh "${CPU_VM}"

# Pull latest branch and start pipeline in a tmux session
log "Pulling branch '${BRANCH}' and starting pipeline on ${CPU_VM}..."
ssh_cmd "${CPU_VM}" -- bash -s <<EOF
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"

git fetch --all --prune
git checkout "${BRANCH}"
git pull origin "${BRANCH}"

# Kill any stale pipeline tmux session
tmux kill-session -t pipeline 2>/dev/null || true

# Start fresh tmux session: run pipeline (no shutdown — orchestrator handles that)
tmux new-session -d -s pipeline "
    source \"\$HOME/miniforge3/etc/profile.d/conda.sh\" 2>/dev/null || true
    conda activate snakemake
    snakemake \\
        --cores \$(nproc) \\
        --use-conda \\
        --rerun-triggers mtime \\
        --configfile config/config.yaml \\
        2>&1 | tee pipeline.log
    echo 'Pipeline finished.'
"

echo "Pipeline started in tmux session 'pipeline'."
EOF

# Poll pipeline.log for completion — VM stays running so we can copy directly
log "Polling pipeline.log for completion (this may take 1-3 hours)..."
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
    # Also check if snakemake failed (VM still running but pipeline errored)
    FAILED="$(ssh_cmd "${CPU_VM}" --command \
        "grep -c 'Error\|Exiting because' \$HOME/splice-neoepitope-pipeline/pipeline.log 2>/dev/null || echo 0" \
        2>/dev/null | tail -1)"
    if [[ "${FAILED}" -ge 1 ]]; then
        echo "ERROR: Pipeline appears to have failed. Check pipeline.log on ${CPU_VM}." >&2
        exit 1
    fi
    # Show last 3 lines of pipeline.log so the user can see what's happening
    PROGRESS="$(ssh_cmd "${CPU_VM}" --command \
        "tail -3 \$HOME/splice-neoepitope-pipeline/pipeline.log 2>/dev/null" \
        2>/dev/null | grep -v '^$' | sed 's/^/    /')"
    log "  Still running — checking again in 60s..."
    [[ -n "${PROGRESS}" ]] && echo "${PROGRESS}"
    sleep 60
done

# ===========================================================================
# Phase 2 — Copy results from CPU VM, then stop it
# ===========================================================================
log ""
log "=== Phase 2: Copying results from ${CPU_VM} ==="

# VM is still running from Phase 1 — upload results to GCS directly, then stop
log "Uploading results from ${CPU_VM} to gs://${GCS_BUCKET}..."
ssh_cmd "${CPU_VM}" -- bash -s <<EOF
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"
for dir in ${RESULT_DIRS[*]}; do
    if [[ -d "results/\${dir}" ]]; then
        echo "  Uploading results/\${dir}/ ..."
        gcloud storage cp -r "results/\${dir}" "${GCS_PATH}/\${dir}"
    else
        echo "  results/\${dir}/ not found — skipping."
    fi
done
echo "Upload complete."
EOF

log "Stopping ${CPU_VM} to save cost..."
gcloud compute instances stop "${CPU_VM}" --zone="${ZONE}"
log "  ${CPU_VM} stopped."

# ===========================================================================
# Phase 3 — GPU Spot VM: TCRdock structural validation
# ===========================================================================
log ""
log "=== Phase 3: GPU Spot VM — TCRdock ==="

# Create GPU VM if it does not exist
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

# Clone / update repo on GPU VM
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

# Install conda + Snakemake (no reference download needed — just the runner)
log "Installing conda and Snakemake on GPU VM..."
ssh_cmd "${GPU_VM}" -- bash -s <<'EOF'
set -euo pipefail

# Install conda if absent
if ! command -v conda &>/dev/null && [[ ! -f "$HOME/miniforge3/bin/conda" ]]; then
    echo "  Installing Miniforge3..."
    curl -fsSL https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
        -o /tmp/miniforge.sh
    bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"
    rm /tmp/miniforge.sh
    "$HOME/miniforge3/bin/conda" init bash
fi

source "$HOME/miniforge3/etc/profile.d/conda.sh"

# Create snakemake env if absent
if ! conda env list | grep -q "^snakemake "; then
    echo "  Creating snakemake conda environment..."
    conda create -n snakemake -c conda-forge -c bioconda \
        "snakemake>=8.0,<9" python=3.11 -y
fi
echo "  Snakemake environment ready."
EOF

# Run setup_tcrdock_vm.sh (CUDA + TCRdock + AlphaFold params)
log "Running setup_tcrdock_vm.sh on GPU VM (CUDA + TCRdock + AlphaFold params)..."
log "  This downloads ~3.5 GB of AlphaFold params — may take 15-20 min."
ssh_cmd "${GPU_VM}" -- bash -s <<'EOF'
set -euo pipefail
cd "$HOME/splice-neoepitope-pipeline"
bash scripts/setup_tcrdock_vm.sh
EOF

# Download results from GCS directly onto GPU VM
log "Downloading results from gs://${GCS_BUCKET} onto GPU VM..."
ssh_cmd "${GPU_VM}" -- bash -s <<EOF
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"
mkdir -p results
gcloud storage cp -r "${GCS_PATH}/*" results/
echo "Download complete."
EOF

# Run Snakemake with TCRdock overlay — only TCRdock + report rules will run
# (all earlier outputs are already present from the CPU pipeline)
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
    --configfile config/config.yaml \
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
log "        ${GPU_VM}:~/splice-neoepitope-pipeline/results/local/reports \\"
log "        ./tcrdock_report --zone=${ZONE}"
log ""
log "To delete the GPU VM and stop disk charges:"
log "    gcloud compute instances delete ${GPU_VM} --zone=${ZONE} --quiet"
log "================================================================"
