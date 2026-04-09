#!/usr/bin/env bash
# =============================================================================
# run_cloud_gpu.sh — Spin up a GCP Spot GPU VM, run full pipeline, auto-stop
# =============================================================================
#
# What this does:
#   1. Create a Spot n1-standard-4 + NVIDIA T4 VM in europe-west1-b
#   2. Wait for SSH to become available
#   3. Clone / update the repo and check out the current branch
#   4. Run setup_cloud.sh  — installs conda, Snakemake, downloads GRCh38 ref
#   5. Run setup_tcrdock_vm.sh — installs CUDA, TCRdock, AlphaFold v2 params
#   6. Run the pipeline with --configfile config/tcrdock_gpu.yaml overlay
#   7. VM shuts itself down when the pipeline finishes (saves cost)
#
# Requirements (on your local machine):
#   - gcloud CLI authenticated: gcloud auth login && gcloud auth application-default login
#   - Active GCP project set: gcloud config set project <PROJECT_ID>
#   - GPU quota for T4 in europe-west1 (request at console.cloud.google.com if needed)
#
# Usage:
#   bash scripts/run_tcrdock_cloud.sh [--branch <branch>] [--zone <zone>]
#
# Options:
#   --branch <branch>   Git branch to check out on the VM (default: current branch)
#   --zone <zone>       GCP zone (default: europe-west1-b; us-central1 quota exhausted)
#
# VM spec:
#   n1-standard-4 (4 vCPU, 15 GB RAM) + NVIDIA T4 Spot (~$0.14/hr)
#   50 GB boot disk (enough for AlphaFold params + reference + results)
#
# After the run:
#   - Results are in ~/splice-neoepitope-pipeline/results/ on the VM *disk*
#   - The VM is stopped (not deleted) so you can inspect logs and results with:
#       gcloud compute scp splice-tcrdock-spot:~/splice-neoepitope-pipeline/pipeline.log .
#       gcloud compute scp -r splice-tcrdock-spot:~/splice-neoepitope-pipeline/results ./cloud_results
#   - Delete the VM when done to stop disk charges:
#       gcloud compute instances delete splice-tcrdock-spot --zone=<zone> --quiet
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
VM_NAME="splice-tcrdock-spot"
ZONE="europe-west1-b"
MACHINE_TYPE="n1-standard-4"
ACCELERATOR="type=nvidia-tesla-t4,count=1"
DISK_SIZE="50GB"
IMAGE_FAMILY="ubuntu-2204-lts"
IMAGE_PROJECT="ubuntu-os-cloud"
REPO_URL="https://github.com/Jin-HoMLee/splice-neoepitope-pipeline.git"
REPO_DIR='$HOME/splice-neoepitope-pipeline'  # evaluated on the VM

# Default branch = current local branch
BRANCH="$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "main")"

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --branch)
            BRANCH="${2:?--branch requires a value}"; shift 2 ;;
        --zone)
            ZONE="${2:?--zone requires a value}"; shift 2 ;;
        *)
            echo "Unknown argument: $1" >&2
            echo "Usage: $0 [--branch <branch>] [--zone <zone>]" >&2
            exit 1 ;;
    esac
done

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

log "============================================================"
log "TCRdock cloud run"
log "  VM:     ${VM_NAME}"
log "  Zone:   ${ZONE}"
log "  Branch: ${BRANCH}"
log "============================================================"

# ---------------------------------------------------------------------------
# Step 1 — Create the Spot GPU VM (skip if it already exists)
# ---------------------------------------------------------------------------
log "Step 1: Creating VM ${VM_NAME} ..."

if gcloud compute instances describe "${VM_NAME}" --zone="${ZONE}" &>/dev/null; then
    VM_STATUS="$(gcloud compute instances describe "${VM_NAME}" --zone="${ZONE}" \
        --format='value(status)')"
    log "  VM already exists (status: ${VM_STATUS})"
    if [[ "${VM_STATUS}" == "TERMINATED" ]]; then
        log "  VM is stopped — starting it..."
        gcloud compute instances start "${VM_NAME}" --zone="${ZONE}"
    fi
else
    gcloud compute instances create "${VM_NAME}" \
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
    log "  VM created."
fi

# ---------------------------------------------------------------------------
# Step 2 — Wait for SSH
# ---------------------------------------------------------------------------
log "Step 2: Waiting for SSH to become available (up to 3 min)..."

MAX_WAIT=180
ELAPSED=0
until gcloud compute ssh "${VM_NAME}" --zone="${ZONE}" \
        --command="echo ssh-ok" --quiet 2>/dev/null; do
    if [[ ${ELAPSED} -ge ${MAX_WAIT} ]]; then
        echo "ERROR: SSH did not become available after ${MAX_WAIT}s." >&2
        exit 1
    fi
    sleep 10
    ELAPSED=$((ELAPSED + 10))
    log "  Still waiting... (${ELAPSED}s)"
done
log "  SSH ready."

# ---------------------------------------------------------------------------
# Step 3 — Clone / update repo and check out branch
# ---------------------------------------------------------------------------
log "Step 3: Cloning / updating repo (branch: ${BRANCH})..."

# We pass BRANCH and REPO_URL as env vars into the remote shell
gcloud compute ssh "${VM_NAME}" --zone="${ZONE}" -- bash -s <<EOF
set -euo pipefail

REPO_DIR="\$HOME/splice-neoepitope-pipeline"
BRANCH="${BRANCH}"
REPO_URL="${REPO_URL}"

if [[ -d "\${REPO_DIR}/.git" ]]; then
    echo "  Repo already cloned — fetching and checking out \${BRANCH}..."
    git -C "\${REPO_DIR}" fetch --all --prune
    git -C "\${REPO_DIR}" checkout "\${BRANCH}"
    git -C "\${REPO_DIR}" pull origin "\${BRANCH}"
else
    echo "  Cloning \${REPO_URL} (branch: \${BRANCH})..."
    git clone --branch "\${BRANCH}" "\${REPO_URL}" "\${REPO_DIR}"
fi

echo "  Repo ready at \${REPO_DIR} (branch: \$(git -C \${REPO_DIR} rev-parse --abbrev-ref HEAD))"
EOF

# ---------------------------------------------------------------------------
# Step 4 — Run setup_cloud.sh (conda + Snakemake + reference data)
# ---------------------------------------------------------------------------
log "Step 4: Running setup_cloud.sh (conda + Snakemake + GRCh38 reference)..."
log "  This downloads ~900 MB of reference data — may take 20-40 min."

gcloud compute ssh "${VM_NAME}" --zone="${ZONE}" -- bash -s <<'EOF'
set -euo pipefail
cd "$HOME/splice-neoepitope-pipeline"
bash scripts/setup_cloud.sh
EOF

# ---------------------------------------------------------------------------
# Step 5 — Run setup_tcrdock_vm.sh (CUDA + TCRdock + AlphaFold params)
# ---------------------------------------------------------------------------
log "Step 5: Running setup_tcrdock_vm.sh (CUDA, TCRdock, AlphaFold params)..."
log "  AlphaFold params download ~3.5 GB — may take 10-20 min."

gcloud compute ssh "${VM_NAME}" --zone="${ZONE}" -- bash -s <<'EOF'
set -euo pipefail
cd "$HOME/splice-neoepitope-pipeline"
bash scripts/setup_tcrdock_vm.sh
EOF

# ---------------------------------------------------------------------------
# Step 6 — Run pipeline with TCRdock enabled; VM shuts down on completion
# ---------------------------------------------------------------------------
log "Step 6: Running pipeline with TCRdock enabled..."
log "  The VM will shut down automatically when the pipeline finishes."
log "  To monitor progress, open a second terminal and run:"
log "    gcloud compute ssh ${VM_NAME} --zone=${ZONE} -- tail -f splice-neoepitope-pipeline/pipeline.log"

gcloud compute ssh "${VM_NAME}" --zone="${ZONE}" -- bash -s <<'EOF'
set -euo pipefail
cd "$HOME/splice-neoepitope-pipeline"

# Activate the snakemake environment (created by setup_cloud.sh)
eval "$(conda shell.bash hook 2>/dev/null || "$HOME/miniforge3/bin/conda" shell.bash hook)"
conda activate snakemake

# Run pipeline with both config files (TCRdock overlay enables GPU step)
snakemake \
    --cores "$(nproc)" \
    --use-conda \
    --rerun-triggers mtime \
    --configfile config/config.yaml \
    --configfile config/tcrdock_gpu.yaml \
    2>&1 | tee pipeline.log

# Shut down the VM (saves cost — disk is preserved for result retrieval)
echo "Pipeline finished — shutting down VM..."
sudo shutdown -h now
EOF

# ---------------------------------------------------------------------------
# Done — print retrieval instructions
# ---------------------------------------------------------------------------
log "============================================================"
log "Pipeline submitted. VM will stop when finished."
log ""
log "To copy results back to your local machine:"
log "  gcloud compute scp ${VM_NAME}:~/splice-neoepitope-pipeline/pipeline.log . --zone=${ZONE}"
log "  gcloud compute scp -r ${VM_NAME}:~/splice-neoepitope-pipeline/results ./cloud_results --zone=${ZONE}"
log ""
log "To delete the VM and stop disk charges:"
log "  gcloud compute instances delete ${VM_NAME} --zone=${ZONE} --quiet"
log "============================================================"
