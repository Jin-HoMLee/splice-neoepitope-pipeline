#!/usr/bin/env bash
# =============================================================================
# run_cloud_gpu.sh — Pipeline lifecycle: CPU steps → TCRdock GPU → auto-stop
# =============================================================================
#
# Three phases:
#
#   Phase 1 — CPU VM (neoepitope-predict-cpu)
#     Sync prior state from GCS, pull branch, run pipeline steps 1-5
#     (alignment → MHCflurry) in a tmux session, poll until 100% done.
#
#   Phase 2 — GCS handoff
#     Upload results, logs, and Snakemake metadata to GCS; stop CPU VM.
#
#   Phase 3 — GPU Spot VM (pipeline-spot-gpu)
#     Create VM, install conda + Snakemake + CUDA + TCRdock + AlphaFold params,
#     sync results/logs/metadata from GCS, run TCRdock + report rules, upload
#     results/logs/metadata back to GCS, then VM auto-stops.
#
# Modes:
#   --mode prod   Full genome, production config (default)
#   --mode test   chr22, 500K reads — fast validation (~1 hour)
#
# Detached mode (--detach):
#   Launches a tiny e2-micro "orchestrator" VM that runs this script in a
#   tmux session. You can close your laptop — the orchestrator manages all
#   three phases using internal networking, then auto-stops when done.
#
# Requirements:
#   - gcloud CLI authenticated: gcloud auth login
#   - Active project: gcloud config set project <PROJECT_ID>
#   - P100 GPU quota in the target zone (PREEMPTIBLE_NVIDIA_P100_GPUS >= 1)
#
# Usage:
#   bash scripts/run_cloud_gpu.sh --samples config/samples/patient_002.tsv [--mode test|prod] [--branch <branch>] [--zone <zone>] [--detach]
#
# After the run, download the report:
#   gcloud storage cp -r gs://splice-neoepitope-project/results/{patient_id}/reports ./report   # test or prod
#   open report/report.html
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
CPU_VM="neoepitope-predict-cpu"
GPU_VM="pipeline-spot-gpu"
ZONE="europe-west1-b"
CPU_MACHINE_TYPE="n2-highmem-8"  # 64 GB RAM: razers3 (OptiType) peaks at ~36 GB on full RNA-seq FASTQs; n2 preferred over n1 for availability
GPU_MACHINE_TYPE="n1-standard-4"
ACCELERATOR="type=nvidia-tesla-p100,count=1"
DISK_SIZE="100GB"
IMAGE_FAMILY="common-cu128-ubuntu-2204-nvidia-570"  # CUDA 12.8 pre-installed
IMAGE_PROJECT="deeplearning-platform-release"
REPO_URL="https://github.com/Jin-HoMLee/splice-neoepitope-pipeline.git"
ORCH_VM="neoepitope-orchestrator"

GCS_BUCKET="splice-neoepitope-project"

BRANCH="$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "main")"
MODE="prod"
SAMPLES=""
DETACH=false
CLOUD_INTERNAL=false

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --mode)    MODE="${2:?--mode requires test or prod}";      shift 2 ;;
        --branch)  BRANCH="${2:?--branch requires a value}";     shift 2 ;;
        --zone)    ZONE="${2:?--zone requires a value}";         shift 2 ;;
        --samples) SAMPLES="${2:?--samples requires a path}";    shift 2 ;;
        --detach)  DETACH=true;                                   shift ;;
        --_cloud-internal) CLOUD_INTERNAL=true;                   shift ;;
        *)
            echo "Unknown argument: $1" >&2
            echo "Usage: $0 [--mode test|prod] [--branch <branch>] [--zone <zone>] [--detach]" >&2
            exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Mode-dependent settings
# ---------------------------------------------------------------------------
case "${MODE}" in
    test)
        CONFIG_FILE="config/test_config.yaml"
        RESULTS_PATH="results"
        LOGS_PATH="logs"
        GCS_RESULTS_PATH="gs://${GCS_BUCKET}/results"
        GCS_LOGS_PATH="gs://${GCS_BUCKET}/logs"
        PREPARE_DATA_SCRIPT="scripts/prepare_test_data.sh"
        [[ -z "${SAMPLES}" ]] && SAMPLES="config/samples/patient_001_test.tsv"
        ;;
    prod|production)
        CONFIG_FILE="config/config.yaml"
        RESULTS_PATH="results"
        LOGS_PATH="logs"
        GCS_RESULTS_PATH="gs://${GCS_BUCKET}/results"
        GCS_LOGS_PATH="gs://${GCS_BUCKET}/logs"
        PREPARE_DATA_SCRIPT=""
        [[ -z "${SAMPLES}" ]] && { echo "ERROR: --samples required for prod mode (e.g. --samples config/samples/patient_002.tsv)" >&2; exit 1; }
        ;;
    *)
        echo "ERROR: --mode must be 'test' or 'prod', got '${MODE}'" >&2
        exit 1 ;;
esac

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
        # --tunnel-through-iap: connects without a public IP via Google IAP.
        # If you see a "consider installing NumPy" warning, install it into
        # gcloud's own virtualenv: ~/.config/gcloud/virtenv/bin/pip install numpy
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
log "Cloud GPU pipeline run (mode: ${MODE})"
log "  CPU VM:  ${CPU_VM}"
log "  GPU VM:  ${GPU_VM}"
log "  Zone:    ${ZONE}"
log "  Branch:  ${BRANCH}"
log "  Config:  ${CONFIG_FILE}"
log "  Bucket:  gs://${GCS_BUCKET}"
if [[ "${DETACH}" == true ]]; then
    log "  Detach:  yes (orchestrator VM: ${ORCH_VM})"
fi
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
# Ensure tmux and git are available
sudo apt-get update -qq && sudo apt-get install -y -qq tmux git >/dev/null 2>&1 || true
# Pre-generate SSH key so gcloud compute ssh doesn't prompt interactively
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
        2>&1 | tee orchestrator.log
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
    log "When done, results are uploaded to GCS automatically:"
    log "  gcloud storage cp -r ${GCS_RESULTS_PATH}/reports ./tcrdock_report"
    log "  open tcrdock_report/[patient_id]/report.html"
    log ""
    log "The orchestrator VM auto-stops when finished. To delete it:"
    log "  gcloud compute instances delete ${ORCH_VM} --zone=${ZONE} --quiet"
    log "================================================================"
    exit 0
fi

# Stop both VMs on exit (success or failure) to avoid charges.
trap '
    log "Cleanup: ensuring VMs are stopped..."
    gcloud compute instances stop "${CPU_VM}" --zone="${ZONE}" 2>/dev/null || true
    gcloud compute instances stop "${GPU_VM}" --zone="${ZONE}" 2>/dev/null || true
' EXIT

# ---------------------------------------------------------------------------
# Ensure GCS handoff bucket exists
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
        2>/dev/null || \
        log "  Warning: could not set bucket IAM (may already be set)."
else
    log "  Warning: could not determine project number — skipping bucket IAM grant."
    log "  (If uploads fail later, run: gsutil iam ch serviceAccount:<SA>:objectAdmin gs://${GCS_BUCKET})"
fi

# ===========================================================================
# Phase 1 — CPU pipeline (steps 1-5)
# ===========================================================================
log ""
log "=== Phase 1: CPU pipeline (${CONFIG_FILE}) ==="

# ---------------------------------------------------------------------------
# 1a — VM provisioning
# ---------------------------------------------------------------------------
CPU_STATUS="$(vm_status "${CPU_VM}")"
if [[ "${CPU_STATUS}" == "TERMINATED" ]]; then
    log "Starting ${CPU_VM}..."
    gcloud compute instances start "${CPU_VM}" --zone="${ZONE}"
elif [[ "${CPU_STATUS}" == "NOT_FOUND" ]]; then
    log "CPU VM ${CPU_VM} not found — creating it..."
    gcloud compute instances create "${CPU_VM}" \
        --zone="${ZONE}" \
        --machine-type="${CPU_MACHINE_TYPE}" \
        --image-family="ubuntu-2204-lts" \
        --image-project="ubuntu-os-cloud" \
        --boot-disk-size="100GB" \
        --boot-disk-type=pd-ssd \
        --scopes=cloud-platform \
        --metadata=enable-oslogin=FALSE
fi

wait_for_ssh "${CPU_VM}"

# ---------------------------------------------------------------------------
# 1b — Environment setup (conda, Snakemake, reference data)
# ---------------------------------------------------------------------------
log "Running setup_cloud.sh on ${CPU_VM}..."
ssh_cmd "${CPU_VM}" -- bash -s -- --repo-branch "${BRANCH}" --no-next-steps < scripts/setup_cloud.sh

if [[ -n "${PREPARE_DATA_SCRIPT}" ]]; then
    log "Preparing data on ${CPU_VM} (${PREPARE_DATA_SCRIPT})..."
    ssh_cmd "${CPU_VM}" -- bash -s <<REMOTE
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"
bash ${PREPARE_DATA_SCRIPT} --no-next-steps
REMOTE
fi

# ---------------------------------------------------------------------------
# 1c — GCS sync: pull branch and download prior results, logs, and metadata
# ---------------------------------------------------------------------------
log "Pulling branch '${BRANCH}' and syncing prior GCS state on ${CPU_VM}..."
ssh_cmd "${CPU_VM}" -- bash -s <<EOF
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"

git fetch --all --prune
git checkout "${BRANCH}"
git pull origin "${BRANCH}"

# Download previous results, logs, and metadata from GCS so Snakemake can skip
# already-completed steps.  No-ops silently on a first run (bucket is empty).
mkdir -p "${RESULTS_PATH}"
gcloud storage rsync "${GCS_RESULTS_PATH}" "${RESULTS_PATH}" --recursive --preserve-posix 2>/dev/null \
    && echo "Results synced from GCS." \
    || echo "No prior results in GCS — starting fresh."
mkdir -p "${LOGS_PATH}"
gcloud storage rsync "${GCS_LOGS_PATH}" "${LOGS_PATH}" --recursive --preserve-posix 2>/dev/null \
    && echo "Logs synced from GCS." \
    || echo "No prior logs in GCS — starting fresh."
rm -rf .snakemake/metadata
mkdir -p .snakemake/metadata
gcloud storage rsync "gs://${GCS_BUCKET}/.snakemake/metadata" ".snakemake/metadata" --recursive --preserve-posix 2>/dev/null \
    && echo "Snakemake metadata downloaded." \
    || echo "No prior metadata in GCS — starting fresh."
EOF

# ---------------------------------------------------------------------------
# 1d — Pipeline execution
# ---------------------------------------------------------------------------
log "Starting pipeline in tmux session on ${CPU_VM}..."
ssh_cmd "${CPU_VM}" -- bash -s <<EOF
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"

tmux kill-session -t pipeline 2>/dev/null || true

tmux new-session -d -s pipeline "
    source \"\$HOME/miniforge3/etc/profile.d/conda.sh\" 2>/dev/null || true
    conda activate snakemake
    snakemake --unlock --configfile ${CONFIG_FILE} --config samples_tsv=${SAMPLES} 2>/dev/null || true
    snakemake \\
        --cores \$(nproc) \\
        --use-conda \\
        --rerun-triggers mtime \\
        --rerun-incomplete \\
        --configfile ${CONFIG_FILE} \\
        --config samples_tsv=${SAMPLES} \\
        2>&1 | tee pipeline.log
    echo 'Pipeline finished.'
"

echo "Pipeline started in tmux session 'pipeline'."
EOF

# ---------------------------------------------------------------------------
# 1e — Poll for completion
# ---------------------------------------------------------------------------
log "Polling pipeline.log for completion..."
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
# Phase 2 — Upload results, logs, and Snakemake metadata to GCS, stop CPU VM
# ===========================================================================
log ""
log "=== Phase 2: Uploading results, logs, and metadata from ${CPU_VM} ==="

log "Uploading results/, logs/, and .snakemake/metadata/ to gs://${GCS_BUCKET}/ (with --preserve-posix)..."
ssh_cmd "${CPU_VM}" -- bash -s <<EOF
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"
gcloud storage rsync "${RESULTS_PATH}" "${GCS_RESULTS_PATH}" --recursive --preserve-posix
echo "Results upload complete."
gcloud storage rsync "${LOGS_PATH}" "${GCS_LOGS_PATH}" --recursive --preserve-posix
echo "Logs upload complete."
gcloud storage rsync ".snakemake/metadata" "gs://${GCS_BUCKET}/.snakemake/metadata" --recursive
echo "Snakemake metadata upload complete."
EOF

log "Stopping ${CPU_VM} to save cost..."
gcloud compute instances stop "${CPU_VM}" --zone="${ZONE}" || \
    log "  Warning: could not stop ${CPU_VM} — stop it manually to avoid charges."

# ===========================================================================
# Phase 3 — GPU Spot VM: TCRdock structural validation
# ===========================================================================
log ""
log "=== Phase 3: GPU Spot VM — TCRdock ==="

# ---------------------------------------------------------------------------
# 3a — VM provisioning
# ---------------------------------------------------------------------------
GPU_STATUS="$(vm_status "${GPU_VM}")"
if [[ "${GPU_STATUS}" == "NOT_FOUND" ]]; then
    log "Creating GPU Spot VM ${GPU_VM}..."
    gcloud compute instances create "${GPU_VM}" \
        --zone="${ZONE}" \
        --machine-type="${GPU_MACHINE_TYPE}" \
        --accelerator="${ACCELERATOR}" \
        --maintenance-policy=TERMINATE \
        --provisioning-model=STANDARD \
        --instance-termination-action=STOP \
        --image-family="${IMAGE_FAMILY}" \
        --image-project="${IMAGE_PROJECT}" \
        --boot-disk-size="${DISK_SIZE}" \
        --boot-disk-type=pd-ssd \
        --scopes=cloud-platform \
        --metadata=enable-oslogin=FALSE
    log "  GPU VM created."
elif [[ "${GPU_STATUS}" == "TERMINATED" ]]; then
    log "GPU VM ${GPU_VM} exists (stopped) — patching scopes then starting..."
    gcloud compute instances set-service-account "${GPU_VM}" \
        --zone="${ZONE}" \
        --scopes=cloud-platform
    gcloud compute instances start "${GPU_VM}" --zone="${ZONE}"
else
    log "GPU VM ${GPU_VM} already running (status: ${GPU_STATUS})."
fi

wait_for_ssh "${GPU_VM}"

# ---------------------------------------------------------------------------
# 3b — Environment setup (repo, conda, Snakemake, TCRdock deps)
# ---------------------------------------------------------------------------
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

# tmux is used to run the TCRdock pipeline in a detachable session.
# Not present by default on the Deep Learning VM image.
if ! command -v tmux &>/dev/null; then
    sudo apt-get install -y -q tmux
fi

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
ssh_cmd "${GPU_VM}" -- bash -s <<EOF
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"
bash scripts/setup_tcrdock_vm.sh "${CONFIG_FILE}"
EOF

# ---------------------------------------------------------------------------
# 3c — GCS sync: download results, logs, and metadata; verify CPU-phase outputs
# ---------------------------------------------------------------------------
log "Syncing results/, logs/, and .snakemake/metadata/ from GCS onto GPU VM..."
ssh_cmd "${GPU_VM}" -- bash -s <<EOF
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"
# TCRdock Docker outputs are root-owned; wipe before downloading so rsync
# can overwrite them.  GCS is authoritative for the GPU VM — always start fresh.
sudo rm -rf "${RESULTS_PATH}"
mkdir -p "${RESULTS_PATH}"
gcloud storage rsync "${GCS_RESULTS_PATH}" "${RESULTS_PATH}" --recursive --preserve-posix 2>/dev/null \
    && echo "Results synced from GCS." \
    || echo "No prior results in GCS — starting fresh."
sudo rm -rf "${LOGS_PATH}"
mkdir -p "${LOGS_PATH}"
gcloud storage rsync "${GCS_LOGS_PATH}" "${LOGS_PATH}" --recursive --preserve-posix 2>/dev/null \
    && echo "Logs synced from GCS." \
    || echo "No prior logs in GCS — starting fresh."
# Wipe any stale metadata from previous GPU-VM runs before replacing with
# the CPU VM's authoritative version.  Leftover code/params hashes from an
# old run cause false "code changed" triggers and re-run of CPU-phase rules.
rm -rf .snakemake/metadata
mkdir -p .snakemake/metadata
gcloud storage rsync "gs://${GCS_BUCKET}/.snakemake/metadata" ".snakemake/metadata" --recursive --preserve-posix 2>/dev/null \
    && echo "Snakemake metadata downloaded." \
    || echo "No prior metadata in GCS — starting fresh."
EOF

log "Verifying CPU-phase outputs are present on GPU VM..."
ssh_cmd "${GPU_VM}" -- bash -s <<EOF
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"
MHCF="\$(find "${RESULTS_PATH}" -name "mhc_affinity.tsv" | head -1)"
if [[ -z "\${MHCF}" ]]; then
    echo "ERROR: No mhc_affinity.tsv found in ${RESULTS_PATH}/." >&2
    echo "       CPU phase may not have completed, or GCS sync failed." >&2
    exit 1
fi
echo "  Found: \${MHCF}"
EOF

# ---------------------------------------------------------------------------
# 3d — Pipeline execution (MHCflurry sentinel, TCRdock, poll for completion)
# ---------------------------------------------------------------------------
log "Starting TCRdock in tmux session on ${GPU_VM}..."
ssh_cmd "${GPU_VM}" -- bash -s <<EOF
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"

# Ensure MHCflurry models are present before the main TCRdock run.
# Target just the sentinel so it uses the correct python.yaml conda env
# (mhcflurry-downloads is not in the snakemake bootstrap env).
source "\$HOME/miniforge3/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate snakemake
snakemake \
    --cores 1 \
    --use-conda \
    --rerun-triggers mtime \
    --configfile ${CONFIG_FILE} config/tcrdock_gpu.yaml \
    --config samples_tsv=${SAMPLES} \
    -- \
    resources/mhcflurry_models.done
conda deactivate 2>/dev/null || true

# Remove the CPU-generated report so Snakemake regenerates it with the
# TCRdock structural view via generate_report_with_structure.
find ${RESULTS_PATH} -name report.html -delete

tmux kill-session -t tcrdock 2>/dev/null || true
tmux new-session -d -s tcrdock "
    source \"\$HOME/miniforge3/etc/profile.d/conda.sh\" 2>/dev/null || true
    conda activate snakemake
    snakemake --unlock --configfile ${CONFIG_FILE} config/tcrdock_gpu.yaml --config samples_tsv=${SAMPLES} 2>/dev/null || true
    snakemake \\
        --cores \$(nproc) \\
        --use-conda \\
        --rerun-triggers mtime \\
        --rerun-incomplete \\
        --configfile ${CONFIG_FILE} config/tcrdock_gpu.yaml \\
        --config samples_tsv=${SAMPLES} \\
        2>&1 | tee tcrdock.log
    echo 'TCRdock pipeline finished.'
"
echo "TCRdock started in tmux session 'tcrdock'."
EOF

log "Polling tcrdock.log for completion..."
log "  Monitor with:"
log "    gcloud compute ssh ${GPU_VM} --zone=${ZONE} --tunnel-through-iap -- tail -f splice-neoepitope-pipeline/tcrdock.log"
while true; do
    DONE="$(ssh_cmd "${GPU_VM}" --command \
        "grep -cE 'steps \(100%\) done|Nothing to be done' \$HOME/splice-neoepitope-pipeline/tcrdock.log 2>/dev/null || echo 0" \
        2>/dev/null | tail -1)"
    if [[ "${DONE}" -ge 1 ]]; then
        log "  TCRdock finished (100% done detected in log)."
        break
    fi
    FAILED="$(ssh_cmd "${GPU_VM}" --command \
        "grep -c 'Error\|Exiting because' \$HOME/splice-neoepitope-pipeline/tcrdock.log 2>/dev/null || echo 0" \
        2>/dev/null | tail -1)"
    if [[ "${FAILED}" -ge 1 ]]; then
        echo "ERROR: TCRdock appears to have failed. Check tcrdock.log on ${GPU_VM}." >&2
        exit 1
    fi
    PROGRESS="$(ssh_cmd "${GPU_VM}" --command \
        "tail -3 \$HOME/splice-neoepitope-pipeline/tcrdock.log 2>/dev/null" \
        2>/dev/null | grep -v '^$' | sed 's/^/    /')"
    log "  Still running — checking again in 60s..."
    [[ -n "${PROGRESS}" ]] && echo "${PROGRESS}"
    sleep 60
done

# ---------------------------------------------------------------------------
# 3e — GCS upload: results, logs, and metadata back to GCS
# ---------------------------------------------------------------------------
log "Uploading final results, logs, and metadata from ${GPU_VM} to GCS..."
ssh_cmd "${GPU_VM}" -- bash -s <<EOF
set -euo pipefail
cd "\$HOME/splice-neoepitope-pipeline"
gcloud storage rsync "${RESULTS_PATH}" "${GCS_RESULTS_PATH}" --recursive --preserve-posix
echo "Results upload complete."
gcloud storage rsync "${LOGS_PATH}" "${GCS_LOGS_PATH}" --recursive --preserve-posix
echo "Logs upload complete."
gcloud storage rsync ".snakemake/metadata" "gs://${GCS_BUCKET}/.snakemake/metadata" --recursive
echo "Snakemake metadata upload complete."
EOF

# ===========================================================================
# Done
# ===========================================================================
log "================================================================"
log "All phases complete. GPU VM will stop shortly."
log ""
log "To retrieve the report:"
log "  gcloud storage cp -r ${GCS_RESULTS_PATH}/{patient_id}/reports ./report"
log "  open report/report.html"
log ""
log "To delete the GPU VM and stop disk charges:"
log "    gcloud compute instances delete ${GPU_VM} --zone=${ZONE} --quiet"
log "================================================================"
