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
#       [--mode test|prod] [--branch <branch>] [--zone <zone>] [--detach] [--keep-vm] [--on-demand]
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
# Pin to specific image name (not family): family aliases on deeplearning-platform-release
# have all moved to nvidia-580/cu129 builds that ship -open driver variants, which require
# GSP firmware Pascal lacks (P100 PCI ID 10de:15f8 fails with "GPU ... is not supported by
# open nvidia.ko because it does not include the required GPU System Processor"). This image
# predates that shift: CUDA 12.4 toolkit + DLVM install-driver.sh that picks proprietary
# (closed) 550.90.07 for n1 machine types — Pascal-compatible. DEPRECATED but READY; verified
# 2026-05-28 on probe VM. See CLAUDE.md "Infrastructure" + Issue #522.
# If Google deletes this image, fallback is a custom bake — see CLAUDE.md.
IMAGE_NAME="pytorch-latest-cu124-v20250327-ubuntu-2204"
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
# Provisioning model for the pipeline VM, applied at CREATE time only (immutable
# afterwards). SPOT (default) is ~60-70% cheaper but preemptible; --on-demand
# forces STANDARD if SPOT P100 capacity is unavailable. An already-created VM
# keeps its original model — delete + re-run to switch.
PROVISIONING_MODEL="SPOT"

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
        --on-demand) PROVISIONING_MODEL="STANDARD";              shift ;;
        --_cloud-internal) CLOUD_INTERNAL=true;                  shift ;;
        *)
            echo "Unknown argument: $1" >&2
            echo "Usage: $0 [--mode test|prod] [--branch <branch>] [--zone <zone>] [--samples <path>] [--detach] [--keep-vm] [--on-demand]" >&2
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
log "  Provision: ${PROVISIONING_MODEL} (applied only when CREATING the VM; existing VMs keep their original model)"
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
        $([[ "${PROVISIONING_MODEL}" == "STANDARD" ]] && echo "--on-demand") \\
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
    # Protect the fresh bucket immediately (versioning + retention lifecycle, Issue #627)
    # so first-deploy results aren't overwritten unrecoverably before someone runs setup.
    log "Applying result-protection (versioning + lifecycle) to the new bucket..."
    bash "$(dirname "$0")/setup_gcs_archiving.sh" "gs://${GCS_BUCKET}" \
        || log "  WARNING: setup_gcs_archiving.sh failed — run it manually to protect this bucket."
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
# Pre-run snapshot of prior results (Issue #658)
# ---------------------------------------------------------------------------
# Object Versioning protects the upcoming overwrite (90-day noncurrent retention),
# but an explicit, named snapshot of the small funnel/classification files gives a
# live, diffable before-state for before/after comparative analyses — no
# generation-number juggling. Large regenerable artifacts (BAMs, AlphaFold
# intermediates) are left to versioning to avoid per-run archive bloat.
# Sentinel = reports/report.tsv (the run's final output): a prior run that failed
# before the report is treated as "no prior results" — its partial junctions stay
# in versioning rather than the archive (a partial-run snapshot has little value).
PRIOR_RESULTS="${GCS_PATH}/${PATIENT_ID}"
if gcloud storage ls "${PRIOR_RESULTS}/reports/report.tsv" &>/dev/null; then
    SNAP_TS="$(date -u +%Y%m%dT%H%M%SZ)"
    SNAP_DEST="${GCS_PATH}/_archive/${PATIENT_ID}_pre_${SNAP_TS}"
    log "Snapshotting prior headline results -> ${SNAP_DEST}/ (Issue #658)..."
    if snap_out="$(gcloud storage cp -r "${PRIOR_RESULTS}/reports" "${PRIOR_RESULTS}/junctions" "${SNAP_DEST}/" 2>&1)"; then
        log "  Pre-run snapshot complete (reports/ + junctions/)."
    else
        log "  WARNING: pre-run snapshot failed — continuing (Object Versioning still protects the overwrite):"
        log "    gcloud: $(printf '%s\n' "${snap_out}" | tail -1)"
    fi
else
    log "No prior results for ${PATIENT_ID} — skipping pre-run snapshot."
fi

# ---------------------------------------------------------------------------
# Create or start pipeline VM
# ---------------------------------------------------------------------------
VM_STATUS="$(vm_status "${PIPELINE_VM}")"
FRESH_BOOT=false
if [[ "${VM_STATUS}" == "TERMINATED" ]]; then
    log "Starting ${PIPELINE_VM}..."
    [[ "${PROVISIONING_MODEL}" == "SPOT" ]] && log "  Note: provisioning model is fixed at creation — this existing VM keeps its original model. To run on SPOT, delete it first: gcloud compute instances delete ${PIPELINE_VM} --zone=${ZONE} --quiet"
    gcloud compute instances start "${PIPELINE_VM}" --zone="${ZONE}"
elif [[ "${VM_STATUS}" == "NOT_FOUND" ]]; then
    log "Pipeline VM ${PIPELINE_VM} not found — creating it (provisioning: ${PROVISIONING_MODEL})..."
    # Provisioning model is set here and is immutable for the VM's lifetime.
    # SPOT adds --instance-termination-action=STOP so a preemption STOPS the VM
    # (preserving the boot disk + in-flight results); the next run resumes via
    # `instances start` + snakemake --rerun-incomplete. Empty array on STANDARD;
    # the ${arr[@]+...} guard keeps the empty case safe under `set -u` (bash 3.2).
    SCHEDULING_ARGS=()
    if [[ "${PROVISIONING_MODEL}" == "SPOT" ]]; then
        SCHEDULING_ARGS=(--provisioning-model=SPOT --instance-termination-action=STOP)
    fi
    gcloud compute instances create "${PIPELINE_VM}" \
        --zone="${ZONE}" \
        --machine-type="${MACHINE_TYPE}" \
        --accelerator="${ACCELERATOR}" \
        --maintenance-policy=TERMINATE \
        --image="${IMAGE_NAME}" \
        --image-project="${IMAGE_PROJECT}" \
        --boot-disk-size="${DISK_SIZE}" \
        --boot-disk-type=pd-balanced \
        --scopes=cloud-platform \
        --metadata=enable-oslogin=FALSE,install-nvidia-driver=True \
        ${SCHEDULING_ARGS[@]+"${SCHEDULING_ARGS[@]}"}
    log "  VM created. Driver install runs on first boot (~2 min)."
    FRESH_BOOT=true
else
    log "${PIPELINE_VM} already running (status: ${VM_STATUS})."
fi

wait_for_ssh "${PIPELINE_VM}"

# ---------------------------------------------------------------------------
# Verify the NVIDIA driver. The DLVM install-driver.sh runs on first boot
# (triggered by `install-nvidia-driver=True` metadata above) and installs
# proprietary closed driver 550.90.07 via DKMS — Pascal-compatible. The
# install takes ~2 min after SSH becomes ready; retry with backoff to
# cover that window. On subsequent VM starts the driver persists via
# DKMS (kernel-module rebuild on kernel upgrade). See CLAUDE.md
# "Infrastructure" + Issue #522.
# ---------------------------------------------------------------------------
# Extend the retry window on fresh boot — driver init can be slowed by
# concurrent unattended-upgrades holding dpkg locks. 36 × 10s = 6 min on
# fresh boot, 18 × 10s = 3 min on warm restart.
MAX_ATTEMPTS=18
[[ "${FRESH_BOOT}" == true ]] && MAX_ATTEMPTS=36
log "Verifying NVIDIA driver on ${PIPELINE_VM} (allow ~$((MAX_ATTEMPTS * 10 / 60)) min for first-boot install)..."
DRIVER_INFO=""
for attempt in $(seq 1 "${MAX_ATTEMPTS}"); do
    if DRIVER_INFO=$(ssh_cmd "${PIPELINE_VM}" -- nvidia-smi --query-gpu=name,driver_version --format=csv,noheader 2>/dev/null | head -1); then
        if [[ -n "${DRIVER_INFO}" ]]; then break; fi
    fi
    sleep 10
done
if [[ -z "${DRIVER_INFO}" ]]; then
    log "ERROR: nvidia-smi failed on ${PIPELINE_VM} after $((MAX_ATTEMPTS * 10 / 60)) min — driver install did not complete."
    log "  Investigate before retrying:"
    log "    gcloud compute ssh ${PIPELINE_VM} --zone=${ZONE} --tunnel-through-iap -- 'sudo journalctl -u google-startup-scripts --no-pager | tail -50; nvidia-smi'"
    exit 1
fi
log "  Driver OK: ${DRIVER_INFO}"

# Assert Pascal compatibility (script invariant: P100 needs driver < 575).
# `nvidia-smi` returning empty already catches the open-driver case (P100
# fails to initialize), but an explicit version check makes the requirement
# self-documenting and fail-fast if a future DLVM script bypasses the n1
# machine-type guard and installs a 575+ open driver that somehow reports.
DRIVER_MAJOR=$(echo "${DRIVER_INFO}" | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 | cut -d. -f1)
if [[ -n "${DRIVER_MAJOR}" && "${DRIVER_MAJOR}" -ge 575 ]]; then
    log "ERROR: Driver ${DRIVER_MAJOR}.x >= 575 dropped Pascal (SM 6.0) support. P100 will not work."
    log "  Pin a different IMAGE_NAME (see top of script) or bake a custom image with driver < 575."
    exit 1
fi

# ---------------------------------------------------------------------------
# Fresh-boot quiet period. The DLVM image runs install-driver.sh in parallel
# with cloud-init + apt unattended-upgrades on first boot. Holding a long
# SSH session through this window gets the connection killed (sshd is
# restarted when unattended-upgrades bumps openssh-server); even short-poll
# SSH attempts can fail mid-handshake. The simplest reliable mitigation is
# to wait host-side for the unattended-upgrades window to close, ~5 min.
# Skipped on warm restarts (VM was TERMINATED then started) — cloud-init +
# unattended-upgrades only fire on the first boot. Caught 2026-05-28
# runs #1 + #2 — both crashed during step 1/8 of setup_vm.sh and during
# the initial cloud-init wait respectively, before this sleep was added.
# ---------------------------------------------------------------------------
if [[ "${FRESH_BOOT}" == true ]]; then
    log "First-boot quiet period (5 min) — letting cloud-init + unattended-upgrades complete..."
    sleep 300
    log "  Quiet period done; sshd should be stable."
fi

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

# Issue #664: drop any pipeline.log left from a prior successful run on this VM.
# Its "steps (100%) done" line would make the completion poller false-complete
# within ~2s (before the new run's tee truncates the file), uploading stale results.
rm -f pipeline.log

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
    # Issue #669: anchor failure detection on Snakemake's real failure signatures,
    # not a bare 'Error' substring. A bare 'Error' matched benign log content —
    # non-fatal tool stderr (tee'd in via 2>&1), conda/solver messages, INFO prose
    # mentioning "error" — aborting healthy runs. 'Error in rule <name>:' marks a
    # failed job; 'Exiting because a job execution failed' is Snakemake's terminal
    # failure line. Both are emitted only on a genuine failure.
    FAILED="$(ssh_cmd "${PIPELINE_VM}" --command \
        "grep -cE 'Error in rule |Exiting because a job execution failed' \$HOME/splice-neoepitope-pipeline/pipeline.log 2>/dev/null || echo 0" \
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
