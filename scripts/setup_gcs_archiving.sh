#!/usr/bin/env bash
# =============================================================================
# setup_gcs_archiving.sh — protect GCS results from destructive overwrites
# =============================================================================
#
# Issue #627. run_cloud_gpu.sh uploads results with `gcloud storage cp -r`, an
# in-place overwrite. Without object versioning a re-run silently destroys the
# prior results with no recovery. This script enables versioning and applies a
# retention lifecycle so every overwrite leaves a recoverable noncurrent version.
#
# Idempotent: safe to re-run. Reads the policy from scripts/gcs_lifecycle.json.
#
# Lifecycle (see gcs_lifecycle.json):
#   - noncurrent versions under results/<patient>/alignment/  expire after  7 days
#     (regenerable intermediates, ~75% of the footprint)
#   - all other noncurrent versions                           expire after 90 days
#
# MAINTENANCE: GCS lifecycle prefixes cannot wildcard `*/alignment/`, so each new
# patient's `results/<patient>/alignment/` prefix must be added to the 7-day rule
# in gcs_lifecycle.json. Forgetting only costs a little extra storage (old
# alignment versions then fall under the 90-day rule) — never data loss.
#
# Recovery (restore a clobbered object from a noncurrent version):
#   gcloud storage ls --all-versions gs://<bucket>/<path>      # find generation
#   gcloud storage cp gs://<bucket>/<path>#<generation> gs://<bucket>/<path>
#
# Usage:
#   bash scripts/setup_gcs_archiving.sh [gs://bucket-name]
# =============================================================================
set -euo pipefail

BUCKET="${1:-gs://splice-neoepitope-project}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LIFECYCLE_FILE="${SCRIPT_DIR}/gcs_lifecycle.json"

[[ -f "${LIFECYCLE_FILE}" ]] || { echo "ERROR: ${LIFECYCLE_FILE} not found." >&2; exit 1; }

echo "Enabling object versioning on ${BUCKET}..."
gcloud storage buckets update "${BUCKET}" --versioning

echo "Applying retention lifecycle from ${LIFECYCLE_FILE}..."
gcloud storage buckets update "${BUCKET}" --lifecycle-file="${LIFECYCLE_FILE}"

echo "Done. Current config:"
gcloud storage buckets describe "${BUCKET}" \
  --format="yaml(versioning_enabled, lifecycle_config)"
