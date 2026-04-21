# =============================================================================
# Rule module: shared helpers — included first by the Snakefile
# =============================================================================
#
# All included .smk files share one Snakemake namespace, so helpers defined
# here are available everywhere without re-importing or duplicating.
#
# =============================================================================

import csv
import os
from pathlib import Path

# ── Config validation ─────────────────────────────────────────────────────────
if "samples_tsv" not in config:
    raise ValueError(
        "samples_tsv is required. Pass it at runtime:\n"
        "  snakemake --config samples_tsv=config/samples/<patient>.tsv\n"
        "  or set it in a configfile (e.g. config/test_config.yaml)."
    )

# ── Config aliases ────────────────────────────────────────────────────────────
_RES  = config["output"]["results_dir"]
_LOGS = config["output"]["logs"]

# ── Wildcard constraints ──────────────────────────────────────────────────────
wildcard_constraints:
    patient_id="[^/]+",
    sample="[^/]+",

# ── Sample manifest reader ────────────────────────────────────────────────────

def _read_samples_tsv(samples_tsv, patient_id=None):
    """Read rows from the samples TSV, skipping blank/comment lines.

    Args:
        samples_tsv: path to the TSV (string or Path)
        patient_id:  if given, return only rows matching this patient

    Returns a list of dicts (one per sample row).

    Raises ValueError if the TSV contains more than one patient_id and no
    patient_id filter is given — each run must target exactly one patient.
    """
    samples_path = Path(samples_tsv)
    if not samples_path.exists():
        return []
    rows = []
    with samples_path.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pid = (row.get("patient_id") or "").strip()
            if not pid or pid.startswith("#"):
                continue
            row["patient_id"] = pid
            rows.append(row)
    if patient_id is None:
        unique_patients = {r["patient_id"] for r in rows}
        if len(unique_patients) > 1:
            raise ValueError(
                f"samples_tsv must contain exactly one patient_id per run, "
                f"found {len(unique_patients)}: {', '.join(sorted(unique_patients))}.\n"
                "Run the pipeline separately for each patient."
            )
    else:
        rows = [r for r in rows if r["patient_id"] == patient_id]
    return rows


# ── GCS FASTQ path helper ─────────────────────────────────────────────────────
#
# Convention for samples.tsv FASTQ paths:
#   gs://...   — Google Cloud Storage URI; downloaded to data/ by download_fastq
#   https://…  — public HTTPS URL (e.g. Backblaze B2 custom domain); curl download
#   data/...   — local path; must already exist (test data only)

_REMOTE_SCHEMES = ("gs://", "https://")

def _local_fastq(path, patient_id="", sample_id=""):
    """Return the local path for a FASTQ.

    Remote paths (gs://, https://) map to data/{patient_id}/{sample_id}/<filename>
    to avoid collisions when two patients/samples share the same FASTQ filename.
    Local paths are returned unchanged.
    """
    if path and any(path.startswith(s) for s in _REMOTE_SCHEMES):
        return os.path.join("data", patient_id, sample_id, Path(path).name)
    return path
