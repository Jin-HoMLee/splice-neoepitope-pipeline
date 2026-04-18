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


# ── Sample manifest reader ────────────────────────────────────────────────────

def _read_samples_tsv(samples_tsv, patient_id=None):
    """Read rows from the samples TSV, skipping blank/comment lines.

    Args:
        samples_tsv: path to the TSV (string or Path)
        patient_id:  if given, return only rows matching this patient

    Returns a list of dicts (one per sample row).
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
    if patient_id is not None:
        rows = [r for r in rows if r["patient_id"] == patient_id]
    return rows


# ── GCS FASTQ path helper ─────────────────────────────────────────────────────
#
# Convention for samples.tsv FASTQ paths:
#   gs://...  — publicly accessible GCS URI (project staging bucket or open
#               foreign bucket); downloaded to data/ as a temp file before
#               alignment by the download_fastq rule in download.smk
#   data/...  — local path; must already exist (test data only)

def _local_fastq(path, patient_id="", sample_id=""):
    """Return the local path for a FASTQ.

    gs:// paths map to data/{patient_id}/{sample_id}/<filename> to avoid
    collisions when two patients/samples share the same FASTQ filename.
    Local paths are returned unchanged.
    """
    if path and path.startswith("gs://"):
        return os.path.join("data", patient_id, sample_id, Path(path).name)
    return path
