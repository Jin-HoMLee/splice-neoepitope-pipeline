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

# Log-tree convention: per-patient jobs (a {patient_id} wildcard) write to
# logs/<patient_id>/<step>/; shared reference/index/model-build jobs (no
# {patient_id}) write to logs/_shared/<step>/ instead of the top of logs/.
# Namespacing the shared logs keeps every step-name folder at one depth (2)
# and mirrors the results-side segregation (indices/, references/) and the
# research/experiments/_shared/ convention.
_SHARED_LOG = os.path.join(_LOGS, "_shared")

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
        if len(unique_patients) == 0:
            raise ValueError(
                f"samples_tsv '{samples_tsv}' contains no valid sample rows."
            )
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


# ── GTEx pan-tissue blacklist path (Issue #211/#212) ────────────────────────────
# The population-normal junction blacklist lives on GCS in production (a gs://
# reference_bed) and as a small committed/fetched fixture for chr22 tests (a local
# reference_bed set via config/test_config.yaml). The download_gtex_pan_tissue_bed
# rule (download.smk) stages the gs:// object to the local path below; the
# filter_junctions input function (filter_junctions.smk) consumes the SAME path.
# Single-source the derivation here so the download rule's output and the filter
# rule's input cannot drift — a mismatch passes the dry-run clean but fails at
# execute time with a missing-input error (a known dry-run-blind class, CLAUDE.md).

def _gtex_filter_cfg():
    """Return the gtex_filter config block (empty dict when the key is absent)."""
    return config.get("gtex_filter", {}) or {}


def _gtex_blacklist_bed():
    """Local path the GTEx pan-tissue blacklist resolves to, or "" when disabled.

    A gs:// ``reference_bed`` is staged to ``references/gtex/<basename>`` by the
    ``download_gtex_pan_tissue_bed`` rule; a local ``reference_bed`` (the chr22
    test fixture) is consumed in place. Returns "" when the filter is disabled or
    no ``reference_bed`` is set, which makes the filter a no-op (the input function
    then omits the optional input and filter_junctions.py skips GTEx tagging).
    """
    cfg = _gtex_filter_cfg()
    if not cfg.get("enabled", False):
        return ""
    ref = (cfg.get("reference_bed") or "").strip()
    if not ref:
        return ""
    if ref.startswith("gs://"):
        return os.path.join("references", "gtex", os.path.basename(ref))
    return ref
