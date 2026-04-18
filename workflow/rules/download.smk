# =============================================================================
# Rule module: FASTQ download from GCS
# =============================================================================
#
# Registers a download_fastq rule when samples.tsv contains gs:// FASTQ paths.
# Works for both the project staging bucket and publicly accessible foreign
# buckets (e.g. open-access research datasets).
#
# Ingestion from other sources (ENA, SRA, controlled-access) is handled
# outside the pipeline via scripts/ingest_*.sh, which stage FASTQs into a
# GCS bucket before the pipeline runs.
#
# Depends on helpers from common.smk: _read_samples_tsv, _local_fastq
# =============================================================================

# Build a one-time lookup: local relative path → GCS source URI.
# Keyed on the full local path (data/{patient_id}/{sample_id}/{filename}) so
# two samples with identically-named FASTQs cannot collide.
_GCS_FASTQ_MAP = {}
for _row in _read_samples_tsv(config["samples_tsv"]):
    _pid = _row["patient_id"]
    _sid = (_row.get("sample_id") or "").strip()
    for _key in ("fastq1", "fastq2"):
        _path = (_row.get(_key) or "").strip()
        if _path.startswith("gs://"):
            _local = _local_fastq(_path, _pid, _sid)
            _GCS_FASTQ_MAP[_local] = _path


if _GCS_FASTQ_MAP:
    rule download_fastq:
        """Copy a FASTQ from GCS to local data/. Deleted automatically once alignment completes."""
        output:
            fastq=temp("data/{patient_id}/{sample}/{filename}"),
        log:
            os.path.join(_LOGS, "{patient_id}", "download", "{filename}.log"),
        params:
            gcs_path=lambda wildcards: _GCS_FASTQ_MAP[
                os.path.join("data", wildcards.patient_id, wildcards.sample, wildcards.filename)
            ],
        shell:
            "set -euo pipefail && mkdir -p $(dirname {output.fastq}) && gsutil cp {params.gcs_path} {output.fastq} 2>&1 | tee {log}"
