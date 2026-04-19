# =============================================================================
# Rule module: FASTQ download from remote storage (GCS or Backblaze B2)
# =============================================================================
#
# Registers a download_fastq rule when samples.tsv contains remote FASTQ paths.
# Supported URI schemes:
#   gs://   — Google Cloud Storage (gsutil cp)
#   b2://   — Backblaze B2, public buckets (b2 file download, no credentials needed)
#
# Depends on helpers from common.smk: _read_samples_tsv, _local_fastq
# =============================================================================

# Build a one-time lookup: local relative path → remote source URI.
# Keyed on the full local path (data/{patient_id}/{sample_id}/{filename}) so
# two samples with identically-named FASTQs cannot collide.
_REMOTE_FASTQ_MAP = {}
for _row in _read_samples_tsv(config["samples_tsv"]):
    _pid = _row["patient_id"]
    _sid = (_row.get("sample_id") or "").strip()
    for _key in ("fastq1", "fastq2"):
        _path = (_row.get(_key) or "").strip()
        if _path.startswith("gs://") or _path.startswith("b2://"):
            _local = _local_fastq(_path, _pid, _sid)
            _REMOTE_FASTQ_MAP[_local] = _path


if _REMOTE_FASTQ_MAP:
    rule download_fastq:
        """Download a FASTQ from GCS or Backblaze B2 to local data/. Deleted automatically once alignment completes."""
        output:
            fastq=temp("data/{patient_id}/{sample}/{filename}"),
        log:
            os.path.join(_LOGS, "{patient_id}", "download", "{sample}", "{filename}.log"),
        params:
            source_path=lambda wildcards: _REMOTE_FASTQ_MAP[
                os.path.join("data", wildcards.patient_id, wildcards.sample, wildcards.filename)
            ],
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.fastq})
            if [[ "{params.source_path}" == b2://* ]]; then
                b2 file download "{params.source_path}" {output.fastq} 2>&1 | tee {log}
            else
                gsutil cp {params.source_path} {output.fastq} 2>&1 | tee {log}
            fi
            """
