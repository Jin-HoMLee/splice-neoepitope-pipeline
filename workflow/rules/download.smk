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

# Build a one-time lookup: local filename → GCS source URI.
# Example: {"202211_..._R1.fastq.gz": "gs://osteosarc-genomics/.../202211_..._R1.fastq.gz"}
_GCS_FASTQ_MAP = {}
for _row in _read_samples_tsv(config["samples_tsv"]):
    for _key in ("fastq1", "fastq2"):
        _path = (_row.get(_key) or "").strip()
        if _path.startswith("gs://"):
            _GCS_FASTQ_MAP[Path(_local_fastq(_path)).name] = _path


if _GCS_FASTQ_MAP:
    rule download_fastq:
        """Copy a FASTQ from GCS to local data/. Deleted automatically once alignment completes."""
        output:
            fastq=temp("data/{filename}"),
        log:
            os.path.join(OUT["logs"], "download", "{filename}.log"),
        params:
            gcs_path=lambda wildcards: _GCS_FASTQ_MAP[wildcards.filename],
        shell:
            "gsutil cp {params.gcs_path} {output.fastq} 2>&1 | tee {log}"
