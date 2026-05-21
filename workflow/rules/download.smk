# =============================================================================
# Rule module: FASTQ download from remote storage (GCS or public HTTPS)
# =============================================================================
#
# Registers a download_fastq rule when samples.tsv contains remote FASTQ paths.
# Supported URI schemes:
#   gs://    — Google Cloud Storage (gsutil cp)
#   https:// — public HTTPS URL, e.g. Backblaze B2 custom domain (curl)
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
        if any(_path.startswith(s) for s in _REMOTE_SCHEMES):
            _local = _local_fastq(_path, _pid, _sid)
            _REMOTE_FASTQ_MAP[_local] = _path


if _REMOTE_FASTQ_MAP:
    rule download_fastq:
        """Download a FASTQ from GCS or public HTTPS to local data/. Deleted automatically once alignment completes."""
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
            if [[ "{params.source_path}" == https://* ]]; then
                curl --fail -L --no-progress-meter -o {output.fastq} "{params.source_path}" 2>&1 | tee {log}
            else
                gsutil cp {params.source_path} {output.fastq} 2>&1 | tee {log}
            fi
            """


# =============================================================================
# VDJdb release download — Issue #204
# =============================================================================

rule download_vdjdb_release:
    """Download a pinned VDJdb release, verify SHA256, extract vdjdb_full.txt.

    Sentinel-gated for idempotency. Re-runs are no-ops once the sentinel exists.
    """
    output:
        vdjdb_tsv = f"resources/vdjdb/{config['tcrdock']['vdjdb_release']}/vdjdb_full.txt",
        sentinel = f"resources/vdjdb/{config['tcrdock']['vdjdb_release']}/.download.done",
    log:
        f"logs/download/vdjdb_{config['tcrdock']['vdjdb_release']}.log",
    params:
        release = config["tcrdock"]["vdjdb_release"],
        sha256 = config["tcrdock"]["vdjdb_sha256"],
    conda:
        "../envs/vdjdb.yaml"
    shell:
        """
        set -euo pipefail
        DIR=$(dirname {output.sentinel})
        mkdir -p "$DIR"
        ZIP="$DIR/vdjdb-{params.release}.zip"
        URL="https://github.com/antigenomics/vdjdb-db/releases/download/{params.release}/vdjdb-{params.release}.zip"

        echo "Downloading VDJdb {params.release} from $URL" >> {log} 2>&1
        curl -fsSL "$URL" -o "$ZIP" >> {log} 2>&1

        echo "Verifying SHA256..." >> {log} 2>&1
        ACTUAL=$(shasum -a 256 "$ZIP" | awk '{{print $1}}')
        if [ "$ACTUAL" != "{params.sha256}" ]; then
            echo "SHA256 mismatch! expected={params.sha256} actual=$ACTUAL" >> {log} 2>&1
            exit 1
        fi
        echo "SHA256 OK" >> {log} 2>&1

        echo "Extracting..." >> {log} 2>&1
        unzip -o "$ZIP" -d "$DIR/extracted" >> {log} 2>&1
        cp "$DIR/extracted/vdjdb-{params.release}/vdjdb_full.txt" {output.vdjdb_tsv}

        touch {output.sentinel}
        echo "Done." >> {log} 2>&1
        """
