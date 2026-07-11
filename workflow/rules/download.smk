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
        vdjdb_tsv = f"references/vdjdb/{config['tcrdock']['vdjdb_release']}/vdjdb_full.txt",
        sentinel = f"references/vdjdb/{config['tcrdock']['vdjdb_release']}/.download.done",
    log:
        os.path.join(_SHARED_LOG, "download", f"vdjdb_{config['tcrdock']['vdjdb_release']}.log"),
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


# =============================================================================
# GTEx pan-tissue blacklist download — Issue #211/#212
# =============================================================================
# Stage the genome-wide GTEx pan-tissue junction blacklist (Snaptron gtexv2) from
# GCS to a local path for the filter_junctions step. Registered ONLY when the
# filter is enabled AND reference_bed is a gs:// URI — the chr22 test fixture is a
# local path (config/test_config.yaml, prepared by scripts/prepare_test_data.sh)
# and needs no download rule. Relies on HOST gsutil (gcloud SDK): no conda env
# declares gsutil, and this rule only ever fires in production on the VM where
# gcloud exists.

_GTEX_CFG = config.get("gtex_filter", {}) or {}
_GTEX_REF = (_GTEX_CFG.get("reference_bed") or "").strip()
if _GTEX_CFG.get("enabled", False) and _GTEX_REF.startswith("gs://"):

    rule download_gtex_pan_tissue_bed:
        """Download the GTEx pan-tissue junction blacklist BED from GCS (Issue #211/#212).

        Idempotent: with no inputs, Snakemake won't re-run once the local BED exists.
        Optionally SHA256-verified against gtex_filter.reference_bed_sha256 — this
        guards against the silent truncation that bit the #211 genome-wide build (a
        partial copy would otherwise pass downstream unnoticed). Uses host gsutil."""
        output:
            bed=_gtex_blacklist_bed(),
        params:
            src=_GTEX_REF,
            sha256=(_GTEX_CFG.get("reference_bed_sha256") or "").strip(),
        log:
            os.path.join(_SHARED_LOG, "download", "gtex_pan_tissue_bed.log"),
        shell:
            """
            set -euo pipefail
            mkdir -p "$(dirname {output.bed})"
            echo "Downloading GTEx pan-tissue blacklist from {params.src}" > {log} 2>&1
            gsutil cp {params.src} {output.bed} >> {log} 2>&1
            if [ -n "{params.sha256}" ]; then
                echo "Verifying SHA256..." >> {log} 2>&1
                ACTUAL=$(shasum -a 256 {output.bed} | awk '{{print $1}}')
                if [ "$ACTUAL" != "{params.sha256}" ]; then
                    echo "SHA256 mismatch! expected={params.sha256} actual=$ACTUAL" >> {log} 2>&1
                    rm -f {output.bed}
                    exit 1
                fi
                echo "SHA256 OK" >> {log} 2>&1
            fi
            echo "Done." >> {log} 2>&1
            """


# =============================================================================
# IMGT germline download via stitchrdl — Issue #204
# =============================================================================

rule download_rmsk_chrom:
    """Fetch one chromosome's UCSC RepeatMasker track as BED (Issue #919).

    Used to categorize junctions by repeat overlap when evaluating the
    NH-uniqueness prefilter (alignment.uniqueness_filter): the filter claims to
    remove multimapper-driven calls at repeat copies, and this is what lets us
    check that the junctions it removes are in fact repeat-overlapping.

    Not part of the default DAG - request the BED explicitly:

        snakemake --cores 1 -- references/rmsk/hg38/rmsk.chr22.bed

    Lands in references/ (gitignored) rather than resources/test/, because it is
    downloaded reference data, not a committed fixture. Idempotent: no inputs, so
    Snakemake will not re-fetch once the BED exists. No conda env - the fetcher is
    stdlib-only, matching the sibling download rules above. The fetcher fails
    loudly on a truncated UCSC response rather than writing a short BED that would
    silently understate repeat overlap.
    """
    output:
        bed="references/rmsk/{ucsc_genome}/rmsk.{chrom}.bed",
    log:
        os.path.join(_SHARED_LOG, "download", "rmsk_{ucsc_genome}_{chrom}.log"),
    wildcard_constraints:
        ucsc_genome=r"[A-Za-z0-9]+",
        chrom=r"chr[0-9A-Za-z_]+",
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        python workflow/scripts/fetch_rmsk.py \
            --genome {wildcards.ucsc_genome} \
            --chrom {wildcards.chrom} \
            --output {output.bed} \
            > {log} 2>&1
        """


rule download_imgt_germlines:
    """Download IMGT germline reference data via stitchrdl.

    Sentinel-gated for idempotency. stitchrdl pulls latest IMGT release; the
    rule pins by sentinel only (IMGT itself doesn't expose tagged releases —
    see Known limitations in the design spec).
    """
    output:
        sentinel = "references/imgt_germlines/.download.done",
    log:
        os.path.join(_SHARED_LOG, "download", "imgt_germlines.log"),
    conda:
        "../envs/vdjdb.yaml"
    shell:
        """
        set -euo pipefail
        DIR=$(dirname {output.sentinel})
        mkdir -p "$DIR"
        echo "Downloading IMGT germline data via stitchrdl..." >> {log} 2>&1
        # stitchrdl writes to its default cache; we copy into references/ for reproducibility
        stitchrdl --species HUMAN >> {log} 2>&1
        # Discover stitchr's cache dir from python (preferred over hard-coding)
        STITCHR_CACHE=$(python -c "import IMGTgeneDL; from pathlib import Path; print(Path(IMGTgeneDL.__file__).parent / 'data' / 'HUMAN')")
        if [ -d "$STITCHR_CACHE" ]; then
            cp -r "$STITCHR_CACHE" "$DIR/HUMAN"
            echo "Copied IMGT HUMAN data from $STITCHR_CACHE to $DIR/HUMAN" >> {log} 2>&1
        else
            echo "WARNING: stitchr cache not found at $STITCHR_CACHE — reproducibility copy skipped." >> {log} 2>&1
            echo "WARNING: pipeline still works via stitchr's own cache (fetch_vdjdb_panel reads from there);" >> {log} 2>&1
            echo "WARNING: but $DIR/HUMAN will be EMPTY. Investigate IMGTgeneDL package layout if you need the local copy for audit." >> {log} 2>&1
        fi
        touch {output.sentinel}
        echo "Done." >> {log} 2>&1
        """
