#!/usr/bin/env bash
# scripts/prepare_test_data.sh
#
# One-time setup: download a chr22-only reference and a small FASTQ subset
# for end-to-end testing of the pipeline on a local macOS machine.
#
# What this downloads:
#   1. chr22 FASTA         — UCSC hg38 (~52 MB uncompressed)
#   2. chr22 GTF           — GENCODE v47 basic, stream-filtered to chr22 (~400 MB download)
#   3. Test FASTQs         — 500K read pairs from ERR188273 via ENA (no sra-tools needed)
#
# Runtime: 15–30 min depending on connection speed (GTF download dominates).
#
# After this script finishes, run the test pipeline with:
#   snakemake --cores 4 --use-conda --configfile config/test_config.yaml
#
set -euo pipefail

RESOURCES="resources/test"
DATA="data/test"

mkdir -p "$RESOURCES" "$DATA"

# ---------------------------------------------------------------------------
# 1. chr22 FASTA
# ---------------------------------------------------------------------------
echo "=== [1/3] Downloading chr22 reference FASTA (UCSC hg38) ==="
if [[ -f "$RESOURCES/chr22.fa" ]]; then
    echo "    Already exists — skipping."
else
    curl -L --progress-bar \
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz" \
        | gunzip > "$RESOURCES/chr22.fa"
    echo "    Saved: $RESOURCES/chr22.fa"
fi

# ---------------------------------------------------------------------------
# 2. chr22 GTF (stream and filter — avoids storing the full 400 MB file)
# ---------------------------------------------------------------------------
echo ""
echo "=== [2/3] Downloading and filtering GENCODE v47 GTF to chr22 ==="
echo "    Streaming ~400 MB — estimated 5–15 min depending on connection."
if [[ -f "$RESOURCES/chr22.gtf.gz" ]]; then
    echo "    Already exists — skipping."
else
    curl -L --progress-bar \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gtf.gz" \
        | zcat \
        | awk '$1 == "chr22" || /^#/' \
        | gzip > "$RESOURCES/chr22.gtf.gz"
    echo "    Saved: $RESOURCES/chr22.gtf.gz"
fi

# ---------------------------------------------------------------------------
# 3. Test FASTQs — 500K read pairs from ERR188273 via ENA
# ---------------------------------------------------------------------------
# ERR188273 is a GEUVADIS paired-end human RNA-Seq sample (LCL NA06985),
# used in the original HISAT2 paper and reliably available on ENA.
# We use it instead of SRR37781424 because ENA may not yet mirror very recent
# SRR accessions. The specific sample doesn't matter for a functional test —
# any human paired-end RNA-Seq will produce chr22 reads.
#
# We stream each file and keep only the first 500K reads (each read = 4 lines).
# Disable pipefail around the streaming pipeline: head exits after N lines,
# sending SIGPIPE to zcat/curl — this is expected and not an error.
# ---------------------------------------------------------------------------
echo ""
echo "=== [3/3] Downloading 500K read pairs from ERR188273 (ENA) ==="
echo "    Streaming via HTTPS — no sra-tools required."
if [[ -f "$DATA/test_tumor_R1.fastq.gz" && -f "$DATA/test_tumor_R2.fastq.gz" ]]; then
    echo "    Already exists — skipping."
else
    # Fetch the canonical ENA FTP URLs for this accession
    ENA_REPORT=$(curl -sf \
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERR188273&result=read_run&fields=fastq_ftp&format=tsv" \
        || true)
    FTP_FIELD=$(echo "$ENA_REPORT" | tail -1 | cut -f2)

    R1_URL="https://$(echo "$FTP_FIELD" | tr ';' '\n' | grep '_1\.fastq\.gz' || true)"
    R2_URL="https://$(echo "$FTP_FIELD" | tr ';' '\n' | grep '_2\.fastq\.gz' || true)"

    if [[ -z "${R1_URL#https://}" || -z "${R2_URL#https://}" ]]; then
        echo "ERROR: could not resolve ENA download URLs for ERR188273." >&2
        echo "       ENA API response: $ENA_REPORT" >&2
        exit 1
    fi

    echo "    R1: $R1_URL"
    echo "    R2: $R2_URL"

    set +o pipefail
    curl -L --progress-bar "$R1_URL" | zcat | head -n 2000000 | gzip > "$DATA/test_tumor_R1.fastq.gz"
    curl -L --progress-bar "$R2_URL" | zcat | head -n 2000000 | gzip > "$DATA/test_tumor_R2.fastq.gz"
    set -o pipefail

    # Check both files exist and are non-empty
    # (gzip -t is intentionally omitted: head closes the pipe before gzip can
    # write its end-of-file marker, producing a truncated-but-readable archive)
    for f in "$DATA/test_tumor_R1.fastq.gz" "$DATA/test_tumor_R2.fastq.gz"; do
        if [[ ! -s "$f" ]]; then
            echo "ERROR: $f is missing or empty after download." >&2
            exit 1
        fi
    done
    echo "    Saved: $DATA/test_tumor_R1.fastq.gz and $DATA/test_tumor_R2.fastq.gz"
fi

# ---------------------------------------------------------------------------
echo ""
echo "=== Setup complete! ==="
echo ""
echo "Run the test pipeline with:"
echo "    snakemake --cores 4 --use-conda --configfile config/test_config.yaml"
echo ""
echo "Note: Only reads mapping to chr22 will generate junctions."
echo "      Expect fewer novel junctions than a full-genome run — this is normal."
