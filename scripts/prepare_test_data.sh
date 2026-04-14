#!/usr/bin/env bash
# scripts/prepare_test_data.sh
#
# One-time setup: download a chr22-only reference and a small FASTQ subset
# for end-to-end testing of the pipeline on a local macOS machine.
#
# What this downloads:
#   1. chr22 FASTA         — UCSC hg38 (~52 MB uncompressed)
#   2. chr22 GTF           — GENCODE v47 basic, stream-filtered to chr22 (~400 MB download)
#   3. Test FASTQs         — 500K reads each from a matched gastric cancer tumor/normal pair
#                            (SRR9143066 tumor, SRR9143065 normal) via ENA HTTPS
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
# 3. Test FASTQs — matched gastric cancer tumor/normal pair via ENA
# ---------------------------------------------------------------------------
# SRR9143066: gastric cancer surgical section (Primary Tumor)
# SRR9143065: normal stomach tissue adjacent to tumor (Solid Tissue Normal)
# Both are Illumina HiSeq 3000, single-end reads.
# Source: ENA, freely available without controlled-access requirements.
#
# We stream each file and keep only the first 500K reads (4 lines per read).
# Disable pipefail around the streaming pipeline: head exits after N lines,
# sending SIGPIPE to zcat/curl — this is expected and not an error.
# (gzip -t is intentionally omitted: head closes the pipe before gzip can
# write its end-of-file marker, producing a truncated-but-readable archive)
# ---------------------------------------------------------------------------
echo ""
echo "=== [3/3] Downloading 500K reads each from SRR9143066 (tumor) and SRR9143065 (normal) ==="
echo "    Streaming via ENA HTTPS — no sra-tools required."

TUMOR_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/006/SRR9143066/SRR9143066.fastq.gz"
NORMAL_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/005/SRR9143065/SRR9143065.fastq.gz"

# Destination paths — must match fastq1/fastq2 in config/test_samples.tsv
TUMOR_DEST="$DATA/SRR9143066_test.fastq.gz"
NORMAL_DEST="$DATA/SRR9143065_test.fastq.gz"

# Migrate old filenames if present
[[ -f "$DATA/test_tumor.fastq.gz"  && ! -f "$TUMOR_DEST"  ]] && mv "$DATA/test_tumor.fastq.gz"  "$TUMOR_DEST"
[[ -f "$DATA/test_normal.fastq.gz" && ! -f "$NORMAL_DEST" ]] && mv "$DATA/test_normal.fastq.gz" "$NORMAL_DEST"

if [[ -f "$TUMOR_DEST" ]]; then
    echo "    Tumor FASTQ already exists — skipping."
else
    echo "    Downloading tumor (SRR9143066)..."
    set +o pipefail
    curl -L --progress-bar "$TUMOR_URL" | zcat | head -n 2000000 | gzip > "$TUMOR_DEST"
    set -o pipefail
    if [[ ! -s "$TUMOR_DEST" ]]; then
        echo "ERROR: $TUMOR_DEST is missing or empty after download." >&2
        exit 1
    fi
    echo "    Saved: $TUMOR_DEST"
fi

if [[ -f "$NORMAL_DEST" ]]; then
    echo "    Normal FASTQ already exists — skipping."
else
    echo "    Downloading normal (SRR9143065)..."
    set +o pipefail
    curl -L --progress-bar "$NORMAL_URL" | zcat | head -n 2000000 | gzip > "$NORMAL_DEST"
    set -o pipefail
    if [[ ! -s "$NORMAL_DEST" ]]; then
        echo "ERROR: $NORMAL_DEST is missing or empty after download." >&2
        exit 1
    fi
    echo "    Saved: $NORMAL_DEST"
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
