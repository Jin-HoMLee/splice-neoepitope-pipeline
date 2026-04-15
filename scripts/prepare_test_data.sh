#!/usr/bin/env bash
# scripts/prepare_test_data.sh
#
# One-time setup: download a chr22+chr6 reference and a small FASTQ subset
# for end-to-end testing of the pipeline on a local macOS machine.
#
# What this downloads:
#   1. chr22 + chr6 FASTA  — UCSC hg38 (chr22 ~52 MB, chr6 ~58 MB compressed)
#   2. chr22 + chr6 GTF    — GENCODE v47 basic, stream-filtered (~400 MB download)
#   3. Test FASTQs         — 500K reads each from a matched gastric cancer tumor/normal pair
#                            (SRR9143066 tumor, SRR9143065 normal) via ENA HTTPS
#
# chr6 is included so that HLA typing (arcasHLA) can be tested end-to-end:
# the test FASTQs are whole-transcriptome and contain HLA reads, but they can
# only map if chr6 is present in the HISAT2 index.
#
# Runtime: 20–40 min depending on connection speed (GTF download dominates).
#
# After this script finishes, run the test pipeline with:
#   snakemake --cores 4 --use-conda --configfile config/test_config.yaml
#
set -euo pipefail

RESOURCES="resources/test"
DATA="data/test"

mkdir -p "$RESOURCES" "$DATA"

# ---------------------------------------------------------------------------
# 1. chr22 + chr6 FASTA
# ---------------------------------------------------------------------------
echo "=== [1/4] Downloading chr22 + chr6 reference FASTA (UCSC hg38) ==="
if [[ -f "$RESOURCES/chr22.fa" ]]; then
    echo "    chr22 already exists — skipping."
else
    curl -L --progress-bar \
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz" \
        | gunzip > "$RESOURCES/chr22.fa"
    echo "    Saved: $RESOURCES/chr22.fa"
fi

if [[ -f "$RESOURCES/chr6.fa" ]]; then
    echo "    chr6 already exists — skipping."
else
    curl -L --progress-bar \
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr6.fa.gz" \
        | gunzip > "$RESOURCES/chr6.fa"
    echo "    Saved: $RESOURCES/chr6.fa"
fi

# Combine into a single reference FASTA for HISAT2 indexing
if [[ -f "$RESOURCES/chr22_chr6.fa" ]]; then
    echo "    Combined chr22+chr6 FASTA already exists — skipping."
else
    cat "$RESOURCES/chr22.fa" "$RESOURCES/chr6.fa" > "$RESOURCES/chr22_chr6.fa"
    echo "    Saved: $RESOURCES/chr22_chr6.fa"
fi

# ---------------------------------------------------------------------------
# 2. chr22 + chr6 GTF (stream and filter — avoids storing the full 400 MB file)
# ---------------------------------------------------------------------------
echo ""
echo "=== [2/4] Downloading and filtering GENCODE v47 GTF to chr22 + chr6 ==="
echo "    Streaming ~400 MB — estimated 5–15 min depending on connection."
if [[ -f "$RESOURCES/chr22_chr6.gtf.gz" ]]; then
    echo "    Already exists — skipping."
else
    curl -L --progress-bar \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gtf.gz" \
        | zcat \
        | awk '$1 == "chr22" || $1 == "chr6" || /^#/' \
        | gzip > "$RESOURCES/chr22_chr6.gtf.gz"
    echo "    Saved: $RESOURCES/chr22_chr6.gtf.gz"
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
echo "=== [3/4] Downloading 500K reads each from SRR9143066 (tumor) and SRR9143065 (normal) ==="
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
# 4. Invalidate cached HISAT2 index if reference was updated
# ---------------------------------------------------------------------------
echo ""
echo "=== [4/4] Checking HISAT2 index ==="
INDEX_DIR="resources/test/hisat2_index"
# The index references chr22_chr6.fa; if it was built from chr22.fa alone,
# delete it so Snakemake rebuilds it with the combined reference.
if [[ -f "$INDEX_DIR/index.done" ]]; then
    if grep -q "chr22_chr6" "$INDEX_DIR/index.done" 2>/dev/null; then
        echo "    HISAT2 index already built from chr22+chr6 — skipping."
    else
        echo "    Removing stale chr22-only HISAT2 index — will be rebuilt by Snakemake."
        rm -rf "$INDEX_DIR"
    fi
else
    echo "    No existing HISAT2 index found — Snakemake will build it."
fi

# ---------------------------------------------------------------------------
echo ""
echo "=== Setup complete! ==="
echo ""
echo "Run the test pipeline with:"
echo "    snakemake --cores 4 --use-conda --configfile config/test_config.yaml"
echo ""
echo "Note: Reads mapping to chr22 generate novel splice junctions."
echo "      Reads mapping to chr6 enable HLA typing via arcasHLA."
echo "      Expect fewer novel junctions than a full-genome run — this is normal."
