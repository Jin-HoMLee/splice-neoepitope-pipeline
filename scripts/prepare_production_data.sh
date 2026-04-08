#!/usr/bin/env bash
# scripts/prepare_production_data.sh
#
# Download full FASTQs for the matched gastric cancer tumor/normal pair.
#
# Samples:
#   SRR9143066 — Primary Tumor (gastric cancer surgical section, HiSeq 3000, single-end)
#   SRR9143065 — Solid Tissue Normal (adjacent stomach tissue, HiSeq 3000, single-end)
#
# Files are downloaded via ENA HTTPS — no sra-tools or controlled access required.
# The reference genome (GRCh38 FASTA + GENCODE GTF) is downloaded by setup_cloud.sh.
#
# Usage:
#   bash scripts/prepare_production_data.sh
#
# After this script finishes, run the pipeline with:
#   conda activate snakemake
#   snakemake --cores $(nproc) --use-conda 2>&1 | tee pipeline.log ; bash auto_stop.sh
#
set -euo pipefail

DATA="data"
mkdir -p "$DATA"

TUMOR_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/006/SRR9143066/SRR9143066.fastq.gz"
NORMAL_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR914/005/SRR9143065/SRR9143065.fastq.gz"

echo "=== Downloading production FASTQs ==="
echo "    Tumor:  SRR9143066 (~several GB — may take 30–60 min)"
echo "    Normal: SRR9143065 (~several GB — may take 30–60 min)"
echo ""

if [[ -f "$DATA/gastric_tumor.fastq.gz" ]]; then
    echo "    Tumor FASTQ already exists — skipping."
else
    echo "    Downloading tumor (SRR9143066)..."
    curl --fail -L --progress-bar "$TUMOR_URL" -o "$DATA/gastric_tumor.fastq.gz"
    if ! gzip -t "$DATA/gastric_tumor.fastq.gz" 2>/dev/null; then
        echo "ERROR: Tumor FASTQ download appears corrupt." >&2; exit 1
    fi
    echo "    Saved: $DATA/gastric_tumor.fastq.gz"
fi

echo ""

if [[ -f "$DATA/gastric_normal.fastq.gz" ]]; then
    echo "    Normal FASTQ already exists — skipping."
else
    echo "    Downloading normal (SRR9143065)..."
    curl --fail -L --progress-bar "$NORMAL_URL" -o "$DATA/gastric_normal.fastq.gz"
    if ! gzip -t "$DATA/gastric_normal.fastq.gz" 2>/dev/null; then
        echo "ERROR: Normal FASTQ download appears corrupt." >&2; exit 1
    fi
    echo "    Saved: $DATA/gastric_normal.fastq.gz"
fi

echo ""
echo "=== Downloads complete! ==="
echo ""
echo "Run the pipeline with:"
echo "    conda activate snakemake"
echo "    snakemake --cores \$(nproc) --use-conda 2>&1 | tee pipeline.log ; bash auto_stop.sh"
