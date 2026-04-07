#!/usr/bin/env bash
# scripts/prepare_test_data.sh
#
# One-time setup: download a chr22-only reference and a small FASTQ subset
# for end-to-end testing of the pipeline on a local macOS machine.
#
# What this downloads:
#   1. chr22 FASTA         — UCSC hg38 (~52 MB uncompressed)
#   2. chr22 GTF           — GENCODE v47 basic, stream-filtered to chr22 (~400 MB download)
#   3. Test FASTQs         — 500K read pairs from SRR37781424 via fasterq-dump (~200 MB)
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
# 3. Test FASTQs — 500K read pairs from SRR37781424
# ---------------------------------------------------------------------------
echo ""
echo "=== [3/3] Downloading 500K read pairs from SRR37781424 (SRA) ==="
echo "    Requires fasterq-dump (sra-tools). Estimated ~200 MB download."
if [[ -f "$DATA/test_tumor_R1.fastq.gz" ]]; then
    echo "    Already exists — skipping."
else
    fasterq-dump \
        -X 500000 \
        --split-files \
        --threads 4 \
        --outdir "$DATA" \
        --temp "$DATA" \
        SRR37781424

    mv "$DATA/SRR37781424_1.fastq" "$DATA/test_tumor_R1.fastq"
    mv "$DATA/SRR37781424_2.fastq" "$DATA/test_tumor_R2.fastq"
    gzip "$DATA/test_tumor_R1.fastq" "$DATA/test_tumor_R2.fastq"
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
