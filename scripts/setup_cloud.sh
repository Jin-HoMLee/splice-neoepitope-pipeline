#!/usr/bin/env bash
# scripts/setup_cloud.sh
#
# One-time environment setup for a fresh GCP Ubuntu VM.
#
# What this does:
#   1. Installs system dependencies (git, curl, wget)
#   2. Installs Miniforge3 (conda) if not already present
#   3. Creates a 'snakemake' conda environment with Snakemake
#   4. Clones the pipeline repository (if not already present)
#   5. Downloads the full GRCh38 reference FASTA + GENCODE v47 GTF
#
# Usage:
#   Fresh VM (before the repo is cloned — run from your home directory):
#     curl -fsSL https://raw.githubusercontent.com/Jin-HoMLee/splice-neoepitope-pipeline/main/scripts/setup_cloud.sh \
#         -o setup_cloud.sh && bash setup_cloud.sh [--repo-branch <branch>]
#
#   After cloning the repository:
#     bash scripts/setup_cloud.sh [--repo-branch <branch>]
#
# After this script finishes:
#   cd splice-neoepitope-pipeline
#   Edit config/samples.tsv with your sample FASTQ paths
#   conda activate snakemake
#   snakemake --cores $(nproc) --use-conda 2>&1 | tee pipeline.log ; bash auto_stop.sh
#
set -euo pipefail

REPO_URL="https://github.com/Jin-HoMLee/splice-neoepitope-pipeline.git"
REPO_DIR="$HOME/splice-neoepitope-pipeline"
REPO_BRANCH="main"

# Parse optional --repo-branch argument
while [[ $# -gt 0 ]]; do
    case $1 in
        --repo-branch)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --repo-branch requires a value, e.g. --repo-branch main" >&2
                exit 1
            fi
            REPO_BRANCH="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

echo "=== Splice Neoepitope Pipeline — Cloud VM Setup ==="
echo ""

# ---------------------------------------------------------------------------
# 1. System dependencies
# ---------------------------------------------------------------------------
echo "[1/5] Installing system dependencies..."
sudo apt-get update -qq
sudo apt-get install -y git curl wget tmux
echo "    Done."

# ---------------------------------------------------------------------------
# 2. Check / install conda (Miniforge3)
# ---------------------------------------------------------------------------
echo ""
if [[ -f "$HOME/miniforge3/bin/conda" ]]; then
    echo "[2/5] conda already installed at: $HOME/miniforge3"
    echo "      Skipping Miniforge3 installation."
else
    echo "[2/5] Installing Miniforge3..."
    MINIFORGE_INSTALLER="/tmp/Miniforge3.sh"
    curl --fail --location --retry 3 --retry-delay 2 --progress-bar \
        "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" \
        -o "$MINIFORGE_INSTALLER"
    if [[ ! -s "$MINIFORGE_INSTALLER" ]]; then
        echo "Error: Miniforge installer download failed or produced an empty file." >&2
        rm -f "$MINIFORGE_INSTALLER"; exit 1
    fi
    if ! head -n 1 "$MINIFORGE_INSTALLER" | grep -q '^#!/'; then
        echo "Error: Downloaded Miniforge installer does not appear to be a shell script." >&2
        rm -f "$MINIFORGE_INSTALLER"; exit 1
    fi
    bash "$MINIFORGE_INSTALLER" -b -p "$HOME/miniforge3"
    rm "$MINIFORGE_INSTALLER"

    source "$HOME/miniforge3/etc/profile.d/conda.sh"
    "$HOME/miniforge3/bin/conda" init bash
    echo "    Miniforge3 installed at: $HOME/miniforge3"
fi

# Ensure conda is on PATH for the rest of this script
if [[ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
fi

# ---------------------------------------------------------------------------
# 3. Create / update the snakemake environment
# ---------------------------------------------------------------------------
echo ""
if [[ -d "$HOME/miniforge3/envs/snakemake" ]]; then
    echo "[3/5] 'snakemake' environment already exists — skipping creation."
else
    echo "[3/5] Creating 'snakemake' conda environment..."
    conda create -n snakemake -c conda-forge -c bioconda \
        "snakemake>=8.0,<9" python=3.11 -y
    echo "    Created 'snakemake' environment."
fi

# ---------------------------------------------------------------------------
# 4. Clone the repository
# ---------------------------------------------------------------------------
echo ""
if [[ -d "$REPO_DIR/.git" ]]; then
    echo "[4/5] Repository already cloned at $REPO_DIR — pulling latest."
    git -C "$REPO_DIR" fetch origin
    git -C "$REPO_DIR" checkout "$REPO_BRANCH"
    git -C "$REPO_DIR" pull origin "$REPO_BRANCH"
else
    echo "[4/5] Cloning repository (branch: $REPO_BRANCH)..."
    git clone --branch "$REPO_BRANCH" "$REPO_URL" "$REPO_DIR"
fi

# ---------------------------------------------------------------------------
# 5. Download reference data
# ---------------------------------------------------------------------------
echo ""
echo "[5/5] Downloading reference data (this will take 20–60 min)..."

RESOURCES="$REPO_DIR/resources"
mkdir -p "$RESOURCES"

# GRCh38 primary assembly FASTA (uncompressed — required by HISAT2/STAR)
FASTA_GZ="$RESOURCES/GRCh38.primary_assembly.genome.fa.gz"
FASTA="$RESOURCES/GRCh38.primary_assembly.genome.fa"

if [[ -f "$FASTA" ]]; then
    echo "    GRCh38 FASTA already exists — skipping."
else
    echo "    Downloading GRCh38 primary assembly FASTA (~830 MB compressed)..."
    curl --fail -L --progress-bar \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz" \
        -o "$FASTA_GZ"
    echo "    Decompressing FASTA (~3.1 GB uncompressed)..."
    gunzip "$FASTA_GZ"
    if [[ ! -s "$FASTA" ]]; then
        echo "ERROR: FASTA decompression produced an empty file." >&2; exit 1
    fi
    echo "    Saved: $FASTA"
fi

# GENCODE v47 annotation GTF (kept gzipped — pipeline reads with zcat)
GTF="$RESOURCES/gencode.v47.annotation.gtf.gz"

if [[ -f "$GTF" ]]; then
    echo "    GENCODE GTF already exists — skipping."
else
    echo "    Downloading GENCODE v47 annotation GTF (~50 MB compressed)..."
    curl --fail -L --progress-bar \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz" \
        -o "$GTF"
    if ! gzip -t "$GTF" 2>/dev/null; then
        echo "ERROR: GTF download appears corrupt (failed gzip integrity check)." >&2; exit 1
    fi
    echo "    Saved: $GTF"
fi

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
echo ""
echo "=== Setup complete! ==="
echo ""
echo "Next steps:"
echo "  0. Open a new shell (or run: source ~/.bashrc) so that 'conda' is on your PATH"
echo "  1. Edit config/samples.tsv with your sample FASTQ paths and sample types"
echo "  2. Activate the Snakemake environment:"
echo "       conda activate snakemake"
echo "  3. From inside a tmux session, run the pipeline:"
echo "       cd $REPO_DIR"
echo "       snakemake --cores \$(nproc) --use-conda 2>&1 | tee pipeline.log ; bash auto_stop.sh"
echo ""
echo "     The VM will shut down automatically when the pipeline finishes."
echo "     To detach from tmux without stopping the run: Ctrl+B, then D."
