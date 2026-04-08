#!/usr/bin/env bash
# scripts/setup_local.sh
#
# One-time environment setup for local development on macOS.
#
# What this does:
#   1. Installs Miniforge3 (conda) if not already present
#   2. Creates a 'snakemake' conda environment with Snakemake
#   3. Prints next steps
#
# Usage:
#   bash scripts/setup_local.sh
#
# After this script finishes:
#   bash scripts/prepare_test_data.sh          # download chr22 reference + test FASTQs
#   snakemake --cores 4 --use-conda --configfile config/test_config.yaml
#
set -euo pipefail

echo "=== Splice Neoepitope Pipeline — Local Setup ==="
echo ""

# ---------------------------------------------------------------------------
# 1. Check / install conda (Miniforge3)
# ---------------------------------------------------------------------------
if command -v conda &>/dev/null; then
    echo "[1/2] conda already installed at: $(conda info --base)"
    echo "      Skipping Miniforge3 installation."
else
    echo "[1/2] Installing Miniforge3..."
    # Miniforge uses "MacOSX" where uname returns "Darwin"
    _OS=$(uname)
    _ARCH=$(uname -m)
    [[ "$_OS" == "Darwin" ]] && _OS="MacOSX"
    MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-${_OS}-${_ARCH}.sh"
    MINIFORGE_INSTALLER="/tmp/Miniforge3.sh"

    curl --fail --location --retry 3 --retry-delay 2 --progress-bar \
        "$MINIFORGE_URL" -o "$MINIFORGE_INSTALLER"
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

    # Make conda available in this session
    source "$HOME/miniforge3/etc/profile.d/conda.sh"

    # Add to shell profile for future sessions
    "$HOME/miniforge3/bin/conda" init "$(basename "$SHELL")"

    echo "    Miniforge3 installed at: $HOME/miniforge3"
    echo "    Shell profile updated — open a new terminal after setup completes."
fi

# Ensure conda is on PATH for the rest of this script
if [[ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
fi

# ---------------------------------------------------------------------------
# 2. Create / update the snakemake environment
# ---------------------------------------------------------------------------
echo ""
if conda env list | grep -q "^snakemake "; then
    echo "[2/2] 'snakemake' environment already exists — skipping creation."
else
    echo "[2/2] Creating 'snakemake' conda environment..."
    conda create -n snakemake -c conda-forge -c bioconda \
        "snakemake>=8.0,<9" python=3.11 -y
    echo "    Created 'snakemake' environment."
fi

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
echo ""
echo "=== Setup complete! ==="
echo ""
echo "Activate the snakemake environment:"
echo "    conda activate snakemake"
echo ""
echo "Download test data (one-time, ~15–30 min):"
echo "    bash scripts/prepare_test_data.sh"
echo ""
echo "Run the test pipeline (chr22, 500K reads):"
echo "    snakemake --cores 4 --use-conda --configfile config/test_config.yaml"
echo ""
echo "Note (macOS M1/M2, 8 GB RAM):"
echo "  - Use HISAT2 (default) — STAR requires 32 GB RAM and will not run locally."
echo "  - HISAT2 alignment needs ~8 GB; close other apps before running."
