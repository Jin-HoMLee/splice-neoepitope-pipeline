#!/usr/bin/env bash
# Issue #1048 - pVACtools on MHCflurry-only (class I, zero NetMHC). Feasibility spike.
#
# Builds an open-only pVACtools env and runs pvacbind + pvacseq with MHCflurry as the
# sole prediction algorithm. See README.md for the verdict and the gotcha table.
#
# Usage:  bash reproduce.sh [workdir]     (default workdir: ./.spike_run, gitignored)
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${1:-$HERE/.spike_run}"
VENV="$WORK/pvac_venv"
mkdir -p "$WORK"

# --- Gotcha 1: this Python must be built with Tk. pvactools imports turtle -> tkinter
# at package import (pvactools/lib/vector_visualization.py). pyenv 3.11.5 and the
# homebrew 3.11/3.13 builds lack _tkinter; pyenv 3.10.13 has it.
PYBIN="${PYBIN:-$HOME/.pyenv/versions/3.10.13/bin/python}"
"$PYBIN" -c "import _tkinter" 2>/dev/null \
  || { echo "FATAL: $PYBIN lacks _tkinter; pvactools cannot even be imported. Set PYBIN to a Tk-enabled Python." >&2; exit 1; }

# --- Gotcha 2 (macOS): clang++ here cannot find <vector>. The CommandLineTools
# include dir exists but omits the libc++ headers, so point the build at the SDK copy.
# vaxrank -> shellinford is a C++ extension and fails to build without this.
if [[ "$(uname -s)" == "Darwin" ]]; then
  SDK="$(xcrun --show-sdk-path)"
  export SDKROOT="$SDK"
  export CPPFLAGS="-isysroot $SDK -I$SDK/usr/include/c++/v1"
  export CXXFLAGS="$CPPFLAGS"
  export CFLAGS="-isysroot $SDK"
fi

command -v uv >/dev/null || { echo "FATAL: uv not found." >&2; exit 1; }
uv venv --seed --python "$PYBIN" "$VENV"
PY="$VENV/bin/python"

# --- Gotcha 3: install pvactools WITHOUT its declared deps, then supply them minus
# mhcnuggets. pvactools hard-requires mhcnuggets==2.4.1 (JHU NON-COMMERCIAL), which
# hard-requires tensorflow. Dropping it keeps MHCflurry/MHCflurryEL selectable and
# keeps the environment open-only. h5py must then be added explicitly: pvactools
# imports it in lib/output_parser.py but normally gets it via the tensorflow stack.
uv pip install --python "$PY" --no-deps pvactools==7.0.1
uv pip install --python "$PY" \
  'vcfpy==0.13.8' requests 'PyYAML>=5.1' 'biopython==1.77' networkx simanneal 'numpy==1.26.4' \
  'pandas<2.1.0' wget pysam Pillow pymp-pypi mock 'vaxrank>=1.1.0' 'varcode>=1.1.0' \
  testfixtures 'gtfparse==2.0.1' 'pyfaidx>=0.7.1' 'fsspec<=2025.3.0' packaging pyarrow \
  'polars==0.16.18' XlsxWriter openpyxl deepdiff 'scikit-learn>=1.6.0,<1.8.0' imblearn \
  h5py 'setuptools<81'

# --- Gotcha 4: override pvactools' mhcflurry==2.0.6 pin. 2.0.6 lazily imports
# tensorflow at PREDICT time (undeclared runtime dep) so it installs fine and then
# crashes. 2.2.x is torch-backed and is what our pipeline runs. pvactools shells out
# to the `mhcflurry-predict` CLI, whose contract is unchanged, so the override is safe.
# --- Gotcha 5: setuptools<81 above, because mhcflurry still imports pkg_resources.
# Bounded upward on purpose: the mhcflurry version is load-bearing here (2.0.6 is the
# exact thing being dodged), and a future major could change the `mhcflurry-predict`
# CLI contract or its output column names and silently break reproduction.
uv pip install --python "$PY" 'mhcflurry>=2.1,<3'

# --- Gotcha 6: pvactools invokes the BARE name `mhcflurry-predict` via PATH. Calling
# pvacbind by absolute path is not enough; without this the run dies with
# "An error occurred while calling MHCflurry:" and an EMPTY message (swallowed
# FileNotFoundError). Put the venv bin on PATH.
export PATH="$VENV/bin:$PATH"

# Preconditions for the open-only claim: assert the blocked stack is genuinely absent.
for b in netMHCpan netMHC netMHCIIpan netMHCcons netMHCstabpan; do
  command -v "$b" >/dev/null && { echo "FATAL: $b on PATH; this is not an open-only env." >&2; exit 1; }
done
"$PY" -c "import mhcnuggets" 2>/dev/null && { echo "FATAL: mhcnuggets importable." >&2; exit 1; }
"$PY" -c "import tensorflow" 2>/dev/null && { echo "FATAL: tensorflow importable." >&2; exit 1; }
echo "OK: no NetMHC binaries, no mhcnuggets, no tensorflow."

mhcflurry-downloads fetch models_class1_presentation

# --- AC3: pvacbind on the seeded peptide fixture (3 known A*02:01 epitopes).
rm -rf "$WORK/out_pvacbind"
pvacbind run "$HERE/fixtures/peptides.fa" spike1048 "HLA-A*02:01" MHCflurry "$WORK/out_pvacbind" -e1 9

# --- AC2: pvacseq on pVACtools' own VEP-annotated example VCF.
# download_example_data refuses to write into an existing directory, so clear it first.
rm -rf "$WORK/example"
pvacseq download_example_data "$WORK/example"
rm -rf "$WORK/out_pvacseq"
pvacseq run "$WORK/example/pvacseq_example_data/annotated.expression.vcf.gz" \
  HCC1395_TUMOR_DNA "HLA-A*02:01" MHCflurry "$WORK/out_pvacseq" -e1 9 \
  --normal-sample-name HCC1395_NORMAL_DNA

echo
echo "=== committed-vs-fresh comparison ==="
# By default this is ADVISORY: it reports MATCH/DIFFER but does not fail the run, because
# float jitter across MHCflurry/torch builds and platforms is expected and is not a defect.
# A zero exit therefore means "both modules ran open-only", NOT "outputs were byte-identical".
# Set STRICT=1 to make divergence a hard failure (that is the mode used to substantiate the
# byte-identical claim in README.md, on the authoring machine).
divergences=0
for pair in \
  "out_pvacbind/MHC_Class_I/spike1048.MHC_I.filtered.tsv:outputs/pvacbind.MHC_I.filtered.tsv" \
  "out_pvacseq/MHC_Class_I/HCC1395_TUMOR_DNA.MHC_I.filtered.tsv:outputs/pvacseq.MHC_I.filtered.tsv"; do
  fresh="$WORK/${pair%%:*}"; committed="$HERE/${pair##*:}"
  if diff -q "$fresh" "$committed" >/dev/null 2>&1; then
    echo "  MATCH  ${pair##*:}"
  else
    echo "  DIFFER ${pair##*:} (float jitter across MHCflurry/torch builds is expected; inspect before alarm)"
    divergences=$((divergences + 1))
  fi
done

if [[ "${STRICT:-0}" == "1" && "$divergences" -gt 0 ]]; then
  echo "FATAL: STRICT=1 and $divergences output(s) diverged from the committed copies." >&2
  exit 1
fi

echo
echo "GREEN: pvacbind + pvacseq both completed class-I predictions on MHCflurry alone."
if [[ "$divergences" -gt 0 ]]; then
  echo "NOTE: $divergences output(s) differed from the committed copies (advisory; re-run with STRICT=1 to gate on this)."
fi
