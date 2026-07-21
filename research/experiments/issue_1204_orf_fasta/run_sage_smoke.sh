#!/bin/bash
# research/experiments/issue_1204_orf_fasta/run_sage_smoke.sh
# AC-3 recovery smoke: build the search DB (canonical proteome + our junction
# ORF partition), generate a theoretical spectrum for a real junction peptide
# from the current FASTA, run a nonspecific Sage search, and ASSERT that the
# junction peptide is recovered (a check that can fail). Not a CI check - the
# CI-side check is the structural test in workflow/tests/.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ORF_FASTA="${1:?path to emitted junction_orf.fasta}"
# Resolve to absolute before any cd (the cat + generator run after cd $HERE).
ORF_ABS="$(cd "$(dirname "$ORF_FASTA")" && pwd)/$(basename "$ORF_FASTA")"
CANON="$HERE/../../../workflow/tests/data/orf_fasta/canonical_tiny.fasta"

# cd so the config's relative "fasta"/"output_directory" paths resolve here.
cd "$HERE"
# Clean prior run artifacts so a re-run cannot read stale results.
rm -rf sage_out combined.fasta synthetic.mgf

# Competing-target DB: canonical proteome + our junction ORF partition.
cat "$CANON" "$ORF_ABS" > combined.fasta

# Generate a spectrum for a real junction peptide from the current FASTA.
TARGET="$(python3 make_spectrum.py --fasta "$ORF_ABS" --out synthetic.mgf)"
echo "target junction peptide: $TARGET"

# sage must be on PATH (github.com/lazear/sage v0.14.7 release binary).
sage sage_config.json synthetic.mgf

# Recovery assertion: the junction peptide must appear in the PSM table
# (column 2 = peptide). Exact whole-field match.
if tail -n +2 sage_out/results.sage.tsv | cut -f2 | grep -qx "$TARGET"; then
  echo "SAGE SMOKE OK - recovered junction peptide $TARGET"
else
  echo "SAGE SMOKE FAILED - $TARGET not recovered in sage_out/results.sage.tsv" >&2
  exit 1
fi
