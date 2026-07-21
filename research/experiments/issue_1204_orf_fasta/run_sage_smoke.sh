#!/bin/bash
# research/experiments/issue_1204_orf_fasta/run_sage_smoke.sh
# One-time local loadability smoke: does Sage load the junction-partition FASTA
# concatenated with the canonical proteome and complete a nonspecific search?
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ORF_FASTA="${1:?path to emitted junction_orf.fasta}"
CANON="$HERE/../../../workflow/tests/data/orf_fasta/canonical_tiny.fasta"
cat "$CANON" "$ORF_FASTA" > "$HERE/combined.fasta"
# cd so the config's relative "fasta"/"output_directory" paths resolve here,
# regardless of the caller's cwd or Sage's own relative-path convention.
cd "$HERE"
# Clean a stale output dir so a re-run does not read prior results.
rm -rf sage_out
# sage must be on PATH (brew install sage-proteomics, or the github.com/lazear/sage release binary)
sage sage_config.json synthetic.mgf
echo "SAGE SMOKE OK"
