#!/usr/bin/env bash
# Issue #965 leaf A - ASNEO open-only local smoke (macOS arm64, CPU).
# Runs the patched (option-B, no netMHCpan) ASNEO caller on ASNEO's OWN bundled
# test SJ.out.tab, subset to chr22, against the small hg19 chr22 FASTA - no STAR,
# no hand-crafted input. Proves the open-only ASNEO caller runs end to end on
# arm64 CPU and emits junction-derived candidate peptides (the artifact that would
# feed MHCflurry). Reuses the #566 asneo conda env + option-B patch.
#
# Run from repo root:  conda activate asneo  (or rely on `conda run -n asneo`)
#   bash research/experiments/issue_965_open_caller_audit/install/asneo_smoke.sh <WORKDIR>
set -euo pipefail

REPO="$(git rev-parse --show-toplevel)"
WORK="${1:-$REPO/.asneo_smoke_work}"          # scratch (clone + genome + outputs); keep OUT of the repo tree
PATCH="$REPO/research/experiments/issue_566_asneo_crosscheck/apply_optionB_patch.py"
mkdir -p "$WORK"

# 1. ASNEO clone (ships its own test/SRR2660032.SJ.out.tab).
[ -d "$WORK/ASNEO" ] || git clone --depth 1 https://github.com/bm2-lab/ASNEO.git "$WORK/ASNEO"

# 2. hg19 chr22 FASTA (UCSC) - small single-chr genome, fits 8 GB.
if [ ! -f "$WORK/chr22.fa" ]; then
  curl -sL "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz" -o "$WORK/chr22.fa.gz"
  gunzip -f "$WORK/chr22.fa.gz"
fi

# 3. chr22-scale junction table from ASNEO's own test data (a true smoke, not full scale).
awk '$1=="chr22"' "$WORK/ASNEO/test/SRR2660032.SJ.out.tab" > "$WORK/chr22_SJ.out.tab"

# 4. Option-B patch (strip/bypass bundled netMHCpan; stop at candidate peptides).
conda run -n asneo python "$PATCH" "$WORK/ASNEO/ASNEO.py"

# 5. Run the patched caller. Default thresholds yield ~0 peptides at chr22 coverage;
#    relaxed thresholds exercise the full translate -> peptide -> normal-subtract tail.
cd "$WORK/ASNEO"
conda run -n asneo python ASNEO.py -j "$WORK/chr22_SJ.out.tab" -a HLA-A02:01 \
  -l 8,9,10,11 --reads 2 --psi 0.05 -g "$WORK/chr22.fa" -o "$WORK/out"

echo "candidate peptides:"; wc -l "$WORK/out"/*/putative_peptide.txt 2>/dev/null || \
  find "$WORK/out" -name putative_peptide.txt -exec wc -l {} \;
