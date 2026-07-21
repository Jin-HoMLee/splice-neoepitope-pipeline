#!/usr/bin/env bash
# ASNEO runner (#966 AC-1): run the option-B-patched ASNEO caller on a given
# chr22-scale SJ.out.tab, no STAR, against the hg19 chr22 FASTA. Reuses the
# proven #965 asneo_smoke.sh sequence but parametrized on an input junction
# table, and normalizes the output to a deterministic path.
# Usage: bash run_asneo.sh <sj_out_tab> <workdir>
set -euo pipefail
SJ_TAB="${1:?sj_out_tab required}"
WORK="${2:?workdir required}"
REPO="$(git rev-parse --show-toplevel)"
PATCH="$REPO/research/experiments/issue_566_asneo_crosscheck/apply_optionB_patch.py"
mkdir -p "$WORK" "$WORK/out"

# 1. ASNEO clone (out of the repo tree).
[ -d "$WORK/ASNEO" ] || git clone --depth 1 https://github.com/bm2-lab/ASNEO.git "$WORK/ASNEO"

# 2. hg19 chr22 FASTA (UCSC) - small single-chr genome, fits 8 GB.
if [ ! -f "$WORK/chr22.fa" ]; then
  curl -sL "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz" -o "$WORK/chr22.fa.gz"
  gunzip -f "$WORK/chr22.fa.gz"
fi

# 3. Option-B patch (strip/bypass bundled netMHCpan; stop at candidate peptides).
#    The patcher is fail-loud non-idempotent, so skip it on a cached clone.
grep -q "option B" "$WORK/ASNEO/ASNEO.py" || conda run -n asneo python "$PATCH" "$WORK/ASNEO/ASNEO.py"

# 4. Run the patched caller at relaxed chr22 thresholds.
cd "$WORK/ASNEO"
conda run -n asneo python ASNEO.py -j "$SJ_TAB" -a HLA-A02:01 \
  -l 8,9,10,11 --reads 2 --psi 0.05 -g "$WORK/chr22.fa" -o "$WORK/out"

# 5. Normalize the output location: ASNEO writes putative_peptide.txt into a
#    sample subdir of -o. Find it and copy to a deterministic top-level path so
#    the Python side never has to guess (guards the path-discovery no-op).
SRC="$(find "$WORK/out" -name putative_peptide.txt | head -1)"
[ -n "$SRC" ] || { echo "ERROR: no putative_peptide.txt produced under $WORK/out" >&2; exit 1; }
DST="$WORK/out/putative_peptide.txt"
[ "$SRC" = "$DST" ] || cp "$SRC" "$DST"
echo "ASNEO candidate peptides: $(grep -cvE '^\s*(#|$)' "$DST") -> $DST"
