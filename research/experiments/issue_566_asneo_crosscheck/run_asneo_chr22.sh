#!/usr/bin/env bash
# #566 ASNEO cross-check — turnkey chr22 PoC runner (option B, MHCflurry-swap).
#
# Produces ASNEO's junction-derived candidate peptides (putative_peptide.txt) from
# the same chr22 test reads our pipeline uses, on hg19, by RE-ALIGNING to hg19 with
# STAR (ASNEO consumes STAR SJ.out.tab natively — avoids cross-build coordinate liftover).
#
# VM-BOUND: needs STAR (not installed/usable on the arm64 laptop; project policy =
# STAR is VM-only). Run on the project VM (or any host with STAR). The downstream
# MHCflurry concordance vs our pipeline output is the experiment NOTEBOOK step, not here.
#
# Prereqs: conda env `asneo` (research/experiments/issue_566_asneo_crosscheck/asneo_env.yml)
#          + STAR on PATH. Run from the repo root.
set -euo pipefail

EXP="research/experiments/issue_566_asneo_crosscheck"
WORK="${EXP}/outputs/chr22_run"          # gitignored run area
REF="references/hg19/chr22.fa"
STAR_IDX="indices/star_hg19_chr22"
FASTQ="data/test/SRR9143066_test.fastq.gz"   # tumor (SRR9143066), single-end
ASNEO_DIR="${WORK}/ASNEO"
HLA_PLACEHOLDER="HLA-A02:01"             # option B: unused for candidate generation (MHC step bypassed)
THREADS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"

mkdir -p "$WORK" "$(dirname "$REF")" "$STAR_IDX"

# 1. hg19 chr22 reference (UCSC; verified 200 OK 2026-06-23) ----------------------
if [[ ! -s "$REF" ]]; then
  echo "[1/5] Fetching hg19 chr22 (UCSC)..."
  curl -L --progress-bar "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz" \
    | gunzip > "$REF"
else echo "[1/5] hg19 chr22 present — skip."; fi

# 2. STAR hg19 chr22 index (small genome -> reduced --genomeSAindexNbases) --------
if [[ ! -s "${STAR_IDX}/SAindex" ]]; then
  echo "[2/5] Building STAR hg19 chr22 index..."
  STAR --runMode genomeGenerate --genomeDir "$STAR_IDX" --genomeFastaFiles "$REF" \
       --genomeSAindexNbases 11 --runThreadN "$THREADS"
else echo "[2/5] STAR index present — skip."; fi

# 3. Align chr22 test reads -> SJ.out.tab ----------------------------------------
echo "[3/5] Aligning $FASTQ to hg19 chr22..."
STAR --runMode alignReads --genomeDir "$STAR_IDX" --readFilesIn "$FASTQ" \
     --readFilesCommand zcat --runThreadN "$THREADS" \
     --outSAMtype BAM Unsorted --outFileNamePrefix "${WORK}/chr22_"
SJ="${WORK}/chr22_SJ.out.tab"
test -s "$SJ" || { echo "ERROR: $SJ not produced"; exit 1; }

# 4. Clone + option-B patch ASNEO ------------------------------------------------
echo "[4/5] Cloning + patching ASNEO (option B)..."
[[ -d "$ASNEO_DIR" ]] || git clone -q --depth 1 https://github.com/bm2-lab/ASNEO.git "$ASNEO_DIR"
python "${EXP}/apply_optionB_patch.py" "${ASNEO_DIR}/ASNEO.py"
# NB: src/software.tar.gz (bundled NetMHCpan/NetCTLpan/pepmatch) is NOT extracted —
# option B uses conda bedtools and bypasses the MHC binaries entirely.

# 5. Run patched ASNEO -> candidate peptides -------------------------------------
echo "[5/5] Running patched ASNEO -> putative_peptide.txt..."
python "${ASNEO_DIR}/ASNEO.py" -j "$SJ" -a "$HLA_PLACEHOLDER" -g "$REF" -o "${WORK}/asneo_out"
echo "DONE -> ${WORK}/asneo_out/putative_peptide.txt"
echo "Next: notebook.ipynb — score these + our pipeline peptides via MHCflurry; compute concordance."
