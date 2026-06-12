#!/usr/bin/env bash
#
# star_flag_sweep.sh — STAR sensitivity-flag benchmark for Issue #411.
#
# Benchmarks the deferred Nature-paper STAR flags ONE AT A TIME against the
# current production baseline (workflow/rules/alignment.smk star_align), on the
# patient_002 tumor RNA-seq sample. Each variant is measured on:
#   - total SJ.out.tab junction count
#   - % annotated (SJ.out.tab col 6 == 1) vs novel (col 6 == 0)
#   - junctions with >=1 uniquely-mapping read (col 7 > 0)
#   - wall-clock runtime
#   - multi-mapper utilization (Log.final.out: % mapped to multiple loci,
#     and % mapped to too-many loci — the cap's discard bucket)
#
# Design notes:
#   * CPU-only. STAR here emits SJ.out.tab only (--outSAMtype None), no BAM,
#     no GPU. Run on the (warm) production VM with the P100 idle — Option B.
#   * The production STAR index (indices/star/, built WITH the GTF) PERSISTS on
#     the VM disk (it is not a Snakemake temp()); it is reused, not rebuilt.
#   * The patient_002 FASTQs ARE temp() in the pipeline and get cleaned up after
#     alignment — so they are re-downloaded once from the public B2 bucket and
#     reused across every variant.
#   * The --sjdbGTFfile-removal variant needs a SEPARATE no-GTF index. It builds
#     to ${NOGTF_INDEX_DIR} (default indices/star_nogtf/) so it does NOT clobber
#     the production indices/star/ cache (which is not auto-invalidated — a
#     documented gotcha).
#   * Idempotent: a variant whose SJ.out.tab already exists is skipped, so an
#     interrupted sweep resumes cheaply.
#   * --outSAMtype BAM ... (issue #411 last bullet) is intentionally NOT swept:
#     it is only needed if a downstream rule consumes a BAM, which the
#     junction-only pipeline does not. Out of scope for this sensitivity sweep.
#
# Usage:
#   bash scripts/star_flag_sweep.sh            # full sweep (2-pass, faithful)
#   TWOPASS=None bash scripts/star_flag_sweep.sh   # single-pass quick look (~half the time)
#   THREADS=16 bash scripts/star_flag_sweep.sh
#
# Run inside tmux on the VM; it is a multi-hour unattended job.

set -euo pipefail

# ── Config (env-overridable) ─────────────────────────────────────────────────
THREADS="${THREADS:-$(nproc)}"
INDEX_DIR="${INDEX_DIR:-indices/star}"
NOGTF_INDEX_DIR="${NOGTF_INDEX_DIR:-indices/star_nogtf}"
GENOME="${GENOME:-references/GRCh38.primary_assembly.genome.fa}"
GTF="${GTF:-references/gencode.v47.annotation.gtf.gz}"
OUTROOT="${OUTROOT:-results/star_flag_sweep}"
TWOPASS="${TWOPASS:-Basic}"   # Basic = production-faithful; None = single-pass quick look

# patient_002 tumor (BG003082_T0), paired-end, from config/samples/patient_002.tsv
FASTQ_DIR="${FASTQ_DIR:-data/patient_002/BG003082_T0}"
FASTQ1_URL="https://b2.osteosarc.com/rna-seq/fastq/bostongene_2022/202211_bostongene_tumor_rna_BG003082_R1.fastq.gz"
FASTQ2_URL="https://b2.osteosarc.com/rna-seq/fastq/bostongene_2022/202211_bostongene_tumor_rna_BG003082_R2.fastq.gz"
FASTQ1="${FASTQ_DIR}/$(basename "$FASTQ1_URL")"
FASTQ2="${FASTQ_DIR}/$(basename "$FASTQ2_URL")"

# Baseline align flags = the production star_align command (alignment.smk).
BASELINE_FLAGS=(
  --outSAMtype None
  --twopassMode "$TWOPASS"
  --limitSjdbInsertNsj 2000000
)

# Variants: label → extra flag(s) layered on top of the baseline, ONE knob each.
# (sjdbGTFfile-removal is handled separately — it is an index change, not an
#  align flag.)
VARIANT_LABELS=(
  matchNmin_0.33
  scoreMin_0.33
  multimap_20
  intronMax_500k
  matesGap_1M
  sjOverhang_1
)
declare -A VARIANT_FLAGS=(
  [matchNmin_0.33]="--outFilterMatchNminOverLread 0.33"
  [scoreMin_0.33]="--outFilterScoreMinOverLread 0.33"
  [multimap_20]="--outFilterMultimapNmax 20"
  [intronMax_500k]="--alignIntronMax 500000"
  [matesGap_1M]="--alignMatesGapMax 1000000"
  [sjOverhang_1]="--alignSJDBoverhangMin 1"
)

mkdir -p "$OUTROOT"
RESULTS_TSV="${OUTROOT}/sweep_metrics.tsv"
printf 'variant\tflag\ttotal_sj\tannotated\tnovel\tpct_annotated\tuniq_supported\tmultiloci_pct\ttoomany_pct\truntime_s\n' > "$RESULTS_TSV"

log()  { printf '\n\033[1;34m[sweep %(%H:%M:%S)T]\033[0m %s\n' -1 "$*"; }

# ── Step 0: FASTQs (download once, reuse) ────────────────────────────────────
# Atomic download to .part then mv, so an interrupted/partial transfer is never
# mistaken for a complete cached file by the -s check below. --http1.1 avoids
# B2's HTTP/2 "stream not closed cleanly" flakiness (curl 92); --retry-all-errors
# retries transport-level failures (curl 92 is not an HTTP status code, so a
# plain --retry would not catch it).
mkdir -p "$FASTQ_DIR"
for pair in "$FASTQ1_URL|$FASTQ1" "$FASTQ2_URL|$FASTQ2"; do
  url="${pair%%|*}"; dst="${pair##*|}"
  if [[ -s "$dst" ]]; then
    log "FASTQ cached: $dst"
  else
    log "Downloading $(basename "$dst") from B2 ..."
    curl --fail -L --http1.1 --retry 5 --retry-delay 10 --retry-all-errors \
         --no-progress-meter -o "${dst}.part" "$url"
    mv "${dst}.part" "$dst"
  fi
done

# ── Step 1: baseline index (reuse warm prod index; build only if absent) ─────
if [[ -f "${INDEX_DIR}/index.done" || -f "${INDEX_DIR}/SAindex" ]]; then
  log "Reusing existing STAR index: $INDEX_DIR"
else
  log "Building STAR index (with GTF) → $INDEX_DIR (this is the ~32 GB-RAM step)"
  mkdir -p "$INDEX_DIR"
  GTF_FILE="$GTF"
  if [[ "$GTF" == *.gz ]]; then
    GTF_FILE="$(mktemp)"; trap 'rm -f "$GTF_FILE"' EXIT
    gunzip -c "$GTF" > "$GTF_FILE"
  fi
  STAR --runMode genomeGenerate --runThreadN "$THREADS" \
       --genomeDir "$INDEX_DIR" --genomeFastaFiles "$GENOME" \
       --sjdbGTFfile "$GTF_FILE" --sjdbOverhang 100
  touch "${INDEX_DIR}/index.done"
fi

# ── Helper: run one STAR alignment + record metrics ──────────────────────────
run_variant() {
  local label="$1" index_dir="$2"; shift 2
  local extra=("$@")
  local outdir="${OUTROOT}/${label}"
  local prefix="${outdir}/star_"
  local sjtab="${prefix}SJ.out.tab"
  local logfinal="${prefix}Log.final.out"

  if [[ -s "$sjtab" ]]; then
    log "SKIP ${label} (SJ.out.tab already present)"
  else
    mkdir -p "$outdir"
    log "RUN  ${label}   extra: ${extra[*]:-<none>}"
    local t0 t1; t0=$(date +%s)
    # Tolerate a single variant's failure: log it and skip, but keep the sweep
    # going so one bad flag does not lose the whole unattended run. (STAR in an
    # if-condition is exempt from set -e.)
    if STAR --runMode alignReads --runThreadN "$THREADS" \
         --genomeDir "$index_dir" \
         --readFilesIn "$FASTQ1" "$FASTQ2" --readFilesCommand zcat \
         --outFileNamePrefix "$prefix" \
         "${BASELINE_FLAGS[@]}" "${extra[@]}"; then
      t1=$(date +%s); echo $((t1 - t0)) > "${outdir}/runtime_s.txt"
    else
      log "!! ${label} FAILED (STAR nonzero) — skipping this variant, continuing sweep"
      return 0
    fi
  fi

  # Metrics (guard: only emit a row if an SJ.out.tab actually exists)
  [[ -s "$sjtab" ]] || { log "!! ${label}: no SJ.out.tab — no metrics row"; return 0; }
  local secs total annot novel pct uniq multi toomany flagstr
  secs=$(cat "${outdir}/runtime_s.txt" 2>/dev/null || echo NA)
  total=$(wc -l < "$sjtab")
  annot=$(awk '$6==1' "$sjtab" | wc -l)
  novel=$(awk '$6==0' "$sjtab" | wc -l)
  pct=$(awk -v a="$annot" -v t="$total" 'BEGIN{printf "%.1f", t?100*a/t:0}')
  uniq=$(awk '$7>0' "$sjtab" | wc -l)
  multi=$(grep "% of reads mapped to multiple loci" "$logfinal" | awk '{print $NF}' | tr -d '%')
  toomany=$(grep "% of reads mapped to too many loci" "$logfinal" | awk '{print $NF}' | tr -d '%')
  flagstr="${extra[*]:-<baseline>}"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$label" "$flagstr" "$total" "$annot" "$novel" "$pct" "$uniq" "${multi:-NA}" "${toomany:-NA}" "$secs" \
    >> "$RESULTS_TSV"
}

# ── Step 2: baseline ─────────────────────────────────────────────────────────
# Reuse a config-identical existing production STAR run as the baseline column
# (set REUSE_BASELINE_DIR) rather than recomputing it. The existing run MUST be
# the same prod star_align config: same index, same FASTQs, same 2-pass flags,
# same STAR build. We stage its raw artifacts into the baseline outdir; run_variant
# then sees the SJ.out.tab as "already present" and only computes the metrics row.
if [[ -n "${REUSE_BASELINE_DIR:-}" ]]; then
  if [[ -s "${REUSE_BASELINE_DIR}/star_SJ.out.tab" && -s "${REUSE_BASELINE_DIR}/star_Log.final.out" ]]; then
    log "Reusing existing baseline from $REUSE_BASELINE_DIR (no recompute)"
    mkdir -p "${OUTROOT}/baseline"
    cp -f "${REUSE_BASELINE_DIR}/star_SJ.out.tab"    "${OUTROOT}/baseline/star_SJ.out.tab"
    cp -f "${REUSE_BASELINE_DIR}/star_Log.final.out" "${OUTROOT}/baseline/star_Log.final.out"
    echo "reused" > "${OUTROOT}/baseline/runtime_s.txt"
  else
    log "!! REUSE_BASELINE_DIR set but artifacts missing — falling back to a fresh baseline run"
  fi
fi
run_variant "baseline" "$INDEX_DIR"

# ── Step 3: one-flag-at-a-time variants (1–6) ────────────────────────────────
for label in "${VARIANT_LABELS[@]}"; do
  # shellcheck disable=SC2206
  extra=(${VARIANT_FLAGS[$label]})
  run_variant "$label" "$INDEX_DIR" "${extra[@]}"
done

# ── Step 4: --sjdbGTFfile removal (separate no-GTF index) ─────────────────────
if [[ -f "${NOGTF_INDEX_DIR}/index.done" || -f "${NOGTF_INDEX_DIR}/SAindex" ]]; then
  log "Reusing existing no-GTF index: $NOGTF_INDEX_DIR"
else
  log "Building no-GTF STAR index → $NOGTF_INDEX_DIR (separate from prod cache)"
  mkdir -p "$NOGTF_INDEX_DIR"
  # No --sjdbOverhang here: STAR rejects >0 overhang when generating a genome
  # without annotations (--sjdbGTFfile / --sjdbFileChrStartEnd). Overhang only
  # has meaning relative to an annotated splice-junction database, which this
  # index deliberately omits.
  if STAR --runMode genomeGenerate --runThreadN "$THREADS" \
       --genomeDir "$NOGTF_INDEX_DIR" --genomeFastaFiles "$GENOME"; then
    touch "${NOGTF_INDEX_DIR}/index.done"
  else
    log "!! no-GTF index build FAILED — skipping noGTF_index variant"
  fi
fi
# NOTE: with no sjdb in the index, EVERY junction is reported as novel (col 6 == 0)
# by construction — % annotated is meaningless for this variant; compare on total
# and unique-supported counts.
if [[ -f "${NOGTF_INDEX_DIR}/SAindex" ]]; then
  run_variant "noGTF_index" "$NOGTF_INDEX_DIR"
else
  log "!! skipping noGTF_index (no usable index)"
fi

# ── Step 5: summary table (markdown, ready to paste into developer.md) ────────
log "Sweep complete. Metrics → $RESULTS_TSV"
echo
column -t -s $'\t' "$RESULTS_TSV"
echo
{
  echo "| variant | total SJ | %annot | uniq-supported | multi-loci % | too-many % | runtime s |"
  echo "|---|---|---|---|---|---|---|"
  tail -n +2 "$RESULTS_TSV" | awk -F'\t' \
    '{printf "| %s | %s | %s | %s | %s | %s | %s |\n",$1,$3,$6,$7,$8,$9,$10}'
} | tee "${OUTROOT}/summary_table.md"
