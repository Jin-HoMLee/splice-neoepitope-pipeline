#!/usr/bin/env bash
# scripts/pm/arc_labels.sh
# Idempotently create the arc:* (narrative throughline) and arc-phase:* (focus
# marker) labels. Re-runnable: --force updates colour/description if the label
# already exists. Spec: docs/superpowers/specs/2026-06-05-arc-work-structuring-design.md
set -euo pipefail
REPO="Jin-HoMLee/splice-neoepitope-pipeline"

# name | hex colour (no '#') | description
ARC_LABELS=(
  "arc:aligner-junctions|0F766E|Arc: aligner & junction extraction (STAR/HISAT2 verification, SJ/annotation correctness)"
  "arc:junction-filtering|047857|Arc: junction filtering & tumor-specificity (GTEx + matched-normal + AlphaGenome, integrity)"
  "arc:variant-cohort|15803D|Arc: variant calling & cohort expansion (somatic calling, SpliceAI/MMSplice, new patients)"
  "arc:scoring-tcr-pmhc|B45309|Arc: neoepitope scoring & TCR-pMHC modeling (MHCflurry calibration, scorers, structures)"
  "arc:results-reporting|A16207|Arc: results, reporting & manuscript (HLA panel runs, report layer, decks, writeup)"
  "arc:cloud-reproducibility|1D4ED8|Arc: cloud execution & reproducibility (Batch, run registry, GPU, run_cloud_gpu.sh)"
  "arc:memory-methodology|7C3AED|Arc: memory & methodology (MEMORY.md slimming/audit, MM role, persona framing)"
  "arc:board-governance|9D174D|Arc: board governance & enforcement (Kanban governance, hooks/guards, PM-tooling evals)"
)
PHASE_LABELS=(
  "arc-phase:active|2DA44E|Arc focus: actively pulled now (cap 3)"
  "arc-phase:next|D4A72C|Arc focus: queued next"
  "arc-phase:later|8C959F|Arc focus: parked"
)

create() {
  local spec="$1" name color desc
  IFS='|' read -r name color desc <<<"$spec"
  gh label create "$name" --repo "$REPO" --color "$color" --description "$desc" --force
}
for spec in "${ARC_LABELS[@]}" "${PHASE_LABELS[@]}"; do create "$spec"; done
echo "Created/updated ${#ARC_LABELS[@]} arc + ${#PHASE_LABELS[@]} arc-phase labels on $REPO."
