#!/usr/bin/env bash
# scripts/pm/apply_arc_labels.sh
# Reconcile arc-phase:* labels to the arc taxonomy (Issue #973, model A').
#
# SOURCE-OF-TRUTH MODEL (A', Issue #973):
#   - MEMBERSHIP (issue -> arc) is LABEL-authoritative: an issue's arc:<slug>
#     label (set at triage) is the truth. This script NEVER adds or removes an
#     arc:<slug> label.
#   - The manifest (arc_taxonomy.tsv) is authoritative ONLY for the arc ROSTER and
#     each arc's slate PHASE. One line per arc: "<arc_slug> <phase> [description]".
#     It carries NO member lists -- those drifted by construction (the bug #973 fixed).
#   - arc-phase:<phase> labels are DERIVED here: for each manifest arc, every OPEN
#     issue carrying arc:<slug> gets its arc-phase set to the arc's manifest phase,
#     and any stale arc-phase:* is removed. Never hand-set an arc-phase label.
#
# MODES:
#   (default)  apply     -- reconcile arc-phase labels on the live board.
#   --check    read-only -- report drift and exit 2 if any, else exit 0. Wire into
#              the arc-review ritual / a CI guard to catch drift early.
#
# DRIFT surfaced (both modes print it; --check exits 2 on any, mutates nothing):
#   - phase mismatch : an open arc member whose arc-phase != its arc's manifest phase.
#   - multi-arc      : >1 arc:<slug> label -> ambiguous phase; skipped, never auto-set.
#   - unknown-slug   : a defined arc:<slug> label absent from the manifest roster
#                      (a retired arc, or a new arc not yet added to the manifest).
#
# bash 3.2-compatible (macOS default): no associative arrays, no mapfile.
# Prereq: run scripts/pm/arc_labels.sh first (the labels must exist in the repo).
set -euo pipefail
REPO="Jin-HoMLee/splice-neoepitope-pipeline"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CHECK=0
MANIFEST=""
for a in "$@"; do
  case "$a" in
    --check) CHECK=1 ;;
    -*) echo "unknown flag: $a" >&2; exit 2 ;;
    *) MANIFEST="$a" ;;
  esac
done
MANIFEST="${MANIFEST:-$SCRIPT_DIR/arc_taxonomy.tsv}"
[[ -f "$MANIFEST" ]] || { echo "manifest not found: $MANIFEST" >&2; exit 1; }

mode="apply"; [[ $CHECK -eq 1 ]] && mode="check (read-only)"
echo "arc taxonomy reconcile [$mode] <- $MANIFEST"
drift=0
failed=0

# --- roster of valid arc slugs (for the unknown-slug guard) ---
roster=" "
while read -r arc phase rest || [[ -n "${arc:-}" ]]; do
  [[ -z "${arc:-}" || "${arc:0:1}" == "#" ]] && continue
  roster="${roster}${arc} "
done < "$MANIFEST"

# --- unknown-slug guard: any defined arc:* label not in the roster ---
while IFS= read -r lbl; do
  [[ -z "$lbl" ]] && continue
  case "$roster" in
    *" $lbl "*) : ;;
    *) echo "  DRIFT unknown-slug: label $lbl is defined but absent from the manifest roster"; drift=1 ;;
  esac
done < <(gh label list --repo "$REPO" --limit 200 --json name \
          --jq '.[].name | select(startswith("arc:"))')

# --- per-arc phase reconcile (membership discovered live via labels) ---
while read -r arc phase rest || [[ -n "${arc:-}" ]]; do
  [[ -z "${arc:-}" || "${arc:0:1}" == "#" ]] && continue
  want_phase="arc-phase:$phase"
  members=$(gh issue list --repo "$REPO" --label "$arc" --state open \
              --json number --jq '.[].number')
  for n in $members; do
    labels=$(gh issue view "$n" --repo "$REPO" --json labels \
              --jq '.labels[].name | select(startswith("arc:") or startswith("arc-phase:"))')

    n_arcs=$(printf '%s\n' "$labels" | grep -c '^arc:' || true)
    if [[ "${n_arcs:-0}" -gt 1 ]]; then
      arcs_csv=$(printf '%s\n' "$labels" | grep '^arc:' | paste -sd, -)
      echo "  DRIFT multi-arc: #$n carries $n_arcs arc labels ($arcs_csv) -> ambiguous phase, skipping (resolve at arc review)"
      drift=1
      continue
    fi

    cur_phases=$(printf '%s\n' "$labels" | grep '^arc-phase:' || true)
    stale=""
    have_want=0
    while IFS= read -r p; do
      [[ -z "$p" ]] && continue
      if [[ "$p" == "$want_phase" ]]; then have_want=1; else stale="$stale $p"; fi
    done <<EOF
$cur_phases
EOF

    if [[ -n "$stale" || $have_want -eq 0 ]]; then
      drift=1
      cur_csv=$(printf '%s\n' "$cur_phases" | grep . | paste -sd, - || true)
      if [[ $CHECK -eq 1 ]]; then
        echo "  DRIFT phase: #$n ($arc) has [${cur_csv:-none}], wants $want_phase"
      else
        remove_args=()
        for p in $stale; do remove_args+=(--remove-label "$p"); done
        echo "  fix #$n -> $want_phase${stale:+ (removing:$stale)}"
        # Tolerate a transient gh/network failure on one issue: log it and keep
        # going so a single 5xx can't abort the whole batch. The reconcile is
        # idempotent, so a re-run picks up any that failed (tracked in `failed`).
        if ! gh issue edit "$n" --repo "$REPO" --add-label "$want_phase" \
              ${remove_args[@]+"${remove_args[@]}"}; then
          echo "  WARN: edit #$n failed (transient?); re-run to retry"
          failed=$((failed + 1))
        fi
      fi
    fi
  done
done < "$MANIFEST"

if [[ $CHECK -eq 1 ]]; then
  [[ $drift -ne 0 ]] && { echo "arc taxonomy: DRIFT found (see above)."; exit 2; }
  echo "arc taxonomy: clean - arc-phase labels agree with the manifest."
  exit 0
fi
if [[ $failed -ne 0 ]]; then
  echo "Done with $failed transient edit failure(s) - re-run to retry (idempotent)." >&2
  exit 3
fi
echo "Done reconciling arc-phase labels from $MANIFEST."
