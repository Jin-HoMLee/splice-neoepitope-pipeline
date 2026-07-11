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
#   - unknown-slug   : a defined arc:<slug> label absent from the manifest roster.
#
# ADVISORY NOTE (printed, never affects the exit code):
#   - cap            : more arcs marked `active` than ACTIVE_SLATE_GUIDELINE (default 3).
#                      Softened from a hard gate in Issue #1102: `active` was measured
#                      not to predict throughput, and a WIP limit belongs on work in
#                      flight, not on categories. Real WIP guards: check_ready_queue.sh.
#
# MANIFEST DEFECTS fail fast (exit 2, both modes): a phase outside {active,next,later}
# is a typo that would make every `gh issue edit --add-label arc-phase:<typo>` fail,
# so it is rejected up front rather than surfacing as a spurious "transient" retry.
#
# RESILIENCE: only the MUTATING call (gh issue edit) is retry-tolerant -- a transient
# 5xx on one issue logs + continues (exit 3 signals a retry) so a single failure can't
# abort the batch. READ calls (gh issue list / label list) fail loud via `set -e` by
# design: a truncated/failed read must not be mistaken for "clean".
#
# bash 3.2-compatible (macOS default): no associative arrays, no mapfile.
# Prereq: run scripts/pm/arc_labels.sh first (the labels must exist in the repo).
set -euo pipefail
REPO="Jin-HoMLee/splice-neoepitope-pipeline"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Advisory focus-slate guideline (Issue #1102). Not a gate: exceeding it prints a
# NOTE and never affects the exit code. Tunable for a deliberate wide slate.
# A non-numeric override would make `[[ n -gt $G ]]` emit a cryptic bash error and
# silently skip the NOTE (it does NOT abort -- an `if` condition is exempt from
# set -e). Fall back to the default rather than let an advisory knob misbehave.
ACTIVE_SLATE_GUIDELINE="${ACTIVE_SLATE_GUIDELINE:-3}"
if ! [[ "$ACTIVE_SLATE_GUIDELINE" =~ ^[0-9]+$ ]]; then
  echo "  WARN: ACTIVE_SLATE_GUIDELINE='$ACTIVE_SLATE_GUIDELINE' is not a non-negative integer; using 3" >&2
  ACTIVE_SLATE_GUIDELINE=3
fi

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
seen_multi=" "
seen_parent=" "

# --- manifest pass: roster + phase-vocab validation + advisory active-slate note ---
roster=" "
active_count=0
while read -r arc phase rest || [[ -n "${arc:-}" ]]; do
  [[ -z "${arc:-}" || "${arc:0:1}" == "#" ]] && continue
  case "$phase" in
    active|next|later) : ;;
    *) echo "  MANIFEST ERROR: arc $arc has invalid phase '$phase' (expected active|next|later)" >&2; exit 2 ;;
  esac
  roster="${roster}${arc} "
  [[ "$phase" == "active" ]] && active_count=$((active_count + 1))
done < "$MANIFEST"

if [[ $active_count -gt $ACTIVE_SLATE_GUIDELINE ]]; then
  echo "  NOTE cap: $active_count arcs marked 'active' (focus-slate guideline is ~$ACTIVE_SLATE_GUIDELINE; advisory)"
fi

# --- unknown-slug guard: any defined arc:* label not in the roster ---
while IFS= read -r lbl; do
  [[ -z "$lbl" ]] && continue
  case "$roster" in
    *" $lbl "*) : ;;
    *) echo "  DRIFT unknown-slug: label $lbl is defined but absent from the manifest roster"; drift=1 ;;
  esac
done < <(gh label list --repo "$REPO" --limit 200 --json name \
          --jq '.[].name | select(startswith("arc:"))')

# --- parent set (Issue #1103) ---
# Parents/epics are parked in the `Epic` Status and are never pulled, so an
# arc-phase focus marker is meaningless on them: they carry NONE. And because a
# parent legitimately spans themes (#1036 spans scoring-tcr-pmhc +
# immunogenicity-benchmark), multi-arc is LEGAL at the parent tier - the
# one-arc rule binds leaves only. Fetched once (one paginated query), not per
# issue, so the reconcile stays O(arcs) reads.
PARENTS=" "
while IFS= read -r pn; do
  [[ -n "${pn:-}" ]] && PARENTS="${PARENTS}${pn} "
done < <(gh api graphql --paginate -f query='
  query($endCursor: String) {
    repository(owner: "Jin-HoMLee", name: "splice-neoepitope-pipeline") {
      issues(states: OPEN, first: 100, after: $endCursor) {
        pageInfo { hasNextPage endCursor }
        nodes { number subIssuesSummary { total } }
      }
    }
  }' --jq '.data.repository.issues.nodes[] | select(.subIssuesSummary.total > 0) | .number')

# --- per-arc phase reconcile (membership discovered live via labels) ---
# One list call per arc returns each member's number + its arc/arc-phase labels
# inline (O(arcs) reads, not O(members)); --limit 1000 avoids gh's default 30-row
# truncation, which would silently hide drift on high-volume arcs.
while read -r arc phase rest || [[ -n "${arc:-}" ]]; do
  [[ -z "${arc:-}" || "${arc:0:1}" == "#" ]] && continue
  want_phase="arc-phase:$phase"
  while IFS=$'\t' read -r n lbls || [[ -n "${n:-}" ]]; do
    [[ -z "${n:-}" ]] && continue

    n_arcs=0
    cur_phases=""
    arcs_list=""
    for l in $lbls; do
      case "$l" in
        arc-phase:*) cur_phases="$cur_phases $l" ;;
        arc:*) n_arcs=$((n_arcs + 1)); arcs_list="$arcs_list $l" ;;
      esac
    done

    # --- parent tier (#1103): multi-arc is legal, and NO arc-phase is carried ---
    # A parent is parked in `Epic` Status and never pulled, so the focus marker is
    # meaningless on it. Multi-arc on a parent is therefore not ambiguous - it is
    # simply unphased. The one-arc rule binds LEAVES.
    if [[ "$PARENTS" == *" $n "* ]]; then
      if [[ -n "$cur_phases" ]]; then
        case "$seen_parent" in
          *" $n "*) : ;;  # already handled via another of its arcs
          *)
            drift=1
            seen_parent="${seen_parent}${n} "
            cur_csv="${cur_phases# }"; cur_csv="${cur_csv// /,}"
            if [[ $CHECK -eq 1 ]]; then
              echo "  DRIFT parent-phase: #$n has [$cur_csv], wants none (parents carry no arc-phase)"
            else
              remove_args=()
              for p in $cur_phases; do remove_args+=(--remove-label "$p"); done
              echo "  fix #$n -> strip arc-phase (parent; removing:$cur_phases)"
              if ! gh issue edit "$n" --repo "$REPO" "${remove_args[@]}"; then
                echo "  WARN: edit #$n failed (transient?); re-run to retry"
                failed=$((failed + 1))
              fi
            fi
            ;;
        esac
      fi
      continue
    fi

    # --- leaf tier: exactly one arc ---
    if [[ $n_arcs -gt 1 ]]; then
      drift=1
      case "$seen_multi" in
        *" $n "*) : ;;  # already reported for another arc
        *)
          arcs_csv="${arcs_list# }"; arcs_csv="${arcs_csv// /,}"
          echo "  DRIFT multi-arc: #$n carries $n_arcs arc labels ($arcs_csv) -> one arc per leaf (multi-arc is legal only on a parent, #1103)"
          seen_multi="${seen_multi}${n} "
          ;;
      esac
      continue
    fi

    stale=""
    have_want=0
    for p in $cur_phases; do
      if [[ "$p" == "$want_phase" ]]; then have_want=1; else stale="$stale $p"; fi
    done

    if [[ -n "$stale" || $have_want -eq 0 ]]; then
      drift=1
      if [[ $CHECK -eq 1 ]]; then
        cur_csv="${cur_phases# }"; cur_csv="${cur_csv// /,}"
        echo "  DRIFT phase: #$n ($arc) has [${cur_csv:-none}], wants $want_phase"
      else
        remove_args=()
        for p in $stale; do remove_args+=(--remove-label "$p"); done
        echo "  fix #$n -> $want_phase${stale:+ (removing:$stale)}"
        if ! gh issue edit "$n" --repo "$REPO" --add-label "$want_phase" \
              ${remove_args[@]+"${remove_args[@]}"}; then
          echo "  WARN: edit #$n failed (transient?); re-run to retry"
          failed=$((failed + 1))
        fi
      fi
    fi
  done < <(gh issue list --repo "$REPO" --label "$arc" --state open --limit 1000 \
             --json number,labels \
             --jq '.[] | "\(.number)\t" + ([.labels[].name | select(startswith("arc:") or startswith("arc-phase:"))] | join(" "))')
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
