#!/usr/bin/env bash
# scripts/pm/apply_arc_labels.sh
# Re-sync arc:* + arc-phase:* labels on issues to match scripts/pm/arc_taxonomy.tsv.
#
# True re-sync (Issue #689), not add-only: each manifest line lists one arc, one
# arc-phase, and the issues under it. For every listed issue this:
#   1. removes any arc:* / arc-phase:* label NOT matching the manifest pair, then
#   2. adds the manifest-correct arc + arc-phase.
# So after an arc-review slate change (a phase edit, a re-tag, a split/merge in
# arc_taxonomy.tsv) a re-run leaves each LISTED issue with exactly its manifest
# pair — the TSV is authoritative for every issue it lists. (Scope note: an issue
# silently *removed* from the manifest keeps its old arc labels; de-labeling a
# dropped issue would need a separate full `gh issue list --label arc:*` sweep,
# out of scope here.) Idempotent: a no-op manifest re-applies the same labels and
# removes nothing.
#
# Prereq: run scripts/pm/arc_labels.sh first (the labels must exist in the repo).
set -euo pipefail
REPO="Jin-HoMLee/splice-neoepitope-pipeline"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${1:-$SCRIPT_DIR/arc_taxonomy.tsv}"
[[ -f "$MANIFEST" ]] || { echo "manifest not found: $MANIFEST" >&2; exit 1; }

while read -r arc phase rest || [[ -n "${arc:-}" ]]; do
  [[ -z "${arc:-}" || "${arc:0:1}" == "#" ]] && continue
  want_phase="arc-phase:$phase"
  for n in $rest; do
    # Collect the issue's current arc:* / arc-phase:* labels, then flag for removal
    # any that aren't the manifest-correct pair. read-loop (not mapfile) + guarded
    # array expansion keep this bash 3.2-compatible (macOS default shell).
    remove_args=()
    while IFS= read -r lbl; do
      [[ -z "$lbl" ]] && continue
      [[ "$lbl" == "$arc" || "$lbl" == "$want_phase" ]] && continue
      remove_args+=(--remove-label "$lbl")
    done < <(gh issue view "$n" --repo "$REPO" --json labels \
              --jq '.labels[].name | select(startswith("arc:") or startswith("arc-phase:"))')

    if [[ ${#remove_args[@]} -gt 0 ]]; then
      echo "  #$n -> $arc + $want_phase (removing stale: ${remove_args[*]//--remove-label/})"
    else
      echo "  #$n -> $arc + $want_phase"
    fi
    gh issue edit "$n" --repo "$REPO" \
      --add-label "$arc" --add-label "$want_phase" \
      ${remove_args[@]+"${remove_args[@]}"}
  done
done < "$MANIFEST"
echo "Done re-syncing arc taxonomy from $MANIFEST."
