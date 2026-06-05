#!/usr/bin/env bash
# scripts/pm/apply_arc_labels.sh
# Apply arc:* + arc-phase:* labels to issues per scripts/pm/arc_taxonomy.tsv.
# Idempotent: gh add-label of an already-present label is a no-op.
# Prereq: run scripts/pm/arc_labels.sh first (labels must exist).
set -euo pipefail
REPO="Jin-HoMLee/splice-neoepitope-pipeline"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${1:-$SCRIPT_DIR/arc_taxonomy.tsv}"
[[ -f "$MANIFEST" ]] || { echo "manifest not found: $MANIFEST" >&2; exit 1; }

while read -r arc phase rest; do
  [[ -z "${arc:-}" || "${arc:0:1}" == "#" ]] && continue
  for n in $rest; do
    echo "  #$n -> $arc + arc-phase:$phase"
    gh issue edit "$n" --repo "$REPO" \
      --add-label "$arc" --add-label "arc-phase:$phase"
  done
done < "$MANIFEST"
echo "Done applying arc taxonomy from $MANIFEST."
