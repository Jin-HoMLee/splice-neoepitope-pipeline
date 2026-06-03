#!/usr/bin/env bash
# scripts/check_ready_queue.sh
#
# Surfaces a starved Ready queue: counts open issues at board Status
# `Ready` (the committed pull-queue) and emits a [REPLENISH] nudge when
# the count falls below THRESHOLD. Intended for the PM morning-routine
# "Ready-queue replenishment" phase (Phase 2.8).
#
# Why: under late-commitment Kanban the failure mode shifts from a bloated
# Backlog to a starved Ready queue — too few committed, refined items for
# Dev/Sci to pull. Phase 2.6 (check_milestone_health.sh) watches milestone
# deadlines; this is its mirror, watching committed pull-queue depth. The
# nudge pulls a Backlog -> Ready commitment sweep forward off the bi-weekly
# cadence (see shared/feedback_board_hygiene.md).
#
# Usage:
#   bash scripts/check_ready_queue.sh [--threshold <count>]
#
# Env:
#   READY_QUEUE_THRESHOLD  override default threshold (3)
#
# Exit codes:
#   0 — Ready queue healthy (count >= threshold)
#   2 — Ready queue below threshold (replenishment nudge)
#   1 — usage / runtime error

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THRESHOLD="${READY_QUEUE_THRESHOLD:-3}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threshold) [[ -n "${2:-}" ]] || { echo "--threshold requires a value" >&2; exit 1; }; THRESHOLD="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,25p' "$0" | sed 's|^# \{0,1\}||'
            exit 0
            ;;
        *) echo "unknown arg: $1" >&2; exit 1 ;;
    esac
done

[[ "$THRESHOLD" =~ ^[0-9]+$ ]] || { echo "threshold must be a non-negative integer, got: $THRESHOLD" >&2; exit 1; }

# board_open_items.py paginates the ProjectV2 board (#9) and filters to open
# items at Status=Ready; its progress chatter goes to stderr (dropped here),
# the JSON array to stdout.
READY_JSON="$(python3 "$SCRIPT_DIR/board_open_items.py" --status Ready --json 2>/dev/null)" || {
    echo "error: board_open_items.py failed (gh auth / network?)" >&2
    exit 1
}

COUNT="$(printf '%s' "$READY_JSON" | jq 'length' 2>/dev/null)" || {
    echo "error: could not parse Ready-queue JSON" >&2
    exit 1
}

if [[ "$COUNT" -lt "$THRESHOLD" ]]; then
    echo "[REPLENISH] Ready queue at ${COUNT} (< threshold ${THRESHOLD}) — pull a Backlog → Ready commitment sweep (see shared/feedback_board_hygiene.md)."
    exit 2
fi

echo "Ready queue: healthy — ${COUNT} committed items ready to pull (threshold ${THRESHOLD})."
exit 0
