#!/usr/bin/env bash
# scripts/check_ready_queue.sh
#
# Surfaces an under-stocked Ready queue using a proactive per-role floor plus a
# total cap, and routes a shortfall to the right remedy (commit vs groom).
#
#   - PER-ROLE PROACTIVE FLOOR (default 5) for PM / Scientist / Developer.
#     The floor is a TARGET shelf depth: keep ~5 DoR-ready items committed per
#     role so that when a role sits down there is always a curated shortlist to
#     pull, without paying the commitment-decision tax mid-session. It is
#     proactive (kept stocked ahead of demand), NOT gated on current
#     consumption: the buffer's whole value is being stocked before demand
#     arrives. Memory Manager is EXCLUDED (memory-curation pulls cross-repo /
#     opportunistically, not from a maintained board Ready buffer).
#   - TOTAL CAP (default 18 = floor x 3 roles + 3 headroom) - a WIP limit on the
#     commitment buffer so Ready doesn't over-deepen and inflate lead time.
#
# Shortfall routing (Issue #902 facet 1, refined 2026-06-30). A role below floor
# is NOT always a "commit now" alarm - it depends on whether committable work
# exists:
#   - [REPLENISH role] : short AND the role has triaged Backlog candidates ->
#     commit the highest-priority DoR-ready ones. If none actually meet the
#     Definition of Ready, that is a grooming gap - refine/file, do NOT stuff
#     low-value work into Ready just to hit the number. Arc-aware (Issue #931):
#     the nudge flags how many candidates are on an active arc (arc-phase:active)
#     and says to prefer those, so re-selection follows the active slate without
#     relying on memory (the arc axis lives in shared/feedback_arc_review.md).
#   - [GROOMING-GAP role] : short AND the role has NO Backlog candidates at all
#     -> nothing is committable; the remedy is intake/grooming, not commitment.
#
# Why this shape: the original #754 fixed floor-5 fired a chronic "you're short"
# nag that could not always be honestly satisfied (post-ship lulls, an
# un-groomed backlog), which created pressure to stuff junk to silence it. The
# fix is not a smaller buffer (5 is a good shelf from the consuming side) nor a
# consumption gate (that deletes the stocked-shelf benefit) - it is to make the
# shortfall interpretable and never demand a commitment that isn't there.
# HONEST LIMIT: DoR-readiness ("scope is clear") is human judgment, not a field,
# so the script cannot tell an un-groomed backlog from a ready one. It surfaces
# Backlog candidate counts as context; the human applies DoR on the REPLENISH
# branch. In-progress counts are reported as demand context only (not gating).
#
# Intended for the PM morning-routine Replenishment beat (commitment half).
#
# Usage:
#   bash scripts/check_ready_queue.sh [--floor <count>] [--cap <count>]
#
# Env:
#   READY_QUEUE_FLOOR      override per-role floor (default 5)
#   READY_QUEUE_CAP        override total cap (default 18)
#   BOARD_ITEMS_JSON_FILE  read the open-items JSON array from this file instead
#                          of calling board_open_items.py (test seam). The array
#                          holds objects with at least `.status` and `.labels`
#                          across all statuses; the script filters Ready /
#                          In progress / Backlog (the test fixtures also carry
#                          `.number`).
#
# Exit codes:
#   0 - healthy (every role at/above floor AND total below cap)
#   2 - needs attention (a role below floor [REPLENISH or GROOMING-GAP], or
#       total at/over cap). Distinguish via stdout: [REPLENISH] (commit DoR-ready
#       work) vs [GROOMING-GAP] (groom/intake) vs [CAP] (hold new commitments).
#   1 - usage / runtime error

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FLOOR="${READY_QUEUE_FLOOR:-5}"
CAP="${READY_QUEUE_CAP:-18}"

# Roles subject to the per-role floor. MM is intentionally excluded (#754).
FLOOR_ROLES=(pm scientist developer)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --floor) [[ -n "${2:-}" ]] || { echo "--floor requires a value" >&2; exit 1; }; FLOOR="$2"; shift 2 ;;
        --cap)   [[ -n "${2:-}" ]] || { echo "--cap requires a value" >&2; exit 1; }; CAP="$2"; shift 2 ;;
        -h|--help)
            awk 'NR>1 && /^#/ {sub(/^# ?/, ""); print; next} NR>1 {exit}' "$0"
            exit 0
            ;;
        *) echo "unknown arg: $1" >&2; exit 1 ;;
    esac
done

[[ "$FLOOR" =~ ^[0-9]+$ ]] || { echo "floor must be a non-negative integer, got: $FLOOR" >&2; exit 1; }
[[ "$CAP" =~ ^[0-9]+$ ]]   || { echo "cap must be a non-negative integer, got: $CAP" >&2; exit 1; }

# Source the open-items JSON: a fixture file (test seam) when set, else the live
# board via board_open_items.py (progress chatter on stderr is dropped; the JSON
# array lands on stdout). We fetch ALL open items and filter by status here, so
# the Ready buffer, the In-progress demand context, and the Backlog candidate
# pool all come from one snapshot.
if [[ -n "${BOARD_ITEMS_JSON_FILE:-}" ]]; then
    ITEMS_JSON="$(cat "$BOARD_ITEMS_JSON_FILE")" || {
        echo "error: could not read BOARD_ITEMS_JSON_FILE=$BOARD_ITEMS_JSON_FILE" >&2
        exit 1
    }
else
    ITEMS_JSON="$(python3 "$SCRIPT_DIR/board_open_items.py" --json 2>/dev/null)" || {
        echo "error: board_open_items.py failed (gh auth / network?)" >&2
        exit 1
    }
fi

# Split into Ready / In-progress / Backlog arrays once.
READY_JSON="$(printf '%s' "$ITEMS_JSON" | jq '[.[] | select(.status == "Ready")]' 2>/dev/null)" || {
    echo "error: could not parse open-items JSON" >&2
    exit 1
}
INPROG_JSON="$(printf '%s' "$ITEMS_JSON" | jq '[.[] | select(.status == "In progress")]' 2>/dev/null)"
BACKLOG_JSON="$(printf '%s' "$ITEMS_JSON" | jq '[.[] | select(.status == "Backlog")]' 2>/dev/null)"

TOTAL="$(printf '%s' "$READY_JSON" | jq 'length')"

needs_attention=0
ready_breakdown=""
inprog_breakdown=""
backlog_breakdown=""

# Per-role proactive floor check. An item counts toward every role label it
# carries (a role:pm + role:memory_manager item counts toward pm), so we test
# label membership rather than a single role field.
for role in "${FLOOR_ROLES[@]}"; do
    ready_count="$(printf '%s' "$READY_JSON" | jq --arg r "role:$role" '[.[] | select(.labels | index($r))] | length')"
    inprog_count="$(printf '%s' "$INPROG_JSON" | jq --arg r "role:$role" '[.[] | select(.labels | index($r))] | length')"
    backlog_count="$(printf '%s' "$BACKLOG_JSON" | jq --arg r "role:$role" '[.[] | select(.labels | index($r))] | length')"
    # Active-arc Backlog candidates for this role (carry arc-phase:active). The
    # daily gate is arc-aware: on a REPLENISH shortfall we flag active-arc
    # candidates so re-selection follows the active slate without relying on
    # someone remembering the lens at replenishment (Issue #931). The arc axis /
    # active slate is defined in shared/feedback_arc_review.md.
    active_count="$(printf '%s' "$BACKLOG_JSON" | jq --arg r "role:$role" '[.[] | select(.labels | index($r)) | select(.labels | index("arc-phase:active"))] | length')"
    ready_breakdown+="${role}=${ready_count} "
    inprog_breakdown+="${role}=${inprog_count} "
    backlog_breakdown+="${role}=${backlog_count} "
    if [[ "$ready_count" -lt "$FLOOR" ]]; then
        if [[ "$backlog_count" -ge 1 ]]; then
            if [[ "$active_count" -ge 1 ]]; then
                arc_hint=", ${active_count} on an active arc (arc-phase:active) - prefer these"
            else
                arc_hint=", none on an active arc - commit the top DoR-ready pool item, or flag whether a next-phase arc should be promoted"
            fi
            echo "[REPLENISH ${role}: ${ready_count} < ${FLOOR}] - ${backlog_count} Backlog candidate(s)${arc_hint}; commit the highest-priority DoR-ready ones. If none meet DoR, it's a grooming gap - refine/file, do not stuff (see shared/feedback_board_hygiene.md)."
        else
            echo "[GROOMING-GAP ${role}: ${ready_count} < ${FLOOR}] - 0 Backlog candidates; nothing committable for this role. File/intake, do not stuff."
        fi
        needs_attention=1
    fi
done

# Total cap check (soft WIP limit on the commitment buffer).
if [[ "$TOTAL" -ge "$CAP" ]]; then
    echo "[CAP] Ready at ${TOTAL} (>= ${CAP}) - hold new commitments; the Ready buffer is at its WIP limit."
    needs_attention=1
fi

# Always print the full per-role breakdown (Ready buffer + In-progress demand
# context + Backlog candidate pool) so every run is self-documenting.
status_line="Ready by role: ${ready_breakdown%% } | In progress: ${inprog_breakdown%% } | Backlog candidates: ${backlog_breakdown%% } | total Ready ${TOTAL} (floor ${FLOOR}/role, cap ${CAP})"

if [[ "$needs_attention" -eq 1 ]]; then
    echo "${status_line} - needs attention."
    exit 2
fi

echo "${status_line} - healthy."
exit 0
