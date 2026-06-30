#!/usr/bin/env bash
# scripts/check_ready_queue.sh
#
# Surfaces an unhealthy Ready queue using a demand-aware per-role trigger plus
# a total cap.
#
#   - PER-ROLE DEMAND-AWARE FLOOR (default 3) for PM / Scientist / Developer.
#     A role's [REPLENISH] fires only when that role is actually *consuming* -
#     i.e. it has >= 1 item In progress - AND its Ready count is below the
#     floor. A role idle on a quiet/post-ship board (0 In progress) does NOT
#     nag, even at 0 Ready: with no demand drawing the buffer down, thin Ready
#     is the correct late-commitment posture (Last Responsible Moment). Memory
#     Manager is EXCLUDED (memory-curation pulls cross-repo / opportunistically,
#     not from a maintained board Ready buffer).
#   - TOTAL CAP (default 12 = floor x 3 roles + 3 headroom) - a WIP limit on
#     the commitment buffer so Ready doesn't over-deepen and inflate lead time.
#     The cap is demand-independent: an over-deep buffer is over-commitment
#     regardless of consumption.
#
# Why demand-aware (Issue #902 facet 1, web-grounded 2026-06-30): a Kanban
# replenishment "minimum order point" should be sized from consumption x lead
# time and replenished demand-pull ("stock arriving just before it runs out"),
# never as a guessed fixed number. The buffer exists only to keep the
# downstream (In progress) from starving - with no demand there is nothing to
# starve. We have no flow data yet, so the faithful translation is to gate the
# trigger on actual draw-down and keep the number to one WIP slate (3), then
# retune from flow data later. Replaces the #754 fixed floor-5/cap-18, which
# fired false REPLENISH pressure on quiet post-ship boards (2026-06-29 and
# 2026-06-30 incidents).
#
# The floor is a replenish *trigger*, not a force-commit quota: below floor for
# a consuming role, pull that role's highest-priority DoR-ready Backlog; if
# fewer than the floor are DoR-ready, commit what's ready and flag a grooming
# gap - never force low-value work into Ready just to hit the number. The cap is
# a soft WIP limit: at/over cap, hold new commitments.
#
# Intended for the PM morning-routine Replenishment beat (commitment half).
#
# Usage:
#   bash scripts/check_ready_queue.sh [--floor <count>] [--cap <count>]
#
# Env:
#   READY_QUEUE_FLOOR      override per-role floor (default 3)
#   READY_QUEUE_CAP        override total cap (default 12)
#   BOARD_ITEMS_JSON_FILE  read the open-items JSON array from this file instead
#                          of calling board_open_items.py (test seam). The array
#                          holds objects with at least `.status` and `.labels`
#                          across all statuses; the script filters Ready /
#                          In progress (the test fixtures also carry `.number`).
#
# Exit codes:
#   0 - healthy (no consuming role below floor AND total below cap)
#   2 - needs attention (a consuming role below floor, or total at/over cap)
#       Both conditions share exit 2 - distinguish via stdout:
#       `[REPLENISH <role>: K < F]` (commit more) vs `[CAP] …` (hold).
#   1 - usage / runtime error

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FLOOR="${READY_QUEUE_FLOOR:-3}"
CAP="${READY_QUEUE_CAP:-12}"

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
# the Ready buffer and the In-progress demand signal come from one snapshot.
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

# Split into Ready and In-progress arrays once.
READY_JSON="$(printf '%s' "$ITEMS_JSON" | jq '[.[] | select(.status == "Ready")]' 2>/dev/null)" || {
    echo "error: could not parse open-items JSON" >&2
    exit 1
}
INPROG_JSON="$(printf '%s' "$ITEMS_JSON" | jq '[.[] | select(.status == "In progress")]' 2>/dev/null)"

TOTAL="$(printf '%s' "$READY_JSON" | jq 'length')"

needs_attention=0
ready_breakdown=""
inprog_breakdown=""

# Per-role demand-aware floor check. An item counts toward every role label it
# carries (a role:pm + role:memory_manager item counts toward pm), so we test
# label membership rather than a single role field.
for role in "${FLOOR_ROLES[@]}"; do
    ready_count="$(printf '%s' "$READY_JSON" | jq --arg r "role:$role" '[.[] | select(.labels | index($r))] | length')"
    inprog_count="$(printf '%s' "$INPROG_JSON" | jq --arg r "role:$role" '[.[] | select(.labels | index($r))] | length')"
    ready_breakdown+="${role}=${ready_count} "
    inprog_breakdown+="${role}=${inprog_count} "
    if [[ "$ready_count" -lt "$FLOOR" ]]; then
        if [[ "$inprog_count" -ge 1 ]]; then
            echo "[REPLENISH ${role}: ${ready_count} < ${FLOOR}] - consuming (${inprog_count} In progress); pull ${role}'s highest-priority DoR-ready Backlog (see shared/feedback_board_hygiene.md)."
            needs_attention=1
        else
            # Below floor but no demand: quiet board. Holding thin Ready is the
            # correct late-commitment posture, not a replenish trigger.
            echo "[quiet ${role}: ${ready_count} Ready, 0 In progress] - no active demand; holding (Last Responsible Moment)."
        fi
    fi
done

# Total cap check (soft WIP limit on the commitment buffer, demand-independent).
if [[ "$TOTAL" -ge "$CAP" ]]; then
    echo "[CAP] Ready at ${TOTAL} (>= ${CAP}) - hold new commitments; the Ready buffer is at its WIP limit."
    needs_attention=1
fi

# Always print the full per-role breakdown (Ready buffer + In-progress demand)
# so every run is self-documenting, not just the attention path.
status_line="Ready by role: ${ready_breakdown%% } | In progress by role: ${inprog_breakdown%% } | total Ready ${TOTAL} (floor ${FLOOR}/role, cap ${CAP})"

if [[ "$needs_attention" -eq 1 ]]; then
    echo "${status_line} - needs attention."
    exit 2
fi

echo "${status_line} - healthy."
exit 0
