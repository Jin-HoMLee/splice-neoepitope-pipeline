#!/usr/bin/env bash
# scripts/check_ready_queue.sh
#
# Surfaces an unhealthy Ready queue using a two-part gate (Issue #754):
#   - PER-ROLE FLOOR (default 5) for PM / Scientist / Developer — a Kanban
#     "minimum order point" that prevents role-starvation. Memory Manager is
#     EXCLUDED (memory-curation pulls cross-repo / opportunistically, not from
#     a maintained board Ready buffer).
#   - TOTAL CAP (default 18 = floor x 3 roles + 3 headroom) — a WIP limit on
#     the commitment buffer so Ready doesn't over-deepen and inflate lead time.
#     The headroom above the summed floors (15) leaves a reachable "healthy"
#     band (total 15-17 = all floors met, under cap). A zero-headroom cap=15
#     could only ever read REPLENISH or CAP, never healthy (#754 review).
#
# The floor is a replenish *trigger*, not a force-commit quota: below floor
# for a role, pull that role's highest-priority DoR-ready Backlog; if fewer
# than the floor are DoR-ready, commit what's ready and flag a grooming gap —
# never force low-value work into Ready just to hit the number. The cap is a
# soft WIP limit: at/over cap, hold new commitments.
#
# Replaces the prior single count-only floor-of-3, which was blind to per-role
# distribution and priority mix (it skipped a warranted commitment on
# 2026-06-16 while PM had 0 / Dev had 1 pullable Ready item).
#
# Intended for the PM morning-routine Replenishment beat (commitment half).
#
# Usage:
#   bash scripts/check_ready_queue.sh [--floor <count>] [--cap <count>]
#
# Env:
#   READY_QUEUE_FLOOR      override per-role floor (default 5)
#   READY_QUEUE_CAP        override total cap (default 18)
#   READY_QUEUE_JSON_FILE  read the Ready-items JSON array from this file
#                          instead of calling board_open_items.py (test seam)
#
# Exit codes:
#   0 — healthy (every role at/above floor AND total below cap)
#   2 — needs attention (a role below floor, or total at/over cap)
#       Both conditions share exit 2 — distinguish via stdout:
#       `[REPLENISH <role>: K < F]` (commit more) vs `[CAP] …` (hold).
#   1 — usage / runtime error

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

# Source the Ready-items JSON: a fixture file (test seam) when set, else the
# live board via board_open_items.py (progress chatter on stderr is dropped;
# the JSON array lands on stdout).
if [[ -n "${READY_QUEUE_JSON_FILE:-}" ]]; then
    READY_JSON="$(cat "$READY_QUEUE_JSON_FILE")" || {
        echo "error: could not read READY_QUEUE_JSON_FILE=$READY_QUEUE_JSON_FILE" >&2
        exit 1
    }
else
    READY_JSON="$(python3 "$SCRIPT_DIR/board_open_items.py" --status Ready --json 2>/dev/null)" || {
        echo "error: board_open_items.py failed (gh auth / network?)" >&2
        exit 1
    }
fi

TOTAL="$(printf '%s' "$READY_JSON" | jq 'length' 2>/dev/null)" || {
    echo "error: could not parse Ready-queue JSON" >&2
    exit 1
}

needs_attention=0
breakdown=""

# Per-role floor check. An item counts toward every role label it carries
# (a role:pm + role:memory_manager item counts toward pm), so we test label
# membership rather than a single role field.
for role in "${FLOOR_ROLES[@]}"; do
    count="$(printf '%s' "$READY_JSON" | jq --arg r "role:$role" '[.[] | select(.labels | index($r))] | length')" || {
        echo "error: could not count Ready items for role:$role" >&2
        exit 1
    }
    breakdown+="${role}=${count} "
    if [[ "$count" -lt "$FLOOR" ]]; then
        echo "[REPLENISH ${role}: ${count} < ${FLOOR}] — pull ${role}'s highest-priority DoR-ready Backlog (see shared/feedback_board_hygiene.md)."
        needs_attention=1
    fi
done

# Total cap check (soft WIP limit on the commitment buffer).
if [[ "$TOTAL" -ge "$CAP" ]]; then
    echo "[CAP] Ready at ${TOTAL} (>= ${CAP}) — hold new commitments; the Ready buffer is at its WIP limit."
    needs_attention=1
fi

# Always print the full per-role breakdown so every run is self-documenting
# (not just the healthy path) — roles at/above floor are otherwise silent.
if [[ "$needs_attention" -eq 1 ]]; then
    echo "Ready by role: ${breakdown%% } — total ${TOTAL} (floor ${FLOOR}/role, cap ${CAP}) — needs attention."
    exit 2
fi

echo "Ready by role: ${breakdown%% } — total ${TOTAL} (floor ${FLOOR}/role, cap ${CAP}) — healthy."
exit 0
