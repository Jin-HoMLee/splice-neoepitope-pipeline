#!/usr/bin/env bash
# scripts/check_ready_queue.sh
#
# Surfaces an under-stocked Ready queue using a proactive per-role floor plus a
# total cap, and routes a shortfall to the right remedy (commit vs groom).
#
#   - PER-ROLE PROACTIVE FLOOR (default 5) for PM / Scientist / Developer / MM.
#     The floor is a TARGET shelf depth: keep ~5 DoR-ready items committed per
#     role so that when a role sits down there is always a curated shortlist to
#     pull, without paying the commitment-decision tax mid-session. It is
#     proactive (kept stocked ahead of demand), NOT gated on current
#     consumption: the buffer's whole value is being stocked before demand
#     arrives. Memory Manager is held to the SAME floor as every role (Jin-Ho
#     2026-07-04, Issue #1006): the pre-#902 MM floor-exemption (#705) assumed a
#     fixed floor would nag a bursty/blocked MM lane into stuffing junk, but #902
#     facet-1 already solved stuffing for all roles - a thin lane reads
#     GROOMING-GAP (groom, don't stuff), so a genuinely-empty MM lane is honest,
#     not a violation. The web-canonical rule is to start every workstream at a
#     uniform limit and only scale by capacity AFTER flow data shows a mismatch;
#     the exemption was a pre-emptive scale with no such data, so it is retired.
#   - TOTAL CAP (default 23 = floor x 4 roles + 3 headroom) - a WIP limit on the
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
# NOT-PULLABLE READY ITEMS (Issue #1248 shipped the first cut; Issue #1294
# generalised it). The floor counts Ready *depth*, not *pullable* depth, so an
# Issue that cannot actually be worked ("build only if it slips a 2nd time", an
# unsettled decision fork, a date gate) must not count as healthy stock - Issue
# #841 was swept into Ready twice this way, Issue #929 sat there with unsettled
# decision-fork ACs. Each Ready item now carries a derived `not_pullable` reason
# from board_open_items.py, computed by the pullability predicate
# (scripts/pm/pullability.py) over natively-owned STRUCTURED sources - GitHub
# `blockedBy`, the `needs-design` / `trigger-gated` labels, and the `Start date`
# field - NOT the Issue body. Such items are BOTH excluded from the pullable
# floor count AND surfaced as [NOT-PULLABLE role], because either alone leaves a
# silent failure: exclude-only hides them, surface-only keeps the queue reading
# healthy.
#
# ONE DECISION HOME. The gate reason comes from one predicate consulting one
# authoritative source per cause; nothing is copied into a second store, so no
# state can drift. The prose scan (scripts/pm/not_pullable.py -> propose_label)
# is demoted to a low-recall PROPOSER that only SUGGESTS a label at the
# Backlog->Ready commitment; a human applies it, and the label - not the prose -
# is what the predicate reads. So a phrasing the proposer misses is a non-event
# once the label is set. RESIDUAL LIMIT: a gate not yet marked with a label or
# date (an un-backfilled Issue) is not seen until its marker is applied; the
# capability-prerequisite axis ("gated on sequence access", Issue #817) is a
# separate axis this predicate does not model.
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
#   bash scripts/check_ready_queue.sh [--floor <count>] [--cap <count>] [--review-limit <count>]
#
# Env:
#   READY_QUEUE_FLOOR      override per-role floor (default 5)
#   READY_QUEUE_CAP        override total cap (default 23)
#   REVIEW_WIP_LIMIT       override review-column WIP limit (default 10)
#   BOARD_ITEMS_JSON_FILE  read the open-items JSON array from this file instead
#                          of calling board_open_items.py (test seam). The array
#                          holds objects with at least `.status` and `.labels`
#                          across all statuses; the script filters Ready /
#                          In progress / Backlog / review-column PRs (the test
#                          fixtures carry `.number` + `.kind`).
#
# Exit codes:
#   0 - healthy (every role at/above floor, total below cap, review columns
#       below WIP limit)
#   2 - needs attention (a role below floor [REPLENISH or GROOMING-GAP], a
#       body-gated item in Ready [NOT-PULLABLE], total at/over cap [CAP], or
#       review columns at/over limit [REVIEW-DEBT]).
#       Distinguish via stdout tag: [REPLENISH] (commit DoR-ready work) vs
#       [GROOMING-GAP] (groom/intake) vs [NOT-PULLABLE] (decommit it to Backlog
#       or resolve the gate) vs [CAP] (hold new commitments) vs
#       [REVIEW-DEBT] (human review bandwidth is the binding constraint).
#   1 - usage / runtime error

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FLOOR="${READY_QUEUE_FLOOR:-5}"
CAP="${READY_QUEUE_CAP:-23}"
REVIEW_LIMIT="${REVIEW_WIP_LIMIT:-10}"

# Roles subject to the per-role floor. MM is held to the same floor as every
# role (Jin-Ho 2026-07-04, Issue #1006, retiring the #705/#754 MM exemption).
FLOOR_ROLES=(pm scientist developer memory_manager)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --floor)         [[ -n "${2:-}" ]] || { echo "--floor requires a value" >&2; exit 1; }; FLOOR="$2"; shift 2 ;;
        --cap)           [[ -n "${2:-}" ]] || { echo "--cap requires a value" >&2; exit 1; }; CAP="$2"; shift 2 ;;
        --review-limit)  [[ -n "${2:-}" ]] || { echo "--review-limit requires a value" >&2; exit 1; }; REVIEW_LIMIT="$2"; shift 2 ;;
        -h|--help)
            awk 'NR>1 && /^#/ {sub(/^# ?/, ""); print; next} NR>1 {exit}' "$0"
            exit 0
            ;;
        *) echo "unknown arg: $1" >&2; exit 1 ;;
    esac
done

[[ "$FLOOR" =~ ^[0-9]+$ ]]        || { echo "floor must be a non-negative integer, got: $FLOOR" >&2; exit 1; }
[[ "$CAP" =~ ^[0-9]+$ ]]          || { echo "cap must be a non-negative integer, got: $CAP" >&2; exit 1; }
[[ "$REVIEW_LIMIT" =~ ^[0-9]+$ ]] || { echo "review-limit must be a non-negative integer, got: $REVIEW_LIMIT" >&2; exit 1; }

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
# PR cards only — a PR and its linked Issue both sit in review columns
# (PR at Ready for review, Issue at In review), so counting cards would
# double each review unit. Filtering to kind=="PR" gives the real PRs-
# awaiting-review count.
REVIEW_JSON="$(printf '%s' "$ITEMS_JSON" | jq '[.[] | select((.status == "Ready for review" or .status == "In review") and .kind == "PR")]' 2>/dev/null)"

# Cap counts ALL Ready cards (gated included); the floor counts only pullable
# ones. This asymmetry is the documented Kanban WIP-limit convention - blocked
# items still consume a WIP slot, so they count against the cap, but they cannot
# satisfy the Ready floor because they are not workable. Not a local choice.
TOTAL="$(printf '%s' "$READY_JSON" | jq 'length')"

needs_attention=0
ready_breakdown=""
inprog_breakdown=""
backlog_breakdown=""

# Per-role proactive floor check. An item counts toward every role label it
# carries (a role:pm + role:memory_manager item counts toward pm), so we test
# label membership rather than a single role field.
for role in "${FLOOR_ROLES[@]}"; do
    # Split this role's Ready items into pullable vs body-gated (Issue #1248).
    # `not_pullable` is a derived reason string from board_open_items.py, or
    # null/absent when the item is workable. An absent key reads as pullable, so
    # older fixtures and non-Issue cards are unaffected.
    role_ready_json="$(printf '%s' "$READY_JSON" | jq --arg r "role:$role" '[.[] | select(.labels | index($r))]')"
    ready_count="$(printf '%s' "$role_ready_json" | jq '[.[] | select(.not_pullable == null)] | length')"
    gated_json="$(printf '%s' "$role_ready_json" | jq '[.[] | select(.not_pullable != null)]')"
    gated_count="$(printf '%s' "$gated_json" | jq 'length')"
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

    # Surface body-gated Ready items explicitly. Excluding them from the count
    # alone would fix the false-healthy read but make them INVISIBLE, which just
    # trades one silent failure for another; surfacing alone would leave the
    # queue still reading healthy, which is the original bug. Both halves are
    # needed, so Issue #1248's "exclude OR surface" is resolved as "exclude AND
    # surface". A gated item sitting in Ready is a miscommit, and the remedy
    # (decommit it to Backlog) is always honestly available, so this flags
    # attention without the un-satisfiable-nag pathology described above.
    if [[ "$gated_count" -ge 1 ]]; then
        gated_detail="$(printf '%s' "$gated_json" | jq -r '[.[] | "#\(.number) (\(.not_pullable))"] | join(", ")')"
        echo "[NOT-PULLABLE ${role}: ${gated_count}] - in Ready but body-gated, excluded from the pullable count: ${gated_detail}. Decommit to Backlog, or resolve the gate (see shared/feedback_board_hygiene.md)."
        needs_attention=1
    fi
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

# Review-column WIP limit (advisory, keyed to human review bandwidth).
# Counts PRs (not cards — a PR and its linked Issue both sit in review
# columns, so card-count double-counts). Default 10 (top of the 5-10
# PRs/reviewer band from governance research;
# agent_team_governance_research_2026-07.md §7); tunable via --review-limit
# or REVIEW_WIP_LIMIT env. Advisory not blocking (house style).
REVIEW_TOTAL="$(printf '%s' "$REVIEW_JSON" | jq 'length')"
if [[ "$REVIEW_TOTAL" -ge "$REVIEW_LIMIT" ]]; then
    echo "[REVIEW-DEBT] PRs awaiting review at ${REVIEW_TOTAL} (>= ${REVIEW_LIMIT}) - human review is the binding throughput constraint; merge some PRs or slow down dispatch."
    needs_attention=1
fi

# Always print the full per-role breakdown (Ready buffer + In-progress demand
# context + Backlog candidate pool + review column load) so every run is
# self-documenting.
status_line="Ready by role: ${ready_breakdown%% } | In progress: ${inprog_breakdown%% } | Backlog candidates: ${backlog_breakdown%% } | Review PRs: ${REVIEW_TOTAL} | total Ready ${TOTAL} (floor ${FLOOR}/role, cap ${CAP}, review limit ${REVIEW_LIMIT})"

if [[ "$needs_attention" -eq 1 ]]; then
    echo "${status_line} - needs attention."
    exit 2
fi

echo "${status_line} - healthy."
exit 0
