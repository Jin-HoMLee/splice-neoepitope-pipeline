#!/usr/bin/env bash
# scripts/check_milestone_health.sh
#
# Surfaces open GitHub milestones that are overdue OR imminent (within
# THRESHOLD days of due_on). Intended for the PM morning-routine
# "Milestone health" phase.
#
# Why: per-Issue AC audits catch unticked boxes, but the milestone-level
# deadline is a separate axis that has no automated watch. Caught
# 2026-05-20 when `i1 - S4 - EDA - Junction Filtering Observability`
# went 5 days overdue without any session noticing.
#
# Usage:
#   bash scripts/check_milestone_health.sh [--threshold <days>] [--repo <owner/repo>]
#
# Env:
#   MILESTONE_HEALTH_THRESHOLD_DAYS  override default threshold (7)
#   REPO                             override default repo
#
# Exit codes:
#   0 — all clear OR only imminent (no overdue)
#   2 — at least one open milestone is overdue
#   1 — usage / runtime error

set -euo pipefail

REPO="${REPO:-Jin-HoMLee/splice-neoepitope-pipeline}"
THRESHOLD="${MILESTONE_HEALTH_THRESHOLD_DAYS:-7}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threshold) [[ -n "${2:-}" ]] || { echo "--threshold requires a value" >&2; exit 1; }; THRESHOLD="$2"; shift 2 ;;
        --repo) [[ -n "${2:-}" ]] || { echo "--repo requires a value" >&2; exit 1; }; REPO="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,22p' "$0" | sed 's|^# \{0,1\}||'
            exit 0
            ;;
        *) echo "unknown arg: $1" >&2; exit 1 ;;
    esac
done

[[ "$THRESHOLD" =~ ^[0-9]+$ ]] || { echo "threshold must be a non-negative integer, got: $THRESHOLD" >&2; exit 1; }

# NOTE: fetches first 100 open milestones only; GitHub's milestones endpoint
# paginates silently past per_page=100. If this repo ever exceeds 100 open
# milestones, swap to `gh api --paginate ... | jq --slurp 'add | ...'`.
FILTERED=$(gh api "repos/${REPO}/milestones?state=open&per_page=100" | jq --argjson threshold "$THRESHOLD" '
    [
        .[]
        | select(.due_on != null)
        | . + {days_until: (((.due_on | fromdateiso8601) - now) / 86400 | floor)}
        | select(.days_until <= $threshold)
    ]
    | sort_by(.days_until)
')

COUNT=$(echo "$FILTERED" | jq 'length')

if [[ "$COUNT" -eq 0 ]]; then
    echo "Milestone health: all clear — no overdue or imminent open milestones (threshold ${THRESHOLD}d)."
    exit 0
fi

printf "%-65s %-12s %-14s %5s %6s   %s\n" "TITLE" "DUE_ON" "DAYS" "OPEN" "CLOSED" "URL"
printf '%.0s-' {1..65}; printf ' '; printf '%.0s-' {1..12}; printf ' '; printf '%.0s-' {1..14}; printf ' '; printf '%.0s-' {1..5}; printf ' '; printf '%.0s-' {1..6}; printf '   '; printf '%.0s-' {1..40}; printf '\n'

echo "$FILTERED" | jq -r '
    .[] | [
        (.title | .[0:65]),
        (.due_on | sub("T.*$"; "")),
        (.days_until | tostring),
        (.open_issues | tostring),
        (.closed_issues | tostring),
        .html_url
    ] | @tsv
' | while IFS=$'\t' read -r title due days open closed url; do
    if [[ "$days" -lt 0 ]]; then
        days_display="${days#-}d overdue"
    elif [[ "$days" -eq 0 ]]; then
        days_display="due today"
    else
        days_display="${days}d left"
    fi
    printf "%-65s %-12s %-14s %5s %6s   %s\n" "$title" "$due" "$days_display" "$open" "$closed" "$url"
done

OVERDUE=$(echo "$FILTERED" | jq '[.[] | select(.days_until < 0)] | length')
if [[ "$OVERDUE" -gt 0 ]]; then
    exit 2
fi
exit 0
