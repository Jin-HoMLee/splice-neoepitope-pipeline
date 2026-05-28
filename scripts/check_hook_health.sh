#!/usr/bin/env bash
# Summarize fires logged to .claude/hook_fires.jsonl per hook.
# Operator-facing tool for hook-promotion decisions (Issue #453).
set -euo pipefail

LOG_PATH=".claude/hook_fires.jsonl"

# Helper function to get dock Issue for a hook.
# Duplicated from HOOK_CONFIG in .claude/hooks/recheck_dispatch.py — keep in sync (3 entries).
get_dock_issue() {
    case "$1" in
        target_sync_check)
            echo "Issue #454"
            ;;
        recheck_milestone)
            echo "(no dock Issue filed yet)"
            ;;
        recheck_parent_status)
            echo "(no dock Issue filed yet)"
            ;;
        *)
            echo "(unknown hook)"
            ;;
    esac
}

HOOKS=(target_sync_check recheck_milestone recheck_parent_status)

echo "Hook health summary (${LOG_PATH})"
echo ""

if [[ ! -f "${LOG_PATH}" ]]; then
    echo "No fires logged yet."
    echo ""
    echo "Wired hooks:"
    for hook in "${HOOKS[@]}"; do
        echo "  - ${hook}  (dock: $(get_dock_issue "$hook"))"
    done
    exit 0
fi

total=0
for hook in "${HOOKS[@]}"; do
    # Filter JSONL lines for this hook (one JSON object per line).
    # jq-safe wrapper: skip parse errors line-by-line, count valid matches.
    count=$({ jq -R 'fromjson? | select(.hook == $h) | .hook' --arg h "${hook}" "${LOG_PATH}" 2>/dev/null || true; } | wc -l)
    count=$((count))  # Convert to numeric (removes leading spaces)
    echo "${hook}"
    echo "  Fires:  ${count}"
    if [[ "${count}" -gt 0 ]]; then
        first=$({ jq -rR 'fromjson? | select(.hook == $h) | .ts' --arg h "${hook}" "${LOG_PATH}" 2>/dev/null || true; } | head -1 || echo "")
        first=${first:-}
        last=$({ jq -rR 'fromjson? | select(.hook == $h) | .ts' --arg h "${hook}" "${LOG_PATH}" 2>/dev/null || true; } | tail -1 || echo "")
        last=${last:-}
        echo "  First:  ${first}"
        echo "  Last:   ${last}"
    fi
    echo "  Dock:   $(get_dock_issue "$hook")"
    echo ""
    total=$((total + count))
done

echo "Total: ${total} fires across ${#HOOKS[@]} hooks"
