#!/usr/bin/env bash
# scripts/audit_and_merge.sh
#
# Closure-ritual gate before `gh pr merge`. Refuses to merge if any `- [ ]`
# remains on the PR body Test plan OR on any Issue in closingIssuesReferences
# (under that Issue's Acceptance criteria heading).
#
# Why: the closure-ritual rule has broken 4× in 10 days despite being inlined
# into shared MEMORY.md. Memory of a declarative rule cannot reliably survive
# the action-distance from session start to merge. This script enforces at
# the exact moment of action.
#
# Usage:
#   bash scripts/audit_and_merge.sh <PR_NUMBER> [--squash|--merge|--rebase] [--delete-branch|--no-delete-branch]
#
# Examples:
#   bash scripts/audit_and_merge.sh 387
#   bash scripts/audit_and_merge.sh 387 --rebase --no-delete-branch
#
# Exit codes:
#   0 — merged successfully
#   1 — audit failed (unticked boxes); the unticked lines are printed
#   2 — usage error

set -euo pipefail

REPO="${REPO:-Jin-HoMLee/splice-neoepitope-pipeline}"

usage() {
    echo "Usage: $0 <PR_NUMBER> [--squash|--merge|--rebase] [--delete-branch|--no-delete-branch]" >&2
    exit 2
}

[[ $# -ge 1 ]] || usage
PR="$1"; shift
[[ "$PR" =~ ^[0-9]+$ ]] || usage

MERGE_TYPE="--squash"
DELETE_FLAG="--delete-branch"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --squash|--merge|--rebase) MERGE_TYPE="$1" ;;
        --delete-branch)            DELETE_FLAG="--delete-branch" ;;
        --no-delete-branch)         DELETE_FLAG="" ;;
        *) echo "Unknown flag: $1" >&2; usage ;;
    esac
    shift
done

# Extract unticked `- [ ]` lines under a specific `## <Heading>` markdown section.
# Stops at the next `## ` heading. Handles multiple occurrences of the heading.
unticked_under() {
    local heading="$1"
    awk -v h="^## ${heading}([[:space:]]|$)" '
        $0 ~ h    { in_section=1; next }
        /^## /    { in_section=0 }
        in_section && /^- \[ \]/ { print }
    '
}

count_boxes_under() {
    local heading="$1"
    awk -v h="^## ${heading}([[:space:]]|$)" '
        $0 ~ h    { in_section=1; next }
        /^## /    { in_section=0 }
        in_section && /^- \[[ xX]\]/ { c++ }
        END { print c+0 }
    '
}

FAILED=0

PR_BODY=$(gh pr view "$PR" --repo "$REPO" --json body --jq .body)
TEST_PLAN_GAPS=$(printf '%s\n' "$PR_BODY" | unticked_under "Test plan")
TEST_PLAN_TOTAL=$(printf '%s\n' "$PR_BODY" | count_boxes_under "Test plan")
if [[ -n "$TEST_PLAN_GAPS" ]]; then
    echo "✗ PR #${PR} Test plan has unticked boxes:" >&2
    printf '%s\n' "$TEST_PLAN_GAPS" | sed 's/^/    /' >&2
    FAILED=1
elif ! printf '%s\n' "$PR_BODY" | grep -q '^## Test plan'; then
    echo "⚠ PR #${PR} has no '## Test plan' section (audit skipped)." >&2
fi

LINKED_ISSUES=$(gh pr view "$PR" --repo "$REPO" --json closingIssuesReferences --jq '.closingIssuesReferences[].number')
LINKED_COUNT=0
AC_TOTAL=0
for ISSUE in $LINKED_ISSUES; do
    LINKED_COUNT=$((LINKED_COUNT + 1))
    ISSUE_BODY=$(gh issue view "$ISSUE" --repo "$REPO" --json body --jq .body)
    AC_GAPS=$(printf '%s\n' "$ISSUE_BODY" | unticked_under "Acceptance criteria")
    AC_TOTAL=$((AC_TOTAL + $(printf '%s\n' "$ISSUE_BODY" | count_boxes_under "Acceptance criteria")))
    if [[ -n "$AC_GAPS" ]]; then
        echo "✗ Issue #${ISSUE} Acceptance criteria has unticked boxes:" >&2
        printf '%s\n' "$AC_GAPS" | sed 's/^/    /' >&2
        FAILED=1
    fi
done

[[ "$FAILED" -eq 1 ]] && exit 1

MERGE_ARGS=("$MERGE_TYPE")
[[ -n "$DELETE_FLAG" ]] && MERGE_ARGS+=("$DELETE_FLAG")
gh pr merge "$PR" --repo "$REPO" "${MERGE_ARGS[@]}"

echo "✓ PR #${PR} merged (${TEST_PLAN_TOTAL} test-plan boxes ticked, ${AC_TOTAL} AC boxes across ${LINKED_COUNT} linked issues ticked)."
