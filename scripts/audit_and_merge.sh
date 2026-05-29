#!/usr/bin/env bash
# scripts/audit_and_merge.sh
#
# Closure-ritual gate before `gh pr merge`. Refuses to merge if any of:
#   1. `- [ ]` remains on the PR body Test plan, OR
#   2. `- [ ]` remains under any linked Issue's Acceptance criteria, OR
#   3. Any linked Issue is missing a `**Priority rationale:**` line, OR
#   4. The would-be squash body (PR title + body + commit messages) contains a
#      closing keyword (`close|fix|resolve` + `#N`) targeting an Issue OUTSIDE
#      closingIssuesReferences — a silent unintended auto-close on merge.
#
# Why: declarative rules of the form "always include X" reliably drift across
# session boundaries even when inlined into shared MEMORY.md. The gate enforces
# at the moment of action (`gh pr merge`), so the rule cannot be forgotten
# regardless of how much context has accumulated since session start.
#
# Closure-ritual gate (1+2) shipped via Issue #357 after 4× drift in 10 days.
# Priority-rationale gate (3) added via Issue #481 after live drift in the
# Issue #264 closure flow. Stray-closing-keyword gate (4) added via Issue #559
# after PR #543 auto-closed epic Issue #538 via a commit-body keyword that the
# closingIssuesReferences API does not surface.
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
#   1 — audit failed (unticked boxes or missing priority rationale); the gaps are printed
#   2 — usage error

set -euo pipefail

REPO="${REPO:-Jin-HoMLee/splice-neoepitope-pipeline}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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
PR_RATIONALE_OK=0
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
    # Priority-rationale gate: case-insensitive substring grep. A line like
    # "**Priority rationale:** (deferred to #X)" passes — the deferral phrasing
    # still contains the keyword, no special-case needed.
    if ! printf '%s\n' "$ISSUE_BODY" | grep -qi "priority rationale"; then
        echo "✗ Issue #${ISSUE} is missing a '**Priority rationale:**' line." >&2
        echo "    Add a one-sentence rationale near the bottom of the body explaining the P0/P1/P2/P3 choice." >&2
        echo "    To intentionally defer, write e.g. '**Priority rationale:** (deferred to #X)' — the keyword alone passes." >&2
        FAILED=1
    else
        PR_RATIONALE_OK=$((PR_RATIONALE_OK + 1))
    fi
done

# Gate 4: stray closing-keyword in the would-be squash body (Issue #559). The
# detector scans PR title + body + every commit message for a `close|fix|resolve`
# + #N whose N is NOT in closingIssuesReferences. Fails OPEN (non-blocking) on a
# missing interpreter or gh error — a hiccup must not block a legitimate merge.
PYTHON="$(command -v python3 || command -v python || true)"
if [[ -z "$PYTHON" ]]; then
    echo "⚠ stray-closer check skipped (no python on PATH)." >&2
elif ! STRAY_OUT=$(REPO="$REPO" "$PYTHON" "$SCRIPT_DIR/../tools/ci/stray_closers.py" "$PR" 2>&1); then
    printf '%s\n' "$STRAY_OUT" >&2
    FAILED=1
elif [[ -n "$STRAY_OUT" ]]; then
    printf '%s\n' "$STRAY_OUT" >&2   # fail-open warning (e.g. gh error); non-blocking
fi

[[ "$FAILED" -eq 1 ]] && exit 1

MERGE_ARGS=("$MERGE_TYPE")
[[ -n "$DELETE_FLAG" ]] && MERGE_ARGS+=("$DELETE_FLAG")
gh pr merge "$PR" --repo "$REPO" "${MERGE_ARGS[@]}"

echo "✓ PR #${PR} merged (${TEST_PLAN_TOTAL} test-plan boxes ticked, ${AC_TOTAL} AC boxes ticked + ${PR_RATIONALE_OK}/${LINKED_COUNT} priority rationales present, no stray closing-keyword)."
