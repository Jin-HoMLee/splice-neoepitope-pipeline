#!/usr/bin/env bash
# scripts/audit_and_merge.sh
#
# Closure-ritual gate before `gh pr merge`. Refuses to merge if any of:
#   1. `- [ ]` remains on the PR body Test plan, OR
#   2. `- [ ]` remains under any linked Issue's Acceptance criteria, OR
#   3. Any linked Issue is missing a `**Priority rationale:**` line, OR
#   4. The would-be squash body (PR title + body + commit messages) contains a
#      closing keyword (`close|fix|resolve` + `#N`) targeting an Issue OUTSIDE
#      closingIssuesReferences — a silent unintended auto-close on merge, OR
#   5. A role-tagged linked Issue lacks a `## <merge-date>` lab-notebook entry in
#      research/lab_notebook/<role>.md referencing the PR/Issue (Issue #409).
# After those pass, a bot-review-offer gate (6) ensures a `@-claude review` was
# offered on the PR before merging (interactive prompt, or block + --skip-review-offer
# carve-out when non-interactive).
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
# closingIssuesReferences API does not surface. Lab-notebook gate (5) added via
# Issue #409 to move the post-hoc closure-audit notebook check to merge time
# (prevention over post-merge cleanup comment). Bot-review-offer gate (6) added
# via Issue #443 after PR #441 + PR #442 merged without the review offer.
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
#   1 — audit failed (unticked boxes, missing priority rationale, or a role-tagged
#       linked Issue missing its lab-notebook entry); OR the bot-review-offer gate
#       blocked (no review offered: cancelled interactively, or non-interactive with
#       no --skip-review-offer). The gaps/reason are printed
#   2 — usage error

set -euo pipefail

REPO="${REPO:-Jin-HoMLee/splice-neoepitope-pipeline}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    echo "Usage: $0 <PR_NUMBER> [--squash|--merge|--rebase] [--delete-branch|--no-delete-branch] [--skip-review-offer]" >&2
    exit 2
}

[[ $# -ge 1 ]] || usage
PR="$1"; shift
[[ "$PR" =~ ^[0-9]+$ ]] || usage

MERGE_TYPE="--squash"
DELETE_FLAG="--delete-branch"
SKIP_REVIEW_OFFER=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --squash|--merge|--rebase) MERGE_TYPE="$1" ;;
        --delete-branch)            DELETE_FLAG="--delete-branch" ;;
        --no-delete-branch)         DELETE_FLAG="" ;;
        --skip-review-offer)        SKIP_REVIEW_OFFER=1 ;;
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

# Gate 5: lab-notebook entry for role-tagged linked Issues (Issue #409). Moves
# the post-hoc closure-audit notebook check (tools/ci/closure_audit.py, run by
# the closure-audit workflow AFTER merge) to merge time, so a role-tagged PR
# cannot merge without its `## <today>` lab-notebook entry — prevention, not a
# post-merge cleanup comment. Single-sourced via tools/ci/lab_notebook_gate.py
# (a thin wrapper around closure_audit.audit_pr_pre_merge). The wrapper exits 1
# on a real gap, 0 when clear/exempt/skip-marked, and fails OPEN (exit 0 +
# warning) on a gh error — so a hiccup never blocks a legitimate merge. The entry
# is read from the working tree, so run this from the PR branch where the entry
# was written. Routine-ship bypass: `<!-- skip-lab-notebook: routine -->` in the
# PR body (same marker the post-hoc bot honors).
if [[ -z "$PYTHON" ]]; then
    echo "⚠ lab-notebook check skipped (no python on PATH)." >&2
elif ! NB_OUT=$("$PYTHON" "$SCRIPT_DIR/../tools/ci/lab_notebook_gate.py" "$PR" 2>&1); then
    printf '%s\n' "$NB_OUT" >&2
    FAILED=1
elif [[ -n "$NB_OUT" ]]; then
    printf '%s\n' "$NB_OUT" >&2   # fail-open warning (e.g. gh error); non-blocking
fi

[[ "$FAILED" -eq 1 ]] && exit 1

# Bot-review-offer gate (Issue #443). The closure-ritual checks have passed, so
# the PR is otherwise mergeable — now ensure a bot review (@-claude review) was
# offered on it before merging. Two same-morning slips (PR #441 + PR #442 merged
# without the offer) crossed the memory→mechanism escalation threshold; both were
# Claude-driven (non-interactive) merges, so this gate must bite in BOTH the
# interactive and non-interactive cases. Deterministic detection lives in
# tools/ci/bot_review_offer.py (prints OFFERED / NOT_OFFERED, fails open to
# OFFERED on a gh error). The rule's "skip for trivial PRs" carve-out keeps the
# judgment call with the operator (interactive option b / non-interactive
# --skip-review-offer flag) rather than coding a triviality heuristic.
#
#   --skip-review-offer set   → explicit trivial-PR bypass; proceed.
#   no python on PATH         → fail open; proceed (mirrors stray-closer gate).
#   review already offered    → proceed.
#   not offered, interactive  → prompt (a) offer now / (b) skip / (c) cancel.
#   not offered, non-interact → BLOCK (exit 1) with guidance — do NOT silently
#                               merge; this is the exact agent-side slip the gate
#                               exists to catch.
if [[ "$SKIP_REVIEW_OFFER" -eq 1 ]]; then
    echo "→ Bot-review-offer gate bypassed (--skip-review-offer; trivial PR)." >&2
elif [[ -z "$PYTHON" ]]; then
    # $PYTHON is resolved once in gate 4's setup and reused here — keep gate 4
    # before this block if either is ever refactored.
    echo "⚠ bot-review-offer check skipped (no python on PATH)." >&2
else
    REVIEW_STATUS="$("$PYTHON" "$SCRIPT_DIR/../tools/ci/bot_review_offer.py" "$PR")" || REVIEW_STATUS="OFFERED"
    if [[ "$REVIEW_STATUS" == "NOT_OFFERED" ]]; then
        if [[ -t 0 && -t 1 ]]; then
            echo "" >&2
            echo "No bot-review trigger (@-claude review) found on PR #${PR}. Offer one before merging?" >&2
            echo "  (a) offer review now — post the trigger comment, then merge" >&2
            echo "  (b) skip — trivial PR (typo / news_log / single-line doc fix)" >&2
            echo "  (c) cancel merge" >&2
            read -r -p "Choose [a/b/c]: " BR_CHOICE || BR_CHOICE=""
            case "$BR_CHOICE" in
                a|A)
                    gh pr comment "$PR" --repo "$REPO" --body "@claude review"
                    echo "✓ Posted review trigger on PR #${PR}." >&2
                    ;;
                b|B)
                    echo "→ Proceeding without a review offer (trivial PR)." >&2
                    ;;
                *)
                    echo "✗ Merge cancelled — no review offered." >&2
                    exit 1
                    ;;
            esac
        else
            echo "✗ PR #${PR} has no bot-review trigger (@-claude review), and this merge is" >&2
            echo "  running non-interactively (cannot prompt for a/b/c)." >&2
            echo "  → Offer a review first (ask the user, then post the trigger) and re-run, OR" >&2
            echo "  → pass --skip-review-offer for a genuinely trivial PR." >&2
            exit 1
        fi
    fi
fi

MERGE_ARGS=("$MERGE_TYPE")
[[ -n "$DELETE_FLAG" ]] && MERGE_ARGS+=("$DELETE_FLAG")
gh pr merge "$PR" --repo "$REPO" "${MERGE_ARGS[@]}"

echo "✓ PR #${PR} merged (${TEST_PLAN_TOTAL} test-plan boxes ticked, ${AC_TOTAL} AC boxes ticked + ${PR_RATIONALE_OK}/${LINKED_COUNT} priority rationales present, no stray closing-keyword, lab-notebook entry verified)."
