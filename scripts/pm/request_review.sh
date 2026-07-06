#!/usr/bin/env bash
#
# request_review.sh - request a bot code-review on a PR AND advance its board
# cards to "In review", deterministically, in one command.
#
# Why: the "In review" flip is normally done by the post_gh_pr_review_request.py
# PostToolUse hook (Issue #996). That hook is version-fragile - it silently did
# not fire for PR #1057 / PR #1058 (2026-07-06), stranding the linked Issues at
# "Ready for review", and /hooks confirmed the hooks were registered yet not
# executing. A critical board move must not depend on a hook as its sole path.
# This wrapper does the same move deterministically; the hook stays as the
# automatic path and this is the belt-and-suspenders. Issue #1063; sibling of
# set_status.sh / new_branch.sh / audit_and_merge.sh.
#
# Usage:
#   scripts/pm/request_review.sh <pr#-or-url> [--flip-only]
#
#     (default)     post "@claude review" on the PR + flip the PR card and every
#                   linked Issue card (via closingIssuesReferences) to In review.
#     --flip-only   flip the board only; do NOT post a comment. Use when the
#                   review was already requested (avoids a double bot-review).
#
# Behaviour:
#   - flips the PR's own card AND each linked Issue's card on board #9,
#   - skips a card already In review / Done / parked in Epic (idempotent, never
#     overwrites a terminal/epic state) - mirrors the hook's SKIP_STATUSES,
#   - fails SAFE per card: a gh error on one card is reported and the rest still
#     run; the comment post is best-effort likewise. Exit 1 if anything failed.
#
# Board #9 IDs: if a mutation 404s they may have been regenerated - re-query per
# set_status.sh's header recipe (the canonical Status option-id map lives there).
#
# NOTE: -e is deliberately OFF - per-card fail-safe is handled explicitly so one
# bad card never aborts the rest.
set -uo pipefail

PROJECT_NUMBER=9
PROJECT_ID="PVT_kwHOB17eGc4BSomP"
STATUS_FIELD_ID="PVTSSF_lAHOB17eGc4BSomPzhAHFf8"
IN_REVIEW_OPTION="df73e18b"
REVIEW_TRIGGER="@claude review"

usage() {
  cat >&2 <<EOF
Usage: ${0##*/} <pr#-or-url> [--flip-only]

  (default)    post "$REVIEW_TRIGGER" on the PR + flip the PR card and every
               linked Issue card to In review on board #$PROJECT_NUMBER.
  --flip-only  flip the board only; do not post a comment (review already asked).
EOF
  exit 2
}

PR=""
FLIP_ONLY=0
for arg in "$@"; do
  case "$arg" in
    --flip-only) FLIP_ONLY=1 ;;
    -h|--help)   usage ;;
    -*)          echo "Error: unknown flag '$arg'." >&2; usage ;;
    *)           if [[ -n "$PR" ]]; then echo "Error: multiple PR refs given." >&2; usage; fi; PR="$arg" ;;
  esac
done
[[ -n "$PR" ]] || usage

# owner/repo of the current repo (cards resolve from the cwd repo, like set_status.sh).
OWNER_REPO="$(gh repo view --json owner,name --jq '.owner.login + " " + .name' 2>&1)" \
  || { echo "Error: could not resolve the current repo (gh unauthenticated / not a repo):" >&2; printf '%s\n' "$OWNER_REPO" >&2; exit 1; }
read -r OWNER REPO <<<"$OWNER_REPO"

# Resolve the PR: its number + linked Issue numbers (closingIssuesReferences).
PR_JSON="$(gh pr view "$PR" --json number,closingIssuesReferences 2>&1)" \
  || { echo "Error: could not resolve PR '$PR' in $OWNER/$REPO:" >&2; printf '%s\n' "$PR_JSON" >&2; exit 1; }
PR_NUMBER="$(printf '%s' "$PR_JSON" | jq -r '.number')"

# bash 3.2 (macOS system bash): no mapfile - collect linked issue numbers by loop.
LINKED=()
while IFS= read -r n; do
  [[ -n "$n" ]] && LINKED+=("$n")
done < <(printf '%s' "$PR_JSON" | jq -r '.closingIssuesReferences[]?.number')

# flip_card <number> <label>: resolve the board-#9 item id + current Status for
# an Issue OR PR by number and set it to In review, unless it is already there /
# Done / Epic. Echoes the outcome. Returns 1 on a hard gh failure (aggregated).
flip_card() {
  local number="$1" label="$2" result item_id current
  if ! result="$(
    gh api graphql -f query="
      query {
        repository(owner: \"$OWNER\", name: \"$REPO\") {
          issueOrPullRequest(number: $number) {
            ... on Issue {
              projectItems(first: 10) { nodes { id project { number }
                fieldValues(first: 20) { nodes {
                  ... on ProjectV2ItemFieldSingleSelectValue {
                    name field { ... on ProjectV2FieldCommon { name } } } } } } }
            }
            ... on PullRequest {
              projectItems(first: 10) { nodes { id project { number }
                fieldValues(first: 20) { nodes {
                  ... on ProjectV2ItemFieldSingleSelectValue {
                    name field { ... on ProjectV2FieldCommon { name } } } } } } }
            }
          }
        }
      }" --jq "
        .data.repository.issueOrPullRequest.projectItems.nodes[]
        | select(.project.number == $PROJECT_NUMBER)
        | [.id, ((.fieldValues.nodes[]? | select(.field.name == \"Status\") | .name) // \"\")]
        | @tsv
    " 2>&1
  )"; then
    echo "  ! $label: board query failed: $result" >&2
    return 1
  fi
  if [[ -z "$result" ]]; then
    echo "  - $label: no card on board #$PROJECT_NUMBER, skipped."
    return 0
  fi
  IFS=$'\t' read -r item_id current <<<"$result"
  case "$current" in
    "In review")   echo "  = $label: already In review, no-op." ; return 0 ;;
    "Done"|"Epic") echo "  = $label: $current (terminal/epic), left as-is." ; return 0 ;;
  esac
  # $p/$i/$f/$o below are GraphQL variables (passed via -f), not shell vars, so
  # the single quotes are intentional - suppress SC2016's expand-in-quotes hint.
  # shellcheck disable=SC2016
  if ! gh api graphql -f query='
    mutation($p:ID!, $i:ID!, $f:ID!, $o:String!) {
      updateProjectV2ItemFieldValue(input:{
        projectId:$p, itemId:$i, fieldId:$f, value:{ singleSelectOptionId:$o }
      }) { projectV2Item { id } }
    }' -f p="$PROJECT_ID" -f i="$item_id" -f f="$STATUS_FIELD_ID" -f o="$IN_REVIEW_OPTION" >/dev/null 2>&1; then
    echo "  ! $label: failed to set In review (board option IDs may be regenerated)." >&2
    return 1
  fi
  echo "  + $label: '${current:-<none>}' -> In review."
  return 0
}

FAILED=0

# 1. Request the review (best-effort) unless --flip-only.
if [[ "$FLIP_ONLY" -eq 0 ]]; then
  if gh pr comment "$PR" --body "$REVIEW_TRIGGER" >/dev/null 2>&1; then
    echo "Requested review on PR #$PR_NUMBER ('$REVIEW_TRIGGER' posted)."
  else
    echo "! Failed to post '$REVIEW_TRIGGER' on PR #$PR_NUMBER." >&2
    FAILED=1
  fi
else
  echo "--flip-only: not posting a comment on PR #$PR_NUMBER."
fi

# 2. Flip the PR card + each linked Issue card to In review.
flip_card "$PR_NUMBER" "PR #$PR_NUMBER" || FAILED=1
if [[ ${#LINKED[@]} -eq 0 ]]; then
  echo "  (PR #$PR_NUMBER has no linked Issues via closingIssuesReferences.)"
else
  for iss in "${LINKED[@]}"; do
    flip_card "$iss" "Issue #$iss" || FAILED=1
  done
fi

exit "$FAILED"
