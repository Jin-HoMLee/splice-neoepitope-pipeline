#!/usr/bin/env bash
#
# request_review.sh - request a bot code-review on a PR AND deterministically
# advance its board card + every linked Issue card to "In review", in ONE
# command, independent of the version-fragile PostToolUse hook.
#
# Why (Issue #1063): the `post_gh_pr_review_request.py` PostToolUse hook is the
# automatic path that flips a PR's linked Issue to "In review" when `@claude
# review` is posted. On 2026-07-06 it silently did NOT fire for PR #1057 + #1058
# (both linked Issues had to be flipped by hand) - not a code bug (it fired
# cleanly 4x on 2026-07-04) but a Claude-Code-version-dependent hook-execution
# lapse, the exact failure class CLAUDE.md documents as having flip-flopped
# across builds. A critical board move must not depend on a version-fragile hook
# as its SOLE path. This wrapper is the belt-and-suspenders the flow can rely on;
# the hook stays as the automatic path (this is idempotent, so both firing is a
# no-op).
#
# Sibling of scripts/pm/set_status.sh (Issue-card status wrapper),
# scripts/new_branch.sh, and scripts/audit_and_merge.sh.
#
# Usage:
#   scripts/pm/request_review.sh <pr#> [--flip-only]
#
#   <pr#>          the PR number (in the current repo).
#   --flip-only    flip the board cards only; do NOT post the `@claude review`
#                  comment. Use when the review was already requested (a
#                  re-flip), so a second bot review is not triggered.
#
# Examples:
#   scripts/pm/request_review.sh 1063               # post trigger + flip cards
#   scripts/pm/request_review.sh 1063 --flip-only   # flip cards only
#
# Behaviour:
#   - posts `@claude review` on the PR (unless --flip-only),
#   - flips the PR card AND every linked Issue card (resolved via the PR's
#     closingIssuesReferences) on board #9 to "In review" (option df73e18b),
#   - idempotent: a card already "In review" is a no-op,
#   - skips a card parked in "Epic" (a parent) or terminal "Done" (never moved),
#   - if the PR itself is not yet on board #9 (the sibling post_gh_pr_create hook
#     may also have failed to fire), it is added first, then flipped,
#   - per-card fail-safe: a `gh` error on one card is reported and does NOT abort
#     the remaining cards; each step prints a clear one-line outcome,
#   - unknown args / non-numeric PR -> usage message, exit 2.
#
# Exit codes:
#   0 - trigger posted (if requested) and all cards flipped or safely skipped,
#   1 - one or more steps failed (comment post, a board query, or a card flip);
#       every other step still ran. The failing step(s) are printed,
#   2 - usage error.
#
# The IDs below are for board #9 (mirrored from set_status.sh / the board hooks -
# a user-level project, so they are stable + repo-independent). If a mutation
# 404s they may have been regenerated; re-query per set_status.sh's header.
set -uo pipefail

PROJECT_ID="PVT_kwHOB17eGc4BSomP"
PROJECT_NUMBER=9
PROJECT_OWNER="Jin-HoMLee"
STATUS_FIELD_ID="PVTSSF_lAHOB17eGc4BSomPzhAHFf8"
IN_REVIEW_OPTION="df73e18b"

# Only PRs in these repos feed board #9 (mirrors the board hooks' TRACKED_REPOS).
# Used to gate the add-to-board fallback so an unrelated repo's PR is never
# boarded onto #9 by accident.
is_tracked() {
  case "$1/$2" in
    "$PROJECT_OWNER/splice-neoepitope-pipeline") return 0 ;;
    "$PROJECT_OWNER/claude-personas-splice-neoepitope-pipeline") return 0 ;;
    *) return 1 ;;
  esac
}

usage() {
  cat >&2 <<EOF
Usage: ${0##*/} <pr#> [--flip-only]

  <pr#>          PR number in the current repo.
  --flip-only    flip the board cards only; do not post the '@claude review'
                 comment (for a review already requested - avoids a double
                 bot-review).

Requests a bot review on the PR and flips the PR card + every linked Issue card
on board #9 to 'In review'. Idempotent; skips Epic/Done cards; per-card fail-safe.
EOF
  exit 2
}

# --- args --------------------------------------------------------------------

PR=""
FLIP_ONLY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --flip-only) FLIP_ONLY=1 ;;
    -h|--help)   usage ;;
    -*)          echo "Error: unknown flag '$1'." >&2; usage ;;
    *)
      if [[ -n "$PR" ]]; then echo "Error: unexpected extra argument '$1'." >&2; usage; fi
      PR="$1" ;;
  esac
  shift
done

[[ -n "$PR" ]] || usage
[[ "$PR" =~ ^[0-9]+$ ]] || { echo "Error: pr# must be numeric, got '$PR'." >&2; usage; }

# --- board helpers (mirror set_status.sh) ------------------------------------

# resolve_card <owner> <repo> <number>
# Prints one TAB-separated line "ITEM_ID<TAB>CURRENT_STATUS<TAB>TITLE" when the
# Issue-or-PR #number has a card on board #9, else nothing. `issueOrPullRequest`
# resolves either a PR (the PR card) or an Issue (a linked-issue card) with one
# query shape - the number space is shared per repo, so it is unambiguous.
# Returns gh's exit status (non-zero => query failed, distinct from "no card").
resolve_card() {
  local owner="$1" repo="$2" num="$3"
  gh api graphql -f query="
    query {
      repository(owner: \"$owner\", name: \"$repo\") {
        issueOrPullRequest(number: $num) {
          ... on Issue {
            title
            projectItems(first: 10) { nodes { id project { number }
              fieldValues(first: 20) { nodes {
                ... on ProjectV2ItemFieldSingleSelectValue {
                  name field { ... on ProjectV2FieldCommon { name } } } } } } }
          }
          ... on PullRequest {
            title
            projectItems(first: 10) { nodes { id project { number }
              fieldValues(first: 20) { nodes {
                ... on ProjectV2ItemFieldSingleSelectValue {
                  name field { ... on ProjectV2FieldCommon { name } } } } } } }
          }
        }
      }
    }" --jq "
      .data.repository.issueOrPullRequest as \$n
      | \$n.projectItems.nodes[]?
      | select(.project.number == $PROJECT_NUMBER)
      | [.id, ((.fieldValues.nodes[]? | select(.field.name == \"Status\") | .name) // \"\"), (\$n.title // \"\")]
      | @tsv
  "
}

# add_to_board <content-url> -> prints the new/existing project-item id.
# `gh project item-add` no-ops server-side if the content is already on the board
# (dedup by content) and returns the existing item either way, so it is safe.
add_to_board() {
  local out
  out="$(gh project item-add "$PROJECT_NUMBER" --owner "$PROJECT_OWNER" --url "$1" --format json)" || return 1
  jq -r '.id' <<<"$out"
}

# set_in_review <item-id> -> move the card to "In review". Same mutation as
# set_status.sh / the board hooks.
set_in_review() {
  gh api graphql -f query='
    mutation($p:ID!, $i:ID!, $f:ID!, $o:String!) {
      updateProjectV2ItemFieldValue(input:{
        projectId:$p, itemId:$i, fieldId:$f,
        value:{ singleSelectOptionId:$o }
      }) { projectV2Item { id } }
    }' -f p="$PROJECT_ID" -f i="$1" -f f="$STATUS_FIELD_ID" -f o="$IN_REVIEW_OPTION" >/dev/null
}

# flip_card <owner> <repo> <number> <kind> <content-url|-> <add-if-missing 0|1>
# Resolve the card, apply the skip rules (In review / Epic / Done), and flip to
# "In review" otherwise. Prints exactly one outcome line. Returns 0 on a flip or
# a legitimate skip; 1 on an error (so the caller can mark overall failure). Every
# gh call is guarded, so a hiccup on one card never aborts the run.
flip_card() {
  local owner="$1" repo="$2" num="$3" kind="$4" content_url="$5" add="$6"
  local result item current

  if ! result="$(resolve_card "$owner" "$repo" "$num" 2>&1)"; then
    echo "warn  $kind #$num: board query failed: $result" >&2
    return 1
  fi
  # Third @tsv field (title) is read into `_` - unused here, but it must be
  # consumed so a title with spaces does not spill into `current`.
  IFS=$'\t' read -r item current _ <<<"$result"

  if [[ -z "${item:-}" ]]; then
    if [[ "$add" == "1" && "$content_url" != "-" ]]; then
      if ! item="$(add_to_board "$content_url" 2>&1)"; then
        echo "warn  $kind #$num: not on board #9 and item-add failed: $item" >&2
        return 1
      fi
      current=""
      echo "add   $kind #$num: added to board #9"
    else
      echo "warn  $kind #$num: no card on board #9 (skipped)" >&2
      return 1
    fi
  fi

  case "${current:-}" in
    "In review") echo "noop  $kind #$num: already 'In review'"; return 0 ;;
    "Epic")      echo "skip  $kind #$num: parked in 'Epic' (not moved)"; return 0 ;;
    "Done")      echo "skip  $kind #$num: terminal 'Done' (not moved)"; return 0 ;;
  esac

  local mut_err
  if ! mut_err="$(set_in_review "$item" 2>&1)"; then
    echo "warn  $kind #$num: flip to 'In review' failed: $mut_err" >&2
    return 1
  fi
  echo "flip  $kind #$num: '${current:-<none>}' -> 'In review'"
  return 0
}

# --- resolve the PR ----------------------------------------------------------

if ! PR_JSON="$(gh pr view "$PR" --json url,number,title,closingIssuesReferences 2>&1)"; then
  echo "Error: could not resolve PR #$PR (not a PR in this repo / gh error):" >&2
  printf '%s\n' "$PR_JSON" >&2
  exit 1
fi

PR_URL="$(jq -r '.url // ""' <<<"$PR_JSON")"
PR_TITLE="$(jq -r '.title // ""' <<<"$PR_JSON")"
# newline-separated linked Issue numbers (closingIssuesReferences is same-repo).
LINKED="$(jq -r '.closingIssuesReferences[]?.number // empty' <<<"$PR_JSON")"

if [[ -z "$PR_URL" ]]; then
  echo "Error: PR #$PR view returned no url (unexpected gh output)." >&2
  exit 1
fi

# owner/repo from the PR's own URL (authoritative), so linked-Issue lookups use
# the PR's repo rather than assuming cwd.
REST="${PR_URL#https://github.com/}"
OWNER="${REST%%/*}"
REST2="${REST#*/}"
REPO="${REST2%%/*}"

echo "PR #$PR ($OWNER/$REPO): ${PR_TITLE:-?}"

OVERALL_RC=0

# --- 1. request the review (unless --flip-only) ------------------------------

if [[ "$FLIP_ONLY" -eq 1 ]]; then
  echo "--flip-only: not posting '@claude review'."
else
  if comment_err="$(gh pr comment "$PR" --repo "$OWNER/$REPO" --body "@claude review" 2>&1)"; then
    echo "post  PR #$PR: posted '@claude review'"
  else
    echo "warn  PR #$PR: failed to post '@claude review': $comment_err" >&2
    OVERALL_RC=1
  fi
fi

# --- 2. flip the PR card (add to board first if the create hook never did) ----

ADD_PR=0
if is_tracked "$OWNER" "$REPO"; then ADD_PR=1; fi
flip_card "$OWNER" "$REPO" "$PR" "PR" "$PR_URL" "$ADD_PR" || OVERALL_RC=1

# --- 3. flip every linked Issue card -----------------------------------------

if [[ -z "$LINKED" ]]; then
  echo "note  PR #$PR has no linked Issues (closingIssuesReferences empty)."
else
  while IFS= read -r ISSUE; do
    [[ -n "$ISSUE" ]] || continue
    # Issues are never auto-added to the board here (add-if-missing = 0): an
    # unboarded linked Issue is an anomaly to surface, not to silently create.
    flip_card "$OWNER" "$REPO" "$ISSUE" "Issue" "-" 0 || OVERALL_RC=1
  done <<<"$LINKED"
fi

exit "$OVERALL_RC"
