#!/usr/bin/env bash
#
# set_status.sh - move an Issue's card on board #9 ("JH M Lee Lab", user project
# under Jin-HoMLee) to a named Status.
#
# A thin, idempotent DRY wrapper over the raw
#   gh api graphql updateProjectV2ItemFieldValue
# mutation. Before this script the mutation appeared only inside the board hooks
# (post_gh_pr_create.py, post_gh_pr_review_request.py), so every manual status
# move meant hand-supplying the project / field / option IDs - the most
# error-prone step in the flow ("query, don't guess IDs"). The Status-name ->
# option-id map lives here, in ONE place; callers reference this script.
#
# Usage:
#   scripts/pm/set_status.sh <issue#> <status-name>
#
#   <status-name> one of:
#     Backlog | Ready | In progress | Ready for review | In review | Done | Epic
#
# Examples:
#   scripts/pm/set_status.sh 1024 "In progress"
#   scripts/pm/set_status.sh 1024 "Ready for review"
#
# Behaviour:
#   - resolves the Issue's project item id on board #9 itself,
#   - unknown status name / bad args           -> usage message, exit 2,
#   - Issue has no card on board #9            -> loud error, exit 1,
#   - card already in the target status        -> no-op, exit 0.
#
# The IDs below are for board #9. If a mutation 404s, they may have been
# regenerated (updateProjectV2Field regenerates option IDs) - re-query with:
#   gh api graphql -f query='query{ user(login:"Jin-HoMLee"){ projectV2(number:9){
#     id field(name:"Status"){ ... on ProjectV2SingleSelectField {
#       id options { id name } } } } }'
set -euo pipefail

PROJECT_ID="PVT_kwHOB17eGc4BSomP"
PROJECT_NUMBER=9
STATUS_FIELD_ID="PVTSSF_lAHOB17eGc4BSomPzhAHFf8"

# Canonical Status-name -> single-select option-id map for board #9.
# Single source of truth (mirrored, historically, in the quick-win-burndown
# skill appendix and the board hook constants). A `case` rather than an
# associative array so the script runs on macOS system bash 3.2 (no `declare -A`).
status_option_id() {
  case "$1" in
    "Backlog")          echo "f75ad846" ;;
    "Ready")            echo "61e4505c" ;;
    "In progress")      echo "47fc9ee4" ;;
    "Ready for review") echo "8bf9192f" ;;
    "In review")        echo "df73e18b" ;;
    "Done")             echo "98236657" ;;
    "Epic")             echo "9f872564" ;;
    *)                  return 1 ;;
  esac
}

usage() {
  cat >&2 <<EOF
Usage: ${0##*/} <issue#> <status-name>

  <status-name> one of:
    Backlog | Ready | In progress | Ready for review | In review | Done | Epic

Moves the Issue's card on board #9 to <status-name>.
Idempotent: already-in-target is a no-op.
EOF
  exit 2
}

[[ $# -eq 2 ]] || usage
ISSUE="$1"
STATUS="$2"

[[ "$ISSUE" =~ ^[0-9]+$ ]] || { echo "Error: issue# must be numeric, got '$ISSUE'." >&2; usage; }

OPTION_ID="$(status_option_id "$STATUS")" || { echo "Error: unknown status name '$STATUS'." >&2; usage; }

# owner/repo of the current repo (the Issue's card resolves from its own repo).
read -r OWNER REPO < <(gh repo view --json owner,name --jq '.owner.login + " " + .name')

# Resolve the Issue's board-#9 project item id + its current Status name.
# @tsv keeps the two fields tab-separated so a status with a space
# ("In progress") survives the read.
IFS=$'\t' read -r ITEM_ID CURRENT < <(
  gh api graphql -f query="
    query {
      repository(owner: \"$OWNER\", name: \"$REPO\") {
        issue(number: $ISSUE) {
          projectItems(first: 10) {
            nodes {
              id
              project { number }
              fieldValues(first: 20) {
                nodes {
                  ... on ProjectV2ItemFieldSingleSelectValue {
                    name
                    field { ... on ProjectV2FieldCommon { name } }
                  }
                }
              }
            }
          }
        }
      }
    }" --jq "
      .data.repository.issue.projectItems.nodes[]
      | select(.project.number == $PROJECT_NUMBER)
      | [.id, ((.fieldValues.nodes[]? | select(.field.name == \"Status\") | .name) // \"\")]
      | @tsv
  " 2>/dev/null
) || true

if [[ -z "${ITEM_ID:-}" ]]; then
  echo "Error: Issue #$ISSUE has no card on board #$PROJECT_NUMBER (nothing to move)." >&2
  exit 1
fi

if [[ "${CURRENT:-}" == "$STATUS" ]]; then
  echo "Issue #$ISSUE already in '$STATUS' - no-op."
  exit 0
fi

gh api graphql -f query='
  mutation($p:ID!, $i:ID!, $f:ID!, $o:String!) {
    updateProjectV2ItemFieldValue(input:{
      projectId:$p, itemId:$i, fieldId:$f,
      value:{ singleSelectOptionId:$o }
    }) { projectV2Item { id } }
  }' -f p="$PROJECT_ID" -f i="$ITEM_ID" -f f="$STATUS_FIELD_ID" -f o="$OPTION_ID" >/dev/null

echo "Issue #$ISSUE: '${CURRENT:-<none>}' -> '$STATUS'"
