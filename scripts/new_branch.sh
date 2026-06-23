#!/usr/bin/env bash
# new_branch.sh — create the canonical <type>/<role>/issue-N-slug branch.
# Wraps `gh issue develop` (preserves the Issue<->branch Development-panel link).
# Spec: docs/superpowers/specs/2026-06-22-new-branch-helper-design.md (Issue #578)
set -euo pipefail

REPO_OWNER="Jin-HoMLee"
REPO_NAME="splice-neoepitope-pipeline"

die() { echo "new_branch.sh: $*" >&2; exit 1; }

# --- pure helpers (no I/O) --------------------------------------------------

# Leading Conventional-Commit type of an Issue title: "feat(scripts): x" -> "feat".
parse_type_from_title() {
  local title="$1"
  printf '%s' "$title" | sed -n 's/^\([a-zA-Z]\{1,\}\)\((.*)\)\{0,1\}:.*/\1/p' | tr '[:upper:]' '[:lower:]'
}

# Sanitize a slug: lowercase; spaces/underscores -> -; drop non [a-z0-9-];
# collapse repeated -; trim leading/trailing -.
sanitize_slug() {
  local s="$1"
  s=$(printf '%s' "$s" | tr '[:upper:]' '[:lower:]' | tr ' _' '-')
  s=$(printf '%s' "$s" | tr -cd 'a-z0-9-')
  s=$(printf '%s' "$s" | sed -E 's/-+/-/g; s/^-//; s/-$//')
  printf '%s' "$s"
}

# assemble_branch <type> <role> <issue#> <slug>
assemble_branch() {
  printf '%s/%s/issue-%s-%s' "$1" "$2" "$3" "$4"
}

# --- effectful ---------------------------------------------------------------

# fetch_issue <N> -> prints "title\trole1,role2\tsubTotal"
fetch_issue() {
  local n="$1" out
  # $number is a GraphQL variable (bound via -F number=), intentionally not a
  # shell var — the single-quoted query must stay unexpanded.
  # shellcheck disable=SC2016
  out=$(gh api graphql -F number="$n" -f query='
    query($number: Int!) {
      repository(owner: "'"$REPO_OWNER"'", name: "'"$REPO_NAME"'") {
        issue(number: $number) {
          title
          labels(first: 50) { nodes { name } }
          subIssuesSummary { total }
        }
      }
    }')
  # Null-tolerant: a non-existent issue returns issue:null (exit 0) — default to
  # empty/0 so the caller can detect "not found" instead of jq aborting on null.
  local title roles sub
  title=$(printf '%s' "$out" | jq -r '.data.repository.issue.title // ""')
  roles=$(printf '%s' "$out" | jq -r '[((.data.repository.issue.labels.nodes // [])[].name) | select(startswith("role:")) | ltrimstr("role:")] | join(",")')
  sub=$(printf '%s' "$out" | jq -r '.data.repository.issue.subIssuesSummary.total // 0')
  # \037 (Unit Separator) — non-whitespace so empty fields are preserved by read
  # (a TAB delimiter collapses consecutive separators, eating an empty role field).
  printf '%s\037%s\037%s\n' "$title" "$roles" "$sub"
}

# --- main --------------------------------------------------------------------

# require_value <flag> <value> -> dies if value is missing or itself a flag,
# so `--type --dry-run` can't swallow --dry-run (which would clear the dry-run
# guard and attempt a real branch-create).
require_value() {
  case "${2:-}" in
    ""|-*) die "$1 requires a value" ;;
  esac
}

# pick_role <override> <comma-list> -> echoes the single role or dies.
pick_role() {
  local override="$1" list="$2"
  if [ -n "$override" ]; then printf '%s' "$override"; return; fi
  case "$list" in
    "") die "no role label on issue; pass --role <pm|scientist|developer>" ;;
    *,*) die "multiple roles: ${list//,/, } — pass --role <one>" ;;
    *) printf '%s' "$list" ;;
  esac
}

DRY_RUN=0
NO_ISSUE=0
TYPE_OVERRIDE=""
ROLE_OVERRIDE=""
POSITIONAL=()
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) DRY_RUN=1; shift ;;
    --no-issue) NO_ISSUE=1; shift ;;
    --type) require_value --type "${2:-}"; TYPE_OVERRIDE="$2"; shift 2 ;;
    --role) require_value --role "${2:-}"; ROLE_OVERRIDE="$2"; shift 2 ;;
    *) POSITIONAL+=("$1"); shift ;;
  esac
done

# --- issueless fallback: --no-issue <type> <role> <slug> --------------------
if [ "$NO_ISSUE" -eq 1 ]; then
  # type+role are positional here; --type/--role would silently shift a slot.
  [ -z "$TYPE_OVERRIDE" ] && [ -z "$ROLE_OVERRIDE" ] \
    || die "--type/--role are not valid with --no-issue (type and role are positional)"
  TYPE_RAW="${POSITIONAL[0]:-}"
  ROLE_RAW="${POSITIONAL[1]:-}"
  SLUG_RAW="${POSITIONAL[2]:-}"
  { [ -n "$TYPE_RAW" ] && [ -n "$ROLE_RAW" ] && [ -n "$SLUG_RAW" ]; } \
    || die "usage: new_branch.sh --no-issue <type> <role> <short-slug> [--dry-run]"
  # sanitize all three for a canonical, valid git ref (the issue-linked path
  # gets a lowercased title-type + label-role for free; do the same here).
  TYPE=$(sanitize_slug "$TYPE_RAW")
  ROLE=$(sanitize_slug "$ROLE_RAW")
  SLUG=$(sanitize_slug "$SLUG_RAW")
  { [ -n "$TYPE" ] && [ -n "$ROLE" ] && [ -n "$SLUG" ]; } \
    || die "type, role, and slug must each be non-empty after sanitization"
  BRANCH="$TYPE/$ROLE/$SLUG"
  if [ "$DRY_RUN" -eq 1 ]; then printf '%s\n' "$BRANCH"; exit 0; fi
  git fetch origin --quiet
  git checkout -b "$BRANCH" origin/main   # fresh base (not local HEAD)
  exit 0
fi

# --- issue-linked path ------------------------------------------------------
ISSUE="${POSITIONAL[0]:-}"
SLUG_RAW="${POSITIONAL[1]:-}"

[ -n "$ISSUE" ] || die "usage: new_branch.sh <issue#> <short-slug> [--type T] [--role R] [--dry-run]"
case "$ISSUE" in (*[!0-9]*|'') die "issue number must be a positive integer";; esac

# Explicit capture so a failed fetch or a non-existent issue gives a clean error
# (not a raw jq trace): a missing issue returns issue:null → empty title here.
DATA=$(fetch_issue "$ISSUE") || die "could not fetch issue #$ISSUE (gh error?)"
IFS=$'\037' read -r TITLE ROLES SUBTOTAL <<<"$DATA"
[ -n "$TITLE" ] || die "issue #$ISSUE not found"

# parent guard (replicates the gh-issue-develop PreToolUse hook, which can't see
# the nested gh call): refuse epics — branching off them auto-closes the parent.
if [ "${SUBTOTAL:-0}" -gt 0 ]; then
  die "issue #$ISSUE is a parent/epic ($SUBTOTAL sub-issues) — branch off a leaf sub-issue, or file a closure sub-issue under it"
fi

# type
TYPE="$TYPE_OVERRIDE"
[ -n "$TYPE" ] || TYPE=$(parse_type_from_title "$TITLE")
[ -n "$TYPE" ] || die "could not derive type from title; pass --type"

# role
ROLE=$(pick_role "$ROLE_OVERRIDE" "$ROLES")

# slug
[ -n "$SLUG_RAW" ] || die "missing slug for issue #$ISSUE ($TITLE); suggested: new_branch.sh $ISSUE <your-slug>"
SLUG=$(sanitize_slug "$SLUG_RAW")
[ -n "$SLUG" ] || die "slug is empty after sanitization"

BRANCH=$(assemble_branch "$TYPE" "$ROLE" "$ISSUE" "$SLUG")

if [ "$DRY_RUN" -eq 1 ]; then
  printf '%s\n' "$BRANCH"
  exit 0
fi

# gh issue develop --base main branches server-side from fresh remote main, so
# no local `git fetch` is needed here (it wouldn't affect the server-side base).
gh issue develop "$ISSUE" --name "$BRANCH" --base main --checkout
