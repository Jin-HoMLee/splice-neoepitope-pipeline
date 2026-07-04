#!/usr/bin/env bash
# Poll a PR for a landed bot code-review and print its verdict (Issue #864).
#
# Launch this as a BACKGROUND task right after posting `@claude review`, so the
# harness re-invokes the agent when it exits and the verdict is surfaced without
# the user having to notice the review landed. It never merges - relaying the
# verdict (and any blocking findings) is the whole job.
#
# Usage:
#   poll_bot_review.sh <PR> [--since ISO8601] [--timeout-min N] [--interval-sec N]
#
# Options:
#   --since        Watermark; only a review newer than this counts. Defaults to
#                  "now" (UTC), so a review already on the PR is ignored. Pass the
#                  `@claude review` trigger comment's createdAt for a precise mark.
#   --timeout-min  Give up after N minutes (default 25); exit 3 = notify-and-stop,
#                  no auto re-ping (avoids double-review noise).
#   --interval-sec Poll cadence in seconds (default 30).
#
# Exit codes: 0 = review landed (verdict on stdout); 2 = usage error;
#             3 = timed out with no landed review;
#             4 = gave up after too many consecutive gh/jq failures (check auth /
#                 that the PR exists) - a persistent error surfaced loudly rather
#                 than masquerading as a clean timeout.
#
# Detection lives in the sibling match_review.jq so it is unit-tested independently
# of this polling shell (see tests/). The `gh ... | jq -f` form is deliberate: an
# earlier inline `gh --jq --arg` poller was rejected by gh ("accepts at most 1
# arg") and silently timed out every run - piping raw JSON to real jq avoids that
# whole failure class.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MATCH_JQ="$SCRIPT_DIR/match_review.jq"

PR=""
SINCE=""
TIMEOUT_MIN=25
INTERVAL_SEC=30

die() { echo "poll_bot_review: $1" >&2; exit 2; }

# Reject a missing value or one that looks like the next flag, so
# `--since --timeout-min 5` errors instead of silently taking "--timeout-min" as
# the watermark (which, via jq's lexicographic compare, sorts below every real
# timestamp and would false-match everything). (PR #993 review, nit 1.)
need_val() { case "${2:-}" in ""|-*) die "option $1 requires a value";; esac; }

while [ $# -gt 0 ]; do
  case "$1" in
    --since)        need_val "$1" "${2:-}"; SINCE="$2"; shift 2 ;;
    --timeout-min)  need_val "$1" "${2:-}"; TIMEOUT_MIN="$2"; shift 2 ;;
    --interval-sec) need_val "$1" "${2:-}"; INTERVAL_SEC="$2"; shift 2 ;;
    -h|--help)      awk 'NR==1{next} /^set -euo/{exit} /^#/{sub(/^# ?/,"");print}' "${BASH_SOURCE[0]}"; exit 0 ;;
    -*)             die "unknown option: $1" ;;
    *)              [ -z "$PR" ] || die "unexpected extra argument: $1"; PR="$1"; shift ;;
  esac
done

[ -n "$PR" ] || die "missing required <PR> argument"
[[ "$PR" =~ ^[0-9]+$ ]] || die "PR must be a number, got: $PR"
[[ "$TIMEOUT_MIN" =~ ^[0-9]+$ ]] || die "--timeout-min must be a non-negative integer, got: $TIMEOUT_MIN"
[[ "$INTERVAL_SEC" =~ ^[0-9]+$ ]] || die "--interval-sec must be a non-negative integer, got: $INTERVAL_SEC"
[ -f "$MATCH_JQ" ] || die "detection predicate not found at $MATCH_JQ"

if [ -z "$SINCE" ]; then
  SINCE="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
fi

# Give up only after this many *consecutive* failed polls, so a transient gh blip
# (network, secondary rate-limit, a 5xx) does not abort an unattended background
# watch, while a persistent error (bad auth, deleted PR) still surfaces loudly
# instead of masquerading as a clean timeout. (PR #993 review, finding 1.)
MAX_CONSECUTIVE_FAILURES=5

deadline=$(( $(date +%s) + TIMEOUT_MIN * 60 ))
fails=0

while true; do
  # Pipe raw JSON to real jq (NOT `gh --jq`) so `jq --arg` works. The `if` guards
  # the assignment so a non-zero pipeline (gh blip, under `set -euo pipefail`)
  # does NOT abort the script - we count it and keep polling to the timeout.
  if verdict="$(gh pr view "$PR" --json comments | jq -r --arg w "$SINCE" -f "$MATCH_JQ")"; then
    fails=0
  else
    verdict=""
    fails=$(( fails + 1 ))
    if [ "$fails" -ge "$MAX_CONSECUTIVE_FAILURES" ]; then
      echo "poll_bot_review: $fails consecutive gh/jq failures polling PR #$PR; giving up (check auth / that the PR exists)." >&2
      exit 4
    fi
  fi
  if [ -n "$verdict" ]; then
    echo "Bot review landed on PR #$PR:"
    echo "$verdict"
    exit 0
  fi
  if [ "$(date +%s)" -ge "$deadline" ]; then
    echo "poll_bot_review: no bot review on PR #$PR within ${TIMEOUT_MIN}m (watermark ${SINCE}); stopping." >&2
    exit 3
  fi
  sleep "$INTERVAL_SEC"
done
