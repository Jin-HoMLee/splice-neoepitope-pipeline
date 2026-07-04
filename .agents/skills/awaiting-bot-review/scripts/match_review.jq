# Detection predicate for a landed bot code-review (Issue #864).
#
# Input: the JSON from `gh pr view <PR> --json comments`
#   { "comments": [ { "author": {"login": "..."}, "body": "...", "createdAt": "..." }, ... ] }
# Arg:   $w = watermark, an ISO-8601 UTC timestamp (e.g. the `@claude review`
#        trigger comment's createdAt). Only comments strictly newer than $w count,
#        so a PRIOR review on the same PR never re-matches.
#
# A review has LANDED when a comment is:
#   - authored by the bot (login == "claude"),
#   - newer than the watermark (createdAt > $w), and
#   - in the finished (not interim "working") state -> body contains "Claude finished".
#
# ISO-8601 UTC strings (identical `...Z` shape) compare correctly with `>`
# lexicographically, so no date parsing is needed.
#
# Output: the newest matching comment's body, or nothing (empty) when none has
# landed yet. An empty result is the caller's "keep polling" signal.
# The bot-review comment author login is exactly "claude" (NOT "claude[bot]" or
# "github-actions[bot]") - empirically verified 2026-07-04 against real landed
# reviews on PRs #985 and #993 (`gh pr view N --json comments --jq
# '[.comments[]|select(.body|test("Claude finished")).author.login]|unique'` -> ["claude"]).
# Pinned so a future reader doesn't "fix" it to a [bot] form and silently break
# detection. If a future GitHub-App migration changes it, re-verify and update here.
[
  .comments[]
  | select(
      .author.login == "claude"
      and .createdAt > $w
      and (.body | contains("Claude finished"))
    )
]
| sort_by(.createdAt)
| last
| if . == null then empty else .body end
