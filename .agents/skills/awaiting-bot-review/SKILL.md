---
name: awaiting-bot-review
description: Use immediately after posting the `@claude review` trigger on a PR (or the moment you are about to request a bot code-review) to await the verdict hands-free. Launches a background poller that re-invokes you when the bot's review lands and surfaces its verdict + any blocking findings, so you never wait for the user to say "the review is back". Triggers on requesting a bot/`@claude` PR review, "await the review", "watch for the review", "poll for the review verdict". Never merges - relaying the verdict is the whole job.
---

# awaiting-bot-review

After you post `@claude review` on a PR, the GitHub Action takes ~5-7 minutes to
post its verdict. Without this skill the user has to notice it landed and prompt
you. This skill makes the wait hands-free: launch a background poller on the
review request, and the harness re-invokes you when it exits so you can relay the
verdict yourself.

## When this fires

- You just posted (or are about to post) the `@claude review` trigger on a PR.
- The user asks you to await / watch / poll for a bot review verdict.

Scope is **bot PR reviews only** - this is not a general "wait for any CI/any
comment" tool.

## Procedure

1. **Capture a watermark** - the `createdAt` of the trigger comment you just
   posted (so a *prior* review on the same PR can't false-match). Get it with:

   ```bash
   gh pr view <PR> --json comments \
     --jq '[.comments[] | select(.body|test("@claude review"))] | last | .createdAt'
   ```

   (This `gh --jq` form is safe *because* it passes no `--arg`; the poller uses
   `gh ... | jq -f` precisely because `gh --jq --arg` is the trap - do not "fix"
   one to match the other.)

   If you can't get it cleanly, omit `--since`; the poller defaults the watermark
   to "now", which is safe (it just ignores anything already on the PR).

2. **Launch the bundled poller as a background task** (this is what makes the
   wait hands-free - a hook can only *remind*, it cannot deliver the async
   "review landed" re-invocation; only a harness-tracked background task can):

   ```bash
   .agents/skills/awaiting-bot-review/scripts/poll_bot_review.sh <PR> --since <watermark>
   ```

   Run it with the Bash tool's background mode (`run_in_background: true`). Keep
   working; the harness re-invokes you when the script exits.

3. **On exit, relay - never act.**
   - Exit `0`: the verdict body is on stdout. Summarize it for the user: the
     overall call (LGTM / requests-changes) and any **blocking** findings, with
     your take on each. **Do not merge** - merging stays a user-gated decision.
   - Exit `3`: timed out (default 25 min) with no landed review. Tell the user it
     hasn't come back yet; do not auto re-ping (that risks a double review). They
     can ask you to relaunch.
   - Exit `2`: a usage error in how the poller was invoked - fix the args and
     relaunch.

## The bundled script

`scripts/poll_bot_review.sh <PR> [--since ISO8601] [--timeout-min N] [--interval-sec N]`

Detection (in `scripts/match_review.jq`, unit-tested independently): a comment
authored by `claude`, newer than the watermark, whose body contains
`"Claude finished"` (the finished state, not the interim "working" state).

**Why a tested bundled script, not inline polling:** an inline poller written on
PR #862 used `gh ... --jq` with jq's `--arg`, which `gh` rejects
(`accepts at most 1 arg`), so it errored every iteration and timed out instead of
catching the review. The script pipes raw JSON to real `jq -f`, and
`tools/ci/test_poll_bot_review.py` locks in both the predicate and the arg
parsing so that whole failure class can't regress. Re-deriving the poll inline
each time invites the bug back - always use this script.
