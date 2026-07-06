---
name: quick-win-burndown
description: Use when the user asks to run a quick-win burn-down, clear quick wins off the board, or bring a batch of small role-tagged Issues to the merge gate together. Selects premise-verified quick wins for one role, brings each through branch/build/test/PR/review to the merge gate (stopping BEFORE merge, so they pile under In review for one final human look), loops until the field is exhausted, and surfaces the boundary case instead of grinding marginal items. Triggers on "any quick wins?", "burn down the quick wins", "clear the board", "pile them under In review", "bring quick wins to the gate", or any request to work through several small Issues autonomously and stop short of merging.
---

# quick-win-burndown

Clear a role's easy board work in one pass: repeatedly pick the next real quick win, take it all the way to the merge command but not through it, and leave it parked at **In review** for the human to glance at and merge. Keep going until nothing clean is left, then stop and say so.

The mechanical PR flow here is already covered by other tools (`new_branch.sh`, `audit_and_merge.sh`, the `awaiting-bot-review` skill, the merge-gate cadence). What this skill actually earns its keep on is the two hard parts: **selecting** which board items are genuinely quick wins (most that look like it aren't), and **knowing when to stop** instead of churning marginal work to hit a number.

## When this fires

- The user asks for "quick wins", a "burn-down", to "clear the board", or to batch small Issues to the gate.
- Implicitly: any "keep bringing me small wins until you run out" request scoped to a role.

Default the role to the one this workspace implements (its `role:` label). Ask only if it's genuinely ambiguous.

## The loop

```
select the candidate field  ->  pre-flight triage  ->  pick the next clean-go  ->  bring it to the gate (stop before merge)  ->  repeat
                                 (batch the human-input                 |
                                  questions once, up front)             +--> when the next-best isn't a clean quick win: STOP and surface it
```

Normally the commitment to *start* each Issue is the human's, and the merge-gate cadence only governs how far to run *once started*. The burn-down request is what changes that: it's a **standing grant to keep starting** quick wins without asking each time. That grant is only safe because of two guardrails - the per-candidate freshness ratification (Step 1) that stops you starting stale work, and the stop condition (Step 3) that ends the grant the moment the field stops being clean. Without those, "keep starting things autonomously" is exactly the anti-pattern the start-is-human-gated rule guards against.

## Step 1 - Select the candidate field (the part worth thinking about)

Scan the board **role-scoped and paginated** - a single-page query silently truncates because the board is Done-first sorted:

```bash
python scripts/board_open_items.py --role <role> --status Backlog --json
python scripts/board_open_items.py --role <role> --status Ready   --json
```

Keep only **XS / S** items (a quick win is small by definition; an M can qualify only if it's genuinely mechanical).

**Prefer the Ready queue; treat a Backlog pull as a commitment act.** Routine warm-up / "what's next" pulls come from the committed pool and never from Backlog ([[morning-routine-shared-backbone]] pull discipline) - stocking Ready with junk to feed a warm-up distorts priority. A burn-down is the *sanctioned exception*: the user has explicitly asked to clear quick wins, which authorizes committing Backlog items. But each Backlog pull is still a commitment, so apply Definition-of-Ready judgment before you pull one, and never stuff low-value work just to keep the loop fed.

Then **ratify each candidate before you branch** - do not trust the board's Status/Priority, which record where the Issue was *triaged*, not whether it's *currently actionable*. This is the canonical `[[issue-freshness-check-before-start]]` gate, not a burn-down invention; run its full 7-check pass (`gh issue view <N> --json body,state,comments,milestone,labels` + a look at linked/parent Issues). The checks and the drop/keep they produce:

| Freshness check | Fires when | Do | Example seen |
|---|---|---|---|
| 1 blockers | `blockedBy` / body blocker still open | drop (blocked) | #725 (blocked on #719); a *gated* AC "sequenced behind" unstarted work is the same shape - #842 |
| 2 supersession | a newer Issue overlaps / "supersedes #N" | drop, close-as-superseded | - |
| 3 stale/contradictory ACs | ACs appended post-filing conflict with originals (e.g. a tool switch) | drop or strip the dead ACs | #378 (HISAT2 vs STAR) |
| 4 parent/linked state | parent already closed COMPLETED, orphaning a verify task | drop (moot) | - |
| 5 priority inversion | a superseding Issue sits at *lower* priority | flag as hygiene, reconcile | - |
| 6 premise/target exists | the file/path/tool/infra it operates on is gone | drop (stale-premise) - `ls`/`git log` the target, don't trust the body | #183 (GCS; GCP decommissioned), #445 (retired `news_log.md`) |
| 7 already in progress | a parked branch / open PR / prior note has it 80-90% done | **resume, don't re-open greenfield** - grep `post-it.md`, `git worktree list`, `gh pr list --search <N>` | - |

Two burn-down-specific drops the freshness checks don't name, but that still disqualify a *quick win*:

- **not-a-PR** - verification-only or research-only; no code artifact to park at In review (#800, needs a live event that hasn't happened). A pipeline **run** ("execute the pipeline") is the same: it produces results, not a mergeable diff (#193).
- **not actually quick** - well-specified but wide (100+ occurrences, touches a cross-script contract, or collides with an open PR). #677 (~142-occurrence rename) is the archetype; it is a real task, just not a *quick* one, so it belongs in Step 3's stop condition, not the loop.

What survives is the real quick-win field. Rank: Ready before Backlog, then `arc-phase:active` first, then smallest, then freshest. This ranked field is the input to Step 1.5.

## Step 1.5 - Pre-flight triage (front-load the human-input decisions)

Before starting *any* item, read **all** of Step 1's surviving candidates' bodies in one pass and bucket each. This is the Definition-of-Ready / refinement judgment pulled to the **front** of the burn-down instead of discovered item-by-item mid-loop - the whole point is that a decision knowable from the issue body gets made once, up front, not when the loop stalls on it.

- **Clean-go** - passes freshness and the approach is unambiguous from the body. Runs autonomously in the loop.
- **Needs-one-answer** - a single decision blocks it *and that decision is answerable straight from the body* (a go/no-go, a which-of-two-approaches, a scope call). Do not start it silently, and do not drop it - it just needs one answer.
- **Drop** - genuinely needs design, is too wide, or is blocked. Name it with a one-line reason; don't debate it here (that is Step 3's boundary-case surfacing).

Then **surface every needs-one-answer question together, as one upfront batch**, and let the user clear them in a single pass - an `AskUserQuestion` batch when they are mutually-exclusive picks, a short numbered list when a freeform answer is likely (a batched `AskUserQuestion` loses all typed text if the user closes without answering every tab, per [[ask-user-question-batching]]). Fold the answered items into clean-go. The loop (Steps 2-3) then runs end-to-end through clean-go + cleared items with **no mid-loop input stops** for this knowable-upfront class.

**Caveat - only the knowable-upfront class.** This front-loads decisions readable from the issue body. Some problems only surface once you are inside the code (an API that misbehaves, a hidden coupling), so pre-flight triage **shrinks** mid-loop interruptions, it does not eliminate them. That residual is expected: handle it with the loop's normal freshness/stop logic, not as a triage failure. Establishing case (2026-07-04): a burn-down ran #784 clean to the gate, then stalled at #1002, whose marker-write approach (a memory-instructed stamp vs a hook signal) was decidable straight from the body - exactly the question an upfront batch front-loads.

## Step 2 - Bring one item to the merge gate

This is the per-item pipeline, and it is exactly the `[[autonomy-merge-gate-cadence]]` quick-win flow - run the whole reversible chain unsupervised and **stop at the merge command**, the one irreversible outward-facing act, which stays with the human. Don't over-checkpoint: stopping at PR-open to hand back the review or lab-notebook as separate asks is stopping *too early*, inside the autonomous zone.

1. **Start clean + branch.** `git checkout main` for a predictable starting point - do **not** `git pull origin main` here: `scripts/new_branch.sh` bases the branch off *remote* main server-side (`gh issue develop --base main`), so a local pull is redundant, and in the fan-out loop (Step 3) it would merge main *into* the previous item's still-checked-out feature branch. Then `scripts/new_branch.sh <issue#> <short-slug>` (never hand-roll the branch name; the helper preserves the Issue<->branch link and refuses epics).
2. **Move the Issue to In progress** on board #9 before opening the branch (the pre-PR step the auto-hooks don't cover): `scripts/pm/set_status.sh <issue#> "In progress"` (resolves the item id itself, idempotent). Raw-graphql fallback + IDs in the appendix.
3. **Build + verify.** Make the change; run the real verification, not just the tests a dry-run would pass. For CI/YAML or docs changes there may be no unit test - drive the actual behavior instead (run the validator green *and* red; build the venv and run pytest against it). Restore any file you mutate for a red-path test from a backup copy, never `git checkout` (it wipes uncommitted work).
4. **Commit, push, PR.** One logical change; PR body carries a ticked Test plan; tick the Issue's `## Acceptance criteria` boxes once truly satisfied.
5. **Request the bot review** (`@-claude review`, hyphenated only to dodge the mention guard - the literal trigger is `@claude review`) and **launch the `awaiting-bot-review` skill** so the wait is hands-free. Do not hand-roll a poll loop.
6. **Address the review**, reply with the fix SHA, then write the lab-notebook entry (it goes *after* review, before merge, per the closure ritual).
7. **STOP.** Do not run `gh pr merge` or `audit_and_merge.sh`. Print the exact merge command for the user and leave the card at In review. Requesting the review already flipped it there via the post-review hook.

The card now sits at In review, closure-ritual-clean, one command from done.

## Step 3 - Loop, and stop deliberately

Go back to Step 1's surviving field and take the next one. Two items in flight at once is fine and good - start the next build while the previous review bakes (see fan-out below).

**The stop condition is the whole point.** Stop the autonomous loop and hand back to the user the moment the next-best remaining candidate is not a clean quick win - i.e. it's only reachable by relaxing "quick" (a wide refactor), it needs a decision only the user can make (a go/no-go, a rescope), or it collides with open work. Do **not** stuff marginal or stale items into the loop to keep it running; an idle quick-win field is the correct, honest end state, and the board's own philosophy is to prune the option pool, not to grind it. When you stop, report: what's parked at In review (with merge commands), and the one boundary-case item with your read on why it's not in scope - as a question, not an action.

## Gotchas (all seen live)

- **Lab-notebook merge conflict on parallel same-role PRs.** Every PR inserts a `### <time>` entry at the same top-of-today anchor in `research/lab_notebook/<role>.md`, so the *second* PR the user merges will conflict there. Flag this up front in your closing report; offer to rebase the second after the first lands. It's trivial to resolve (keep both entries) but surprising if unannounced.
- **Bot-review bake-time.** A review takes ~5-7 min. Don't idle-poll: use `awaiting-bot-review`, and meanwhile start the next item's build. Never chain `sleep`s in Bash (the harness blocks it and the 2-min shell cap kills long polls anyway).
- **Branch-switching shows stale file contents.** Hopping between two in-flight branches makes the working tree show each branch's version; a "file modified since read" error usually just means you switched branches, not that work was lost - your commits are safe on their own branch. Prefer finishing one item before switching, or use worktrees for true parallelism.
- **Board Status field IDs must be queried, not guessed** (see appendix for the current ones).
- **`CLAUDE.md` is a symlink to `AGENTS.md`** - edit the real target; a direct write to the symlink is refused.

## Appendix - board #9 IDs (verified 2026-07-04)

Query these fresh if a mutation 404s (`updateProjectV2Field` regenerates option IDs):

- Project: `PVT_kwHOB17eGc4BSomP` (user project #9, `Jin-HoMLee`)
- Status field: `PVTSSF_lAHOB17eGc4BSomPzhAHFf8`
- Options: Backlog `f75ad846` · Ready `61e4505c` · In progress `47fc9ee4` · Ready for review `8bf9192f` · In review `df73e18b` · Done `98236657` · Epic `9f872564`

Set an Issue's status - **prefer the wrapper**, which resolves the project item id itself, maps the status name from the single canonical place, and is idempotent:

```bash
scripts/pm/set_status.sh <issue#> "In progress"
```

Raw-graphql fallback (get the project item id from `issue(number:N){ projectItems }`), only if the wrapper is unavailable:

```bash
gh api graphql -f query='mutation($item:ID!){ updateProjectV2ItemFieldValue(input:{
  projectId:"PVT_kwHOB17eGc4BSomP", itemId:$item,
  fieldId:"PVTSSF_lAHOB17eGc4BSomPzhAHFf8",
  value:{singleSelectOptionId:"47fc9ee4"}}){ projectV2Item{ id } } }' -f item="<ITEM_ID>"
```
