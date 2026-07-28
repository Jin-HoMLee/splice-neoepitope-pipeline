# Pullability predicate - design (Issue #1294, PR1)

**Status:** approved 2026-07-24 (brainstorm with Jin-Ho).
**Issue:** [#1294](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1294) - "one pullability predicate over natively-owned sources, replacing four independent detectors".
**Scope of this spec:** the predicate-first first PR only. Two pieces are deliberately deferred to sequenced follow-ups (see the Deferred section).

## Problem

Whether a Ready-queue Issue can actually be worked ("pulled") is decided today in several independent places.
The Backlog-to-Ready floor guard (`scripts/check_ready_queue.sh`) counts Ready *depth*, and consumes a `not_pullable` reason computed by a prose scanner (`scripts/pm/not_pullable.py`).
Native `blockedBy` state is checked elsewhere again (`scripts/pm/scan_unblocked.py`, `scripts/pm/scan_prose_deps.py`).
Date gates and external-event triggers have no structured home at all - they live only in Issue prose and in a role's private post-it.

Scattered coverage of one rule is the **shotgun surgery** smell: one reason to change ("what makes an Issue un-pullable?") forces edits across many files, and the characteristic failure is *missing one*, producing a silent inconsistency.
The scanner's first real use (the 2026-07-23 Replenishment pass over 17 Backlog candidates) produced three distinct misses in bodies written by the same author who wrote the patterns, because English has unbounded ways to phrase "we have not decided yet" and a regex set does not converge.

## What is NOT the fix

**Not one big field.**
Single Source of Truth means each data *element* is mastered in exactly one place, not that everything lives in one place.
A dependency is genuinely owned by GitHub's `blockedBy`; a date is genuinely a date.
Copying those into one custom "blocked reason" field would create duplicates that drift - trading a coverage bug for a silent consistency bug, which is worse.
So the sources stay plural.

**Not more regexes.**
Growing the prose pattern set on every near-miss is the denylist widen that never terminates ([Issue #1150](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1150) posture).

## The fix: consolidate the decision, not the state

One **pullability predicate**: a single function answering "can this Issue be pulled now?", consulting each natively-owned source, returning a reason from one enumerated taxonomy, or `None` if pullable.
Sources stay where they are mastered; the *decision* gets exactly one home.
Everything that used to decide for itself (`check_ready_queue.sh`, the floor guard) stops deciding and just asks the predicate.

## Architecture

### The predicate module - `scripts/pm/pullability.py`

A pure function over an Issue's already-fetched structured fields.
No network and no stored state, mirroring the existing `not_pullable.py` discipline - any run re-derives the answer, and it is unit-testable in isolation.

```
assess(item) -> Optional[Reason]
```

`item` carries the fields `board_open_items.py` already fetches for each card: its label list, its open `blockedBy` edges, and its `Start date` field value.

The `Reason` values are one enumerated taxonomy (a small set of module-level constants), each mapped to exactly one authoritative source, checked in a fixed order:

| Reason string | Sole authoritative source (read, never copied) |
|---|---|
| `blocked-by-issue: #N` | any *open* GitHub `blockedBy` edge on the Issue |
| `needs-design` | the `needs-design` label is present |
| `trigger-gated` | the `trigger-gated` label is present (new label, this PR) |
| `date-gated: <date>` | the `Start date` field is set to a future date |

`assess` returns the first matching reason, or `None`.
Adding a fifth cause later is one enum entry plus one source read, in this one file - the AC's "one edit in one place".

Ordering rationale: `blocked-by-issue` first because a hard dependency subsumes any softer gate; the label causes next; `date-gated` last because a date gate is the weakest ("not yet" rather than "not workable").

### Prose scan demoted to a proposer - `scripts/pm/not_pullable.py`

`not_pullable.py` is repurposed from an authority into a **low-recall convenience proposer**.
It no longer feeds the floor guard's decision.
At the Backlog-to-Ready transition it *suggests* a label for a human to apply:
a decision-shaped body suggests `needs-design`; a trigger-marker body suggests `trigger-gated`.
Its module docstring states explicitly that low recall is **intended, not a defect**: the label is the authority, so a missed phrasing is a non-event, and there is never a reason to grow its pattern set to chase completeness.

The four known misses become its regression fixtures, documenting *accepted* proposer behavior:
[Issue #1111](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1111) and [Issue #353](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/353) (noun-form decision ACs),
[Issue #1241](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1241) (checkbox-prefixed marker, folded from [Issue #1299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1299)),
[Issue #876](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/876) (no body gate at all).
A fixture that the proposer misses is asserted as an accepted miss, not a failure.

### Wiring - minimal blast radius

`scripts/board_open_items.py` computes the per-item `not_pullable` field from `pullability.assess(...)` instead of `scan_not_pullable(body)`.
The JSON **key name `not_pullable` is kept**, so `scripts/check_ready_queue.sh` needs no logic change - it already splits Ready into `not_pullable == null` (pullable) versus non-null (gated).
The reason string it prints now comes from the taxonomy.

`check_ready_queue.sh`'s cap comment is upgraded to cite the WIP-limit convention it matches - blocked items count against the WIP cap while the floor excludes them - rather than reading as a local choice (an AC of the Issue).

## One-time backfill (cutover parity)

The predicate only excludes an Issue once its label or date is actually set, and the proposer is deliberately low-recall.
So at cutover the structured markers must be **backfilled** onto the currently-known gated Issues, or the queue would regress: items the prose scan catches today by text would leak into the pullable count.

Backfill targets are the standing gated set (each **freshness-checked before labeling**, not bulk-applied):
[Issue #841](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/841) (`trigger-gated`),
[Issue #876](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/876) (`trigger-gated`),
[Issue #929](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/929) (`needs-design`),
[Issue #1247](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1247) (`Start date` ~2026-08-03),
[Issue #817](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/817) (`Start date` 2026-12-01),
[Issue #585](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/585) (condition-gated - re-read to classify).
The decision-shaped Backlog items [Issue #1111](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1111) and [Issue #353](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/353) get `needs-design` if the freshness read confirms the fork is still open.

## Testing

- `pullability.assess` unit tests: one matched-pair per reason (a gated fixture returns the reason; an identical-but-ungated fixture returns `None`), plus the ordering case (an Issue that is both blocked and dated returns `blocked-by-issue`).
- The proposer's four known-miss fixtures assert *accepted* behavior (miss recorded as accepted, hit recorded as a suggestion).
- A cross-module test that `board_open_items.py` emits `not_pullable` from the predicate (swap a fixture's label and the field flips), so the wiring cannot silently regress.
- Live smoke: run `check_ready_queue.sh` post-backfill and confirm the known-gated Issues are excluded from the pullable count and appear in the `[NOT-PULLABLE ...]` line - a matched-pair control (before backfill they leak in; after, they do not).

## Deferred to sequenced follow-ups

Per the predicate-first scope decision:

- **[Issue #876](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/876) as full structured Issue state** - writing its revival gate into the Issue body as structured state and removing the [Issue #1248](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1248) "permanent known miss" note.
  PR1 gives #876 a `trigger-gated` label (so the predicate excludes it), which is enough for the queue; the deeper "eliminate the unqueryable location" work is its own commit-point.
- **The DoR-text edit** in personas `shared/feedback_board_hygiene.md` (an unresolved-decision AC is not Ready) - MM-routed, lands as an MM commit, tracked with [Issue #1260](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1260)'s cross-repo shape.
- **Any new board *field*** - not needed: date gates reuse the existing `Start date` DATE field.

## Notes and residual risks

- **[Issue #715](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/715) constraint is satisfied, not open.**
  Confirmed 2026-07-24: `Jin-HoMLee` is a User account, so `issueType` returns NONE and Issue Types are unavailable.
  That is *why* the structured surfaces are labels plus a Projects date field, not a Type field.
- **`Start date` reuse.**
  Reusing `Start date` for "earliest-pullable" needs no new field, but carries a minor risk of a future collision with a roadmap "Start" semantic.
  Accepted for now; if a roadmap use of `Start date` ever lands, a dedicated "Not before" date field is the clean split.
- **`trigger-gated` label creation** is a plain `gh label create`, not the destructive single-select-option regeneration path (the `updateProjectV2Field` options-are-destructive rule) - labels are independent GitHub objects, not Projects single-select options.
</content>
</invoke>
