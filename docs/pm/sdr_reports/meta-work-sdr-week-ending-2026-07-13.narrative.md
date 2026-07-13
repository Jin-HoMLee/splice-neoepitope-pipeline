<!-- Author-owned narrative for Meta-work SDR - week ending 2026-07-13. Sections 3/4/5 only.
     The script regenerates the HTML from this file + fresh board data;
     it never overwrites this sidecar once it exists. -->

## Deliverables (Review layer)

**25 delivered, evenly split across the three lanes** (PM 10, Developer 9, Scientist 8).
The full list is in the Inventory appendix; the week has a single dominant theme and it is worth naming rather than enumerating.

**The theme: mechanisms that shipped inert, and the checks that could not see it.**
This week we did not mostly build new capability.
We mostly discovered that capability we believed we had was not running, and then built the ability to notice.

- **The auto-review hook had never fired.** [PR #1124](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/1124) shipped the first-pass bot-review auto-request ([Issue #1073](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1073)), the highest-leverage rung of [epic #1072](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1072), because human review bandwidth is our binding throughput constraint. It then emerged that its matcher never matched a heredoc-created PR, which is how essentially every PR here is opened, so the mechanism was inert from birth ([PR #1131](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/1131) / [Issue #1130](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1130)). Developer's own lab-notebook title for it is the honest one: *the hook I shipped never ran.* It fired for the first time on 2026-07-13.
- **A hook mirrored the Issue but not the PR** ([Issue #1108](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1108)), stranding PRs in the wrong review column.
- **The SDR's own trend was wrong, and so was the reference I filed against it** ([Issue #1099](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1099) / [PR #1128](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/1128)). The buggy script was closer to the truth than the yardstick I measured it with.
- **This report gained an arrival axis** ([Issue #811](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/811) / [PR #1143](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/1143)), splitting throughput into committed versus unplanned, and corrected the [#902](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/902) facet-3 claim that priority already carried that signal. It does not: urgency and arrival are orthogonal.
- **The coordinator autonomy envelope closed** ([Issue #1074](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1074) / [PR #1087](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/1087)), giving us the 3-tier reversibility ladder that gates the epic-1072 dispatch rungs.
- On the science side, the **NH-uniqueness prefilter** landed opt-in ([Issue #919](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/919) / [PR #1113](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/1113)) and surfaced that the **STAR path already pays the recall cost silently today**, which reframed the aligner fork ([Issue #1112](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1112)) from a preference into a correctness question. **pVACtools** was confirmed MHCflurry-only, class I, zero NetMHC ([Issue #1048](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1048)), and the **two-resolution registry design** landed ([Issue #1084](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1084)).

## Descoped (closed NOT_PLANNED)

Two, both legitimate supersessions with the need routed onward.
Neither is silently dropped scope.

- [#1071](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1071) fix(scientist): harden #680 registry derive scripts. **Why descoped:** superseded, folded whole into [#1076](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1076) (registry tooling hardening), where both items are now acceptance criteria alongside the single-source-of-truth refactor. No work was lost; no work had been started.
- [#451](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/451) Consider unifying `workflow/tests/.venv` + `research/.venv` into a repo-root `.venv`. **Why descoped:** the project settled the **opposite** and wrote it down normatively (`CLAUDE.md`, "Python environments": *no root-level `.venv` exists or should exist*). Closed by the Friday stale-Issue self-review, which is the cadence working exactly as intended: a stale option was actively pruned rather than left to age forever.

**Worth noting about both:** the closure-audit bot posted false-positive gap comments on each. That is a known defect ([#1137](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1137), n=3), not a real closure-ritual violation. A NOT_PLANNED close ships no work, so its acceptance criteria are *supposed* to stay unticked, and the bot cannot currently tell that apart from a genuine gap.

## Carried-forward & routing

None. Every issue opened in the window that was going to close, closed.

## Retrospective (process/health)

**Flow is healthy. Throughput is not the constraint, and the numbers say so plainly.**

25 delivered this week against a 4-week trend of **3 / 12 / 70 / 25**.
The 70 is the 2026-06-30 governance-sweep burst and is not a rate; roughly **25 per week is the sustainable line**, and we hit it.
Work in progress is nearly empty across every lane (PM 0, Scientist 1, Developer 0, Memory Manager 0), Ready sits at 17 against a cap of 23, and the review column holds 6 PRs against a limit of 10.
Nothing is jammed, nothing is starved, and no WIP guard is tripped.

**The one number that is drifting is median cycle time: 0.0 -> 1.8 -> 2.0 -> 2.9 days.**
It is still short in absolute terms, so this is a watch item and not yet a finding.
The plausible benign explanation is composition: as the trivially-fast governance items get pruned out of the pool, what remains is simply bigger work.
The plausible malign one is that review is starting to queue.
**Next week's number decides which**, and this is the first SDR positioned to tell them apart, because it is the first one riding an exact trend series ([#1099](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1099)).

**The arrival split is UNUSABLE this week, and I want that stated louder than the number itself.**

The report says **0% unplanned**, in every week, including the burst.
That is **not a measurement.** The `unplanned` label has been applied to **zero issues repo-wide**, ever, because it shipped four hours ago.
So the honest reading of that row is **"no data"**, and the report currently **cannot tell the difference** between a genuine 0% and a marker nobody has used.
It guards the empty-window case (which correctly returns "no data") and does not guard this one, so an unused marker prints a confident, clean, actionable `0.0%`.

That is the exact failure this metric was built to avoid, reproduced one layer up, in the code I merged today.
It is [#1144](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1144), and this week is its first live instance.
**Do not read the arrival axis, do not retune WIP against it, until #1144 lands.**

**The real retrospective finding, and it is a pattern, not an incident.**

Four separate things this month were **shipped, tested, green, and inert**: the auto-review hook ([#1130](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1130)), the review-column mirror ([#1108](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1108)), the SDR trend ([#1099](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1099)), and now the arrival marker ([#1144](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1144)).
In every case the unit tests passed, because the unit tests exercised the function and **not the trigger**.
A mechanism that never fires has no failing test and no surviving mutant; it produces a clean log, and a clean log reads exactly like good behavior.

This is now well enough understood to be structural rather than anecdotal, and it is being answered rather than admired: [epic #1135](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1135) carries it, with **[#1140](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1140) (hook-liveness contract tests: drive every hook's real trigger)** as the direct answer and [#1141](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1141) (mutation testing) as the complement.
The sequencing matters: **mutation testing is structurally blind to the biggest damage class here**, so #1140 is the one that pays.

**A second pattern, filed but worth naming here: our governance instruments do not implement our governance model.**
Three now, independently found: the SDR calls lead time "cycle time" because it cannot see the `Backlog -> Ready` commitment act ([#1138](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1138)); the dispatch digest under-counts multi-role Ready items and invents phantom shortfalls ([#1139](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1139), which came within one beat of triggering a false replenishment); and the arrival marker cannot see that same commitment act ([#1144](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1144)).
The board has a commitment point, and **not one of our tools can read it.**
A fourth instance stops this being three bugs and makes it an architectural finding about how the PM tooling reads the board.

**Actions out of this review:**

1. **Do not retune WIP limits this week.** The input the retune is supposed to consume ([#902](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/902)) does not exist yet. Retuning against a fabricated 0% would be worse than not retuning.
2. **[#1140](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1140) is the priority mechanism**, ahead of #1141.
3. **Watch median cycle time.** One more week of climb makes it a finding.
4. **Scientist Ready is at 4 of 5** after two honest decommits ([#601](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/601), [#817](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/817)) removed unstartable work from the shelf. Restock at the next replenishment; the shortfall is real but the number is finally truthful.
