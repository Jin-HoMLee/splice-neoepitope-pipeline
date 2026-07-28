<!-- Author-owned narrative for Meta-work SDR - week ending 2026-07-28. Sections 3/4/5 only.
     The script regenerates the HTML from this file + fresh board data;
     it never overwrites this sidecar once it exists. -->

## Deliverables (Review layer)

One descope, everything else delivered (counts in the summary block above - this narrative deliberately does not restate them, since the report recomputes from live board data on every render and a hard-coded figure here would silently drift out of agreement with it).

The week had an unusually clear spine: **almost everything that shipped was about making a silent failure loud.** That was not planned as a theme, and it is worth naming because it says something about where the system currently hurts.

**The commitment-hygiene chain (PM).** The week's largest piece closed a line of work that had been running for a fortnight. [Issue #1248](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1248) taught the Ready-floor guard to exclude body-gated Issues from the pullable count, and [Issue #1294](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1294) then consolidated four scattered "can this be worked?" checks into a single predicate over natively-owned sources (GitHub `blockedBy`, `needs-design` / `trigger-gated` labels, and the `Start date` field), with one enumerated reason taxonomy and the prose body-scan demoted to a deliberately low-recall proposer. The standalone bug fix [Issue #1299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1299) was folded into that work rather than shipped on its own, which is the descope recorded below.

The pattern under all three: gates that lived in prose, in a private post-it, or in a human's memory were moved into structured state a predicate can read. The failure being fixed was never "the gate was wrong", it was "the gate was invisible to every automated check".

**Board and tooling correctness (Developer).** Three of the five Developer items are the same shape. [Issue #1151](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1151) replaced a board mutation that could only ever report success with an assert-before-write resolver. [Issue #1242](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1242) added a lint for the zsh word-split anti-pattern, where a loop iterates once over an entire blob and silently does nothing. [Issue #1221](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1221) removed a recurring manual merge-conflict step by giving the lab notebooks a `merge=union` driver. Each targets a defect that previously presented as a clean run.

**One genuine science fix.** [Issue #1278](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1278) corrected minus-strand junction contigs that were swapping their exon arms and fabricating the splice boundary. This is the only item this week that changed pipeline output rather than process, and it is a correctness bug that produced plausible-looking wrong answers.

**Research and decision deliverables (Scientist, PM).** [Issue #1284](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1284) settled whether the literature-mining route for functional splice positives is exhausted, and concluded it is, redirecting effort toward the negative class. [Issue #601](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/601) delivered the AF3-class structure-backend re-evaluation verdict. On the PM side, [Issue #1280](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1280) (competitive positioning) and [Issue #265](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/265) (multi-agent autonomy SOTA) both closed as research artifacts, and [Issue #250](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/250) closed out the decision-telemetry mechanism with its first real signal.

## Descoped (closed NOT_PLANNED)

- [Issue #1299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1299) `fix(pm): not_pullable misses a trigger marker inside a checkbox AC` - **descoped because the fix became the wrong shape, not because the bug was unreal.** The bug reproduces exactly as filed. A best-practice cross-check on 2026-07-23 reframed [Issue #1294](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1294) from "route each cause to its native home" to "one predicate over natively-owned sources, with the body scan demoted to a low-recall proposer". Under that reframe, sharpening a regex in the demoted scanner is work against the grain of the decision. **Routed:** the behavior is now covered by the [Issue #1294](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1294) predicate reading a structured label, and the case survives as a regression fixture rather than a shipped fix. No residual need.

## Carried-forward & routing

**Nothing carried forward this week.** Every item that entered the window closed within it, so there is no carry to route.

This is worth one sentence of caution rather than celebration: a zero-carry week is as easily a symptom of low intake as of high completion, and this week's intake was low (see the throughput note below). It is not evidence of health on its own.

## Retrospective (process/health)

Five findings. Four are actionable and carry a named board carrier; one is explicitly an observation with no carrier, recorded as such rather than quietly dropped.

**1. The closure ritual holds in this repo and does not hold in the personas repo.** The 2026-07-28 closure audit scanned all 31 Issues closed across both repos since 2026-07-21. Five closed `COMPLETED` with unticked acceptance criteria: two here, three in personas. The split matters - the project repo's `scripts/audit_and_merge.sh` blocks a merge on unticked boxes, and the personas merge path appears to have no equivalent, which would explain the asymmetry better than discipline does.
**Carrier:** [Discussion #1311](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/discussions/1311) (Memory Manager, asking whether to tick or whether personas deliberately runs a lighter ritual and my audit should be scoped out).

**2. An unticked box was hiding a real defect, which is the argument for the ritual.** [Issue #1156](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1156) closed `COMPLETED` with 5/5 criteria unticked and no closing PR. Verifying each against live files found four genuinely met and one not: a full paginated sweep showed **2 of 18 open parents still carrying a Priority** the convention forbids. Both were cleared and re-read to confirm. The criterion that failed is the one that explicitly said *"verified by reading the board field back, not by trusting the mutation's exit code"* - so the instruction was correct and simply was not followed, on the very Issue that wrote it down.
**Carrier:** fixed inline this session; verification recorded on the Issue. No further carrier needed.

**3. `check_milestone_health.sh` is structurally blind to an entire class.** Line 50 filters on `select(.due_on != null)`, so a milestone with no due date is dropped from the input set and can never be reported however complete it is. The live instance, `i6 - S5 - Modeling`, had been finished and invisible; it surfaced only because I listed milestones directly while investigating a different one. The check has been passing cleanly every day while unable to see the class.
**Carrier:** [Issue #1312](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1312).

**4. This report is 11 days late, and the reason is structural.** The Friday-cleanup beat carries three weekly obligations - branch hygiene, the per-role stale-Issue review, and this SDR - and it hangs off the morning routine. Last Friday's session opened on a resume greeting, correctly ran the light resume routine, and dropped all three silently. A weekly obligation whose trigger is a property of one session's opening message is not really weekly.
**Carrier:** [Issue #1314](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1314).

**5. Throughput is down sharply and I am declining to call it a problem - observation only, no carrier.** See the weekly trend table above: delivered-per-week has fallen by roughly three quarters from the early-July peak, while median cycle time has risen from under a day to around three. Read naively that is a collapse. I do not think the data supports that reading yet, for three reasons: the earlier windows include large quick-win burn-down batches which inflate counts with small items; this window's last three days had zero closures because they were not worked; and the composition shifted toward fewer, larger items ([Issue #1294](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1294) alone was a multi-day build), which raises median cycle time without indicating a flow problem. Declaring a capacity trend from four points across a known composition change would be over-reading. **Recorded as an observation to re-check next week** rather than routed to a carrier; if the next window also lands near 14 with a rising median, that is a third point on a consistent composition and worth an Issue then.

## Routing summary

Four actionable findings, four named carriers ([Discussion #1311](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/discussions/1311), inline fix, [Issue #1312](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1312), [Issue #1314](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1314)). One observation explicitly marked `observation only, no carrier` with the reasoning stated. Nothing in this report is left un-routed.
