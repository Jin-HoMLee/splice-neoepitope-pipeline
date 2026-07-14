# ADR-0003: Converge on STAR as the single aligner, run locally via a linux-64 container

- **Status:** Accepted
- **Date:** 2026-07-14 (decision date; ratified on merge of the PR closing [Issue #1112](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1112))
- **Deciders:** Jin-Ho (steer), Scientist (science sign-off), Developer (runnability evidence), PM (commitment)

## Context

The pipeline supports two aligners (`alignment.aligner: star | hisat2`) and has treated them as interchangeable backends.
A 2026-07-10 investigation spun out of [Issue #919](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/919) showed that they are **not** interchangeable, and that the stated reason for keeping HISAT2 as the local aligner does not survive measurement.

**Finding 1 - the two paths already disagree about what a junction *is*.**
The HISAT2 path (`bed12_to_junctions.py`, via regtools) counts **all** reads crossing a junction, multimappers included.
The STAR path (`star_sj_to_junctions.py:112`) reads `SJ.out.tab` **column 7 only** (uniquely-mapping reads) and never column 8, so a junction supported only by multimappers is discarded outright.
Flipping `alignment.aligner` therefore quietly changes which candidates reach neoepitope prediction.

**Finding 2 - "STAR needs more than 8 GB, unusable locally" is false at chromosome scale.**
Measured on the M1 (8 GB): a chr22 STAR index build peaks at **730 MB in 14 s**; aligning 500k reads peaks at ~230-380 MB.
The ~30 GB figure is the whole-genome suffix array, not chr22. RAM was never the local blocker.

**Finding 3 - the real local blocker is a broken bioconda macOS STAR build.**
Both the native `osx-arm64` and the Rosetta `osx-64` builds install cleanly and build the index correctly, but on alignment they report **0 input reads** for every FASTQ, including a synthetic read cut from the chr22 reference itself (`end of input stream, nextChar=-1` on the first byte).
Index generation works; the read parser returns instant EOF.
Reproduces on two independent architecture builds and is input-independent, so it is the build, not our data or invocation.
Secondarily, `workflow/envs/star.yaml` pins `star=2.7.10b`, which has **no osx-arm64 build at all** (bioconda's only arm64 STAR is 2.7.11b), so the STAR env is unsolvable on Apple Silicon as written.

Together these remove the RAM argument for the two-aligner split and raise the cost of maintaining it.

## Decision

**Converge on STAR as the single aligner, and run it locally via a `linux-64` container.**

- **STAR is the defensible convergence target on the science axis.** It is the field-standard aligner a neoantigen-pipeline reviewer expects, and its output carries the uniqueness distinction natively (`SJ.out.tab` cols 7 and 8).
  HISAT2's one remaining virtue is a native laptop binary that works, which is a developer-convenience argument, not a science one.
- **The container is the build that works, and it is identical to production.** It eliminates the two-path divergence rather than documenting it, and costs no new infrastructure concept: we already run TCRdock this way ([ADR-0001](0001-tcrdock-via-docker.md)).
  It also makes the `star=2.7.10b` arm64-unsolvable pin moot, since the pin resolves on `linux-64`.
- RAM stays trivial at chromosome scale (Finding 2), so the local story never depends on the RAM claim again.

**Binding condition on the science sign-off (Scientist, 2026-07-14).**
The Scientist signed off on **the aligner choice only** and explicitly **declined to sign off on the uniqueness policy**, on the grounds that no one has yet characterized the multimapped spliced reads (n=338 on the chr22 tumor) and that asserting "unique-reads-only is correct" today would be an *a priori* claim dressed as a verdict.

That distinction is load-bearing here, because **converging on STAR would otherwise decide the uniqueness question by construction**: with `star_sj_to_junctions.py:112` reading col 7 only, adopting STAR as the sole aligner silently adopts unique-reads-only as pipeline policy - never by a decision, only as an artifact of which column a script happens to read.
Note this is *our script's* choice, not STAR's: col 8 carries the multi-mapping counts and is right there.

Therefore, as a condition of this ADR:

- The col-7-only behavior **must become an explicit, configurable policy**, not an implementation detail. Carried by [Issue #1118](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1118).
- Its **default is set by [Issue #1122](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1122)'s verdict** (what the multimapped reads actually are), not by this ADR.
- **Aligner convergence must not be allowed to launder an undecided scientific question into a settled one.**

## Consequences

- HISAT2 stops being a co-equal backend. Its removal is not immediate: the containerized STAR path must first run the chr22 test end-to-end locally ([Issue #1162](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1162)), which is the gate on retiring the second path.
- Local development gains a container dependency it did not have. This is the accepted cost of a single code path that is identical to production; it is the same trade already made for TCRdock.
- The uniqueness semantics stay **explicitly open** until [Issue #1122](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1122) reports and [Issue #1118](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1118) lands the knob. Until then, a reader must not treat "STAR is the aligner" as implying "multimapper-supported junctions are correctly dropped".
- The CLAUDE.md claim that STAR is unusable locally is corrected to distinguish whole-genome (RAM-bound) from chromosome-scale (fine), naming the broken macOS build as the actual blocker.

## Alternatives considered

- **B. Converge on STAR, compiled from source for arm64.**
  Rejected: a native binary, but we own the maintenance of a compile toolchain and must still validate the arm64 build numerically against the container. The container gets the same result with no compile to own.

- **C. Keep both aligners, make the divergence explicit.**
  Rejected: it preserves two definitions of "a junction exists" in a pipeline whose output is a therapeutic target list, and pays maintenance for both. The status quo is not the safe option here; it is the incoherent one.

- **D. Converge on HISAT2.**
  Rejected, and actively argued against by the Scientist: it trades the field-standard aligner for local convenience on a pipeline whose credibility rests on methodological orthodoxy.
