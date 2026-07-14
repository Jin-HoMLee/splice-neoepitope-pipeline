# Issue #1122 - what ARE the multimapped spliced reads, and what should we do with them?

**Issue:** [#1122](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1122)
**Companion evidence base:** [`issue_919_nh_uniqueness_filter/`](../issue_919_nh_uniqueness_filter/) (Developer built the instrument and the A/B; this is the interpretation and the decision)
**Status:** complete.
**Decision:** the NH-uniqueness filter stays **default OFF** - and, as currently *applied*, it is **structurally disqualified** for a matched tumor/normal design at any index and any depth.

## The answer, in one paragraph

`NH` is a property of an individual **read** (how many places its sequence could go), not of a **junction**. Two libraries sampling the same biological junction draw different reads and therefore get different `NH` profiles, by chance. A junction-level gate built on `NH`, applied **independently to the tumor and the normal**, erodes the two arms independently - and a matched subtraction cannot survive having its two arms independently eroded. That is not a tunable weakness; it is a category error about what `NH` measures, and it is why the filter must not be turned on, rather than merely why the A/B was unconvincing.

## The demonstration

The **IGLJ3 -> IGLC3** junction (`chr22:22905026:22906341:+`) is the canonical J-to-C join in a rearranged immunoglobulin lambda light chain - a real, functional antibody transcript from plasma cells in the tissue. It is present in **both** the tumor and the matched normal.

| | supporting reads, filter OFF | filter ON |
|---|---|---|
| **tumor** | 4 (all `NH=1`) | **4** |
| **normal** | 5 (one `NH=1`, three `NH=3`, one `NH=2`) | **1** |

The filter leaves the tumor untouched and destroys the normal. Normal support falls below `min_normal_reads: 2`, the junction stops counting as "seen in normal", and its tumor copy is **promoted from `normal_shared` to `tumor_exclusive`**. A junction with five reads of matched-normal support becomes a tumor-specific neoepitope candidate *because a precision filter was switched on*.

The mechanism is visible at read level in the two BAMs: a read at `pos=22905015` with CIGAR `11M1316N39M` - the **same alignment of the same junction** - is `NH=1` in the tumor and `NH=2` in the normal.

**Scope, stated honestly.** Across the 91 junctions detected in *both* libraries, the filter changes the tumor:normal read ratio for **4**, and exactly **one** became a false `tumor_exclusive`. This is an **existence proof, not a rate** - 91 shared junctions is far too small to estimate a frequency, and n=1 is not a percentage. It is sufficient anyway, because the cost function is asymmetric: a false `tumor_exclusive` is a candidate *therapeutic target*. Notably, **two of the four asymmetries are immunoglobulin junctions** (IGLJ3-IGLC3, IGLC5-IGLC6) - a striking concentration, given how small a slice of chr22 the IG lambda locus is, and exactly the paralog-rich domain (IG/HLA/TCR) this pipeline targets.

## Why this settles AC-4 without the whole-genome run

Developer's load-bearing finding is correct: **chr22 cannot measure the filter's precision.** `NH` is index-relative, so the population the filter exists to catch - a read whose true locus is off-chr22, landing on chr22's single copy with `NH=1` - is invisible here. No refinement of the A/B fixes that, and [#1095](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1095) would be needed to measure how often the filter *helps*.

But we do not need to measure how often it helps. The defect above is **structural, not statistical**: it does not depend on the index, the chromosome, or the depth. A filter that corrupts the matched comparison is disqualified regardless of its hit rate. **AC-4 answer: yes, the fixture can settle this. #1095 is not a prerequisite.**

## The field default, and why its justification does not transfer

Unique-only **is** the field standard and it is not cargo cult. [LeafCutter's documentation](https://davidaknowles.github.io/leafcutter/articles/Usage.html) states its "most restrictive filter is the requirement that reads considered be uniquely mapped"; STAR institutionalises the same split in `SJ.out.tab` (col 7 unique / col 8 multi).

The default is calibrated for **quantification across replicates**: LeafCutter and the differential-splicing literature compare *ratios* within junction clusters across many samples, where per-read multimapping noise is applied to both conditions in expectation and **averaged over replicates**. A misplaced read costs a slightly wrong PSI.

Our estimand is different: a **single-pair binary subtraction**, one tumor and one normal, with a hard `min_normal_reads: 2` gate. Nothing is averaged, the noise does not cancel, and it **flips a binary call**. A misplaced read costs a false therapeutic target.

So: unique-only is right for its home problem and wrong for ours. The paralogy of HLA/IG/TCR determines *where* it fails; the **single-pair design** determines *that* it fails. This is a real differentiation point - the splice-neoepitope field inherits its tooling from splice *quantification*, and inherits filter defaults whose statistical justification does not survive the change of estimand. [Splicing neoepitope prediction is sensitive to methodological differences](https://doi.org/10.1101/2025.09.10.674685) (2025) already shows junction 9-mer sets diverging 69-85% across similar pipelines on parameter choices alone; this is a concrete mechanism for one such divergence.

## Two things found on the way

### The real depth gate is a silent floating mean -> [#1161](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1161)

`filter_junctions.py:400` keeps a tumor junction only if `reads > mean(reads)` for that file. On this fixture the mean is **1.205**, so it collapses to `reads >= 2` and removes **1,638 of 1,872 junctions (87.5%)** - versus the NH filter's 14%. The threshold is a function of the library's own depth, not a decision; on a deeper library it silently becomes `reads >= 4` or higher. It is also what makes the pipeline **non-monotonic**: because the mean is recomputed from whatever survives, an upstream filter that only *removes reads* can *add* a downstream candidate. Filed separately as [#1161](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1161).

This also corrects #919's headline. Its "16.2% of the tumor_exclusive candidate set (251/1,550)" is measured on the **raw** junction set. Production does not build candidates from the raw set - it mean-filters first, so the real candidate set is **151**, and the filter's real cost is **43/151 = 28.5%**. Measuring against the raw set understated the damage by nearly 2x, because 86% of what the NH filter removes is single-read and would have been dropped by the mean filter anyway.

### Two probes that do not work - recorded so they are not re-attempted

- **The splice-motif probe is void.** "Are the removed junctions less canonical?" looks like the ideal index-independent test. It cannot fire: **all 1,872 junctions are GT-AG** (999 `GT..AG` on `+`, 873 `CT..AC` on `-`, zero exceptions), because HISAT2's `XS` tag is motif-derived and `regtools -s XS` consumes it. The motif is **already an inclusion gate upstream**, so a motif-based probe restates the gate and can only return a null. (Checked separately: all 2,274 spliced alignments carry `XS`, so regtools discards nothing - the absence of GC-AG / AT-AC junctions is HISAT2's own, not a loss we introduce.) This does *not* sink [#1116](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1116): Portcullis-style **anchor-vs-intron Hamming** is a different signal and remains viable.
- **The NH filter does not clean up [#1147](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1147)'s low-complexity artifacts.** Scoring each junction's translated contig (24 nt flanking each side) for Shannon entropy: the filter removes **0 of 267** low-complexity contigs and leaves all 9 behind, so it would very slightly *concentrate* them. **#1122 is paralogy at the alignment stage; #1147 is low-complexity at the translation stage.** They are orthogonal and need separate filters.

## Developer's findings, independently checked (AC-5)

Re-derived from the committed `issue_919` outputs, the two chr22 BAMs, and an annotated-intron set built from the GENCODE chr22 GTF (independent of the BED Developer used):

- **The stratified null reproduces exactly.** Unstratified enrichment 1.13x; within the unannotated pool 0.98x. I recover his exact partition - lost/annotated 11, lost/unannotated 256, retained/annotated 276, retained/unannotated 1,329 - from a *different* annotation source. A genuine cross-validation, and the strongest single result in his record.
- **The mismatch table reconciles; only its header is wrong.** Rows sum to 2,277 against 2,274 spliced alignments. The unit is **junction instances**, not alignments: 15 reads carry two introns and are counted twice. Re-tallied per instance, `NH=1` = 1,939 (= 435 + 1,504, exactly) and `NH>1` = 350 (= 338 + 12 annotated).
- **Both retractions were correct.** `NH>1` is not a dirtier population than `NH=1`; the mismatch profiles are indistinguishable. **The multimapped spliced reads are largely real RNA, not misalignment junk** - which is AC-1's answer, and the opposite of the filter's founding premise. Concretely, on chr22 they are immunoglobulin J-C and C-C junctions, chr22-internal paralogs (PI4KA/PI4KAP2, the GGT family, SEC14L4/L6), and a bulk of single-read junctions statistically indistinguishable from the retained pool.
- **The "chr22 inverts the filter" argument is correct** and is the most valuable thing in his record.

## Consequence for #1112 / #1118

`star_sj_to_junctions.py:112` reads `SJ.out.tab` **column 7 only** (uniquely-mapped reads), discarding column 8 (multimappers). That **is** the unique-only semantic. So **converging on STAR silently adopts, by construction, the exact policy this Issue disqualifies** - never as a decision, only as an artifact of which column a script reads. Column 8 is right there. This must become an explicit, configurable policy, and on the evidence here its default should be **count all reads**.

## Reproducing

Everything below runs from **committed inputs only** (the `issue_919` A/B outputs plus `resources/test/chr22.gtf.gz`) - no BAM, no network:

```bash
research/.venv/bin/python \
  research/experiments/issue_1122_multimapped_reads/nh_matched_subtraction.py \
  --json research/experiments/issue_1122_multimapped_reads/outputs/results.json
```

The read-level `NH` evidence (the `11M1316N39M` read that is `NH=1` in tumor and `NH=2` in normal) needs the chr22 BAMs, which are gitignored. Regenerate them with the chr22 test pipeline, then:

```bash
samtools view results/patient_001_test/alignment/SRR9143066_test/SRR9143066_test.bam \
  chr22:22904000-22907000 | awk '$6 ~ /N/'
```

## Outputs

| file | what it is |
|---|---|
| `nh_matched_subtraction.py` | the analysis: the mean-filter profile, the end-to-end candidate cost, and the matched-arm asymmetry |
| `outputs/results.json` | machine-readable results, including the full asymmetric-junction list |

## Cross-experiment deps

- **Reads from:** [`issue_919_nh_uniqueness_filter/outputs/`](../issue_919_nh_uniqueness_filter/outputs/) (Developer's committed A/B junction sets) and `resources/test/chr22.gtf.gz`.
- **Feeds:** [#1118](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1118) (the aligner-coherence defect - this decides its default), [#1112](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1112) (the STAR col-7 trap), [#1116](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1116) (which probe is viable), [#1161](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1161) (the mean-reads gate).
