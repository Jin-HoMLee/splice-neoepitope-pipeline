# SRA-accession fetch to gzipped FASTQ for ENA-unmirrored cohorts

Design for [Issue #1296](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1296).
Author: Developer.
Date: 2026-07-28.

> **Revision note.** The first draft of this spec designed the whole tool from scratch and never checked prior art.
> A best-practice pass found [kingfisher](https://github.com/wwood/kingfisher-download), an actively-maintained tool that already implements the same fetch-and-convert chain.
> This version adopts it and keeps only the wrapper that is genuinely ours.
> The superseded from-scratch reasoning is retained in "Rejected alternatives" rather than deleted, because the constraint analysis there is what makes the adoption case legible.

## Problem

`scripts/prepare_production_data.sh` acquires FASTQs by hardcoded ENA HTTPS URL, and its own header states the assumption: "Files are downloaded via ENA HTTPS - no sra-tools or controlled access required."
That assumption fails for the Courcelles CRC cohort needed by the [Issue #1176](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1176) back-half.

Re-verified live on 2026-07-28, independently of the Issue's 2026-07-23 checks:

- ENA reports 26 runs for `PRJNA1372708` with **0** having `fastq_ftp` and **0** having `submitted_ftp`.
- The NCBI SDL API returns exactly one file type for these runs: `sra`. There are no original submitted FASTQs to shortcut to, so conversion is unavoidable.
- The public ODP mirror is reachable unauthenticated: `SRR36274715` is 6,120,422,876 bytes and `SRR36274703` is 6,475,790,169 bytes, both HTTP 200.
- Library layout is `PAIRED`, strategy `RNA-Seq`, platform `ILLUMINA`.

This is a tooling gap, not an access-permission problem.
The gap will recur, since the ENA-has-metadata-but-no-bytes case is not specific to this cohort.

## Prior art

The section the first draft was missing.

[**kingfisher**](https://github.com/wwood/kingfisher-download) procures sequence files from ENA, NCBI SRA, AWS and GCP, with an explicit fallback chain: try ENA `.fastq.gz`, else fetch `.sra` from the AWS Open Data Program and convert, else fall back to NCBI `prefetch`.
That is the same chain this Issue describes, including the ODP path our cohort needs.

Verified via `gh api` rather than from search summaries:

| Property | Value (checked 2026-07-28) |
|---|---|
| Last push | 2026-07-26 (two days before this spec) |
| Latest release | v0.5.0, 2026-04-10 |
| Stars / open issues | 312 / 10 |
| Archived | No |
| License | GPL-3.0 |
| Packaging | bioconda, **noarch**, v0.5.0 |

Two capabilities decide the matter:

- `--check-md5sums` is supported **specifically for the `aws-http` method**, which is the ODP path our accessions require.
- Its most recently merged PR is "pipe-sracat-into-pigz": streaming `.sra` extraction piped directly into gzip, with no uncompressed intermediate.

That second point eliminates the disk constraint that drove the entire first draft (see "Rejected alternatives").

**Licensing.** This repo is MIT; kingfisher is GPL-3.0.
Invoking a CLI as a subprocess does not create a derivative work, so there is no license interaction.
We must not vendor or link its source.

## Decisions

### D1. Adopt kingfisher as the engine; own only the wrapper

Build-versus-adopt resolves to adopt.
Kingfisher already implements the download, the source fallback chain, MD5 verification on our exact path, and streaming conversion.
Re-implementing those is work whose best possible outcome is parity with a tool that is maintained by someone else and was updated two days ago.

What stays ours is what kingfisher has no opinion about, because it is specific to this repo:

- Output naming and layout that a `config/samples/<patient>.tsv` samplesheet consumes unchanged.
- The provenance record, which doubles as the idempotency sentinel (D4).
- The skip-if-complete check (AC5).
- The preflight free-space check (D5).
- A short pointer note on `prepare_production_data.sh` (AC6).

This shrinks the surface we own to roughly the assertion-and-bookkeeping layer, which is precisely the part this Issue's acceptance criteria are actually about.

### D2. A thin shell launcher plus a tested Python core

Unchanged from the first draft, and the adoption in D1 strengthens it: with fetching delegated, what remains is almost entirely bookkeeping and assertion, which is data manipulation rather than utility-calling.

The [Google Shell Style Guide](https://google.github.io/styleguide/shellguide.html) sanctions shell when "mostly calling other utilities and doing relatively little data manipulation", and directs a rewrite past ~100 lines or non-straightforward control flow.
Orchestration (checking for a Docker daemon, starting colima) stays in shell.
Everything else is Python with pytest coverage.

Two project-local facts reinforce it: local shell scripts must be bash-3.2-safe on macOS, and [Issue #1043](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1043) exists because untested shell has accumulated here.

### D3. Pin kingfisher and the container by digest

Container images are pinned by **digest**, not tag.
A tag is a mutable pointer that can silently resolve to different bytes; a digest is content-addressed and cannot.
This is standard guidance for reproducible scientific workflows.

kingfisher is pinned to an explicit version (`0.5.0` at time of writing), not to `latest`.

The image is built from bioconda rather than pulling `ncbi/sra-tools` directly, since kingfisher brings its own extraction dependencies.
Note that had we used `ncbi/sra-tools` directly, Docker Hub's `latest` resolves to **3.4.1**, which `AGENTS.md` already records as segfaulting ("use version 3.1.1 on GCP VMs - newer versions (3.4.x) have a segfault bug").
Delegating the sra-tools dependency to kingfisher's own pinning removes that trap from our surface, and is a further argument for D1.

### D4. The provenance record is the idempotency sentinel

AC5 requires that an already-completed run is skipped, which needs a completion marker regardless.
A zero-byte `.done` file carries no information, so the marker instead *is* the provenance record.
Existence means complete; contents answer what produced the data.

Recorded fields: accession, resolved source method, upstream MD5, kingfisher version, image digest, exact command, UTC timestamp, output paths with sizes and read counts, and the subset cap if any.

Reproducibility guidance (for example [QIIME 2 Provenance Replay](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1011676)) is consistent that a checksum alone is necessary but not sufficient, and that tool versions, parameters and environment must also be recorded.

### D5. Verification is subset-first, with the full pull deferred

Streaming extraction removes the scratch-space problem but not the output-size problem.
Working from NCBI's sizing (final FASTQs approximately 7x the accession), each accession yields roughly 13 GiB gzipped, so both together land about **26 GiB against 39 GiB free**.
Peak during the second run, holding the first run's output plus the second `.sra` plus the second output, is near **32 GiB**.

So the preflight free-space check remains load-bearing, and dropping each `.sra` as soon as its conversion succeeds remains required rather than tidy.

The mechanism is verified on a bounded read subset of **both** real accessions, exercising every step in minutes and a few GiB.
The full-data pull runs separately and unattended.

### D6. Viability is proven before the wrapper is written

D1 rests on documentation and repository metadata, not on kingfisher having run here.
So the **first** implementation step is a bounded viability check: run kingfisher against `SRR36274715` in the linux-64 container, read-capped, and confirm it produces valid paired gzipped FASTQ via the ODP path with MD5 checking on.

If that fails, the fallback is the from-scratch design preserved in "Rejected alternatives", amended with the `fastq-dump` mitigations recorded there.
Committing to adoption before this check would be exactly the "verified by documentation" error this project keeps having to correct.

## Architecture

### `scripts/fetch_sra_fastq.sh`

A launcher, targeted well under 100 lines.

- Verify the `docker` CLI exists and the daemon is reachable, starting colima if not, reusing the pattern proven in `scripts/run_local_linux64.sh`.
- Verify the repo is under `$HOME`, since colima mounts only `$HOME` and `/tmp/colima`, and a bind mount of any other path silently yields an empty directory inside the container.
- Exec the Python core, forwarding arguments unchanged.

### `scripts/sra_fetch.py`

Tests in `scripts/tests/test_sra_fetch.py`, one of the five directories the CI pytest job enumerates.
This is `scripts/`, not `workflow/scripts/`: nothing here is invoked by a Snakemake rule.

```
already_complete(acc, outdir)   -> bool
preflight_disk(acc, outdir)     -> None   (raises if insufficient)
run_kingfisher(acc, …)          -> FetchResult
normalise_outputs(result, …)    -> OutputPaths   (samplesheet-consumable names)
assert_integrity(paths)         -> None   (raises on any failure)
write_provenance(...)           -> Path
```

## Data flow, per accession

1. **Skip check.** If `data/<ACC>.provenance.json` exists, report and return.
2. **Preflight disk.** Refuse to start if free space is insufficient, stating required versus available.
3. **Fetch and convert** by invoking kingfisher in the pinned container, with `--check-md5sums` and the ODP method, writing gzipped FASTQ.
4. **Normalise outputs** to `data/<ACC>_1.fastq.gz` and `data/<ACC>_2.fastq.gz`.
5. **Assert integrity**: `gzip -t` on each mate, and equal read counts between mates.
6. **Drop the `.sra`** and any scratch.
7. **Write the provenance sentinel last**, so it exists only if every prior step passed.

### Output layout

Flat in `data/`, matching `prepare_production_data.sh`'s convention, so the paths drop into a samplesheet's `fastq1` and `fastq2` columns with no renaming (AC4).
Scratch lives under `data/.sra_cache/`.
A `SINGLE`-layout accession emits `data/<ACC>.fastq.gz` and records the layout rather than asserting a pairing that never existed.

## Error handling

Every integrity check fails closed.
A download that cannot be verified is treated as corrupt, never as acceptable.
This follows the project's fail-safe-not-fail-open rule: a gate must never false-PASS.

| # | Failure | Behavior |
|---|---|---|
| E1 | Accession unknown, or no source has bytes | Fatal, accession named, worded so a genuinely-unavailable accession does not read as a bug |
| E2 | kingfisher MD5 check fails | Fatal, partial output removed so a re-run cannot resume onto a poisoned file |
| E3 | kingfisher exits non-zero | Fatal, its stderr surfaced verbatim rather than summarised |
| E4 | `docker` absent or daemon unreachable | Caught in the launcher before any download begins, with install guidance |
| E5 | Expected mate file missing after conversion | Fatal, naming which mate |
| E6 | Mate counts unequal | Fatal, **both counts printed** |
| E7 | Insufficient disk at preflight | Refuses to start, states required versus available |
| E8 | Extra or orphan output file | Warning only, recorded in provenance |

## Testing

### Unit tests, `scripts/tests/test_sra_fetch.py`

The load-bearing ones are matched pairs: identical inputs, one variable flipped, opposite expected outcomes.
A test that can only pass is not a check.

| Control | Falsifier |
|---|---|
| Equal mate counts pass | Unequal mate counts fail |
| Sentinel present skips without invoking kingfisher | Sentinel absent proceeds |
| Successful run writes a sentinel | Failed fetch writes **no** sentinel |
| Sufficient disk proceeds | Insufficient disk refuses |
| Both mates present pass | Missing mate fails, naming it |

kingfisher itself is stubbed at the subprocess boundary, so these stay fast and offline.
Per current Python guidance, HTTP-level interactions (where any remain) are stubbed with `responses` rather than hand-rolled mocks.

### Live subset smoke

Both real accessions, read-capped: real outputs, real `gzip -t`, real equal mate counts, real provenance written, and a real second invocation that skips.

Unit stubs validate our *model* of kingfisher, not kingfisher itself, so the live smoke is required by the project's live-integration-smoke rule for any change touching an external boundary.

## Acceptance criteria disposition

| AC | Disposition |
|---|---|
| 1 - fetch to gzipped FASTQ without host sra-tools | Met (sra-tools lives in the container, via kingfisher) |
| 2 - verified end to end on both accessions | **Mechanism half met** (subset, both accessions). Full-data half deferred to a carrier Issue |
| 3 - integrity asserted, not assumed | Met, and exceeded: upstream MD5 in addition to `gzip -t` and mate-count equality |
| 4 - samplesheet-consumable layout | Met |
| 5 - idempotent and resumable | Met |
| 6 - note on `prepare_production_data.sh` | Met |

AC2's full-data half is carved into a follow-up Issue, with a native `blockedBy` edge from [Issue #1176](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1176) onto it.
[Issue #1176](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1176) stays correctly blocked on real data, rather than appearing unblocked merely because a tool now exists.

## Rejected alternatives

### A. Build the whole tool from scratch (the first draft)

Rejected on prior art.
Retained because its constraint analysis is what makes D1 legible, and because it is the fallback if D6's viability check fails.

That design downloaded the `.sra` directly from ODP, verified the SDL-provided MD5 itself, and converted with **`fastq-dump --split-3 --gzip`**.
It chose `fastq-dump` over `fasterq-dump` because, per the [NCBI wiki](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump), `fasterq-dump` needs roughly **17x the accession size** during conversion (about 104 GiB for a 6.12 GiB run) against 39 GiB free, and offers no `--gzip`.

**The correctness tax that choice carried, and which the first draft failed to weigh:** `fastq-dump --split-3` gives both mates identical read IDs unless `--readids` is passed, and `--readids` in turn breaks BWA.
Worse, default read-length filtering discards short reads while keeping their mates, silently desynchronising the pair unless `-M 0` is passed.
There is also an open sra-tools issue on improper pairing for BAM-submitted entries.

So this path trades a disk problem for a correctness problem, in a pipeline whose premise is junction accuracy.
If it is ever revived, it must carry `--split-3 -M 0 --gzip` with no `--readids`, and the mate-count assertion becomes the falsifier for residual desync rather than a formality.

### B. `fasterq-dump` directly

Rejected: infeasible on this host at ~17x the accession size, per the sizing above.

### C. A Snakemake rule

Rejected: acquisition happens before the DAG, the Issue scopes this alongside `prepare_production_data.sh`, and it would drag a multi-GiB network fetch into rule execution.

## Out of scope

- The genome-wide pipeline run itself.
- Samplesheet authoring, the normal-filter design decision, the MS search, subset FDR, entrapment, and the known-answer control. All Scientist.
- Any change to `prepare_production_data.sh` beyond the AC6 pointer note.
- Vendoring, patching or linking kingfisher source. It is invoked as a CLI, which is also what keeps GPL-3.0 and MIT non-interacting.

## Resolved during implementation

- The kingfisher container image digest, recorded in the PR and in provenance.
- Which download method kingfisher actually selects for these accessions, confirmed by D6's viability check rather than assumed to be `aws-http`.
- The precise read-cap value for the subset smoke.
