# SRA-accession fetch to gzipped FASTQ for ENA-unmirrored cohorts

Design for [Issue #1296](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1296).
Author: Developer.
Date: 2026-07-28.

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

## Decisions

### D1. A thin shell launcher plus a tested Python core

The [Google Shell Style Guide](https://google.github.io/styleguide/shellguide.html) says to rewrite in a structured language past ~100 lines or non-straightforward control flow, but also that shell is appropriate when "mostly calling other utilities and doing relatively little data manipulation."
This tool has two halves that fall on opposite sides of that line.

Orchestration (checking for a Docker daemon, starting colima, invoking the converter) is utility-calling and stays in shell.
The assertion layer (parsing SDL JSON, comparing MD5s, computing and comparing mate counts, managing sentinels) is data manipulation and moves to Python, where it can be unit-tested.

Two project-local facts reinforce the split.
Local shell scripts must be bash-3.2-safe on macOS, which constrains a script of this size.
[Issue #1043](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1043) exists because untested shell has accumulated in this repo.

Three of this Issue's six acceptance criteria are about *proving* the output is correct, and a pure-bash implementation gives no way to prove the prover works.

### D2. `fastq-dump --split-3 --gzip`, not `fasterq-dump`

Per the [NCBI sra-tools wiki](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump), `fasterq-dump` produces final FASTQs about 7x the accession size and needs scratch of about 1.5x the FASTQ size, for a total of roughly **17x the accession** during conversion.
For a 6.12 GiB accession that is about 104 GiB.
The development host has 39 GiB free, so this path is not tight, it is impossible.

The same source confirms there is no `--gzip` option in `fasterq-dump`, so the uncompressed intermediate cannot be avoided there.

`fastq-dump --split-3 --gzip` writes compressed output directly, keeping peak usage near the accession size plus the compressed output.
It is meaningfully slower.
That trade is acceptable because the full data pull is deferred and unattended (see D5), so wall-clock is the cheap thing to spend.

`fastq-dump` is the older tool but is not, as far as could be established, officially deprecated.

### D3. Pin the converter image by digest

Container images are pinned by **digest** (`ncbi/sra-tools@sha256:...`), not by tag.
A tag is a mutable pointer and can silently resolve to different bytes over time; a digest is content-addressed and cannot.
This is standard guidance for reproducible scientific workflows.

The tag to resolve the digest *from* is chosen by evidence, not inheritance.
`AGENTS.md` records "use version 3.1.1 on GCP VMs - newer versions (3.4.x) have a segfault bug", and Docker Hub's `latest` currently resolves to 3.4.1, so accepting the default would walk directly into a documented failure.
Available tags are `latest`, `3.4.1`, `3.3.0`, `3.2.1`, `3.1.0`, `3.0.1`.
Note that no `3.1.1` image exists, only `3.1.0`.

Because verification runs on a bounded subset (D5) and therefore takes minutes, the candidate tag is tested directly against our own data rather than trusted on the strength of a GCP-era note.
Starting candidate is `3.3.0`.
Whatever passes is pinned by its digest and the digest is recorded.

NCBI publishes both `linux/amd64` and `linux/arm64` builds.
The `AGENTS.md` warning about sra-tools on macOS arm64 concerns *conda on macOS*, which is a different environment from NCBI's own Linux build in a container, so the native arm64 image is tried first (no Rosetta translation, materially faster on Apple Silicon) with `linux/amd64` as the guaranteed fallback.

### D4. The provenance record is the idempotency sentinel

AC5 requires that an already-completed run is skipped, which needs a completion marker regardless.
A zero-byte `.done` file carries no information, so the marker instead *is* the provenance record.
Existence means complete; contents answer what produced the data.

This adds no artifact beyond the one idempotency already required.
Reproducibility guidance (for example [QIIME 2 Provenance Replay](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1011676)) is consistent that a checksum alone is necessary but not sufficient, and that tool versions, parameters, and environment must also be recorded.
That matters here specifically because a version is being pinned around a segfault.

### D5. Verification is subset-first, with the full pull deferred

The full pull is roughly 12.6 GiB of download and, per D2's sizing, about 26 GiB of gzipped output across the two accessions.
Sequentially, deleting each `.sra` immediately after its conversion, the second run peaks near **32 GiB against 39 GiB free**.
It fits, with under 7 GiB of headroom.

The mechanism is therefore verified on a bounded read subset of **both** real accessions, which exercises every step in minutes and a few GiB.
The full-data pull runs separately and unattended.

The narrow headroom is why the preflight free-space check (E8 below) is load-bearing rather than a nicety, and why dropping the `.sra` the moment conversion succeeds is required rather than tidy.

## Architecture

Two units with one clear boundary.

### `scripts/fetch_sra_fastq.sh`

A launcher, targeted well under 100 lines.
Its entire job is establishing that the run *can* proceed, then handing off.

- Verify the `docker` CLI exists and the daemon is reachable, starting colima if not, reusing the pattern already proven in `scripts/run_local_linux64.sh`.
- Verify the repo is under `$HOME`, since colima mounts only `$HOME` and `/tmp/colima`, and a bind mount of any other path silently yields an empty directory inside the container.
- Exec the Python core, forwarding arguments unchanged.

No parsing, no assertions, no arrays.

### `scripts/sra_fetch.py`

Everything that can be wrong.
Tests live in `scripts/tests/test_sra_fetch.py`, which is one of the five directories the CI pytest job enumerates, so they actually run.

This is `scripts/`, not `workflow/scripts/`: nothing here is invoked by a Snakemake rule, and `workflow/scripts/` is reserved for rule-invoked code.

Public surface:

```
resolve_accession(acc)        -> AccessionMeta   (url, size, md5, layout)
already_complete(acc, outdir) -> bool
preflight_disk(meta, outdir)  -> None            (raises if insufficient)
download_sra(meta, cachedir)  -> Path            (resumable, MD5-verified)
convert(sra_path, outdir, …)  -> ConversionResult
assert_integrity(result)      -> None            (raises on any failure)
write_provenance(...)         -> Path
```

Each is independently testable and communicates through plain data.

## Data flow, per accession

1. **Resolve.** Query the NCBI SDL API for the download URL, exact size, and MD5. Query ENA for `library_layout`.
2. **Skip check.** If `data/<ACC>.provenance.json` exists, report and return. This is the first step, so a re-run costs one API call.
3. **Preflight disk.** Compute the requirement from the SDL-reported size and refuse to start if free space is insufficient, stating required versus available.
4. **Download** the `.sra` to `data/.sra_cache/<ACC>/`, resumable via HTTP range requests, with at most 3 attempts and exponential backoff between them.
5. **Verify MD5** against the SDL value.
6. **Convert** inside the digest-pinned `ncbi/sra-tools` container.
7. **Assert integrity** (see below).
8. **Drop the `.sra`**, freeing roughly 6 GiB before the next accession begins.
9. **Write the provenance sentinel last**, so it exists only if every prior step passed.

### Output layout

`data/<ACC>_1.fastq.gz` and `data/<ACC>_2.fastq.gz`, flat in `data/`, matching the convention of `prepare_production_data.sh`'s outputs.
These paths drop directly into a `config/samples/<patient>.tsv` samplesheet's `fastq1` and `fastq2` columns with no renaming, satisfying AC4.

The `.sra` scratch and cache live under `data/.sra_cache/`, so the machinery does not clutter the data directory.

A `SINGLE`-layout accession emits `data/<ACC>.fastq.gz` and records the layout, rather than asserting a pairing that never existed.

## Error handling

Every integrity check fails closed.
A download that cannot be verified is treated as corrupt, never as acceptable.
This follows the project's fail-safe-not-fail-open rule: a gate must never false-PASS.

| # | Failure | Behavior |
|---|---|---|
| E1 | SDL API unreachable, or accession unknown | Fatal, accession named |
| E2 | SDL returns no `sra` file type | Fatal, worded so a genuinely-unavailable accession does not read as a bug |
| E3 | MD5 mismatch | Fatal, **and the bad `.sra` is deleted** so a re-run cannot resume onto a poisoned file |
| E4 | Download interrupted | Resumable via range requests, at most 3 attempts with backoff; no sentinel was written, so a re-run continues |
| E5 | `docker` absent or daemon unreachable | Caught in the launcher before anything downloads, with install guidance |
| E6 | Converter exits non-zero | Fatal, partial outputs cleaned |
| E7 | Mate counts unequal | Fatal, **both counts printed** |
| E8 | Insufficient disk at preflight | Refuses to start, states required versus available |
| E9 | `--split-3` emits an orphan third file | Warning only. This is legitimate, and it is recorded in provenance |

## Testing

### Unit tests, `scripts/tests/test_sra_fetch.py`

The load-bearing ones are matched pairs: identical inputs, one variable flipped, opposite expected outcomes.
A test that can only pass is not a check.

| Control | Falsifier |
|---|---|
| Correct MD5 passes | Corrupted MD5 fails **and removes the file** |
| Equal mate counts pass | Unequal mate counts fail |
| Sentinel present skips without touching the network | Sentinel absent proceeds |
| Successful run writes a sentinel | Failed conversion writes **no** sentinel |
| Sufficient disk proceeds | Insufficient disk refuses |

Plus: SDL JSON parsing (happy path, missing `sra` type, HTTP error), and emitted paths matching samplesheet expectations.

### Live subset smoke

Both real accessions, read-capped: real outputs, real `gzip -t`, real equal mate counts, real provenance written, and a real second invocation that skips.

Unit fixtures validate our *model* of the SDL API, not the API itself, so the live smoke is required by the project's live-integration-smoke rule for any change that touches an external boundary.

## Acceptance criteria disposition

| AC | Disposition |
|---|---|
| 1 - fetch to gzipped FASTQ without host sra-tools | Met |
| 2 - verified end to end on both accessions | **Mechanism half met** (subset, both accessions). Full-data half deferred to a carrier Issue |
| 3 - integrity asserted, not assumed | Met, and exceeded: upstream MD5 in addition to `gzip -t` and mate-count equality |
| 4 - samplesheet-consumable layout | Met |
| 5 - idempotent and resumable | Met |
| 6 - note on `prepare_production_data.sh` | Met |

AC2's full-data half is carved into a follow-up Issue, with a native `blockedBy` edge from [Issue #1176](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1176) onto it.
[Issue #1176](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1176) stays correctly blocked on real data, rather than appearing unblocked merely because a tool now exists.

## Out of scope

- The genome-wide pipeline run itself.
- Samplesheet authoring, the normal-filter design decision, the MS search, subset FDR, entrapment, and the known-answer control. All Scientist.
- Any change to `prepare_production_data.sh` beyond the AC6 pointer note. The two scripts sit alongside each other; neither replaces the other.

## Resolved during implementation

- The exact image digest, and the tag it was resolved from, decided by the subset run and recorded in the PR.
- Whether the `linux/arm64` image works, with `linux/amd64` as the fallback.
- The precise read-cap value for the subset smoke.
