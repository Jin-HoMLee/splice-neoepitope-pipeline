# #966 AC-1: ASNEO end-to-end runner in the caller-benchmark harness

**Status:** design, approved 2026-07-21.
**Issue:** [#966](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/966) (benchmark harness skeleton, leaf B of [#679](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/679)), AC-1.
**PR in flight:** [#1251](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/1251) (AC-2/3/5 landed; this design completes AC-1).

## Goal

Satisfy #966 AC-1: "Harness runs >=1 open caller end-to-end on a chr22-scale input."
Today the harness only *ingests* a hand-produced caller TSV.
This design adds a *run* step so the harness itself invokes a caller, produces its output, and collects it into the common schema.

## Decision summary

- **Caller = ASNEO** (Apache-2.0, `github.com/bm2-lab/ASNEO`), open-only via the [#566](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/566) option-B patch.
  It is a true end-to-end splice-neoantigen caller, its chr22 smoke already passes on arm64 (800 candidate peptides), and the `asneo` conda env exists on this box.
  splice2neo was the alternative but it is a junction-to-peptide library (our own stack's component, not an external caller), was only ever smoked on its own hg19 toy fixtures, and needs an hg19->hg38 liftover to touch our chr22 junctions.
- **ASNEO output is peptide-only by construction.**
  Its stop-point file `putative_peptide.txt` (`ASNEO.py:246-253`) is `peps - norm_peps`: a normal-subtracted *set* of k-mer peptide strings.
  The set operation is over bare strings, so the link from each peptide back to its source junction / novel isoform is discarded there, and the mapping is not even 1:1 (chr22 smoke: 60 isoforms -> 800 de-duplicated peptides).
- **Schema shape = an explicit `record_level` discriminator, not silent nullable fields.**
  Peptide-only callers emit `record_level="peptide"` with null junction-level fields; junction-anchored callers emit `record_level="junction"` and must carry a `junction_id`/`strand`.
  A legal-combination validator enforces this.
  Chosen over "just make the fields nullable" (the anti-pattern the schema-design literature warns against) and over two separate dataclasses (correct but heavier for a flat-TSV research harness).
- **Junction traceability for ASNEO is deferred to a new follow-up, [#1258](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1258)** (patch ASNEO to emit the source junction pre-subtraction), NOT to #1100 (which is the unrelated registry-side canonical `junction_id` scheme).
  Peptide-only ASNEO is consistent with how the project already plans to compare it: the #566 cross-check design is peptide-set (Jaccard) concordance, MHC-agnostic, precisely because ASNEO is peptide-only.

## Architecture

Three units, mirroring the existing adapter pattern, plus one schema change.

### 1. `runners/asneo.py` - the run step

`run_asneo(sj_tab_path, workdir) -> peptide_txt_path`.

- Invokes patched (option-B, no netMHCpan) ASNEO in the `asneo` conda env on a chr22-scale `SJ.out.tab`, returns the path to the produced `putative_peptide.txt`.
- Implementation: a Python function that shells out via `conda run -n asneo python ASNEO.py -j <sj_tab> ...` on a clone patched by #566's `apply_optionB_patch.py`, reusing the exact step sequence proven in the #965 `asneo_smoke.sh` (clone -> patch -> run at the relaxed chr22 thresholds). It is a subprocess wrapper, not a reimplementation of ASNEO.
- Two-phase, fail-loud: a validation phase (env present, patch anchors match exactly once, input `SJ.out.tab` well-formed) precedes the run; a post-run check asserts a non-empty `putative_peptide.txt` was produced.

### 2. `adapters/asneo.py` - the ingest step

`parse_asneo(path) -> list[CommonRecord]`.

- One record per candidate peptide: `caller="asneo"`, `peptide`, `genome_build="hg19"`, `event_type="junction"`, `record_level="peptide"`, `junction_id=None`, `strand=None`.
- `provenance = {source: "asneo", sj_tab_digest, thresholds}` so the input junction set is recoverable even though per-peptide linkage is not.
- Two-phase: validate the file exists and is the expected bare-peptide-per-line shape (skip comment/blank lines) before building records.

### 3. Registry - `RUNNERS` beside `ADAPTERS`

Add a `RUNNERS` dict in `collect.py` parallel to `ADAPTERS`, so registering a caller is a one-line entry in each: a runner (how to invoke it) and an adapter (how to parse it).
This extends the AC-5 extension point to cover the run step too.

### 4. Schema change - `common_schema.py`

- Add `record_level: str` (values `"junction"` | `"peptide"`); make `junction_id: Optional[str] = None` and `strand: Optional[str] = None`.
- Add a `validate(record)` legal-combination check: `record_level=="junction"` requires non-null `junction_id` and `strand`; `record_level=="peptide"` requires both null.
  Raise loudly on an illegal combination (the #1237 faceted legal-combination-matrix pattern).
- `record_level` is orthogonal to the existing `event_type` (which is splice-event *kind*, e.g. junction vs exon-skip, not record *granularity*).
- Update the splice2neo adapter to set `record_level="junction"` (it already fills `junction_id`/`strand`).

## Data flow

```
chr22 SJ.out.tab (committed fixture)
   -> runners.run_asneo (asneo conda env, option-B patch)
   -> putative_peptide.txt (peptide-only)
   -> adapters.parse_asneo -> list[CommonRecord] (record_level="peptide")
   -> collect merge + validate
   -> records_to_tsv -> outputs/asneo_chr22_unified.tsv   (ticks AC-1)
```

## CLI / invocation (AC-3)

Extend `collect.py` with a run-then-ingest mode alongside the existing ingest-only mode:

- `python collect.py --run asneo:<sj_tab> --out unified.tsv` - run the caller, then ingest.
- `python collect.py --input asneo:<putative_peptide.txt> --out unified.tsv` - ingest a pre-run file (unchanged).

Documented in the experiment README with the exact chr22 invocation.

## Inputs and artifacts (AC-4)

- Commit a chr22-subset `SJ.out.tab` (ASNEO's bundled `test/SRR2660032.SJ.out.tab` filtered to chr22, ~few hundred KB) as the deterministic AC-1 input.
- The hg19 chr22 FASTA (~12 MB) and the ASNEO clone stay in a gitignored scratch dir.
- No artifact exceeds 100 MB, so AC-4's `data_manifest.yaml` requirement stays vacuous (documented as such).

## Error handling

- Runner fails loud if: the `asneo` env is missing, the option-B patch anchors do not match exactly once, the input `SJ.out.tab` is malformed, or 0 peptides are produced at the given thresholds.
- Adapter fails loud on an unreadable / wrong-shape output file.
- `validate()` fails loud on any illegal `record_level` / field combination.

## Testing

- Unit: `parse_asneo` on a small peptide-only fixture -> records with `record_level="peptide"`, null junction/strand, comment/blank lines skipped.
- Unit: `validate()` matched-pair - a `record_level="junction"` record with a null `junction_id` raises; the same record with a `junction_id` passes (red-before / green-after).
- Live integration (per the live-integration-smoke rule): the actual chr22 ASNEO run recorded in the experiment, producing `outputs/asneo_chr22_unified.tsv`. This is what ticks AC-1.

## Downstream limitation (on the record)

- With `junction_id=None`, ASNEO records join **peptide-level** concordance (the #566 Jaccard design) but **not junction-level** detection concordance ([#968](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/968)) until [#1258](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1258) recovers the linkage.
- This is consistent with the existing plan (ASNEO was always going to be compared peptide-set-wise), and AC-1 ("harness runs a caller end-to-end") is fully met regardless.
- ASNEO does not *detect* junctions itself - it is fed a STAR `SJ.out.tab`, so the junction-detection signal lives in that input, not in ASNEO's peptide output.

## Best-practice grounding

- Explicit discriminator over nullable-for-all-options: "make illegal states unrepresentable" / tagged-variant guidance ([Inside.java DoP](https://inside.java/2024/06/03/dop-v1-1-illegal-states/), [F# for Fun and Profit](https://fsharpforfunandprofit.com/posts/designing-with-types-making-illegal-states-unrepresentable/)); internally consistent with our own #1237 faceted decision.
- Per-source adapter normalization into one common schema + per-caller conda env for reproducibility: [nf-core/variantbenchmarking](https://github.com/nf-core/variantbenchmarking).
- Self-registering registry for pluggable tools: standard plugin-registry pattern ([GeeksforGeeks](https://www.geeksforgeeks.org/system-design/registry-pattern/)).

## Out of scope (YAGNI)

- Manifest-driven registry entries (a machine-readable per-caller description of input type / emitted fields / genome build) would strengthen AC-5, but AC-1 does not need it. Noted, not built.
- ASNEO junction-linkage recovery: [#1258](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1258).
- splice2neo runner, SNAF (Linux/Docker-only): later increments.

## AC mapping

| AC | This design |
|---|---|
| AC-1 (run a caller end-to-end on chr22) | ASNEO runner + the recorded chr22 run producing `asneo_chr22_unified.tsv` |
| AC-2 (common schema) | extended with `record_level`; landed core in #1251 |
| AC-3 (CLI entry point) | `collect.py --run caller:sj_tab` + README |
| AC-4 (data_manifest for >100 MB) | vacuous (no >100 MB artifact); documented |
| AC-5 (low-friction extension point) | `RUNNERS` + `ADAPTERS` one-line register per caller |
