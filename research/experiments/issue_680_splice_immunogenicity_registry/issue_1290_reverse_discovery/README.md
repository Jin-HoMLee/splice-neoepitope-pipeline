# Issue #1290 - reverse-discovery: sequence-remap IEDB assayed epitopes to splice neojunctions

Parent: [Issue #680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680) (splice-immunogenicity registry) · This sub-experiment: [Issue #1290](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1290)

**Status:** stage 1 of 4 complete (the IEDB pull). Stages 2 to 4 not started.

## Goal

All prior registry curation is *forward* (find splice neoantigens, then test immunogenicity) or annotation-based DB mining ([Issue #734](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/734): search IEDB for splice-*titled* references).
The annotation-based mine is the **floor**: it can only find epitopes whose authors already knew were splice-derived, and #734 confirmed splice epitopes are routinely filed as ordinary protein peptides.
This experiment is the **ceiling**: take every T-cell-assayed epitope regardless of annotation, remap the peptide *sequence* to the genome, and test junction-crossing independently.

The prize is the **negatives**.
The registry's binding constraint is its hard-negative class, and IEDB stores T-cell-assay negatives in bulk.

## Stages

| Stage | What | State |
|---|---|---|
| 1 | Pull the human-source T-cell-assayed epitope set (positive + negative) from the IQ-API | **done** |
| 2 | Peptide to source transcript, anchored on `parent_source_antigen` | not started |
| 3 | GENCODE gate: annotated junction means normal self-peptide and is discarded; only an unannotated/aberrant junction is a candidate | not started |
| 4 | Two-gate verify, dedup against `registry.tsv`, fold survivors with provenance | not started |

Stage 3 is the load-bearing one.
A peptide spanning an *annotated* junction is a normal self-peptide, not a neoantigen, so without that gate the yield number is meaningless.

## Stage 1 output

Regenerate (about 1 minute, no credentials needed):

```bash
research/.venv/bin/python \
  research/experiments/issue_680_splice_immunogenicity_registry/issue_1290_reverse_discovery/fetch_iedb_tcell_assays.py
```

| File | Committed | Notes |
|---|---|---|
| `fetch_iedb_tcell_assays.py` | yes | the regenerator |
| `outputs/iedb_human_source_tcell_assays.provenance.json` | yes | counts, column list, sha256 of the TSV |
| `outputs/iedb_human_source_tcell_assays.tsv` | **no** (gitignored, 15 MB) | regenerate with the command above |

The TSV is gitignored because it is network-derived and therefore fails the "offline-regenerable" tie-breaker in the 10-100 MB size band of [`docs/research_artifact_conventions.md`](../../../../docs/research_artifact_conventions.md).
It differs from the [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) AlphaGenome precedent in a way that matters: the IEDB IQ-API is **public, unauthenticated, and unmetered**, so the committed regenerator *is* the fetch path for a fresh clone, and no R2 mirror is needed to keep the artifact reachable.

## Numbers as of the 2026-07-28 pull

| Quantity | Count |
|---|---|
| Human-source T-cell assays | 66,358 |
| Carrying a `parent_source_antigen_name` anchor | 66,353 (99.99%) |
| Positive (all four `Positive*` grades) | 33,674 |
| Negative | 32,684 |

These reproduce the 2026-07-22 feasibility probe exactly.

## API facts worth not re-deriving

- Endpoint is PostgREST at `https://query-api.iedb.org/tcell_search`, public and unauthenticated.
- `offset` **without** `order` is refused with HTTP 400. The API declines to page inconsistently rather than doing it silently, so the puller always sends `order=tcell_id.asc`.
- A page is capped server-side at 10,000 rows, so the human-source pool is 7 pages.
- Exact counts come from a `Prefer: count=exact` request, read off the `Content-Range` header.
- The `qualitative_measure` vocabulary for this pool is exactly five values (`Positive`, `Positive-High`, `Positive-Intermediate`, `Positive-Low`, `Negative`). They sum to 66,358, which equals the unfiltered pool total, so the enumeration is exhaustive and there are no nulls. `normalise_outcome()` raises on anything else rather than defaulting, because silently bucketing a new IEDB category would move rows between the two classes.

## Known open risk

Yield is unproven.
Stage 1 establishes only that the *inputs* exist in quantity; it says nothing about how many reverse-map to real neojunctions.
An empty result remains a live outcome, and AC3 explicitly allows documenting that with numbers.

## Tests

`../tests/test_fetch_iedb_tcell_assays.py` (16 tests, run with `research/.venv/bin/python -m pytest`).
Every integrity check is paired with a case that must fail, so no check can only ever confirm.
