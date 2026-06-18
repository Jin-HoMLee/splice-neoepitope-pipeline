# Issue #299 — Kwok et al. PSR_GTEx validation (NEJ_GNAS / NEJ_RPL22)

**Issue:** [#299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/299) — re-derive `PSR_GTEx` for Kwok et al. NEJ_GNAS / NEJ_RPL22 to validate the DISCUSSION threshold-tradeoff
**Status:** complete (Issue closed)
**Manuscript link:** supports the Kwok threshold-tradeoff subsection in [PR #300](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/300) DISCUSSION

## Goal

Re-derive `PSR_GTEx` for the two validated public neojunctions in [Kwok et al. *Nature* 2025](https://doi.org/10.1038/s41586-024-08552-0) — `NEJ_GNAS` and `NEJ_RPL22` — and test whether Kwok's `<1%` GTEx threshold was the discriminator for these targets, or whether our stricter `min_read_count: 1` filter would also have caught them. Replaces a population-level claim in the manuscript with a concrete target-specific check.

## Inputs

- Snaptron GTEx endpoint (queried live in-notebook; hg19 for GNAS, hg38 for RPL22 region per notebook §§).
- Kwok et al. neojunction definitions (donor/acceptor shifts) transcribed from the paper + SSNIP repo.

## Outputs

- None persisted — analysis is self-contained; results (PSR_GTEx values, per-sample tallies, interpretation) live in the notebook. The `outputs/` slot is kept empty (`.gitkeep`) per the `research/experiments/` convention layout.

## Result

- 🎯 **NEJ_RPL22** — tradeoff IS real and target-specific: detected in 1/9,662 GTEx samples (1 read), `PSR_GTEx = 0.0000%` (below Kwok's `<1%`, so Kwok keeps it) but our `min_read_count: 1` rejects it. A documented case of a validated public neoepitope our pipeline would reject (sample-membership caveat in notebook §9).
- 🟡 **NEJ_GNAS** — tradeoff does not apply: no matching acceptor shift detectable in GTEx; both filters behave the same.

## Caveats

- **Snaptron GTEx ≠ Kwok's exact 9,166-sample subset** — robustness of the RPL22 call depends on whether the single positive sample is in Kwok's subset (notebook §9).
- **Coordinate/assembly care** — GNAS queried on hg19, RPL22 region on hg38; see per-section notes.
