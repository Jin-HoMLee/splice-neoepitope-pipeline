# Issue #1204 - junction-ORF search-DB artifacts

This experiment directory holds two artifacts for the junction-ORF search-DB work:

- **Sage recovery smoke** (this file, below) - the AC-3 live verification (Issue #1204).
- **Slide deck** `slides.qmd` - the experiment-tier deck (Issue #1262). Render with
  `quarto render slides.qmd --to revealjs`; figures regenerate via
  `figures/_regenerate_figures.py` (reads `outputs/orf_run_summary.json` and the cached
  real rulegraph `figures/_rulegraph.dot`). Rendered `slides.html` is gitignored.

## Sage junction-peptide recovery smoke

One-time pre-merge verification that the emitted ORF-stretch FASTA is consumable
by a nonspecific MS search AND that a real junction peptide is recovered from it
(AC-3). NOT run in CI; the CI-side check is the structural test in
`workflow/tests/test_orf_fasta_from_contigs.py`.

## What it proves

The concatenated FASTA (canonical proteome + our junction ORF partition) is the
search *database*, not a sample: a nonspecific Sage search digests it in silico
and matches candidate peptides against observed spectra. The canonical proteome
is included as a *competing target* (per the #1176 run4 brief) so a junction
peptide only wins when it out-scores any canonical explanation - it is not a
pre-exclusion filter.

The smoke exercises the full path end to end: build the DB, generate a
theoretical spectrum for a **real junction peptide taken from the current
FASTA**, run the search, and assert the junction peptide is recovered. The
assertion can fail (a bogus peptide is not recovered), so a green run is
evidence, not ceremony.

## Files (committed)

- `sage_config.json` - minimal nonspecific Sage config (no-enzyme cleavage,
  8-11 aa peptide length, decoys on).
- `make_spectrum.py` - picks a genuine junction ORF-stretch peptide (9-11 aa)
  from the emitted FASTA and writes a theoretical MS2 spectrum (precursor +
  full b/y ion ladder from monoisotopic masses). Deriving the target from the
  current FASTA keeps the smoke valid as the junction set evolves - no
  hand-computed masses, no hardcoded peptide that could drift out of the DB.
- `run_sage_smoke.sh` - builds `combined.fasta`, generates the spectrum, runs
  `sage`, and asserts recovery of the junction peptide.

Run artifacts (`combined.fasta`, `synthetic.mgf`, `sage_out/`, logs) are
regenerated on every run and gitignored.

## Run

    conda activate snakemake
    snakemake --cores 4 --use-conda --configfile config/test_config.yaml -- results/patient_001_test/ms_search_db/junction_orf.fasta
    bash research/experiments/issue_1204_orf_fasta/run_sage_smoke.sh \
        results/patient_001_test/ms_search_db/junction_orf.fasta

Sage itself is not installed by this repo's conda envs; install it once locally
before running the smoke (the v0.14.7 aarch64-apple-darwin release binary from
`github.com/lazear/sage`, on PATH as `sage`).

## Pass criteria

`SAGE SMOKE OK - recovered junction peptide <SEQ>` and exit 0: Sage loaded and
digested the concatenated FASTA, searched the generated spectrum, and the
junction peptide appears in `sage_out/results.sage.tsv` attributed to its
junction protein ID. (Sage reports "0 at 1% FDR" on a single spectrum - FDR is
undefined for n=1 - but the PSM itself is at rank 1; statistical significance
against real spectra is #1176 AC-6, not this smoke.)

This is a manual, one-time local check. Run it once pre-merge and paste the
`SAGE SMOKE OK` line into the PR Test plan. It is intentionally not wired into
CI - Sage is a separate binary this repo does not manage, and the structural
correctness of the emitted FASTA (headers, alphabet, no stop/X, no duplicate IDs
against the canonical proteome) is already covered by
`workflow/tests/test_orf_fasta_from_contigs.py`.
