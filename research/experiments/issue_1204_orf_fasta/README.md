# Issue #1204 - Sage loadability smoke

One-time pre-merge verification that the emitted ORF-stretch FASTA is consumable
by a nonspecific MS search (AC-3).
NOT run in CI; the CI-side check is the structural test in `workflow/tests/test_orf_fasta_from_contigs.py`.

## Files

- `sage_config.json` - minimal nonspecific Sage config (no-enzyme cleavage,
  8-11 aa peptide length, decoys on).
- `synthetic.mgf` - one synthetic MS2 spectrum. Precursor and fragment m/z
  values correspond to `MKLPQRSTV`, a 9-mer substring of the canonical
  decoy/target sequence in `workflow/tests/data/orf_fasta/canonical_tiny.fasta`,
  so a hit is plausible but not required (see Pass criteria).
- `run_sage_smoke.sh` - concatenates the canonical proteome with an emitted
  `junction_orf.fasta` into `combined.fasta`, then invokes `sage` against it.

## Run

    conda activate snakemake
    snakemake --cores 4 --use-conda --configfile config/test_config.yaml -- results/patient_001_test/ms_search_db/junction_orf.fasta
    bash research/experiments/issue_1204_orf_fasta/run_sage_smoke.sh \
        results/patient_001_test/ms_search_db/junction_orf.fasta

Sage itself is not installed by this repo's conda envs; install it once
locally before running the smoke (e.g. `brew install sage-proteomics`, or a
release binary from `github.com/lazear/sage`, on PATH as `sage`).

## Pass criteria

Sage exits 0, writes `sage_out/`, and logs the loaded target+decoy count
without a FASTA parse error. Recovering the synthetic peptide is a bonus, not
required; the deliverable is loadability of the concatenated FASTA by a
nonspecific search, not a confirmed hit (real hits against real spectra are
#1176 AC-6).

This is a manual, one-time local check (see Task 5 / Step 6 in the design
doc): run it once pre-merge and paste the `SAGE SMOKE OK` line into the PR
Test plan. It is intentionally not wired into CI - Sage is a separate binary
this repo does not manage, and the structural correctness of the emitted
FASTA (headers, alphabet, no stop/X, no duplicate IDs against the canonical
proteome) is already covered by
`workflow/tests/test_orf_fasta_from_contigs.py`.
