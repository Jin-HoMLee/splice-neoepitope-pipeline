"""Build HISAT2 read-input arguments (single-end vs paired-end).

Used by `alignment.smk` to translate a sample's fastq paths into the HISAT2
read-input flag form: `-U <r1>` for single-end, `-1 <r1> -2 <r2>` for
paired-end. Extracted from the rule's inline shell so the single-end /
paired-end selection is unit-testable (see `test_hisat2_command.py`),
mirroring the `strandness.py` helper. Paired-end has no local end-to-end
coverage otherwise (the chr22 smoke fixture is single-end only) and its only
prior coverage (the GCP patient_002 runs) is gone with the decommissioned
infra - see Issue #962.
"""


def build_read_args(fastq1, fastq2):
    """Return the HISAT2 read-input args for a sample.

    Paired-end (`-1 <r1> -2 <r2>`) when a non-empty `fastq2` is given, else
    single-end (`-U <r1>`). `fastq2` may be `None`, an empty/whitespace-only
    string, or a path; whitespace-only is treated as single-end, matching the
    samples.tsv convention in `strandness.get_strandness_from_row`.
    """
    r1 = str(fastq1 or "").strip()
    r2 = str(fastq2 or "").strip()
    if not r1:
        raise ValueError("fastq1 is required for HISAT2 alignment (got empty).")
    if r2:
        return f"-1 {r1} -2 {r2}"
    return f"-U {r1}"
