"""Build STAR's --readFilesIn argument (single-end vs paired-end).

Used by `alignment.smk` to translate a sample's fastq paths into STAR's
`--readFilesIn` value: `<r1>` for single-end, `<r1> <r2>` for paired-end.
Extracted from the `star_align` rule's inline shell so the single-end /
paired-end selection is unit-testable (see `test_star_command.py`), mirroring
`hisat2_command.build_read_args`. The one shape difference from the HISAT2 helper
is that STAR takes a bare space-separated file list, not `-U` / `-1 / -2` flags.
Consistency follow-up to Issue #962 (which extracted the HISAT2 side); the STAR
branch is behavior-correct today, this only removes the HISAT2/STAR asymmetry.
"""


def build_read_files_in(fastq1, fastq2):
    """Return STAR's --readFilesIn argument for a sample.

    Paired-end (`<r1> <r2>`) when a non-empty `fastq2` is given, else single-end
    (`<r1>`). `fastq2` may be `None`, an empty/whitespace-only string, or a path;
    whitespace-only is treated as single-end, matching `build_read_args` and the
    samples.tsv convention in `strandness.get_strandness_from_row`.
    """
    r1 = str(fastq1 or "").strip()
    r2 = str(fastq2 or "").strip()
    if not r1:
        raise ValueError("fastq1 is required for STAR alignment (got empty).")
    if r2:
        return f"{r1} {r2}"
    return r1
