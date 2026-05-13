"""Map biological strandness to HISAT2 --rna-strandness flag values.

Used by alignment.smk to translate per-sample `strandness` column values
(`unstranded` / `forward` / `reverse`) into the SE/PE flag pair that HISAT2
expects (`F`/`R` for single-end, `FR`/`RF` for paired-end).
"""

_FORWARD = {"forward"}
_REVERSE = {"reverse"}


def get_strandness_flag(strandness, is_paired_end):
    """Return the HISAT2 --rna-strandness value for a (strandness, SE/PE) pair.

    Returns the empty string when unstranded, unrecognized, or missing — the
    caller suppresses the --rna-strandness flag entirely in that case.
    """
    normalized = (strandness or "").strip().lower()
    if normalized in _FORWARD:
        return "FR" if is_paired_end else "F"
    if normalized in _REVERSE:
        return "RF" if is_paired_end else "R"
    return ""
