"""Map biological strandness to HISAT2 --rna-strandness flag values.

Used by alignment.smk to translate per-sample `strandness` column values
(`unstranded` / `forward` / `reverse`) into the SE/PE flag pair that HISAT2
expects (`F`/`R` for single-end, `FR`/`RF` for paired-end).
"""

_VALID = {"", "unstranded", "forward", "reverse"}


def get_strandness_flag(strandness, is_paired_end):
    """Return the HISAT2 --rna-strandness value for a (strandness, SE/PE) pair.

    Missing values (`None`, empty string, whitespace) and `unstranded` return
    `""` — the caller then omits the `--rna-strandness` flag entirely. Any
    non-empty unrecognized string raises `ValueError`: a typo like `forwrd`
    would otherwise silently produce unstranded alignment for that sample.
    """
    normalized = (strandness or "").strip().lower()
    if normalized not in _VALID:
        raise ValueError(
            f"Unrecognized strandness value {strandness!r}. "
            "Expected: unstranded | forward | reverse (or empty/missing)."
        )
    if normalized == "forward":
        return "FR" if is_paired_end else "F"
    if normalized == "reverse":
        return "RF" if is_paired_end else "R"
    return ""


def get_strandness_from_row(row):
    """Resolve the HISAT2 strandness flag for a samples.tsv row dict.

    Detects SE vs PE from the presence of a non-empty `fastq2` field and
    delegates to `get_strandness_flag`.
    """
    is_pe = bool((row.get("fastq2") or "").strip())
    return get_strandness_flag(row.get("strandness"), is_pe)
