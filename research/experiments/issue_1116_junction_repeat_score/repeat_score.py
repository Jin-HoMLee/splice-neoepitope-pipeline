"""Anchor-vs-intron repeat-embedding score for splice junctions.

[Issue #1116](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1116).

The Portcullis (GigaScience 2018) repeat test, as a pure function.

## What it measures

For a junction, compare each **exonic anchor** against the intronic sequence
adjacent to the *opposite* splice site:

    5' side:  exon_left[-k:]   vs  intron[-k:]   -> hamming_5p
    3' side:  exon_right[:k]   vs  intron[:k]    -> hamming_3p

If an anchor closely matches the intronic sequence at the other end of the
intron, then the aligner had **no sequence basis** to place that anchor outside
the intron rather than inside it. The junction may be an artifact of a repeat
rather than a real splice event.

**Low Hamming distance = suspicious.** This inverts the usual "bigger is worse"
intuition, so it is stated in every docstring here and pinned in the tests.

## Why this score, and not `NH`

`NH` (the multimapper count) is **index-relative**: on a single-chromosome index
a genome-wide multimapper looks unique, which is what made the #919 chr22 A/B
degenerate and unfixable without a whole-genome run we cannot afford.

This score consults **only the reference sequence and the junction coordinates**.
No index, no BAM, no `NH` tag. So it is computable on the chr22 fixture we
already have, and it means the same thing there as it would genome-wide.

It is also the *shape* of signal that
[Issue #1122](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1122)
permits: a **junction-level score** that can be carried into the tumor/normal
comparison, rather than a per-read gate applied to each arm independently. The
per-arm gate is what manufactured a false `tumor_exclusive` candidate and got
the NH filter disqualified.

## What it is NOT

It is not a verdict. It is a feature. Nothing here decides whether a junction is
admitted - that decision needs the evaluation in `evaluate_chr22.py` first, and
then a Scientist call.
"""

from dataclasses import dataclass

# Portcullis compares a window of exonic anchor against the opposite intron end.
# 10 is short enough to stay inside the ~8 bp minimum anchor that regtools/STAR
# require, and long enough that a chance match is unlikely (4^-10).
DEFAULT_ANCHOR = 10


class JunctionTooCloseToContigEnd(ValueError):
    """Not enough flanking reference sequence to build a full-length anchor.

    Raised rather than silently truncating the anchor. A truncated anchor would
    be scored over fewer positions than a full one, biasing junctions near
    contig ends toward LOW (suspicious) distances - a systematic artifact that
    would be indistinguishable from a real finding.
    """


def hamming(a: str, b: str) -> int:
    """Number of differing positions between two equal-length sequences.

    Case-insensitive: the reference soft-masks repeats in lowercase, and repeats
    are precisely the regions this score exists to detect - treating `a` != `A`
    would make every soft-masked repeat look maximally *un*-repeat-like, which is
    exactly backwards.

    `N` is unknown, not a wildcard: `N` vs `N` counts as a **mismatch**. A
    generous `N` would drag N-rich regions toward a low (suspicious) distance and
    flag their junctions as artifacts on no evidence.
    """
    if len(a) != len(b):
        raise ValueError(f"hamming needs equal lengths, got {len(a)} and {len(b)}")
    a = a.upper()
    b = b.upper()
    return sum(1 for x, y in zip(a, b) if x != y or x == "N")


@dataclass(frozen=True)
class RepeatScore:
    """Repeat-embedding distances for one junction. Lower = more suspicious."""

    hamming_5p: int
    hamming_3p: int

    @property
    def min_hamming(self) -> int:
        """The more suspicious of the two sides.

        `min`, not mean: an artifact needs only **one** ambiguous anchor to be
        placeable somewhere else. Averaging would let a clean side mask a
        perfectly repeat-embedded one.
        """
        return min(self.hamming_5p, self.hamming_3p)


def repeat_embedding_score(
    reference: str,
    intron_start: int,
    intron_end: int,
    anchor: int = DEFAULT_ANCHOR,
) -> RepeatScore:
    """Score one junction for repeat embedding.

    Args:
        reference:    Contig sequence (the whole chromosome).
        intron_start: First intron base, **1-based inclusive** (STAR `SJ.out.tab` col 2).
        intron_end:   Last intron base, **1-based inclusive** (STAR `SJ.out.tab` col 3).
        anchor:       Anchor window length, in bases.

    Returns:
        RepeatScore. **Low distances mean suspicious** (repeat-embedded).

    Raises:
        ValueError: if the intron is shorter than the anchor - the two compared
            windows would then overlap inside the intron, so the two distances
            stop being independent measurements and the score is meaningless
            rather than merely noisy.
        JunctionTooCloseToContigEnd: if there is not enough flanking sequence for
            a full-length anchor on either side.

    Strand is deliberately not a parameter: Hamming distance between two genomic
    windows is strand-agnostic, and reverse-complementing both windows leaves the
    per-position comparison unchanged.
    """
    intron_len = intron_end - intron_start + 1
    if intron_len < anchor:
        raise ValueError(
            f"intron of {intron_len} bp is shorter than the {anchor} bp anchor; "
            "the compared windows would overlap and the two distances would not "
            "be independent"
        )

    # Convert to 0-based half-open.
    i0 = intron_start - 1
    i1 = intron_end  # exclusive

    left_anchor_start = i0 - anchor
    right_anchor_end = i1 + anchor
    if left_anchor_start < 0 or right_anchor_end > len(reference):
        raise JunctionTooCloseToContigEnd(
            f"junction at {intron_start}-{intron_end} lacks {anchor} bp of flanking "
            f"sequence within a contig of {len(reference)} bp"
        )

    left_anchor = reference[left_anchor_start:i0]       # exonic, 5' of the donor
    right_anchor = reference[i1:right_anchor_end]       # exonic, 3' of the acceptor
    intron_5p = reference[i0 : i0 + anchor]             # first bases of the intron
    intron_3p = reference[i1 - anchor : i1]             # last bases of the intron

    # Each anchor is compared against the intron end at the OPPOSITE splice site:
    # that is the placement the aligner could have chosen instead.
    return RepeatScore(
        hamming_5p=hamming(left_anchor, intron_3p),
        hamming_3p=hamming(right_anchor, intron_5p),
    )
