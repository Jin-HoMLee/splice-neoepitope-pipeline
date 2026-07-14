"""Tests for the anchor-vs-intron repeat-embedding score.

[Issue #1116](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1116).

The score is the Portcullis "is this splice site embedded in a repeat" test
(GigaScience 2018): compare each exonic **anchor** against the intronic sequence
adjacent to the *opposite* splice site. If they match closely, the aligner could
have placed that anchor inside the intron just as well as outside it, so the
junction may be an alignment artifact rather than a real splice event.

    5' side:   ...exonL[-k:] | INTRON ...................... intron[-k:] | exonR...
               \\___________/                               \\__________/
                left anchor          compared against        intron 3' end

    3' side:   ...exonL | intron[:k] ..................... INTRON | exonR[:k]...
                         \\_________/                                \\________/
                        intron 5' start   compared against          right anchor

**Low Hamming distance = suspicious** (repeat-embedded). This is the inverse of
the usual intuition, so the tests below pin the direction explicitly.

Why this test and not the NH gate: it is a **pure sequence test on the junction**
and never consults the alignment index, so it is immune to the false-unique
artifact that made the #919 chr22 A/B degenerate (`NH` is index-relative: on a
single-chromosome index, genome-wide multimappers look unique). And per
[Issue #1122](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1122)
it is the *shape* of signal we are allowed to use - a junction-level score, not a
per-arm per-read gate.
"""

import pytest

from repeat_score import (
    JunctionTooCloseToContigEnd,
    hamming,
    repeat_embedding_score,
)


class TestHamming:
    def test_counts_mismatching_positions(self):
        assert hamming("ACGT", "ACGT") == 0
        assert hamming("ACGT", "ACGA") == 1
        assert hamming("ACGT", "TGCA") == 4

    def test_is_case_insensitive(self):
        """Reference FASTA soft-masks repeats in lowercase.

        Treating `a` != `A` would make every soft-masked repeat look like a
        perfect mismatch - i.e. maximally *un*-repeat-like. That would inv&ert
        the score exactly where it matters most.
        """
        assert hamming("acgt", "ACGT") == 0
        assert hamming("AcGt", "AcGa") == 1

    def test_rejects_unequal_lengths(self):
        with pytest.raises(ValueError):
            hamming("ACGT", "ACG")

    def test_n_base_counts_as_a_mismatch(self):
        """`N` is unknown, not "matches anything" - it must not lower the score.

        A generous `N` would make N-rich regions look repeat-embedded (low
        distance) and get their junctions flagged as artifacts.
        """
        assert hamming("ACGT", "ACGN") == 1
        assert hamming("NNNN", "NNNN") == 4


class TestRepeatEmbeddingScore:
    """Intron coords are 1-based inclusive, matching STAR `SJ.out.tab` cols 2-3."""

    def test_perfect_repeat_scores_zero_on_both_sides(self):
        """The pathological case: anchors identical to the opposite intron end.

        Reference is built so the left anchor (`AAAACCCCGG`) is byte-identical to
        the intron's last 10 bases, and the right anchor to the intron's first 10.
        An aligner has no sequence basis to prefer this junction over a shifted
        one, so both distances are 0 - maximally suspicious.
        """
        left_anchor = "AAAACCCCGG"
        right_anchor = "TTTTGGGGAA"
        intron = right_anchor + "CATCATCATCAT" + left_anchor
        ref = "GC" + left_anchor + intron + right_anchor + "GC"

        intron_start = len("GC" + left_anchor) + 1          # 1-based first intron base
        intron_end = intron_start + len(intron) - 1         # 1-based last intron base

        score = repeat_embedding_score(ref, intron_start, intron_end, anchor=10)

        assert score.hamming_5p == 0
        assert score.hamming_3p == 0
        assert score.min_hamming == 0

    def test_unrelated_sequence_scores_high(self):
        """A genuine junction: anchors bear no resemblance to the intron ends."""
        left_anchor = "AAAAAAAAAA"
        right_anchor = "CCCCCCCCCC"
        intron = "GGGGGGGGGG" + "ATATATATAT" + "TTTTTTTTTT"
        ref = "GC" + left_anchor + intron + right_anchor + "GC"

        intron_start = len("GC" + left_anchor) + 1
        intron_end = intron_start + len(intron) - 1

        score = repeat_embedding_score(ref, intron_start, intron_end, anchor=10)

        # left anchor AAAA... vs intron 3' end TTTT... -> all 10 differ
        assert score.hamming_5p == 10
        # right anchor CCCC... vs intron 5' start GGGG... -> all 10 differ
        assert score.hamming_3p == 10
        assert score.min_hamming == 10

    def test_min_hamming_takes_the_more_suspicious_side(self):
        """One repeat-embedded side is enough to make a junction doubtful.

        `min` (not mean) is the summary, because an artifact needs only ONE
        ambiguous anchor to be placeable elsewhere.
        """
        left_anchor = "AAAAAAAAAA"
        right_anchor = "CCCCCCCCCC"
        # intron 5' start == right anchor exactly -> hamming_3p == 0
        # intron 3' end   == unrelated            -> hamming_5p == 10
        intron = right_anchor + "ATATATATAT" + "TTTTTTTTTT"
        ref = "GC" + left_anchor + intron + right_anchor + "GC"

        intron_start = len("GC" + left_anchor) + 1
        intron_end = intron_start + len(intron) - 1

        score = repeat_embedding_score(ref, intron_start, intron_end, anchor=10)

        assert score.hamming_3p == 0
        assert score.hamming_5p == 10
        assert score.min_hamming == 0, "one embedded side is enough"

    def test_soft_masked_reference_scores_the_same_as_uppercase(self):
        """Repeat regions arrive soft-masked; the score must not depend on case."""
        left_anchor = "AAAACCCCGG"
        right_anchor = "TTTTGGGGAA"
        intron = right_anchor + "CATCATCATCAT" + left_anchor
        ref = "GC" + left_anchor + intron + right_anchor + "GC"

        intron_start = len("GC" + left_anchor) + 1
        intron_end = intron_start + len(intron) - 1

        upper = repeat_embedding_score(ref, intron_start, intron_end, anchor=10)
        lower = repeat_embedding_score(ref.lower(), intron_start, intron_end, anchor=10)

        assert (upper.hamming_5p, upper.hamming_3p) == (lower.hamming_5p, lower.hamming_3p)

    def test_intron_shorter_than_the_anchor_is_rejected(self):
        """A short intron would make the two compared windows OVERLAP.

        With intron length < anchor, `intron[:k]` and `intron[-k:]` share bases,
        so the two distances stop being independent and the score is meaningless
        rather than merely noisy. Refuse instead of returning a number.
        """
        ref = "AAAAAAAAAA" + "GTAG" + "CCCCCCCCCC"
        with pytest.raises(ValueError, match="intron"):
            repeat_embedding_score(ref, 11, 14, anchor=10)

    def test_junction_too_close_to_contig_start_is_rejected(self):
        """No room for a full left anchor -> refuse, never silently truncate.

        A truncated anchor would be compared against a full-length window, and
        the Hamming distance would be computed over fewer positions - producing a
        systematically LOWER (more suspicious) score for junctions near contig
        ends. A silent bias toward flagging edge junctions is worse than a gap.
        """
        left = "AAA"  # only 3 bases before the intron, anchor wants 10
        intron = "GGGGGGGGGGGGGGGGGGGG"
        ref = left + intron + "CCCCCCCCCC"
        with pytest.raises(JunctionTooCloseToContigEnd):
            repeat_embedding_score(ref, len(left) + 1, len(left) + len(intron), anchor=10)

    def test_junction_too_close_to_contig_end_is_rejected(self):
        right = "CCC"  # only 3 bases after the intron
        intron = "GGGGGGGGGGGGGGGGGGGG"
        ref = "AAAAAAAAAA" + intron + right
        with pytest.raises(JunctionTooCloseToContigEnd):
            repeat_embedding_score(ref, 11, 10 + len(intron), anchor=10)
