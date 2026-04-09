"""Tests for translate_peptides.py — junction-spanning 9-mer extraction."""

from translate_peptides import extract_spanning_9mers

# Default config: upstream_nt=26, so min_start=2, max_start=23
UPSTREAM_NT = 26


def _no_stop_contig(length: int = 50) -> str:
    """Return a contig of given length with no stop codons in any frame.

    'GCC' * n: frame0=Ala, frame1=Pro (CCG), frame2=Arg (CGC) — all non-stop.
    """
    return ("GCC" * (length // 3 + 1))[:length]


class TestExtractSpanning9mers:
    def test_returns_list_of_tuples(self):
        seq = _no_stop_contig(50)
        result = extract_spanning_9mers(seq, upstream_nt=UPSTREAM_NT)
        assert isinstance(result, list)
        for item in result:
            start_nt, peptide = item
            assert isinstance(start_nt, int)
            assert isinstance(peptide, str)
            assert len(peptide) == 9

    def test_all_starts_within_valid_range(self):
        # min_start=2, max_start=23
        seq = _no_stop_contig(50)
        result = extract_spanning_9mers(seq, upstream_nt=UPSTREAM_NT)
        for start_nt, _ in result:
            assert 2 <= start_nt <= 23, f"start_nt={start_nt} outside [2, 23]"

    def test_frame0_start_zero_excluded(self):
        # Frame 0 first window starts at 0 < min_start=2 → excluded
        seq = _no_stop_contig(50)
        result = extract_spanning_9mers(seq, upstream_nt=UPSTREAM_NT, reading_frames=[0])
        starts = [s for s, _ in result]
        assert 0 not in starts

    def test_frame2_start_two_included(self):
        # Frame 2 first window starts at 2 == min_start=2 → included
        seq = _no_stop_contig(50)
        result = extract_spanning_9mers(seq, upstream_nt=UPSTREAM_NT, reading_frames=[2])
        starts = [s for s, _ in result]
        assert 2 in starts

    def test_window_start_24_excluded(self):
        # start=24 > max_start=23 → excluded in all frames
        seq = _no_stop_contig(50)
        result = extract_spanning_9mers(seq, upstream_nt=UPSTREAM_NT)
        starts = [s for s, _ in result]
        assert 24 not in starts

    def test_stop_codon_excluded(self):
        # Insert a stop codon (TAA) starting at nt 6 (frame 0, window index 2).
        # Any 27 nt window containing nt 6-8 is discarded for frame 0.
        seq = list(_no_stop_contig(50))
        seq[6], seq[7], seq[8] = "T", "A", "A"  # TAA at nt 6-8
        result = extract_spanning_9mers("".join(seq), upstream_nt=UPSTREAM_NT, reading_frames=[0])
        for _, peptide in result:
            assert "*" not in peptide

    def test_x_ambiguous_excluded(self):
        seq = list(_no_stop_contig(50))
        # Write 'N' at nt 6-8 in frame 0 → will translate to X
        seq[6], seq[7], seq[8] = "N", "N", "N"
        result = extract_spanning_9mers("".join(seq), upstream_nt=UPSTREAM_NT, reading_frames=[0])
        for _, peptide in result:
            assert "X" not in peptide

    def test_short_contig_no_results(self):
        # Contig shorter than 27 nt (min window size) → no results
        result = extract_spanning_9mers("GCC" * 8, upstream_nt=UPSTREAM_NT)
        assert result == []

    def test_empty_sequence_returns_empty(self):
        result = extract_spanning_9mers("", upstream_nt=UPSTREAM_NT)
        assert result == []

    def test_single_frame_subset(self):
        seq = _no_stop_contig(50)
        all_frames = extract_spanning_9mers(seq, upstream_nt=UPSTREAM_NT, reading_frames=[0, 1, 2])
        frame0_only = extract_spanning_9mers(seq, upstream_nt=UPSTREAM_NT, reading_frames=[0])
        # frame0_only should be a strict subset
        assert set(frame0_only).issubset(set(all_frames))
        assert len(frame0_only) < len(all_frames)

    def test_known_peptide_sequence(self):
        # 50 nt of "AAA": all Lys (K) in frame 0.
        # Valid starts in frame 0: 3, 6, 9, 12, 15, 18, 21.
        seq = "A" * 50
        result = extract_spanning_9mers(seq, upstream_nt=UPSTREAM_NT, reading_frames=[0])
        assert all(peptide == "K" * 9 for _, peptide in result)
        starts = sorted(s for s, _ in result)
        assert starts == [3, 6, 9, 12, 15, 18, 21]
