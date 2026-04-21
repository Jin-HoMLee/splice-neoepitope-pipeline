"""Tests for translate_peptides.py — junction-spanning peptide extraction."""

from translate_peptides import extract_spanning_peptides

# Default config: upstream_nt=27, contig_length=54
UPSTREAM_NT = 27
CONTIG_LEN = 54


def _no_stop_contig(length: int = CONTIG_LEN) -> str:
    """Return a contig of given length with no stop codons in any frame.

    'GCC' * n: frame0=Ala, frame1=Pro (CCG), frame2=Arg (CGC) — all non-stop.
    """
    return ("GCC" * (length // 3 + 1))[:length]


class TestExtractSpanningPeptides:
    def test_returns_list_of_tuples(self):
        seq = _no_stop_contig()
        result = extract_spanning_peptides(seq, upstream_nt=UPSTREAM_NT, peptide_lengths=[9])
        assert isinstance(result, list)
        for item in result:
            start_nt, peptide = item
            assert isinstance(start_nt, int)
            assert isinstance(peptide, str)
            assert len(peptide) == 9

    def test_9mer_starts_within_valid_range(self):
        # For 9-mers: min_start=3, max_start=24
        seq = _no_stop_contig()
        result = extract_spanning_peptides(seq, upstream_nt=UPSTREAM_NT, peptide_lengths=[9])
        for start_nt, _ in result:
            assert 3 <= start_nt <= 24, f"start_nt={start_nt} outside [3, 24]"

    def test_8mer_starts_within_valid_range(self):
        # For 8-mers: min_start=27-21=6, max_start=24
        seq = _no_stop_contig()
        result = extract_spanning_peptides(seq, upstream_nt=UPSTREAM_NT, peptide_lengths=[8])
        for start_nt, _ in result:
            assert 6 <= start_nt <= 24, f"start_nt={start_nt} outside [6, 24]"

    def test_10mer_starts_within_valid_range(self):
        # For 10-mers: min_start=27-27=0, max_start=24
        seq = _no_stop_contig()
        result = extract_spanning_peptides(seq, upstream_nt=UPSTREAM_NT, peptide_lengths=[10])
        for start_nt, _ in result:
            assert 0 <= start_nt <= 24, f"start_nt={start_nt} outside [0, 24]"

    def test_peptide_length_matches_requested(self):
        seq = _no_stop_contig()
        for length in [8, 9, 10]:
            result = extract_spanning_peptides(seq, upstream_nt=UPSTREAM_NT, peptide_lengths=[length])
            for _, peptide in result:
                assert len(peptide) == length

    def test_multi_length_contains_all_lengths(self):
        seq = _no_stop_contig()
        result = extract_spanning_peptides(seq, upstream_nt=UPSTREAM_NT, peptide_lengths=[8, 9, 10])
        lengths_found = {len(p) for _, p in result}
        assert lengths_found == {8, 9, 10}

    def test_frame0_start_zero_excluded_for_9mer(self):
        # Frame 0 first window starts at 0 < min_start=3 for 9-mers → excluded
        seq = _no_stop_contig()
        result = extract_spanning_peptides(seq, upstream_nt=UPSTREAM_NT, peptide_lengths=[9], reading_frames=[0])
        starts = [s for s, _ in result]
        assert 0 not in starts

    def test_frame0_start_zero_included_for_10mer(self):
        # For 10-mers: min_start=0, so start=0 in frame 0 → included
        seq = _no_stop_contig()
        result = extract_spanning_peptides(seq, upstream_nt=UPSTREAM_NT, peptide_lengths=[10], reading_frames=[0])
        starts = [s for s, _ in result]
        assert 0 in starts

    def test_stop_codon_excluded(self):
        # Insert a stop codon (TAA) starting at nt 6 in frame 0.
        # Any window containing nt 6-8 is discarded for frame 0.
        seq = list(_no_stop_contig())
        seq[6], seq[7], seq[8] = "T", "A", "A"
        result = extract_spanning_peptides("".join(seq), upstream_nt=UPSTREAM_NT, peptide_lengths=[9], reading_frames=[0])
        for _, peptide in result:
            assert "*" not in peptide

    def test_x_ambiguous_excluded(self):
        seq = list(_no_stop_contig())
        seq[6], seq[7], seq[8] = "N", "N", "N"
        result = extract_spanning_peptides("".join(seq), upstream_nt=UPSTREAM_NT, peptide_lengths=[9], reading_frames=[0])
        for _, peptide in result:
            assert "X" not in peptide

    def test_short_contig_no_results(self):
        # Contig shorter than 24 nt (min 8-mer window) → no results
        result = extract_spanning_peptides("GCC" * 7, upstream_nt=UPSTREAM_NT)
        assert result == []

    def test_empty_sequence_returns_empty(self):
        result = extract_spanning_peptides("", upstream_nt=UPSTREAM_NT)
        assert result == []

    def test_single_frame_subset(self):
        seq = _no_stop_contig()
        all_frames = extract_spanning_peptides(seq, upstream_nt=UPSTREAM_NT, peptide_lengths=[9], reading_frames=[0, 1, 2])
        frame0_only = extract_spanning_peptides(seq, upstream_nt=UPSTREAM_NT, peptide_lengths=[9], reading_frames=[0])
        assert set(frame0_only).issubset(set(all_frames))
        assert len(frame0_only) < len(all_frames)

    def test_known_9mer_starts_frame0(self):
        # 54 nt of "AAA": all Lys (K) in frame 0.
        # For 9-mers with upstream_nt=27: valid frame-0 starts in [3, 24]: 3,6,...,24
        seq = "A" * CONTIG_LEN
        result = extract_spanning_peptides(seq, upstream_nt=UPSTREAM_NT, peptide_lengths=[9], reading_frames=[0])
        assert all(peptide == "K" * 9 for _, peptide in result)
        starts = sorted(s for s, _ in result)
        assert starts == [3, 6, 9, 12, 15, 18, 21, 24]

    def test_known_10mer_starts_frame0(self):
        # For 10-mers with upstream_nt=27: valid frame-0 starts in [0, 24]: 0,3,...,24
        seq = "A" * CONTIG_LEN
        result = extract_spanning_peptides(seq, upstream_nt=UPSTREAM_NT, peptide_lengths=[10], reading_frames=[0])
        assert all(peptide == "K" * 10 for _, peptide in result)
        starts = sorted(s for s, _ in result)
        assert starts == [0, 3, 6, 9, 12, 15, 18, 21, 24]

    def test_known_8mer_starts_frame0(self):
        # For 8-mers with upstream_nt=27: valid frame-0 starts in [6, 24]: 6,9,...,24
        seq = "A" * CONTIG_LEN
        result = extract_spanning_peptides(seq, upstream_nt=UPSTREAM_NT, peptide_lengths=[8], reading_frames=[0])
        assert all(peptide == "K" * 8 for _, peptide in result)
        starts = sorted(s for s, _ in result)
        assert starts == [6, 9, 12, 15, 18, 21, 24]
