"""Tests for translate_peptides.py — in-silico translation logic."""

from translate_peptides import translate_contig


class TestTranslateContig:
    def test_frame0_basic(self):
        # ATG = Met, GGG = Gly, ... simple 50 nt contig, no stop codon
        # Use a known sequence: 48 nt → 16 aa in frame 0
        seq = "ATG" * 16  # 48 nt, all Met
        results = translate_contig(seq, reading_frames=[0], peptide_length=16)
        assert len(results) == 1
        frame, peptide = results[0]
        assert frame == 0
        assert peptide == "M" * 16

    def test_all_three_frames(self):
        # GCC (Ala) repeated: no stop codons in any reading frame
        # Frame 0: GCC GCC... = Ala; Frame 1: CCG CCG... = Pro; Frame 2: CGC CGC... = Arg
        seq = "GCC" * 17  # 51 nt
        results = translate_contig(seq, reading_frames=[0, 1, 2], peptide_length=16)
        frames_returned = {r[0] for r in results}
        assert frames_returned == {0, 1, 2}

    def test_stop_codon_truncates_peptide(self):
        # TAA is a stop codon; peptide should be truncated before it
        # Frame 0: ATG ATG TAA ... → "MM" (2 aa) — too short, filtered out
        # Use a longer prefix: 8× ATG then TAA
        seq = "ATG" * 8 + "TAA" + "ATG" * 5
        results = translate_contig(seq, reading_frames=[0], peptide_length=16)
        assert len(results) == 1
        frame, peptide = results[0]
        assert "*" not in peptide
        assert len(peptide) == 8

    def test_immediate_stop_codon_excluded(self):
        # TAA at position 0 in frame 0 → no peptide for that frame
        seq = "TAA" + "ATG" * 20
        results = translate_contig(seq, reading_frames=[0], peptide_length=16)
        assert len(results) == 0

    def test_short_peptide_filtered_out(self):
        # Only 6 aa before stop → below _MIN_PEPTIDE_LENGTH (8), excluded
        seq = "ATG" * 6 + "TAA" + "ATG" * 10
        results = translate_contig(seq, reading_frames=[0], peptide_length=16)
        assert len(results) == 0

    def test_peptide_trimmed_to_target_length(self):
        seq = "ATG" * 20  # 60 nt → 20 aa in frame 0, trimmed to 16
        results = translate_contig(seq, reading_frames=[0], peptide_length=16)
        assert len(results) == 1
        assert len(results[0][1]) == 16

    def test_empty_sequence_returns_nothing(self):
        results = translate_contig("", reading_frames=[0, 1, 2], peptide_length=16)
        assert results == []
