"""Tests for run_mhcflurry.py — 9-mer generation and junction-spanning filter."""

from run_mhcflurry import _generate_9mers, _parse_frame_offset


class TestParseFrameOffset:
    def test_frame1(self):
        header = "chr22:100:200:+|chr22:100-200:+|Primary Tumor|frame1"
        assert _parse_frame_offset(header) == 0

    def test_frame2(self):
        header = "chr22:100:200:+|chr22:100-200:+|Primary Tumor|frame2"
        assert _parse_frame_offset(header) == 1

    def test_frame3(self):
        header = "chr22:100:200:+|chr22:100-200:+|Primary Tumor|frame3"
        assert _parse_frame_offset(header) == 2

    def test_no_frame_token(self):
        header = "some_header_without_frame"
        assert _parse_frame_offset(header) is None

    def test_invalid_frame_token_returns_none(self):
        for bad in ("frame0", "frame4", "frame10"):
            header = f"chr22:100:200:+|{bad}"
            assert _parse_frame_offset(header) is None, f"expected None for {bad}"


class TestGenerate9mers:
    # A 16-mer with no stop codons or invalid characters
    SEQ = "ACDEFGHIKLMNPQRS"  # 16 standard amino acids

    def test_no_filter_returns_all_8_positions(self):
        # Without frame_offset, all 8 windows (16 - 9 + 1) are returned
        result = _generate_9mers(self.SEQ, window_size=9, frame_offset=None)
        assert len(result) == 8

    def test_positions_are_1_indexed(self):
        result = _generate_9mers(self.SEQ, window_size=9, frame_offset=None)
        positions = [pos for pos, _ in result]
        assert positions == [1, 2, 3, 4, 5, 6, 7, 8]

    def test_junction_filter_frame0_drops_position1_and_last(self):
        # frame_offset=0, upstream_nt=26: valid start_nt range is [2, 23]
        # position i (0-indexed): start_nt = 0 + i*3
        # i=0 → start_nt=0  (< 2) → excluded
        # i=1 → start_nt=3  ✓
        # i=7 → start_nt=21 ✓
        # i=8 would be start_nt=24 (> 23) → excluded, but 16-mer only has 8 windows
        result = _generate_9mers(self.SEQ, window_size=9, frame_offset=0, upstream_nt=26)
        positions = [pos for pos, _ in result]
        assert 1 not in positions  # start_nt=0, purely upstream
        assert 2 in positions      # start_nt=3, spans junction

    def test_junction_filter_frame1_drops_position1(self):
        # frame_offset=1: start_nt = 1 + i*3
        # i=0 → start_nt=1 (< 2) → excluded
        # i=1 → start_nt=4 ✓
        result = _generate_9mers(self.SEQ, window_size=9, frame_offset=1, upstream_nt=26)
        positions = [pos for pos, _ in result]
        assert 1 not in positions
        assert 2 in positions

    def test_junction_filter_frame2_includes_position1(self):
        # frame_offset=2: start_nt = 2 + i*3
        # i=0 → start_nt=2, which equals min_start=2 → included
        result = _generate_9mers(self.SEQ, window_size=9, frame_offset=2, upstream_nt=26)
        positions = [pos for pos, _ in result]
        assert 1 in positions

    def test_stop_codon_excluded(self):
        seq = "ACDEFGHIK*MNPQRS"
        result = _generate_9mers(seq, window_size=9, frame_offset=None)
        nmers = [nmer for _, nmer in result]
        assert all("*" not in nmer for nmer in nmers)

    def test_x_excluded(self):
        seq = "ACDEFGHIKXMNPQRS"
        result = _generate_9mers(seq, window_size=9, frame_offset=None)
        nmers = [nmer for _, nmer in result]
        assert all("X" not in nmer for nmer in nmers)
