"""Tests for assemble_contigs.py — contig assembly helpers."""

import pandas as pd

from assemble_contigs import (
    _build_downstream_bed,
    _build_upstream_bed,
    _has_soft_clip,
    _parse_fasta,
    _write_zero_stats,
    assemble_contigs,
)


class TestHasSoftClip:
    def test_all_uppercase_returns_false(self):
        assert _has_soft_clip("ATCGATCGATCG") is False

    def test_lowercase_bases_return_true(self):
        assert _has_soft_clip("ATCGatcg") is True

    def test_all_lowercase_returns_true(self):
        assert _has_soft_clip("atcgatcg") is True

    def test_empty_string_returns_false(self):
        assert _has_soft_clip("") is False


class TestBuildUpstreamBed:
    def _make_junctions(self):
        return pd.DataFrame({
            "chrom": ["chr22"],
            "start": [200],   # 0-based junction start
            "end":   [300],
            "strand": ["+"],
            "junction_id": ["chr22:201:300:+"],
        })

    def test_upstream_end_equals_junction_start(self):
        df = self._make_junctions()
        bed = _build_upstream_bed(df, upstream_nt=26)
        assert bed.iloc[0]["end"] == 200

    def test_upstream_start_is_junction_start_minus_nt(self):
        df = self._make_junctions()
        bed = _build_upstream_bed(df, upstream_nt=26)
        assert bed.iloc[0]["start"] == 174

    def test_upstream_start_clipped_at_zero(self):
        df = pd.DataFrame({
            "chrom": ["chr22"], "start": [10], "end": [200],
            "strand": ["+"], "junction_id": ["chr22:11:200:+"],
        })
        bed = _build_upstream_bed(df, upstream_nt=26)
        assert bed.iloc[0]["start"] == 0


class TestBuildDownstreamBed:
    def _make_junctions(self):
        return pd.DataFrame({
            "chrom": ["chr22"],
            "start": [200],
            "end":   [300],   # 0-based junction end
            "strand": ["+"],
            "junction_id": ["chr22:201:300:+"],
        })

    def test_downstream_start_equals_junction_end(self):
        df = self._make_junctions()
        bed = _build_downstream_bed(df, downstream_nt=24)
        assert bed.iloc[0]["start"] == 300

    def test_downstream_end_is_junction_end_plus_nt(self):
        df = self._make_junctions()
        bed = _build_downstream_bed(df, downstream_nt=24)
        assert bed.iloc[0]["end"] == 324


class TestParseFasta:
    def test_parses_single_record(self, tmp_path):
        fa = tmp_path / "test.fa"
        fa.write_text(">header1\nATCGATCG\n")
        result = _parse_fasta(fa)
        assert result == {"header1": "ATCGATCG"}

    def test_parses_multiple_records(self, tmp_path):
        fa = tmp_path / "test.fa"
        fa.write_text(">seq1\nAAAA\n>seq2\nCCCC\n")
        result = _parse_fasta(fa)
        assert result["seq1"] == "AAAA"
        assert result["seq2"] == "CCCC"

    def test_strips_bedtools_coordinate_suffix(self, tmp_path):
        # bedtools getfasta -name appends ::chrom:start-end to the header
        fa = tmp_path / "test.fa"
        fa.write_text(">chr22:201:300:+::chr22:174-200\nATCG\n")
        result = _parse_fasta(fa)
        assert "chr22:201:300:+" in result

    def test_empty_file_returns_empty_dict(self, tmp_path):
        fa = tmp_path / "empty.fa"
        fa.write_text("")
        assert _parse_fasta(fa) == {}


class TestEmptyPipelineEarlyReturns:
    """Issue #215 follow-up: when the pipeline produces zero candidates, the
    early-return paths must still write a zero-count stats TSV so the
    cross-step aggregator does not fail with missing input.
    """

    _NOVEL_COLUMNS = [
        "junction_id", "chrom", "start", "end", "strand",
        "junction_origin", "sample_id", "sample_type",
    ]

    def _expected_zero_stats(self):
        return pd.DataFrame(
            [
                {"category": "contigs_written", "count": 0},
                {"category": "skipped_softclip", "count": 0},
                {"category": "skipped_length", "count": 0},
            ],
            columns=["category", "count"],
        )

    def test_write_zero_stats_helper_emits_three_rows(self, tmp_path):
        out = tmp_path / "stats.tsv"
        _write_zero_stats(out)
        got = pd.read_csv(out, sep="\t")
        pd.testing.assert_frame_equal(got, self._expected_zero_stats())

    def test_write_zero_stats_helper_no_path_is_noop(self, tmp_path):
        # No exception, no file created
        _write_zero_stats(None)
        assert list(tmp_path.iterdir()) == []

    def test_empty_novel_junctions_writes_zero_stats(self, tmp_path):
        novel_tsv = tmp_path / "novel.tsv"
        pd.DataFrame(columns=self._NOVEL_COLUMNS).to_csv(
            novel_tsv, sep="\t", index=False
        )
        fasta_out = tmp_path / "contigs.fa"
        stats_out = tmp_path / "contig_stats.tsv"

        assemble_contigs(
            novel_junctions_tsv=novel_tsv,
            genome_fasta=tmp_path / "fake.fa",  # never read in early-return path
            output_fasta=fasta_out,
            stats_output_path=stats_out,
        )

        assert fasta_out.exists()
        got = pd.read_csv(stats_out, sep="\t")
        pd.testing.assert_frame_equal(got, self._expected_zero_stats())

    def test_only_normal_shared_writes_zero_stats(self, tmp_path):
        novel_tsv = tmp_path / "novel.tsv"
        pd.DataFrame({
            "junction_id": ["j1", "j2"],
            "chrom": ["chr22", "chr22"],
            "start": [100, 200],
            "end": [200, 300],
            "strand": ["+", "+"],
            "junction_origin": ["normal_shared", "normal_shared"],
            "sample_id": ["s1", "s1"],
            "sample_type": ["tumor", "tumor"],
        }).to_csv(novel_tsv, sep="\t", index=False)
        fasta_out = tmp_path / "contigs.fa"
        stats_out = tmp_path / "contig_stats.tsv"

        assemble_contigs(
            novel_junctions_tsv=novel_tsv,
            genome_fasta=tmp_path / "fake.fa",
            output_fasta=fasta_out,
            stats_output_path=stats_out,
        )

        assert fasta_out.exists()
        got = pd.read_csv(stats_out, sep="\t")
        pd.testing.assert_frame_equal(got, self._expected_zero_stats())


class TestSoftMaskedReferenceAssembles:
    """A soft-masked (lower-case) genome reference must NOT cause contigs to be
    dropped. Lower-case in a genome FASTA is repeat soft-masking, not an
    alignment soft-clip; contigs here come from ``bedtools getfasta`` on the
    reference, so they never contain alignment soft-clips. Production (GENCODE)
    is upper-case, but a UCSC-style soft-masked reference (e.g. the chr22 test
    fixture) was silently dropping every contig. The assembled contig must be
    emitted and upper-cased.
    """

    def _make_tumor_exclusive_junction(self, path):
        pd.DataFrame({
            "junction_id": ["chr22:201:300:+"],
            "chrom": ["chr22"], "start": [200], "end": [300],
            "strand": ["+"], "junction_origin": ["tumor_exclusive"],
            "sample_id": ["s1"], "sample_type": ["tumor"],
        }).to_csv(path, sep="\t", index=False)

    def test_softmasked_sequence_still_emits_uppercase_contig(self, tmp_path, monkeypatch):
        novel_tsv = tmp_path / "novel.tsv"
        self._make_tumor_exclusive_junction(novel_tsv)

        # Stub the external bedtools getfasta (unavailable in the test venv) to
        # emit soft-masked (lower-case) sequence, as a UCSC-masked reference would.
        def fake_getfasta(bed_path, genome_fasta, out_fa):
            with open(out_fa, "w") as fh:
                fh.write(">chr22:201:300:+\n" + "acgt" * 100 + "\n")
        monkeypatch.setattr("assemble_contigs._run_bedtools_getfasta", fake_getfasta)

        fasta_out = tmp_path / "contigs.fa"
        assemble_contigs(
            novel_junctions_tsv=novel_tsv,
            genome_fasta=tmp_path / "fake.fa",
            output_fasta=fasta_out,
            stats_output_path=tmp_path / "stats.tsv",
        )

        contigs = fasta_out.read_text()
        assert ">chr22:201:300:+" in contigs, "soft-masked contig was dropped"
        seq_lines = [ln for ln in contigs.splitlines() if ln and not ln.startswith(">")]
        assert seq_lines, "no contig sequence emitted"
        assert all(ln.isupper() for ln in seq_lines), "contig not upper-cased"


class TestMinusStrandArmOrder:
    """Issue #1278: a minus-strand junction contig must concatenate the two exon
    arms in TRANSCRIPT order, not genomic order. ``bedtools getfasta -s`` reverse-
    complements each arm but does NOT reorder them, so on the minus strand the
    correct 5'->3' contig is downstream-arm + upstream-arm. The old code emitted
    upstream + downstream unconditionally, fabricating the junction boundary for
    every minus-strand junction (~half of all junctions), invisible to CI (dry-run
    never executes) and to the other tests here (all strand '+').
    """

    def _contig_for_strand(self, tmp_path, monkeypatch, strand):
        junction_id = f"chr22:201:300:{strand}"
        novel_tsv = tmp_path / "novel.tsv"
        pd.DataFrame({
            "junction_id": [junction_id],
            "chrom": ["chr22"], "start": [200], "end": [300],
            "strand": [strand], "junction_origin": ["tumor_exclusive"],
            "sample_id": ["s1"], "sample_type": ["tumor"],
        }).to_csv(novel_tsv, sep="\t", index=False)

        # Stub getfasta: upstream arm = all-A, downstream arm = all-C, so the
        # ORDER of the two arms in the emitted contig is unambiguous. Each marker
        # stands in for the getfasta -s output (already transcript-oriented per
        # arm); only the concatenation order is under test.
        def fake_getfasta(bed_path, genome_fasta, out_fa, strand=True):
            base = "A" if "upstream" in str(bed_path) else "C"
            with open(out_fa, "w") as fh:
                fh.write(f">{junction_id}\n{base * 100}\n")
        monkeypatch.setattr("assemble_contigs._run_bedtools_getfasta", fake_getfasta)

        fasta_out = tmp_path / "contigs.fa"
        assemble_contigs(
            novel_junctions_tsv=novel_tsv,
            genome_fasta=tmp_path / "fake.fa",
            output_fasta=fasta_out,
            stats_output_path=tmp_path / "stats.tsv",
        )
        seqs = [ln for ln in fasta_out.read_text().splitlines()
                if ln and not ln.startswith(">")]
        assert seqs, "no contig emitted"
        return seqs[0]

    def test_minus_strand_orders_downstream_arm_first(self, tmp_path, monkeypatch):
        # Transcript 5'->3' on the minus strand is downstream-arm (C) then
        # upstream-arm (A). The bug emits A*27 + C*27 (genomic order) instead.
        contig = self._contig_for_strand(tmp_path, monkeypatch, "-")
        assert contig == "C" * 27 + "A" * 27, (
            "minus-strand contig must be downstream+upstream (transcript order); "
            f"got {contig[:6]}...{contig[-6:]}"
        )

    def test_plus_strand_orders_upstream_arm_first(self, tmp_path, monkeypatch):
        # Guard: the fix must NOT unconditionally swap. Plus strand keeps
        # upstream-arm (A) then downstream-arm (C).
        contig = self._contig_for_strand(tmp_path, monkeypatch, "+")
        assert contig == "A" * 27 + "C" * 27


class TestMinusStrandIntegrationThroughTranslate:
    """Issue #1278 integration guard: a minus-strand junction carried through
    assemble_contigs THEN translate_peptides must yield the biologically-correct
    junction-spanning peptide, not just correct arm order.

    This is the guard the unit test cannot give. translate_peptides translates
    forward frames only and assumes the contig is already in mRNA sense with the
    breakpoint at ``upstream_nt`` - the exact assumption that made the bug silent.
    Every prior fixture across both scripts used strand '+', so this whole class
    of error (wrong peptides from a fabricated minus-strand junction) was untested
    end to end. Recommended by the PR #1286 bot review.
    """

    def test_minus_strand_yields_transcript_correct_peptide(self, tmp_path, monkeypatch):
        from translate_peptides import extract_spanning_peptides

        junction_id = "chr22:201:300:-"
        novel_tsv = tmp_path / "novel.tsv"
        pd.DataFrame({
            "junction_id": [junction_id],
            "chrom": ["chr22"], "start": [200], "end": [300],
            "strand": ["-"], "junction_origin": ["tumor_exclusive"],
            "sample_id": ["s1"], "sample_type": ["tumor"],
        }).to_csv(novel_tsv, sep="\t", index=False)

        # getfasta -s returns each arm already in transcript sense. Codon-clean
        # markers make the arm ORDER legible in the translated peptide:
        #   upstream arm   = GCT*9 -> Ala ('A')
        #   downstream arm = GGT*9 -> Gly ('G')
        # Correct minus-strand contig is dn+up, so frame-0 translation is G*9 then
        # A*9 and junction-spanning 8-mers read G...A. The buggy up+dn order would
        # translate to A*9 then G*9 (A...G) - a different, wrong peptide set.
        def fake_getfasta(bed_path, genome_fasta, out_fa, strand=True):
            arm = "GCT" * 9 if "upstream" in str(bed_path) else "GGT" * 9
            with open(out_fa, "w") as fh:
                fh.write(f">{junction_id}\n{arm}\n")
        monkeypatch.setattr("assemble_contigs._run_bedtools_getfasta", fake_getfasta)

        fasta_out = tmp_path / "contigs.fa"
        assemble_contigs(
            novel_junctions_tsv=novel_tsv,
            genome_fasta=tmp_path / "fake.fa",
            output_fasta=fasta_out,
            stats_output_path=tmp_path / "stats.tsv",
        )
        contig = [ln for ln in fasta_out.read_text().splitlines()
                  if ln and not ln.startswith(">")][0]

        peptides = {
            pep for _, pep in extract_spanning_peptides(
                contig, upstream_nt=27, peptide_lengths=[8, 9, 10], reading_frames=[0]
            )
        }

        # Transcript order (dn+up) -> a G-then-A spanning peptide.
        assert "GGGGGGGA" in peptides, (
            "minus-strand junction did not translate to the transcript-correct "
            f"spanning peptide; got {sorted(peptides)}"
        )
        # The buggy genomic order (up+dn) would produce an A-then-G peptide; it
        # must never appear from a correctly-oriented minus-strand contig.
        assert "AAAAAAAG" not in peptides, (
            "minus-strand arms concatenated in genomic (buggy) order"
        )
