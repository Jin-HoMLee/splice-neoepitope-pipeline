"""Tests for blastp_filter.py — filtering translated 9-mers against the human proteome."""

import csv
from pathlib import Path
from unittest.mock import patch

from blastp_filter import _parse_exact_matches, _write_peptide_fasta, blastp_filter


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fake_blastp(hits_content: str):
    """Return a _run_blastp replacement that writes hits_content to hits_path."""
    def _inner(fasta_path, db_prefix, hits_path, threads):
        Path(hits_path).write_text(hits_content)
    return _inner


def _make_peptides_tsv(tmp_path: Path, rows: list[dict]) -> Path:
    tsv = tmp_path / "peptides.tsv"
    with tsv.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["contig_key", "start_nt", "peptide"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)
    return tsv


def _read_tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------------------
# _write_peptide_fasta
# ---------------------------------------------------------------------------

class TestWritePeptideFasta:
    def test_basic(self, tmp_path):
        peptides = [
            {"contig_key": "c1", "start_nt": "0", "peptide": "ACDEFGHIK"},
            {"contig_key": "c2", "start_nt": "3", "peptide": "LMNPQRSTV"},
        ]
        fasta = tmp_path / "out.fasta"
        _write_peptide_fasta(peptides, fasta)
        lines = fasta.read_text().splitlines()
        assert lines == [">ACDEFGHIK", "ACDEFGHIK", ">LMNPQRSTV", "LMNPQRSTV"]

    def test_deduplicates_identical_sequences(self, tmp_path):
        peptides = [
            {"contig_key": "c1", "start_nt": "0", "peptide": "ACDEFGHIK"},
            {"contig_key": "c2", "start_nt": "3", "peptide": "ACDEFGHIK"},
            {"contig_key": "c3", "start_nt": "6", "peptide": "LMNPQRSTV"},
        ]
        fasta = tmp_path / "out.fasta"
        _write_peptide_fasta(peptides, fasta)
        lines = fasta.read_text().splitlines()
        assert lines.count(">ACDEFGHIK") == 1
        assert len(lines) == 4  # 2 unique peptides × 2 lines each

    def test_empty_input_writes_empty_file(self, tmp_path):
        fasta = tmp_path / "out.fasta"
        _write_peptide_fasta([], fasta)
        assert fasta.read_text() == ""


# ---------------------------------------------------------------------------
# _parse_exact_matches
# ---------------------------------------------------------------------------

class TestParseExactMatches:
    def test_100_identity_hit_included(self, tmp_path):
        hits = tmp_path / "hits.tsv"
        hits.write_text("ACDEFGHIK\tsp|P12345|PROT_HUMAN\t100.000\n")
        assert _parse_exact_matches(hits) == {"ACDEFGHIK": "sp|P12345|PROT_HUMAN"}

    def test_near_100_identity_excluded(self, tmp_path):
        hits = tmp_path / "hits.tsv"
        hits.write_text("ACDEFGHIK\tsp|P12345|PROT_HUMAN\t99.999\n")
        assert _parse_exact_matches(hits) == {}

    def test_partial_identity_excluded(self, tmp_path):
        hits = tmp_path / "hits.tsv"
        hits.write_text("ACDEFGHIK\tsp|P12345|PROT_HUMAN\t88.889\n")
        assert _parse_exact_matches(hits) == {}

    def test_99_identity_excluded(self, tmp_path):
        hits = tmp_path / "hits.tsv"
        hits.write_text("ACDEFGHIK\tsp|P12345|PROT_HUMAN\t99.000\n")
        assert _parse_exact_matches(hits) == {}

    def test_empty_file_returns_empty_dict(self, tmp_path):
        hits = tmp_path / "hits.tsv"
        hits.write_text("")
        assert _parse_exact_matches(hits) == {}

    def test_missing_file_returns_empty_dict(self, tmp_path):
        assert _parse_exact_matches(tmp_path / "nonexistent.tsv") == {}

    def test_keeps_first_accession_for_duplicate_peptide(self, tmp_path):
        hits = tmp_path / "hits.tsv"
        hits.write_text(
            "ACDEFGHIK\tsp|P11111|FIRST_HUMAN\t100.000\n"
            "ACDEFGHIK\tsp|P22222|SECOND_HUMAN\t100.000\n"
        )
        result = _parse_exact_matches(hits)
        assert result == {"ACDEFGHIK": "sp|P11111|FIRST_HUMAN"}

    def test_mixed_identity_only_exact_retained(self, tmp_path):
        hits = tmp_path / "hits.tsv"
        hits.write_text(
            "SELFPEPAA\tsp|P11111|SELF_HUMAN\t100.000\n"
            "PARTIALMATCH\tsp|P22222|PART_HUMAN\t95.000\n"
            "ANOTHSELF\tsp|P33333|ANOT_HUMAN\t100.000\n"
        )
        result = _parse_exact_matches(hits)
        assert set(result.keys()) == {"SELFPEPAA", "ANOTHSELF"}
        assert "PARTIALMATCH" not in result


# ---------------------------------------------------------------------------
# blastp_filter (mocking _run_blastp to avoid needing a real BLAST database)
# ---------------------------------------------------------------------------

class TestBlastpFilter:
    def test_splits_novel_and_excluded(self, tmp_path):
        peptides_tsv = _make_peptides_tsv(tmp_path, [
            {"contig_key": "c1", "start_nt": "0", "peptide": "SELFPEPAA"},
            {"contig_key": "c2", "start_nt": "3", "peptide": "NOVELPEPB"},
        ])
        hits = "SELFPEPAA\tsp|P11111|SELF_HUMAN\t100.000\n"

        with patch("blastp_filter._run_blastp", _fake_blastp(hits)):
            blastp_filter(
                peptides_tsv=peptides_tsv,
                db_prefix="fake_db",
                novel_tsv=tmp_path / "novel.tsv",
                excluded_tsv=tmp_path / "excluded.tsv",
                hits_path=tmp_path / "hits.tsv",
            )

        novel = _read_tsv(tmp_path / "novel.tsv")
        excluded = _read_tsv(tmp_path / "excluded.tsv")
        assert len(novel) == 1
        assert novel[0]["peptide"] == "NOVELPEPB"
        assert len(excluded) == 1
        assert excluded[0]["peptide"] == "SELFPEPAA"
        assert excluded[0]["matched_accession"] == "sp|P11111|SELF_HUMAN"

    def test_all_novel_when_no_hits(self, tmp_path):
        peptides_tsv = _make_peptides_tsv(tmp_path, [
            {"contig_key": "c1", "start_nt": "0", "peptide": "NOVELPEPB"},
        ])
        with patch("blastp_filter._run_blastp", _fake_blastp("")):
            blastp_filter(
                peptides_tsv=peptides_tsv,
                db_prefix="fake_db",
                novel_tsv=tmp_path / "novel.tsv",
                excluded_tsv=tmp_path / "excluded.tsv",
                hits_path=tmp_path / "hits.tsv",
            )
        assert len(_read_tsv(tmp_path / "novel.tsv")) == 1
        assert len(_read_tsv(tmp_path / "excluded.tsv")) == 0

    def test_all_excluded(self, tmp_path):
        peptides_tsv = _make_peptides_tsv(tmp_path, [
            {"contig_key": "c1", "start_nt": "0", "peptide": "SELFPEPAA"},
        ])
        hits = "SELFPEPAA\tsp|P11111|SELF_HUMAN\t100.000\n"
        with patch("blastp_filter._run_blastp", _fake_blastp(hits)):
            blastp_filter(
                peptides_tsv=peptides_tsv,
                db_prefix="fake_db",
                novel_tsv=tmp_path / "novel.tsv",
                excluded_tsv=tmp_path / "excluded.tsv",
                hits_path=tmp_path / "hits.tsv",
            )
        assert len(_read_tsv(tmp_path / "novel.tsv")) == 0
        assert len(_read_tsv(tmp_path / "excluded.tsv")) == 1

    def test_output_files_have_correct_headers(self, tmp_path):
        peptides_tsv = _make_peptides_tsv(tmp_path, [
            {"contig_key": "c1", "start_nt": "0", "peptide": "SELFPEPAA"},
            {"contig_key": "c2", "start_nt": "3", "peptide": "NOVELPEPB"},
        ])
        hits = "SELFPEPAA\tsp|P11111|SELF_HUMAN\t100.000\n"
        with patch("blastp_filter._run_blastp", _fake_blastp(hits)):
            blastp_filter(
                peptides_tsv=peptides_tsv,
                db_prefix="fake_db",
                novel_tsv=tmp_path / "novel.tsv",
                excluded_tsv=tmp_path / "excluded.tsv",
                hits_path=tmp_path / "hits.tsv",
            )
        novel_header = (tmp_path / "novel.tsv").read_text().splitlines()[0]
        excl_header = (tmp_path / "excluded.tsv").read_text().splitlines()[0]
        assert novel_header == "contig_key\tstart_nt\tpeptide"
        assert excl_header == "contig_key\tstart_nt\tpeptide\tmatched_accession"

    def test_same_peptide_from_two_contigs_both_excluded(self, tmp_path):
        """The same peptide sequence arising from two different junctions should
        both be excluded if the sequence matches the proteome."""
        peptides_tsv = _make_peptides_tsv(tmp_path, [
            {"contig_key": "c1", "start_nt": "0", "peptide": "SELFPEPAA"},
            {"contig_key": "c2", "start_nt": "3", "peptide": "SELFPEPAA"},
        ])
        hits = "SELFPEPAA\tsp|P11111|SELF_HUMAN\t100.000\n"
        with patch("blastp_filter._run_blastp", _fake_blastp(hits)):
            blastp_filter(
                peptides_tsv=peptides_tsv,
                db_prefix="fake_db",
                novel_tsv=tmp_path / "novel.tsv",
                excluded_tsv=tmp_path / "excluded.tsv",
                hits_path=tmp_path / "hits.tsv",
            )
        excluded = _read_tsv(tmp_path / "excluded.tsv")
        assert len(excluded) == 2
        assert all(row["peptide"] == "SELFPEPAA" for row in excluded)
        assert excluded[0]["contig_key"] == "c1"
        assert excluded[1]["contig_key"] == "c2"
