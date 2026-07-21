# workflow/tests/test_orf_fasta_from_contigs.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from orf_fasta_from_contigs import (
    crossing_orf_stretch,
    _orf_header,
    emit_orf_fasta,
    _parse_fasta,
)


def test_clean_crossing_stretch_frame0():
    # 20 codons, no stop: one stretch spanning nt [0,60), crosses 30.
    contig = "".join(["AAA"] * 20)  # all Lys, no stop
    assert crossing_orf_stretch(contig, 0) == "K" * 20


def test_stop_at_breakpoint_yields_none():
    # Stop codon occupies codon index 10 -> nt [30,33). Upstream run ends at
    # nt 30 (entirely upstream), downstream run starts at nt 33: neither crosses.
    contig = "".join(["AAA"] * 10 + ["TAA"] + ["AAA"] * 9)
    assert crossing_orf_stretch(contig, 0) is None


def test_upstream_only_stretch_dropped():
    # Stop at codon 5 (nt 15) then run to end: the surviving run nt [18,60)
    # starts downstream? No - it starts at 18 < 30 < 60, so it DOES cross.
    # Construct an upstream-only run: stop at codon 12 (nt 36), first run
    # nt [0,36) crosses; make first run end BEFORE 30 with a stop at codon 8.
    contig = "".join(["AAA"] * 8 + ["TAA"] + ["AAA"] * 11)
    # first run nt [0,24) entirely upstream (24 < 30) -> dropped;
    # second run nt [27,60) crosses (27 < 30 < 60) -> returned.
    assert crossing_orf_stretch(contig, 0) == "K" * 11


def test_frame_offset_shifts_breakpoint_math():
    # frame 1: codon i covers nt [1+3i, 1+3i+3); 19 codons over nt [1,58).
    contig = "G" + "".join(["AAA"] * 19)
    assert crossing_orf_stretch(contig, 1) == "K" * 19


def test_orf_header_reformat():
    h = "JUNC1|chr7:100-200:+|tumor"
    assert _orf_header(h, 2) == "JUNC1|2|chr7:100-200:+"


def test_orf_header_rejects_malformed():
    import pytest

    with pytest.raises(ValueError, match="malformed contig header"):
        _orf_header("no-pipe-here", 0)


def _write_fasta(path, records):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")


def test_emit_drops_x_and_short_and_writes_headers(tmp_path):
    contigs = tmp_path / "contigs.fa"
    out = tmp_path / "orf.fa"
    # J1: clean 20-Lys crossing stretch in frame 0.
    clean = "".join(["AAA"] * 20)
    # J2: contains an ambiguous base (N) inside a crossing run -> 'X' -> dropped.
    ambig = "".join(["AAA"] * 9) + "ANA" + "".join(["AAA"] * 10)
    _write_fasta(contigs, [
        ("J1|chr1:10-20:+|tumor", clean),
        ("J2|chr2:30-40:-|tumor", ambig),
    ])
    emit_orf_fasta(contigs, out, flank_nt=30, min_peptide_len=8)
    text = out.read_text()
    assert ">J1|0|chr1:10-20:+" in text
    assert "K" * 20 in text
    assert "J2|" not in text  # X-containing stretch dropped
    # headers unique
    headers = [ln for ln in text.splitlines() if ln.startswith(">")]
    assert len(headers) == len(set(headers))


def test_emit_drops_short_crossing_stretch(tmp_path):
    contigs = tmp_path / "contigs.fa"
    out = tmp_path / "orf.fa"
    # J1 frame 0 translates to KKKKKK*KKKK*KKKKKKKK (20 codons):
    #   codons 0-5   (nt [0,18))  -> run "KKKKKK", entirely upstream, dropped
    #   codon 6      (nt [18,21)) -> stop
    #   codons 7-10  (nt [21,33)) -> run "KKKK", crosses breakpoint nt 30
    #                                 (21 < 30 < 33) but is only 4 aa
    #   codon 11     (nt [33,36)) -> stop
    #   codons 12-19 (nt [36,60)) -> run "KKKKKKKK", entirely downstream
    # So the ONLY crossing stretch is the 4-aa run, below min_peptide_len=8
    # -> dropped.
    short_crossing = (
        "".join(["AAA"] * 6) + "TAA" + "".join(["AAA"] * 4) + "TAA"
        + "".join(["AAA"] * 8)
    )
    assert len(short_crossing) == 60
    # J2 frame 0: clean 20-Lys crossing stretch (matched-pair control - same
    # min_peptide_len=8 must NOT drop this long stretch, else the assertion
    # on J1 would be vacuous).
    long_crossing = "".join(["AAA"] * 20)
    _write_fasta(contigs, [
        ("J1|chr1:10-20:+|tumor", short_crossing),
        ("J2|chr2:30-40:-|tumor", long_crossing),
    ])
    emit_orf_fasta(contigs, out, flank_nt=30, min_peptide_len=8)
    text = out.read_text()
    assert ">J1|0|chr1:10-20:+" not in text  # short crossing stretch dropped
    assert ">J2|0|chr2:30-40:-" in text  # long crossing stretch kept
    assert "K" * 20 in text


def test_emit_skips_wrong_length_contig(tmp_path):
    contigs = tmp_path / "contigs.fa"
    out = tmp_path / "orf.fa"
    wrong_length = "A" * 45  # != 2 * flank_nt (60) -> must be skipped
    valid = "".join(["AAA"] * 20)  # 60 nt, clean frame-0 crossing stretch
    _write_fasta(contigs, [
        ("J1|chr1:10-20:+|tumor", wrong_length),
        ("J2|chr2:30-40:-|tumor", valid),
    ])
    emit_orf_fasta(contigs, out, flank_nt=30, min_peptide_len=8)
    text = out.read_text()
    assert "J1|" not in text  # wrong-length contig skipped, no record at all
    assert ">J2|0|chr2:30-40:-" in text  # valid contig still written


_VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def test_multichr_capability_and_structural_validity(tmp_path):
    data = Path(__file__).parent / "data" / "orf_fasta"
    out = tmp_path / "orf.fa"
    emit_orf_fasta(data / "contigs_multichr.fa", out, flank_nt=30, min_peptide_len=8)
    records = _parse_fasta(out)  # reuse the emitter's parser

    # AC-1: entries from all three chromosomes (not chr22-scoped).
    chroms = {h.split("|")[2].split(":")[0] for h, _ in records}
    assert {"chr1", "chr7", "chrX"} <= chroms

    # Structural: unique headers, valid alphabet, no stop/X.
    headers = [h for h, _ in records]
    assert len(headers) == len(set(headers))
    for _, seq in records:
        assert set(seq) <= _VALID_AA

    # Concatenation with the canonical proteome yields no duplicate IDs.
    canon = _parse_fasta(data / "canonical_tiny.fasta")
    all_ids = [h.split()[0] for h, _ in records + canon]
    assert len(all_ids) == len(set(all_ids))
