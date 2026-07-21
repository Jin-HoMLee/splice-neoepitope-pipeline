# workflow/tests/test_orf_fasta_from_contigs.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from orf_fasta_from_contigs import crossing_orf_stretch


def _codon_seq(aas, frame=0, upstream_codons=10):
    # Build a 60-nt contig whose frame-`frame` translation is `aas`.
    # One unambiguous codon per residue; breakpoint at nt 30.
    table = {"M": "ATG", "K": "AAA", "L": "CTG", "F": "TTT", "*": "TAA", "P": "CCC"}
    body = "".join(table[a] for a in aas)
    pad = "G" * frame
    seq = (pad + body)
    return (seq + "C" * 60)[:60]


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


from orf_fasta_from_contigs import _orf_header, emit_orf_fasta


def test_orf_header_reformat():
    h = "JUNC1|chr7:100-200:+|tumor"
    assert _orf_header(h, 2) == "JUNC1|2|chr7:100-200:+"


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
