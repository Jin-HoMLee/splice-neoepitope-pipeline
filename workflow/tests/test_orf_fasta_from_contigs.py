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
