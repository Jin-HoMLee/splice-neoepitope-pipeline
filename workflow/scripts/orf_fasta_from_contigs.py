#!/usr/bin/env python3
"""orf_fasta_from_contigs.py - emit breakpoint-crossing ORF stretches as FASTA.

For each junction contig (30 nt upstream + 30 nt downstream of the breakpoint),
translate all three reading frames, split each translation at stop codons, and
keep the single stop-free ORF stretch that crosses the breakpoint (nt 30). Each
kept stretch is written as a FASTA mini-protein for a nonspecific MS search.
"""

import argparse
import logging
from pathlib import Path

from Bio.Seq import Seq

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

_CODON = 3


def crossing_orf_stretch(contig_seq: str, frame: int, breakpoint_nt: int = 30) -> str | None:
    """Return the amino-acid string of the stop-free stretch crossing the
    breakpoint in this frame, or None.

    Codon i (0-based) covers nt [frame + 3i, frame + 3i + 3). A stretch
    spanning codons [c0, c1) covers nt [frame + 3*c0, frame + 3*c1) and
    crosses iff frame + 3*c0 < breakpoint_nt < frame + 3*c1.
    """
    usable = contig_seq[frame:]
    usable = usable[: (len(usable) // _CODON) * _CODON]
    aa = str(Seq(usable).translate())
    c0 = 0
    for i, residue in enumerate(list(aa) + ["*"]):  # sentinel flushes last run
        if residue == "*":
            if i > c0:
                nt_start = frame + _CODON * c0
                nt_end = frame + _CODON * i
                if nt_start < breakpoint_nt < nt_end:
                    return aa[c0:i]
            c0 = i + 1
    return None
