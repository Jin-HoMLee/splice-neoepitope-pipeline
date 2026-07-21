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


def _parse_fasta(fasta_path):
    records = []
    header = None
    parts = []
    with Path(fasta_path).open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(parts)))
                header = line[1:]
                parts = []
            else:
                parts.append(line)
    if header is not None:
        records.append((header, "".join(parts)))
    return records


def _orf_header(contig_header, frame):
    # contig header: {junction_id}|{chrom}:{start}-{end}:{strand}|{sample_type}
    fields = contig_header.split("|")
    junction_id, coords = fields[0], fields[1]
    return f"{junction_id}|{frame}|{coords}"


def emit_orf_fasta(contigs_fasta, output_fasta, flank_nt=30, min_peptide_len=8,
                   frames=(0, 1, 2)):
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    records = _parse_fasta(contigs_fasta)

    n_written = n_drop_x = n_drop_short = n_none = 0
    with output_fasta.open("w") as out:
        for header, seq in records:
            seq = seq.upper()
            if len(seq) != 2 * flank_nt:
                log.warning("Skipping %s: length %d != %d", header, len(seq),
                            2 * flank_nt)
                continue
            for frame in frames:
                stretch = crossing_orf_stretch(seq, frame, breakpoint_nt=flank_nt)
                if stretch is None:
                    n_none += 1
                    continue
                if "X" in stretch:
                    n_drop_x += 1
                    continue
                if len(stretch) < min_peptide_len:
                    n_drop_short += 1
                    continue
                out.write(f">{_orf_header(header, frame)}\n{stretch}\n")
                n_written += 1

    log.info(
        "ORF FASTA: %d written, %d dropped (X), %d dropped (<%d aa), "
        "%d frames with no crossing stretch",
        n_written, n_drop_x, n_drop_short, min_peptide_len, n_none,
    )


def _snakemake_main():
    log_file = snakemake.log[0]  # type: ignore[name-defined]  # noqa: F821
    logging.getLogger().addHandler(logging.FileHandler(log_file))
    emit_orf_fasta(
        contigs_fasta=snakemake.input.contigs_fasta,  # type: ignore[name-defined]  # noqa: F821
        output_fasta=snakemake.output.orf_fasta,  # type: ignore[name-defined]  # noqa: F821
        flank_nt=snakemake.params.flank_nt,  # type: ignore[name-defined]  # noqa: F821
        min_peptide_len=snakemake.params.min_peptide_len,  # type: ignore[name-defined]  # noqa: F821
    )


def _cli_main():
    p = argparse.ArgumentParser(description="Emit breakpoint-crossing ORF-stretch FASTA.")
    p.add_argument("--contigs-fasta", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--flank-nt", type=int, default=30)
    p.add_argument("--min-peptide-len", type=int, default=8)
    args = p.parse_args()
    emit_orf_fasta(args.contigs_fasta, args.output, args.flank_nt, args.min_peptide_len)


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]  # noqa: F821
        _snakemake_main()
    except NameError:
        _cli_main()
