#!/usr/bin/env python3
"""regenerate.py — (re)build the hermetic annotate-flag canary fixture (Issue #377).

Emits three tiny, committed files the ``annotate-flag-canary`` CI job consumes:

- ``genome.fa``      — a synthetic single-contig FASTA (``chrT``), long enough to
                       cover every fixture coordinate. Sequence content is
                       irrelevant: ``regtools junctions annotate``'s
                       ``known_junction`` flag is annotation-based, not motif-based
                       (verified empirically — a non-canonical motif still scores
                       ``known_junction = 1`` when the junction is annotated).
- ``annotation.gtf`` — synthetic multi-exon genes. Their consecutive exon
                       boundaries are the "known" junctions, by the *same* rule
                       both ``build_reference_junctions.py`` (home-rolled source of
                       truth) and ``regtools junctions annotate`` use.
- ``junctions.bed12`` — the anchor-outer BED12 that ``regtools junctions extract``
                       would emit, for every known junction plus a few fabricated
                       *novel* ones (a real donor paired with an off-annotation
                       acceptor). Novel junctions exercise the ``known_junction = 0``
                       direction so the canary tests both, not just the annotated case.

Why hermetic/synthetic rather than a chr22 slice: the canary's job is to catch a
**coordinate-semantics** desync between the home-rolled flag and regtools (the
Issue #370 anchor-outer bug class). That is a property of the coordinate math,
not of real biology — so a tiny self-consistent fixture tests it precisely, with
no 400 MB GENCODE download and no gitignored chr22 data in CI.

Anchor sizes are varied across junctions so the fixture exercises the
``blockSizes``/``blockStarts`` math non-trivially (a hardcoded-anchor assumption
would break here).

Run from the repo root (snakemake env: has the deps + a built regtools to self-verify):

    conda activate snakemake
    python resources/test/annotate_canary/regenerate.py
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]

CHROM = "chrT"
EXON_LEN = 100
# Deterministic gene layout: (n_exons, intron_len, strand) cycled across genes.
GENE_SPECS = [
    (3, 200, "+"),
    (2, 300, "-"),
    (4, 150, "+"),
    (3, 250, "-"),
    (2, 180, "+"),
    (4, 220, "-"),
    (3, 160, "+"),
    (2, 280, "-"),
]
GENE_GAP = 600              # intergenic gap between genes
ANCHORS = (25, 50, 75, 100, 147)
NOVEL_ACCEPTOR_SHIFT = 777  # push the acceptor off annotation -> novel pair


def _build_genes():
    """Return (genes, known_junctions) with deterministic coordinates.

    genes: list of (gene_id, transcript_id, strand, [(exon_start_1based, exon_end_1based)])
    known_junctions: list of (strand, donor_0based, acceptor_0based_exclusive)
    """
    genes = []
    known = []
    cursor = 1001  # 1-based genomic start of the first exon
    for i, (n_exons, intron_len, strand) in enumerate(GENE_SPECS):
        exons = []
        pos = cursor
        for _ in range(n_exons):
            exons.append((pos, pos + EXON_LEN - 1))
            pos += EXON_LEN + intron_len
        gid, tid = f"g{i}", f"t{i}"
        genes.append((gid, tid, strand, exons))
        # Consecutive exon boundaries -> junctions (0-based half-open).
        for (s0, e0), (s1, _e1) in zip(exons, exons[1:]):
            donor = e0            # 1-based inclusive end == 0-based exclusive end
            acceptor = s1 - 1     # 1-based start -> 0-based start
            known.append((strand, donor, acceptor))
        cursor = pos + GENE_GAP
    return genes, known, cursor


def _bed12_line(donor, acceptor, strand, anchor, name):
    chrom_start = donor - anchor
    chrom_end = acceptor + anchor
    return (
        f"{CHROM}\t{chrom_start}\t{chrom_end}\t{name}\t10\t{strand}\t"
        f"{chrom_start}\t{chrom_end}\t255,0,0\t2\t{anchor},{anchor}\t0,{acceptor - chrom_start}"
    )


def build(out_dir: Path = HERE):
    genes, known, end_coord = _build_genes()
    genome_len = end_coord + max(ANCHORS) + 500

    # genome.fa — content irrelevant to the flag; repeating ACGT for determinism.
    seq = ("ACGT" * (genome_len // 4 + 1))[:genome_len]
    fa = out_dir / "genome.fa"
    with fa.open("w") as fh:
        fh.write(f">{CHROM}\n")
        for i in range(0, genome_len, 60):
            fh.write(seq[i:i + 60] + "\n")

    # annotation.gtf — gene/transcript/exon records.
    gtf = out_dir / "annotation.gtf"
    with gtf.open("w") as fh:
        for gid, tid, strand, exons in genes:
            gstart, gend = exons[0][0], exons[-1][1]
            tattr = f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{gid}";'
            fh.write(f'{CHROM}\tsynthetic\tgene\t{gstart}\t{gend}\t.\t{strand}\t.\t'
                     f'gene_id "{gid}"; gene_name "{gid}";\n')
            fh.write(f'{CHROM}\tsynthetic\ttranscript\t{gstart}\t{gend}\t.\t{strand}\t.\t{tattr}\n')
            for s, e in exons:
                fh.write(f'{CHROM}\tsynthetic\texon\t{s}\t{e}\t.\t{strand}\t.\t{tattr}\n')

    # junctions.bed12 — known + fabricated-novel, varied anchors.
    bed12 = out_dir / "junctions.bed12"
    novel = [(st, d, a + NOVEL_ACCEPTOR_SHIFT) for (st, d, a) in known[::4]]
    with bed12.open("w") as fh:
        for i, (strand, donor, acceptor) in enumerate(known):
            fh.write(_bed12_line(donor, acceptor, strand,
                                 ANCHORS[i % len(ANCHORS)], f"KNOWN{i:03d}") + "\n")
        for i, (strand, donor, acceptor) in enumerate(novel):
            fh.write(_bed12_line(donor, acceptor, strand,
                                 ANCHORS[i % len(ANCHORS)], f"NOVEL{i:03d}") + "\n")

    print(f"Wrote {len(genes)} genes, {len(known)} known + {len(novel)} novel junctions")
    print(f"  {fa}\n  {gtf}\n  {bed12}")
    return fa, gtf, bed12, len(known), len(novel)


def _self_verify(fa, gtf, bed12):
    """Run the real crosscheck end-to-end if a regtools binary is reachable.

    Builds the reference BED from the same GTF (the home-rolled source of truth),
    then asserts the crosscheck passes — so a broken fixture is caught at
    regeneration time, not in CI.
    """
    regtools = shutil.which("regtools")
    if not regtools:
        print("regtools not on PATH — skipping self-verify (CI will run it).")
        return
    sys.path.insert(0, str(REPO_ROOT / "workflow" / "scripts"))
    from build_reference_junctions import extract_junctions, write_bed
    import crosscheck_annotate_flag as cc

    with tempfile.TemporaryDirectory() as td:
        ref_bed = Path(td) / "reference_junctions.bed"
        write_bed(extract_junctions(gtf), ref_bed)
        rc = cc.main([
            "--bed12", str(bed12), "--reference-fasta", str(fa),
            "--gtf", str(gtf), "--reference-bed", str(ref_bed),
            "--regtools-bin", regtools,
        ])
    if rc != 0:
        raise SystemExit("Self-verify FAILED — fixture does not cross-check clean.")
    print("Self-verify OK: fixture cross-checks at 100% agreement.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-verify", action="store_true",
                        help="skip the regtools self-verify step")
    args = parser.parse_args()
    fa, gtf, bed12, _, _ = build()
    if not args.no_verify:
        _self_verify(fa, gtf, bed12)
