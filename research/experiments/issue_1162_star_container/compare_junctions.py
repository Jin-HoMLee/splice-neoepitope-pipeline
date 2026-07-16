#!/usr/bin/env python3
"""Compare junction outputs for Issue #1162 AC4.

Two comparisons:

1. aligner (STAR vs HISAT2), knob held fixed - characterizes the residual
   difference and lets us attribute it to the aligner algorithm rather than the
   read-support semantic.
2. knob (uniqueness off vs on), aligner held fixed - the matched-pair control.
   Prediction (AC4): with the knob OFF, junctions carrying only multi-mapped
   support must be present and per-junction counts must include multi-mappers;
   turning the knob ON must drop exactly those. If the two runs are identical,
   the knob is not wired into that path - and that is the finding.

Inputs are the classified `novel_junctions.tsv` (10-col) and/or the raw
`raw_junctions.tsv` (2-col: junction_id<TAB>reads). Set membership is keyed on
`junction_id`; for novel_junctions we restrict to `junction_origin ==
tumor_exclusive` (the set that reaches neoepitope prediction) unless --all-origins.
"""
import argparse
import sys
from pathlib import Path


def load_novel(path: Path, all_origins: bool) -> dict[str, int]:
    """junction_id -> mapped_reads, restricted to tumor_exclusive by default."""
    out: dict[str, int] = {}
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            if not all_origins and f[idx["junction_origin"]] != "tumor_exclusive":
                continue
            out[f[idx["junction_id"]]] = int(float(f[idx["mapped_reads"]]))
    return out


def load_raw(path: Path) -> dict[str, int]:
    """junction_id -> reads from the 2-col raw_junctions.tsv."""
    out: dict[str, int] = {}
    with path.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            out[f[0]] = int(float(f[1]))
    return out


def report(label_a: str, a: dict[str, int], label_b: str, b: dict[str, int]) -> None:
    sa, sb = set(a), set(b)
    both = sa & sb
    only_a = sa - sb
    only_b = sb - sa
    print(f"\n=== {label_a}  vs  {label_b} ===")
    print(f"  {label_a:<24} : {len(sa)} junctions")
    print(f"  {label_b:<24} : {len(sb)} junctions")
    print(f"  shared                   : {len(both)}")
    print(f"  {label_a}-only           : {len(only_a)}")
    print(f"  {label_b}-only           : {len(only_b)}")
    if sa == sb:
        # count deltas on the shared set (matters for the knob matched-pair)
        diffs = [(j, a[j], b[j]) for j in both if a[j] != b[j]]
        print(f"  identical membership; per-junction count diffs: {len(diffs)}")
        for j, va, vb in sorted(diffs)[:10]:
            print(f"    {j}: {label_a}={va} {label_b}={vb}")
    else:
        for j in sorted(only_a)[:10]:
            print(f"    only {label_a}: {j} (reads={a[j]})")
        for j in sorted(only_b)[:10]:
            print(f"    only {label_b}: {j} (reads={b[j]})")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--a", required=True, help="TSV (novel_junctions or raw_junctions)")
    p.add_argument("--b", required=True, help="TSV to compare against")
    p.add_argument("--label-a", default="A")
    p.add_argument("--label-b", default="B")
    p.add_argument("--raw", action="store_true", help="inputs are 2-col raw_junctions.tsv")
    p.add_argument("--all-origins", action="store_true",
                   help="novel_junctions: include all origins, not just tumor_exclusive")
    args = p.parse_args()

    loader = load_raw if args.raw else (lambda path: load_novel(path, args.all_origins))
    a = loader(Path(args.a))
    b = loader(Path(args.b))
    report(args.label_a, a, args.label_b, b)
    return 0


if __name__ == "__main__":
    sys.exit(main())
