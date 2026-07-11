#!/usr/bin/env python3
"""junction_repeat_overlap.py - Compare two junction sets and categorize the difference by repeat overlap.

Built to evaluate the NH-uniqueness prefilter (Issue #919), and reusable for the
whole-genome re-run that Issue #1095 exists for.

The filter's claim is specific: HISAT2 emits one arbitrary copy of a multimapped
read, so spliced reads over repeats produce junction calls at arbitrary repeat
copies, and dropping ``NH>1`` alignments should remove those. That claim is
falsifiable. If the junctions the filter removes are *not* enriched for repeat
overlap relative to the ones it keeps, the filter is removing something else -
low-coverage junctions, say - and the mechanism story is wrong even if the
counts look plausible.

So this compares a junction set produced with the filter off (A) against one
produced with it on (B), and reports repeat overlap for the junctions **lost**
(in A, not B) against the junctions **retained** (in both), as an enrichment
ratio. Junctions **gained** (in B, not A) are reported too: the filter can only
remove alignments, so a gain is a red flag, not a feature.

A junction is called repeat-overlapping when its **donor or acceptor splice
site** falls inside a RepeatMasker interval. The splice sites are the right
probe, not the whole intron: an intron spanning a repeat in its middle is
unremarkable (introns are repeat-dense), whereas a splice site *inside* a repeat
copy is the signature of the misplaced-read failure mode.

Junction IDs are the ``chrom:donor:acceptor:strand`` form emitted by
``bed12_to_junctions.py``, where ``donor`` is 1-based and ``acceptor`` is
0-based-exclusive, so the intron is ``[donor - 1, acceptor)`` in 0-based
half-open coordinates.

Usage:
    python workflow/scripts/junction_repeat_overlap.py \\
        --junctions-off results/.../raw_junctions.tsv \\
        --junctions-on  results/.../raw_junctions.filtered.tsv \\
        --rmsk references/rmsk/hg38/rmsk.chr22.bed \\
        --label "tumor (SRR9143066)"
"""

import argparse
import bisect
import logging
from collections import Counter
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


class Repeat(NamedTuple):
    start: int          # 0-based inclusive
    end: int            # 0-based exclusive
    name: str
    rep_class: str
    rep_family: str


class Junction(NamedTuple):
    chrom: str
    start: int          # 0-based inclusive intron start (the donor site)
    end: int            # 0-based exclusive intron end
    strand: str
    reads: int

    @property
    def junction_id(self) -> str:
        return f"{self.chrom}:{self.start + 1}:{self.end}:{self.strand}"

    @property
    def donor_site(self) -> int:
        """0-based position of the first intronic base at the donor."""
        return self.start

    @property
    def acceptor_site(self) -> int:
        """0-based position of the last intronic base at the acceptor."""
        return self.end - 1


class RepeatIndex:
    """Per-chromosome interval index over a RepeatMasker BED.

    Repeats can nest and overlap, so a point lookup cannot stop at the first hit
    found by bisect; it scans back over any interval that could still span the
    query. The scan-back window is the longest repeat on that chromosome, which
    keeps the lookup bounded without assuming the intervals are disjoint.
    """

    def __init__(self, by_chrom: Dict[str, List[Repeat]]):
        self._by_chrom = {}
        for chrom, repeats in by_chrom.items():
            ordered = sorted(repeats, key=lambda r: r.start)
            starts = [r.start for r in ordered]
            longest = max((r.end - r.start for r in ordered), default=0)
            self._by_chrom[chrom] = (ordered, starts, longest)

    @property
    def chromosomes(self) -> List[str]:
        return sorted(self._by_chrom)

    def __len__(self) -> int:
        return sum(len(v[0]) for v in self._by_chrom.values())

    def at(self, chrom: str, position: int) -> List[Repeat]:
        """Every repeat spanning `position` (0-based) on `chrom`."""
        entry = self._by_chrom.get(chrom)
        if entry is None:
            return []
        ordered, starts, longest = entry
        hits = []
        i = bisect.bisect_right(starts, position)
        j = i - 1
        while j >= 0 and starts[j] >= position - longest:
            repeat = ordered[j]
            if repeat.start <= position < repeat.end:
                hits.append(repeat)
            j -= 1
        return hits


def load_rmsk(path: Path) -> RepeatIndex:
    """Load a BED6+2 RepeatMasker file as written by fetch_rmsk.py."""
    by_chrom: Dict[str, List[Repeat]] = {}
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(
                    f"{path}: expected BED6+2 (chrom start end name score strand "
                    f"repClass repFamily), got {len(fields)} columns"
                )
            by_chrom.setdefault(fields[0], []).append(
                Repeat(
                    start=int(fields[1]),
                    end=int(fields[2]),
                    name=fields[3],
                    rep_class=fields[6],
                    rep_family=fields[7],
                )
            )
    return RepeatIndex(by_chrom)


def parse_junction_id(junction_id: str, reads: int = 0) -> Junction:
    """Parse `chrom:donor:acceptor:strand` into 0-based half-open coordinates.

    Mirrors filter_junctions._parse_junction_id: donor is 1-based (so the intron
    start is donor - 1), acceptor is already 0-based exclusive.
    """
    parts = junction_id.split(":")
    if len(parts) != 4:
        raise ValueError(f"malformed junction id: {junction_id!r}")
    chrom, donor, acceptor, strand = parts
    return Junction(
        chrom=chrom,
        start=int(donor) - 1,
        end=int(acceptor),
        strand=strand,
        reads=reads,
    )


def load_junctions(path: Path) -> Dict[str, Junction]:
    """Load a raw_junctions.tsv (`junction_id<TAB>reads`), keyed by junction id."""
    junctions: Dict[str, Junction] = {}
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            reads = int(fields[1]) if len(fields) > 1 else 0
            junction = parse_junction_id(fields[0], reads)
            junctions[junction.junction_id] = junction
    return junctions


def classify(junction: Junction, index: RepeatIndex) -> Tuple[List[Repeat], List[Repeat]]:
    """Repeats spanning the donor site and the acceptor site, respectively."""
    return index.at(junction.chrom, junction.donor_site), index.at(
        junction.chrom, junction.acceptor_site
    )


def overlaps_repeat(junction: Junction, index: RepeatIndex) -> bool:
    """True when either splice site falls inside a repeat."""
    donor, acceptor = classify(junction, index)
    return bool(donor or acceptor)


def repeat_classes(junction: Junction, index: RepeatIndex) -> List[str]:
    """Distinct repeat classes touching either splice site."""
    donor, acceptor = classify(junction, index)
    return sorted({r.rep_class for r in donor + acceptor})


class Comparison(NamedTuple):
    lost: List[Junction]        # in OFF, not in ON - what the filter removed
    retained: List[Junction]    # in both
    gained: List[Junction]      # in ON, not in OFF - should be empty


def compare(off: Dict[str, Junction], on: Dict[str, Junction]) -> Comparison:
    lost = [j for jid, j in off.items() if jid not in on]
    retained = [j for jid, j in off.items() if jid in on]
    gained = [j for jid, j in on.items() if jid not in off]
    return Comparison(lost=lost, retained=retained, gained=gained)


def repeat_rate(junctions: Sequence[Junction], index: RepeatIndex) -> Tuple[int, int, float]:
    """(n_overlapping, n_total, fraction) for a junction set."""
    total = len(junctions)
    hits = sum(1 for j in junctions if overlaps_repeat(j, index))
    return hits, total, (hits / total if total else 0.0)


def class_counter(junctions: Sequence[Junction], index: RepeatIndex) -> Counter:
    counts: Counter = Counter()
    for junction in junctions:
        classes = repeat_classes(junction, index)
        if not classes:
            counts["(none)"] += 1
        for rep_class in classes:
            counts[rep_class] += 1
    return counts


def _pct(hits: int, total: int) -> str:
    return f"{hits}/{total} ({100 * hits / total:.1f}%)" if total else f"0/0 (n/a)"


def format_report(comparison: Comparison, index: RepeatIndex, label: str) -> str:
    """Markdown report: does the filter remove repeat-overlapping junctions?"""
    lost_hits, lost_n, lost_frac = repeat_rate(comparison.lost, index)
    kept_hits, kept_n, kept_frac = repeat_rate(comparison.retained, index)
    enrichment = (lost_frac / kept_frac) if kept_frac else float("inf")

    lines = [
        f"## Junction/repeat overlap - {label}",
        "",
        f"| set | junctions | splice site in a repeat |",
        f"|---|---|---|",
        f"| lost to the filter | {lost_n} | {_pct(lost_hits, lost_n)} |",
        f"| retained | {kept_n} | {_pct(kept_hits, kept_n)} |",
        f"| gained (should be 0) | {len(comparison.gained)} | - |",
        "",
        f"**Enrichment of repeat overlap among lost junctions: {enrichment:.2f}x**",
        "",
    ]

    single_read_lost = sum(1 for j in comparison.lost if j.reads == 1)
    if comparison.lost:
        lines.append(
            f"Single-read junctions among the lost: "
            f"{_pct(single_read_lost, len(comparison.lost))}"
        )
        lines.append("")

    lost_classes = class_counter(comparison.lost, index)
    kept_classes = class_counter(comparison.retained, index)
    lines.extend([
        "| repeat class at splice site | lost | retained |",
        "|---|---|---|",
    ])
    for rep_class, count in lost_classes.most_common():
        lines.append(f"| {rep_class} | {count} | {kept_classes.get(rep_class, 0)} |")
    lines.append("")
    return "\n".join(lines)


def write_tsv(comparison: Comparison, index: RepeatIndex, path: Path) -> None:
    """Per-junction categorization, so the report's numbers can be audited."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("junction_id\tfate\treads\tdonor_in_repeat\tacceptor_in_repeat\trepeat_classes\n")
        for fate, junctions in (
            ("lost", comparison.lost),
            ("retained", comparison.retained),
            ("gained", comparison.gained),
        ):
            for junction in sorted(junctions, key=lambda j: (j.chrom, j.start)):
                donor, acceptor = classify(junction, index)
                handle.write(
                    f"{junction.junction_id}\t{fate}\t{junction.reads}\t"
                    f"{bool(donor)}\t{bool(acceptor)}\t"
                    f"{','.join(repeat_classes(junction, index)) or '.'}\n"
                )


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--junctions-off", required=True, type=Path, help="raw_junctions.tsv, filter OFF")
    parser.add_argument("--junctions-on", required=True, type=Path, help="raw_junctions.tsv, filter ON")
    parser.add_argument("--rmsk", required=True, type=Path, help="RepeatMasker BED from fetch_rmsk.py")
    parser.add_argument("--label", default="chr22", help="Label for the report heading")
    parser.add_argument("--output-tsv", type=Path, help="Optional per-junction categorization TSV")
    args = parser.parse_args(argv)

    index = load_rmsk(args.rmsk)
    log.info("Loaded %d repeats across %s", len(index), ", ".join(index.chromosomes))

    off = load_junctions(args.junctions_off)
    on = load_junctions(args.junctions_on)
    log.info("Junctions: %d (filter off) vs %d (filter on)", len(off), len(on))

    comparison = compare(off, on)
    if comparison.gained:
        # The filter only ever removes alignments, so a gained junction means the
        # comparison is not measuring what it claims to.
        log.warning(
            "%d junctions appear ONLY with the filter on - the filter cannot add "
            "alignments, so this needs explaining before the report is trusted",
            len(comparison.gained),
        )

    print(format_report(comparison, index, args.label))

    if args.output_tsv:
        write_tsv(comparison, index, args.output_tsv)
        log.info("Wrote per-junction categorization to %s", args.output_tsv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
