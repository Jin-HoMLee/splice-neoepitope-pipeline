#!/usr/bin/env python3
"""crosscheck_annotate_flag.py — CI canary: home-rolled ``annotated`` flag vs.
``regtools junctions annotate``.

[Issue #370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370)
was a semantic-misread bug — BED12 coords that looked syntactically fine but
encoded the wrong thing (anchor outers vs. intron donor/acceptor). Unit tests of
the parser alone catch *that* bug but not *the next one of the same family*
(a swapped BED12 column, an off-by-one, a regtools version bump that changes
BED12 semantics). This canary adds a second, independent source of truth:
regtools' own annotator, which agrees with its BED12 semantics by construction.

Two flags are computed for the **same** input BED12 and compared:

1. **regtools**: ``regtools junctions annotate`` → ``known_junction`` column.
2. **home-rolled**: ``bed12_to_junctions`` → ``filter_junctions._parse_junction_id``
   → membership in the GENCODE-derived reference junction BED.

Both derive their *coordinates* by independent logic, so a coordinate-semantics
desync (the #370 bug class) makes the two flags disagree — and the canary fails
loudly instead of silently shipping wrong predictions.

Coordinate facts (verified empirically, regtools 1.0.0):
- annotate ``start`` = intron donor (0-based); annotate ``end`` =
  acceptor-exclusive **+ 1**. The home-rolled key (``chrom``, ``start``,
  ``end_exclusive``, ``strand``) is recovered as ``(chrom, start, end - 1, strand)``.
- ``known_junction`` is column 14 (1-indexed): ``1`` = annotated, ``0`` = novel.
- **A gzipped GTF makes regtools annotate silently emit zero rows** — the GTF is
  gunzipped to a temp file before the subprocess call.

See: https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/377
"""

import argparse
import contextlib
import gzip
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

# Importable both as a workflow script (conftest adds workflow/scripts to path
# for tests) and as a CLI run from the repo root in CI.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bed12_to_junctions import convert_bed12_to_junctions  # noqa: E402
from filter_junctions import _load_reference_junctions, _parse_junction_id  # noqa: E402

Key = tuple  # (chrom: str, start: int, end_exclusive: int, strand: str)

# Documented fixed column positions, used only if the regtools header is absent.
_FALLBACK_COLS = {"chrom": 0, "start": 1, "end": 2, "strand": 5, "known_junction": 13}


def parse_regtools_annotate(lines: Iterable[str]) -> dict:
    """Parse ``regtools junctions annotate`` output into ``{key: known_junction}``.

    ``key`` is ``(chrom, start, end - 1, strand)`` — the ``end - 1`` undoes
    regtools' acceptor-exclusive **+ 1** so the key matches the home-rolled one.
    Column positions are resolved from the header row by name (robust to a
    regtools version reordering columns); if no header is seen, documented fixed
    positions are used. Blank and ``#`` comment lines are skipped.
    """
    flags: dict = {}
    cols: Optional[dict] = None
    for raw in lines:
        line = raw.rstrip("\n")
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if cols is None and "known_junction" in fields and "chrom" in fields:
            cols = {name: i for i, name in enumerate(fields)}
            continue
        if cols is None:
            cols = dict(_FALLBACK_COLS)  # no header — parse this line as data
        try:
            chrom = fields[cols["chrom"]]
            start = int(fields[cols["start"]])
            end = int(fields[cols["end"]])
            strand = fields[cols["strand"]]
            known = int(fields[cols["known_junction"]])
        except (IndexError, ValueError):
            continue
        flags[(chrom, start, end - 1, strand)] = known
    return flags


def homerolled_flags(junction_ids: Iterable[str], reference: frozenset) -> dict:
    """Map each parseable ``chr:donor:acceptor:strand`` id → reference membership.

    ``1`` if the parsed junction is in the reference (annotated), else ``0``.
    Unparseable ids are skipped (not represented in the result).
    """
    flags: dict = {}
    for jid in junction_ids:
        parsed = _parse_junction_id(jid)
        if parsed is None:
            continue
        flags[parsed] = 1 if parsed in reference else 0
    return flags


@dataclass
class CrosscheckResult:
    n_regtools: int
    n_homerolled: int
    n_compared: int
    n_agree: int
    agreement: float
    coverage: float
    disagreements: list = field(default_factory=list)


def crosscheck(regtools: dict, homerolled: dict) -> CrosscheckResult:
    """Compare the two flag maps over their shared keys.

    ``agreement`` = fraction of shared keys whose flags match.
    ``coverage``  = shared keys / the larger of the two maps — guards against a
    total coordinate desync producing a tiny (vacuously 100%-agreeing)
    intersection. Both are ``0.0`` when there is nothing to compare, so an empty
    intersection reads as *failure*, never as healthy agreement.
    """
    common = regtools.keys() & homerolled.keys()
    n_compared = len(common)
    n_agree = 0
    disagreements: list = []
    for key in common:
        if regtools[key] == homerolled[key]:
            n_agree += 1
        else:
            disagreements.append((key, regtools[key], homerolled[key]))
    disagreements.sort()
    agreement = n_agree / n_compared if n_compared else 0.0
    denom = max(len(regtools), len(homerolled))
    coverage = n_compared / denom if denom else 0.0
    return CrosscheckResult(
        n_regtools=len(regtools),
        n_homerolled=len(homerolled),
        n_compared=n_compared,
        n_agree=n_agree,
        agreement=agreement,
        coverage=coverage,
        disagreements=disagreements,
    )


def run_regtools_annotate(
    bed12, reference_fasta, gtf, regtools_bin: str = "regtools"
) -> list:
    """Run ``regtools junctions annotate`` and return its stdout lines.

    A ``.gz`` GTF is gunzipped to a temp file first — regtools silently emits
    zero rows when handed a gzipped GTF.
    """
    with contextlib.ExitStack() as stack:
        gtf_path = str(gtf)
        if gtf_path.endswith(".gz"):
            tmp = stack.enter_context(
                tempfile.NamedTemporaryFile("wb", suffix=".gtf", delete=True)
            )
            with gzip.open(gtf, "rb") as src:
                shutil.copyfileobj(src, tmp)
            tmp.flush()
            gtf_path = tmp.name
        proc = subprocess.run(
            [regtools_bin, "junctions", "annotate", str(bed12),
             str(reference_fasta), gtf_path],
            capture_output=True, text=True, check=True,
        )
    return proc.stdout.splitlines()


def compute_homerolled(bed12, reference_bed) -> dict:
    """Run the production home-rolled path on ``bed12`` and return its flag map."""
    with tempfile.TemporaryDirectory() as td:
        tsv = Path(td) / "junctions.tsv"
        convert_bed12_to_junctions(bed12, tsv)
        ids = [
            line.split("\t", 1)[0]
            for line in tsv.read_text().splitlines()
            if line.strip()
        ]
    reference = _load_reference_junctions(reference_bed)
    return homerolled_flags(ids, reference)


def _print_report(result: CrosscheckResult, threshold: float) -> None:
    print("=== annotate-flag cross-check (Issue #377) ===")
    print(f"regtools junctions:    {result.n_regtools}")
    print(f"home-rolled junctions: {result.n_homerolled}")
    print(f"compared (shared keys): {result.n_compared}")
    print(f"agreement: {result.agreement:.4f}  (threshold {threshold})")
    print(f"coverage:  {result.coverage:.4f}  (threshold {threshold})")
    if result.disagreements:
        print(f"disagreements ({len(result.disagreements)}), first 20:")
        for key, rt, hr in result.disagreements[:20]:
            print(f"  {key[0]}:{key[1]}:{key[2]}:{key[3]}  regtools={rt} home-rolled={hr}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Cross-check the home-rolled annotated flag against "
        "regtools junctions annotate (CI canary, Issue #377)."
    )
    parser.add_argument("--bed12", required=True, help="regtools BED12 fixture")
    parser.add_argument("--reference-fasta", required=True, help="genome FASTA")
    parser.add_argument("--gtf", required=True, help="GENCODE GTF (.gtf or .gtf.gz)")
    parser.add_argument(
        "--reference-bed", required=True,
        help="GENCODE-derived reference junction BED (home-rolled source of truth)",
    )
    parser.add_argument("--regtools-bin", default="regtools", help="regtools binary")
    parser.add_argument(
        "--threshold", type=float, default=0.99,
        help="minimum agreement AND coverage fraction (default 0.99)",
    )
    args = parser.parse_args(argv)

    regtools = parse_regtools_annotate(
        run_regtools_annotate(
            args.bed12, args.reference_fasta, args.gtf, args.regtools_bin
        )
    )
    homerolled = compute_homerolled(args.bed12, args.reference_bed)
    result = crosscheck(regtools, homerolled)
    _print_report(result, args.threshold)

    ok = (
        result.n_compared > 0
        and result.agreement >= args.threshold
        and result.coverage >= args.threshold
    )
    if not ok:
        print(
            "FAIL: home-rolled annotated flag diverges from regtools "
            "junctions annotate — see disagreements above. This is the "
            "Issue #370 bug class (anchor outers vs. intron coords).",
            file=sys.stderr,
        )
        return 1
    print("OK: home-rolled annotated flag agrees with regtools.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
