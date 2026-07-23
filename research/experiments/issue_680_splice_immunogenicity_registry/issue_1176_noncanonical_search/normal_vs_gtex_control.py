#!/usr/bin/env python3
"""Matched-pair control: what does a matched normal remove that GTEx cannot?

Issue #1176 needs a normal-filter decision for the Courcelles CRC cohort
(PXD071022 / GSE312236), where every one of the 26 RNA-seq runs is a tumor
biopsy. With no matched normal, `classify_junctions` labels every unannotated
tumor junction `tumor_exclusive` (with a warning) and the GTEx pan-tissue
population filter becomes the sole specificity control.

This script measures the cost of that substitution on data where we DO have a
matched normal (the chr22 gastric pair), by running the pipeline's own
`classify_junctions` twice with exactly one variable flipped:

  Run A: tumor + matched normal + GTEx   (our normal configuration)
  Run B: tumor + GTEx only               (the Courcelles configuration)

Everything else - tumor file, reference BED, GTEx blacklist, thresholds - is
byte-identical between runs, so the delta is attributable to the normal arm
alone.

FALSIFIER (this check can fail, and would have): if run B equals run A, the
matched normal contributes nothing beyond the population filter. A non-zero
delta is the number of patient-private non-tumor junctions that a population
reference structurally cannot see, and which therefore become false
`tumor_exclusive` candidates.

Prerequisite: the GTEx chr22 blacklist is gitignored and regenerable. Build it
with (note --region, NOT --restrict-chrom alone - see README):

    python workflow/scripts/build_gtex_pan_tissue_ref.py \\
      --output-bed resources/test/gtex_gtexv2_pan_tissue_junctions.chr22.bed \\
      --output-qc  resources/test/gtex_gtexv2_pan_tissue_junctions.chr22.qc.json \\
      --min-samples 1 --region chr22 --restrict-chrom chr22

Run from the repo root:  python <this script>
"""
import csv
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(
    subprocess.run(["git", "rev-parse", "--show-toplevel"],
                   capture_output=True, text=True, check=True).stdout.strip()
)
HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
sys.path.insert(0, str(REPO / "workflow" / "scripts"))

from filter_junctions import classify_junctions  # noqa: E402

TUMOR = REPO / "results/patient_001_test/alignment/SRR9143066_test/junctions.tsv"
NORMAL = REPO / "results/patient_001_test/alignment/SRR9143065_test/junctions.tsv"
REF_BED = REPO / "resources/test/chr22_reference_junctions.bed"
GTEX_BED = REPO / "resources/test/gtex_gtexv2_pan_tissue_junctions.chr22.bed"
GTF = REPO / "resources/test/chr22.gtf.gz"

FUNNEL = ["junctions_raw", "mean_reads_filtered", "annotated_discarded",
          "normal_shared", "gtex_pantissue_shared", "tumor_exclusive"]


def _write_manifest(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["file_id", "sample_type"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def _run(tag, junction_files, manifest_rows):
    man, out, stats = OUT / f"manifest_{tag}.tsv", OUT / f"classified_{tag}.tsv", OUT / f"stats_{tag}.tsv"
    _write_manifest(man, manifest_rows)
    classify_junctions(
        junction_files=[str(p) for p in junction_files],
        manifest_path=str(man),
        reference_bed=str(REF_BED),
        output_path=str(out),
        stats_output_path=str(stats),
        min_normal_reads=2,
        gencode_gtf=str(GTF),
        gtex_bed=str(GTEX_BED),
    )
    counts = {}
    with open(stats) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            # stats values are emitted as floats; the funnel is integral
            counts[row["category"]] = int(float(row["count"]))
    with open(out) as fh:
        rows = {r["junction_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    return counts, rows


def main():
    missing = [p for p in (TUMOR, NORMAL, REF_BED, GTEX_BED, GTF) if not p.exists()]
    if missing:
        print("FATAL: missing inputs:", *(f"\n  {p}" for p in missing))
        print("\nThe GTEx blacklist is gitignored/regenerable - see the module docstring.")
        sys.exit(1)
    OUT.mkdir(exist_ok=True)

    a, rows_a = _run("A_with_normal", [TUMOR, NORMAL],
                     [{"file_id": "SRR9143066_test", "sample_type": "Primary Tumor"},
                      {"file_id": "SRR9143065_test", "sample_type": "Solid Tissue Normal"}])
    b, rows_b = _run("B_no_normal", [TUMOR],
                     [{"file_id": "SRR9143066_test", "sample_type": "Primary Tumor"}])

    sel = lambda rows, origin: {k for k, v in rows.items() if v["junction_origin"] == origin}
    ns_a = sel(rows_a, "normal_shared")
    te_a, te_b = sel(rows_a, "tumor_exclusive"), sel(rows_b, "tumor_exclusive")
    gx_b = sel(rows_b, "gtex_pantissue_shared")
    recovered, lost = sorted(ns_a & gx_b), sorted(ns_a & te_b)

    print(f"{'category':<24}{'A (matched normal)':>20}{'B (GTEx only)':>16}")
    for k in FUNNEL:
        print(f"  {k:<22}{a.get(k, '-'):>20}{b.get(k, '-'):>16}")

    # The funnel must reconcile in both arms, or the comparison is meaningless.
    for tag, c in (("A", a), ("B", b)):
        total = sum(c[k] for k in FUNNEL[1:])
        assert total == c["junctions_raw"], f"{tag} funnel does not reconcile: {total} != {c['junctions_raw']}"
    print("\nboth funnels reconcile against junctions_raw")

    delta = len(te_b) - len(te_a)
    print(f"\nnormal_shared with a matched normal : {len(ns_a)}")
    print(f"  also blacklisted by GTEx          : {len(recovered)}")
    print(f"  LOST (population filter blind)    : {len(lost)}")
    print(f"tumor_exclusive  A={len(te_a)}  B={len(te_b)}  delta=+{delta} "
          f"({delta / len(te_a) * 100:.1f}% inflation)")

    if delta == 0:
        print("\nFALSIFIER FIRED: the matched normal adds nothing beyond GTEx.")
    else:
        print("\nCaught ONLY by the matched normal (patient-private, population-invisible):")
        for j in lost:
            print(f"  {j}  reads={rows_b[j]['mapped_reads']}")

    json.dump({"A": a, "B": b, "delta_tumor_exclusive": delta,
               "normal_shared_a": len(ns_a), "recovered_by_gtex": recovered,
               "lost_to_tumor_exclusive": lost},
              open(OUT / "normal_vs_gtex.json", "w"), indent=2, sort_keys=True)
    print(f"\nwrote {OUT / 'normal_vs_gtex.json'}")


if __name__ == "__main__":
    main()
