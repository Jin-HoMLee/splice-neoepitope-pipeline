#!/usr/bin/env python3
"""Recompute Issue #1161 depth sensitivity from committed chr22 artifacts.

This is intentionally standard-library-only. It uses the committed raw junction
TSVs and the legacy mean-filter classifications from Issue #1162. Thresholds at
or above the legacy effective threshold are exact because every legacy-retained
row is available in the classified TSV. A threshold below that point is not
reported: its discarded rows lack the reference/GTEx classification provenance
needed for an exact candidate count.
"""
from __future__ import annotations

import csv
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "research" / "experiments" / "issue_1162_star_container"
OUT = Path(__file__).parent / "outputs" / "depth_sensitivity.tsv"

PANELS = {
    "hisat2_chr22": {
        "raw": BASE / "hisat2_off" / "SRR9143066_tumor_raw_junctions.tsv",
        "normal": BASE / "hisat2_off" / "SRR9143065_normal_raw_junctions.tsv",
        "classified": BASE / "hisat2_off" / "novel_junctions.tsv",
        "legacy_threshold": 2,
    },
    "star_chr22": {
        "raw": BASE / "star_off" / "SRR9143066_tumor_raw_junctions.tsv",
        "normal": BASE / "star_off" / "SRR9143065_normal_raw_junctions.tsv",
        "classified": BASE / "star_off" / "novel_junctions.tsv",
        "legacy_threshold": 3,
    },
}

ORIGINS = ("normal_shared", "gtex_pantissue_shared", "tumor_exclusive")


def read_raw(path: Path) -> dict[str, float]:
    rows: dict[str, float] = {}
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) >= 2:
            rows[fields[0]] = float(fields[1])
    return rows


def read_classified(path: Path) -> dict[str, tuple[float, str]]:
    rows: dict[str, tuple[float, str]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows[row["junction_id"]] = (float(row["mapped_reads"]), row["junction_origin"])
    return rows


def classify_at_threshold(
    classified: dict[str, tuple[float, str]],
    normal: dict[str, float],
    threshold: int,
) -> dict[str, int]:
    counts = {origin: 0 for origin in ORIGINS}
    for junction_id, (reads, legacy_origin) in classified.items():
        if reads < threshold:
            continue
        if legacy_origin == "normal_shared" and normal.get(junction_id, 0) < threshold:
            origin = "tumor_exclusive"
        elif legacy_origin == "normal_shared":
            origin = "normal_shared"
        elif legacy_origin == "gtex_pantissue_shared":
            origin = "gtex_pantissue_shared"
        elif legacy_origin == "tumor_exclusive":
            origin = "tumor_exclusive"
        else:
            origin = "annotated_discarded"
        counts[origin] += 1
    return counts


def emit_panel(panel: str, spec: dict[str, object], writer: csv.DictWriter) -> None:
    raw = read_raw(spec["raw"])  # type: ignore[arg-type]
    normal = read_raw(spec["normal"])  # type: ignore[arg-type]
    classified = read_classified(spec["classified"])  # type: ignore[arg-type]
    legacy_threshold = int(spec["legacy_threshold"])
    mean_reads = statistics.mean(raw.values())

    # The current output is the legacy mean-filter baseline, whose origin labels
    # used the historical two-read normal floor rather than a shared threshold.
    baseline = {origin: 0 for origin in ORIGINS}
    for _junction_id, (_reads, origin) in classified.items():
        if origin in baseline:
            baseline[origin] += 1
    legacy_below = sum(reads < legacy_threshold for reads in raw.values())
    legacy_retained = len(raw) - legacy_below
    writer.writerow({
        "panel": panel,
        "mode": "legacy_mean",
        "threshold": legacy_threshold,
        "raw_junctions": len(raw),
        "mean_reads": f"{mean_reads:.6f}",
        "effective_threshold": legacy_threshold,
        "below_min_reads": legacy_below,
        "retained_junctions": legacy_retained,
        "classified_rows": sum(baseline.values()),
        "annotated_discarded": legacy_retained - sum(baseline.values()),
        **baseline,
        "predicted_presenters": "not_available",
    })

    for threshold in (5, 10, 20):
        if threshold < legacy_threshold:
            continue
        counts = classify_at_threshold(classified, normal, threshold)
        below_min = sum(reads < threshold for reads in raw.values())
        retained = len(raw) - below_min
        writer.writerow({
            "panel": panel,
            "mode": "explicit_minimum",
            "threshold": threshold,
            "raw_junctions": len(raw),
            "mean_reads": f"{mean_reads:.6f}",
            "effective_threshold": threshold,
            "below_min_reads": below_min,
            "retained_junctions": retained,
            "classified_rows": sum(counts.values()),
            "annotated_discarded": retained - sum(counts.values()),
            **counts,
            "predicted_presenters": "not_available",
        })


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "panel", "mode", "threshold", "raw_junctions", "mean_reads",
        "effective_threshold", "below_min_reads", "retained_junctions",
        "classified_rows", "annotated_discarded", *ORIGINS, "predicted_presenters",
    ]
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for panel, spec in PANELS.items():
            emit_panel(panel, spec, writer)
    print(OUT)


if __name__ == "__main__":
    main()
