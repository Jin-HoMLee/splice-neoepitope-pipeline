#!/usr/bin/env python3
"""alphagenome_spike.py — Live API smoke test for Issue #223.

Verifies that the AlphaGenome API:
  1. Authenticates with the key in $ALPHAGENOME_API_KEY (or .env).
  2. Returns OutputType.SPLICE_JUNCTIONS for a small chr22 reference interval.
  3. Returns variant-vs-reference junction outputs for a fake variant context.
  4. Surfaces wall-clock latency for both calls.

This is a one-shot probe — does not write any artefacts the pipeline depends on.
Output goes to stdout so it can be pasted into the Issue #223 thread.
"""

import os
import sys
import time


def load_env(path=".env"):
    """Mirror research/scripts/zotero_add.py — minimal .env loader, stdlib only."""
    if not os.path.exists(path):
        return
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                key, _, value = line.partition("=")
                os.environ.setdefault(key.strip(), value.strip())


def summarise_output(output, label):
    """Print a compact summary of an Output object's splice junction track."""
    print(f"\n  {label}:")
    print(f"    output type: {type(output).__name__}")
    sj = getattr(output, "splice_junctions", None)
    if sj is None:
        print("    splice_junctions: <attribute missing>")
        return
    print(f"    splice_junctions type: {type(sj).__name__}")
    # JunctionData exposes .values (np.ndarray), .metadata (pd.DataFrame),
    # .interval (Interval). For values + interval the shape is the useful
    # summary; for metadata we surface column names directly because
    # DataFrame.shape is also populated (e.g. (367, 8)) but less informative
    # than the column inventory.
    for attr in ("values", "interval"):
        val = getattr(sj, attr, None)
        if val is None:
            continue
        shape = getattr(val, "shape", None)
        if shape is not None:
            print(f"    {attr}.shape: {shape}")
        else:
            print(f"    {attr}: {val!r}")
    md = getattr(sj, "metadata", None)
    if md is not None:
        cols = list(md.columns)
        print(f"    metadata: DataFrame with {len(md)} rows, "
              f"columns={cols[:8]}{'…' if len(cols) > 8 else ''}")


def main() -> int:
    load_env()
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        sys.exit("Error: ALPHAGENOME_API_KEY must be set in .env or environment.")

    # Late imports so a missing-SDK error doesn't obscure the env-var check above.
    from alphagenome.data import genome
    from alphagenome.models import dna_client
    from alphagenome.models.dna_output import OutputType

    # AlphaGenome only accepts these four specific input lengths:
    # 16384, 131072, 524288, 1048576. Other powers of 2 (e.g. 32768, 65536,
    # 262144) raise ValueError. We use 131072 (~128 kb) — smallest size that
    # comfortably spans a typical multi-exon gene.
    INTERVAL_WIDTH = 131_072
    INTERVAL_START = 42_000_000
    interval = genome.Interval(
        chromosome="chr22",
        start=INTERVAL_START,
        end=INTERVAL_START + INTERVAL_WIDTH,
    )
    # Synthetic SNV inside the interval — does not need to be real for the probe;
    # purpose is to verify the variant API path works and returns ref/alt outputs.
    variant = genome.Variant(
        chromosome="chr22",
        position=INTERVAL_START + INTERVAL_WIDTH // 2,
        reference_bases="A",
        alternate_bases="G",
    )

    print(f"Creating client …")
    model = dna_client.create(api_key)
    print(f"  {type(model).__name__} ready.")
    print(f"  Using SPLICE_JUNCTIONS output type "
          f"({OutputType.SPLICE_JUNCTIONS.name}).")

    print(f"\nCall 1 — predict_interval on {interval}")
    t0 = time.perf_counter()
    interval_output = model.predict_interval(
        interval=interval,
        requested_outputs=[OutputType.SPLICE_JUNCTIONS],
        ontology_terms=None,
    )
    elapsed1 = time.perf_counter() - t0
    print(f"  wall-clock: {elapsed1:.2f}s")
    summarise_output(interval_output, "interval-only result")

    print(f"\nCall 2 — predict_variant on {interval} with {variant}")
    t0 = time.perf_counter()
    variant_output = model.predict_variant(
        interval=interval,
        variant=variant,
        requested_outputs=[OutputType.SPLICE_JUNCTIONS],
        ontology_terms=None,
    )
    elapsed2 = time.perf_counter() - t0
    print(f"  wall-clock: {elapsed2:.2f}s")
    # VariantOutput exposes .reference and .alternate Output objects;
    # show both so we can confirm the diff is observable.
    if hasattr(variant_output, "reference"):
        summarise_output(variant_output.reference, "variant.reference")
    if hasattr(variant_output, "alternate"):
        summarise_output(variant_output.alternate, "variant.alternate")
    if not hasattr(variant_output, "reference"):
        # Fallback if the API shape differs — dump non-private attrs.
        print(f"  VariantOutput attrs: "
              f"{[a for a in dir(variant_output) if not a.startswith('_')]}")

    # Inspect track metadata so Scientist can plan tissue filters for #224/#225.
    sj = interval_output.splice_junctions
    md = getattr(sj, "metadata", None)
    if md is not None and hasattr(md, "columns"):
        print(f"\nTrack metadata: {len(md)} tracks, columns={list(md.columns)}")
        if "data_source" in md.columns:
            print(f"  data_source counts: "
                  f"{md['data_source'].value_counts().to_dict()}")
        if "Assay title" in md.columns:
            print(f"  Assay title counts (top 5): "
                  f"{md['Assay title'].value_counts().head(5).to_dict()}")

    print(f"\nDone. Latency: predict_interval={elapsed1:.2f}s, "
          f"predict_variant={elapsed2:.2f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
