# AlphaGenome Primer — Reading the Output

Reference for interpreting AlphaGenome predictions in the context of this
pipeline (predicted-normal junction filter, [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203)
and follow-ups). Written after the [Issue #223](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/223)
spike — working code lives in [`scripts/alphagenome_spike.py`](../scripts/alphagenome_spike.py).

The numeric examples below come from one live API call against the spike
script's chr22 test region. Numbers are illustrative; rerun the spike for
fresh values.

---

## Table of Contents

1. [What AlphaGenome is](#what-alphagenome-is)
2. [Tracks — the core concept](#tracks--the-core-concept)
3. [Output shape](#output-shape)
4. [API methods — query shapes, not different assays](#api-methods--query-shapes-not-different-assays)
5. [`data_source` is training-time provenance, not runtime](#data_source-is-training-time-provenance-not-runtime)
6. [Per-track interpretation](#per-track-interpretation)
7. [Operational details](#operational-details)
8. [For the predicted-normal filter](#for-the-predicted-normal-filter)

---

## What AlphaGenome is

AlphaGenome is DeepMind's sequence-to-regulation foundation model: it takes
a DNA sequence as input and emits predicted functional-genomics measurements
(RNA-seq, ATAC-seq, ChIP-seq, splice junctions, contact maps, etc.) in a
single forward pass. **Free for non-commercial use** via Python SDK and
hosted API; key obtained at <https://deepmind.google.com/science/alphagenome>.

End-to-end means: there is no internal step where AlphaGenome simulates raw
RNA-seq reads and re-analyses them. The model maps DNA → predicted assay
summary statistic directly. **No internal aligner, no FASTQ generation, no
runtime database lookups** — everything the model "knows" is baked into
its trained weights.

---

## Tracks — the core concept

A **track** is one output channel of the model. Each track corresponds to
one specific (biosample × assay × consortium) combination the model was
calibrated against during training.

For example, one row of the metadata DataFrame:

| metadata field | example value |
|---|---|
| `name` | `junction_CL:0000084 polyA plus RNA-seq` |
| `ontology_curie` | `CL:0000084` (Cell Ontology — T-cell) |
| `biosample_name` | `T-cell` |
| `biosample_type` | `primary_cell` |
| `biosample_life_stage` | `adult` |
| `gtex_tissue` | (empty for ENCODE-derived tracks) |
| `data_source` | `encode` |
| `Assay title` | `polyA plus RNA-seq` |

When we ask for `OutputType.SPLICE_JUNCTIONS`, we get **367 tracks** in
total: 313 from ENCODE, 54 from GTEx; 195 trained on total RNA-seq, 172 on
polyA-plus RNA-seq (counts from the spike output, current as of 2026-05-02).

**Tracks are NOT sub-models.** AlphaGenome is one model with a shared
backbone (the DNA encoder) and 367 lightweight output heads. All heads were
trained jointly via multi-task learning. What the model learned about gene
structure from ENCODE training data informs its GTEx predictions too —
they share the encoder. That's why one inference call returns predictions
for all 367 tracks at once in ~1–3 seconds.

### Glossary terms used here

- **Biosample** — consortium jargon for "the biological material the
  experiment was performed on." Covers primary cells, immortalised cell
  lines, in-vitro-differentiated cells, and bulk tissue samples. Broader
  than "tissue."
- **Assay** — what experimental measurement the track predicts: RNA-seq,
  ATAC-seq, ChIP-seq, splice junctions, etc. Set by `OutputType` in the SDK.
- **Consortium** — which lab/project produced the original measurements:
  ENCODE (cell lines + primary cells + tissues), GTEx (postmortem human
  tissues — our preferred source for "normal").
- **RNA-seq protocols in the inventory:**
  - *polyA plus RNA-seq* — captures only mature, polyadenylated mRNA
  - *total RNA-seq* — captures all RNA species (mRNA + ncRNA + pre-mRNA)
  Different protocols sample different RNA populations → different
  junction signal.

---

## Output shape

For `requested_outputs=[OutputType.SPLICE_JUNCTIONS]`, the SDK returns an
`Output` (or `VariantOutput`) whose `.splice_junctions` is a `JunctionData`
object with:

- `.values` — the predicted-signal matrix, shape `(n_junctions, n_tracks)`
- `.metadata` — DataFrame describing the tracks (one row per column of
  `.values`), shape `(n_tracks, 8)`
- `.interval` — the genomic region the prediction is for

Concrete shape for the spike's chr22:42M-42.13M call to `predict_interval`:

```
.values: (6084, 367)        — 6084 candidate junctions × 367 tracks
.metadata: (367, 8)         — 367 track descriptors

                    track 0     track 1    ...   track 366
                    (T-cell     (liver           (lung
                    polyA)      total)           polyA)
junction 0          0.04        0.92      ...    0.01
junction 1          0.81        0.03      ...    0.78
junction 2          0.00        0.51      ...    0.00
...
junction 6083       0.65        0.05      ...    0.41
```

**Critical: junctions are SHARED across tracks.** The model emits one set
of candidate junctions (the rows) and provides per-track predicted signal
for each (the columns). A junction with high signal in track A (T-cell)
and low signal in track B (liver) is the model expressing tissue
specificity. **Tissue specificity is a property of the score distribution
within a column, not of the junction set.**

---

## API methods — query shapes, not different assays

The SDK exposes several methods that look like they predict different
things. They don't — what is predicted is controlled by `requested_outputs`,
which is orthogonal to the API method. **You can request any `OutputType`
through any API method.**

| Method | Input | Returns |
|---|---|---|
| `predict_interval(interval)` | a genome `Interval` | `Output` — fetches reference, predicts on it |
| `predict_sequence(sequence, interval=...)` | a raw DNA **string** + optional interval | `Output` — predicts on the string you provide |
| `predict_variant(interval, variant)` | interval + one specified variant | `VariantOutput` with `.reference` and `.alternate` |
| `predict_intervals(intervals)` | sequence of intervals | list of `Output`s (parallel batch) |
| `predict_variants(intervals, variants)` | sequences of intervals + variants | list of `VariantOutput`s (parallel batch) |

`predict_variant` does NOT generate variants — you specify exactly one
variant as input, and the function returns predictions for both the
reference and the variant-applied DNA. Same for `predict_variants`: it
processes N user-specified variants independently in parallel, returning
N separate `VariantOutput`s — it is **NOT** a multi-variant compositor.

### Empirical: `predict_variant` filters to variant-affected junctions

In the spike, `predict_interval` on a 131kb chr22 region returned **6084
junctions**. `predict_variant` on the same interval (with one synthetic
SNV) returned only **503 junctions** in both `.reference` and `.alternate`.

We verified empirically that all 503 junctions show ref ≠ alt in at least
one track (none have identical values across all 367 tracks). So
`predict_variant` is filtering its output to **only the junctions where
the variant has measurable effect somewhere in the (junctions × tracks)
matrix** — not geographic windowing, not random sub-sampling.

The inclusion threshold appears very lenient: in the spike, mean |ref − alt|
across the 503 included junctions is ~0.00014. Effectively any non-zero
variant effect, however small, qualifies a junction for output. This means
the 503 number is **"all junctions with any model-detectable response to
the variant"** — not a high-confidence "this variant matters here" set.
The meaningful biological signal lives in the *magnitude* of ref-vs-alt
deltas, not the count of returned junctions.

---

## `data_source` is training-time provenance, not runtime

The `data_source` column in track metadata (`encode` / `gtex`) describes
**which consortium produced the training data** that the corresponding
output head was calibrated against — not a runtime input.

At inference, AlphaGenome does NOT consult ENCODE or GTEx databases. The
model performs a pure DNA-sequence forward pass; the consortium-trained
patterns live in the model weights. No internet access to consortium
archives is needed beyond the API call itself.

Why this matters for our use case: GTEx profiled splicing in **healthy
human tissues** from postmortem donors, which is exactly the "normal"
baseline the predicted-normal filter wants. ENCODE's track inventory
includes many cancer-derived cell lines (HeLa, K562) — useful for general
regulatory inference but a poor proxy for "normal patient tissue."
Filtering to `data_source == 'gtex'` gives us 54 tracks aligned with our
question.

---

## Per-track interpretation

For one track at one junction, the `.values[junction_i, track_j]` cell
gives the model's predicted signal under the following counterfactual:

> *If* the input DNA were transcribed in cells of biosample X, *and* you
> ran assay Y on the resulting RNA, *then* the measured signal at junction
> J would be approximately this value.

The reference DNA is fixed across tracks — what varies is which cellular
context (which output head) you're reading from. Tissue specificity is
encoded in the model's weights (different output heads for different
biosamples), not in the input.

### Worked example

```python
md = output.splice_junctions.metadata
liver_polya_idx = md[
    (md["biosample_name"] == "liver")
    & (md["Assay title"] == "polyA plus RNA-seq")
    & (md["data_source"] == "gtex")
].index[0]

values = output.splice_junctions.values
print(values[42, liver_polya_idx])   # → 0.81
```

Reading: "junction 42 (some donor → acceptor pair in chr22:42M-42.13M) is
predicted to show signal ~0.81 under liver polyA-plus RNA-seq from GTEx."
Higher value indicates a more confident prediction that the junction
would be observed in that biological context.

---

## Operational details

- **Required input lengths.** The model only accepts these region widths:
  16384, 131072, 524288, 1048576 (powers-of-2-ish, up to 1Mb). Other
  lengths raise `ValueError`. Pad or trim accordingly.
- **Per-call latency.** ~1–3s wall-clock for a 131kb interval. Cost scales
  modestly with input length thanks to the shared backbone.
- **Rate limits.** Intentionally not documented; the team's posture per
  the [community FAQ](https://www.alphagenomecommunity.com/t/what-are-the-alphagenome-api-limits/673)
  is "increase parallel workers (default 5) until you see
  `RESOURCE_EXHAUSTED`, then back off." `predict_intervals` and
  `predict_variants` accept a `max_workers` parameter for batching.
- **SDK install.** `pip install alphagenome` (PyPI; not source build).
- **API key.** Set `ALPHAGENOME_API_KEY` in `.env`. Loader pattern matches
  [`research/scripts/zotero_add.py`](../research/scripts/zotero_add.py)'s
  `load_env`.

Minimal usage:

```python
from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models.dna_output import OutputType

model = dna_client.create(API_KEY)
interval = genome.Interval(chromosome="chr22", start=42_000_000, end=42_131_072)
output = model.predict_interval(
    interval=interval,
    requested_outputs=[OutputType.SPLICE_JUNCTIONS],
    ontology_terms=None,
)
junctions = output.splice_junctions.values   # shape (n_junctions, 367)
metadata = output.splice_junctions.metadata  # DataFrame describing the tracks
```

See [`scripts/alphagenome_spike.py`](../scripts/alphagenome_spike.py) for a
runnable end-to-end example with both interval and variant calls.

---

## For the predicted-normal filter

Use case: filter tumor candidate junctions by checking whether they're
predicted to have non-trivial signal in *any* healthy human tissue
([Issue #212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212)).

Sketch:

```python
# 1. Restrict to GTEx tracks (healthy human tissues).
md = output.splice_junctions.metadata
gtex_cols = md.index[md["data_source"] == "gtex"]   # 54 tracks
gtex_signals = output.splice_junctions.values[:, gtex_cols]  # (n_junctions, 54)

# 2. Aggregate across tissues — pan-tissue "any signal" check.
max_signal_per_junction = gtex_signals.max(axis=1)   # (n_junctions,)

# 3. Threshold (TBD by Scientist; coarse default to start with).
predicted_normal = max_signal_per_junction > THRESHOLD
```

Reading: a tumor junction matching one of the model's `predicted_normal`
candidates is plausibly normal splicing → filter out. A tumor junction the
model predicts absent across all 54 GTEx tracks (low signal everywhere) is
a better tumor-specific candidate.

The threshold is the open design question for [Issue #212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212)
and [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224)
— Scientist's call after the validation experiments in [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203).

For the experiment-specific tool choice (which API method to use, what to
benchmark against), see Scientist's design in [#203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203)
and the implementation in [#224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224)
— this primer is intentionally scoped to "what AlphaGenome does," not "how
to design experiments around it."
