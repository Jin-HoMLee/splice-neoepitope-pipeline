# Issue 681 - reanalyze public HLA-ligandome MS for splice-junction-derived neoepitopes

Assembles the AC-1 splice-neoepitope panel and intersects it against public HLA-ligandome atlases.

```bash
research/.venv/bin/python research/experiments/issue_681_ligandome_ms_reanalysis/build_panel.py
research/.venv/bin/python research/experiments/issue_681_ligandome_ms_reanalysis/recover_presented.py
```

## The panel (89,415 unique peptides)

| stratum | n unique | source |
|---|---|---|
| `ours` | 2,156 | our pipeline, patients 001 + 002 (`mhc_presentation.tsv`, fetched from R2) |
| `snaf_pred` | 87,258 | SNAF-T predicted splice neoantigens: Data S1 (TCGA-SKCM), Data S1 (Van Allen), Data S3 (ovarian) |
| `snaf_ms` | 4,715 | SNAF peptides **detected in HLA-ligandome MS** (Data S2, melanoma cohort) |

All strata are canonical-alphabet 9/10-mers. `ours` is disjoint from the SNAF strata (different tumor
types: gastric vs melanoma/ovarian). `snaf_ms` is almost entirely contained in `snaf_pred` (4,714/4,715).

## Why `snaf_ms` is the load-bearing stratum

It is not extra coverage - it is the **ceiling control**. Those peptides are splice-junction-derived
*and* known to be MS-detectable, because SNAF detected them in real immunopeptidomics. So they answer
the question that a plain null cannot: **can a reference atlas report this peptide class at all?**

This matters because these atlases were built by searching MS spectra against a canonical proteome,
which raises the possibility that a junction-spanning peptide is *structurally* unreportable - in which
case every null we score is guaranteed by the method rather than observed in the biology. It is not the
case: the atlases return 79 and 23 hits against the SNAF panel, so the ceiling is real, just low.

Yesterday's control (`HLA Ligand Atlas n SysteMHC = 171`) proved only that **exact string matching works
across sources**. That is a control on the machinery, not on the ceiling. See
`feedback_search_key_must_not_mirror_the_gate.md`: check the CEILING, not just the floor.

## Results

**The axis that decides everything is the reference's SEARCH DATABASE, not its tissue.** This was not
obvious going in - the natural framing is benign-vs-tumor - and it only became visible once caAtlas was
added as a third reference.

| reference | tissue | search database | n (class-I) |
|---|---|---|---|
| HLA Ligand Atlas | **benign** | canonical proteome | 95,318 |
| caAtlas | **tumor** | **100% UniProt-mapped** (measured: 557,450/557,450 records) | 195,217 |
| SysteMHC non-UniProt | **tumor** | **non-canonical by construction** | 8,911 |

| panel | vs HLA Ligand Atlas (canonical) | vs caAtlas (canonical) | vs SysteMHC non-UniProt (**non-canonical**) |
|---|---|---|---|
| `snaf_pred` (87,258) | 79 (0.091%) | 188 (0.215%) | 29 (0.033%) |
| `snaf_ms` - **MS-confirmed** (4,715) | **0** (exp. 4.3) | **0** (exp. 10.1) | **4** (0.085%) |
| `ours` (2,156) | 0 (exp. 2.0) | 0 (exp. 4.6) | 0 (exp. 0.7) |

**Read the middle row.** These are 4,715 splice-junction peptides that real immunopeptidomics *did*
detect - they exist, they are presented, SNAF found them in melanoma MS. They score **zero against both
canonical-search references** (Fisher p = 0.035 and p = 8.8e-05 against the predicted-set base rate) and
are found **only** in the one non-canonical reference. That is not biology. That is the reference's
search database deciding in advance what it is able to report: a peptide spanning a novel junction is
not in UniProt, so a UniProt-restricted search cannot return it, however abundantly it is presented.

**Consequence for AC-2, which is the real finding of this experiment:** an exact-match intersection
against a canonical-search atlas has a **structural ceiling of ~zero for genuinely novel
splice-junction peptides**. Every hit such an atlas *does* return (the 79, the 188) is by definition a
peptide that is **not novel** - a predicted "neoantigen" whose sequence coincides with a canonical
protein. Those hits are worth having (they are tumor-specificity failures of the prediction step, see
below), but they are the opposite of what AC-2 was reaching for. The only productive reference class is
one built from a **non-canonical / open search** - SysteMHC's non-UniProt set today, or re-searching raw
spectra ourselves, which is exactly what the AC-10 Comet stretch goal proposes.

1. **Our zero is underpowered, not a specificity pass.** At SNAF's empirical base rate (0.091%) our
   2,156-peptide panel expects 1.95 hits against the benign atlas, and `P(0 | 1.95) = 0.14`. A zero is
   fully consistent with our `tumor_exclusive` filter being *exactly as leaky as SNAF's*. The panel is
   ~40x too small to tell. (This retracts the "specificity pass" reading first posted on the Issue,
   2026-07-13.) Against caAtlas the expected count is 4.6 and `P(0) = 0.010`, which *looks* significant -
   but read it with the ceiling in mind: it says our peptides coincide with canonical sequences *less*
   often than SNAF's do, which is a statement about novelty, not about specificity. Some of that may
   simply be the low-complexity junk in our panel ([Issue #1147](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1147)) matching nothing anywhere.
2. **The benign-vs-tumor dissociation I first reported does not survive the third reference.** I had read
   `snaf_ms` depletion-in-benign + enrichment-in-tumor as a specificity contrast. caAtlas is *tumor* and
   yields **0** - so the contrast was never benign-vs-tumor. It is canonical-vs-non-canonical, and the
   enrichment leg was never significant anyway (2.6x, p = 0.087, n = 4). Superseded.
3. **The 79 (and 188) canonical-atlas hits are tumor-specificity failures of the prediction step:**
   predicted "neoantigens" that are in fact presented on healthy tissue, or that are simply not novel.
   Genuine ligand-like sequences (mean Shannon entropy 2.67 bits), not low-complexity artifacts. A natural
   feed for the AC-9 decoy seed.

## AC-4 - the recovered presentation tier (29 rows), and its null

`recover_presented.py` -> `outputs/presented_untested_rows.tsv`, `outputs/recovery_summary.json`.

Every row is conditioned on the one reference that *can* report this peptide class (SysteMHC's
non-UniProt set, 8,911 class-I peptides). The rows land as `label=untested` /
`tier=presentation-prevalence` / `assay_context=prevalence_only` and **must never be pooled into the
functional positive set** (AC-5, AC-6): they are evidence of *presentation*, not of immunogenicity.
Nobody T-cell-tested these.

**The control is the content.** An exact 8-11mer match between an 87k panel and an 8.9k atlas is not
self-evidently meaningful - class-I peptides are short, drawn from a 20-letter alphabet, and strongly
composition-biased at anchor positions, so some collisions are expected from composition alone. Two
independent decoys (200 replicates, seed 681):

| panel | n | hits | decoy(shuffle) mean | p | decoy(reverse) |
|---|---|---|---|---|---|
| `snaf_pred` | 87,258 | **29** (0.033%) | 0.010 | 0.005 | 3 |
| `snaf_ms` (MS-confirmed) | 4,715 | **4** (0.085%) | 0.000 | 0.005 | 0 |
| `ours` | 2,156 | 0 | 0.000 | 1.000 | 0 |

`shuffle` permutes each panel peptide's residues (preserves length **and** composition exactly, destroys
only sequence order) - the strict null. `reverse` is the conventional proteomics decoy and is weaker
(it preserves order structure). The observed hits survive both, but note the reverse decoy recovers
**3**, so the honest read is ~29 hits against a background of ~3, not against zero. `p = 0.005` is the
floor of a 200-replicate empirical test (the add-one estimator `(ge+1)/(replicates+1)` = 1/201), not a
precise p-value.

**The seed only pins the null if the RNG is consumed in a fixed order.** `decoy_null` iterates
`sorted(panel)`, never the set itself: Python hash-randomizes `str` per process, so set-iteration order
varies run to run and the decoy stream drifted even with `seed=681` - which is exactly how an earlier
draft of this README quoted a `shuffle` mean of 0.01 while the JSON it was derived from said 0.02. Two
runs under different `PYTHONHASHSEED` values now produce byte-identical output. (Caught in review of
[PR #1172](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/1172).)

**Our own panel contributes nothing to this tier, and that is a power statement, not a result.** At
SNAF's observed rate against this reference (0.033%), our 2,156 peptides expect **0.72** hits. A zero is
what you get from a panel this size whether our junctions are presented or not. Same shape as the
retracted "specificity pass" - the honest label is *underpowered*.

## AC-7 - per-allele coverage: the skew is NOT demonstrably broken

The tempting headline was available and is wrong. The recovered set spans **25 distinct alleles over 28
allele-assignments**, with A\*02:01 at just **1/28 (3.6%)** - which looks like the A\*02:01 dominance of
classical immunopeptidomics being broken. It is not, for three reasons, and the third is fatal:

1. **The reference's own A\*02:01 share is 6.3%** (188/2,980 unique `(peptide, allele)` pairs). At that
   background the recovered set *expects* **1.77** A\*02:01 assignments and we observe **1**. There is no
   depletion here - the recovered distribution is statistically indistinguishable from the reference's.
   (Background and recovered are counted in the **same unit** - unique `(peptide, allele)` pairs, not raw
   SysteMHC records - or the expectation would compare two different quantities.)
2. **The reference is allele-diverse by construction.** SysteMHC leans on mono-allelic cell-line data,
   which spreads across 103 distinct alleles by design. Allele breadth is a property of how the atlas was
   *built*, not of splice peptides. And 28 assignments spread over a 103-allele reference is ~28 draws
   from a 103-category multinomial, so *almost every cell is expected to be 0 or 1*: the "1 each" table
   is what pure sampling noise looks like, and carries essentially no allele-specific signal.
3. **Most of the reference carries no allele at all** (only 2,325 of its peptides have one; the rest are
   `NA` or `unclassified`), and **11 of our 29 recovered peptides carry none**. The coverage table is
   computed on a minority of a minority.

So a per-allele table here mostly **restates the reference's allele composition** - the same failure shape
as `feedback_search_key_must_not_mirror_the_gate.md` (a key that mirrors the gate can only confirm it).
**AC-7's answer is therefore: we cannot state that the A\*02:01 skew is broken.** Answering AC-7's actual
question needs allele assignments from a re-search we control (AC-10), not from a borrowed column.

## SpliceMutr publishes no peptide set

AC-1 originally read "our pipeline + SNAF/SpliceMutr public outputs". The SpliceMutr half of that premise
is false, verified against four independent artifacts (not one failed fetch): the data-availability
statement names only recount3 junctions + code; the GitHub repo is code; the journal supplements are a
subtype key and a Wilcoxon table, with zero peptide-shaped tokens in the docx; and the 228 MB Zenodo
archive is `.git` history plus scripts, including the `extract_peptides.R` that *would* produce peptides,
with none of its output. SpliceMutr's published output is a per-sample *splicing antigenicity score*.
Recovering peptides would mean running their pipeline over TCGA recount3 junctions + HLA genotypes: a
compute project, not an ingestion.

## Data provenance (`data/` - gitignored, rebuild with `bash fetch_data.sh`)

`data/` is covered by the repo-wide `data/` ignore rule (`.gitignore:231`), which is the right
convention here: `outputs/` is the committed artifact, `data/` is re-fetchable input. `fetch_data.sh`
rebuilds it.

| file | source |
|---|---|
| `hla_ligand_atlas_aggregated.tsv.gz` | HLA Ligand Atlas `rel/2020.12/aggregated.tsv.gz` - 223,246 peptides, **benign** tissue (95,318 in the class-I 8-11mer band) |
| `systemhc_nonuniprot.html` | SysteMHC Atlas v2.0 `/Non-UniProt` page - 20,000 rows, 15,130 unique peptides (8,911 class-I), **tumor** samples |
| `patient_00{1,2}_mhc.tsv` | our `results/patient_00{1,2}/predictions/mhc_presentation.tsv`, from R2 |
| `hgnc_symbols.txt` | HGNC complete set, filtered to **protein-coding + Approved** - the gene list driving the caAtlas harvest |
| *(by hand)* `scitranslmed.ade2886_data_s1_to_s15.xlsx` | SNAF supplementary workbook, DOI `10.1126/scitranslmed.ade2886`. **Not fetched and not committed** (29 MB, publisher's file). `fetch_data.sh` documents where to put it; both scripts fail without it. |

**The caAtlas provenance chain** (`harvest_caatlas.py`), because the committed `.gz` is the AC-2
ceiling-control input and half the headline finding rests on it:

```
hgnc_symbols.txt (19,273 protein-coding genes)
  -> harvest_caatlas.py            gene-by-gene fetch, resumable
  -> outputs/caatlas_raw.jsonl     557,450 peptide rows (gitignored, 60 MB)
  -> harvest_caatlas.py --distill-only
  -> outputs/caatlas_peptides.txt.gz   195,217 unique class-I peptides (committed, 0.9 MB)
```

The `.gz` is committed deliberately (under the 10 MB band, and network-derived from a live site so it
is not offline-regenerable); the raw JSONL is not (60 MB, single consumer). The distill step is written
with `mtime=0` so re-running it on identical input yields an identical file rather than a spurious diff.
**Verified 2026-07-14:** re-distilling the raw harvest reproduces the committed peptide set exactly.

**Two traps in the SysteMHC scrape, both of which produced wrong numbers on the first pass (2026-07-13):**

1. **A truncated download looks like a complete one.** The first fetch cut off mid-`<td>` at 3.2 MB and
   silently yielded 3,870 peptides instead of 8,911. The complete file is 7,400,953 bytes and ends with
   `</html>`. Check the terminator.
2. **Do not sweep every `<td>`.** The peptide is column 2 of a fixed 10-column row. An indiscriminate
   cell sweep pulls metadata words that happen to be valid amino-acid strings into the peptide set -
   observed: `METASTATIC` and `PANCREAS`.

The served table is exactly **20,000 rows**, a round number that is very likely a server-side cap, so
even 15,130 is a lower bound on what SysteMHC exposes - and still far below the **78,959** non-UniProt
peptides claimed in Table 1 of the paper. That gap is a data-availability finding in its own right (AC-3).

SNAF's supplementary workbook is read from local Zotero storage (item `TZGPRIFK`, attachment `WK4DHT6M`)
rather than committed: it is 29 MB and is the publisher's file, not ours to redistribute.

## Background: proteasome-spliced peptides are NOT what we study

[`explainers/reverse_cis_splicing.html`](explainers/reverse_cis_splicing.html) - open in a browser.

A recurring confusion, and one that is baked into at least one reference set we may use (NeoAtlas-Protein
is a proteasome-spliced catalogue): **proteasome-catalyzed peptide splicing (PCPS) is post-translational
peptide ligation**, encoded nowhere in the transcriptome. Our splice-junction-derived neoepitopes are
transcript-encoded. The explainer walks the four ligation topologies and shows why a *forward cis*-spliced
product is sequence-identical to an exon-skipping junction peptide (so MS alone cannot separate them),
while a *reverse cis* product cannot be encoded by any transcript at all. Load-bearing consequence for us:
our peptides have RNA-seq as an orthogonal corroborating assay; PCPS claims have no nucleic-acid check
available, which is why that field's false-discovery dispute is unresolved.

## Still open

- **AC-10 (the main event, promoted from stretch 2026-07-14):** re-search raw spectra against a
  non-canonical database we build from our own junction panel. Everything above converges on this: a
  borrowed atlas either *cannot* report our peptide class (canonical search: structural ceiling ~zero)
  or reports it through a column we do not control (SysteMHC: 17% allele-assigned, 83% of the claimed
  peptides never released). AC-2 and AC-7 both terminate in "you need your own search."
- AC-9 decoy seed: the 79 + 188 canonical-atlas hits (predicted "neoantigens" presented on healthy tissue).
- AC-3 SysteMHC splice-origin fraction, and the ~74k claimed-but-unreleased peptides (author request drafted).
- AC-8 validator enforcement of the `tier`-not-`label` filter discipline, now that 29 presentation rows exist.
