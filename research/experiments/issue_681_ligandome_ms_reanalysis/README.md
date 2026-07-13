# Issue 681 - reanalyze public HLA-ligandome MS for splice-junction-derived neoepitopes

Assembles the AC-1 splice-neoepitope panel and intersects it against public HLA-ligandome atlases.

```bash
research/.venv/bin/python research/experiments/issue_681_ligandome_ms_reanalysis/build_panel.py
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

The two references differ in the way that matters: **HLA Ligand Atlas is benign tissue**, while the
**SysteMHC non-UniProt set is tumor** (Lymphoma 11,020 rows, Metastatic 4,580, Lung cancer 647,
Melanoma 387). So they are not two samples of one thing - they are the negative and positive halves of
a specificity contrast.

| panel | vs HLA Ligand Atlas (**benign**) | vs SysteMHC non-UniProt (**tumor**) |
|---|---|---|
| `snaf_pred` (87,258) | 79 (0.091%) | 29 (0.033%) |
| `snaf_ms` (4,715) | **0** (expected 4.3) | **4** (0.085%) |
| `ours` (2,156) | 0 (expected 2.0) | 0 (expected 0.7) |

1. **Our zero is underpowered, not a specificity pass.** At SNAF's empirical base rate (0.091%) our
   2,156-peptide panel expects 1.95 hits against the benign atlas, and `P(0 | 1.95) = 0.14`. A zero is
   fully consistent with our `tumor_exclusive` filter being *exactly as leaky as SNAF's*. The panel is
   ~40x too small to tell. (This retracts the "specificity pass" reading first posted on the Issue,
   2026-07-13.)
2. **One leg of the dissociation holds, one does not.** MS-confirmed splice antigens are **depleted** in
   the benign atlas (0 vs 4.3 expected, Fisher p = 0.035). They are enriched 2.6x in the tumor reference,
   but at **p = 0.087 that leg is not significant** (n = 4 hits). Direction matches the biology; the
   evidence does not yet carry it. Do not report this as a double dissociation.
3. **The 79 benign-atlas hits are tumor-specificity failures of the prediction step:** predicted
   "neoantigens" that are in fact presented on healthy tissue. Genuine ligand-like sequences (mean Shannon
   entropy 2.67 bits), not low-complexity artifacts. A natural feed for the AC-9 decoy seed.

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

## Still open

- caAtlas (tumor immunopeptidome) - per-gene harvest, the reference most likely to carry real hits.
- AC-4 tier rows: the 23 SysteMHC hits are the presentation-prevalence candidates.
