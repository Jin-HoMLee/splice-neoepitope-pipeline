# ORF-stretch FASTA emitter (genome-wide junction-spanning MS search DB) - design

Issue: [#1204](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1204) (`feat(translation)`), prerequisite for [#1176](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1176) AC-6.
Date: 2026-07-21.
Author: Developer.

## Goal

Produce a FASTA of the candidate "mini-proteins" that tumor-specific splice junctions would translate to, so a nonspecific (no-enzyme) mass-spectrometry search engine can search real immunopeptidome spectra for junction-derived peptides.
This FASTA is the missing input for the [#1176](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1176) re-search; it delivers **evidence of presentation**, not a binding prediction.

## Why a new emitter (not the existing k-mer path)

`translate_peptides.py` was built for MHCflurry: it slides fixed-length 8/9/10-mer windows across the breakpoint and writes TSV rows.
An MS engine needs the opposite - variable-length protein sequences (stop-to-stop ORF stretches) in FASTA, because the engine does its own in-silico digestion.
Feeding it pre-cut k-mers defeats the nonspecific search and silently drops 11-mers and any longer peptide the spectra contain.
Two mismatches: (1) fixed k-mers vs. whole ORF stretches, TSV vs. FASTA; (2) 27-nt flanks cap peptides at 10-mers, but class-I ligands run to 11 (need 30-nt flanks).

## Architecture

Reuse the existing two-stage split; replace only the second stage.

1. **Reuse `assemble_contigs.py`** at 30-nt flanks (60-nt contig) via a decoupled invocation.
   Gives the `tumor_exclusive` filter, junction dedup, `bedtools getfasta`, soft-mask handling, strand handling, and the coord-carrying header for free.
2. **New `workflow/scripts/orf_fasta_from_contigs.py`** (+ a new Snakemake rule; both siblings are rule-invoked).
   Reads the 60-nt contig FASTA, translates all three frames of the full contig, splits at stop codons, keeps the ORF stretch that **crosses the breakpoint**, and writes it as a FASTA mini-protein.

### ORF-stretch extraction (the load-bearing algorithm)

For a 60-nt contig, the breakpoint boundary sits at nt 30 (0-based): upstream = `[0,30)`, downstream = `[30,60)`.
For each frame `f` in `{0,1,2}`:

- Translate `contig[f:]` codon by codon; amino acid `aa[i]` maps to nt span `[f+3i, f+3i+3)`.
- Split the AA string into maximal stop-free runs (ORF stretches); record each run's nt span `[f+3*c_start, f+3*c_end)`.
- **Keep a stretch iff `nt_start < 30 < nt_end`** - it straddles the breakpoint, so it contains junction-novel sequence.
  A stretch entirely upstream (`nt_end <= 30`) or entirely downstream (`nt_start >= 30`) is canonical and redundant with the concatenated proteome, so it is dropped.
- There is at most one crossing stretch per (contig, frame): a single breakpoint sits inside at most one stop-free run.
  If a stop falls across the breakpoint region, that frame emits nothing.
- **Drop a crossing stretch that contains `X`** (ambiguous residue - not a valid MS DB entry) and any stretch shorter than `min_peptide_len` (default 8, the minimum class-I ligand).
  Log direct/dropped-X/dropped-short counts at INFO.

### Header format

`>{junction_id}|{frame}|{chrom}:{start}-{end}:{strand}`, parsed from the contig header (`{junc_id}|{chrom}:{start}-{end}:{strand}|{sample_type}`, `sample_type` dropped).
Headers stay traceable to `junction_id`/coords for the [#1176](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1176) AC-6 entrapment + RNA-seq cross-check.

## Decisions

- **Fork 1 - flank width: decoupled (A).**
  The ORF path uses its own 30-nt flank param; `translation.peptide_lengths` stays `[8,9,10]`, so the production MHCflurry DAG and its fixtures are untouched.
  Realized at the DAG level as a separate 30-nt contig-assembly invocation feeding the emitter (a distinct contig output, not the 27-nt one the MHC path consumes).
- **Fork 2 - AC-3 verification: both.**
  A structural check in CI (valid FASTA parse, unique headers, valid AA alphabet, clean concatenation with a canonical proteome) **plus** a one-time local **Sage** smoke (MIT, native arm64-macOS, CPU-only) with `cleave_at=""` against a tiny synthetic mzML, recorded in the Test plan per the live-integration-smoke rule.
  Proves real loadability at $0 without the Courcelles spectra (those are AC-6).
- **Optional target, not default DAG.**
  The ORF FASTA is an opt-in target under a new config block (`ms_search_db: {flank_nt: 30, min_peptide_len: 8}`), not in `rule all`.
- **Promotable internals.**
  `translate` + `stop-split` are written as clean, consumer-agnostic functions so the emitter can later be promoted to the canonical peptide-generation intermediate ([#1255](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1255)) as a wiring change, not a rewrite.

## Genome-wide capability (AC-1)

No code is chr22-bound: `assemble_contigs.py` and the emitter already take a `--genome-fasta` and a junction TSV, so pointing them at a genome-wide reference is the capability.
AC-1 is demonstrated with a **multi-chromosome fixture** (junctions on several chromosomes against a multi-chr reference), proving the path is not chr22-scoped - a cheap, falsifiable capability check, not a full genome-wide alignment run.
The *actual* genome-wide DB stays RNA-seq-bounded to the cohort (Courcelles `GSE312236`) and is built at [#1176](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1176) AC-6, not here.

## Testing

- **Unit** (`workflow/tests/`): ORF-stretch extraction with matched-pair falsifiers - clean crossing ORF -> one stretch/frame; stop at the breakpoint -> nothing; upstream-only -> dropped; `X` -> dropped; min-len filter; header format; minus-strand contig.
- **Structural CI check**: run the emitter on the multi-chr fixture, assert unique headers, valid AA alphabet (no `*`/`X`), and duplicate-free concatenation with a small canonical proteome FASTA.
- **Live Sage smoke** (one-time, pre-merge, Test plan): concatenate the emitted junction FASTA + a tiny canonical proteome, run Sage nonspecific against a synthetic mzML, confirm the search completes (exit 0). Fixtures + recipe under `research/experiments/issue_1204_orf_fasta/`.

## Out of scope

- The actual cohort DB assembly (canonical + cRAP + decoys) and the MS re-search - [#1176](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1176) AC-6.
- Promoting the emitter to the canonical peptide-generation intermediate - [#1255](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1255), trigger-gated on MHC-II.
- Any change to the production MHCflurry k-mer path.
