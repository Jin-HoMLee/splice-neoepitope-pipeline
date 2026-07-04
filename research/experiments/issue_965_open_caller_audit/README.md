# Issue #965 - Open-only splice-caller install + license audit (benchmark leaf A)

**Parent epic:** [#679](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/679) (first open head-to-head benchmark of splice-neoantigen callers), arc `immunogenicity-benchmark`.
**Leaf A** of the [#964](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/discussions/964) decomposition; feeds the harness-skeleton leaf B ([#966](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/966)).

**Goal.** Install and smoke-run each open splice-neoantigen caller under the **open-only** constraint, standardizing the MHC-presentation step on **MHCflurry (Apache-2.0)** in place of the blocked NetMHCpan, and produce a per-tool openness/license note plus a **go/no-go list** of which callers actually run open-only.

## Two axes, not one

The audit separates two questions that the issue phrases as one, because for these tools they diverge sharply:

1. **Open-only viability** - is the tool's own license permissive/redistributable, and can its *core* run avoid every paid / non-redistributable / academic-only dependency once presentation is MHCflurry? This is the question the epic needs answered to pick a caller set.
2. **Local-arm64 smoke feasibility** - does it actually run on our current keep-alive target (macOS arm64 / Apple M1 / 8 GB / CPU-only / **no Docker daemon**), where **STAR is VM-only** (its human index build exceeds 8 GB)?

A tool can be a clean open-only GO yet be un-smokeable locally (SNAF: MIT but Linux-x86/Docker-only). Conversely no tool that is a licensing NO-GO becomes viable by being locally installable. Leaf B (the harness) consumes the **open-only GO** set; **where** each runs (local arm64 vs free-GPU/Linux) is a scheduling detail carried in the table below.

## Go/no-go table (headline)

| Caller | License (source-verified) | Presentation step | Non-redistributable binary | **Open-only** | Local arm64 CPU smoke | Runs where |
|--------|---------------------------|-------------------|----------------------------|:-------------:|:---------------------:|------------|
| **splice2neo** | MIT (TRON gGmbH) | none - stops at peptides; feeds our MHCflurry | none | **GO** | **YES** (toy fixtures) | local arm64 |
| **SNAF** | MIT (F. Li) | native `binding_method='MHCflurry'` flag | none | **GO** | NO (Py3.7/TF2.3 linux-x86; amd64 AltAnalyze Docker) | Linux / free-GPU |
| **ASNEO** | Apache-2.0 (bm2-lab) | bundled NetMHCpan bypassed (#566 patch) -> MHCflurry | yes, **stripped/bypassed** by [#566](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/566) | **GO** | PARTIAL (caller runs on a pre-made `SJ.out.tab`; STAR front-end VM-bound) | local arm64 (caller) / VM (alignment) |
| **NeoSplice** | Apache-2.0 code, but vendors NetMHCpan-4.0 + NetMHCIIpan-3.2; README says "academic and non-profit" | hard `--netMHCpan_path`, fixed `.xls` nM parser | **yes, vendored** | **NO-GO** as-is | NO (Python 2.7 EOL, no arm64) | needs a Py3 fork + MHCflurry shim |
| **SINE** | GPL-3.0 (Guo Lab UCSD) | NetMHCpan **+ NetMHCIIpan** load-bearing (PHBR class I **and** II) | none bundled | **NO-GO** as-is | NO (Trinity-in-Singularity; netMHCpan x86) | class-I-only re-plumb on Linux |
| **NeoHunter** | **MIT** (own) + Apache-2.0 (bundled ASNEO) - issue's "academic-only" label is **wrong** | pervasive NetMHCpan **+ NetMHCstabpan** (no MHCflurry equivalent for stability) | **yes, vendored** netMHCpan-4.0 x86 | **EXCLUDE** (operational) | NO (STAR >30 GB index, hg19 tens-of-GB refs, x86) | - (wraps ASNEO anyway) |
| **REAL-neo / SPLICE-neo** | no software released; paper is CC BY-NC 4.0 | hard-wraps NetMHC/NetMHCpan/NetMHCIIpan suite | N/A | **NO-GO** (no installable code) | N/A | - (not distributed) |

**Open-only GO set (feeds leaf B): `splice2neo`, `SNAF`, `ASNEO`** - 3 of 7.
**NO-GO / EXCLUDE: `NeoSplice`, `SINE`, `NeoHunter`, `REAL-neo/SPLICE-neo`** - 4 of 7.

## Corrections to the issue / epic framing (source-verified)

Three assumptions in [#965](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/965) / [#679](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/679) did not survive the audit and are corrected here:

- **NeoHunter is NOT "free academic use."** Its own code is **MIT** (the license detector missed it because the file is misspelled `LISCENCE.md`) and the bundled ASNEO module is Apache-2.0 - both permissive. The disqualifier is **operational** (vendored netMHCpan-4.0 x86 binary; hard-wired NetMHCstabpan, which MHCflurry cannot replace; hg19 + tens-of-GB references + STAR's >30 GB index), not licensing. NeoHunter's splicing branch is literally a wrapper around ASNEO, so standalone ASNEO is the better comparator.
- **`splice2neo` is a library, not an end-to-end caller.** It ingests an upstream junction caller's output (LeafCutter / SplAdder / RegTools / SpliceAI / ...) and emits peptide/context sequences; it has **no** MHC step. So it slots into the benchmark as a *junction-to-peptide component* to compare against our own contig/translation stage, paired with our MHCflurry behind it - not as a standalone pipeline peer. Default build is hg19 (needs `liftOver` to our hg38).
- **`REAL-neo / SPLICE-neo` has no released code.** SPLICE-neo is a module inside Mayo Clinic's in-house REAL-neo framework (Wickland et al., JITC 2024); the paper ships no code/Zenodo/GitHub link and hard-wraps the NetMHC suite. There is nothing to install or license-audit. The "R spliceneo/splice-neo" the issue also names resolves to `splice2neo` (audited above).

Also noted for leaf B: NeoSplice uses **STAR + msbwt/msbwt-is (Burrows-Wheeler)**, not MOSAIK/bwa as the issue text guessed.

## Per-tool notes

Each fact below was verified against the tool's actual `LICENSE` / README / env files (fetched 2026-07-04), not recall. Full per-tool fact sheets: [`audit_notes.md`](audit_notes.md). Install recipes: [`install/`](install/).

### splice2neo - GO, locally smokeable
`TRON-Bioinformatics/splice2neo` v0.6.14, R >= 4.4, **MIT**. No NetMHCpan, no bundled binaries, no MHC step (stops at peptides). Bioconductor stack (GenomicRanges/Biostrings/...), builds CPU-only on arm64. Ships toy fixtures (`toy_junc_df`, `toy_transcripts`, `toy_cds`) so a real smoke needs no external caller. Reference: `BSgenome.Hsapiens.UCSC.hg19` (~700 MB). Install recipe: [`install/splice2neo_install.R`](install/splice2neo_install.R). Smoke: [`install/splice2neo_smoke.R`](install/splice2neo_smoke.R).

### SNAF - GO on license, Linux/Docker for the run
`frankligy/SNAF`, Python, **MIT**. MHC step is a first-class flag: `snaf.initialize(..., binding_method='MHCflurry')` calls MHCflurry in-process with zero NetMHCpan dependency (both engines register under the same internal attribute, so downstream is engine-agnostic). No bundled binaries. Blocked locally by hard Linux-x86 pins (Python 3.7.12, TensorFlow 2.3.0, PyMC3) and an **amd64-only AltAnalyze Docker image** for junction extraction - and we have no Docker daemon. 2.72 GB reference bundle. Runs on Linux / free-GPU only; smoke there needs >=2 BAMs (AltAnalyze is cohort-level). Install recipe (Linux): [`install/snaf_install.md`](install/snaf_install.md).

### ASNEO - GO, largely solved by #566
`bm2-lab/ASNEO`, Python, **Apache-2.0** code. Vendors NetMHCpan-4.0/NetCTLpan-1.1/pepmatch (non-redistributable, x86) inside `src/software.tar.gz`, but [#566](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/566) already ships a tested `apply_optionB_patch.py` that bypasses all of them (stops ASNEO at its normal-subtracted candidate peptides, repoints bedtools to PATH so the tarball is never extracted) plus an arm64-validated `asneo_env.yml` (python 3.8, biopython 1.79). The remaining open-only work is only the MHCflurry scoring of `putative_peptide.txt`. Input is a STAR `SJ.out.tab` (hg19), so the *front-end alignment* is VM-bound, but the **caller-proper runs CPU-only on arm64** given a junction table. See [`install/asneo_notes.md`](install/asneo_notes.md) and the #566 experiment dir.

### NeoSplice - NO-GO as-is
`pirl-unc/NeoSplice`, Python **2.7** (EOL), **Apache-2.0** code but **vendors** netMHCpan-4.0 + netMHCIIpan-3.2 Linux-x86 binaries and its README narrows use to "academic and non-profit." NetMHCpan is a hard code dependency: `run_netMHCpan()` shells out with a required `--netMHCpan_path`, and step 9 parses NetMHCpan's `.xls` at fixed nM-affinity column offsets - swapping to MHCflurry means rewriting both the producer call and the parser (and re-expressing the 500 nM affinity cutoffs as presentation percentiles). No osx-arm64 path (2018-pinned C-extension stack). Viable only as a **Python-3 fork with an MHCflurry shim**, run in emulated Linux Docker - a port, not a drop-in.

### SINE - NO-GO as-is
`GuoLabUCSD/SINE` (Splice Isoform Neoantigen Evaluator), Python 3 + Bash, **GPL-3.0**. No bundled binaries (clean repo), but NetMHCpan **and NetMHCIIpan** are load-bearing: it computes PHBR from both class-I and class-II `.xls` output. MHCflurry (Apache-2.0) can replace **class I only** - it has **no class-II predictor**, and class-II binding is a core selling point of the tool - so a pure open-only run is impossible without dropping half the method and rewriting the predictor/PHBR wrapper. Second blocker: de-novo assembly runs **Trinity via Singularity/Apptainer** (Linux containers) on our no-Docker arm64 box. GPL copyleft also constrains any code we adapt.

### NeoHunter - EXCLUDE
`XuegongLab/NeoHunter`, Python, **MIT** (own) + Apache-2.0 (bundled ASNEO). A genuine multi-source caller (SNV/indel/fusion/**splicing**), with the splicing branch delegating to bundled ASNEO. Excluded on operations, not license: vendored netMHCpan-4.0 x86 binary, pervasive NetMHCpan **plus NetMHCstabpan** (stability - no MHCflurry equivalent, so that filter can't be swapped, only dropped, changing ranking), hg19-only with tens-of-GB references (GATK bundle, STAR-Fusion CTAT), and STAR's >30 GB index build. None fit the arm64 / 8 GB / CPU-only / MHCflurry target. Its splicing signal is just ASNEO, which we cover directly.

### REAL-neo / SPLICE-neo - NO-GO (no code)
Mayo Clinic in-house framework (Wickland et al., *J Immunother Cancer* 2024). No public repository, no software license, no code-availability statement - the paper wraps the paid NetMHC/NetMHCpan/NetMHCIIpan suite via pVACtools. Nothing installable.

## Openness documented end-to-end

- **Presentation:** MHCflurry 2.2.0 (Apache-2.0), already in the `snakemake` env; the open standard replacing NetMHCpan across every GO caller.
- **GO callers:** splice2neo (MIT), SNAF (MIT), ASNEO (Apache-2.0). No non-redistributable artifact is installed or shipped for any of them (ASNEO's bundled binaries are bypassed, never extracted).
- **Data:** the leaf-B/C smoke data (chr22 fixtures; SG-NEx AWS Open Data; recount3/Snaptron; GENCODE) are all open - access verified in [#964](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/discussions/964).

## Status + what feeds leaf B

- **splice2neo: installed + smoked on arm64 (PASS).** `install/splice2neo_install.R` builds the MIT package + Bioconductor deps + hg19 BSgenome from source on this M1 (see the toolchain gotcha below); `install/splice2neo_smoke.R` on the bundled toy fixtures gives **17/17 junctions -> a context sequence -> 15 mutated proteins -> 14 junction `peptide_context` neoepitope candidates (5 frame-shift; 3 junctions yield no in-frame peptide)** - the exact artifact that feeds MHCflurry. The committed output slice is [`outputs/splice2neo_smoke_out.tsv`](outputs/splice2neo_smoke_out.tsv); the `INSTALL_OK` / `SMOKE_OK` console logs are gitignored (`*.log`) and kept locally.
- **ASNEO:** open-only path validated (env built + patch tested) in [#566](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/566); the caller-proper runs CPU-only on arm64 given a pre-made hg19 `SJ.out.tab`. A local toy-`SJ.out.tab` run + the MHCflurry-scoring notebook are the remaining local proofs (fast follow-up, #566 territory).
- **SNAF:** open-only GO confirmed; its smoke is deferred to the Linux/free-GPU leaf (no local Docker daemon, no arm64 build).
- Leaf B (harness) should target the **3 GO callers**, scheduling SNAF's runs on Linux/free-GPU and running splice2neo (+ the ASNEO caller-proper) locally.

### arm64 toolchain gotcha (reproducibility note)

On this box (homebrew R 4.6.0, Bioconductor 3.23), **no CRAN/Bioconductor arm64 binaries exist yet**, so the whole R stack builds from source - and the default compile line omits `-isysroot`, so Apple clang cannot find libc++ headers and every C++ source package dies with `fatal error: 'cmath' file not found`. `splice2neo_install.R` fixes this by pointing `CPLUS_INCLUDE_PATH` at the active SDK's `usr/include/c++/v1`. A more portable alternative recipe is a conda-forge R env (prebuilt osx-arm64 binaries, no compilation) - preferred if this build is ever reproduced on a different toolchain.

## Files

| File | Role |
|------|------|
| `README.md` | this audit + go/no-go table |
| `audit_notes.md` | full per-tool source-verified fact sheets |
| `install/splice2neo_install.R` | splice2neo + Bioconductor deps + hg19 BSgenome (arm64) |
| `install/splice2neo_smoke.R` | toy-fixture smoke |
| `install/snaf_install.md` | SNAF Linux/Docker recipe (open-only, MHCflurry flag) |
| `install/asneo_notes.md` | ASNEO open-only recipe pointer (#566 patch + env) |
| `outputs/` | committed smoke output slice (`splice2neo_smoke_out.tsv`); build/smoke `*.log` are gitignored, kept locally |
