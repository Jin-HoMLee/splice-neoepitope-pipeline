# Issue #965 - full per-tool fact sheets (source-verified 2026-07-04)

Every license / dependency claim below was verified against the tool's actual `LICENSE`, README, and env/build files (fetched 2026-07-04), not recall. Verdicts are summarized in [`README.md`](README.md).

---

## splice2neo - GO (library, locally smokeable)

- **Repo:** `github.com/TRON-Bioinformatics/splice2neo` v0.6.14 (GitHub-only; not CRAN/Bioconductor). Paper: Lang et al. 2024, *Bioinformatics Advances* 4(1):vbae080 (PMID 38863673).
- **Language:** R package, `Depends: R (>= 4.4)`. Heavy Bioconductor `Imports` (IRanges, GenomicRanges, GenomicFeatures, Biostrings, rtracklayer, liftOver, BSgenome.Hsapiens.UCSC.hg19) + tidyverse.
- **License: MIT.** `DESCRIPTION`: `License: MIT + file LICENSE`; `LICENSE`: `YEAR: 2021 / COPYRIGHT HOLDER: TRON gGmbH`; `LICENSE.md` = verbatim MIT ("...to deal in the Software without restriction..."). GitHub API shows `NOASSERTION` only due to the two-file layout; text is unambiguously MIT.
- **Scope:** a junction -> peptide **library**, not a caller. Ingests upstream junction/effect predictor output (SpliceAI, CI-SpliceAI, MMsplice, Pangolin; LeafCutter, RegTools, SplAdder, IRFinder, SUPPA2, STAR, StringTie), normalizes to a unified junction format, modifies the transcript, translates to peptide/context sequence. Does **not** call junctions itself and does **no** MHC prediction.
- **NetMHCpan:** none - `gh search code` for `netmhc`/`mhcflurry` = 0 hits. MHCflurry is bolted on downstream (feed its peptide strings to our existing MHCflurry stage); nothing to swap inside the tool.
- **Bundled non-redistributable binaries:** none. Pure R + Bioconductor.
- **Reference data:** `BSgenome.Hsapiens.UCSC.hg19` (~700 MB). Default build is hg19 -> needs `liftOver` to our hg38. Ships toy fixtures (`toy_junc_df`, `toy_transcripts`, `toy_cds`).
- **arm64 CPU:** feasible - standard CRAN/Bioconductor, no GPU/x86. (On this box, R 4.6.0 / Bioconductor 3.23 have no prebuilt arm64 binaries, so deps compile from source - slow but clean.)
- **Smoke:** `add_context_seq(toy_junc_df, transcripts=toy_transcripts, size=400, bsg=hg19)` on bundled toy data - no external data needed. See [`install/splice2neo_smoke.R`](install/splice2neo_smoke.R).
- **Verdict: GO** (library, not a caller). Slots into the benchmark as a junction->peptide component vs our contig/translation stage.

---

## SNAF - GO on license; Linux/Docker for the run

- **Repo:** `github.com/frankligy/SNAF`. Paper: Li et al., *Sci Transl Med* 2024, DOI 10.1126/scitranslmed.ade2886.
- **Language:** Python, pinned **3.7.12**. Stack: TensorFlow 2.3.0, MHCflurry 2.0.5, PyMC3 3.11.2 / theano-pymc.
- **License: MIT.** `LICENSE` (raw, main): "Copyright (c) 2021 Guangyuan(Frank) Li ... Permission is hereby granted, free of charge ... without restriction ...". No non-commercial clause. (Optional external netMHCpan/TMHMM carry their own academic licenses but are avoidable.)
- **NetMHCpan:** optional. `binding_prediction(self, hlas, binding_method=None)` selects `'netMHCpan'` (external subprocess, needs a user `software_path`) vs `'MHCflurry'` (in-process Python, the documented default). Both register under the same internal `netMHCpan_el` attribute -> downstream is engine-agnostic. **Swap = a config flag** (`initialize(..., binding_method='MHCflurry')`), no code patch.
- **Bundled non-redistributable binaries:** none (verified via git tree). Large blobs are open reference/data + SNAF's own DeepImmuno CNN weights.
- **Reference data:** 2.72 GB bundle (`altanalyze.org/SNAF/download.tar.gz`, Content-Length 2721440081). **AltAnalyze** (junction extraction from BAM) ships only as a **linux/amd64 Docker image** `frankligy123/altanalyze:0.7.0.1`.
- **arm64 CPU:** no native path - Python 3.7 + TF 2.3.0 are Linux-x86 only (no osx-arm64 conda 3.7; no arm64 TF 2.3), and AltAnalyze is amd64-container-only. Runs only under Docker/amd64 emulation or on Linux. This box has **no Docker daemon**.
- **Smoke:** no tiny fixture ships; AltAnalyze needs >=2 BAMs (cohort-level). Minimal proof = run AltAnalyze on 2 small BAMs (e.g. our chr22 SRR9143066/65) -> `counts.original.pruned.txt` -> `snaf.initialize(..., binding_method='MHCflurry')`. (Genome-wide junction UIDs mean a chr22 subset proves plumbing, not recall.)
- **Verdict: GO** (open-only) / run on Linux or free-GPU. Recipe: [`install/snaf_install.md`](install/snaf_install.md).

---

## ASNEO - GO; open-only path already engineered by #566

- **Repo:** `github.com/bm2-lab/ASNEO` (Aging 2020 authors; `zzb23/ASNEO` is the original mirror). Apache-2.0. Full line-level map in the [#566 experiment README](../issue_566_asneo_crosscheck/README.md).
- **License: Apache-2.0** (LICENSE at repo root). Permissive; covers ASNEO's Python code only, NOT the bundled third-party binaries.
- **NetMHCpan + bundled-binary trap:** `src/software.tar.gz` unpacks `netMHCpan-4.0`, `netCTLpan-1.1`, `pepmatch_db_x86_64` (all non-redistributable, x86, DTU academic EULA). **[#566](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/566) already resolved this:** `apply_optionB_patch.py` (tested, string-anchored) stops ASNEO at its normal-subtracted candidate peptides (`putative_peptide.txt`, `ASNEO.py:252`), repoints bedtools to PATH so `software.tar.gz` is never extracted, and skips `ProcessNeo`. The bundled binaries are fully bypassed. Remaining open-only work = MHCflurry-score `putative_peptide.txt` (the not-yet-written concordance notebook).
- **Reference data:** hg19/GRCh37 FASTA, ASNEO's own `data/iRefSeq.bed`, a normal proteome. Input = STAR `SJ.out.tab` (+ optional BAM); invocation `python ASNEO.py -j <SJ.out.tab> -a <HLA> -g hg19.fa -o <outdir>`.
- **arm64 CPU:** `asneo_env.yml` (python 3.8, biopython 1.79 pinned <1.80, bedtools 2.31.1) solves + import-validates on osx-arm64 (#566 built it locally; the `asneo` conda env exists on this box). The bundled binaries are irrelevant (bypassed). **Real blocker = STAR** for the `SJ.out.tab` (whole-genome index >8 GB; STAR is VM-only per project policy). A chr22-only index can fit 8 GB, but the clean laptop path is to feed a pre-generated `SJ.out.tab` directly to patched ASNEO with `-j`, which runs CPU-only on arm64.
- **Smoke: PASS on arm64.** `install/asneo_smoke.sh` runs the patched caller on ASNEO's own bundled `test/SRR2660032.SJ.out.tab` subset to chr22 (no STAR): all 11 stages end to end, **800 junction-derived candidate peptides** at relaxed thresholds (0 at default - low chr22 coverage). Only the STAR front-end (fresh `SJ.out.tab` from FASTQs) stays VM-bound.
- **Verdict: GO** (open-only solved by #566; caller-proper now smoked locally). Remaining: MHCflurry scoring of `putative_peptide.txt` (#848). Recipe: [`install/asneo_notes.md`](install/asneo_notes.md) + [`install/asneo_smoke.sh`](install/asneo_smoke.sh).

---

## NeoSplice - NO-GO as-is

- **Repo:** `github.com/pirl-unc/NeoSplice` (licensed/maintained fork; `Benjamin-Vincent-Lab/NeoSplice` redirects here). Paper: DOI 10.1093/bioadv/vbac032.
- **Language:** Python **2.7** (EOL) - ~11 CLI scripts glued by an example Nextflow.
- **License: Apache-2.0** code (verbatim LICENSE), **but** README says "NeoSplice is free for academic and non-profit use" and the repo vendors non-redistributable binaries -> the **distributed bundle is not commercially redistributable**.
- **NetMHCpan:** hard code dependency. `kmer_graph_inference.py:80 run_netMHCpan()` shells out with a **required** `--netMHCpan_path`; `SV_summarization.py` parses NetMHCpan `.xls` at fixed nM-affinity column offsets and bins by nM (<=50/150/500). Swap = rewrite the producer call **and** the `.xls` parser, and re-express nM cutoffs as presentation percentiles - a code patch, not a flag.
- **Bundled non-redistributable binaries:** yes - `stage/netMHCpan-4.0-docker/Linux_x86_64/...` and `stage/netMHCIIpan-3.2-docker/Linux_x86_64/...` (COPY'd by the Dockerfile). Their presence in an Apache-2.0 repo is itself a license conflict.
- **Deps:** 2018-pinned C-extension stack (pysam 0.14.1, numpy 1.16.0, scipy 1.2.0), msbwt/msbwt-is (BWT), samtools 1.7, STAR (user supplies BAMs). Bundled `Reference_peptidome/` ~292 MB (GRCh38 8-11mers, regenerable). (Correction to issue: STAR + msbwt, not MOSAIK/bwa.)
- **arm64 CPU:** no path - Py2.7 + 2018 C-ext deps have no osx-arm64 wheels; msbwt-is is Linux C; bundled NetMHC binaries are Linux-x86. Only route is emulated Linux Docker, where `samtools sort -m 15G` + whole-transcriptome BWT strain 8 GB.
- **Verdict: NO-GO** as-is. GO-with-caveat only as a Python-3 fork with an MHCflurry shim replacing `run_netMHCpan()` + its parser, run in emulated Linux.

---

## SINE - NO-GO as-is

- **Repo:** `github.com/GuoLabUCSD/SINE` (Splice Isoform Neoantigen Evaluator; uses OutSplice upstream). Paper: *Int J Mol Sci* 2025, 26(1):205, DOI 10.3390/ijms26010205 (PMC11720059). (The Zotero candidate IS this tool.)
- **Language:** Python 3 (tested 3.10.13) + Bash orchestration.
- **License: GPL-3.0** (LICENSE header "GNU GENERAL PUBLIC LICENSE Version 3"; `gh api .../license` -> `spdx_id: GPL-3.0`). Strong copyleft - constrains any code we copy/adapt. `xls_parse.py` also adapts Rachel Marty's PyPresent.
- **NetMHCpan:** load-bearing for **both** class I and class II. `supplemental/prepare_netmhcpan.py` shells out to `netMHCpan` and `netMHCIIpan`, parses `.xls` percentiles, computes PHBR. **MHCflurry swap is partial only:** MHCflurry can replace class-I but **has no class-II predictor**, and class-II binding is a core selling point - the class-II half cannot be swapped. Needs rewriting the predictor wrapper + `calc_phbr.py`.
- **Bundled non-redistributable binaries:** none in-repo (clean); external tools user-installed.
- **Deps:** pandas 1.5.3, pysam 0.22.0, biopython 1.81, numpy 1.26.0 (all arm64-fine) + bedtools, samtools, **Trinity (via Singularity/Apptainer)**, STAR (upstream), NetMHCpan/NetMHCIIpan 4.1 (academic-gated). Reference: Ensembl CDS FASTA + GTF (r111).
- **arm64 CPU:** poor - Python layer is fine, but Trinity runs in a Linux container (no Docker daemon here) and NetMHCpan ships Linux/Darwin-x86 only.
- **Verdict: NO-GO** as-is (class-II NetMHCpan lock + Trinity-in-container). GO-with-caveat only for a class-I-only, re-plumbed MHCflurry variant on Linux.

---

## NeoHunter - EXCLUDE (operational, not license)

- **Repo:** `github.com/XuegongLab/NeoHunter`. Paper: Ma et al., *Quantitative Biology* 2024, DOI 10.1002/qub2.28 (PMC12806199).
- **Language:** Python scripts orchestrated by `NeoHunter.py` (subprocess-shells to external tools). Splicing branch runs a bundled ASNEO in a `conda run -n ASNEO_env`.
- **License:** **MIT** for NeoHunter's own code - the license detector missed it because the file is misspelled `LISCENCE.md` ("MIT License / Copyright (c) 2023 XuegongLab / ... without restriction ..."); bundled ASNEO is Apache-2.0. **So the issue's "free academic use" label is wrong - both are permissive/redistributable.**
- **Splice relevance:** yes - multi-source (`-at snv,indel,fusion,splicing`), splicing via bundled ASNEO (STAR -> `ASNEO.py -j SJ.out.tab`). A legitimate splice comparator in principle, but its splice signal IS ASNEO.
- **NetMHCpan:** pervasive and hard-coupled - NetMHCpan-4.1 **and** NetMHCstabpan-1.0 required (`config.yaml` paths); `parse_netMHC.py` parses fixed-width stdout in two code paths (NeoHunter core + ASNEO). **NetMHCstabpan (stability) has no MHCflurry equivalent** - that filter can only be dropped, changing ranking.
- **Bundled non-redistributable binaries:** yes - `software/ASNEO/src/netMHCpan-4.0/Linux_x86_64/bin/netMHCpan` (+ data), plus more x86 binaries in `software.tar.gz`.
- **Deps/refs:** GATK 4.2.6.1, MiXCR, STAR 2.7.8a, STAR-Fusion + CTAT `GRCh37_gencode_v19` (tens of GB), VEP r105 cache, OptiType, etc. **hg19-only.**
- **arm64 CPU:** not feasible - vendored x86 binaries, NetMHCpan Linux/x86, STAR >30 GB index build (banned locally), tens-of-GB hg19/CTAT downloads on 8 GB.
- **Verdict: EXCLUDE.** License is fine; the operational lock-in (netMHCpan + NetMHCstabpan + hg19 + STAR + x86) is decisive, and its splice branch is just ASNEO, which we cover directly.

---

## REAL-neo / SPLICE-neo - NO-GO (no released code)

- **Identity:** SPLICE-neo is a module inside Mayo Clinic's in-house **REAL-neo** framework (Wickland et al., *J Immunother Cancer* 2024, DOI 10.1136/jitc-2024-008988, PMID 38754917). One project, not two tools.
- **Repo:** none found - the article ships no GitHub/GitLab/Zenodo/code-availability link (data statement points only to TCGA GDC); `gh search repos/code` returns nothing. **Not distributed as installable software.**
- **License:** no software license exists; the article is CC BY-NC 4.0 (non-commercial), applying to the paper, not code.
- **MHC:** hard-wraps six class-I predictors (NetMHC/NetMHCpan/SMM/SMMPMBEC/PickPocket/MHCflurry) + three class-II (NetMHCII/NetMHCIIpan) via pVACtools - so it depends on the paid NetMHC suite. MHCflurry is already one of its six, but the framework itself is unobtainable, so "swap" is moot.
- **Verdict: NO-GO** (nothing installable). The "R spliceneo/splice-neo" the issue also names resolves to `splice2neo` (audited above).
