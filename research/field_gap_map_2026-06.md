# Field gap map — splice-neoantigen → immunotherapy chain

**Snapshot date:** 2026-06-05 · **Scope:** open gaps across the full cancer-neoantigen → immunotherapy translational chain that a small, non-academic, open-source, independent group can realistically work on.

> **How this was produced.** Two fan-out research passes (web search → source fetch → adversarial verification → synthesis). Pass 1 (8 verified findings, stages 1–2) and Pass 2 (25 findings surviving 3-vote perspective-diverse verification, stages 3–4 + licensing/data-openness). Each gap was checked by three independent skeptics asking: *is it still open (not saturated)?*, *are the data **and** tools truly free+open?*, *is it realistic on one P100 with no wet lab?* Confidence reflects how those votes landed. **Several load-bearing sources are 2024–2026 preprints** — treated as reproduction targets, not settled facts. This is a point-in-time snapshot; the field moves fast.

**Group constraints assumed throughout:** local CPU-only core (M1, 8 GB) + opportunistic free GPU (Kaggle / Colab); forward paid path = RunPod L4/Ada + Cloudflare R2 (see [`docs/migration_runbook.md`](docs/migration_runbook.md)); no clusters, no large-scale pretraining; **no wet lab**; **public + free data and tools only** (open-source code, open weights, free databases — no paywalled / dbGaP-EGA-controlled / commercial-licensed / proprietary-weight resources). Strengths: software engineering, pipelines, public-data reanalysis, benchmarking. Existing niche: splice-junction-derived neoantigen prediction via a Snakemake pipeline (HISAT2/STAR → regtools/`SJ.out.tab` junctions → GENCODE + matched-normal filtering → MHCflurry presentation → TCRdock structure).

> **Compute-posture correction (2026-06-26, [#854](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/854)).** This map was written 2026-06-05 assuming a single NVIDIA P100 (16 GB) + modest CPU VMs on GCP. That GCP/P100 stack was decommissioned on 2026-06-26 to a $0-budget posture; the live compute tier is now the local CPU-only core plus opportunistic free GPU, with RunPod L4/Ada as the funded-revival target (the archived P100 operational recipe lives at [`docs/legacy/gcp_p100_setup.md`](docs/legacy/gcp_p100_setup.md)). Read every "P100" / "single-GPU" / "16 GB" feasibility framing below (methodology note, §2.6, §3.6, §3.8, §3.9, the shortlist) against this current tier. In practice per-gap compute feasibility now **splits into two lanes**: **local-CPU-feasible** (detection, reanalysis, benchmarking, most presentation work - runs on the M1 core) vs **paid-GPU-gated** (AF3-class / Boltz-2 structure work - deferred to an opportunistic free GPU or the funded RunPod path). The Pascal-specific engineering questions (§3.8/§3.9) keep archival value but are now gated on the paid-GPU lane rather than an always-on P100.

**Workability lanes:** **M** = Methods & tools · **B** = Benchmarks & datasets · **R** = Reanalysis & biology · **Rep** = Reproducibility & validation.
**Confidence:** ✓✓ strong (3/3 skeptics held) · ✓ solid · ~ caveated/weak. *Vote splits appear in two forms: `n/3` = n of 3 skeptics held the gap open (supplement pass); `n-m` = n confirmed / m refuted (deep-research pass).*

---

## Contents

- [0. The MHC-predictor licensing wall](#0-the-mhc-predictor-licensing-wall) — the cross-cutting constraint; **read first**
- [1. Stage 1: Detection](#1-stage-1-detection) — 6 entries · your niche, least saturated
- [2. Stage 2: Presentation](#2-stage-2-presentation) — 7 entries
- [3. Stage 3: TCR-pMHC recognition and structure](#3-stage-3-tcr-pmhc-recognition-and-structure) — 10 entries
- [4. Stage 4: Clinical and immunogenicity validation](#4-stage-4-clinical-and-immunogenicity-validation) — 4 entries
- [5. Cross-cutting: what's actually reachable under open-only](#5-cross-cutting-whats-actually-reachable-under-open-only) — 2 entries
- [6. Best bets (ranked)](#6-best-bets-ranked) — the shortlist
- [7. Read-before-you-cite caveats](#7-read-before-you-cite-caveats)
- [8. Method and provenance](#8-method-and-provenance)

---

## 0. The MHC-predictor licensing wall

The single most consequential cross-cutting finding. It is invisible to well-funded academic groups because it doesn't bind them — but it binds you, and that asymmetry is itself the opportunity.

- **NetMHCpan / NetMHCIIpan are BLOCKED for a non-academic group.** Academic-only, non-commercial, binary-only, no-redistribution ("Any use ... which results in any form of commercialization is not allowed"; "If you are not a member of a publicly funded Academic/Education/Research Institution you must obtain a commercial license"). Yet they are the de-facto MHC-binding step in nearly every published splice-neoantigen pipeline — NeoSplice (requires NetMHCpan 4.0), ASNEO (ships it bundled), NeoHunter (requires NetMHCpan 4.1 + NetMHCstabpan), SPLICE-neo/REAL-neo (6-tool consensus dominated by NetMHC variants).
  - Sources: <https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/>, <https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/>
- **No fully-open MHC class-II presentation predictor exists at all.** NetMHCIIpan-4.3 and MixMHC2pred are both non-commercial / binary-only. There is no Apache/MIT/GPL class-II option to drop in.
- **Pep2Vec is a license trap.** Ships an open-*looking* weight binary (`pep2vec.bin`), but `LICENSE.txt` restricts use to publicly-funded academic institutions, prohibits modification and redistribution, and **explicitly forbids publishing benchmark results without Genentech's written consent** — which directly blocks your core activity. BigMHC similarly carries a custom JHU non-commercial license.
  - Sources: <https://github.com/Genentech/Pep2Vec>, <https://github.com/KarchinLab/bigmhc/blob/master/LICENSE>
- **ASNEO illegally bundles the NetMHCpan binary** inside its repo (`src/software.tar.gz`), so its "Apache-2.0" badge misrepresents the effective license of the installed pipeline. Unmaintained since 2020.
  - Source: <https://github.com/bm2-lab/ASNEO>
- **MHCflurry (Apache-2.0, already in your stack) is the *one* fully-open class-I substitute.** <https://github.com/openvax/mhcflurry/blob/master/LICENSE>
- **Detection-side parallel:** **OpenSpliceAI** (GPL-3.0, own freshly-trained OSAI-MANE weights) cleanly replaces the license-encumbered SpliceAI (PolyForm code + CC-BY-NC weights) and Pangolin (NC-provenance weights), and runs whole-genome on a single GPU.
  - Sources: <https://github.com/Kuanhao-Chao/OpenSpliceAI>, <https://elifesciences.org/articles/107454>

**Implication:** you are *structurally* one of the few groups that must — and can therefore credibly — build and benchmark the open-only alternative. That reframes "we can't use the standard tools" into "we own the open re-tooling lane."

---

## 1. Stage 1: Detection

*Your niche; least saturated.*

### 1.1 Intron-retention (IR) neoantigen detection is genuinely unsolved · M, B · Med–Large · ✓✓
**Gap.** No published filter combination reaches FDR < 0.50 for IR, so the leading open tool (splice2neo, TRON/BioNTech) excludes IR entirely. Short-read IR callers show > 47% non-reproducible calls and F1 < 0.26 vs long-read ground truth; 2026 IR-neoantigen work still concedes false positives.
**Public resources.** splice2neo (open R + EasyQuant); IRFinder-S / IRFinder2; long-read truth from SG-NEx / HCC1395 (see §5).
**Open/free status.** Fully open; the blocker is scientific difficulty, not access.
**Evidence.** Lang et al. 2024, *Bioinformatics Advances* vbae080 ("For IR events no filter combination led to an estimated FDR lower than 0.50"; IR excluded). Broseus et al. 2022, *Genome Biology* (F1 < 0.26 across 8 tools; 47.7% single-tool calls; Fleiss κ = 0.113).
**Sources.** <https://academic.oup.com/bioinformaticsadvances/article/4/1/vbae080/7684965> · <https://github.com/TRON-Bioinformatics/splice2neo> · <https://pmc.ncbi.nlm.nih.gov/articles/PMC9652823/>

### 1.2 No single splicing modality is complete — combined open detector missing · M, B · Med–Large · ✓
**Gap.** DNA-splice-site-mutation and de-novo-RNA-splicing detection capture largely non-overlapping event sets: across TCGA BRCA/LUAD/LUSC/LIHC, only 3,501 single-exon-skipping events are shared vs 5,165 (28.7%) DNA-only and 9,371 (52%) RNA-only (~19% overlap). SPLICE-neo uniquely handles >2-exon skipping (26–33% of skipping events); splice2neo integrates SpliceAI/MMSplice/Pangolin with LeafCutter/SplAdder/RegTools/IRFinder — but no head-to-head modality integration + evaluation exists.
**Open/free status.** Open. **Caveat:** the overlap figures are splicing *events*, not validated neoantigens (this is why verification split 2-1).
**Evidence.** Wickland et al. 2024 (SPLICE-neo, *J Immunother Cancer*, PMC11097882): "essential to profile aberrant RNA splicing events from both DNA and RNA mutational sources."
**Sources.** <https://pmc.ncbi.nlm.nih.gov/articles/PMC11097882/> · <https://academic.oup.com/bioinformaticsadvances/article/4/1/vbae080/7684965>

### 1.3 No head-to-head benchmark of the open splice-neoantigen callers · B, Rep · Medium · ✓✓
**Gap.** SNAF, NeoHunter (wraps ASNEO on STAR `SJ.out.tab`), SPLICE-neo/REAL-neo, splice2neo, SINE, NeoSplice all exist independently with no common dataset, ground truth, or metrics. Building the first public benchmark + reproducibility harness is squarely your software/benchmarking strength.
**Public resources.** SNAF (PyPI/GitHub/Zenodo; GTEx+TCGA normal filter); NeoHunter (XuegongLab); ASNEO (Apache-2.0 *code*, bm2-lab); public TCGA RNA-seq via recount3.
**Open/free status.** Mostly open. **Caveats:** NeoHunter is stated "free academic use" (may not be OSI-open); SNAF's full-cohort-scale claim split 1-2 in verification; ASNEO bundles NetMHCpan (see §0).
**Sources.** <https://pmc.ncbi.nlm.nih.gov/articles/PMC11517820/> · <https://pmc.ncbi.nlm.nih.gov/articles/PMC12806199/> · <https://github.com/XuegongLab/NeoHunter> · <https://github.com/bm2-lab/ASNEO>

### 1.4 Short-read IR-calling unreliability is unreproduced on cancer lines · Rep, B · Low–Med · ✓✓
**Gap.** Broseus et al. 2022 (none of 8 tools F1 > 0.26 vs matched long reads) was done on only two **non-cancer** specimens (whole blood, iPSC). It has not been reproduced or extended to cancer cell lines — yet IR is a claimed neoantigen source. A load-bearing silent-failure mode for any short-read splice pipeline, including yours.
**Public resources.** Original SRA: SRP065930 (HX1 blood, PacBio Iso-Seq + short read), SRP098984 (iPSC); extend on SG-NEx cancer lines; tools IRFinder-S, superintronic, iREAD, KMA, IntEREst + rMATS/MAJIQ/SUPPA2 (all FOSS); original scripts MIT on GitHub/Zenodo.
**Open/free status.** Fully open. No controlled-access dependency.
**Evidence.** Broseus et al., *Genome Biology* 2022 (PMC9652823 / s13059-022-02789-6).
**Sources.** <https://pmc.ncbi.nlm.nih.gov/articles/PMC9652823/> · <https://github.com/GoekeLab/sg-nex-data> · <https://ascopubs.org/doi/10.1200/CCI.21.00124>

### 1.5 (Enabler) An open long-read IR/junction ground-truth benchmark is buildable from SG-NEx + HCC1395 · B, Rep · Medium · ✓
**Gap/enabler.** No standardized open IR/splice-junction ground-truth benchmark is tied to **cancer** lines and pegged to long-read truth. SG-NEx and SEQC2/HCC1395 each ship matched long-read (ONT direct-RNA/cDNA) **and** short-read RNA-seq for the same cancer specimens — exactly the substrate to score a short-read junction/IR caller. Your own HISAT2/regtools and STAR `SJ.out.tab` calls have never been scored against long-read truth on cancer lines.
**Public resources.** SG-NEx (GoekeLab; 5+ cancer lines HCT116/HepG2/A549/MCF7/K562; AWS Open Data `s3://sg-nex-data`); SEQC2/HCC1395 + HCC1395BL (FDA reference); LongBench; ENCODE4 long-read (264 PacBio libs); LRGASP. Tools: IsoQuant, FLAIR, IRFinder-S, superintronic, regtools, STAR.
**Open/free status.** Fully open (AWS Open Data, no controlled access). **Caveat:** the richest *patient* long-read resource (Long-Read POG, 189 tumors) is EGA-controlled — so a patient-level (vs cell-line) benchmark is *not* open.
**Evidence.** Chen et al., *Nature Methods* 2025 (s41592-025-02623-4, SG-NEx).
**Sources.** <https://github.com/GoekeLab/sg-nex-data> · <https://www.nature.com/articles/s41592-025-02623-4> · <https://www.nature.com/articles/s41592-024-02298-3>

### 1.6 (Enabler) OpenSpliceAI adds a fully-open de-novo / variant splice-scoring lane · M, Rep, B · Medium · ✓
**Gap.** Your pipeline detects junctions empirically (HISAT2/STAR) but has no in-silico splice-strength scoring lane to corroborate/prioritize junctions or predict splice-disrupting variant effects. OpenSpliceAI (GPL-3.0, own OSAI-MANE weights) + splice2neo (MIT) is a fully-open path to add one — and an OpenSpliceAI-vs-SpliceAI(CC-BY-NC weights) agreement benchmark is a publishable open-vs-encumbered reproducibility result.
**Sources.** <https://github.com/Kuanhao-Chao/OpenSpliceAI> · <https://elifesciences.org/articles/107454> · <https://github.com/TRON-Bioinformatics/splice2neo>

---

## 2. Stage 2: Presentation

### 2.1 The NetMHC licensing wall + no open class-II predictor (benchmark gap) · B, Rep · Medium · ✓✓
**Gap.** Quantify what published splice-neoantigen pipelines lose when NetMHCpan is swapped for MHCflurry (class-I), and document the complete absence of a fully-open class-II presentation predictor. See §0 for the full licensing picture.
**Public resources.** MHCflurry 2.x (Apache-2.0, open weights, already in stack); IEDB + HLA Ligand Atlas eluted-ligand data for evaluation; open SNAF/NeoSplice/ASNEO peptide outputs as input candidate sets.
**Open/free status.** The *evaluation* is fully open; the tools being evaluated are the constraint — which is the finding.
**Sources.** <https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/> · <https://github.com/openvax/mhcflurry> · <https://pmc.ncbi.nlm.nih.gov/articles/PMC11097882/>

### 2.2 MHC-II predictors biased toward DR over data-scarce DP/DQ · M, B · Medium · ✓✓
**Gap.** Class-II tools predict DR better than DP/DQ because DR data dominates training sets; still open as of Dec 2025. An independent IEDB-derived benchmark (~67,061 binders / ~67,163 non-binders across 20 allotypes — 9 DR, 1 DQ, 10 DP) already exists to extend on the under-served DP/DQ side. NetMHCIIpan-4.3 is the open *baseline* (free web server + downloadable training data; 675,364 EL positives, 142 class-II molecules) — but its standalone *code* is academic-license-restricted.
**Evidence.** Faisal et al. 2024, *Front Immunol* (11-tool MHC-II benchmark); Nilsson et al. (NetMHCIIpan-4.3, *Sci Adv* adj6367).
**Sources.** <https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1293706/full> · <https://www.science.org/doi/10.1126/sciadv.adj6367>

### 2.3 No systematic reanalysis of public HLA-ligandome MS for splice neoepitopes · R, B · Medium · ✓✓
**Gap.** Every MS confirmation of splice neoantigens so far used lab-internal MS or a few bundled cohorts. Nobody has run the inverse, purely-computational reanalysis: take a fixed predicted splice-neoepitope panel and ask whether those exact sequences appear in the already-search-completed peptide tables of HLA Ligand Atlas, SysteMHC v2.0, and caAtlas. **IEAtlas reanalyzed 2,747 public MS samples for non-canonical epitopes but covered only lncRNA / pseudogene / UTR / uORF-dORF — *not* splice junctions**, leaving the splice category a named hole in an otherwise saturated landscape.
**Public resources.** HLA Ligand Atlas (PRIDE PXD020186/PXD024837); SysteMHC Atlas v2.0 (~1.0M class-I + 1.1M class-II peptides); caAtlas (364,283 peptides / 43 cancer datasets); predicted panels from SNAF/SpliceMutr/NeoSplice or your own pipeline.
**Open/free status.** Fully open. **Caveat:** rigorous spectral re-confirmation benefits from FragPipe/MSFragger (free non-commercial only) — Comet is the fully-open fallback.
**Evidence.** IEAtlas, *Nucleic Acids Res* 2023 (D409); NeoSplice (PMC9154024); SNAF (PMC11517820).
**Sources.** <https://academic.oup.com/nar/article/51/D1/D409/6696849> · <https://hla-ligand-atlas.org/search> · <https://academic.oup.com/nar/article/52/D1/D1062/7449490> · <https://github.com/frankligy/SNAF>

### 2.4 SysteMHC v2.0's 78,959 unclassified non-UniProt peptides are an unmined splice substrate · R, B · Low · ✓
**Gap.** SysteMHC Atlas v2.0 reports 78,959 non-UniProt peptides it explicitly does **not** classify by origin ("an exhaustive classification ... will be made available in a future version"). Nobody has asked how many map to tumor splice junctions vs SNV/fusion/ncORF/proteasomal-cis-splicing. Your pipeline already produces a GENCODE+normal-filtered splice-junction peptide reference to intersect against them — a download + sequence-mapping task, no MS re-search.
**Public resources.** SysteMHC v2.0 non-UniProt set + per-sample FASTA; GENCODE; your translator; SPIsnake (to subtract the proteasomal-splicing confounder).
**Evidence.** SysteMHC Atlas v2.0, *Nucleic Acids Res* 2024 (D1062 / PMC10767952).
**Sources.** <https://academic.oup.com/nar/article/52/D1/D1062/7449490> · <https://pmc.ncbi.nlm.nih.gov/articles/PMC10767952/> · <https://github.com/QuantSysBio/SPIsnake>

### 2.5 Public immunopeptidomics MS contamination — no open contaminant tool/benchmark · B, M, Rep · Medium · ~
**Gap.** Two coupled gaps. (1) **Generalizability:** Pep2Vec's headline ~5% contaminant rate is from one corpus; earlier work found 12–40% tryptic contaminants (>42% for HLA-B*57:01/B*35:01) — highly allele- and dataset-dependent, so a single "5%" hides large variance and is unvalidated across independent PRIDE/MassIVE sets. (2) **No open tool:** Pep2Vec is binary-only + academic-only + no-benchmark-publish (§0), so no open, reproducible contaminant-detection model or standardized benchmark corpus exists. Building an open baseline (tryptic-terminus + length + binding-motif + retention-time heuristics) + a labeled benchmark from public PRIDE sets is a 1-2 person methods+benchmark task.
**Open/free status.** Substrate open (PRIDE/MassIVE, MhcVizPipe, IEDB); the SOTA tool (Pep2Vec) is the hard exclusion that *makes* this a gap. **Caveat:** the 5% figure is single-source / self-reported / non-peer-reviewed — a hypothesis to test, not a field rate.
**Sources.** <https://www.biorxiv.org/content/10.1101/2024.10.14.618255v1> · <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8977642/> · <https://pmc.ncbi.nlm.nih.gov/articles/PMC8717601/>

### 2.6 No open, splice-junction-specific FDR-control for MS validation · Rep, M · Medium · ~
**Gap.** Searching spectra against a custom junction-peptide database inflates the search space and over-estimates FDR at fixed thresholds, yet there's no published splice-junction-specific FDR-control protocol quantifying how many "validated" splice neoepitopes survive properly-stratified target-decoy vs naive search. The general machinery (Sequoia + SPIsnake, stratified target-decoy) is open but HPC/Slurm-oriented and built for genomic-ORF / proteasomal-cis-spliced peptides, not RNA splice junctions at 1-2-person/P100 scale.
**Open/free status.** Substrate open (nf-core/MHCquant MIT, Comet open, PRIDE/MassIVE raw, SPIsnake). **Caveat:** MSFragger engine is free non-commercial only (Comet is the open fallback); free-vote split 1/3 in verification on tooling openness.
**Critical terminology guard.** This is **RNA splice-junction-derived** (cis at the transcript level), NOT **proteasomal cis-spliced** peptides — the "Are there indeed spliced peptides?" / Neo-Fusion debate (13–45% → ~1–6%) is about the proteasomal kind and must not be conflated.
**Sources.** <https://pmc.ncbi.nlm.nih.gov/articles/PMC12397870/> · <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8724635/> · <https://github.com/nf-core/mhcquant> · <https://github.com/QuantSysBio/SPIsnake>

### 2.7 Recent presentation SOTA are reproduction targets (but Pep2Vec is blocked) · B, Rep · Medium · ✓
**Gap.** A 2025 benchmark of 17 peptide-HLA-I tools (>290,000 peptides / 44 HLA-I alleles) documents persistent tool-performance variability, limited interpretability, and dataset-quality bias. Pep2Vec claims higher average precision than MHCflurry/NetMHCpan — but it's a license trap (§0). The publishable open move is to benchmark the *fully-open* models (MHCflurry) on splice neoepitopes, explicitly excluding Pep2Vec/BigMHC/NetMHC so results are publishable.
**Caveat.** Both are non-peer-reviewed preprints with self-reported SOTA. Three pass-1 "best-tool" ranking claims were *refuted* in verification (see §7) — do not assert a single best predictor.
**Sources.** <https://www.biorxiv.org/content/10.1101/2025.04.10.648169v1.full> · <https://www.biorxiv.org/content/10.1101/2024.10.14.618255.full.pdf>

---

## 3. Stage 3: TCR-pMHC recognition and structure

### 3.1 No splice-neoepitope slice in any TCR-immunogenicity generalization benchmark · B, Rep · Medium · ✓✓
**Gap.** Every standardized 2024–2026 TCR-pMHC generalization benchmark (ePytope-TCR, IMMREP23, the 2026 arXiv generalization study) is dominated by viral epitopes; splice-junction neoepitopes are absent as an evaluation category. ePytope-TCR shows prediction accuracy correlates with epitope prevalence (Pearson r = 0.69) and all methods drop to ~random for rare epitopes. Nobody has assembled a "splice-derived neoepitope held-out" slice or measured how badly NetTCR-2.2/MixTCRpred/pMTnet do on splice peptides — and you already generate the candidate peptides.
**Public resources.** ePytope-TCR (GitHub+Zenodo, 21 models wrapped); IMMREP23 (GitHub/Kaggle); predictors NetTCR-2.2, MixTCRpred, pMTnet, epiTCR, PanPep, ERGO-II, SCEPTR (all open). Splice ground truth in open-access supplements: SNAF (5 IFN-γ-validated), Kwok et al. *Nature* 2024 (8 validated TCRs, HLA-A*02:01).
**Open/free status.** Tools+harnesses fully open. **Caveat:** splice-TCR ground truth is scattered in supplements (~13–20 triples total) — curation effort, and statistically thin (frame as a *stress-test probe*, not a powered benchmark).
**Sources.** <https://pmc.ncbi.nlm.nih.gov/articles/PMC12366652/> · <https://arxiv.org/html/2606.04994> · <https://pmc.ncbi.nlm.nih.gov/articles/PMC11903331/>

### 3.2 No audit of the splice/non-canonical fraction of TCR-specificity databases · R, B · **Low** (days) · ✓✓
**Gap.** "TCR databases are viral-skewed" is asserted but never quantified for splice. What fraction of VDJdb + IEDB + McPAS-TCR antigen entries are (a) cancer neoantigens at all, (b) specifically splice-junction-derived / non-canonical-ORF / frameshift vs point-mutation? This determines whether any supervised TCR model *could* have learned splice features, and gives you a citable "your recognition step is out-of-distribution" denominator. Pure public-data exercise.
**Public resources.** VDJdb, IEDB T-cell assay export, McPAS-TCR (all free TSV); a unified 190,670-TCR / 2,313-epitope DB already exists (caokai1073 GitHub); classify against SNAF/neojunction catalogs; GENCODE for canonical subtraction.
**Open/free status.** Fully open. No controlled-access component.
**Sources.** <https://pmc.ncbi.nlm.nih.gov/articles/PMC12366652/> · <https://github.com/caokai1073/Papers-for-TCR-antigen-prediction>

### 3.3 CDR3-pLDDT-as-affinity-signal is AF3-only and unreproduced on open weights · Rep, R, M · Medium · ✓✓
**Gap.** A Nov 2025 benchmark claims TCR CDR3-region pLDDT tracks mutation-induced binding-affinity direction in >80% of cases and a CDR3-pLDDT reranking lifts AF3 Top-1 by 2.9–4.3% (beating ipTM/ranking-score). But it's single-study, AF3-specific (AF3 weights are academic-license + request-gated), so the field's most actionable structure-confidence signal is **unreproduced on the open-weight models you can actually run** (Boltz-2, TCRdock, TCRmodel2) — and never tested against actual T-cell immunogenicity (IEDB) rather than in-vitro affinity. A focused, fully-open reproduction + reanalysis that plugs straight into your TCRdock stage.
**Open/free status.** Replication side fully open (Boltz-2 MIT, TCRdock open, data free). The original model's gating is *why* the open reproduction matters.
**Sources.** [Benchmarking TCR-pMHC structure prediction — a unified evaluation and CDR3-based functional insights (biorxiv 2025.11.30.691400)](https://www.biorxiv.org/content/10.64898/2025.11.30.691400v1.full) · <https://github.com/phbradley/TCRdock> · <https://github.com/jwohlwend/boltz>
> Note: bioRxiv DOIs minted from late 2025 use the `10.64898` prefix (not the legacy `10.1101`); both `10.64898` URLs in §3.3/§3.4 were verified to resolve 2026-06-05.

### 3.4 Structure-confidence (ipTM/PAE) isn't calibrated to TCR binding — no open decoy-controlled TCR-pMHC benchmark · Rep, B, R · Medium · ✓
**Gap.** Multiple 2025–2026 papers show internal confidence scores generate geometrically plausible TCR-pMHC complexes but fail to discriminate true binders from mispaired non-binders. The rigorous thousands-of-decoys controlled study was done for nanobody-antigen (106 cognate + 11,342 shuffled across AF3/Boltz-2/Chai-1) but **not** for TCR-pMHC. A direct port — cognate STCRDab/TCR3d structures vs systematically shuffled TCR↔pMHC pairs, scoring whether Boltz-2/TCRdock/Chai-1 ipTM ranks cognate highest, and whether *any* structural confidence predicts IEDB-measured immunogenicity beyond MHCflurry presentation — is wide open.
**Public resources.** STCRDab (weekly PDB), TCR3d 2.0 (free, includes affinities + crossing/incident angles); IEDB/VDJdb/McPAS/10x for labels + negatives; Boltz-1/2 (MIT), TCRdock (open) as scorers; STAG (open negative-set template).
**Open/free status.** Fully open. **Caveat:** AF3 itself is academic-license + request-gated — benchmark the Boltz-2/Chai-1/TCRdock triad and cite AF3 as reference.
**Sources.** [Structural Plausibility Without Binding Specificity — limits of AI-based antibody-antigen confidence scores (biorxiv 2026.03.02.709004)](https://www.biorxiv.org/content/10.64898/2026.03.02.709004v1) · <https://tcr3d.ibbr.umd.edu/> · <https://github.com/jwohlwend/boltz>

### 3.5 Your pipeline (and SNAF) lack a TCR-repertoire recognition step · M · Medium · ✓
**Gap.** Splice-discovery pipelines stop at MHC presentation + a peptide-intrinsic immunogenicity score (DeepImmuno). None ask "does this patient's TCR repertoire contain a clonotype likely to recognize this splice pMHC?" SNAF authors explicitly flag the missing TCR-level model. A lightweight add-on — run open pan-epitope predictors over the patient repertoire × candidate splice pMHC, emit a repertoire-coverage score — is shippable, *if* framed as a calibrated, OOD-flagged add-on (not a confident binder caller; see §3.1/§3.6 caveats).
**Public resources.** MixTCRpred, SCEPTR (CPU-runnable), NetTCR-2.2, PanPep, ERGO-II (open); repertoire extraction via **TRUST4 (MIT)** — prefer over MiXCR (free-academic only) to stay genuinely open for a non-academic group; ImmuneCODE/10x reference repertoires.
**Sources.** <https://pmc.ncbi.nlm.nih.gov/articles/PMC11517820/> · <https://pmc.ncbi.nlm.nih.gov/articles/PMC12366652/>

### 3.6 Decoy/negative-sampling sensitivity untested on the cancer/splice setting · Rep, M · Medium · ✓
**Gap.** Reported TCR-immunogenicity AUCs swing wildly with how negatives are constructed (within-pair shuffle vs background repertoire vs hard negatives), but published stress-tests are on viral benchmarks. Re-scoring the same splice-neoepitope positives under multiple standardized negative-sampling schemes and reporting the AUC swing would expose silent failure modes before you trust any score.
**Public resources.** Pitfalls protocols (bioRxiv 2023.04.06.535863 + reply); ImmuneCODE (free, Azure Open Dataset) + 10x donor TCRs as decoy background; predictors as above (small models, P100/CPU-fine).
**Sources.** <https://www.biorxiv.org/content/10.1101/2023.04.06.535863.full.pdf> · <https://learn.microsoft.com/en-us/azure/open-datasets/dataset-immunecode>

### 3.7 Class-II TCR-pMHC structure prediction is the universal weak spot · B, Rep · Low–Med · ✓
**Gap.** Every 2025–2026 benchmark is class-I-dominated and names class II as the unsolved residual but never isolates it (the JCIM study had 6 class-II vs 21 class-I; original TCRdock benchmark was 8 class-I only). A dedicated class-II-only open benchmark — all class-II TCR-pMHC in STCRDab/TCR3d, temporally held out past each tool's training cutoff, scored across TCRdock/Boltz-2/TCRmodel2 — doesn't exist. Relevant if the pipeline extends to CD4+ neoepitopes.
**Open/free status.** Fully open on TCRdock+TCRmodel2+Boltz-2. **Caveat:** tFold-TCR is PolyForm-Noncommercial (usable for non-commercial benchmarking, not redistributable into a product); AF3 academic-only.
**Sources.** <https://pubmed.ncbi.nlm.nih.gov/40512975/> · <https://academic.oup.com/nar/article/53/D1/D604/7779351> · <https://pmc.ncbi.nlm.nih.gov/articles/PMC9859041/>

### 3.8 (Hardware moat) Single-P100/16 GB feasibility of AF3-class TCR-pMHC models is unverified · Rep, M · Low–Med · ✓
> **Posture note ([#854](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/854)):** the P100 is gone (see the compute-posture correction up top). This "does it fit in 16 GB" question is now a **paid-GPU-lane** item - reframe it against the opportunistic free-GPU tier (Kaggle / Colab) or the funded RunPod L4/Ada target; the Pascal-specific cu126 workarounds keep archival value ([`docs/legacy/gcp_p100_setup.md`](docs/legacy/gcp_p100_setup.md)).
**Gap.** AF3-class models (AF3, Boltz-1/2, Chai-1) OOM on 16 GB at default settings; a full TCR-pMHC complex (~700–900 residues) is near the edge. LMI4Boltz (Oct 2025, MIT) claims +100%/+67% token limits via in-place updates + chunking — but tested only on 24 GB, **never on Pascal/P100**. A documented "does Boltz-2 + LMI4Boltz run on a single 16 GB Pascal GPU, at what complex-size ceiling, with which kernel/precision workarounds" recipe is a pure engineering+reproducibility win nobody with A100s bothers to write — and your existing cu126/Pascal expertise (per CLAUDE.md) transfers directly.
**Open/free status.** Fully open (LMI4Boltz MIT, Boltz-2 MIT). Risk is technical (Pascal arch), not legal.
**Sources.** <https://github.com/tlitfin/lmi4boltz> · <https://www.biorxiv.org/content/10.1101/2025.10.29.684571v1> · <https://github.com/jwohlwend/boltz>

### 3.9 Boltz-2 runnability audit + open TCR-pMHC structure-scoring lane · Rep, B, M · Low–Med · ✓
> **Posture note ([#854](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/854)):** the "two single-P100 constraints" below are now **paid-GPU-lane** questions - the 16 GB co-fit ceiling and Pascal SM 6.0 wheel issues apply to whichever GPU the revival actually uses (opportunistic free tier or RunPod L4/Ada), not an always-on P100.
**Gap.** Boltz-2 (MIT code **and** weights — uniquely unencumbered vs AF3) is a candidate structure-scoring step alongside TCRdock, but two single-P100 constraints are under-documented: (1) structure ~11 GB + affinity 7–8 GB likely won't co-fit on 16 GB (needs splitting/offload); (2) Pascal SM 6.0 wheel compatibility (your documented cu126 pin issue). A reproducible runnability note + a Boltz-2-vs-TCRdock structural-confidence benchmark on your splice neoepitopes — both fully open.
**Sources.** <https://github.com/jwohlwend/boltz> · <https://boltz.bio/boltz2> · <https://github.com/phbradley/TCRdock>

### 3.10 TCRdock register/orientation generalization test · Rep, B · Low–Med · ~ (weak)
**Gap.** TCRdock relies on 12 hand-clustered docking templates + AF2 monomer params and documents failures on divergent binding modes (reversed orientation, displaced footprints) — exactly the out-of-distribution case splice neoepitopes represent. Whether newer template-free open models (Boltz-2, tFold-TCR) escape this ceiling on unusual registers is untested. **Verification flagged this weak (0/3 still-open votes)** — skeptics felt template/orientation generalization is being addressed by the newer models already; treat as low-priority / fold into §3.4.
**Sources.** <https://pmc.ncbi.nlm.nih.gov/articles/PMC9859041/> · <https://github.com/TencentAI4S/tfold>

---

## 4. Stage 4: Clinical and immunogenicity validation

### 4.1 No public benchmark links *predicted* splice-neoantigens to *measured* T-cell immunogenicity · B, M, Rep · Medium · ✓✓
**Gap.** No standardized benchmark pairs computationally-predicted splice neoantigens with experimentally-measured T-cell immunogenicity, the way TESLA does for SNV/indel. **TESLA explicitly excluded splice isoforms**, so its 608 validated peptides can't ground-truth a splice pipeline. Every splice tool ran its own one-off validation on a handful of peptides — no shared harness, no negative set, no held-out split. Assemble the first open harness: pull splice/aberrant-junction/frameshift-from-splicing entries from CEDAR + IEDB + NEPdb + dbPepNeo2.0, reconstruct each peptide's junction, label positive/negative, publish a fixed FASTA+label benchmark + scoring script. Pure curation+software.
**Public resources.** CEDAR (289,028 cancer epitopes / 87,903 T-cell assays); NEPdb (>17,000 validated); dbPepNeo2.0 (801 HC neoantigens + 630 reactive TCRβ); TESLA set as contrast.
**Open/free status.** Fully open. **Caveat:** splice entries are sparse in these DBs (mostly SNV/MS-ligandome) — the benchmark will be small, but that scarcity *is* the publishable point.
**Evidence.** TESLA: Wells et al., *Cell* 2020 (splice isoforms excluded; 608 peptides). SNAF validated 5; Kwok et al. *Nature* 2024 confirmed 2 (GNAS, RPL22).
**Sources.** <https://www.cell.com/cell/fulltext/S0092-8674(20)31156-9> · <https://academic.oup.com/nar/article/51/D1/D845/6761984> · <https://pmc.ncbi.nlm.nih.gov/articles/PMC8078594/> · <https://pmc.ncbi.nlm.nih.gov/articles/PMC9043652/>

### 4.2 No open meta-analysis reconciles the contradictory splice-burden vs ICI/survival associations · R, Rep, B · Medium · ✓✓
**Gap.** The clinical claim "splice-neoantigen burden predicts ICI response/survival" is genuinely unsettled and directionally contradictory, with no harmonized open reanalysis across cohorts using one pipeline. SNAF: high splice-neoantigen burden trended toward *poor* OS in untreated TCGA-SKCM yet *improved* OS under ICB in Van Allen (underpowered). IR-neoantigen load is *favorable* in PDAC but *unfavorable* in multiple myeloma. SpliceMutr: splicing antigenicity higher in ICI responders (CM-038/Riaz, P = 7e-5). SINE: immunogenic splice events lost in responders post-nivolumab — but each with a different tool, normal filter, and burden definition. Run **one** pipeline uniformly across open-tier ICI cohorts and publish the first apples-to-apples meta-analysis.
**Public resources (openness audited).** **Open (full FASTQ):** Hugo (GSE78220), Riaz/CM-038 (GSE91061), Gide (ENA PRJEB23709), Auslander (GEO GSE115821). **Processed-only:** IMvigor210 (raw reads EGA-controlled; processed counts free via R package). **Controlled (exclude or junction-summary only):** Van Allen, Liu, raw TCGA (dbGaP) — but TCGA+GTEx **junction counts** are free via recount3/Snaptron.
**Sources.** <https://pmc.ncbi.nlm.nih.gov/articles/PMC11517820/> · <https://pmc.ncbi.nlm.nih.gov/articles/PMC11648103/> · <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11720059/>

### 4.3 Splice-burden definitions are non-reproducible across tools on identical samples · Rep, B, R · Low–Med · ✓
**Gap.** Before any clinical-association claim can be trusted, the burden *metric* must be reproducible — and it almost certainly isn't. SNAF, SpliceMutr, SINE, NeoSplice, splice2neo each define "splicing-derived neoantigen burden" differently (per-junction MHC-bound count vs cohort-level differential-usage score vs lost-event count), use different normal filters (matched vs within-cohort vs GTEx ~2,500-sample panel vs none) and different junction callers. No study has run ≥2 on the **same** RNA-seq and measured concordance of burden values or patient rankings. If two tools disagree on who is "high burden," every downstream survival curve is suspect.
**Open/free status.** All tools GitHub-open. **Caveat:** SpliceMutr notes "non-academic use needs a license" — anchor on SNAF/SINE (unrestricted), treat SpliceMutr as literature reference. Test cohorts (Hugo/Riaz) GEO-open; recount3 junction layer open.
**Sources.** <https://github.com/frankligy/SNAF> · <https://github.com/FertigLab/splicemutr> · <https://github.com/GuoLabUCSD/SINE>

### 4.4 (Stage 1→4 reanalysis) Splice neoantigens are shared/recurrent/public — validatable without a wet lab · R · Medium · ✓ (efficacy debated)
**Gap/opportunity.** Splice neoantigens are substantially shared and recurrent — up to 90% of melanoma patients; 940 T-cell neoantigens predicted in >15% of TCGA patients; SLC45A2 in 212/472 TCGA-SKCM; recurrent GNAS/RPL22 neojunctions with isolated recognizing TCRs (GNAS tumor-wide across glioma/mesothelioma/prostate/liver). You can cross-match predicted candidates against IEDB positive T-cell assays, VDJdb/10x public TCR data, and PRIDE/MassIVE immunopeptidomics — pure reanalysis.
**Caveat.** The *existence* of public recurrent neojunctions is solid (3-0); downstream **therapeutic efficacy** is genuinely debated (some cohorts null; one functional screen found 2/4 candidates reactive). The "spatial conservation of GNAS" claim split 2-1; some multi-region cohorts are partly controlled-access.
**Sources.** <https://pmc.ncbi.nlm.nih.gov/articles/PMC11517820/> · <https://www.nature.com/articles/s41586-024-08552-0>

---

## 5. Cross-cutting: what's actually reachable under open-only

### 5.1 Recurrence/conservation reanalysis is open-tier executable ONLY via junction summaries (recount3/Snaptron), not raw TCGA · R, Rep, B · Medium · ~
The landmark recurrent-neojunction paper (Kahles/tumor-wide, *Nature* 2025) used raw TCGA FASTQ (~4,900 samples) + 9,166 GTEx FASTQ — both **dbGaP-controlled** (phs000178), so that exact pipeline is *not* reproducible from raw data by an open-only group; GDC harmonized data does not expose per-sample splice-junction-quantification files. **The open path is summary-level:** recount3/Snaptron redistribute precomputed junction counts (tcgav2 ~32M junctions / ~11K TCGA samples; gtexv2 ~33M / ~19K GTEx) — de-identified aggregate counts, free, no login. So the recurrence + cancer-specificity *call* is open-tier executable, but **only at junction-summary granularity** — per-read evidence, allele-specific phasing, and patient HLA genotype (for full neoepitope calling) stay cell-line-only or controlled. *(Verification split 1/3 on "still open" — skeptics felt this is more a feasibility ceiling than a novel gap; framed accordingly.)*
**Sources.** <https://snaptron.cs.jhu.edu/data.html> · <https://link.springer.com/article/10.1186/s13059-021-02533-6> · <https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/467>

### 5.2 A standardized open GTEx pan-normal junction-frequency filter panel is assemblable but doesn't exist · B, M, R · Low–Med · ✓
Your pipeline filters tumor junctions against a single matched normal (or labels everything `tumor_exclusive`). The field standard filters candidate neojunctions against a **large pan-normal panel (GTEx, expressed in <1% of normals)** for cancer-specificity. No open, splice-neoepitope-ready precomputed pan-normal junction-frequency panel exists — every group rebuilds it. Snaptron gtexv2 makes a junction→normal-tissue-prevalence table assemblable openly and packageable as a reusable resource (fits the `resources/` convention). This gates the move from per-patient matched-normal filtering to population-level cancer-specificity.
**Sources.** <https://snaptron.cs.jhu.edu/data.html> · <https://pmc.ncbi.nlm.nih.gov/articles/PMC11903331/>

---

## 6. Best bets (ranked)

Where open × unsaturated × niche-aligned × P100-tractable all line up. (§refs point to the detailed entries above.)

1. **Flagship — first open splice-neoantigen caller benchmark + burden-concordance study.** §1.3 + §1.5 + §4.3. SNAF / NeoHunter / SPLICE-neo / splice2neo / SINE head-to-head on shared public data — SG-NEx cancer lines for long-read detection truth, recount3 for cohort-scale burden concordance. Methods+Benchmarks+Reproducibility; your exact strength; nobody owns it. **✓✓**
2. **Predicted-splice → measured-immunogenicity benchmark.** §4.1. First open FASTA+label harness from CEDAR/IEDB/NEPdb/dbPepNeo2.0; TESLA left splice out. No GPU; scarcity of ground truth is the point. **✓✓**
3. **Splice peptides in the public immunopeptidome.** §2.3 + §2.4. Intersect your junction-peptide translator against HLA Ligand Atlas / SysteMHC v2.0 (the 78,959 unclassified) / caAtlas — the one non-canonical class IEAtlas skipped. **✓✓**
4. **Turn the licensing wall into a paper.** §0 + §2.1. Benchmark MHCflurry vs the blocked NetMHC stack on splice neoepitopes; document the class-II void; flag the ASNEO binary landmine. Differentiated to open groups; de-risks your own stack. **✓✓**
5. **Harmonized meta-analysis: splice-burden vs ICI response.** §4.2. One pipeline across open ICI cohorts (Hugo/Riaz/Gide) to sharpen a contradictory literature. **✓✓**
6. **TCR DB-composition audit → splice-recognition layer.** §3.2 → §3.1/§3.5. The audit (~days) is the low-risk warm-up that motivates building the missing TCR-recognition step (TRUST4 + open pan-epitope predictor). **✓✓ → ✓**
7. **TCR-pMHC structure-confidence calibration on open weights.** §3.3 + §3.4. Reproduce the AF3-only CDR3-pLDDT claim on Boltz-2/TCRdock + test immunogenicity signal over MHCflurry; plugs into your TCRdock stage. **✓✓**
8. **The hard, high-payoff one — intron-retention neoepitope detection.** §1.1 + §1.4. Genuinely unsolved, squarely your niche; higher science risk — attempt *after* the SG-NEx benchmark (best-bet 1) gives you ground truth. **✓✓, Med–Large.**

*Plus the hardware-moat freebie:* §3.8 — a "does Boltz-2 + LMI4Boltz run on a 16 GB Pascal P100" recipe; trivial given your cu126/Pascal expertise; nobody else writes it.

---

## 7. Read-before-you-cite caveats

- **Therapeutic efficacy of splice neoantigens is genuinely debated.** The *existence* of shared/recurrent public neojunctions is solid (3-0); the ICI-response link is not. Hedge.
- **Ground-truth N is tiny** — ~13–20 validated splice-TCR triples across three papers; SNAF validated 5, Kwok et al. 2. Frame stage-3/4 benchmarks as *stress-test probes*, not powered evaluations.
- **Several load-bearing sources are 2024–2026 preprints** (Pep2Vec, the 17-tool MHC-I benchmark, CDR3-pLDDT, LMI4Boltz) — which is *why* they're reproduction targets; don't quote their numbers as settled. bioRxiv full-text returned HTTP 403 to automated fetch — figures were triangulated via search index, not independently re-derived.
- **"Splice peptides are just peptides."** Predictors are sequence-level and source-agnostic, so the *expected* benchmark result is that splice peptides land in the already-known "rare/unseen → near-random" regime. The novelty is being first to *measure it on an open splice slice*, not a surprise failure mode.
- **Terminology guard:** RNA splice-junction-derived ≠ proteasomal cis-spliced peptides (§2.6). Do not conflate.

### Claims that did NOT survive verification (do not cite)
- "SPLICE-neo applied to 11,892 TCGA tumors + 680 normals" — refuted 1-2 (scale claim unverified).
- "STMHCpan / BigMHC are the top MHC-I architectures" — refuted 1-2 (avoid asserting a single best tool).
- "MixMHC2pred / NetMHCIIpan-4.1 are the best open MHC-II predictors (AUC 0.95+)" — refuted 1-2.
- "recount3/Snaptron burden-vs-immune-correlate reanalysis is an under-exploited gap" — **killed** (refuted 2-3): it's already *feasible* (open junction layer + open Thorsson immune tables), so treat it as low-hanging fruit you can simply *do*, not a novel opening.

---

## 8. Method and provenance

Produced 2026-06-05 via two adversarially-verified web-research passes (deep-research harness + a targeted supplement workflow): 5 + 6 search angles, ~44 sources fetched, ~130 claims extracted, 50 claims verified by 3-vote panels. Raw verification transcripts are in the session workflow outputs. Confidence flags reflect skeptic vote splits, not author certainty. Refresh before relying on any preprint-derived figure.
