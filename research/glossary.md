# Glossary

Project-relevant abbreviations and acronyms. The pipeline mixes biology, ML, bioinformatics, and cloud infrastructure - entries here aim to make onboarding (or returning to the project after time away) easier than re-asking each time.

**Format:** `**ABBREV** - Expansion. One-line context. *Domain: bio | ml | cloud | pipeline | stats.*`

**Scope:** project-relevant terms only. Skip generic web/programming acronyms (HTTP, REST, etc.).

---

## A

**AC** - Acceptance Criteria. The checklist in an issue body defining what "done" means for that issue. Per the closure ritual, all AC boxes must be ticked (`- [x]`) or comment-deferred before closing. *Domain: project.*

**ADR** - Architecture Decision Record. A short, dated in-repo document capturing one significant technical decision: context, the decision (including what is explicitly *not* done), alternatives, and consequences/revisit-trigger. Records *why* a choice was made so it isn't silently re-litigated. Lives under `docs/adr/` (being stood up in [Issue #777](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/777)). *Domain: project.*

**AF-Multimer** - AlphaFold-Multimer. DeepMind extension of AlphaFold (Evans et al. 2021) trained on protein complexes; predicts both monomeric structure and inter-chain interactions. Used by TCRdock as the structural backbone. *Domain: bio.*

**AF_confidence** - AlphaFold confidence score. Combines pTM (predicted Template Modeling) and ipTM (interface predicted TM); typically `0.8·ipTM + 0.2·pTM`. Standard global-quality metric for AF-Multimer model selection. *Domain: bio.*

**AF3** - AlphaFold 3 (DeepMind, Abramson et al. *Nature* 2024). Successor to AlphaFold 2 / AF-Multimer; co-folding architecture predicting protein + ligand + nucleic acid + ion complexes in one pass via diffusion-based atom placement. Best-overall TCR-pMHC structure predictor in [Lu et al. 2025](https://www.biorxiv.org/content/10.64898/2025.11.30.691400v1.full) (median DockQ 0.636 class I / 0.679 class II on 70 unseen complexes), motivating [Issue #316](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/316) as the TCRdock-backend modernization track. *Domain: bio.*

**AFDB** - AlphaFold Database (alphafold.ebi.ac.uk). DeepMind/EBI's public repo of ~200M predicted protein structures. *Domain: bio.*

**AP** - Access Point. In HTCondor topology, the machine where jobs are submitted, files staged, and conda envs built. Distinct from the EP (execution point) when there is no shared FS; the AP→EP env-deployment gap is the same problem cloud executors (Google Batch) face. *Domain: cloud.*

## C

**CAPRI** - Critical Assessment of PRedicted Interactions. Community blind-prediction challenge for protein-protein / peptide-MHC docking (since 2001); CAPRI quality bands are the standard yardstick for docking models - DockQ is the continuous reformulation of these bands. Bands map to DockQ ranges: **HQ** (high quality, DockQ ≥ 0.80), **MQ** (medium quality, 0.49 ≤ DockQ < 0.80), **AQ** (acceptable quality, 0.23 ≤ DockQ < 0.49), **incorrect** (DockQ < 0.23). *Domain: bio.*

**CAR-T** - Chimeric Antigen Receptor T-cell therapy. Adoptive cell therapy where patient T cells are transduced with a synthetic chimeric receptor (scFv + CD3ζ + costim domains); recognizes surface antigen directly, MHC-independent. FDA-approved CD19/BCMA-targeting products for B-cell malignancies and myeloma. Distinct from TCR-T (native MHC-restricted TCR). *Domain: bio.*

**CDR3** - Complementarity-Determining Region 3. Hypervariable loop on each TCR chain (α and β) that makes direct contact with the peptide in pMHC; the primary specificity-encoding region of the TCR. Average `pLDDT` across CDR3 residues used as a TCR-pMHC docking quality signal (Lu et al. 2025, [Issue #316](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/316)). *Domain: bio.*

**CHTC** - Center for High Throughput Computing (UW-Madison). Develops and maintains HTCondor; took over maintenance of Snakemake's HTCondor executor plugin at the 2026 Snakemake Hackathon. *Domain: cloud.*

**CoREST** - Corepressor of REST (RE1-Silencing Transcription factor). Chromatin co-repressor complex (LSD1-HDAC1-CoREST core scaffold); originally characterised as a neuronal-gene silencer, recently shown to also regulate RNA-processing genes. Pharmacological inhibition by **corin** in melanoma cell lines disrupts spliceosome activity → induces immunogenic splice-derived neoantigens, reactivates anti-PD-1 in cold tumors ([Fisher et al. 2025](https://insight.jci.org/articles/view/190287)). *Domain: bio.*

**CR** - Complete Remission (a.k.a. Complete Response). Oncologic response criterion: full disappearance of detectable disease following therapy; distinguished from PR (partial response), SD (stable disease), PD (progressive disease) in RECIST 1.1 solid-tumor criteria. The weakest TNBC vaccine responder in Sahin et al. 2026 achieved CR on subsequent anti-PD-1 rescue. *Domain: bio.*

**CRC** - Colorectal Cancer. Common GI adenocarcinoma; clinically subtyped by microsatellite status (MSI vs MSS) which drives ICI responsiveness - MSI-high tumors respond, MSS broadly does not. Recurrent test bed for splice-derived TSA discovery via proteogenomics (e.g. [Hayer et al. 2026, *Mol Cell Proteomics*](https://www.mcponline.org/article/S1535-9476(26)00077-0/fulltext)). *Domain: bio.*

**CVE** - Common Vulnerabilities and Exposures. Standardised identifier system for publicly disclosed security flaws (CVE-YYYY-NNNN format); maintained by MITRE, used industry-wide for tracking known vulns. *Domain: security.*

## D

**DAG** - Directed Acyclic Graph. Snakemake's compiled job-dependency graph; computed from rule wildcards + I/O files, traversed for scheduling. Rendered via `scripts/visualize_dag.sh` after a dry-run. *Domain: pipeline.*

**DDA** - Data-Dependent Acquisition. Mass-spec acquisition mode where the instrument picks the top-N most intense MS1 precursors per cycle and fragments only those. Stochastic and intensity-biased - low-abundance peptides (including most neoantigens) are routinely missed across replicate runs. Historical default; superseded by DIA-MS for immunopeptidomics. *Domain: bio.*

**DIA-MS** - Data-Independent Acquisition Mass Spectrometry. Acquisition mode where the instrument fragments every precursor in a stepped m/z window, regardless of intensity - no peptide is "missed", but resulting MS2 spectra are multiplexed and require a reference spectral library for deconvolution. Public libraries cover canonical proteomes; patient-specific neoantigens need custom libraries (e.g. Pepyrus, [Manakongtreecheep et al. 2026](https://www.nature.com/articles/s41587-026-03003-9)). *Domain: bio.*

**DLT** - Dose-Limiting Toxicity. Pre-defined adverse-event severity threshold (typically CTCAE grade ≥3) in Phase 1 dose-escalation trials; reaching DLT caps further escalation. "No DLTs" is the standard Phase 1 safety endpoint reported in neoantigen-vaccine trials (e.g. GNOS-PV01, [Garfinkle et al. 2026](https://www.nature.com/articles/s43018-026-01163-w)). *Domain: bio.*

**DockQ** - Docking Quality. Continuous 0-1 score combining Fnat (fraction of native contacts), iRMSD (interface RMSD), and LRMSD (ligand RMSD); maps onto CAPRI bands; used to evaluate predicted protein-protein / TCR-pMHC complex structures against experimental ground truth. *Domain: bio.*

## E

**ENCODE** - Encyclopedia of DNA Elements. NIH consortium (since 2003) cataloging functional genome features via standardized high-throughput assays (RNA-seq, ChIP-seq, ATAC-seq, Hi-C, CAGE, bisulfite-seq, etc.). "ENCODE-style assays" = this standard set of functional genomics readouts, regardless of which project produced the data. *Domain: bio.*

**EP** - Execution Point. In HTCondor topology, the machine where a job actually runs. In no-shared-FS topologies the conda env built on the AP must be transported to (or rebuilt on) the EP - the same constraint cloud workers face. *Domain: cloud.*

**ESM-IF1** - Evolutionary Scale Modeling Inverse Folding (v1). Meta AI model (Hsu et al. 2022) trained on 12M protein structures; given a backbone, infers compatible sequences. Used as a structure-aware embedding source by structure-informed TCR-pMHC scorers (e.g. NetTCR-struc). *Domain: ml.*

## F

**FDR** - False Discovery Rate. Expected proportion of false positives among reported discoveries; standard control in MS peptide identification (typically 1% peptide-spectrum match FDR + 1% peptide FDR via target-decoy estimation) and in genomics multiple-testing corrections (Benjamini-Hochberg). *Domain: stats.*

**FEP** - Free Energy Perturbation. Physics-based binding-affinity prediction via alchemical molecular-dynamics simulations: gradually morph one ligand into another (or into nothing) while tracking the free-energy change with thermodynamic integration. Accuracy gold-standard for small-molecule lead-optimization in pharma, but computationally expensive (CPU-hours to GPU-days per ligand-pair). The "1000× faster than FEP" speed claim in the Boltz-2 paper ([Issue #188](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/188)) refers to this small-molecule benchmark - not a TCRdock-style structure-prediction speedup. *Domain: bio.*

**FOXA2** - Forkhead Box A2. Pioneer transcription factor that binds nucleosomal DNA at enhancer / promoter elements; canonical regulator of endoderm specification (gut, liver, pancreas) and adult hepatocyte / β-cell identity. In PDAC, dysregulated FOXA2 binding drives alternative promoter usage that generates neoantigen-encoding tumor-specific transcripts (neoTSTs) - flagged as a major non-mutation NA source ([Zhao et al. 2026 NeoAPP](https://www.biorxiv.org/content/10.64898/2026.02.10.705024v1.full)). *Domain: bio.*

## G

**GATK** - Genome Analysis Toolkit (Broad Institute). Canonical somatic/germline variant calling suite; reference for the Panel of Normals (PoN) concept reused by neoepitope filtering. *Domain: bio.*

**GBM** - Glioblastoma (Multiforme). WHO grade 4 astrocytoma; most aggressive primary adult brain tumor (~15-mo median OS post standard-of-care surgery + radiation + temozolomide). MGMT-unmethylated subset has the worst prognosis (resistant to temozolomide); subject of the GNOS-PV01 personalized DNA vaccine Phase 1 ([Garfinkle et al. 2026](https://www.nature.com/articles/s43018-026-01163-w)). *Domain: bio.*

**GCP** - Google Cloud Platform. Formerly hosted this pipeline's VMs and storage bucket (zone `europe-west4-a`); the GCP stack was decommissioned 2026-06-26 ([Issue #854](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/854)) to stop charges before the free-trial lapse. Current compute is a local CPU-only core; the funded revival path targets a paid GPU provider (RunPod, see `docs/migration_runbook.md`). *Domain: cloud.*

**GCS** - Google Cloud Storage. GCP's object store; formerly held this pipeline's results, logs, and reference data in `gs://splice-neoepitope-project`. Decommissioned with the GCP stack 2026-06-26 ([Issue #854](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/854)) after project data was preserved on Cloudflare R2 - see **[R2](#r)**. *Domain: cloud.*

**GKE** - Google Kubernetes Engine. Managed Kubernetes on GCP. *Domain: cloud.*

**GNOS-PV01** - Personalized DNA neoantigen vaccine (Geneos Therapeutics + WashU Medicine). Encodes up to 40 patient-specific neoantigens on a DNA backbone, delivered intramuscularly with electroporation; Phase 1 in MGMT-unmethylated GBM ([Garfinkle et al., Nat Cancer 2026](https://www.nature.com/articles/s43018-026-01163-w)). Adds **DNA** to the platform diversity census alongside peptide / mRNA / DC. *Domain: bio.*

**GPS** - Genotype Presentation Score. Project-specific scoring formula combining per-allele MHCflurry presentation percentiles into a single rank for the patient's genotype; primary ranking key for top neoepitope candidates. *Domain: pipeline.*

**GTEx** - Genotype-Tissue Expression project. NIH consortium with healthy-tissue RNA-seq across ~50 tissues from ~1,000 donors; canonical reference for tissue-matched normal expression panels. *Domain: bio.*

**GVP-GNN** - Geometric Vector Perceptron Graph Neural Network. Architecture (Jing et al. 2021) for learning on protein 3D structures; respects rotational/translational equivariance natively. Used in NetTCR-struc and similar structure-aware TCR-pMHC scorers. *Domain: ml.*

## H

**HCC** - Hepatocellular Carcinoma. Most common primary liver cancer; typically arises on a background of chronic liver injury (HBV/HCV, cirrhosis, NASH). A major test bed for transcriptome-derived (non-mutation) neoantigen discovery - Lin et al. 2025 mined ~60 neoTSTs/patient from novel splice junctions, intron retention, and TE activation across 1,013 HCC patients ([bioRxiv 2025](https://www.biorxiv.org/content/10.64898/2025.12.07.692877v1)). *Domain: bio.*

**HLA** - Human Leukocyte Antigen. Human MHC class I/II proteins (HLA-A/B/C class I; HLA-DR/DP/DQ class II); patient typing via OptiType drives MHCflurry's per-allele predictions. *Domain: bio.*

**HLA-LOH** - HLA Loss of Heterozygosity. Tumor immune-escape mechanism: somatic deletion or copy-neutral loss of one HLA allele, narrowing the peptide-presentation repertoire and shielding the tumor from neoantigen-specific T cells. Detected from WES (e.g. LOHHLA); candidates predicted on a lost allele are no longer presented. *Domain: bio.*

**HPC** - High-Performance Computing. Cluster compute systems with job schedulers (SLURM, LSF, HTCondor); distinct from cloud "on-demand" compute models (GCP, AWS Batch). The pipeline core currently runs on a local CPU-only host (neither cloud nor HPC) post-GCP-decommission; the funded GPU-revival path targets cloud (RunPod), not HPC. *Domain: cloud.*

**HTCondor** - High-Throughput Condor. Distributed batch scheduler from UW-Madison's CHTC; dominant in physics/astronomy (e.g. CERN ATLAS). Supports no-shared-FS topologies natively - relevant precedent for cloud-executor design. *Domain: cloud.*

## I

**IC50** - half-maximal Inhibitory Concentration. Concentration (in nM) at which a peptide displaces 50% of a reference ligand from MHC; classic affinity metric, retained as informational column alongside the percentile-based presentation score. *Domain: bio.*

**ICI** - Immune Checkpoint Inhibitor. Class of monoclonal antibodies blocking inhibitory T-cell receptors (anti-CTLA-4, anti-PD-1, anti-PD-L1) to release brakes on anti-tumor immunity; standard-of-care across many solid tumors. "Immune cold" tumors (low T-cell infiltrate) often fail to respond - motivating combination strategies including neoantigen vaccines and pharmacologic mis-splicing inducers (corin; [Fisher et al. 2025](https://insight.jci.org/articles/view/190287)). *Domain: bio.*

**IP-MS** - Immunoprecipitation Mass Spectrometry. Core workflow of immunopeptidomics: pan-MHC or allele-specific antibodies affinity-purify HLA-peptide complexes from cell or tissue lysate, peptides are acid-eluted and identified by LC-MS/MS (DDA or DIA-MS). The empirical complement to in-silico binding prediction; ground truth for "is this peptide actually presented?" ([Hayer et al. 2026](https://www.mcponline.org/article/S1535-9476(26)00077-0/fulltext); MHC1-TIP, Pepyrus). *Domain: bio.*

**IR** - Intron Retention. Failure of the spliceosome to excise an intron from pre-mRNA, producing a mature transcript that retains intronic sequence. Elevated in many tumors vs matched normals; IR-derived peptides span exon/intron junctions and act as a non-mutation neoantigen class - ~30% of IR-predicted CRC epitopes elicit measurable CD8+ T cell responses in functional assays ([Manoharan et al. 2026](https://www.nature.com/articles/s41598-026-43687-2)). Distinct from canonical AS-NAs (exon skipping, cryptic splice sites) tracked by SNAF / SpliceMutr. *Domain: bio.*

**IV** - Intravenous. Drug delivery route via direct venous injection; chosen for personalized mRNA cancer vaccines to deliver LNP-mRNA payloads systemically (Sahin et al. 2026; Rojas et al. 2023). *Domain: bio.*

## K

**KIR** - Killer-cell Immunoglobulin-like Receptor. Family of NK-cell receptors whose dominant ligands are HLA-C alleles; HLA-C-KIR engagement inhibits NK cytotoxicity under "missing-self" surveillance. Relevant here because HLA-C is studied primarily through this NK biology rather than T-cell presentation, which explains why TCR repertoire databases (VDJdb, McPAS-TCR) have thin paired α/β coverage for C-locus alleles. See [Issue #204](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/204) empirical coverage check. *Domain: bio.*

## L

**LNP-mRNA** - Lipid Nanoparticle-formulated mRNA. Delivery format: ionizable-lipid encapsulation protects mRNA from RNases and mediates dendritic-cell uptake. Foundation of COVID-19 mRNA vaccines (BioNTech/Moderna) and current personalized cancer-vaccine trials (autogene cevumeran; Sahin et al. 2026). *Domain: bio.*

**LSD1-HDAC1-CoREST** - Core chromatin co-repressor complex: lysine-specific demethylase 1 (LSD1, KDM1A) + histone deacetylase 1 (HDAC1) + CoREST scaffold. Coordinates H3K4 demethylation with H3 deacetylation for transcriptional repression; the molecular target of the small-molecule inhibitor **corin** ([Fisher et al. 2025](https://insight.jci.org/articles/view/190287)). Disruption reroutes the complex away from RNA-processing genes → mis-splicing → immunogenic neopeptides. *Domain: bio.*

**LSF** - Load Sharing Facility (IBM Spectrum). Commercial HPC batch scheduler; used at DKFZ Heidelberg among others. Snakemake has a dedicated LSF executor plugin. *Domain: cloud.*

**LSP** - Language Server Protocol. Microsoft spec (2016) decoupling code-intelligence backends (completion, hover, go-to-definition, diagnostics) from editor frontends - one server, many editors. Claude Code plugins can register LSP servers alongside MCP servers; `/plugin Discover` previews these before install (Claude Code 2.1.145). *Domain: cloud.*

## M

**MCP** - Model Context Protocol. Anthropic's open spec (late 2024) for AI agents to discover and invoke external tools / data sources via per-purpose MCP servers; the standardized way to expose resources to LLMs without bespoke per-tool integration. *"USB-C for AI agents."* Used by Claude Code, the GitHub MCP server, filesystem servers, etc. *Domain: ml.*

**MGMT** - O⁶-methylguanine-DNA methyltransferase. DNA repair enzyme; promoter methylation silences the gene and sensitizes glioma to alkylating chemotherapy (temozolomide). **MGMT-unmethylated** GBM has worse prognosis (temozolomide-resistant) - the poor-response subset targeted by the GNOS-PV01 vaccine trial ([Garfinkle et al. 2026](https://www.nature.com/articles/s43018-026-01163-w)). *Domain: bio.*

**MHC** - Major Histocompatibility Complex. Cell-surface proteins that present peptides to T cells; class I (all nucleated cells, presents endogenous peptides) drives this pipeline. Human MHC = HLA. *Domain: bio.*

**MHC1-TIP** - MHC class I 1-Tube Immunopeptidomics. Low-input, single-tube MHC-I ligandome workflow ([Dollinger et al. 2026, *Comm Bio*](https://www.nature.com/articles/s42003-026-09570-6)); scales to cell lines, patient-derived organoids, and sub-mg ex-vivo tumor fragments - replaces traditional workflows requiring hundreds of millions of cells. Primary RCC application revealed widespread **intratumoral heterogeneity in antigen presentation that is poorly correlated with source protein expression** - relevant to the RNA-Seq-abundance ≠ surface-presentation framing in this pipeline. *Domain: bio.*

**MLP** - Multi-Layer Perceptron. Feedforward neural network of fully-connected layers with nonlinear activations - the generic "dense network" building block. In ImmunoStruct, a biochemical-feature MLP is fused with sequence and structure modules to predict pMHC immunogenicity ([Issue #610](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/610)). *Domain: ml.*

**mRNA prime + peptide-target boost** - Specific heterologous prime-boost format inside the Prime-Target regimen: an mRNA-encoded antigen primes the T cell response, a synthetic peptide of the same antigen amplifies it. Heterologous formats consistently outperform same-format boosts in vaccinology - each modality engages distinct APC subsets and avoids anti-vector immunity that dampens homologous boosters. *Domain: bio.*

**MSA** - Multiple Sequence Alignment. Stack of homologous protein (or DNA) sequences aligned by position; columns reveal residue-level evolutionary conservation and covariation, the signal that drove AlphaFold's accuracy breakthrough (coevolving residues = likely 3D contacts). The standard upstream input to co-folding / structure-prediction models (AlphaFold, AF-Multimer, Boltz, Chai-1, RosettaFold-AllAtom) - built per-query via jackhmmer / HHblits against UniRef + BFD / MGnify. *Domain: bio.*

**MSI** - Microsatellite Instability. Tumor phenotype caused by defective DNA mismatch repair (MMR); hypermutated, high TMB, characteristically responsive to ICI. MSI-high is a tissue-agnostic FDA-approved indication for anti-PD-1 (pembrolizumab). Non-canonical splicing contributes substantively to the CRC immunopeptidome in this class ([Hayer et al. 2026](https://www.mcponline.org/article/S1535-9476(26)00077-0/fulltext)). *Domain: bio.*

**MSS** - Microsatellite Stable. The MMR-proficient majority of solid tumors; lower TMB, broadly ICI-resistant. The challenge case for neoantigen vaccines - fewer canonical somatic mutations motivates non-canonical (splice, RNA-edit, retroelement) TSA discovery; splice-derived antigens are detectable in MSS CRC at rates comparable to MSI ([Hayer et al. 2026](https://www.mcponline.org/article/S1535-9476(26)00077-0/fulltext)). *Domain: bio.*

## N

**NJ** - Neojunction. Tumor-specific or recurrent splice junction absent from canonical normal-tissue annotation; the upstream substrate from which splice neoepitopes are translated. Notation introduced by Nejo et al. (*Nat. Med.* 2023). Kwok et al. additionally use **NEJ** (neoepitope-encoding junction) for the validated subset yielding a presented peptide (e.g. NEJ<sub>GNAS</sub>, NEJ<sub>RPL22</sub>). *Domain: bio.*

**neoTSTs** - Neoantigen-encoding Tumor-Specific Transcripts. Aberrant tumor mRNAs that, when translated, yield peptides absent from normal-tissue proteomes - a broad parent class spanning splice variants, intron retention, TE-fusion transcripts, and alternative-promoter isoforms. PDAC reference scale: median ~351 neoantigens / 56 neoTSTs per sample, vastly exceeding SNV-derived NAs in cold tumors ([Zhao et al. 2026 NeoAPP](https://www.biorxiv.org/content/10.64898/2026.02.10.705024v1.full)). *Domain: bio.*

## O

**OGS** - Open Grid Scheduler. Open-source community fork of SGE; `qsub`-compatible. *Domain: cloud.*

**OKR** - Objectives and Key Results. Goal-setting framework popularised by Intel/Google: one qualitative *Objective* + 3-5 quantitative *Key Results*, usually quarterly. Common in software shops; we use S1-S7 lifecycle stages + iteration capacity instead. *Domain: project.*

**OOD** - Out-of-distribution. Inputs whose features fall in low-density regions of the model's training distribution; predictions on OOD inputs are typically over-confident yet unreliable since the model never learned to constrain them. Detected operationally via density estimators, embedding-space distance to training samples, or ensemble disagreement. Relevant for splice-junction-spanning peptides - systematically under-represented in canonical-proteome MHC training data. *Domain: ml.*

**OTEL** - OpenTelemetry. Vendor-neutral observability framework under the Cloud Native Computing Foundation (CNCF): specs + libraries for distributed traces, metrics, and logs, exportable to backends like Honeycomb, Jaeger, Tempo. Claude Code 2.1.145 added `agent_id` and `parent_agent_id` to `claude_code.tool` spans so nested agent calls show correct trace parenting (relevant if we ever wire trace export from PM/Sci/Dev sessions). *Domain: cloud.*

## P

**PAE** - Predicted Aligned Error. AlphaFold/TCRdock per-residue-pair error estimate (in Å); summarised globally to gauge inter-domain confidence (e.g. relative MHC↔TCR docking confidence). *Domain: bio.*

**PD-1** - Programmed cell Death protein 1 (CD279). Inhibitory T-cell coreceptor; engages PD-L1/PD-L2 on tumor or APCs to dampen T-cell activity. Anti-PD-1 antibodies (pembrolizumab, nivolumab) are the workhorse ICI class; recurrent combination partner in neoantigen-vaccine and splice-modulator trials (e.g. TN2008 + anti-PD-1 in resistant models, [Lin et al. 2025](https://www.nature.com/articles/s41392-024-02118-2)). *Domain: bio.*

**PDB** - Protein Data Bank ([rcsb.org](https://www.rcsb.org)). Public archive of experimentally determined 3D structures (~220k entries); training source for AlphaFold / Boltz / HERMES; reference ground truth for DockQ scoring. *Domain: bio.*

**pLDDT** - predicted Local Distance Difference Test. AlphaFold/TCRdock per-residue confidence score (0-100); per-residue version of the global `AF_confidence`. >90 = high, >70 = confident, <50 = low confidence. CDR3-region average used as a TCR-pMHC docking quality flag (Lu et al. 2025, [Issue #316](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/316)). *Domain: bio.*

**PLM** - Protein Language Model. Transformer trained on raw protein sequences (e.g. ESM-2/3, ProtT5, ProGen) to learn residue-level embeddings that capture structural and functional context - analogous to NLP language models but on amino-acid alphabets. PLM-based structure predictors (e.g. ESMFold, OmegaFold) skip the MSA step and embed single sequences directly; faster than MSA-based methods (AlphaFold, Boltz) but typically less accurate on complexes - Lu et al. 2025 benchmark classifies TCR-pMHC predictors as MSA-based / PLM-based / docking-based. *Domain: ml.*

**pMHC** - peptide-MHC. The complex of a peptide bound in the MHC groove; the molecular target a TCR recognises. "TCR-pMHC interaction" = the central event in adaptive immune recognition. *Domain: bio.*

**PoN** - Panel of Normals. Aggregated reference set of normal samples (originally GATK's somatic variant calling concept); used to filter recurrent germline variation / artefacts that single matched normals miss. *Domain: bio.*

**Prime-Target** - Neoantigen-vaccine regimen pairing an mRNA *prime* with a synthetic peptide *target* boost of the same antigen. The combined heterologous format produces stronger systemic T cell immunity inside immunologically "cold" tumors than either modality alone in mouse models ([Prime-Target preprint, bioRxiv 2026](https://www.biorxiv.org/content/10.64898/2026.01.13.699214v1)). Manuscript DISCUSSION hook on downstream use of our predicted AS-neoepitopes; contrasted against single-modality regimens. *Domain: bio.*

**PSI** - Percent Spliced In. Fraction of transcripts that include a given alternative exon/junction (vs skip); the standard quantitative readout of splicing in RNA-seq analyses. *Domain: bio.*

**PSR** - Positive Sample Rate. Percentage of samples in a cohort expressing a splice junction at a read frequency ≥1% relative to the canonical junction (Kwok et al., *Nature* 2025); thresholded as `PSR_TCGA ≥ 10%` (tumor-recurrence) and `PSR_GTEx < 1%` (normal-exclusion). Distinct from **PSI** (within-sample inclusion fraction): PSR is a binary cross-sample detection rate at a fixed relative-frequency threshold; PSI is a continuous within-sample metric. *Domain: bio.*

## R

**R2** - Cloudflare R2. S3-compatible object store; the current home for this pipeline's results, logs, and reference data after the GCP/**[GCS](#g)** decommission (2026-06-26, [Issue #854](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/854)). Accessed via the s3v4 API with `region=auto`. *Domain: cloud.*

**RCE** - Remote Code Execution. Vulnerability class where an attacker runs arbitrary code on a remote system, typically via input-validation or memory-safety bugs (e.g. CVE-2026-3854's crafted git push trigger). *Domain: security.*

**RMSD** - Root Mean Square Deviation. Average distance (in Å) between corresponding atoms after optimal superposition; the standard structural-similarity metric. Variants used in DockQ: backbone-only, all-atom, interface-only (iRMSD), ligand-only (LRMSD). *Domain: bio.*

**RunPod** - GPU cloud provider (per-second billing on L4 / Ada-class GPUs). The forward paid-provider target for reviving GPU workloads (TCRdock structural validation, GPU MHCflurry) post-GCP; migration plan in `docs/migration_runbook.md`. *Domain: cloud.*

## S

**SGE** - Sun Grid Engine. Classic HPC job scheduler from Sun Microsystems (now Oracle); ancestor of UGE and OGS. Snakemake covers the SGE/UGE/OGS family via `snakemake-executor-plugin-sge`. *Domain: cloud.*

**single-modality regimen** - Neoantigen vaccine using one delivery format throughout (mRNA-only, peptide-only, or DNA-only). The default in Phase 1 trials so far (autogene cevumeran / Sahin et al. 2026 TNBC, GNOS-PV01 GBM); contrasted against heterologous prime-boost approaches like Prime-Target where format diversity within a single regimen broadens the T cell response. *Domain: bio.*

**SLURM** - Simple Linux Utility for Resource Management. The dominant HPC job scheduler in academia. *Domain: cloud.*

**SRSF1** - Serine/arginine-Rich Splicing Factor 1. Canonical SR-protein family member; binds exonic splicing enhancers, regulates alternative 5′ splice-site selection. Oncogenic across multiple cancers; small-molecule inhibitor **TN2008** boosts anti-PD-1 response in resistant models via a dual hit - restoring CD8+ T-cell glycolysis + reprogramming tumor c-Jun/c-myc/JunB transcription ([Lin et al. 2025, *Sig Transduct Target Ther*](https://www.nature.com/articles/s41392-024-02118-2)). *Domain: bio.*

## T

**TCR** - T-Cell Receptor. Heterodimeric (α/β) receptor on T cells that recognises peptide-MHC complexes; modelled here via TCRdock for the top neoepitope candidate. *Domain: bio.*

**TCR-T** - TCR-engineered T-cell therapy. Adoptive cell therapy where autologous patient T cells are transduced with a tumor-reactive TCR (typically discovered in HLA-matched healthy donors to bypass tumor-induced tolerance). MHC-restricted; distinct from CAR-T (synthetic chimeric receptor) and TIL therapy (non-engineered expanded patient T cells). *Domain: bio.*

**TE** - Transposable Element. Repetitive mobile DNA (LINE, SINE, LTR / HERV families) comprising ~45% of the human genome; mostly silenced by methylation in normal tissue. Cancer-associated demethylation reactivates TEs, producing chimeric TE-gene transcripts that encode tumor-specific peptides - an emerging non-mutation neoantigen source alongside splice variants and IR ([Zhao et al. 2026 NeoAPP](https://www.biorxiv.org/content/10.64898/2026.02.10.705024v1.full)). *Domain: bio.*

**TF** - Transcription Factor. DNA-binding regulatory protein controlling gene expression; binding sites are one of the modalities ENCODE-style assays (ChIP-seq) and predictors like AlphaGenome map. *Domain: bio.*

**TIL** - Tumor-Infiltrating Lymphocyte. T cell isolated from a tumor biopsy; expanded ex vivo and re-infused as TIL therapy (FDA-approved 2024 for advanced melanoma as lifileucel/Amtagvi). Distinct from TCR-T (engineered TCR) and CAR-T (synthetic chimeric receptor). *Domain: bio.*

**TNBC** - Triple-Negative Breast Cancer. Breast cancer subtype lacking ER, PR, and HER2 expression; lacks targeted hormonal / HER2 therapy options. ~10-15% of breast cancers; relatively higher tumor mutation burden makes it a recurrent personalized-vaccine trial target (Sahin et al., *Nature* 2026). *Domain: bio.*

**TPM** - Transcripts Per Million. Normalised RNA expression unit; reads per kb of transcript scaled so per-sample values sum to 10⁶, making cross-sample comparisons direct (unlike RPKM/FPKM). Conventional thresholds: ~1 TPM = expressed, <0.5 TPM = below detection floor (functionally silent in that tissue), >10 TPM = moderately-highly expressed. Used as a tier-1 normal-tissue QC for public-neoantigen prioritisation (e.g. Zhang et al. 2026: <0.5 TPM in normal tissues). *Domain: bio.*

**TSA** - Tumor-Specific Antigen. Peptide present on tumor cells but absent from normal tissue; the strict subset of tumor antigens whose immunogenic recognition is unlikely to drive autoimmunity. Distinct from **TAA** (tumor-associated antigen, also expressed in some normal lineages). Splice neoepitopes our pipeline targets are the splice-junction TSA class; junction-spanning TSAs are MS-detectable in MSI and MSS CRC ([Hayer et al. 2026](https://www.mcponline.org/article/S1535-9476(26)00077-0/fulltext)). *Domain: bio.*

## U

**UGE** - Univa Grid Engine. Commercial fork of SGE (Univa, acquired by Altair 2020); drop-in `qsub`-compatible. *Domain: cloud.*

## V

**VAE** - Variational Autoencoder. Generative neural network that encodes inputs into a probabilistic latent space (regularized toward a prior) and decodes samples back, learning a smooth low-dimensional representation. In ImmunoStruct a pMHC-sequence VAE supplies the sequence-embedding branch, fused with a structure graph-transformer + biochemical MLP for immunogenicity prediction ([Issue #610](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/610)). *Domain: ml.*

**V(D)J** - Variable / Diversity / Joining segment recombination. Somatic DNA rearrangement at TCR and immunoglobulin loci during lymphocyte development: one V gene segment, optionally one D segment (β/δ chains and IgH only), and one J segment are spliced together to form the variable domain. Junctional diversity is amplified by exonuclease trimming and TdT-mediated N-nucleotide additions at the V-D and D-J joints - the source of CDR3 hypervariability and the reason CDR3 (V-D-J joint) is orders of magnitude more diverse than germline-encoded CDR1/CDR2 (within a single V segment). Theoretical TCR diversity ~10¹⁵-10²⁰; functional repertoire per individual ~10⁷-10⁸. *Domain: bio.*

## W

**WES** - Whole Exome Sequencing. Sequencing restricted to protein-coding regions (~1-2% of the genome); cheaper than WGS, sufficient for variant calls in coding regions. *Domain: bio.*

**WGS** - Whole Genome Sequencing. Sequencing of the entire genome (~3 Gbp human); alternative input for upstream filtering when matched-normal RNA-seq is missing. *Domain: bio.*
