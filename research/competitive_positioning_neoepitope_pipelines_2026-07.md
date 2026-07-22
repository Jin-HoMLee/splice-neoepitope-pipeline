# Competitive positioning scan: splice-neoepitope pipelines (2026-07)

**Author:** PM (deep-research fan-out, adversarially verified).
**Date:** 2026-07-22.
**Method:** 6 search angles, 24 sources fetched, 115 claims extracted, top 25 put through 3-vote adversarial verification (23 confirmed, 2 refuted, 0 unverified).
**Question:** How do comparable cancer-neoepitope pipelines published or substantially updated in 2024-2026 structure themselves, and where does our pipeline sit competitively?

## Headline

Our pipeline sits in a crowded but genuinely differentiable niche.
Two direct splice-neoantigen tools bracket it, and two broader multi-source suites are ahead of us on presentation breadth.
The things that make us distinct are real and mostly unmatched: matched-normal *junction-level* filtering, TCRdock structural TCR-pMHC validation, and a *measured*-immunogenicity registry/calibrator.
The things we are behind on are also real: class-I-only presentation, single-source (splice only), and local-CPU chr22 scale.

## The comparable landscape

| Tool | Year / status | Source(s) | Presentation | Immunogenicity | TCR-pMHC structural | Tumor-specificity | Packaging |
|---|---|---|---|---|---|---|---|
| **SNAF** (closest peer) | Sci Transl Med 2024, maintained (v0.7.0) | Alt splicing (junction-first) + B-cell surface-protein antigens | Class I (NetMHCpan / MHCflurry) | DeepImmuno (sequence ML) + BayesTS tumor-specificity ranker | **No** (has a non-structural sequence/ML TCR module) | Large normal panel: GTEx >2500 samples / >54 tissues + TCGA paratumor | Python workflow, PyPI |
| **splice2neo** | Bioinformatics Adv 2024, maintained | Alt splicing only (5 event classes incl. intron retention, exitrons) | **None** (stops at peptide/ORF annotation) | None | None | Links junctions to causal somatic SNV/indel | R package |
| **NeoSplice** | Bioinformatics Adv 2022, older | Alt splicing only | Class I only (NetMHCpan 4.0, <500 nM) | MS immunopeptidome presentation only, no functional assay | None | - | repo |
| **pVACtools / pVACseq / pVACview** | Genome Med 2024 + v6/v7 2026, actively maintained | **Multi-source**: SNV/indel, fusions (pVACfuse), splicing (pVACsplice) | **Class I AND II** (MHCflurryEL, NetMHCpanEL, BigMHC-EL; NetMHCIIpanEL) | BigMHC_IM + DeepImmuno (off-the-shelf) | None | cis-splicing SNV-driven variants (not our RNA-seq junction detection) | Docker, field-standard |
| **OmniNeo** (strongest 2025 competitor) | Front Immunol 2025 | **Multi-source**: SNV/indel/frameshift/fusion/non-coding | **Class I AND II** (NetMHCpan 4.1 + NetMHCIIpan 4.0) | OmniNeo-CNN (self-reported AUC 0.88) | **Yes** - TCRmodel2 + PyMOL (demo case, not core step) | - | Nextflow |
| **NeoPredPipe** | BMC Bioinf 2019, older | SNV/indel only | Class I focus (netMHCpan) | Luksza recognition potential (A x R, computational) | None | - | Python |
| **ImmuneApp** | Nat Commun 2024 | (predictor building block, not a pipeline) | Class I only (CNN+LSTM+attention, 349,650 ligands) | AP sub-model (presentation) | None | - | web server + repo |

## Where we sit

### Genuine differentiators (medium confidence - cross-tool inference)

- **Matched-normal junction-LEVEL filtering.** Peers use large normal-tissue *reference panels* (SNAF's GTEx/TCGA), not paired per-patient junction filtering. This is a real design-axis difference, not just a parameter choice.
- **Presentation-likelihood scoring** (MHCflurry Class1PresentationPredictor) - in line with the field's better tools, not a differentiator by itself but not a gap either.
- **TCRdock structural TCR-pMHC validation.** Absent in SNAF, splice2neo, NeoSplice, and pVACtools. **Only OmniNeo has a structural step** (via TCRmodel2, not TCRdock), and even there it is framed as a demonstration case, not a mandatory core step. This is our strongest single differentiator.
- **A *measured*-immunogenicity registry/calibrator** (predicted -> measured T-cell immunogenicity). Unmatched: every competitor uses *pre-trained off-the-shelf* immunogenicity models (DeepImmuno, BigMHC_IM, OmniNeo-CNN). Calibrating on measured T-cell data is a step past the field - but note this is still "being built" on our side, so the claim is aspirational until it lands.

### Real gaps versus 2024-2026 trends (medium confidence)

- **Class-I-only presentation** while peers increasingly do class I AND II (pVACtools NetMHCIIpanEL, OmniNeo NetMHCIIpan 4.0). This is the clearest "behind the trend" gap, and it lines up exactly with the class-II direction Discussion #949 already proposed.
- **Single-source (splice only)** while the leading general pipelines are multi-source (pVACtools: SNV/indel/fusion/splice; OmniNeo: five source classes). Also directly the #949 "sources axis" broadening.
- **Local-CPU chr22 scale** while comparators ship container/workflow-managed genome-scale packaging (OmniNeo Nextflow, pVACtools Docker). Note: our Snakemake+conda packaging is itself reproducible-competitive; the gap is *scale*, not reproducibility discipline.

## Caveats (read before citing)

- **Two claims were refuted on re-vote and must not be relied on:** (1) that pVACtools' pVACsplice uses RegTools junctions from paired WES/RNA the way we do (1-2 refute) - pVACsplice targets cis-splicing SNV-driven variants, a *different mechanism*; and (2) that NeoPredPipe's presentation is definitively class-I-only via netMHCpan on 8-10mers (0-3 refute) - treat its class specifics cautiously.
- **Self-reported metrics:** OmniNeo-CNN's AUC 0.88 is from the tool's own paper on an internal held-out split, not external validation.
- **Split-vote (2-1) claims** resting partly on repo docs rather than paper body: pVACtools multi-predictor count, NeoSplice class-I designation.
- **The two positioning findings are cross-tool inferences (medium confidence)**, not single-source facts, and depend partly on our own self-described current state.
- **Coverage gap:** this scan verified SNAF, splice2neo, NeoSplice, pVACtools, OmniNeo, NeoPredPipe, and ImmuneApp. It did **not** systematically cover IRIS, ASNEO, antigen.garnish, MuPeXI, or NeoFuse - a follow-up sweep should close these, especially IRIS (splicing) and antigen.garnish (class-II/immunogenicity).

## Open questions worth a follow-up

1. Does our matched-normal junction-level filtering actually out-specify SNAF's large-panel GTEx/TCGA approach, or does per-patient pairing under-power on low-support junctions? (empirical, needs a benchmark)
2. Does OmniNeo's TCRmodel2 out/under-perform our TCRdock/AlphaFold path on TCR-pMHC accuracy - should we benchmark against or adopt TCRmodel2?
3. What is the realistic accuracy delta between a measured-immunogenicity calibrator and the off-the-shelf models (DeepImmuno, BigMHC_IM, OmniNeo-CNN)? This is what justifies the calibrator differentiator claim.
4. Do the uncovered tools (IRIS, ASNEO, antigen.garnish, MuPeXI, NeoFuse) change the gap picture, especially on class-II?

## Primary sources

- SNAF - https://www.science.org/doi/10.1126/scitranslmed.ade2886 , https://github.com/frankligy/SNAF
- splice2neo - https://github.com/TRON-Bioinformatics/splice2neo (Bioinformatics Advances, PMID 38863673)
- NeoSplice - https://academic.oup.com/bioinformaticsadvances/article/2/1/vbac032/6581739 , https://github.com/Benjamin-Vincent-Lab/NeoSplice
- pVACtools v6 - https://arxiv.org/abs/2606.26659 ; Genome Medicine 2024 - https://link.springer.com/article/10.1186/s13073-024-01384-7 ; https://github.com/griffithlab/pVACtools
- OmniNeo - https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2025.1727642/full
- NeoPredPipe - https://www.biorxiv.org/content/10.1101/409839.full.pdf , https://github.com/MathOnco/NeoPredPipe
- ImmuneApp - https://www.nature.com/articles/s41467-024-53296-0 , https://github.com/bsml320/ImmuneApp
