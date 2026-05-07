# Lab Notebook — Scientist

Per-role lab notebook for Scientist sessions. Started 2026-05-06 as part of the lab-notebook split (see `lab_notebook.md` freeze note for rationale).

Format and rules unchanged from the unified notebook — see `shared/feedback_lab_notebook.md`. New date sections at the TOP; new time sections at the TOP of the date block; entries are immutable once committed.

---

## 2026-05-07

### 16:19 UTC — Editor: Scientist

#### [Sub-Issue #280](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/280) DISCUSSION — public-vs-personalised splice neoantigen axis (Kwok et al.)

**Author-name correction.** Yesterday's frozen [`lab_notebook.md`](research/lab_notebook.md) entry at lines 23 + 25 attributes [the *Nature* 2025 paper](https://www.nature.com/articles/s41586-024-08552-0) to "Pan et al." — the correct first author is **Darwin W. Kwok**. Same paper, same DOI (`10.1038/s41586-024-08552-0`), Zotero key `5ZT8KC8X`. The stale attribution propagated into [Issue #280](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/280)'s title + body and `research/news_log.md` line 23 before catching it during today's verification work. Per the lab notebook immutability rule (`shared/feedback_lab_notebook.md`), the frozen-snapshot entries are not editable; this note serves as the corrected reference for future lookups. Fixed today: [Issue #280](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/280) title + body (Pan → Kwok, plus expanded ACs covering the threshold-tradeoff scope), and `research/news_log.md` line 37 (rolled into this PR per yesterday's option-2 decision).

**Verification path that surfaced the correction.** The user pushed back on yesterday's first-pass framing — *"biologically I don't know if junctions that are tumor-recurrent across patients AND show up sporadically in normal tissue are really relevant as neoepitopes? Do Pan et al. show that such junctions result in potent neoepitopes?"* — driving a structured re-read of the paper. WebFetch on PubMed by DOI (39972144) returned first-author Darwin W. Kwok; Zotero search by `Kwok` returned the existing entry already in the collection. The correction triggered the immediate news_log + Issue body fixes today.

**Placement.** New 3rd subsection `### Kwok et al.: public neoepitopes from recurrent splicing — the off-the-shelf end of the axis` inserted into [DISCUSSIONS.md](research/manuscript/DISCUSSIONS.md) `## Comparison to related neoantigen-prediction tools` after the existing SpliceMutr (line 557) and ENEO (line 585) subsections. Three paragraphs: (1) what Kwok et al. did — recurrent *GNAS* / *RPL22* splicing across glioma / mesothelioma / prostate / liver, spatially conserved *GNAS* neojunction across multi-region biopsies, TCRs isolated from healthy donors with HLA-dependent peptide-dose-dependent killing of endogenously expressing tumor cell lines (Fig 5, Extended Data Fig 9-10); (2) public-vs-personalised axis — Kwok-class for off-the-shelf shared TCR targets across patients with recurrent NJ + matching HLA, our pipeline for patient-private candidates covering the residual majority; (3) GTEx threshold tradeoff — Kwok `PSR_GTEx < 1%` vs our `min_read_count: 1` (= `PSR_GTEx = 0%`) at the population level, with verifiability of NeoA<sub>GNAS</sub> / NeoA<sub>RPL22</sub> specifically deferred to follow-up [Issue #299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/299).

**Threshold-tradeoff framing — option (a), population-level only.** The user asked the sharper question — *"Were the GNAS or RPL22 neoepitopes derived from splice junctions with non-zero (but <1%) GTEx normals? Because only then can we maybe start to justify the <1% GTEx normals threshold, right?"* — and we exhausted public sources looking for the answer. None of the published material (main text, Extended Data figures 1-10, Springer supplementary MOESM1-5, [`dakwok/SSNIP`](https://github.com/dakwok/SSNIP) GitHub repo, peer review rebuttal MOESM5) reports the specific PSR<sub>GTEx</sub> values for the validated NEJ<sub>GNAS</sub> and NEJ<sub>RPL22</sub>. Extended Data Fig. 1D scatter plots show the public NJ population (colored dots) clustered near `PSR_GTEx = 0` but at the 0–1 axis scale visually indistinguishable from `PSR_GTEx ≈ 0.005`. The DISCUSSION subsection therefore states the tradeoff at the population level (*"at least some of Kwok et al.'s ~789 public NJ pool — which by their inclusion criterion spans `0% ≤ PSR_GTEx < 1%` — would be excluded by this pipeline's filter"*) without claiming the specific GNAS/RPL22 verdict. Verification deferred to follow-up [Issue #299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/299) (lightweight Snaptron query against GTEx v8, ~1–2 hours focused work, no pipeline rerun needed).

**Tooling note.** Added `openpyxl>=3.1` to [research/requirements.txt](research/requirements.txt) for the supplementary-table scan (verifying MOESM1/3/4 contents — none had the per-NJ table). Rolled into this PR.

**HTML-comment TODO marker.** Added `<!-- TODO(#299): re-derive PSR_GTEx ... -->` after the new subsection in [DISCUSSIONS.md](research/manuscript/DISCUSSIONS.md) so [Issue #299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/299)'s outcome has a clear insertion point in the manuscript prose without cluttering the rendered output.

---

## 2026-05-06

### 13:30 UTC — Editor: Scientist

#### [Sub-Issue #268](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/268) INTRODUCTION — Rojas 2023 foundational personalized neoantigen vaccine reference

**Placement.** Single sentence appended to the existing clinical-translation context block in section `## Cancer Neoepitopes` of [INTRODUCTION.md](research/manuscript/INTRODUCTION.md), inserted after the line "...personalized cancer vaccines and adoptive T-cell therapies." That line was the natural anchor — the sole pre-existing clinical-translation statement in INTRODUCTION, previously uncited. Considered Section 4 (HLA Genotype) where Sahin et al. 2017 + Ott et al. 2017 already cluster, and rejected — those are cited there for the *vaccine slot count* argument (10–20 candidates per formulation), a different rhetorical use; piling Rojas in would muddle the load-bearing function of the Section 4 citation. Considered a new paragraph and rejected — one foundational reference doesn't justify its own paragraph; INTRODUCTION should anchor and DISCUSSION should develop.

**Framing.** Rojas et al., *Nature* 2023, 618:144-150 (Zotero `UIN9DIUP`) — Phase 1 autogene cevumeran trial in PDAC, 8/16 vaccine-induced T-cell responders, mRNA-LNP synthesis with up to 20 patient-specific neoantigens. The cited fact ("8 of 16 patients developed vaccine-induced T-cell responses") is the single cleanest summary of the trial's primary endpoint; no need to also cite the recurrence-free survival trend (covered in DISCUSSION via [PR #260](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/260)) or the 6-year follow-up (AACR 2026, surfaced [2026-05-03 14:07 UTC] in news_log; deferred to a future DISCUSSION update if/when published as a full paper).

**Cross-section consistency.** INTRODUCTION now cites Rojas 2023 as foundational; DISCUSSION already develops the autogene cevumeran clinical context and adds the 6-year/long-lived T-cell durability follow-ups. Sahin et al. 2017 + Ott et al. 2017 remain in Section 4 only, doing different work (vaccine slot count). No duplication.

**Citation style.** Author-year inline `(Rojas et al., *Nature* 2023)` matches the existing INTRODUCTION convention used by Yewdell & Bennink 1999, Sahin et al. 2017, Ott et al. 2017, Cai et al. 2026 (added in [PR #282](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/282)), Alam et al. 2023. Reference-list finalisation deferred to [Sub-Issue #272](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/272).

---
