# Lab Notebook — Scientist

Per-role lab notebook for Scientist sessions. Started 2026-05-06 as part of the lab-notebook split (see `lab_notebook.md` freeze note for rationale).

Format and rules unchanged from the unified notebook — see `shared/feedback_lab_notebook.md`. New date sections at the TOP; new time sections at the TOP of the date block; entries are immutable once committed.

---

## 2026-05-06

### 13:30 UTC — Editor: Scientist

#### [Sub-Issue #268](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/268) INTRODUCTION — Rojas 2023 foundational personalized neoantigen vaccine reference

**Placement.** Single sentence appended to the existing clinical-translation context block in section `## Cancer Neoepitopes` of [INTRODUCTION.md](research/manuscript/INTRODUCTION.md), inserted after the line "...personalized cancer vaccines and adoptive T-cell therapies." That line was the natural anchor — the sole pre-existing clinical-translation statement in INTRODUCTION, previously uncited. Considered Section 4 (HLA Genotype) where Sahin et al. 2017 + Ott et al. 2017 already cluster, and rejected — those are cited there for the *vaccine slot count* argument (10–20 candidates per formulation), a different rhetorical use; piling Rojas in would muddle the load-bearing function of the Section 4 citation. Considered a new paragraph and rejected — one foundational reference doesn't justify its own paragraph; INTRODUCTION should anchor and DISCUSSION should develop.

**Framing.** Rojas et al., *Nature* 2023, 618:144-150 (Zotero `UIN9DIUP`) — Phase 1 autogene cevumeran trial in PDAC, 8/16 vaccine-induced T-cell responders, mRNA-LNP synthesis with up to 20 patient-specific neoantigens. The cited fact ("8 of 16 patients developed vaccine-induced T-cell responses") is the single cleanest summary of the trial's primary endpoint; no need to also cite the recurrence-free survival trend (covered in DISCUSSION via [PR #260](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/260)) or the 6-year follow-up (AACR 2026, surfaced [2026-05-03 14:07 UTC] in news_log; deferred to a future DISCUSSION update if/when published as a full paper).

**Cross-section consistency.** INTRODUCTION now cites Rojas 2023 as foundational; DISCUSSION already develops the autogene cevumeran clinical context and adds the 6-year/long-lived T-cell durability follow-ups. Sahin et al. 2017 + Ott et al. 2017 remain in Section 4 only, doing different work (vaccine slot count). No duplication.

**Citation style.** Author-year inline `(Rojas et al., *Nature* 2023)` matches the existing INTRODUCTION convention used by Yewdell & Bennink 1999, Sahin et al. 2017, Ott et al. 2017, Cai et al. 2026 (added in [PR #282](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/282)), Alam et al. 2023. Reference-list finalisation deferred to [Sub-Issue #272](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/272).

---
