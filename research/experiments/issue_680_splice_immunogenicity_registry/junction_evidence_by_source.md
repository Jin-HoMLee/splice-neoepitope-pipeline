# Junction evidence by source - splice immunogenicity registry (#735)

For each of the 10 sources in `registry.tsv`, this file records what junction unit the source
published for its peptides, the resulting grade on the ladder, and the PROVENANCE location that
backs each call.

No-inference rule: grades are assigned only from evidence explicitly quoted in PROVENANCE.md.
Coordinates or identifiers not recorded in PROVENANCE are not inferred from peptide sequence,
from the source type, or from knowledge that the source is "junction-based" in principle.

---

## Grade ladder (summary)

| Grade | Assigned when... |
|---|---|
| `coords` | Source published a genomic donor/acceptor coordinate (chr:start-end) quoted in PROVENANCE |
| `event-id` | Source published a transcript or splice-event identifier (not full coords) quoted in PROVENANCE |
| `gene-mechanism` | Only gene name + canonical splice mechanism available; no junction-level identifier |
| `none` | No recoverable junction; reason recorded in `notes` column |

---

## 1. SNAF - Li et al. 2024, Sci Transl Med (eade2886)

**Junction unit published:** SNAF Data S1 contains per-peptide junction records with two levels of
detail: (a) a SNAF event notation in the format `GENE:EXON_POSITION-EXON` (one endpoint embedded in
an event label) and (b) a full genomic coordinate pair `chr:start-end(strand)`.

**Grade split (3 distinct sub-grades across 6 rows):**

| Peptide | Grade | Junction evidence | PROVENANCE location |
|---|---|---|---|
| `IPDSQGNDI` | `coords` | `chr5:33954504-33963931(-)` | PROVENANCE: "Data S1-TCGA sheet ... chr5:33954504-33963931(-)" |
| `STESITATL` | `coords` | `chr12:55956189-55956949(-)` | PROVENANCE: "Data S1 row ... chr12:55956189-55956949(-)" |
| `IIDNQEPVF` | `event-id` | `CDH19:E12.1-E13.2_66509195` | PROVENANCE: "sm.pdf Supp Fig7B, verbatim 'CDH19:E12.1-E13.2_66509195'" |
| `RLLGTEFQT` | `gene-mechanism` | SLC45A2, alt-5' splice site | PROVENANCE: "Fig5A legend, verbatim 'RLLGTEFQT (RLL)'" - no coord/event-id quoted |
| `FQTRRAMTL` | `gene-mechanism` | SLC45A2, alt-5' splice site | PROVENANCE: "Fig5A/B legend" only - no coord/event-id quoted |
| `TEFQTRRAM` | `gene-mechanism` | SLC45A2, alt-5' splice site | PROVENANCE: "main text + sm.pdf Supp Fig5B" only - no coord/event-id quoted |

**Rationale for IIDNQEPVF = event-id (not coords):** The SNAF event notation
`CDH19:E12.1-E13.2_66509195` encodes a single junction endpoint (66509195) within a
gene:exon-pair label. A full genomic coordinate requires both donor and acceptor positions
(as shown for IPDSQGNDI and STESITATL with their two-coordinate `chr:start-end` form).
Since PROVENANCE quotes only the one-coordinate SNAF event label - not a two-coordinate
chr:start-end - `event-id` is the conservative assignment.

**Caveat on the 2 `coords` rows (`IPDSQGNDI`, `STESITATL`):** the `coords` grade is correct - the genomic coordinate pair was published in SNAF Data S1, so the no-inference rule holds for the junction annotation.
But both rows carry `provenance_grade=inferred` with an open `VERIFY` caveat on the **peptide-identity** match (the exact sequence was matched from a 3-letter figure prefix, `IPD` / `STE`, to the unique Data S1 row, not read verbatim).
The junction coordinate is solid; the peptide-to-coordinate linkage inherits that sequence-identity caveat, so do not over-trust `coords` as fully verified end-to-end until the peptide identity is confirmed against the synthesized-peptide list.

---

## 2. Kim et al. 2025 - mis-splicing neoantigens in SF-mutant leukemias

**Junction unit published:** Figure labels and MS spectra identifying peptides by gene and
splice mechanism (junction-derived, intron retention, poison exon, minor intron). No genomic
coordinates or junction event IDs are quoted in PROVENANCE for any row.

**Grade: `gene-mechanism` for all 5 rows.** Gene names and splice mechanisms are known from
figure annotations (e.g. "CLK3 #1", "RHOT2 #5", "C16orf70", "ATG3", "MYO1F"), but PROVENANCE
does not quote any junction-level identifier.

PROVENANCE location: Kim Fig S1B, Fig 2H, Fig 2C, Fig S2F, Fig S2H/I; Zotero LX6DMXTL;
all `agent-local`.

---

## 3. Xiong et al. 2025 - RCAN1-4, glioblastoma

**Junction unit published:** The RCAN1-4 isoform is a named transcript whose exon4/exon5
junction is the specific splice junction. PROVENANCE cites "SJ exon4/exon5" and Fig3A structure
showing the exon-exon junction (NSDIFSESETR-AKFESLFRTY boundary). The task brief explicitly
lists "RCAN1-4 exon4/exon5 SJ" as an `event-id` example.

**Grade split (1 event-id, 1 none):**

| Peptide | Grade | Junction evidence | PROVENANCE location |
|---|---|---|---|
| `IFSESETRAKF` | `event-id` | `RCAN1-4:exon4/exon5` | PROVENANCE: "Fig3A structure '...SJ exon4/exon5'"; Zotero JJZL6D84 |
| `VFVDGLCRAKF` | `none` | - | RCAN1-1_77-87 constitutive control - not a splice junction; fails splice gate |

**Rationale for VFVDGLCRAKF = none:** This is a constitutive exon peptide used as a same-locus
control. It is not splice-derived (`splice_mechanism_canonical = not_splice`), so junction grading
is not applicable - recorded as `none` with reason in notes.

---

## 4. Fisher et al. 2026 - CoREST

**Junction unit published:** None. PROVENANCE explicitly notes "genes not named in source". All
three peptides (RSQGWLFLR, AEHAHRVPL, VELEDHVML) are identified only from Fig7I/J/K peptide
tables and ELISpot panels, with no gene names and no junction coordinates or event IDs.

**Grade: `none` for all 3 rows.** Without a gene name, the originating junction is unrecoverable.
The reason is appended to `notes` for each row.

PROVENANCE location: Fisher Fig7I peptide table + Fig7J/K ELISpot, p15; Zotero WLSZLMA4;
all `agent-local`.

---

## 5. IRIS - Pan et al., PNAS

**Junction unit published:** PROVENANCE records gene names and the label "alt-splicing" for each
peptide, plus TCR clone identifiers (JPTCR-238, JPTCR-45/47/50/56, JPTCR-13, JPTCR-52). TCR
clone IDs are T-cell receptor identifiers, not splice junction event IDs. No IRIS junction event
IDs (the identifiers that IRIS as a platform would assign to specific intron-retention or
alternative-splicing events) are quoted in PROVENANCE.

**Grade: `gene-mechanism` for all 4 rows.** Gene names (CLASP1, SCAMP3, MAN2C1, AP2A1) and the
"alt-splicing" mechanism are known, but no junction event ID is quoted in PROVENANCE. Upgrading
to `event-id` would require a PROVENANCE-quoted IRIS event ID, which is not present (source is
`agent-web` with a VERIFY flag).

PROVENANCE location: IRIS Dataset S3c; NOT in Zotero - VERIFY in supplement; all `agent-web`.

---

## 6. POSTN-203 study

**Junction unit published:** PROVENANCE quotes the transcript identifier `ENST00000379747` for
POSTN isoform 203. This is a formal Ensembl transcript accession identifying the specific POSTN
splice isoform. The task brief lists "a POSTN ENST... transcript" as an `event-id` example.

**Grade: `event-id` for the 1 row.** junction_id = `ENST00000379747`.

PROVENANCE location: "POSTN-203 (ENST00000379747)"; agent-web; NOT in Zotero - VERIFY.

---

## 7. Kwok et al. 2024/2025 - Nature (s41586-024-08552-0)

**Junction unit published:** Gene names (GNAS, RPL22) and junction mechanisms (frameshift
neojunction, in-frame neojunction) are described in figures. PROVENANCE records that exact peptide
sequences are NOT published (figure-locked in MS peak plots or generic schematics). No genomic
coordinates or junction event IDs are quoted in PROVENANCE.

**Grade: `gene-mechanism` for both rows.** Gene and mechanism are known; the sequences are
`unpublished` (per provenance_grade column), but the junction type is described. Since genes ARE
named and junctions ARE characterized (recurrent public neojunctions), `gene-mechanism` is correct
rather than `none` (cf. Fisher rows where even the gene is unnamed).

PROVENANCE location: Kwok Fig4b-h, Fig5d, ExtData Fig6e/f, ExtData Fig9b; seq NOT legible;
Zotero EA2NMFNZ; `unpublished`.

---

## 8. Manoharan et al. 2026 - intron retention, colorectal cancer (IR-CRC)

**Junction unit published:** Table 1 lists gene names (IR1-IR11), HLA alleles, and identifies
each peptide as intron-retention derived. No genomic coordinates or intron-retention event IDs
(e.g. intron number, chromosome position) are quoted in PROVENANCE.

**Grade: `gene-mechanism` for all 11 rows.** Gene name + intron retention mechanism is the
maximum information available from PROVENANCE.

PROVENANCE location: Manoharan 2026 PMC13096189 Table 1 (IR1-IR11); 2 independent PMC reads;
Zotero ZAT8678F; all `agent-web`.

---

## 9. Merlotti et al. 2023 - exon-TE junctions, NSCLC

**Junction unit published:** Data file S7 "Characteristics of predicted pJETs used for in-vitro
assays" assigns each peptide a numeric junction ID (e.g. 3228-1, 321, 704). These are Merlotti's
pJET (predicted Junction Exon-TE) identifiers. The task brief explicitly lists "a Merlotti
pJET/TE-junction ID" as an `event-id` example.

**Grade: `event-id` for all 10 rows.** All 10 pJET IDs are quoted in PROVENANCE.

| Peptide | junction_id |
|---|---|
| LLGETKVYV | Merlotti2023-S7:3228-1 |
| LLDRFGYHV | Merlotti2023-S7:321 |
| YLWTTFFPL | Merlotti2023-S7:704 |
| SLMQSGSPV | Merlotti2023-S7:299 |
| AILPKANTV | Merlotti2023-S7:4111 |
| YLPYFLKSL | Merlotti2023-S7:2242 |
| ILSGYGPCV | Merlotti2023-S7:2930-1 |
| ILANLPPAL | Merlotti2023-S7:1825 |
| RLLHLESFL | Merlotti2023-S7:1375 |
| VLMWTMAHL | Merlotti2023-S7:4188 |

PROVENANCE location: Merlotti 2023 Data file S7 (ID NNN) for each; Zotero 5ZADNVCB (local supp
xlsx); all `direct`.

---

## 10. Bigot et al. 2021 - SF3B1-mutant uveal melanoma

**Grade split (event-id for 5 S8-validated; gene-mechanism for 30 tetramer-only):**

### 10a. S8-junction-validated (5 rows) - `event-id`

Supplementary Table S8 provides explicit alternative-mRNA-junction validation for 5 peptides:
it shows the aberrant splicing event creates a frameshift (+N nucleotides) with NMD prediction.
The task brief explicitly lists "a Bigot Table-S8 alternative-mRNA-junction" as an `event-id`
example. The A2:N panel peptide identifiers (from Table S2, cross-confirmed in PROVENANCE) are
the event IDs.

| Peptide | gene | junction_id | S8 evidence |
|---|---|---|---|
| LLIRWQHFL | NF1 | Bigot2021-S8:A2:14 | alt mRNA +14, NMD+ |
| AALPILFQV | ATP8B2 | Bigot2021-S8:A2:17 | alt mRNA +20, NMD+ |
| ALLLQLFTL | USP39 | Bigot2021-S8:A2:18 | alt mRNA +40, NMD+ |
| ALLPGLPAA | NET1 | Bigot2021-S8:A2:26 | alt mRNA +22 |
| RLPGVLPRA | MAPK8IP2 | Bigot2021-S8:A2:37 | alt mRNA +16 |

PROVENANCE location: "S2 panel (A2:N) seq + S8 junction + IEDB ref1040678 assays"; Zotero
DJAS2BJ2; `direct`.

### 10b. Tetramer-only (30 rows) - `gene-mechanism`

The other 30 IEDB-positive panel peptides have patient ex-vivo tetramer detection only. Their
splice origin rests on the A2:N panel design (all SF3B1-aberrant-splicing-predicted) and IEDB
source-gene curation. PROVENANCE notes explicitly: "no per-peptide effector/junction" and
"splice origin via panel design + IEDB gene" for these rows. No per-peptide junction event IDs
or coordinates are quoted.

PROVENANCE location: S2 panel (A2:N) seq + IEDB ref1040678 (tetramer); Zotero DJAS2BJ2; `direct`.

---

## Grade tally (all 79 rows)

| Grade | Count | Sources |
|---|---|---|
| `gene-mechanism` | 55 | SNAF x3, Kim x5, IRIS x4, Kwok x2, Manoharan x11, Bigot-tetramer x30 |
| `event-id` | 18 | SNAF-IIDNQEPVF x1, Xiong-RCAN1-4 x1, POSTN x1, Merlotti x10, Bigot-S8 x5 |
| `coords` | 2 | SNAF-IPDSQGNDI, SNAF-STESITATL |
| `none` | 4 | Xiong-RCAN1-1 (constitutive), Fisher x3 (gene unnamed) |
| **Total** | **79** | |

---

## Ambiguous calls - documented

1. **IIDNQEPVF (SNAF):** The SNAF event notation `CDH19:E12.1-E13.2_66509195` contains one
   genomic position (66509195). One could argue this partial coordinate means `coords`. Conservatively
   assigned `event-id` because a junction coordinate requires BOTH donor and acceptor positions;
   PROVENANCE does not quote the second position, and the two-coordinate form (`chr:start-end`) is
   what PROVENANCE explicitly quotes for the other two SNAF rows with genomic evidence.

2. **IRIS rows:** The IRIS paper is a junction-level study with formal IRIS event IDs. However,
   PROVENANCE records only gene names and TCR clone IDs for these 4 rows - not IRIS event IDs.
   Conservatively assigned `gene-mechanism` to reflect only what PROVENANCE quotes. These are
   `agent-web` with a VERIFY flag; upgrading to `event-id` would require a verified PROVENANCE
   re-read confirming the event ID.

3. **Kwok rows (GNAS, RPL22):** Sequences are unpublished. The junction type is described (public
   recurrent neojunction) and genes are named. Assigned `gene-mechanism` (not `none`) because
   gene+mechanism is recoverable even though sequences are not. If the junction had no gene
   name (like Fisher), it would be `none`.
