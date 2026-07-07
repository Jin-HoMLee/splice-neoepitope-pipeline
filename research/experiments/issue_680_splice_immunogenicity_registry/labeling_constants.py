#!/usr/bin/env python3
"""Shared controlled vocabularies + keyword rules for the #735 labeling scheme.

Single source of truth imported by both validate_registry.py and
derive_evidence_strength.py, so the derivation rule and the validator that
enforces it cannot drift. (Previously these were duplicated byte-for-byte in
both scripts and kept in sync only by a manual checklist item.)
"""

# Controlled vocabularies for the two #735 columns.
GRADES = {"coords", "event-id", "gene-mechanism", "none"}
STRENGTHS = {"strong", "weak", "hard", "soft", "na"}

# Controlled vocabulary for the #823 assay_context column. Captures *which
# immunological system* produced the functional readout, so a scoring run can
# weight rows by assay realism (patient ex-vivo detection > healthy-donor IVS
# priming > engineered-TCR readout). Assigned source-keyed from the curated
# rationale in PROVENANCE.md; `unspecified` is the honest value for rows whose
# held provenance records the assay but NOT its T-cell source (no guessing -
# verify-against-source rule). See derive_assay_context.py for the rule.
ASSAY_CONTEXTS = {
    "patient_exvivo",     # patient PBMC/blood ex-vivo tetramer+ or functional
    "patient_til",        # patient tumor-infiltrating (or draining-LN) lymphocytes
    "healthy_donor_ivs",  # healthy-donor in-vitro-sensitized (IVS) T cells
    "cloned_tcr",         # engineered/cloned-TCR functional readout (no primary patient/donor detection)
    "prevalence_only",    # population prevalence / presentation, no per-peptide T-cell assay
    "unspecified",        # functional assay reported but T-cell source not determinable from held provenance
    "na",                 # not applicable (constitutive control / non-tested sequence)
}

# Controlled vocabulary for the #1001 venue_type column. Records the publication
# venue class of each row's source so peer-reviewed evidence is queryable and
# auditable at a glance (GRADE / systematic-review practice down-weights
# non-peer-reviewed sources). Assigned source-keyed from the curated venue audit
# in PROVENANCE.md. `db_recovered` is reserved for a row recovered from a public
# deposit (IEDB / CEDAR) whose underlying publication venue is ambiguous - none
# exist today. There is deliberately NO in-vocabulary "unknown": a source not in
# the venue map derives to the out-of-vocab sentinel `unclassified`, which the
# validator rejects, so a newly-folded source cannot slip in venue-unmarked
# (esp. a preprint). See derive_venue_type.py for the rule.
VENUE_TYPES = {
    "journal",       # peer-reviewed journal article
    "preprint",      # preprint server (bioRxiv / medRxiv / arXiv), not yet peer-reviewed
    "db_recovered",  # recovered from a public deposit (IEDB/CEDAR); underlying venue ambiguous
}
VENUE_UNCLASSIFIED = "unclassified"  # out-of-vocab sentinel: forces classification of a new source

# Source-keyed venue map (#1001): matched on a lowercased substring of the
# `source` column (same convention as derive_assay_context.py). This is the
# single source of truth for both the derivation (derive_venue_type.py reads
# `venue`) and the offline-optional Zotero cross-check in validate_registry.py
# (which reads `doi`). All 15 current sources are peer-reviewed journals (14 per
# the 2026-07-04 audit + COL6A3/STM verified 2026-07-07); DOIs verified against
# Zotero collection Z38GTJNW.
VENUE_BY_SOURCE_SUBSTR = {
    "bigot":     {"venue": "journal", "doi": "10.1158/2159-8290.cd-20-0555"},   # Cancer Discovery
    "kim 2025":  {"venue": "journal", "doi": "10.1016/j.cell.2025.03.047"},     # Cell
    "manoharan": {"venue": "journal", "doi": "10.1038/s41598-026-43687-2"},     # Scientific Reports
    "merlotti":  {"venue": "journal", "doi": "10.1126/sciimmunol.abm6359"},     # Science Immunology
    "snaf":      {"venue": "journal", "doi": "10.1126/scitranslmed.ade2886"},   # Sci Transl Med
    "long-read": {"venue": "journal", "doi": "10.1158/2326-6066.cir-23-0083"},  # Cancer Immunology Research
    "iris":      {"venue": "journal", "doi": "10.1073/pnas.2221116120"},        # PNAS
    "fisher":    {"venue": "journal", "doi": "10.1172/jci.insight.190287"},     # JCI Insight
    "kwok":      {"venue": "journal", "doi": "10.1038/s41586-024-08552-0"},     # Nature
    "xiong":     {"venue": "journal", "doi": "10.1038/s41423-025-01360-0"},     # Cellular & Molecular Immunology (GBM)
    "postn":     {"venue": "journal", "doi": "10.1038/s41435-025-00326-6"},     # Genes & Immunity
    "col6a3":    {"venue": "journal", "doi": "10.1126/scitranslmed.abo6135"},   # Science Translational Medicine (Kim/Immatics tumor-stroma)
}

# Zotero collection that mirrors the registry's sources, and the itemType ->
# venue_type mapping the cross-check uses. A Zotero itemType absent here (e.g.
# conferencePaper) is surfaced by the cross-check, not silently accepted.
ZOTERO_COLLECTION = "Z38GTJNW"
ZOTERO_ITEMTYPE_TO_VENUE = {
    "journalArticle": "journal",
    "preprint": "preprint",
}

# Preprint-server markers (#1001 review finding 1). The source-keyed derivation
# matches on study substring, so it CANNOT self-distinguish a preprint from the
# journal version of an *already-mapped* study (a future "SNAF ... bioRxiv" source
# would first-match "snaf" -> journal). This keyword set closes that hole: the
# validator fails any row whose `source` names a preprint server but whose
# venue_type is not `preprint`. Matched as a lowercased substring of `source`
# only (not `notes`, which may legitimately mention a preprint for a journal row).
PREPRINT_MARKERS = ("biorxiv", "medrxiv", "arxiv", "preprint", "research square",
                    "ssrn", "chemrxiv", "osf.io")

# Effector vs detection-only keyword sets (the positive evidence_strength rule):
# a positive readout naming any EFFECTOR term is `strong`; one with only a
# DETECTION term is `weak`. Mass-spec (`MS`) is presentation, not T-cell
# function, so it is deliberately absent from EFFECTOR.
EFFECTOR = ("ifn", "elispot", "cytotox", "granzyme", "gzmb", "cd107", "cd137",
            "degranul", "killing", "ldh", "caspase", "incucyte", "tcr", "tnf",
            "activation", "in vivo")
DETECTION = ("tetramer", "dextramer", "multimer")

# IVS = in-vitro sensitization. A negative tested in healthy-donor IVS is a
# `soft` negative (failed to prime != intrinsically non-immunogenic); a
# negative with a functional assay on presented antigen and no IVS context is
# `hard`. Used by both the derivation rule and the validator cross-check.
IVS_MARKER = "ivs"
