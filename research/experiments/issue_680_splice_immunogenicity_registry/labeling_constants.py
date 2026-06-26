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
