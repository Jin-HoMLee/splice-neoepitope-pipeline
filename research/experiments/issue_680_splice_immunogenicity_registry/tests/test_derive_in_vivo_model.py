"""Unit tests for the in_vivo_model derivation (#1120).

Every other in-vivo test goes through `validate_registry.violations()` or pins the live
snapshot, which left the **no-guess `unspecified` branch untested by anything**: it is
reached only when a row has an in-vivo marker but its `source` is not in
`MODEL_BY_SOURCE_SUBSTR`, and no live row is in that state. That branch *is* the PR's
core principle - never infer a model type from a keyword - so leaving it unexercised
meant the one rule most likely to be quietly regressed had no guard at all.
(PR #1186 bot review, finding 1.)
"""
import pandas as pd
import pytest

from derive_in_vivo_model import has_in_vivo_readout, in_vivo_model


def _row(readout, source="Xiong 2025 (GBM)"):
    return pd.Series({"readout": readout, "source": source})


# --- has_in_vivo_readout: the gate the whole column hangs off ----------------------


@pytest.mark.parametrize("readout", [
    "IFN-g + LDH killing + in vivo + TCR-T",
    "engineered-TCR IFN-g + in-vivo tumor control",   # hyphenated spelling
    "IFN-g + IN VIVO",                                # case-insensitive
])
def test_marker_is_detected(readout):
    assert has_in_vivo_readout(readout)


@pytest.mark.parametrize("readout", [
    "IFN-g (41% CD8+) + MS",
    "tetramer+ ex-vivo",
    "",
])
def test_absent_marker_is_not_detected(readout):
    assert not has_in_vivo_readout(readout)


# --- in_vivo_model: the three branches --------------------------------------------


def test_known_source_with_an_in_vivo_readout_derives_its_model():
    """Verified first-hand (NSG xenograft), so it is source-keyed, not inferred."""
    assert in_vivo_model(_row("IFN-g + in vivo", "Xiong 2025 (GBM)")) == "xenograft"
    assert in_vivo_model(_row("IFN-g + in vivo", "Kim GB 2022 STM (COL6A3 tumor-stroma)")) == "xenograft"


def test_unknown_source_with_an_in_vivo_readout_is_unspecified_not_guessed():
    """**The no-guess branch, and the reason this file exists.**

    A future source reports an in-vivo experiment but we have not established the model
    first-hand. Xenograft is the *likely* answer and that is exactly why it must not be
    the *derived* one: xenograft-vs-syngeneic is a claim about whether an intact immune
    system was involved, and getting it wrong silently inverts what the row means.
    """
    assert in_vivo_model(_row("IFN-g + cytotoxicity + in vivo", "Some New Study 2027")) == "unspecified"


def test_no_marker_means_no_model_even_for_an_in_vivo_source():
    """The marker-first gate. Keying on `source` alone would stamp a model onto every
    row of an in-vivo source - including Xiong's constitutive `VFVDGLCRAKF` control,
    which was never put in a mouse."""
    assert in_vivo_model(_row("RCAN1-1 constitutive control", "Xiong 2025 (GBM)")) == "none"
    assert in_vivo_model(_row("IFN-g (41% CD8+) + MS", "SNAF (Li 2024, Sci Transl Med)")) == "none"
