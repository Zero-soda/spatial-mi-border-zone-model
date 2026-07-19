from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from sstba.scoring import score_modules


def test_zscore_module_means_are_frozen_and_oriented() -> None:
    expression = pd.DataFrame(
        {
            "G1": [1.0, 2.0, 3.0],
            "G2": [3.0, 2.0, 1.0],
            "G3": [0.0, 1.0, 4.0],
        },
        index=["s1", "s2", "s3"],
    )
    signatures = {"mechanical": ["G1", "G3"], "scar": ["G2"]}

    scores, coverage = score_modules(expression, signatures)

    assert list(scores) == ["mechanical", "scar"]
    assert coverage["mechanical"]["detected"] == 2
    assert coverage["mechanical"]["missing_genes"] == []
    assert np.isclose(scores["mechanical"].mean(), 0.0)
    assert scores.loc["s3", "mechanical"] > scores.loc["s1", "mechanical"]


def test_gene_symbols_are_case_insensitive_and_duplicates_are_rejected() -> None:
    expression = pd.DataFrame(
        {"Postn": [1.0, 3.0], "Cthrc1": [0.0, 4.0]},
        index=["s1", "s2"],
    )
    scores, coverage = score_modules(
        expression,
        {"scar": ["POSTN", "cthrc1", "POSTN"]},
    )

    assert coverage["scar"]["requested"] == 2
    assert coverage["scar"]["detected_genes"] == ["CTHRC1", "POSTN"]
    assert scores.loc["s2", "scar"] > scores.loc["s1", "scar"]


def test_zero_variance_gene_contributes_zero() -> None:
    expression = pd.DataFrame(
        {"G1": [2.0, 2.0, 2.0], "G2": [1.0, 2.0, 3.0]},
        index=["s1", "s2", "s3"],
    )

    scores, _ = score_modules(expression, {"state": ["G1", "G2"]})

    assert np.isfinite(scores.to_numpy()).all()
    assert np.isclose(scores["state"].mean(), 0.0)


def test_missing_module_fails_with_module_name() -> None:
    expression = pd.DataFrame({"G1": [1.0, 2.0]}, index=["s1", "s2"])

    with pytest.raises(ValueError, match="scar"):
        score_modules(expression, {"scar": ["POSTN", "CTHRC1"]})


def test_duplicate_spot_identifiers_are_rejected() -> None:
    expression = pd.DataFrame({"G1": [1.0, 2.0]}, index=["s1", "s1"])

    with pytest.raises(ValueError, match="Duplicate spot"):
        score_modules(expression, {"state": ["G1"]})
