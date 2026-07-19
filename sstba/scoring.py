"""Frozen gene-module scoring for spatial expression matrices."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd


def _normalise_gene_columns(expression: pd.DataFrame) -> pd.DataFrame:
    if expression.index.has_duplicates:
        duplicates = expression.index[expression.index.duplicated()].unique()
        preview = ", ".join(map(str, duplicates[:3]))
        raise ValueError(f"Duplicate spot identifiers are not allowed: {preview}")
    if expression.empty:
        raise ValueError("Expression matrix is empty")

    numeric = expression.apply(pd.to_numeric, errors="coerce")
    numeric.columns = [str(column).strip().upper() for column in numeric.columns]
    if any(not column for column in numeric.columns):
        raise ValueError("Expression matrix contains an empty gene symbol")
    if numeric.columns.has_duplicates:
        numeric = numeric.T.groupby(level=0, sort=False).mean().T
    return numeric


def _gene_zscores(expression: pd.DataFrame) -> pd.DataFrame:
    means = expression.mean(axis=0, skipna=True)
    standard_deviations = expression.std(axis=0, ddof=0, skipna=True)
    safe_standard_deviations = standard_deviations.mask(
        (~np.isfinite(standard_deviations)) | (standard_deviations == 0),
        1.0,
    )
    zscores = expression.subtract(means, axis=1).divide(
        safe_standard_deviations,
        axis=1,
    )
    zero_variance = standard_deviations.index[
        (~np.isfinite(standard_deviations)) | (standard_deviations == 0)
    ]
    if len(zero_variance):
        zscores.loc[:, zero_variance] = 0.0
    return zscores.fillna(0.0)


def score_modules(
    expression: pd.DataFrame,
    signatures: Mapping[str, Sequence[str]],
    *,
    min_detected_genes: int = 1,
) -> tuple[pd.DataFrame, dict[str, dict[str, object]]]:
    """Calculate mean per-gene z scores using frozen module definitions.

    Parameters
    ----------
    expression
        Spots by genes expression matrix.
    signatures
        Mapping from output module name to requested gene symbols.
    min_detected_genes
        Minimum number of genes that must be available for each module.
    """

    if min_detected_genes < 1:
        raise ValueError("min_detected_genes must be at least 1")
    if not signatures:
        raise ValueError("At least one signature module is required")

    numeric = _normalise_gene_columns(expression)
    zscores = _gene_zscores(numeric)
    available = set(zscores.columns)
    score_columns: dict[str, pd.Series] = {}
    coverage: dict[str, dict[str, object]] = {}

    for raw_name, raw_genes in signatures.items():
        name = str(raw_name).strip()
        if not name:
            raise ValueError("Signature module names cannot be empty")
        requested_genes = sorted(
            {
                str(gene).strip().upper()
                for gene in raw_genes
                if str(gene).strip()
            }
        )
        detected_genes = sorted(set(requested_genes) & available)
        missing_genes = sorted(set(requested_genes) - available)
        coverage[name] = {
            "requested": len(requested_genes),
            "detected": len(detected_genes),
            "detected_genes": detected_genes,
            "missing_genes": missing_genes,
        }
        if len(detected_genes) < min_detected_genes:
            raise ValueError(
                f"Signature module '{name}' detected {len(detected_genes)} "
                f"gene(s); at least {min_detected_genes} required"
            )
        score_columns[name] = zscores.loc[:, detected_genes].mean(axis=1)

    return pd.DataFrame(score_columns, index=numeric.index), coverage
