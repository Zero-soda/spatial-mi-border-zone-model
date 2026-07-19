#!/usr/bin/env python3
"""Sparse and ranking utilities for spatial therapeutic prioritization."""

from __future__ import annotations

import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path

from project_paths import project_root
from typing import Iterable, Mapping, Sequence


ROOT = project_root(__file__)
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))

import h5py  # noqa: E402
import numpy as np  # noqa: E402


@dataclass(frozen=True)
class SparseMatrixData:
    """Gene-by-spot CSC matrix with 10x-compatible metadata."""

    data: np.ndarray
    indices: np.ndarray
    indptr: np.ndarray
    shape: tuple[int, int]
    gene_names: list[str]
    barcodes: list[str]
    gene_ids: list[str]

    @classmethod
    def from_dense(
        cls,
        dense: np.ndarray,
        gene_names: Sequence[str] | None = None,
        barcodes: Sequence[str] | None = None,
        gene_ids: Sequence[str] | None = None,
    ) -> "SparseMatrixData":
        arr = np.asarray(dense, dtype=float)
        if arr.ndim != 2:
            raise ValueError("dense matrix must be two-dimensional")
        n_genes, n_spots = arr.shape
        values: list[float] = []
        row_indices: list[int] = []
        indptr = [0]
        for spot_idx in range(n_spots):
            rows = np.flatnonzero(arr[:, spot_idx])
            values.extend(arr[rows, spot_idx].tolist())
            row_indices.extend(rows.tolist())
            indptr.append(len(values))
        names = list(gene_names or [f"gene_{idx}" for idx in range(n_genes)])
        ids = list(gene_ids or names)
        spots = list(barcodes or [f"spot_{idx}" for idx in range(n_spots)])
        if len(names) != n_genes or len(ids) != n_genes or len(spots) != n_spots:
            raise ValueError("metadata dimensions do not match dense matrix")
        return cls(
            data=np.asarray(values, dtype=float),
            indices=np.asarray(row_indices, dtype=np.int64),
            indptr=np.asarray(indptr, dtype=np.int64),
            shape=(n_genes, n_spots),
            gene_names=names,
            barcodes=spots,
            gene_ids=ids,
        )


def _decode(value: object) -> str:
    return value.decode() if isinstance(value, bytes) else str(value)


def load_10x_h5(path: Path) -> SparseMatrixData:
    """Load a filtered 10x matrix without densifying it."""

    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        features = matrix["features"]
        ids_key = "id" if "id" in features else "name"
        shape = tuple(int(value) for value in matrix["shape"][:])
        return SparseMatrixData(
            data=matrix["data"][:].astype(float, copy=False),
            indices=matrix["indices"][:].astype(np.int64, copy=False),
            indptr=matrix["indptr"][:].astype(np.int64, copy=False),
            shape=(shape[0], shape[1]),
            gene_names=[_decode(value) for value in features["name"][:]],
            barcodes=[_decode(value) for value in matrix["barcodes"][:]],
            gene_ids=[_decode(value) for value in features[ids_key][:]],
        )


def normalize_csc_log1p(matrix: SparseMatrixData, scale_factor: float = 10000.0) -> SparseMatrixData:
    """Return log-normalized CSC values while preserving zeros and metadata."""

    data = matrix.data.astype(float, copy=True)
    for spot_idx in range(matrix.shape[1]):
        start = int(matrix.indptr[spot_idx])
        end = int(matrix.indptr[spot_idx + 1])
        total = float(np.sum(data[start:end]))
        if total > 0:
            data[start:end] = np.log1p(data[start:end] / total * scale_factor)
    return SparseMatrixData(
        data=data,
        indices=matrix.indices,
        indptr=matrix.indptr,
        shape=matrix.shape,
        gene_names=matrix.gene_names,
        barcodes=matrix.barcodes,
        gene_ids=matrix.gene_ids,
    )


def subset_spots(matrix: SparseMatrixData, spot_indices: Sequence[int]) -> SparseMatrixData:
    """Select and reorder CSC columns."""

    values: list[float] = []
    rows: list[int] = []
    indptr = [0]
    selected_barcodes: list[str] = []
    for spot_idx in spot_indices:
        start = int(matrix.indptr[spot_idx])
        end = int(matrix.indptr[spot_idx + 1])
        values.extend(matrix.data[start:end].tolist())
        rows.extend(matrix.indices[start:end].tolist())
        indptr.append(len(values))
        selected_barcodes.append(matrix.barcodes[spot_idx])
    return SparseMatrixData(
        data=np.asarray(values, dtype=float),
        indices=np.asarray(rows, dtype=np.int64),
        indptr=np.asarray(indptr, dtype=np.int64),
        shape=(matrix.shape[0], len(spot_indices)),
        gene_names=matrix.gene_names,
        barcodes=selected_barcodes,
        gene_ids=matrix.gene_ids,
    )


def weighted_gene_sum(matrix: SparseMatrixData, weights: np.ndarray) -> np.ndarray:
    """Calculate one weighted sum per gene using sparse CSC columns."""

    weights = np.asarray(weights, dtype=float)
    if len(weights) != matrix.shape[1]:
        raise ValueError("weights length does not match number of spots")
    result = np.zeros(matrix.shape[0], dtype=float)
    for spot_idx, weight in enumerate(weights):
        if not math.isfinite(float(weight)) or weight == 0:
            continue
        start = int(matrix.indptr[spot_idx])
        end = int(matrix.indptr[spot_idx + 1])
        np.add.at(result, matrix.indices[start:end], matrix.data[start:end] * weight)
    return result


def weighted_gene_mean(matrix: SparseMatrixData, weights: np.ndarray) -> np.ndarray:
    """Calculate weighted gene means; weights need not be normalized."""

    weights = np.asarray(weights, dtype=float)
    finite = np.isfinite(weights)
    total = float(np.sum(weights[finite]))
    if total == 0:
        return np.full(matrix.shape[0], np.nan, dtype=float)
    cleaned = np.where(finite, weights, 0.0)
    return weighted_gene_sum(matrix, cleaned) / total


def gene_detection_fraction(matrix: SparseMatrixData, mask: np.ndarray | None = None) -> np.ndarray:
    """Fraction of selected spots with a non-zero value for every gene."""

    selected = np.ones(matrix.shape[1], dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    if len(selected) != matrix.shape[1] or not np.any(selected):
        return np.full(matrix.shape[0], np.nan, dtype=float)
    counts = np.zeros(matrix.shape[0], dtype=float)
    for spot_idx in np.flatnonzero(selected):
        start = int(matrix.indptr[spot_idx])
        end = int(matrix.indptr[spot_idx + 1])
        np.add.at(counts, matrix.indices[start:end], 1.0)
    return counts / float(np.sum(selected))


def gene_pearson(matrix: SparseMatrixData, values: np.ndarray) -> np.ndarray:
    """Pearson correlation of every gene with a spot-level vector."""

    y = np.asarray(values, dtype=float)
    valid = np.isfinite(y)
    n = int(np.sum(valid))
    if len(y) != matrix.shape[1] or n < 3:
        return np.full(matrix.shape[0], np.nan, dtype=float)
    sum_x = np.zeros(matrix.shape[0], dtype=float)
    sum_x2 = np.zeros(matrix.shape[0], dtype=float)
    sum_xy = np.zeros(matrix.shape[0], dtype=float)
    for spot_idx in np.flatnonzero(valid):
        start = int(matrix.indptr[spot_idx])
        end = int(matrix.indptr[spot_idx + 1])
        rows = matrix.indices[start:end]
        data = matrix.data[start:end]
        np.add.at(sum_x, rows, data)
        np.add.at(sum_x2, rows, data * data)
        np.add.at(sum_xy, rows, data * y[spot_idx])
    sum_y = float(np.sum(y[valid]))
    sum_y2 = float(np.sum(y[valid] * y[valid]))
    numerator = sum_xy - sum_x * sum_y / n
    x_ss = sum_x2 - sum_x * sum_x / n
    y_ss = sum_y2 - sum_y * sum_y / n
    denom = np.sqrt(np.maximum(x_ss, 0.0) * max(y_ss, 0.0))
    result = np.full(matrix.shape[0], np.nan, dtype=float)
    np.divide(numerator, denom, out=result, where=denom > 0)
    return result


def gene_linear_slope(matrix: SparseMatrixData, distance: np.ndarray) -> np.ndarray:
    """OLS slope for every gene against a shared signed-distance vector."""

    x = np.asarray(distance, dtype=float)
    valid = np.isfinite(x)
    if len(x) != matrix.shape[1] or np.sum(valid) < 3:
        return np.full(matrix.shape[0], np.nan, dtype=float)
    centered = np.zeros_like(x)
    centered[valid] = x[valid] - float(np.mean(x[valid]))
    denom = float(np.sum(centered[valid] ** 2))
    if denom == 0:
        return np.full(matrix.shape[0], np.nan, dtype=float)
    return weighted_gene_sum(matrix, centered) / denom


def extract_gene_rows(matrix: SparseMatrixData, gene_indices: Sequence[int]) -> np.ndarray:
    """Densify only selected gene rows for plotting or candidate-level tests."""

    lookup = {int(gene_idx): row_idx for row_idx, gene_idx in enumerate(gene_indices)}
    result = np.zeros((len(gene_indices), matrix.shape[1]), dtype=float)
    for spot_idx in range(matrix.shape[1]):
        start = int(matrix.indptr[spot_idx])
        end = int(matrix.indptr[spot_idx + 1])
        for gene_idx, value in zip(matrix.indices[start:end], matrix.data[start:end], strict=True):
            row_idx = lookup.get(int(gene_idx))
            if row_idx is not None:
                result[row_idx, spot_idx] = value
    return result


def build_knn_edges(coords: np.ndarray, k: int = 6) -> list[tuple[int, int]]:
    """Build deterministic symmetrized k-nearest-neighbour edges."""

    coords = np.asarray(coords, dtype=float)
    n = len(coords)
    if n < 2:
        return []
    delta = coords[:, None, :] - coords[None, :, :]
    distances = np.sum(delta * delta, axis=2)
    np.fill_diagonal(distances, np.inf)
    neighbours = np.argsort(distances, axis=1)[:, : min(k, n - 1)]
    return sorted({tuple(sorted((idx, int(other)))) for idx in range(n) for other in neighbours[idx]})


def graph_morans_i(values: np.ndarray, edges: Sequence[tuple[int, int]]) -> float:
    """Moran's I on an undirected graph with binary symmetric weights."""

    arr = np.asarray(values, dtype=float)
    if len(arr) < 3 or not edges or not np.all(np.isfinite(arr)):
        return float("nan")
    centered = arr - float(np.mean(arr))
    denom = float(np.sum(centered * centered))
    if denom == 0:
        return float("nan")
    numerator = float(sum(centered[left] * centered[right] for left, right in edges))
    return float((len(arr) / len(edges)) * numerator / denom)


def filter_one_to_one(records: Iterable[Mapping[str, object]]) -> dict[str, str]:
    """Return confident one-to-one mouse-to-human symbol mappings."""

    result: dict[str, str] = {}
    for record in records:
        if record.get("type") != "ortholog_one2one":
            continue
        try:
            confidence = float(record.get("confidence", 0))
        except (TypeError, ValueError):
            confidence = 0.0
        mouse = str(record.get("mouse", "")).strip()
        human = str(record.get("human", "")).strip()
        if confidence >= 1 and mouse and human:
            result[mouse] = human
    return result


INHIBIT_ACTIONS = {
    "antagonist",
    "blocker",
    "degrader",
    "inhibitor",
    "negative allosteric modulator",
    "suppressor",
}
ACTIVATE_ACTIONS = {
    "activator",
    "agonist",
    "positive allosteric modulator",
    "stimulator",
}


def normalize_action(action: object) -> str:
    text = str(action or "").strip().lower().replace("_", " ")
    if any(term in text for term in INHIBIT_ACTIONS):
        return "inhibit"
    if any(term in text for term in ACTIVATE_ACTIONS):
        return "activate"
    return "unknown"


def independent_source_count(rows: Iterable[Mapping[str, object]]) -> int:
    return len({str(row.get("source", "")).strip() for row in rows if str(row.get("source", "")).strip()})


def assign_tier(row: Mapping[str, object]) -> str:
    mouse = float(row.get("mouse_spatial_score", 0) or 0)
    human = float(row.get("human_support_score", 0) or 0)
    pharmacology = float(row.get("pharmacology_score", 0) or 0)
    penalty = float(row.get("repair_safety_penalty", 0) or 0)
    final_score = mouse + human + pharmacology - penalty
    direction = str(row.get("direction_match", "unknown"))
    source_count = int(float(row.get("pharmacology_source_count", 0) or 0))
    replicate_pass = bool(row.get("replicate_rule_pass", False))
    human_detected = bool(row.get("human_detected", False))
    if penalty >= 12 or direction == "protective_conflict":
        return "caution_protective"
    if (
        final_score >= 65
        and direction == "matched"
        and source_count >= 2
        and replicate_pass
        and human_detected
    ):
        return "tier_a"
    if final_score >= 45:
        return "tier_b"
    return "tier_c"


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(fieldnames), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
