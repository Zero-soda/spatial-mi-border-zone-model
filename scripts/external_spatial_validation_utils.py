#!/usr/bin/env python3
"""Shared readers and fixed-score helpers for external spatial MI validation."""

from __future__ import annotations

import csv
import gzip
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, TextIO

import h5py
import numpy as np

from therapeutic_prioritization_utils import (
    SparseMatrixData,
    normalize_csc_log1p,
    subset_spots,
)


REQUIRED_MANIFEST_FIELDS = (
    "dataset",
    "sample",
    "file_role",
    "url",
    "stage_region",
    "source_citation",
)


@dataclass(frozen=True)
class H5adSpatialData:
    matrix: SparseMatrixData
    coordinates: np.ndarray
    obs_columns: dict[str, list[str]]
    matrix_scale: str


def _decode(value: object) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def validate_manifest_rows(rows: Iterable[Mapping[str, str]]) -> None:
    seen: set[tuple[str, str, str]] = set()
    for row_number, row in enumerate(rows, start=2):
        missing = [field for field in REQUIRED_MANIFEST_FIELDS if not str(row.get(field, "")).strip()]
        if missing:
            raise ValueError(f"manifest row {row_number} missing fields: {', '.join(missing)}")
        if not str(row["url"]).startswith("https://"):
            raise ValueError(f"manifest row {row_number} URL must use HTTPS")
        key = (str(row["dataset"]), str(row["sample"]), str(row["file_role"]))
        if key in seen:
            raise ValueError(f"duplicate manifest key: {key}")
        seen.add(key)


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    validate_manifest_rows(rows)
    return rows


def read_10x_mtx_bundle(matrix_path: Path, features_path: Path, barcodes_path: Path) -> SparseMatrixData:
    gene_ids: list[str] = []
    gene_names: list[str] = []
    with _open_text(features_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if not fields or not fields[0]:
                continue
            gene_ids.append(fields[0])
            gene_names.append(fields[1] if len(fields) > 1 and fields[1] else fields[0])

    with _open_text(barcodes_path) as handle:
        barcodes = [line.strip().split("\t")[0] for line in handle if line.strip()]

    with _open_text(matrix_path) as handle:
        header = handle.readline().strip()
        if not header.startswith("%%MatrixMarket matrix coordinate"):
            raise ValueError(f"unsupported Matrix Market header: {header}")
        dimensions = ""
        for line in handle:
            if line.startswith("%") or not line.strip():
                continue
            dimensions = line
            break
        if not dimensions:
            raise ValueError("Matrix Market dimensions are missing")
        n_genes, n_spots, _ = (int(value) for value in dimensions.split()[:3])
        columns: list[dict[int, float]] = [dict() for _ in range(n_spots)]
        for line in handle:
            if not line.strip() or line.startswith("%"):
                continue
            fields = line.split()
            gene_idx = int(fields[0]) - 1
            spot_idx = int(fields[1]) - 1
            value = float(fields[2])
            if not 0 <= gene_idx < n_genes or not 0 <= spot_idx < n_spots:
                raise ValueError("Matrix Market coordinate exceeds declared shape")
            columns[spot_idx][gene_idx] = columns[spot_idx].get(gene_idx, 0.0) + value

    if n_genes != len(gene_names) or n_spots != len(barcodes):
        raise ValueError(
            f"matrix dimensions {(n_genes, n_spots)} do not match "
            f"features/barcodes {(len(gene_names), len(barcodes))}"
        )

    data: list[float] = []
    indices: list[int] = []
    indptr = [0]
    for column in columns:
        for gene_idx in sorted(column):
            indices.append(gene_idx)
            data.append(column[gene_idx])
        indptr.append(len(data))
    return SparseMatrixData(
        data=np.asarray(data, dtype=float),
        indices=np.asarray(indices, dtype=np.int64),
        indptr=np.asarray(indptr, dtype=np.int64),
        shape=(n_genes, n_spots),
        gene_names=gene_names,
        barcodes=barcodes,
        gene_ids=gene_ids,
    )


def read_tissue_positions(path: Path) -> dict[str, dict[str, float | int]]:
    fields = (
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres",
    )
    with _open_text(path) as handle:
        rows = list(csv.reader(handle))
    if not rows:
        return {}
    if rows[0] and rows[0][0].strip().lower() == "barcode":
        header = [value.strip() for value in rows.pop(0)]
    else:
        header = list(fields)
    result: dict[str, dict[str, float | int]] = {}
    for values in rows:
        if not values:
            continue
        row = dict(zip(header, values, strict=False))
        barcode = row.get("barcode", "").strip()
        if not barcode:
            continue
        result[barcode] = {
            "in_tissue": int(float(row.get("in_tissue", 1))),
            "array_row": int(float(row.get("array_row", 0))),
            "array_col": int(float(row.get("array_col", 0))),
            "pxl_row_in_fullres": float(row.get("pxl_row_in_fullres", math.nan)),
            "pxl_col_in_fullres": float(row.get("pxl_col_in_fullres", math.nan)),
        }
    return result


def _read_h5ad_strings(node: h5py.Dataset | h5py.Group) -> list[str]:
    if isinstance(node, h5py.Dataset):
        values = node[:]
        return [_decode(value) for value in values]
    if "codes" in node and "categories" in node:
        codes = node["codes"][:]
        categories = [_decode(value) for value in node["categories"][:]]
        return [categories[int(code)] if int(code) >= 0 else "" for code in codes]
    raise ValueError(f"unsupported h5ad column encoding at {node.name}")


def _h5ad_index(group: h5py.Group) -> list[str]:
    if "_index" in group:
        return _read_h5ad_strings(group["_index"])
    index_name = group.attrs.get("_index")
    if index_name is not None:
        name = _decode(index_name)
        if name in group:
            return _read_h5ad_strings(group[name])
    raise ValueError(f"h5ad group {group.name} has no readable index")


def _transpose_h5ad_sparse(
    data: np.ndarray,
    indices: np.ndarray,
    indptr: np.ndarray,
    obs_var_shape: tuple[int, int],
    encoding: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_obs, n_var = obs_var_shape
    if encoding == "csr_matrix":
        return data.astype(float, copy=False), indices.astype(np.int64, copy=False), indptr.astype(np.int64, copy=False)
    if encoding != "csc_matrix":
        raise ValueError(f"unsupported h5ad sparse encoding: {encoding}")
    spot_columns: list[dict[int, float]] = [dict() for _ in range(n_obs)]
    for gene_idx in range(n_var):
        start = int(indptr[gene_idx])
        end = int(indptr[gene_idx + 1])
        for spot_idx, value in zip(indices[start:end], data[start:end], strict=True):
            spot_columns[int(spot_idx)][gene_idx] = float(value)
    out_data: list[float] = []
    out_indices: list[int] = []
    out_indptr = [0]
    for column in spot_columns:
        for gene_idx in sorted(column):
            out_indices.append(gene_idx)
            out_data.append(column[gene_idx])
        out_indptr.append(len(out_data))
    return (
        np.asarray(out_data, dtype=float),
        np.asarray(out_indices, dtype=np.int64),
        np.asarray(out_indptr, dtype=np.int64),
    )


def _infer_matrix_scale(data: np.ndarray) -> str:
    if len(data) == 0:
        return "empty"
    sample = data[: min(len(data), 100000)].astype(float, copy=False)
    integer_fraction = float(np.mean(np.isclose(sample, np.round(sample))))
    return "counts" if integer_fraction > 0.995 and float(np.max(sample)) > 10 else "processed"


def load_h5ad_spatial(path: Path) -> H5adSpatialData:
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        matrix_container = (
            handle["raw"]
            if "raw" in handle and "X" in handle["raw"] and "var" in handle["raw"]
            else handle
        )
        var = matrix_container["var"]
        barcodes = _h5ad_index(obs)
        var_index = _h5ad_index(var)
        gene_names = var_index
        for candidate in ("gene_name", "gene_names", "feature_name", "features"):
            if candidate in var:
                values = _read_h5ad_strings(var[candidate])
                if len(values) == len(var_index):
                    gene_names = values
                    break

        x = matrix_container["X"]
        if isinstance(x, h5py.Dataset):
            dense = np.asarray(x[:], dtype=float)
            if dense.shape != (len(barcodes), len(gene_names)):
                raise ValueError("dense h5ad X dimensions do not match obs/var")
            matrix = SparseMatrixData.from_dense(
                dense.T,
                gene_names=gene_names,
                gene_ids=var_index,
                barcodes=barcodes,
            )
        else:
            encoding = _decode(x.attrs.get("encoding-type", "csr_matrix"))
            shape_attr = x.attrs.get("shape")
            if shape_attr is None and "shape" in x:
                shape_attr = x["shape"][:]
            shape = tuple(int(value) for value in shape_attr)
            if shape != (len(barcodes), len(gene_names)):
                raise ValueError(f"h5ad X shape {shape} does not match obs/var")
            data = x["data"][:]
            indices = x["indices"][:]
            indptr = x["indptr"][:]
            out_data, out_indices, out_indptr = _transpose_h5ad_sparse(data, indices, indptr, shape, encoding)
            matrix = SparseMatrixData(
                data=out_data,
                indices=out_indices,
                indptr=out_indptr,
                shape=(shape[1], shape[0]),
                gene_names=gene_names,
                barcodes=barcodes,
                gene_ids=var_index,
            )

        obs_columns: dict[str, list[str]] = {}
        for key, node in obs.items():
            if key == "_index":
                continue
            try:
                values = _read_h5ad_strings(node)
            except (ValueError, TypeError, OSError):
                continue
            if len(values) == len(barcodes):
                obs_columns[key] = values

        coordinates: np.ndarray | None = None
        if "obsm" in handle and "spatial" in handle["obsm"]:
            coordinates = np.asarray(handle["obsm"]["spatial"][:], dtype=float)
        if coordinates is None:
            for left, right in (
                ("pxl_col_in_fullres", "pxl_row_in_fullres"),
                ("array_col", "array_row"),
            ):
                if left in obs and right in obs:
                    coordinates = np.column_stack((obs[left][:], obs[right][:])).astype(float)
                    break
        if coordinates is None or coordinates.shape != (len(barcodes), 2):
            raise ValueError("h5ad has no usable two-dimensional spatial coordinates")

    return H5adSpatialData(
        matrix=matrix,
        coordinates=coordinates,
        obs_columns=obs_columns,
        matrix_scale=_infer_matrix_scale(matrix.data),
    )


def subset_h5ad_spatial(data: H5adSpatialData, spot_indices: Iterable[int]) -> H5adSpatialData:
    indices = np.asarray(list(spot_indices), dtype=int)
    return H5adSpatialData(
        matrix=subset_spots(data.matrix, indices.tolist()),
        coordinates=data.coordinates[indices],
        obs_columns={name: [values[idx] for idx in indices] for name, values in data.obs_columns.items()},
        matrix_scale=data.matrix_scale,
    )


def score_fixed_signatures(
    matrix: SparseMatrixData,
    signatures: Mapping[str, list[str]],
    normalize: bool = True,
) -> tuple[dict[str, np.ndarray], dict[str, int]]:
    working = normalize_csc_log1p(matrix) if normalize else matrix
    symbol_to_indices: dict[str, list[int]] = {}
    for gene_idx, symbol in enumerate(working.gene_names):
        symbol_to_indices.setdefault(symbol.casefold(), []).append(gene_idx)

    module_gene_indices: dict[str, dict[int, str]] = {}
    detected: dict[str, int] = {}
    gene_to_modules: dict[int, list[str]] = {}
    for module, genes in signatures.items():
        unique_symbols = {gene.casefold(): gene for gene in genes}
        mapping: dict[int, str] = {}
        detected_symbols = 0
        for folded, original in unique_symbols.items():
            indices = symbol_to_indices.get(folded, [])
            if not indices:
                continue
            detected_symbols += 1
            for gene_idx in indices:
                mapping[gene_idx] = original
                gene_to_modules.setdefault(gene_idx, []).append(module)
        module_gene_indices[module] = mapping
        detected[module] = detected_symbols

    scores = {module: np.zeros(working.shape[1], dtype=float) for module in signatures}
    for spot_idx in range(working.shape[1]):
        start = int(working.indptr[spot_idx])
        end = int(working.indptr[spot_idx + 1])
        for gene_idx, value in zip(working.indices[start:end], working.data[start:end], strict=True):
            for module in gene_to_modules.get(int(gene_idx), []):
                scores[module][spot_idx] += float(value)
    for module in scores:
        if detected[module]:
            scores[module] /= detected[module]
        else:
            scores[module][:] = np.nan
    return scores, detected


def zscore_vector(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    valid = np.isfinite(arr)
    result = np.full_like(arr, np.nan, dtype=float)
    if not np.any(valid):
        return result
    mean = float(np.mean(arr[valid]))
    sd = float(np.std(arr[valid]))
    result[valid] = 0.0 if sd == 0 else (arr[valid] - mean) / sd
    return result


def _sum_zscores(modules: Mapping[str, np.ndarray], names: tuple[str, ...]) -> np.ndarray:
    available = [zscore_vector(modules[name]) for name in names if name in modules]
    if not available:
        raise ValueError(f"none of the required modules are available: {names}")
    stacked = np.vstack(available)
    with np.errstate(invalid="ignore"):
        return np.nansum(stacked, axis=0)


def composite_framework_scores(modules: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    return {
        "mechanical_border_score": _sum_zscores(
            modules,
            ("CM_BZ1_TRANSITION", "CM_BZ2_MECHANICAL_EDGE"),
        ),
        "immune_fibrotic_activation_score": _sum_zscores(
            modules,
            ("CCR2_IL1B_MYLOID", "TGFB_SIGNALING", "FAP_POSTN_PATHO_FIBROBLAST"),
        ),
        "fibroblast_scar_repair_score": _sum_zscores(
            modules,
            (
                "FAP_POSTN_PATHO_FIBROBLAST",
                "ECM_REMODELING",
                "CTHRC1_REPARATIVE_CF",
                "MYOFIBROBLAST_CONTRACTILE",
            ),
        ),
    }


def raw_framework_scores(modules: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Return unstandardized fixed-signature sums for between-section summaries."""

    return {
        "raw_mechanical_border_score": modules["CM_BZ1_TRANSITION"] + modules["CM_BZ2_MECHANICAL_EDGE"],
        "raw_immune_fibrotic_activation_score": (
            modules["CCR2_IL1B_MYLOID"]
            + modules["TGFB_SIGNALING"]
            + modules["FAP_POSTN_PATHO_FIBROBLAST"]
        ),
        "raw_fibroblast_scar_repair_score": (
            modules["FAP_POSTN_PATHO_FIBROBLAST"]
            + modules["ECM_REMODELING"]
            + modules["CTHRC1_REPARATIVE_CF"]
            + modules["MYOFIBROBLAST_CONTRACTILE"]
        ),
    }


def boundary_transition_index(
    mechanical_border: np.ndarray,
    immune_fibrotic: np.ndarray,
    fibroblast_scar: np.ndarray,
) -> np.ndarray:
    mechanical = zscore_vector(mechanical_border)
    immune = zscore_vector(immune_fibrotic)
    scar = zscore_vector(fibroblast_scar)
    return (immune + scar - mechanical) / 3.0


def rank_biserial_effect(positive: Iterable[float], negative: Iterable[float]) -> float:
    left = np.asarray(list(positive), dtype=float)
    right = np.asarray(list(negative), dtype=float)
    left = left[np.isfinite(left)]
    right = right[np.isfinite(right)]
    if len(left) == 0 or len(right) == 0:
        return float("nan")
    wins = 0
    losses = 0
    for value in left:
        wins += int(np.sum(value > right))
        losses += int(np.sum(value < right))
    return float((wins - losses) / (len(left) * len(right)))


def _linear_r_squared(outcome: np.ndarray, predictors: np.ndarray) -> float:
    y = np.asarray(outcome, dtype=float)
    x = np.asarray(predictors, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    valid = np.isfinite(y) & np.all(np.isfinite(x), axis=1)
    y = y[valid]
    x = x[valid]
    if len(y) < x.shape[1] + 2:
        return float("nan")
    design = np.column_stack((np.ones(len(y)), x))
    coefficients, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    fitted = design @ coefficients
    total = float(np.sum((y - np.mean(y)) ** 2))
    if total == 0:
        return float("nan")
    return float(1.0 - np.sum((y - fitted) ** 2) / total)


def incremental_r_squared(
    outcome: np.ndarray,
    composition_predictors: np.ndarray,
    framework_score: np.ndarray,
) -> dict[str, float]:
    base = np.asarray(composition_predictors, dtype=float)
    if base.ndim == 1:
        base = base[:, None]
    score = np.asarray(framework_score, dtype=float)
    full = np.column_stack((base, score))
    base_r2 = _linear_r_squared(outcome, base)
    full_r2 = _linear_r_squared(outcome, full)
    return {
        "base_r_squared": base_r2,
        "full_r_squared": full_r2,
        "delta_r_squared": full_r2 - base_r2 if np.isfinite(base_r2) and np.isfinite(full_r2) else float("nan"),
    }
