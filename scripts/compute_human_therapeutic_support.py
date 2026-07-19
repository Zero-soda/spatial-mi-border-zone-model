#!/usr/bin/env python3
"""Compute single-section human spatial support for mouse-derived targets."""

from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

from project_paths import project_root
from typing import Mapping


ROOT = project_root(__file__)
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import numpy as np  # noqa: E402

from therapeutic_prioritization_utils import (  # noqa: E402
    build_knn_edges,
    extract_gene_rows,
    gene_detection_fraction,
    gene_pearson,
    load_10x_h5,
    normalize_csc_log1p,
    subset_spots,
    weighted_gene_mean,
    write_tsv,
)


TABLE_DIR = ROOT / "results" / "tables"
MOUSE_SUMMARY = TABLE_DIR / "gse214611_therapeutic_mouse_gene_evidence_summary.tsv"
ORTHOLOGUES = TABLE_DIR / "gse214611_mouse_human_one_to_one_orthologues.tsv"
HUMAN_SPOTS = TABLE_DIR / "gse214611_human_stemi_signature_scores_by_spot.tsv"
HUMAN_H5 = ROOT / "data/raw/gse214611/visium/GSM6613090_V_Human_STEMI/filtered_feature_bc_matrix.h5"
OUT_SUPPORT = TABLE_DIR / "gse214611_therapeutic_human_support.tsv"
OUT_NULL = TABLE_DIR / "gse214611_therapeutic_human_spatial_null_distribution.tsv"

RANDOM_SEED = 20260711
N_PER_STATE = 150
N_MATCHED_NULL = 100
CONTROL_GENES = {"Il1b", "Ccr2", "Ccl2", "Postn", "Cthrc1", "Nppa", "Nppb", "Xirp2", "Flnc", "Tgfb1", "Tgfbr1", "Lox"}

STATE_SCORE_COLUMN = {
    "mechanical_adverse": "mechanical_border_score",
    "immune_fibrotic_adverse": "immune_fibrotic_activation_score",
    "scar_repair_associated": "fibroblast_scar_repair_score",
    "unresolved": "mechanical_border_score",
}


def safe_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def score_human_support(
    detection: float,
    intended_r: float,
    moran_empirical_p: float,
    hotspot_overlap_fraction: float,
) -> float:
    detection_score = 4.0 * max(0.0, min(1.0, detection / 0.10))
    correlation_score = 6.0 * max(0.0, min(1.0, intended_r if math.isfinite(intended_r) else 0.0))
    if moran_empirical_p <= 0.01:
        moran_score = 6.0
    elif moran_empirical_p <= 0.05:
        moran_score = 5.0
    elif moran_empirical_p <= 0.10:
        moran_score = 3.0
    else:
        moran_score = 0.0
    hotspot_score = 4.0 * max(0.0, min(1.0, hotspot_overlap_fraction / 0.50))
    return min(20.0, detection_score + correlation_score + moran_score + hotspot_score)


def select_mouse_candidates(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    selected: dict[str, dict[str, str]] = {}
    by_state: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["replicate_rule_pass"] != "True" or safe_float(row["mouse_spatial_score"]) < 25:
            continue
        by_state[row["primary_state"]].append(row)
    for state in ["mechanical_adverse", "immune_fibrotic_adverse", "scar_repair_associated"]:
        ranked = sorted(by_state[state], key=lambda row: -safe_float(row["mouse_spatial_score"]))
        for row in ranked[:N_PER_STATE]:
            selected[row["gene_id"]] = row
    for row in rows:
        if row["gene_name"] in CONTROL_GENES:
            selected[row["gene_id"]] = row
    return list(selected.values())


def quantile_bins(values: np.ndarray, n_bins: int = 10) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    ranks[order] = np.arange(len(values), dtype=float)
    return np.minimum((ranks / max(len(values), 1) * n_bins).astype(int), n_bins - 1)


def candidate_null_indices(
    candidate_idx: int,
    mean_bins: np.ndarray,
    detection_bins: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    indices = np.arange(len(mean_bins))
    mask = (mean_bins == mean_bins[candidate_idx]) & (detection_bins == detection_bins[candidate_idx])
    pool = indices[mask & (indices != candidate_idx)]
    if len(pool) < N_MATCHED_NULL:
        pool = indices[(mean_bins == mean_bins[candidate_idx]) & (indices != candidate_idx)]
    if len(pool) < N_MATCHED_NULL:
        pool = indices[indices != candidate_idx]
    return rng.choice(pool, size=N_MATCHED_NULL, replace=len(pool) < N_MATCHED_NULL)


def graph_moran_rows(values: np.ndarray, edges: list[tuple[int, int]]) -> np.ndarray:
    n_genes, n_spots = values.shape
    if not edges:
        return np.full(n_genes, np.nan)
    centered = values - np.mean(values, axis=1, keepdims=True)
    denom = np.sum(centered * centered, axis=1)
    left = np.asarray([edge[0] for edge in edges], dtype=int)
    right = np.asarray([edge[1] for edge in edges], dtype=int)
    numerator = np.sum(centered[:, left] * centered[:, right], axis=1)
    result = np.full(n_genes, np.nan)
    np.divide((n_spots / len(edges)) * numerator, denom, out=result, where=denom > 0)
    return result


def hotspot_fraction(gene_values: np.ndarray, score_values: np.ndarray) -> float:
    positive = gene_values > 0
    if np.sum(positive) < 5:
        return 0.0
    gene_threshold = float(np.percentile(gene_values[positive], 90))
    score_threshold = float(np.percentile(score_values, 90))
    gene_hot = positive & (gene_values >= gene_threshold)
    score_hot = score_values >= score_threshold
    return float(np.sum(gene_hot & score_hot) / max(np.sum(gene_hot), 1))


def main() -> None:
    mouse_rows = read_tsv(MOUSE_SUMMARY)
    candidates = select_mouse_candidates(mouse_rows)
    candidate_by_mouse_id = {row["gene_id"]: row for row in candidates}
    orthologues = {
        row["mouse_gene_id"]: row
        for row in read_tsv(ORTHOLOGUES)
        if row["mouse_gene_id"] in candidate_by_mouse_id
    }
    human_rows = read_tsv(HUMAN_SPOTS)
    raw_matrix = load_10x_h5(HUMAN_H5)
    barcode_index = {barcode: idx for idx, barcode in enumerate(raw_matrix.barcodes)}
    aligned_rows = [row for row in human_rows if row["barcode"] in barcode_index]
    matrix = normalize_csc_log1p(subset_spots(raw_matrix, [barcode_index[row["barcode"]] for row in aligned_rows]))
    gene_id_index = {gene_id.split(".")[0]: idx for idx, gene_id in enumerate(matrix.gene_ids)}
    gene_name_index = {name.casefold(): idx for idx, name in enumerate(matrix.gene_names)}
    mapped = []
    for mouse_gene_id, orthologue in orthologues.items():
        human_gene_id = orthologue["human_gene_id"].split(".")[0]
        human_gene_name = orthologue["human_gene_name"]
        human_idx = gene_id_index.get(human_gene_id, gene_name_index.get(human_gene_name.casefold()))
        if human_idx is None:
            continue
        mapped.append((candidate_by_mouse_id[mouse_gene_id], orthologue, human_idx))
    if not mapped:
        raise RuntimeError("No mouse candidates mapped to the human matrix")

    all_mean = weighted_gene_mean(matrix, np.ones(matrix.shape[1], dtype=float))
    all_detection = gene_detection_fraction(matrix)
    correlations = {
        column: gene_pearson(matrix, np.asarray([safe_float(row[column]) for row in aligned_rows], dtype=float))
        for column in set(STATE_SCORE_COLUMN.values())
    }
    coords = np.asarray([[safe_float(row["array_row"]), safe_float(row["array_col"])] for row in aligned_rows])
    edges = build_knn_edges(coords, k=6)
    mean_bins = quantile_bins(all_mean)
    detection_bins = quantile_bins(all_detection)
    rng = np.random.default_rng(RANDOM_SEED)
    null_by_candidate: dict[int, np.ndarray] = {}
    required_indices = {human_idx for _, _, human_idx in mapped}
    for _, _, human_idx in mapped:
        null_indices = candidate_null_indices(human_idx, mean_bins, detection_bins, rng)
        null_by_candidate[human_idx] = null_indices
        required_indices.update(int(idx) for idx in null_indices)
    required = sorted(required_indices)
    dense = extract_gene_rows(matrix, required)
    moran = graph_moran_rows(dense, edges)
    moran_by_index = {gene_idx: moran[row_idx] for row_idx, gene_idx in enumerate(required)}
    dense_by_index = {gene_idx: dense[row_idx] for row_idx, gene_idx in enumerate(required)}

    support_rows = []
    null_rows = []
    for mouse_row, orthologue, human_idx in mapped:
        state = mouse_row["primary_state"]
        score_column = STATE_SCORE_COLUMN[state]
        intended_r = float(correlations[score_column][human_idx])
        observed_moran = float(moran_by_index[human_idx])
        null_indices = null_by_candidate[human_idx]
        null_values = np.asarray([moran_by_index[int(idx)] for idx in null_indices], dtype=float)
        null_values = null_values[np.isfinite(null_values)]
        empirical_p = float((np.sum(null_values >= observed_moran) + 1) / (len(null_values) + 1))
        score_values = np.asarray([safe_float(row[score_column]) for row in aligned_rows], dtype=float)
        hotspot = hotspot_fraction(dense_by_index[human_idx], score_values)
        support_score = score_human_support(
            detection=float(all_detection[human_idx]),
            intended_r=intended_r,
            moran_empirical_p=empirical_p,
            hotspot_overlap_fraction=hotspot,
        )
        support_rows.append(
            {
                "mouse_gene_id": mouse_row["gene_id"],
                "mouse_gene_name": mouse_row["gene_name"],
                "human_gene_id": orthologue["human_gene_id"],
                "human_gene_name": orthologue["human_gene_name"],
                "primary_state": state,
                "mouse_spatial_score": mouse_row["mouse_spatial_score"],
                "human_detected": float(all_detection[human_idx]) > 0,
                "human_detection_fraction": all_detection[human_idx],
                "human_mean_log_expression": all_mean[human_idx],
                "intended_score_column": score_column,
                "intended_score_pearson_r": intended_r,
                "graph_morans_i": observed_moran,
                "matched_null_n": len(null_values),
                "matched_null_moran_p95": np.percentile(null_values, 95) if len(null_values) else np.nan,
                "matched_null_moran_p99": np.percentile(null_values, 99) if len(null_values) else np.nan,
                "moran_empirical_upper_p": empirical_p,
                "hotspot_overlap_fraction": hotspot,
                "human_support_score": support_score,
                "interpretation": "single-section within-human spatial support; not external validation",
            }
        )
        for null_rank, null_idx in enumerate(null_indices, start=1):
            null_rows.append(
                {
                    "human_gene_id": orthologue["human_gene_id"],
                    "human_gene_name": orthologue["human_gene_name"],
                    "null_rank": null_rank,
                    "null_gene_id": matrix.gene_ids[int(null_idx)],
                    "null_gene_name": matrix.gene_names[int(null_idx)],
                    "null_graph_morans_i": moran_by_index[int(null_idx)],
                }
            )
    support_rows.sort(key=lambda row: (-safe_float(row["human_support_score"]), str(row["human_gene_name"])))
    write_tsv(OUT_SUPPORT, support_rows, list(support_rows[0]))
    write_tsv(OUT_NULL, null_rows, list(null_rows[0]))
    print(f"Wrote {OUT_SUPPORT} ({len(support_rows)} mapped candidates)")
    print(f"Wrote {OUT_NULL} ({len(null_rows)} matched-null rows)")


if __name__ == "__main__":
    main()
