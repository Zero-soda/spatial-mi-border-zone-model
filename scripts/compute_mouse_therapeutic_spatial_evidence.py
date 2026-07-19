#!/usr/bin/env python3
"""Compute sample-level all-gene spatial evidence for therapeutic prioritization."""

from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

from project_paths import project_root
from typing import Iterable, Mapping


ROOT = project_root(__file__)
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import numpy as np  # noqa: E402

from batch_map_gse214611_mi_spatial_risk import SAMPLES, locate_visium_dir, sample_stage  # noqa: E402
from therapeutic_prioritization_utils import (  # noqa: E402
    build_knn_edges,
    gene_detection_fraction,
    gene_linear_slope,
    gene_pearson,
    load_10x_h5,
    normalize_csc_log1p,
    subset_spots,
    weighted_gene_mean,
    weighted_gene_sum,
    write_tsv,
)


TABLE_DIR = ROOT / "results" / "tables"
MASTER_TABLE = TABLE_DIR / "gse214611_master_spot_level_source_table.tsv"
OUT_BY_SAMPLE = TABLE_DIR / "gse214611_therapeutic_mouse_gene_evidence_by_sample.tsv"
OUT_SUMMARY = TABLE_DIR / "gse214611_therapeutic_mouse_gene_evidence_summary.tsv"

SCORE_COLUMNS = {
    "mechanical_r": "mechanical_border_score",
    "immune_r": "immune_fibrotic_activation_score",
    "scar_r": "fibroblast_scar_repair_score",
}

CAUTION_GENES = {
    "actc1",
    "atp2a2",
    "cacna1c",
    "cthrc1",
    "flnc",
    "gja1",
    "itgb1",
    "kcnj2",
    "lo x".replace(" ", ""),
    "nppa",
    "nppb",
    "pln",
    "postn",
    "ryr2",
    "scn5a",
    "vcl",
    "xirp1",
    "xirp2",
}


def safe_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def finite_median(values: Iterable[float]) -> float:
    arr = np.asarray(list(values), dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.median(arr)) if len(arr) else float("nan")


def direction_consistency(values: Iterable[float], expected_sign: float) -> float:
    arr = np.asarray(list(values), dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return 0.0
    return float(np.mean(arr * expected_sign > 0))


def leave_one_out_direction(values: Iterable[float], expected_sign: float) -> float:
    arr = np.asarray(list(values), dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) < 2:
        return 0.0
    checks = []
    for idx in range(len(arr)):
        checks.append(float(np.median(np.delete(arr, idx))) * expected_sign > 0)
    return float(np.mean(checks))


def summarize_gene(
    gene_name: str,
    gene_id: str,
    sample_rows: list[Mapping[str, object]],
) -> dict[str, object]:
    """Assign a primary state and summarize replicate-level spatial evidence."""

    day7 = [row for row in sample_rows if row["stage"] == "day7_mi"]
    mechanical_effect = finite_median(safe_float(row["domain3_minus_domain4"]) for row in sample_rows)
    scar_effect = -finite_median(safe_float(row["domain3_minus_domain4"]) for row in day7)
    immune_effect = scar_effect
    mechanical_r = finite_median(safe_float(row["mechanical_r"]) for row in sample_rows)
    scar_r = finite_median(safe_float(row["scar_r"]) for row in day7)
    immune_r = finite_median(safe_float(row["immune_r"]) for row in day7)
    strengths = {
        "mechanical_adverse": max(0.0, mechanical_effect) * max(0.0, mechanical_r),
        "scar_repair_associated": max(0.0, scar_effect) * max(0.0, scar_r),
        "immune_fibrotic_adverse": max(0.0, immune_effect) * max(0.0, immune_r),
    }
    primary_state = max(strengths, key=strengths.get)
    if strengths[primary_state] <= 0:
        primary_state = "unresolved"

    if primary_state == "mechanical_adverse":
        relevant = sample_rows
        expected = 1.0
        intended_r_column = "mechanical_r"
    elif primary_state == "scar_repair_associated":
        relevant = day7
        expected = -1.0
        intended_r_column = "scar_r"
    elif primary_state == "immune_fibrotic_adverse":
        relevant = day7
        expected = -1.0
        intended_r_column = "immune_r"
    else:
        relevant = sample_rows
        expected = 1.0
        intended_r_column = "mechanical_r"

    contrasts = [safe_float(row["domain3_minus_domain4"]) for row in relevant]
    slope_values = [safe_float(row["signed_distance_slope"]) for row in relevant]
    edge_values = [safe_float(row["graph_edge_domain4_minus_domain3"]) for row in relevant]
    intended_values = [safe_float(row[intended_r_column]) for row in relevant]
    detections = [safe_float(row["detection_all"]) for row in relevant]
    replicate_consistency = direction_consistency(contrasts, expected)
    slope_consistency = direction_consistency(slope_values, expected)
    edge_consistency = direction_consistency(edge_values, -expected)
    loso_consistency = leave_one_out_direction(contrasts, expected)
    intended_r = finite_median(intended_values)
    detection = finite_median(detections)
    relevant_effect = finite_median(value * expected for value in contrasts)
    replicate_rule_pass = (
        replicate_consistency >= (5 / 6)
        if primary_state == "mechanical_adverse" and len(relevant) >= 6
        else replicate_consistency >= 1.0
    )
    base_score = (
        12.0 * replicate_consistency
        + 4.0 * slope_consistency
        + 4.0 * edge_consistency
        + 6.0 * max(0.0, min(1.0, intended_r if math.isfinite(intended_r) else 0.0))
        + 3.0 * max(0.0, min(1.0, (detection if math.isfinite(detection) else 0.0) / 0.20))
        + 3.0 * loso_consistency
    )
    return {
        "gene_id": gene_id,
        "gene_name": gene_name,
        "primary_state": primary_state,
        "structural_or_repair_caution": gene_name.casefold() in CAUTION_GENES,
        "n_relevant_samples": len(relevant),
        "relevant_effect_median": relevant_effect,
        "relevant_direction_consistency": replicate_consistency,
        "boundary_slope_consistency": slope_consistency,
        "graph_edge_consistency": edge_consistency,
        "leave_one_out_direction_consistency": loso_consistency,
        "intended_score_correlation_median": intended_r,
        "detection_fraction_median": detection,
        "replicate_rule_pass": replicate_rule_pass,
        "mouse_score_before_effect_percentile": base_score,
        "mechanical_state_strength": strengths["mechanical_adverse"],
        "immune_state_strength": strengths["immune_fibrotic_adverse"],
        "scar_state_strength": strengths["scar_repair_associated"],
    }


def read_master_rows() -> dict[str, list[dict[str, str]]]:
    by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    with MASTER_TABLE.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["species"] == "Mus musculus":
                by_sample[row["sample"]].append(row)
    return by_sample


def find_matrix_h5(visium_dir: Path) -> Path:
    candidates = list(visium_dir.glob("**/filtered_feature_bc_matrix.h5"))
    if (visium_dir / "filtered_feature_bc_matrix.h5").exists():
        candidates.append(visium_dir / "filtered_feature_bc_matrix.h5")
    unique = sorted(set(candidates))
    if len(unique) != 1:
        raise RuntimeError(f"Expected one filtered H5 under {visium_dir}, found {unique}")
    return unique[0]


def graph_edge_weights(rows: list[dict[str, str]], edges: list[tuple[int, int]]) -> tuple[np.ndarray, int]:
    domains = np.asarray([row["author_domain"] for row in rows])
    weights = np.zeros(len(rows), dtype=float)
    cross_edges = []
    for left, right in edges:
        pair = {domains[left], domains[right]}
        if pair != {"3", "4"}:
            continue
        d4 = left if domains[left] == "4" else right
        d3 = right if d4 == left else left
        cross_edges.append((d3, d4))
    if not cross_edges:
        return weights, 0
    for d3, d4 in cross_edges:
        weights[d4] += 1.0
        weights[d3] -= 1.0
    weights /= len(cross_edges)
    return weights, len(cross_edges)


def compute_sample(sample: str, master_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    config = SAMPLES[sample]
    visium_dir = locate_visium_dir(config.unpack_dir)
    matrix = normalize_csc_log1p(load_10x_h5(find_matrix_h5(visium_dir)))
    barcode_index = {barcode: idx for idx, barcode in enumerate(matrix.barcodes)}
    aligned_rows = [row for row in master_rows if row["barcode"] in barcode_index]
    spot_indices = [barcode_index[row["barcode"]] for row in aligned_rows]
    matrix = subset_spots(matrix, spot_indices)
    domains = np.asarray([row["author_domain"] for row in aligned_rows])
    d3_mask = domains == "3"
    d4_mask = domains == "4"
    if np.sum(d3_mask) < 20 or np.sum(d4_mask) < 20:
        raise RuntimeError(f"Insufficient domain 3/4 spots for {sample}")
    d3_mean = weighted_gene_mean(matrix, d3_mask.astype(float))
    d4_mean = weighted_gene_mean(matrix, d4_mask.astype(float))
    detect_d3 = gene_detection_fraction(matrix, d3_mask)
    detect_d4 = gene_detection_fraction(matrix, d4_mask)
    detect_all = gene_detection_fraction(matrix)
    signed_distance = np.asarray(
        [safe_float(row["domain34_signed_distance_fullres_px"]) for row in aligned_rows], dtype=float
    )
    slopes = gene_linear_slope(matrix, signed_distance)
    correlations = {
        output: gene_pearson(matrix, np.asarray([safe_float(row[column]) for row in aligned_rows], dtype=float))
        for output, column in SCORE_COLUMNS.items()
    }
    coords = np.asarray(
        [[safe_float(row["array_row"]), safe_float(row["array_col"])] for row in aligned_rows], dtype=float
    )
    edges = build_knn_edges(coords, k=6)
    edge_weights, n_cross_edges = graph_edge_weights(aligned_rows, edges)
    edge_delta = weighted_gene_sum(matrix, edge_weights) if n_cross_edges else np.full(matrix.shape[0], np.nan)
    rows: list[dict[str, object]] = []
    for gene_idx, (gene_id, gene_name) in enumerate(zip(matrix.gene_ids, matrix.gene_names, strict=True)):
        eligible = max(detect_d3[gene_idx], detect_d4[gene_idx]) >= 0.05
        rows.append(
            {
                "sample": sample,
                "stage": sample_stage(sample),
                "gene_id": gene_id,
                "gene_name": gene_name,
                "n_spots": len(aligned_rows),
                "n_domain3": int(np.sum(d3_mask)),
                "n_domain4": int(np.sum(d4_mask)),
                "n_domain34_graph_edges": n_cross_edges,
                "detection_domain3": detect_d3[gene_idx],
                "detection_domain4": detect_d4[gene_idx],
                "detection_all": detect_all[gene_idx],
                "mean_domain3": d3_mean[gene_idx],
                "mean_domain4": d4_mean[gene_idx],
                "domain3_minus_domain4": d3_mean[gene_idx] - d4_mean[gene_idx],
                "signed_distance_slope": slopes[gene_idx],
                "graph_edge_domain4_minus_domain3": edge_delta[gene_idx],
                "mechanical_r": correlations["mechanical_r"][gene_idx],
                "immune_r": correlations["immune_r"][gene_idx],
                "scar_r": correlations["scar_r"][gene_idx],
                "eligible_detection": eligible,
            }
        )
    print(f"[{sample}] spots={len(aligned_rows)} genes={len(rows)} cross_edges={n_cross_edges}")
    return rows


def percentile_ranks(values: np.ndarray) -> np.ndarray:
    result = np.zeros(len(values), dtype=float)
    finite = np.isfinite(values)
    order = np.argsort(values[finite], kind="mergesort")
    if np.sum(finite) == 1:
        result[finite] = 1.0
        return result
    ranks = np.empty(len(order), dtype=float)
    ranks[order] = np.arange(len(order), dtype=float) / max(len(order) - 1, 1)
    result[finite] = ranks
    return result


def build_summaries(sample_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_gene: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in sample_rows:
        if str(row["eligible_detection"]) not in {"True", "1", "true"} and row["eligible_detection"] is not True:
            continue
        by_gene[(str(row["gene_id"]), str(row["gene_name"]))].append(row)
    summaries = [summarize_gene(gene_name, gene_id, rows) for (gene_id, gene_name), rows in by_gene.items()]
    for state in ["mechanical_adverse", "immune_fibrotic_adverse", "scar_repair_associated"]:
        indices = [idx for idx, row in enumerate(summaries) if row["primary_state"] == state]
        effects = np.asarray([safe_float(summaries[idx]["relevant_effect_median"]) for idx in indices], dtype=float)
        percentiles = percentile_ranks(effects)
        for local_idx, row_idx in enumerate(indices):
            summaries[row_idx]["effect_magnitude_percentile"] = percentiles[local_idx]
    for row in summaries:
        percentile = safe_float(row.get("effect_magnitude_percentile", 0.0))
        if not math.isfinite(percentile):
            percentile = 0.0
        row["mouse_spatial_score"] = min(
            40.0,
            safe_float(row["mouse_score_before_effect_percentile"]) + 8.0 * percentile,
        )
    summaries.sort(key=lambda row: (-safe_float(row["mouse_spatial_score"]), str(row["gene_name"])))
    return summaries


def main() -> None:
    master_by_sample = read_master_rows()
    all_rows: list[dict[str, object]] = []
    for sample in SAMPLES:
        all_rows.extend(compute_sample(sample, master_by_sample[sample]))
    sample_fields = list(all_rows[0])
    write_tsv(OUT_BY_SAMPLE, all_rows, sample_fields)
    summaries = build_summaries(all_rows)
    summary_fields = list(summaries[0])
    write_tsv(OUT_SUMMARY, summaries, summary_fields)
    print(f"Wrote {OUT_BY_SAMPLE}")
    print(f"Wrote {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
