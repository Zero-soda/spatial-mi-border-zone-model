#!/usr/bin/env python3
"""Quantify GSE214611 domain 3/4 spatial contact, distance, and score gradients."""

from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

from project_paths import project_root


ROOT = project_root(__file__)
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))

import numpy as np  # noqa: E402

from batch_map_gse214611_mi_spatial_risk import SAMPLES, sample_stage  # noqa: E402
from map_gse214611_d7_1_spatial_risk import safe_float  # noqa: E402


SCORES_BY_SPOT = ROOT / "results" / "tables" / "gse214611_d3_d7_signature_scores_by_spot.tsv"
OUT_TABLE_DIR = ROOT / "results" / "tables"

SCORE_COLUMNS = [
    "signature_mechanical_border_score",
    "signature_fibroblast_scar_score",
    "signature_fibrotic_risk",
    "prototype_fibrotic_risk",
]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def rows_by_sample(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["sample"]].append(row)
    return grouped


def pairwise_distances(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    delta = left[:, None, :] - right[None, :, :]
    return np.sqrt(np.sum(delta * delta, axis=2))


def median_nearest_neighbor_distance(coords: np.ndarray) -> float:
    distances = pairwise_distances(coords, coords)
    np.fill_diagonal(distances, np.inf)
    return float(np.median(np.min(distances, axis=1)))


def min_distance_to_group(coords: np.ndarray, target_coords: np.ndarray) -> np.ndarray:
    if len(target_coords) == 0:
        return np.full(len(coords), np.nan)
    return np.min(pairwise_distances(coords, target_coords), axis=1)


def nearest_opposite_distances(
    coords: np.ndarray,
    domains: np.ndarray,
    domain_value: int,
    opposite_value: int,
    opposite_coords: np.ndarray,
) -> np.ndarray:
    output = np.full(len(coords), np.nan)
    source_indices = np.where(domains == domain_value)[0]
    if len(source_indices) == 0 or len(opposite_coords) == 0:
        return output
    distances = pairwise_distances(coords[source_indices], opposite_coords)
    output[source_indices] = np.min(distances, axis=1)
    return output


def contact_pair_deltas(
    rows: list[dict[str, str]],
    domain3_indices: np.ndarray,
    domain4_indices: np.ndarray,
    distances_34: np.ndarray,
    threshold: float,
) -> dict[str, float]:
    contact_pairs = np.argwhere(distances_34 <= threshold)
    output = {"contact_pair_count": float(len(contact_pairs))}
    if len(contact_pairs) == 0:
        for column in SCORE_COLUMNS:
            output[f"contact_delta_{column}_domain4_minus_domain3"] = math.nan
        return output

    for column in SCORE_COLUMNS:
        deltas = []
        for domain3_pos, domain4_pos in contact_pairs:
            row3 = rows[int(domain3_indices[int(domain3_pos)])]
            row4 = rows[int(domain4_indices[int(domain4_pos)])]
            deltas.append(safe_float(row4[column]) - safe_float(row3[column]))
        output[f"contact_delta_{column}_domain4_minus_domain3"] = sum(deltas) / len(deltas)
    return output


def mean_for_mask(rows: list[dict[str, str]], mask: np.ndarray, column: str) -> float:
    values = [safe_float(row[column]) for row, keep in zip(rows, mask, strict=True) if keep]
    return sum(values) / len(values) if values else math.nan


def analyze_sample(sample: str, rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], dict[str, str]]:
    coords = np.array(
        [
            [safe_float(row["pxl_col_in_fullres"]), safe_float(row["pxl_row_in_fullres"])]
            for row in rows
        ],
        dtype=float,
    )
    domains = np.array([int(row["annotated"]) for row in rows], dtype=int)
    domain3_mask = domains == 3
    domain4_mask = domains == 4
    domain3_indices = np.where(domain3_mask)[0]
    domain4_indices = np.where(domain4_mask)[0]
    domain3_coords = coords[domain3_mask]
    domain4_coords = coords[domain4_mask]

    median_nn = median_nearest_neighbor_distance(coords)
    contact_threshold = median_nn * 1.35
    neighborhood_radius = median_nn * 2.25

    dist_to_domain3 = min_distance_to_group(coords, domain3_coords)
    dist_to_domain4 = min_distance_to_group(coords, domain4_coords)
    nearest_4_from_3 = nearest_opposite_distances(coords, domains, 3, 4, domain4_coords)
    nearest_3_from_4 = nearest_opposite_distances(coords, domains, 4, 3, domain3_coords)

    if len(domain3_coords) and len(domain4_coords):
        distances_34 = pairwise_distances(domain3_coords, domain4_coords)
        domain3_boundary_local = np.any(distances_34 <= contact_threshold, axis=1)
        domain4_boundary_local = np.any(distances_34 <= contact_threshold, axis=0)
    else:
        distances_34 = np.empty((0, 0))
        domain3_boundary_local = np.array([], dtype=bool)
        domain4_boundary_local = np.array([], dtype=bool)

    is_boundary = np.zeros(len(rows), dtype=bool)
    is_boundary[domain3_indices] = domain3_boundary_local
    is_boundary[domain4_indices] = domain4_boundary_local

    distance_matrix = pairwise_distances(coords, coords)
    neighbor_mask = distance_matrix <= neighborhood_radius

    local_mechanical = np.array([safe_float(row["signature_mechanical_border_score"]) for row in rows])
    local_scar = np.array([safe_float(row["signature_fibroblast_scar_score"]) for row in rows])
    local_risk = np.array([safe_float(row["signature_fibrotic_risk"]) for row in rows])
    local_mechanical_mean = neighbor_mask @ local_mechanical / np.sum(neighbor_mask, axis=1)
    local_scar_mean = neighbor_mask @ local_scar / np.sum(neighbor_mask, axis=1)
    local_risk_mean = neighbor_mask @ local_risk / np.sum(neighbor_mask, axis=1)

    spot_rows = []
    for idx, row in enumerate(rows):
        opposite_distance = ""
        if domains[idx] == 3:
            opposite_distance = f"{nearest_4_from_3[idx]:.6g}"
        elif domains[idx] == 4:
            opposite_distance = f"{nearest_3_from_4[idx]:.6g}"

        spot_rows.append(
            {
                "sample": sample,
                "stage": sample_stage(sample),
                "barcode": row["barcode"],
                "annotated": row["annotated"],
                "array_row": row["array_row"],
                "array_col": row["array_col"],
                "pxl_row_in_fullres": row["pxl_row_in_fullres"],
                "pxl_col_in_fullres": row["pxl_col_in_fullres"],
                "distance_to_domain3_fullres_px": f"{dist_to_domain3[idx]:.6g}",
                "distance_to_domain4_fullres_px": f"{dist_to_domain4[idx]:.6g}",
                "domain34_signed_distance_fullres_px": f"{dist_to_domain4[idx] - dist_to_domain3[idx]:.6g}",
                "nearest_opposite_domain34_distance_fullres_px": opposite_distance,
                "is_domain34_contact_boundary": str(bool(is_boundary[idx])),
                "local_mean_signature_mechanical_border_score": f"{local_mechanical_mean[idx]:.8g}",
                "local_mean_signature_fibroblast_scar_score": f"{local_scar_mean[idx]:.8g}",
                "local_mean_signature_fibrotic_risk": f"{local_risk_mean[idx]:.8g}",
                "local_mechanical_minus_scar_score": f"{local_mechanical_mean[idx] - local_scar_mean[idx]:.8g}",
                "signature_mechanical_border_score": row["signature_mechanical_border_score"],
                "signature_fibroblast_scar_score": row["signature_fibroblast_scar_score"],
                "signature_fibrotic_risk": row["signature_fibrotic_risk"],
                "prototype_fibrotic_risk": row["prototype_fibrotic_risk"],
            }
        )

    contact_deltas = contact_pair_deltas(rows, domain3_indices, domain4_indices, distances_34, contact_threshold)
    summary = {
        "sample": sample,
        "stage": sample_stage(sample),
        "n_spots": str(len(rows)),
        "n_domain3_spots": str(int(np.sum(domain3_mask))),
        "n_domain4_spots": str(int(np.sum(domain4_mask))),
        "median_nearest_neighbor_fullres_px": f"{median_nn:.6g}",
        "contact_threshold_fullres_px": f"{contact_threshold:.6g}",
        "n_domain3_boundary_spots": str(int(np.sum(is_boundary & domain3_mask))),
        "n_domain4_boundary_spots": str(int(np.sum(is_boundary & domain4_mask))),
        "fraction_domain3_boundary": f"{np.mean(is_boundary[domain3_mask]):.6g}",
        "fraction_domain4_boundary": f"{np.mean(is_boundary[domain4_mask]):.6g}",
        "median_domain3_to_domain4_distance_fullres_px": f"{np.nanmedian(nearest_4_from_3):.6g}",
        "median_domain4_to_domain3_distance_fullres_px": f"{np.nanmedian(nearest_3_from_4):.6g}",
    }
    for column in SCORE_COLUMNS:
        mean3 = mean_for_mask(rows, domain3_mask, column)
        mean4 = mean_for_mask(rows, domain4_mask, column)
        mean3_boundary = mean_for_mask(rows, is_boundary & domain3_mask, column)
        mean4_boundary = mean_for_mask(rows, is_boundary & domain4_mask, column)
        summary[f"mean_domain3_{column}"] = f"{mean3:.6g}"
        summary[f"mean_domain4_{column}"] = f"{mean4:.6g}"
        summary[f"mean_domain4_minus_domain3_{column}"] = f"{mean4 - mean3:.6g}"
        summary[f"boundary_mean_domain3_{column}"] = f"{mean3_boundary:.6g}"
        summary[f"boundary_mean_domain4_{column}"] = f"{mean4_boundary:.6g}"
        summary[f"boundary_mean_domain4_minus_domain3_{column}"] = f"{mean4_boundary - mean3_boundary:.6g}"
        summary[f"contact_delta_{column}_domain4_minus_domain3"] = (
            f"{contact_deltas[f'contact_delta_{column}_domain4_minus_domain3']:.6g}"
        )
    summary["contact_pair_count"] = str(int(contact_deltas["contact_pair_count"]))
    return spot_rows, summary


def write_rows(rows: list[dict[str, str]], path: Path) -> None:
    if not rows:
        raise ValueError("No rows to write")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_stage_summary(sample_summaries: list[dict[str, str]], path: Path) -> None:
    numeric_columns = [col for col in sample_summaries[0] if col not in {"sample", "stage"}]
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in sample_summaries:
        grouped[row["stage"]].append(row)

    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["stage", "n_samples", *[f"mean_{col}" for col in numeric_columns]])
        for stage in sorted(grouped):
            stage_rows = grouped[stage]
            values = []
            for col in numeric_columns:
                numeric = [safe_float(row[col]) for row in stage_rows]
                values.append(f"{sum(numeric) / len(numeric):.6g}")
            writer.writerow([stage, len(stage_rows), *values])


def main() -> None:
    rows = read_rows(SCORES_BY_SPOT)
    grouped = rows_by_sample(rows)
    all_spot_rows = []
    sample_summaries = []

    for sample in SAMPLES:
        spot_rows, summary = analyze_sample(sample, grouped[sample])
        all_spot_rows.extend(spot_rows)
        sample_summaries.append(summary)
        print(
            f"[domain34] {sample}: "
            f"D3 boundary {summary['fraction_domain3_boundary']}, "
            f"D4 boundary {summary['fraction_domain4_boundary']}, "
            f"contact pairs {summary['contact_pair_count']}"
        )

    write_rows(
        all_spot_rows,
        OUT_TABLE_DIR / "gse214611_d3_d7_domain34_spatial_features_by_spot.tsv",
    )
    write_rows(
        sample_summaries,
        OUT_TABLE_DIR / "gse214611_d3_d7_domain34_spatial_features_by_sample.tsv",
    )
    write_stage_summary(
        sample_summaries,
        OUT_TABLE_DIR / "gse214611_d3_d7_domain34_spatial_features_by_stage.tsv",
    )
    print(f"[summary] spots: {len(all_spot_rows)}, samples: {len(sample_summaries)}")


if __name__ == "__main__":
    main()
