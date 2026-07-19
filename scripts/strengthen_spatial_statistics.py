#!/usr/bin/env python3
"""Spatial-statistical strengthening analyses for the MI border-zone manuscript."""

from __future__ import annotations

import sys
import csv
from collections import defaultdict
from pathlib import Path

from project_paths import project_root


ROOT = project_root(__file__)
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import numpy as np  # noqa: E402

from spatial_model_utils import (  # noqa: E402
    bootstrap_ci,
    cohens_d,
    domain_label_permutation_margin,
    gearys_c,
    morans_i,
    pairwise_distances,
    read_tsv,
    safe_float,
    signed_distance_slope,
    spatial_block_permutation_margin,
)


TABLE_DIR = ROOT / "results" / "tables"
SPOT = TABLE_DIR / "gse214611_d3_d7_signature_scores_by_spot.tsv"
SAMPLE_DOMAIN = TABLE_DIR / "gse214611_d3_d7_signature_scores_by_sample_domain.tsv"
BOUNDARY_SPOT = TABLE_DIR / "gse214611_d3_d7_domain34_spatial_features_by_spot.tsv"

LABEL_PERMUTATIONS = 5000
BLOCK_PERMUTATIONS = 2000
SPATIAL_BLOCK_GRID_SIZE = 4
RANDOM_SEED = 614611

SCORES = [
    ("signature_mechanical_border_score", "mechanical_border"),
    ("signature_fibroblast_scar_score", "fibroblast_scar_repair"),
    ("signature_fibrotic_risk", "immune_fibrotic_activation_proxy"),
]


def sample_stage(sample: str) -> str:
    return "day3_mi" if sample.startswith("D3_") else "day7_mi"


def format_float(value: float, precision: int = 8) -> str:
    if not np.isfinite(value):
        return "nan"
    return f"{value:.{precision}g}"


def seed_for(prefix: int, sample: str, score_name: str) -> int:
    return prefix + sum(ord(char) for char in f"{sample}:{score_name}")


def rows_by_sample(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["sample"]].append(row)
    return grouped


def median_nearest_neighbor_distance(coords: np.ndarray) -> float:
    valid = np.all(np.isfinite(coords), axis=1)
    coords = coords[valid]
    if len(coords) < 2:
        return float("nan")
    distances = pairwise_distances(coords)
    np.fill_diagonal(distances, np.inf)
    nearest = np.min(distances, axis=1)
    nearest = nearest[np.isfinite(nearest)]
    if len(nearest) == 0:
        return float("nan")
    return float(np.median(nearest))


def percentile_text(values: np.ndarray, percentile: float) -> str:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return "NA"
    return format_float(float(np.percentile(values, percentile)))


def write_lf_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_effect_sizes(rows: list[dict[str, str]]) -> None:
    values_by_key: dict[tuple[str, str, str], dict[str, float]] = defaultdict(dict)
    for row in rows:
        domain = row["annotated"]
        if domain not in {"3", "4"}:
            continue
        stage = sample_stage(row["sample"])
        for score_col, _ in SCORES:
            values_by_key[(stage, row["sample"], domain)][score_col] = safe_float(row[f"mean_{score_col}"])

    out = []
    for stage in ["day3_mi", "day7_mi"]:
        stage_samples = sorted({sample for key_stage, sample, _ in values_by_key if key_stage == stage})
        for score_col, score_name in SCORES:
            paired_samples = [
                sample
                for sample in stage_samples
                if score_col in values_by_key[(stage, sample, "3")]
                and score_col in values_by_key[(stage, sample, "4")]
            ]
            d3 = [values_by_key[(stage, sample, "3")][score_col] for sample in paired_samples]
            d4 = [values_by_key[(stage, sample, "4")][score_col] for sample in paired_samples]
            margins = [left - right for left, right in zip(d3, d4, strict=True)]
            margin_mean, margin_low, margin_high = bootstrap_ci(margins, seed=RANDOM_SEED)
            out.append(
                {
                    "stage": stage,
                    "score": score_name,
                    "contrast": "domain3_minus_domain4",
                    "n_sample_pairs": len(margins),
                    "mean_margin": format_float(margin_mean),
                    "ci95_low": format_float(margin_low),
                    "ci95_high": format_float(margin_high),
                    "cohens_d_domain3_vs_domain4": format_float(cohens_d(d3, d4)),
                    "domain3_values": ";".join(format_float(value, precision=6) for value in d3),
                    "domain4_values": ";".join(format_float(value, precision=6) for value in d4),
                }
            )

    write_lf_tsv(
        TABLE_DIR / "gse214611_stage_domain_effect_sizes.tsv",
        out,
        [
            "stage",
            "score",
            "contrast",
            "n_sample_pairs",
            "mean_margin",
            "ci95_low",
            "ci95_high",
            "cohens_d_domain3_vs_domain4",
            "domain3_values",
            "domain4_values",
        ],
    )


def write_spatial_autocorrelation(spot_rows: list[dict[str, str]]) -> None:
    out = []
    for sample, rows in sorted(rows_by_sample(spot_rows).items()):
        coords = np.array(
            [[safe_float(row["pxl_col_in_fullres"]), safe_float(row["pxl_row_in_fullres"])] for row in rows],
            dtype=float,
        )
        median_nn = median_nearest_neighbor_distance(coords)
        radius = median_nn * 2.25
        for score_col, score_name in SCORES:
            values = np.array([safe_float(row[score_col]) for row in rows], dtype=float)
            out.append(
                {
                    "sample": sample,
                    "stage": sample_stage(sample),
                    "score": score_name,
                    "n_spots": len(rows),
                    "neighbor_radius_fullres_px": format_float(radius),
                    "morans_i": format_float(morans_i(coords, values, radius)),
                    "gearys_c": format_float(gearys_c(coords, values, radius)),
                }
            )

    write_lf_tsv(
        TABLE_DIR / "gse214611_spatial_autocorrelation.tsv",
        out,
        ["sample", "stage", "score", "n_spots", "neighbor_radius_fullres_px", "morans_i", "gearys_c"],
    )


def write_domain_permutation(spot_rows: list[dict[str, str]]) -> None:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in spot_rows:
        if row["annotated"] in {"3", "4"}:
            grouped[(row["sample"], sample_stage(row["sample"]))].append(row)

    out = []
    for (sample, stage), rows in sorted(grouped.items()):
        labels = np.array([row["annotated"] for row in rows])
        coords = np.array(
            [[safe_float(row["pxl_col_in_fullres"]), safe_float(row["pxl_row_in_fullres"])] for row in rows],
            dtype=float,
        )
        for score_col, score_name in SCORES:
            values = np.array([safe_float(row[score_col]) for row in rows], dtype=float)
            observed, p_upper, null = domain_label_permutation_margin(
                values,
                labels,
                positive_label="3",
                negative_label="4",
                n_perm=LABEL_PERMUTATIONS,
                seed=seed_for(RANDOM_SEED, sample, score_name),
            )
            _, block_p_upper, block_null, n_blocks = spatial_block_permutation_margin(
                values,
                labels,
                coords,
                positive_label="3",
                negative_label="4",
                n_perm=BLOCK_PERMUTATIONS,
                seed=seed_for(RANDOM_SEED + 200000, sample, score_name),
                grid_size=SPATIAL_BLOCK_GRID_SIZE,
            )
            out.append(
                {
                    "sample": sample,
                    "stage": stage,
                    "score": score_name,
                    "contrast": "domain3_minus_domain4",
                    "n_permutations": LABEL_PERMUTATIONS,
                    "n_block_permutations": BLOCK_PERMUTATIONS,
                    "n_spatial_blocks": n_blocks,
                    "observed_margin": format_float(observed),
                    "permutation_p_upper": format_float(p_upper),
                    "null_p95": percentile_text(null, 95),
                    "null_p99": percentile_text(null, 99),
                    "block_permutation_p_upper": format_float(block_p_upper),
                    "block_null_p95": percentile_text(block_null, 95),
                    "block_null_p99": percentile_text(block_null, 99),
                }
            )

    write_lf_tsv(
        TABLE_DIR / "gse214611_domain_label_permutation.tsv",
        out,
        [
            "sample",
            "stage",
            "score",
            "contrast",
            "n_permutations",
            "n_block_permutations",
            "n_spatial_blocks",
            "observed_margin",
            "permutation_p_upper",
            "null_p95",
            "null_p99",
            "block_permutation_p_upper",
            "block_null_p95",
            "block_null_p99",
        ],
    )


def signed_distance_definition(rows: list[dict[str, str]]) -> str:
    domain3 = [
        safe_float(row["domain34_signed_distance_fullres_px"])
        for row in rows
        if row["annotated"] == "3"
    ]
    domain4 = [
        safe_float(row["domain34_signed_distance_fullres_px"])
        for row in rows
        if row["annotated"] == "4"
    ]
    median3 = float(np.nanmedian(domain3)) if domain3 else float("nan")
    median4 = float(np.nanmedian(domain4)) if domain4 else float("nan")
    if median3 > 0 and median4 < 0:
        return "domain3_side_positive_domain4_side_negative"
    if median3 < 0 and median4 > 0:
        return "domain3_side_negative_domain4_side_positive"
    return "as_reported_in_domain34_signed_distance_fullres_px"


def write_boundary_distance_gradients(boundary_rows: list[dict[str, str]]) -> None:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in boundary_rows:
        if row["annotated"] in {"3", "4"}:
            grouped[(row["sample"], sample_stage(row["sample"]))].append(row)

    out = []
    for (sample, stage), rows in sorted(grouped.items()):
        distance = np.array([safe_float(row["domain34_signed_distance_fullres_px"]) for row in rows], dtype=float)
        distance_definition = signed_distance_definition(rows)
        for score_col, score_name in SCORES:
            score = np.array([safe_float(row[score_col]) for row in rows], dtype=float)
            slope, intercept = signed_distance_slope(distance, score)
            out.append(
                {
                    "sample": sample,
                    "stage": stage,
                    "score": score_name,
                    "n_domain34_spots": len(rows),
                    "signed_distance_definition": distance_definition,
                    "linear_slope_per_fullres_px": format_float(slope, precision=10),
                    "linear_intercept": format_float(intercept),
                    "slope_per_1000_fullres_px": format_float(slope * 1000),
                }
            )

    write_lf_tsv(
        TABLE_DIR / "gse214611_boundary_distance_gradients.tsv",
        out,
        [
            "sample",
            "stage",
            "score",
            "n_domain34_spots",
            "signed_distance_definition",
            "linear_slope_per_fullres_px",
            "linear_intercept",
            "slope_per_1000_fullres_px",
        ],
    )


def main() -> None:
    spot_rows = read_tsv(SPOT)
    sample_domain_rows = read_tsv(SAMPLE_DOMAIN)
    boundary_rows = read_tsv(BOUNDARY_SPOT)
    write_effect_sizes(sample_domain_rows)
    write_spatial_autocorrelation(spot_rows)
    write_domain_permutation(spot_rows)
    write_boundary_distance_gradients(boundary_rows)
    print("Wrote upgraded spatial-statistics tables")


if __name__ == "__main__":
    main()
