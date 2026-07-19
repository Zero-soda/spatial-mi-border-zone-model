#!/usr/bin/env python3
"""Benchmark the three-output framework against simpler external-study baselines."""

from __future__ import annotations

import csv
from pathlib import Path

from project_paths import project_root

import numpy as np

from external_spatial_validation_utils import (
    incremental_r_squared,
    rank_biserial_effect,
    zscore_vector,
)


ROOT = project_root(__file__)
TABLE_ROOT = ROOT / "results" / "tables" / "external_validation"


COMPARISONS = (
    {
        "comparison": "Yamada prespecified tissue-mean mechanical",
        "dataset": "GSE176092",
        "axis": "mechanical",
        "positive": {"day1"},
        "negative": {"day14"},
        "full": "mean_raw_mechanical_border_score",
        "baseline": "mean_source_bz_baseline",
    },
    {
        "comparison": "Yamada spatial mechanical-to-scar transition",
        "dataset": "GSE176092",
        "axis": "transition",
        "positive": {"day1"},
        "negative": {"day7"},
        "full": "framework_spatial_axis_dominance",
        "baseline": "baseline_spatial_axis_dominance",
    },
    {
        "comparison": "Yamada day7 scar",
        "dataset": "GSE176092",
        "axis": "scar",
        "positive": {"day7"},
        "negative": {"day1"},
        "full": "mean_raw_fibroblast_scar_repair_score",
        "baseline": "mean_generic_fibrosis_baseline",
    },
    {
        "comparison": "Hernandez scar maturation",
        "dataset": "GSE265828",
        "axis": "scar",
        "positive": {"day5"},
        "negative": {"control", "day3"},
        "full": "mean_raw_fibroblast_scar_repair_score",
        "baseline": "mean_generic_fibrosis_baseline",
    },
    {
        "comparison": "Kuppe border-ischaemic mechanical",
        "dataset": "KUPPE2022",
        "axis": "mechanical",
        "positive": {"border", "ischaemic"},
        "negative": {"control", "remote"},
        "full": "mean_raw_mechanical_border_score",
        "baseline": "mean_source_bz_baseline",
    },
    {
        "comparison": "Kuppe fibrotic scar",
        "dataset": "KUPPE2022",
        "axis": "scar",
        "positive": {"fibrotic"},
        "negative": {"control", "remote"},
        "full": "mean_raw_fibroblast_scar_repair_score",
        "baseline": "mean_generic_fibrosis_baseline",
    },
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(row: dict[str, str], column: str) -> float:
    if column == "framework_spatial_axis_dominance":
        return as_float(row, "graph_morans_i_mechanical_border_score") - as_float(
            row, "graph_morans_i_fibroblast_scar_repair_score"
        )
    if column == "baseline_spatial_axis_dominance":
        return as_float(row, "graph_morans_i_source_bz_baseline") - as_float(
            row, "graph_morans_i_generic_fibrosis_baseline"
        )
    try:
        return float(row[column])
    except (KeyError, TypeError, ValueError):
        return float("nan")


def group_label(row: dict[str, str]) -> str:
    return row.get("coarse_region") or row.get("stage_region", "")


def leave_one_out_direction(values: np.ndarray, labels: np.ndarray) -> tuple[int, int, float]:
    tested = 0
    positive = 0
    for omit in range(len(values)):
        keep = np.arange(len(values)) != omit
        left = values[keep & (labels == 1)]
        right = values[keep & (labels == 0)]
        if len(left) == 0 or len(right) == 0:
            continue
        tested += 1
        positive += int(rank_biserial_effect(left, right) > 0)
    return positive, tested, positive / tested if tested else float("nan")


def benchmark_model(
    comparison: dict[str, object],
    rows: list[dict[str, str]],
    model_name: str,
    score_column: str,
) -> dict[str, object]:
    selected = [
        row
        for row in rows
        if group_label(row) in comparison["positive"] or group_label(row) in comparison["negative"]
    ]
    labels = np.asarray([1 if group_label(row) in comparison["positive"] else 0 for row in selected], dtype=int)
    values = np.asarray([as_float(row, score_column) for row in selected], dtype=float)
    finite = np.isfinite(values)
    labels = labels[finite]
    values = values[finite]
    positive_values = values[labels == 1]
    negative_values = values[labels == 0]
    loo_positive, loo_tested, loo_fraction = leave_one_out_direction(values, labels)
    composition = np.column_stack(
        (
            [as_float(row, "mean_cardiomyocyte_composition_surrogate") for row in selected],
            [as_float(row, "mean_fibroblast_composition_surrogate") for row in selected],
        )
    )[finite]
    r2 = incremental_r_squared(labels.astype(float), composition, values)
    return {
        "comparison": comparison["comparison"],
        "dataset": comparison["dataset"],
        "axis": comparison["axis"],
        "model": model_name,
        "score_column": score_column,
        "positive_groups": ";".join(sorted(comparison["positive"])),
        "negative_groups": ";".join(sorted(comparison["negative"])),
        "n_positive_units": len(positive_values),
        "n_negative_units": len(negative_values),
        "positive_mean": float(np.mean(positive_values)) if len(positive_values) else np.nan,
        "negative_mean": float(np.mean(negative_values)) if len(negative_values) else np.nan,
        "mean_difference": (
            float(np.mean(positive_values) - np.mean(negative_values))
            if len(positive_values) and len(negative_values)
            else np.nan
        ),
        "rank_biserial_effect": rank_biserial_effect(positive_values, negative_values),
        "loo_positive": loo_positive,
        "loo_tested": loo_tested,
        "loo_direction_consistency": loo_fraction,
        **r2,
    }


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    mouse_rows = read_tsv(TABLE_ROOT / "external_mouse_section_summary.tsv")
    area_rows = read_tsv(TABLE_ROOT / "external_mouse_author_area_summary.tsv")
    human_rows = read_tsv(TABLE_ROOT / "kuppe_human_section_summary.tsv")
    all_rows = mouse_rows + human_rows
    benchmark_rows: list[dict[str, object]] = []
    for comparison in COMPARISONS:
        dataset_rows = [row for row in all_rows if row["dataset"] == comparison["dataset"]]
        benchmark_rows.append(benchmark_model(comparison, dataset_rows, "three_output_framework", str(comparison["full"])))
        benchmark_rows.append(benchmark_model(comparison, dataset_rows, "simple_published_or_generic_baseline", str(comparison["baseline"])))

    improvement_rows: list[dict[str, object]] = []
    for comparison in COMPARISONS:
        matched = [row for row in benchmark_rows if row["comparison"] == comparison["comparison"]]
        full = next(row for row in matched if row["model"] == "three_output_framework")
        baseline = next(row for row in matched if row["model"] != "three_output_framework")
        improvement_rows.append(
            {
                "comparison": comparison["comparison"],
                "dataset": comparison["dataset"],
                "axis": comparison["axis"],
                "framework_rank_biserial": full["rank_biserial_effect"],
                "baseline_rank_biserial": baseline["rank_biserial_effect"],
                "rank_biserial_gain": float(full["rank_biserial_effect"]) - float(baseline["rank_biserial_effect"]),
                "framework_loo_consistency": full["loo_direction_consistency"],
                "baseline_loo_consistency": baseline["loo_direction_consistency"],
                "framework_delta_r_squared": full["delta_r_squared"],
                "baseline_delta_r_squared": baseline["delta_r_squared"],
            }
        )

    transition_rows: list[dict[str, object]] = []
    for dataset in sorted({row["dataset"] for row in all_rows}):
        dataset_rows = [row for row in all_rows if row["dataset"] == dataset]
        raw_mechanical = np.asarray([as_float(row, "mean_raw_mechanical_border_score") for row in dataset_rows])
        raw_immune = np.asarray([as_float(row, "mean_raw_immune_fibrotic_activation_score") for row in dataset_rows])
        raw_scar = np.asarray([as_float(row, "mean_raw_fibroblast_scar_repair_score") for row in dataset_rows])
        section_transition = (zscore_vector(raw_immune) + zscore_vector(raw_scar) - zscore_vector(raw_mechanical)) / 3.0
        for row, value in zip(dataset_rows, section_transition, strict=True):
            transition_rows.append(
                {
                    "dataset": dataset,
                    "sample": row["sample"],
                    "patient": row.get("patient", ""),
                    "stage_region": row["stage_region"],
                    "coarse_region": group_label(row),
                    "section_boundary_transition_index": float(value),
                    "spot_median_boundary_transition_index": as_float(row, "median_boundary_transition_index"),
                    "spot_p90_boundary_transition_index": as_float(row, "p90_boundary_transition_index"),
                }
            )

    area_gradient_rows: list[dict[str, object]] = []
    for sample in sorted({row["sample"] for row in area_rows}):
        sample_rows = {row["author_area"]: row for row in area_rows if row["sample"] == sample}
        if not {"RZ", "BZ2", "IZ"}.issubset(sample_rows):
            continue
        rz = sample_rows["RZ"]
        bz2 = sample_rows["BZ2"]
        iz = sample_rows["IZ"]
        for axis, positive_area, positive_row, full_column, baseline_column in (
            (
                "mechanical",
                "BZ2",
                bz2,
                "mean_mechanical_border_score",
                "mean_source_bz_baseline",
            ),
            (
                "scar",
                "IZ",
                iz,
                "mean_fibroblast_scar_repair_score",
                "mean_generic_fibrosis_baseline",
            ),
            (
                "transition",
                "IZ",
                iz,
                "mean_boundary_transition_index",
                "mean_boundary_transition_index",
            ),
        ):
            area_gradient_rows.append(
                {
                    "dataset": positive_row["dataset"],
                    "sample": sample,
                    "stage_region": positive_row["stage_region"],
                    "axis": axis,
                    "positive_area": positive_area,
                    "negative_area": "RZ",
                    "framework_or_index_contrast": as_float(positive_row, full_column) - as_float(rz, full_column),
                    "simple_baseline_contrast": as_float(positive_row, baseline_column) - as_float(rz, baseline_column),
                    "framework_direction_supported": (
                        as_float(positive_row, full_column) - as_float(rz, full_column)
                    )
                    > 0,
                }
            )

    write_tsv(TABLE_ROOT / "cross_study_model_benchmark.tsv", benchmark_rows)
    write_tsv(TABLE_ROOT / "cross_study_model_benchmark_gain.tsv", improvement_rows)
    write_tsv(TABLE_ROOT / "cross_study_boundary_transition_index.tsv", transition_rows)
    write_tsv(TABLE_ROOT / "external_author_area_gradient_benchmark.tsv", area_gradient_rows)
    print(
        f"benchmarks={len(benchmark_rows)} comparisons={len(improvement_rows)} "
        f"sections={len(transition_rows)} author_gradients={len(area_gradient_rows)}"
    )


if __name__ == "__main__":
    main()
