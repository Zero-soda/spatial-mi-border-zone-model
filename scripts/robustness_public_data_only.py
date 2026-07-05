#!/usr/bin/env python3
"""Public-data-only robustness checks for the MI border-zone model."""

from __future__ import annotations

import csv
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean

import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[3]
TABLE_DIR = ROOT / "results" / "tables"
FIGURE_DIR = ROOT / "results" / "figures"

MOUSE_SAMPLE_DOMAIN = TABLE_DIR / "gse214611_d3_d7_signature_scores_by_sample_domain.tsv"
MOUSE_SPOT = TABLE_DIR / "gse214611_d3_d7_signature_scores_by_spot.tsv"
BOUNDARY_SAMPLE = TABLE_DIR / "gse214611_d3_d7_domain34_spatial_features_by_sample.tsv"
HUMAN_SPOT = TABLE_DIR / "gse214611_human_stemi_signature_scores_by_spot.tsv"

RANDOM_SEED = 614611
RANDOM_ITERATIONS = 1000
BOUNDARY_THRESHOLD_MULTIPLIERS = [1.10, 1.20, 1.35, 1.50, 1.75]

MECHANICAL_COMPONENTS = ["z_CM_BZ1_TRANSITION", "z_CM_BZ2_MECHANICAL_EDGE"]
IMMUNE_FIBROTIC_COMPONENTS = [
    "z_CCR2_IL1B_MYLOID",
    "z_TGFB_SIGNALING",
    "z_FAP_POSTN_PATHO_FIBROBLAST",
]
SCAR_REPAIR_COMPONENTS = [
    "z_FAP_POSTN_PATHO_FIBROBLAST",
    "z_ECM_REMODELING",
    "z_CTHRC1_REPARATIVE_CF",
    "z_MYOFIBROBLAST_CONTRACTILE",
]

TRUE_AXIS_CONFIGS = [
    {
        "check": "day3_domain3_mechanical_exceeds_domain4",
        "axis": "mechanical_border",
        "components": MECHANICAL_COMPONENTS,
        "stage": "day3_mi",
        "direction": "domain3_minus_domain4",
    },
    {
        "check": "day7_domain3_mechanical_exceeds_domain4",
        "axis": "mechanical_border",
        "components": MECHANICAL_COMPONENTS,
        "stage": "day7_mi",
        "direction": "domain3_minus_domain4",
    },
    {
        "check": "day7_domain4_fibroblast_scar_exceeds_domain3",
        "axis": "fibroblast_scar_repair",
        "components": SCAR_REPAIR_COMPONENTS,
        "stage": "day7_mi",
        "direction": "domain4_minus_domain3",
    },
    {
        "check": "day7_domain4_immune_fibrotic_exceeds_domain3",
        "axis": "immune_fibrotic_activation",
        "components": IMMUNE_FIBROTIC_COMPONENTS,
        "stage": "day7_mi",
        "direction": "domain4_minus_domain3",
    },
]


def safe_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def group_mean(rows: list[dict[str, str]], group_cols: list[str], value_col: str) -> dict[tuple[str, ...], float]:
    sums: dict[tuple[str, ...], float] = defaultdict(float)
    counts: dict[tuple[str, ...], int] = defaultdict(int)
    for row in rows:
        key = tuple(row[col] for col in group_cols)
        sums[key] += safe_float(row[value_col])
        counts[key] += 1
    return {key: sums[key] / counts[key] for key in counts}


def leave_one_sample_robustness(sample_domain_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    samples = sorted({row["sample"] for row in sample_domain_rows})
    rows: list[dict[str, object]] = []
    for omitted in ["none", *samples]:
        kept = [row for row in sample_domain_rows if omitted == "none" or row["sample"] != omitted]
        means = {
            (row["sample"], row["annotated"]): row
            for row in kept
        }

        def stage_values(stage_prefix: str, metric: str, domain_a: str, domain_b: str, direction: str) -> list[float]:
            values: list[float] = []
            for sample in sorted({row["sample"] for row in kept if row["sample"].startswith(stage_prefix)}):
                a = safe_float(means[(sample, domain_a)][metric])
                b = safe_float(means[(sample, domain_b)][metric])
                values.append(a - b if direction == "a_minus_b" else b - a)
            return values

        checks = [
            (
                "day3_domain3_mechanical_exceeds_domain4",
                stage_values("D3_", "mean_signature_mechanical_border_score", "3", "4", "a_minus_b"),
            ),
            (
                "day7_domain3_mechanical_exceeds_domain4",
                stage_values("D7_", "mean_signature_mechanical_border_score", "3", "4", "a_minus_b"),
            ),
            (
                "day7_domain4_scar_repair_exceeds_domain3",
                stage_values("D7_", "mean_signature_fibroblast_scar_score", "3", "4", "b_minus_a"),
            ),
            (
                "day3_domain4_scar_repair_not_yet_dominant",
                [-value for value in stage_values("D3_", "mean_signature_fibroblast_scar_score", "3", "4", "b_minus_a")],
            ),
        ]
        for check, values in checks:
            rows.append(
                {
                    "omitted_sample": omitted,
                    "check": check,
                    "n_remaining_values": len(values),
                    "mean_margin": f"{mean(values):.6g}" if values else "NA",
                    "min_margin": f"{min(values):.6g}" if values else "NA",
                    "direction_preserved": bool(values and min(values) > 0),
                }
            )
    return rows


def add_module_dropout_scores(rows: list[dict[str, str]], components: list[str], out_col: str, drop_component: str | None) -> None:
    active = [component for component in components if component != drop_component]
    for row in rows:
        row[out_col] = f"{sum(safe_float(row[component]) for component in active):.8g}"


def module_dropout_robustness(spot_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    configs = [
        (
            "mechanical_border",
            MECHANICAL_COMPONENTS,
            "domain3_minus_domain4",
            "day3_mi",
            "3",
            "4",
        ),
        (
            "mechanical_border",
            MECHANICAL_COMPONENTS,
            "domain3_minus_domain4",
            "day7_mi",
            "3",
            "4",
        ),
        (
            "fibroblast_scar_repair",
            SCAR_REPAIR_COMPONENTS,
            "domain4_minus_domain3",
            "day7_mi",
            "3",
            "4",
        ),
        (
            "immune_fibrotic_activation",
            IMMUNE_FIBROTIC_COMPONENTS,
            "domain4_minus_domain3",
            "day7_mi",
            "3",
            "4",
        ),
    ]

    out_rows: list[dict[str, object]] = []
    for axis, components, direction, stage, domain3, domain4 in configs:
        for drop_component in [None, *components]:
            working = [dict(row) for row in spot_rows if row["stage"] == stage and row["annotated"] in {domain3, domain4}]
            score_col = f"dropout_{axis}"
            add_module_dropout_scores(working, components, score_col, drop_component)
            sample_domain = group_mean(working, ["sample", "annotated"], score_col)
            margins = []
            for sample in sorted({row["sample"] for row in working}):
                d3 = sample_domain[(sample, domain3)]
                d4 = sample_domain[(sample, domain4)]
                margins.append(d3 - d4 if direction == "domain3_minus_domain4" else d4 - d3)

            out_rows.append(
                {
                    "axis": axis,
                    "stage": stage,
                    "direction_test": direction,
                    "dropped_component": drop_component or "none",
                    "n_samples": len(margins),
                    "mean_margin": f"{mean(margins):.6g}",
                    "min_margin": f"{min(margins):.6g}",
                    "direction_preserved": min(margins) > 0,
                }
            )
    return out_rows


def rank(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    idx = 0
    while idx < len(indexed):
        end = idx
        while end + 1 < len(indexed) and indexed[end + 1][1] == indexed[idx][1]:
            end += 1
        avg_rank = (idx + end + 2) / 2.0
        for pos in range(idx, end + 1):
            ranks[indexed[pos][0]] = avg_rank
        idx = end + 1
    return ranks


def pearson(x: list[float], y: list[float]) -> float:
    mx = mean(x)
    my = mean(y)
    num = sum((a - mx) * (b - my) for a, b in zip(x, y, strict=True))
    den_x = math.sqrt(sum((a - mx) ** 2 for a in x))
    den_y = math.sqrt(sum((b - my) ** 2 for b in y))
    return num / (den_x * den_y) if den_x and den_y else float("nan")


def human_output_separation(human_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    columns = [
        "mechanical_border_score",
        "immune_fibrotic_activation_score",
        "fibroblast_scar_repair_score",
    ]
    out: list[dict[str, object]] = []
    for left_idx, left in enumerate(columns):
        for right in columns[left_idx + 1 :]:
            x = [safe_float(row[left]) for row in human_rows]
            y = [safe_float(row[right]) for row in human_rows]
            cutoff_x = sorted(x)[int(len(x) * 0.90)]
            cutoff_y = sorted(y)[int(len(y) * 0.90)]
            top_x = {idx for idx, value in enumerate(x) if value >= cutoff_x}
            top_y = {idx for idx, value in enumerate(y) if value >= cutoff_y}
            jaccard = len(top_x & top_y) / len(top_x | top_y)
            out.append(
                {
                    "left_score": left,
                    "right_score": right,
                    "pearson_r": f"{pearson(x, y):.6g}",
                    "spearman_r": f"{pearson(rank(x), rank(y)):.6g}",
                    "top_decile_jaccard": f"{jaccard:.6g}",
                    "top_decile_overlap_spots": len(top_x & top_y),
                    "top_decile_union_spots": len(top_x | top_y),
                }
            )
    return out


def component_mean_lookup(spot_rows: list[dict[str, str]], component_cols: list[str]) -> dict[tuple[str, str, str, str], float]:
    sums: dict[tuple[str, str, str, str], float] = defaultdict(float)
    counts: dict[tuple[str, str, str, str], int] = defaultdict(int)
    for row in spot_rows:
        if row["annotated"] not in {"3", "4"}:
            continue
        for component in component_cols:
            key = (row["stage"], row["sample"], row["annotated"], component)
            sums[key] += safe_float(row[component])
            counts[key] += 1
    return {key: sums[key] / counts[key] for key in counts}


def stage_domain_margin(
    component_means: dict[tuple[str, str, str, str], float],
    components: list[str],
    stage: str,
    direction: str,
) -> float:
    samples = sorted({sample for stage_name, sample, _, _ in component_means if stage_name == stage})
    margins = []
    for sample in samples:
        domain3 = sum(component_means[(stage, sample, "3", component)] for component in components)
        domain4 = sum(component_means[(stage, sample, "4", component)] for component in components)
        margins.append(domain3 - domain4 if direction == "domain3_minus_domain4" else domain4 - domain3)
    return mean(margins)


def percentile_from_sorted(values: list[float], q: float) -> float:
    if not values:
        return float("nan")
    idx = min(len(values) - 1, max(0, int(round((len(values) - 1) * q))))
    return values[idx]


def random_signature_negative_control(spot_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    component_cols = [column for column in spot_rows[0] if column.startswith("z_")]
    component_means = component_mean_lookup(spot_rows, component_cols)
    rng = random.Random(RANDOM_SEED)
    out: list[dict[str, object]] = []
    for config in TRUE_AXIS_CONFIGS:
        components = list(config["components"])
        observed = stage_domain_margin(component_means, components, config["stage"], config["direction"])
        random_pool = [component for component in component_cols if component not in components]
        random_margins = []
        example_components = ""
        for iteration in range(RANDOM_ITERATIONS):
            sampled = rng.sample(random_pool, len(components))
            if iteration == 0:
                example_components = ";".join(sampled)
            random_margins.append(stage_domain_margin(component_means, sampled, config["stage"], config["direction"]))
        random_margins.sort()
        random_p95 = percentile_from_sorted(random_margins, 0.95)
        random_p99 = percentile_from_sorted(random_margins, 0.99)
        fraction_at_or_above_observed = sum(value >= observed for value in random_margins) / len(random_margins)
        out.append(
            {
                "check": config["check"],
                "axis": config["axis"],
                "stage": config["stage"],
                "direction_test": config["direction"],
                "true_components": ";".join(components),
                "n_true_components": len(components),
                "n_random_iterations": RANDOM_ITERATIONS,
                "observed_mean_margin": f"{observed:.6g}",
                "random_mean_margin": f"{mean(random_margins):.6g}",
                "random_p95_margin": f"{random_p95:.6g}",
                "random_p99_margin": f"{random_p99:.6g}",
                "fraction_random_ge_observed": f"{fraction_at_or_above_observed:.6g}",
                "observed_exceeds_random_p95": observed > random_p95,
                "observed_exceeds_random_p99": observed > random_p99,
                "example_random_components": example_components,
            }
        )
    return out


def pairwise_distances(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    delta = left[:, None, :] - right[None, :, :]
    return np.sqrt(np.sum(delta * delta, axis=2))


def median_nearest_neighbor_distance(coords: np.ndarray) -> float:
    distances = pairwise_distances(coords, coords)
    np.fill_diagonal(distances, np.inf)
    return float(np.median(np.min(distances, axis=1)))


def boundary_threshold_sensitivity(spot_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in spot_rows:
        grouped[row["sample"]].append(row)

    sample_rows: list[dict[str, object]] = []
    for sample, rows in sorted(grouped.items()):
        coords = np.array(
            [
                [safe_float(row["pxl_col_in_fullres"]), safe_float(row["pxl_row_in_fullres"])]
                for row in rows
            ],
            dtype=float,
        )
        domains = np.array([int(row["annotated"]) for row in rows], dtype=int)
        domain3_indices = np.where(domains == 3)[0]
        domain4_indices = np.where(domains == 4)[0]
        if not len(domain3_indices) or not len(domain4_indices):
            continue
        domain3_coords = coords[domain3_indices]
        domain4_coords = coords[domain4_indices]
        distances_34 = pairwise_distances(domain3_coords, domain4_coords)
        median_nn = median_nearest_neighbor_distance(coords)

        for multiplier in BOUNDARY_THRESHOLD_MULTIPLIERS:
            threshold = median_nn * multiplier
            contact_pairs = np.argwhere(distances_34 <= threshold)
            domain3_boundary = np.any(distances_34 <= threshold, axis=1)
            domain4_boundary = np.any(distances_34 <= threshold, axis=0)
            mech_deltas = []
            scar_deltas = []
            for domain3_pos, domain4_pos in contact_pairs:
                row3 = rows[int(domain3_indices[int(domain3_pos)])]
                row4 = rows[int(domain4_indices[int(domain4_pos)])]
                mech_deltas.append(
                    safe_float(row4["signature_mechanical_border_score"])
                    - safe_float(row3["signature_mechanical_border_score"])
                )
                scar_deltas.append(
                    safe_float(row4["signature_fibroblast_scar_score"])
                    - safe_float(row3["signature_fibroblast_scar_score"])
                )
            sample_rows.append(
                {
                    "sample": sample,
                    "stage": rows[0]["stage"],
                    "threshold_multiplier": f"{multiplier:.2f}",
                    "contact_threshold_fullres_px": f"{threshold:.6g}",
                    "n_contact_pairs": len(contact_pairs),
                    "fraction_domain3_boundary": f"{float(np.mean(domain3_boundary)):.6g}",
                    "fraction_domain4_boundary": f"{float(np.mean(domain4_boundary)):.6g}",
                    "contact_delta_signature_mechanical_border_score_domain4_minus_domain3": (
                        f"{mean(mech_deltas):.6g}" if mech_deltas else "nan"
                    ),
                    "contact_delta_signature_fibroblast_scar_score_domain4_minus_domain3": (
                        f"{mean(scar_deltas):.6g}" if scar_deltas else "nan"
                    ),
                }
            )

    out: list[dict[str, object]] = []
    grouped_stage: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in sample_rows:
        grouped_stage[(str(row["threshold_multiplier"]), str(row["stage"]))].append(row)

    for (threshold_multiplier, stage), rows in sorted(grouped_stage.items()):
        mech_values = [safe_float(str(row["contact_delta_signature_mechanical_border_score_domain4_minus_domain3"])) for row in rows]
        scar_values = [safe_float(str(row["contact_delta_signature_fibroblast_scar_score_domain4_minus_domain3"])) for row in rows]
        expected_scar = "positive" if stage == "day7_mi" else "negative"
        out.append(
            {
                "threshold_multiplier": threshold_multiplier,
                "stage": stage,
                "n_samples": len(rows),
                "mean_contact_pairs": f"{mean(safe_float(str(row['n_contact_pairs'])) for row in rows):.6g}",
                "mean_fraction_domain3_boundary": f"{mean(safe_float(str(row['fraction_domain3_boundary'])) for row in rows):.6g}",
                "mean_fraction_domain4_boundary": f"{mean(safe_float(str(row['fraction_domain4_boundary'])) for row in rows):.6g}",
                "mean_contact_delta_mechanical_domain4_minus_domain3": f"{mean(mech_values):.6g}",
                "mechanical_expected_direction": "negative",
                "mechanical_direction_preserved_all_samples": all(value < 0 for value in mech_values),
                "mean_contact_delta_scar_domain4_minus_domain3": f"{mean(scar_values):.6g}",
                "scar_expected_direction": expected_scar,
                "scar_direction_preserved_all_samples": (
                    all(value > 0 for value in scar_values)
                    if expected_scar == "positive"
                    else all(value < 0 for value in scar_values)
                ),
            }
        )
    return out


def boundary_direction_summary(boundary_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    checks = [
        ("day3_mi", "contact_delta_signature_mechanical_border_score_domain4_minus_domain3", "negative"),
        ("day7_mi", "contact_delta_signature_mechanical_border_score_domain4_minus_domain3", "negative"),
        ("day3_mi", "contact_delta_signature_fibroblast_scar_score_domain4_minus_domain3", "negative"),
        ("day7_mi", "contact_delta_signature_fibroblast_scar_score_domain4_minus_domain3", "positive"),
    ]
    out = []
    for stage, metric, expected in checks:
        values = [safe_float(row[metric]) for row in boundary_rows if row["stage"] == stage]
        if expected == "positive":
            preserved = all(value > 0 for value in values)
        else:
            preserved = all(value < 0 for value in values)
        out.append(
            {
                "stage": stage,
                "metric": metric,
                "expected_direction": expected,
                "n_samples": len(values),
                "mean_value": f"{mean(values):.6g}",
                "min_value": f"{min(values):.6g}",
                "max_value": f"{max(values):.6g}",
                "all_samples_preserved": preserved,
            }
        )
    return out


def write_summary_figure(
    loso_rows: list[dict[str, object]],
    dropout_rows: list[dict[str, object]],
    boundary_rows: list[dict[str, object]],
    human_rows: list[dict[str, object]],
) -> Path:
    width, height = 1400, 820
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    draw.text((40, 30), "Public-data-only robustness summary", fill=(20, 20, 20))

    sections = [
        ("Leave-one-sample-out", loso_rows),
        ("Module dropout", dropout_rows),
        ("Domain 3/4 boundary", boundary_rows),
        ("Human STEMI score separation", human_rows),
    ]
    y = 90
    for title, rows in sections:
        draw.text((40, y), title, fill=(31, 77, 120))
        y += 28
        if title == "Human STEMI score separation":
            for row in rows:
                text = (
                    f"{row['left_score']} vs {row['right_score']}: "
                    f"Spearman r={row['spearman_r']}, top-decile Jaccard={row['top_decile_jaccard']}"
                )
                draw.text((70, y), text, fill=(40, 40, 40))
                y += 24
        else:
            preserved = sum(
                bool(row.get("direction_preserved", row.get("all_samples_preserved", False)))
                for row in rows
            )
            total = len(rows)
            draw.text((70, y), f"Direction-preserved checks: {preserved}/{total}", fill=(40, 40, 40))
            y += 24
            for row in rows[:6]:
                label = row.get("check") or row.get("axis") or row.get("metric")
                margin = row.get("mean_margin", row.get("mean_value", "NA"))
                draw.text((90, y), f"{label}: mean={margin}", fill=(80, 80, 80))
                y += 22
        y += 22

    out_path = FIGURE_DIR / "gse214611_public_only_robustness_summary.png"
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    image.save(out_path)
    return out_path


def main() -> None:
    sample_domain_rows = read_tsv(MOUSE_SAMPLE_DOMAIN)
    spot_rows = read_tsv(MOUSE_SPOT)
    boundary_rows = read_tsv(BOUNDARY_SAMPLE)
    human_rows = read_tsv(HUMAN_SPOT)

    loso = leave_one_sample_robustness(sample_domain_rows)
    dropout = module_dropout_robustness(spot_rows)
    boundary = boundary_direction_summary(boundary_rows)
    human = human_output_separation(human_rows)
    random_control = random_signature_negative_control(spot_rows)
    threshold_sensitivity = boundary_threshold_sensitivity(spot_rows)

    write_tsv(
        TABLE_DIR / "gse214611_loso_stage_domain_robustness.tsv",
        loso,
        ["omitted_sample", "check", "n_remaining_values", "mean_margin", "min_margin", "direction_preserved"],
    )
    write_tsv(
        TABLE_DIR / "gse214611_module_dropout_robustness.tsv",
        dropout,
        ["axis", "stage", "direction_test", "dropped_component", "n_samples", "mean_margin", "min_margin", "direction_preserved"],
    )
    write_tsv(
        TABLE_DIR / "gse214611_boundary_direction_robustness.tsv",
        boundary,
        ["stage", "metric", "expected_direction", "n_samples", "mean_value", "min_value", "max_value", "all_samples_preserved"],
    )
    write_tsv(
        TABLE_DIR / "gse214611_human_stemi_score_separation.tsv",
        human,
        ["left_score", "right_score", "pearson_r", "spearman_r", "top_decile_jaccard", "top_decile_overlap_spots", "top_decile_union_spots"],
    )
    write_tsv(
        TABLE_DIR / "gse214611_random_signature_negative_control.tsv",
        random_control,
        [
            "check",
            "axis",
            "stage",
            "direction_test",
            "true_components",
            "n_true_components",
            "n_random_iterations",
            "observed_mean_margin",
            "random_mean_margin",
            "random_p95_margin",
            "random_p99_margin",
            "fraction_random_ge_observed",
            "observed_exceeds_random_p95",
            "observed_exceeds_random_p99",
            "example_random_components",
        ],
    )
    write_tsv(
        TABLE_DIR / "gse214611_boundary_threshold_sensitivity.tsv",
        threshold_sensitivity,
        [
            "threshold_multiplier",
            "stage",
            "n_samples",
            "mean_contact_pairs",
            "mean_fraction_domain3_boundary",
            "mean_fraction_domain4_boundary",
            "mean_contact_delta_mechanical_domain4_minus_domain3",
            "mechanical_expected_direction",
            "mechanical_direction_preserved_all_samples",
            "mean_contact_delta_scar_domain4_minus_domain3",
            "scar_expected_direction",
            "scar_direction_preserved_all_samples",
        ],
    )
    figure = write_summary_figure(loso, dropout, boundary, human)

    print(f"[loso] {sum(row['direction_preserved'] for row in loso)}/{len(loso)} checks preserved")
    print(f"[module-dropout] {sum(row['direction_preserved'] for row in dropout)}/{len(dropout)} checks preserved")
    print(f"[boundary] {sum(row['all_samples_preserved'] for row in boundary)}/{len(boundary)} checks preserved")
    print(f"[human] {len(human)} pairwise score separation rows")
    print(
        f"[random-control] "
        f"{sum(row['observed_exceeds_random_p95'] for row in random_control)}/{len(random_control)} true axes exceed random p95"
    )
    print(
        f"[threshold] "
        f"{sum(row['mechanical_direction_preserved_all_samples'] and row['scar_direction_preserved_all_samples'] for row in threshold_sensitivity)}/{len(threshold_sensitivity)} stage-threshold checks preserved"
    )
    print(f"[figure] {figure}")


if __name__ == "__main__":
    main()
