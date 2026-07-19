#!/usr/bin/env python3
"""Additional independence checks for the spatial MI manuscript.

This script addresses two remaining reviewer-facing questions without using
new local tissue or outcome data:

1. Do score-derived states remain spatially coherent when author-provided
   domain labels are excluded from state construction?
2. Are the spatial patterns in the single human STEMI section stronger than
   patterns obtained from expression- and detection-matched random genes?

The analyses are sensitivity and falsification checks. They do not constitute
an independent biological or clinical validation cohort.
"""

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
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import numpy as np  # noqa: E402

from batch_map_gse214611_mi_spatial_risk import locate_visium_dir  # noqa: E402
from final_reviewer_upgrade_analyses import (  # noqa: E402
    HUMAN_SAMPLE,
    SIGNATURE_COMPONENTS,
    find_matrix_h5,
    gene_stats,
    load_dense_lognorm,
    matched_control_pool,
    module_mean,
    read_10x_h5,
    signature_indices,
    zscore,
)
from reviewer_rigor_score_state_celltype_audit import (  # noqa: E402
    adjusted_rand_index,
    best_domain_recovery,
    deterministic_kmeans,
    normalized_mutual_information,
)
from spatial_model_utils import read_tsv, safe_float  # noqa: E402


TABLE_DIR = ROOT / "results" / "tables"
SIGNATURES_TSV = Path(__file__).resolve().parents[1] / "config" / "spatial_cardiac_border_zone_signatures.tsv"
MOUSE_SPOT_TABLE = TABLE_DIR / "gse214611_d3_d7_signature_scores_by_spot.tsv"
HUMAN_SPOT_TABLE = TABLE_DIR / "gse214611_human_stemi_signature_scores_by_spot.tsv"

OUT_GRAPH_RECOVERY = TABLE_DIR / "gse214611_graph_smoothed_score_state_recovery.tsv"
OUT_GRAPH_SPOTS = TABLE_DIR / "gse214611_graph_smoothed_score_states_by_spot.tsv"
OUT_HUMAN_NULL_SUMMARY = TABLE_DIR / "gse214611_human_expression_matched_spatial_null.tsv"
OUT_HUMAN_NULL_DISTRIBUTION = TABLE_DIR / "gse214611_human_expression_matched_spatial_null_distribution.tsv"

RANDOM_SEED = 20260710
N_RANDOM = 500
PRIMARY_ALPHA = 0.35
ALPHAS = [0.0, 0.20, PRIMARY_ALPHA, 0.50]

SCORE_COLUMNS = [
    "signature_mechanical_border_score",
    "reviewer_immune_fibrotic_activation_score",
    "signature_fibroblast_scar_score",
]

OUTPUTS = {
    "mechanical_border": "mechanical_border_score",
    "immune_fibrotic_activation": "immune_fibrotic_activation_score",
    "fibroblast_scar_repair": "fibroblast_scar_repair_score",
}


def format_float(value: float, precision: int = 8) -> str:
    if not math.isfinite(value):
        return "nan"
    return f"{value:.{precision}g}"


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_human_signatures(path: Path) -> dict[str, list[str]]:
    signatures: dict[str, list[str]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            signatures[row["signature_id"]] = [
                gene.strip() for gene in row["genes_human"].split(";") if gene.strip()
            ]
    return signatures


def add_immune_score(row: dict[str, str]) -> None:
    value = (
        safe_float(row["z_CCR2_IL1B_MYLOID"])
        + safe_float(row["z_TGFB_SIGNALING"])
        + safe_float(row["z_FAP_POSTN_PATHO_FIBROBLAST"])
    )
    row["reviewer_immune_fibrotic_activation_score"] = format_float(value)


def build_knn_edges(coords: np.ndarray, k: int = 6) -> list[tuple[int, int]]:
    coords = np.asarray(coords, dtype=float)
    distances = np.sum((coords[:, None, :] - coords[None, :, :]) ** 2, axis=2)
    np.fill_diagonal(distances, np.inf)
    nearest = np.argpartition(distances, kth=min(k, len(coords) - 1) - 1, axis=1)[:, :k]
    edges: set[tuple[int, int]] = set()
    for left, neighbours in enumerate(nearest):
        for right in neighbours:
            if left == int(right):
                continue
            edges.add((min(left, int(right)), max(left, int(right))))
    return sorted(edges)


def neighbour_lists(n_spots: int, edges: list[tuple[int, int]]) -> list[list[int]]:
    neighbours = [[] for _ in range(n_spots)]
    for left, right in edges:
        neighbours[left].append(right)
        neighbours[right].append(left)
    return neighbours


def graph_smooth(x: np.ndarray, neighbours: list[list[int]], alpha: float) -> np.ndarray:
    if alpha == 0:
        return x.copy()
    smoothed = np.empty_like(x, dtype=float)
    for idx, local_neighbours in enumerate(neighbours):
        if not local_neighbours:
            smoothed[idx] = x[idx]
            continue
        local_mean = np.mean(x[local_neighbours], axis=0)
        smoothed[idx] = (1.0 - alpha) * x[idx] + alpha * local_mean
    return smoothed


def same_state_edge_fraction(labels: np.ndarray, edges: list[tuple[int, int]]) -> float:
    if not edges:
        return float("nan")
    return float(np.mean([labels[left] == labels[right] for left, right in edges]))


def state_domain_metrics(
    author_domains: list[str], labels: np.ndarray, cluster_id: int, target_domain: str
) -> tuple[float, float, float]:
    member = labels == cluster_id
    target = np.asarray(author_domains) == target_domain
    true_positive = int(np.sum(member & target))
    precision = true_positive / int(np.sum(member)) if np.any(member) else 0.0
    recall = true_positive / int(np.sum(target)) if np.any(target) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return precision, recall, f1


def graph_smoothed_score_states() -> None:
    rows = read_tsv(MOUSE_SPOT_TABLE)
    for row in rows:
        add_immune_score(row)
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["sample"]].append(row)

    summary_rows: list[dict[str, object]] = []
    spot_rows: list[dict[str, object]] = []

    for sample, sample_rows in sorted(grouped.items()):
        x = np.asarray(
            [[safe_float(row[column]) for column in SCORE_COLUMNS] for row in sample_rows],
            dtype=float,
        )
        x = np.column_stack([zscore(x[:, idx]) for idx in range(x.shape[1])])
        coords = np.asarray(
            [[safe_float(row["array_row"]), safe_float(row["array_col"])] for row in sample_rows],
            dtype=float,
        )
        edges = build_knn_edges(coords, k=6)
        neighbours = neighbour_lists(len(sample_rows), edges)
        author_domains = [row["annotated"] for row in sample_rows]

        for alpha in ALPHAS:
            features = graph_smooth(x, neighbours, alpha)
            labels, centres, inertia = deterministic_kmeans(
                features,
                k=4,
                seed=RANDOM_SEED + int(round(alpha * 1000)) + sum(ord(char) for char in sample),
            )
            mechanical_cluster = int(np.argmax(centres[:, 0]))
            scar_cluster = int(np.argmax(centres[:, 2]))
            immune_cluster = int(np.argmax(centres[:, 1]))
            d3 = best_domain_recovery(author_domains, labels, "3")
            d4 = best_domain_recovery(author_domains, labels, "4")
            mech_precision, mech_recall, mech_f1 = state_domain_metrics(
                author_domains, labels, mechanical_cluster, "3"
            )
            scar_precision, scar_recall, scar_f1 = state_domain_metrics(
                author_domains, labels, scar_cluster, "4"
            )
            mech_scar_distinct = mechanical_cluster != scar_cluster
            mech_scar_edges = (
                sum(
                    {int(labels[left]), int(labels[right])} == {mechanical_cluster, scar_cluster}
                    for left, right in edges
                )
                if mech_scar_distinct
                else 0
            )
            summary_rows.append(
                {
                    "sample": sample,
                    "stage": sample_rows[0]["stage"],
                    "graph_smoothing_alpha": alpha,
                    "n_spots": len(sample_rows),
                    "n_graph_edges": len(edges),
                    "ari_vs_author_domains": format_float(adjusted_rand_index(author_domains, labels.tolist())),
                    "nmi_vs_author_domains": format_float(normalized_mutual_information(author_domains, labels.tolist())),
                    "same_state_edge_fraction": format_float(same_state_edge_fraction(labels, edges)),
                    "kmeans_inertia": format_float(inertia),
                    "domain3_best_cluster_f1": format_float(float(d3["f1"])),
                    "domain4_best_cluster_f1": format_float(float(d4["f1"])),
                    "mechanical_high_cluster": mechanical_cluster,
                    "mechanical_high_domain3_precision": format_float(mech_precision),
                    "mechanical_high_domain3_recall": format_float(mech_recall),
                    "mechanical_high_domain3_f1": format_float(mech_f1),
                    "scar_high_cluster": scar_cluster,
                    "scar_high_domain4_precision": format_float(scar_precision),
                    "scar_high_domain4_recall": format_float(scar_recall),
                    "scar_high_domain4_f1": format_float(scar_f1),
                    "immune_high_cluster": immune_cluster,
                    "mechanical_scar_state_separable": int(mech_scar_distinct),
                    "mechanical_scar_graph_edges": mech_scar_edges,
                }
            )

            if alpha == PRIMARY_ALPHA:
                state_names: dict[int, str] = {}
                for cluster_id in range(4):
                    roles = []
                    if cluster_id == mechanical_cluster:
                        roles.append("mechanical")
                    if cluster_id == scar_cluster:
                        roles.append("scar")
                    if cluster_id == immune_cluster:
                        roles.append("immune")
                    state_names[cluster_id] = "/".join(roles) + "-high" if roles else "mixed/low"
                for idx, (row, cluster_id) in enumerate(zip(sample_rows, labels, strict=True)):
                    spot_rows.append(
                        {
                            "sample": sample,
                            "stage": row["stage"],
                            "barcode": row["barcode"],
                            "array_row": row["array_row"],
                            "array_col": row["array_col"],
                            "pxl_row_in_fullres": row["pxl_row_in_fullres"],
                            "pxl_col_in_fullres": row["pxl_col_in_fullres"],
                            "author_domain_audit_only": row["annotated"],
                            "graph_smoothed_cluster_k4": int(cluster_id),
                            "graph_smoothed_state": state_names[int(cluster_id)],
                            "smoothed_mechanical_score": format_float(float(features[idx, 0])),
                            "smoothed_immune_fibrotic_score": format_float(float(features[idx, 1])),
                            "smoothed_fibroblast_scar_score": format_float(float(features[idx, 2])),
                        }
                    )

    write_tsv(
        OUT_GRAPH_RECOVERY,
        summary_rows,
        [
            "sample",
            "stage",
            "graph_smoothing_alpha",
            "n_spots",
            "n_graph_edges",
            "ari_vs_author_domains",
            "nmi_vs_author_domains",
            "same_state_edge_fraction",
            "kmeans_inertia",
            "domain3_best_cluster_f1",
            "domain4_best_cluster_f1",
            "mechanical_high_cluster",
            "mechanical_high_domain3_precision",
            "mechanical_high_domain3_recall",
            "mechanical_high_domain3_f1",
            "scar_high_cluster",
            "scar_high_domain4_precision",
            "scar_high_domain4_recall",
            "scar_high_domain4_f1",
            "immune_high_cluster",
            "mechanical_scar_state_separable",
            "mechanical_scar_graph_edges",
        ],
    )
    write_tsv(
        OUT_GRAPH_SPOTS,
        spot_rows,
        [
            "sample",
            "stage",
            "barcode",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
            "author_domain_audit_only",
            "graph_smoothed_cluster_k4",
            "graph_smoothed_state",
            "smoothed_mechanical_score",
            "smoothed_immune_fibrotic_score",
            "smoothed_fibroblast_scar_score",
        ],
    )


def graph_morans_i(values: np.ndarray, edges: list[tuple[int, int]]) -> float:
    values = np.asarray(values, dtype=float)
    centred = values - np.mean(values)
    denominator = float(np.sum(centred**2))
    if denominator == 0 or not edges:
        return float("nan")
    numerator = 2.0 * sum(centred[left] * centred[right] for left, right in edges)
    weight_sum = 2.0 * len(edges)
    return float((len(values) / weight_sum) * (numerator / denominator))


def hotspot_neighbour_fraction(values: np.ndarray, edges: list[tuple[int, int]], quantile: float = 0.90) -> float:
    threshold = float(np.quantile(values, quantile))
    hotspot = np.asarray(values >= threshold, dtype=bool)
    neighbours = neighbour_lists(len(values), edges)
    hotspot_indices = np.flatnonzero(hotspot)
    if len(hotspot_indices) == 0:
        return float("nan")
    connected = [any(hotspot[other] for other in neighbours[idx]) for idx in hotspot_indices]
    return float(np.mean(connected))


def human_expression_matched_spatial_null() -> None:
    signatures = read_human_signatures(SIGNATURES_TSV)
    visium_dir = locate_visium_dir(HUMAN_SAMPLE.unpack_dir)
    matrix_data = read_10x_h5(find_matrix_h5(visium_dir))
    dense_all = load_dense_lognorm(matrix_data)
    matrix_lookup = {barcode: idx for idx, barcode in enumerate(matrix_data.barcodes)}
    human_rows = read_tsv(HUMAN_SPOT_TABLE)
    keep_indices = np.asarray([matrix_lookup[row["barcode"]] for row in human_rows], dtype=int)
    dense = dense_all[keep_indices]
    coords = np.asarray(
        [[safe_float(row["array_row"]), safe_float(row["array_col"])] for row in human_rows],
        dtype=float,
    )
    edges = build_knn_edges(coords, k=6)
    signature_idx = signature_indices(signatures, matrix_data)
    all_signature_indices = {idx for indices in signature_idx.values() for idx in indices}
    mean_expression, detection = gene_stats(dense_all)

    observed_modules = {
        signature_id: module_mean(dense, indices) for signature_id, indices in signature_idx.items()
    }
    existing_by_barcode = {row["barcode"]: row for row in human_rows}
    summary_rows: list[dict[str, object]] = []
    null_rows: list[dict[str, object]] = []

    for output_name, existing_column in OUTPUTS.items():
        components = SIGNATURE_COMPONENTS[output_name]
        observed = np.sum(
            np.vstack([zscore(observed_modules[component]) for component in components]),
            axis=0,
        )
        existing = np.asarray(
            [safe_float(existing_by_barcode[row["barcode"]][existing_column]) for row in human_rows],
            dtype=float,
        )
        observed_moran = graph_morans_i(observed, edges)
        observed_hotspot = hotspot_neighbour_fraction(observed, edges)
        score_reproduction_r = float(np.corrcoef(observed, existing)[0, 1])
        rng = np.random.default_rng(RANDOM_SEED + sum(ord(char) for char in output_name))
        null_moran: list[float] = []
        null_hotspot: list[float] = []

        for iteration in range(N_RANDOM):
            random_modules: dict[str, np.ndarray] = {}
            for component in components:
                controls = matched_control_pool(
                    signature_idx[component],
                    all_signature_indices,
                    mean_expression,
                    detection,
                    rng,
                    n_per_gene=1,
                )
                random_modules[component] = module_mean(dense, controls)
            random_score = np.sum(
                np.vstack([zscore(random_modules[component]) for component in components]),
                axis=0,
            )
            moran_value = graph_morans_i(random_score, edges)
            hotspot_value = hotspot_neighbour_fraction(random_score, edges)
            null_moran.append(moran_value)
            null_hotspot.append(hotspot_value)
            null_rows.append(
                {
                    "output": output_name,
                    "iteration": iteration + 1,
                    "graph_morans_i": format_float(moran_value),
                    "top_decile_hotspot_neighbour_fraction": format_float(hotspot_value),
                }
            )

        moran_arr = np.asarray(null_moran, dtype=float)
        hotspot_arr = np.asarray(null_hotspot, dtype=float)
        summary_rows.append(
            {
                "sample": "GSM6613090_V_Human_STEMI",
                "output": output_name,
                "n_spots": len(human_rows),
                "n_graph_edges": len(edges),
                "n_expression_matched_random_signatures": N_RANDOM,
                "score_reproduction_pearson_r": format_float(score_reproduction_r),
                "observed_graph_morans_i": format_float(observed_moran),
                "null_morans_i_p95": format_float(float(np.quantile(moran_arr, 0.95))),
                "null_morans_i_p99": format_float(float(np.quantile(moran_arr, 0.99))),
                "moran_empirical_upper_p": format_float(float((np.sum(moran_arr >= observed_moran) + 1) / (N_RANDOM + 1))),
                "observed_top_decile_hotspot_neighbour_fraction": format_float(observed_hotspot),
                "null_hotspot_neighbour_fraction_p95": format_float(float(np.quantile(hotspot_arr, 0.95))),
                "null_hotspot_neighbour_fraction_p99": format_float(float(np.quantile(hotspot_arr, 0.99))),
                "hotspot_empirical_upper_p": format_float(float((np.sum(hotspot_arr >= observed_hotspot) + 1) / (N_RANDOM + 1))),
                "interpretation": "single-section internal falsification; not independent human validation",
            }
        )

    write_tsv(
        OUT_HUMAN_NULL_SUMMARY,
        summary_rows,
        [
            "sample",
            "output",
            "n_spots",
            "n_graph_edges",
            "n_expression_matched_random_signatures",
            "score_reproduction_pearson_r",
            "observed_graph_morans_i",
            "null_morans_i_p95",
            "null_morans_i_p99",
            "moran_empirical_upper_p",
            "observed_top_decile_hotspot_neighbour_fraction",
            "null_hotspot_neighbour_fraction_p95",
            "null_hotspot_neighbour_fraction_p99",
            "hotspot_empirical_upper_p",
            "interpretation",
        ],
    )
    write_tsv(
        OUT_HUMAN_NULL_DISTRIBUTION,
        null_rows,
        [
            "output",
            "iteration",
            "graph_morans_i",
            "top_decile_hotspot_neighbour_fraction",
        ],
    )


def main() -> None:
    graph_smoothed_score_states()
    human_expression_matched_spatial_null()
    print(f"Wrote {OUT_GRAPH_RECOVERY}")
    print(f"Wrote {OUT_GRAPH_SPOTS}")
    print(f"Wrote {OUT_HUMAN_NULL_SUMMARY}")
    print(f"Wrote {OUT_HUMAN_NULL_DISTRIBUTION}")


if __name__ == "__main__":
    main()
