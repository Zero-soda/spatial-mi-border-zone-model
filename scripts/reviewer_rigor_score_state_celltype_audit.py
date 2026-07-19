#!/usr/bin/env python3
"""Reviewer-requested score-state and cell-type marker audits.

This script adds two conservative checks for the public-data-only MI
spatial-state manuscript:

1. Ignore author-provided domain labels during clustering and ask whether
   three-output score space recovers domain-like spatial states.
2. Recalculate domain 3/4 contrasts after residualizing each output against
   canonical cell-type marker scores.
"""

from __future__ import annotations

import csv
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

from project_paths import project_root


ROOT = project_root(__file__)
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import h5py  # noqa: E402
import numpy as np  # noqa: E402

from batch_map_gse214611_mi_spatial_risk import SAMPLES, locate_visium_dir, sample_stage  # noqa: E402
from map_gse214611_d7_1_spatial_risk import sample_slug  # noqa: E402
from spatial_model_utils import bootstrap_ci, read_tsv, safe_float  # noqa: E402


TABLE_DIR = ROOT / "results" / "tables"
SPOT_TABLE = TABLE_DIR / "gse214611_d3_d7_signature_scores_by_spot.tsv"
OUT_STATE_RECOVERY = TABLE_DIR / "gse214611_score_only_state_domain_recovery.tsv"
OUT_CLUSTER_SPOTS = TABLE_DIR / "gse214611_score_only_state_clusters_by_spot.tsv"
OUT_MARKER_SPOTS = TABLE_DIR / "gse214611_celltype_marker_scores_by_spot.tsv"
OUT_ADJUSTED_SAMPLE = TABLE_DIR / "gse214611_celltype_adjusted_domain_contrasts_by_sample.tsv"
OUT_ADJUSTED_SUMMARY = TABLE_DIR / "gse214611_celltype_adjusted_domain_contrasts.tsv"

RANDOM_SEED = 214611

SCORE_COLUMNS = [
    ("mechanical_border", "signature_mechanical_border_score"),
    ("immune_fibrotic_activation", "reviewer_immune_fibrotic_activation_score"),
    ("fibroblast_scar_repair", "signature_fibroblast_scar_score"),
]

SCORE_FEATURE_COLUMNS = [column for _, column in SCORE_COLUMNS]

MARKER_MODULES = {
    "cardiomyocyte_canonical": ["Tnnt2", "Myh6", "Myh7", "Actc1"],
    "resident_fibroblast": ["Dcn", "Lum", "Pdgfra", "Pdgfrl"],
    "myeloid_canonical": ["Ptprc", "Lyz2", "C1qa", "C1qb", "Adgre1"],
    "endothelial_canonical": ["Pecam1", "Kdr", "Vwf"],
    "mural_pericyte_smc": ["Rgs5", "Pdgfrb", "Myh11"],
    "myofibroblast_contractile": ["Acta2", "Tagln", "Myl9", "Tpm2"],
    "activated_scar_fibroblast": ["Postn", "Fap", "Col1a1", "Col1a2"],
}

ADJUSTMENT_MODES = {
    "lineage_marker_covariates": [
        "z_celltype_cardiomyocyte_canonical",
        "z_celltype_resident_fibroblast",
        "z_celltype_myeloid_canonical",
        "z_celltype_endothelial_canonical",
        "z_celltype_mural_pericyte_smc",
    ],
    "expanded_marker_covariates": [
        "z_celltype_cardiomyocyte_canonical",
        "z_celltype_resident_fibroblast",
        "z_celltype_myeloid_canonical",
        "z_celltype_endothelial_canonical",
        "z_celltype_mural_pericyte_smc",
        "z_celltype_myofibroblast_contractile",
        "z_celltype_activated_scar_fibroblast",
    ],
}


def decode(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode()
    return str(value)


def write_lf_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    if not rows:
        raise ValueError(f"No rows to write for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def format_float(value: float, precision: int = 8) -> str:
    if not math.isfinite(value):
        return "nan"
    return f"{value:.{precision}g}"


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    mean = float(np.nanmean(values))
    sd = float(np.nanstd(values))
    if not math.isfinite(sd) or sd == 0:
        sd = 1.0
    return (values - mean) / sd


def add_review_score(row: dict[str, str]) -> None:
    value = (
        safe_float(row["z_CCR2_IL1B_MYLOID"])
        + safe_float(row["z_TGFB_SIGNALING"])
        + safe_float(row["z_FAP_POSTN_PATHO_FIBROBLAST"])
    )
    row["reviewer_immune_fibrotic_activation_score"] = format_float(value)


def rows_by_sample(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["sample"]].append(row)
    return grouped


def comb2(value: int) -> float:
    return value * (value - 1) / 2.0


def adjusted_rand_index(true_labels: list[str], pred_labels: list[int]) -> float:
    contingency: dict[tuple[str, int], int] = Counter(zip(true_labels, pred_labels, strict=True))
    row_counts: Counter[str] = Counter(true_labels)
    col_counts: Counter[int] = Counter(pred_labels)
    n = len(true_labels)
    total_pairs = comb2(n)
    if total_pairs == 0:
        return float("nan")
    sum_comb = sum(comb2(value) for value in contingency.values())
    row_comb = sum(comb2(value) for value in row_counts.values())
    col_comb = sum(comb2(value) for value in col_counts.values())
    expected = row_comb * col_comb / total_pairs
    max_index = 0.5 * (row_comb + col_comb)
    denominator = max_index - expected
    if denominator == 0:
        return 0.0
    return float((sum_comb - expected) / denominator)


def normalized_mutual_information(true_labels: list[str], pred_labels: list[int]) -> float:
    n = len(true_labels)
    if n == 0:
        return float("nan")
    row_counts: Counter[str] = Counter(true_labels)
    col_counts: Counter[int] = Counter(pred_labels)
    joint_counts: Counter[tuple[str, int]] = Counter(zip(true_labels, pred_labels, strict=True))
    mutual_info = 0.0
    for (true_label, pred_label), count in joint_counts.items():
        if count == 0:
            continue
        mutual_info += (count / n) * math.log((count * n) / (row_counts[true_label] * col_counts[pred_label]))
    h_true = -sum((count / n) * math.log(count / n) for count in row_counts.values() if count)
    h_pred = -sum((count / n) * math.log(count / n) for count in col_counts.values() if count)
    denominator = math.sqrt(h_true * h_pred)
    if denominator == 0:
        return 0.0
    return float(mutual_info / denominator)


def kmeans_plus_plus_init(x: np.ndarray, k: int, rng: np.random.Generator) -> np.ndarray:
    n = x.shape[0]
    first = int(rng.integers(0, n))
    centers = [x[first]]
    while len(centers) < k:
        existing = np.vstack(centers)
        distances = ((x[:, None, :] - existing[None, :, :]) ** 2).sum(axis=2).min(axis=1)
        total = float(distances.sum())
        if total <= 0:
            candidate = int(rng.integers(0, n))
        else:
            candidate = int(rng.choice(n, p=distances / total))
        centers.append(x[candidate])
    return np.vstack(centers)


def deterministic_kmeans(x: np.ndarray, k: int, seed: int, n_init: int = 30, max_iter: int = 120) -> tuple[np.ndarray, np.ndarray, float]:
    best_labels: np.ndarray | None = None
    best_centers: np.ndarray | None = None
    best_inertia = float("inf")
    base_rng = np.random.default_rng(seed)
    for _ in range(n_init):
        centers = kmeans_plus_plus_init(x, k, base_rng)
        labels = np.zeros(x.shape[0], dtype=int)
        for _ in range(max_iter):
            distances = ((x[:, None, :] - centers[None, :, :]) ** 2).sum(axis=2)
            new_labels = distances.argmin(axis=1)
            new_centers = centers.copy()
            for cluster_id in range(k):
                member_mask = new_labels == cluster_id
                if member_mask.any():
                    new_centers[cluster_id] = x[member_mask].mean(axis=0)
            if np.array_equal(new_labels, labels):
                centers = new_centers
                labels = new_labels
                break
            centers = new_centers
            labels = new_labels
        inertia = float(((x - centers[labels]) ** 2).sum())
        if inertia < best_inertia:
            best_inertia = inertia
            best_labels = labels.copy()
            best_centers = centers.copy()
    if best_labels is None or best_centers is None:
        raise RuntimeError("k-means failed to initialise")
    return best_labels, best_centers, best_inertia


def label_score_state(center: np.ndarray) -> str:
    names = ["mechanical", "immune_fibrotic", "scar_repair"]
    max_idx = int(np.argmax(center))
    if center[max_idx] < 0.25:
        return "low_or_mixed"
    ordered = sorted(enumerate(center), key=lambda item: item[1], reverse=True)
    if ordered[0][1] - ordered[1][1] < 0.30 and ordered[1][1] > 0.25:
        return f"{names[ordered[0][0]]}_{names[ordered[1][0]]}_mixed"
    return f"{names[max_idx]}_high"


def best_domain_recovery(labels: list[str], clusters: np.ndarray, target_domain: str) -> dict[str, object]:
    total_domain = sum(label == target_domain for label in labels)
    best: dict[str, object] = {
        "cluster": "NA",
        "precision": float("nan"),
        "recall": float("nan"),
        "f1": float("nan"),
        "n_cluster_spots": 0,
        "n_domain_spots": total_domain,
    }
    for cluster_id in sorted(set(clusters.tolist())):
        in_cluster = clusters == cluster_id
        n_cluster = int(in_cluster.sum())
        true_positive = sum(label == target_domain for label, is_member in zip(labels, in_cluster, strict=True) if is_member)
        precision = true_positive / n_cluster if n_cluster else 0.0
        recall = true_positive / total_domain if total_domain else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        if not math.isfinite(float(best["f1"])) or f1 > float(best["f1"]):
            best = {
                "cluster": cluster_id,
                "precision": precision,
                "recall": recall,
                "f1": f1,
                "n_cluster_spots": n_cluster,
                "n_domain_spots": total_domain,
            }
    return best


def write_score_state_recovery(spot_rows: list[dict[str, str]]) -> None:
    recovery_rows = []
    cluster_rows = []
    for sample, rows in sorted(rows_by_sample(spot_rows).items()):
        x = np.array([[safe_float(row[column]) for column in SCORE_FEATURE_COLUMNS] for row in rows], dtype=float)
        x = np.column_stack([zscore(x[:, idx]) for idx in range(x.shape[1])])
        true_labels = [row["annotated"] for row in rows]
        for k in [3, 4]:
            clusters, centers, inertia = deterministic_kmeans(
                x,
                k=k,
                seed=RANDOM_SEED + k * 1000 + sum(ord(char) for char in sample),
            )
            ari = adjusted_rand_index(true_labels, clusters.tolist())
            nmi = normalized_mutual_information(true_labels, clusters.tolist())
            domain3 = best_domain_recovery(true_labels, clusters, "3")
            domain4 = best_domain_recovery(true_labels, clusters, "4")
            state_labels = {idx: label_score_state(centers[idx]) for idx in range(k)}
            recovery_rows.append(
                {
                    "sample": sample,
                    "stage": sample_stage(sample),
                    "n_spots": len(rows),
                    "k_score_states": k,
                    "ari_vs_author_domains": format_float(ari),
                    "nmi_vs_author_domains": format_float(nmi),
                    "inertia": format_float(inertia),
                    "domain3_best_cluster": domain3["cluster"],
                    "domain3_best_cluster_state": state_labels[int(domain3["cluster"])],
                    "domain3_precision": format_float(float(domain3["precision"])),
                    "domain3_recall": format_float(float(domain3["recall"])),
                    "domain3_f1": format_float(float(domain3["f1"])),
                    "domain4_best_cluster": domain4["cluster"],
                    "domain4_best_cluster_state": state_labels[int(domain4["cluster"])],
                    "domain4_precision": format_float(float(domain4["precision"])),
                    "domain4_recall": format_float(float(domain4["recall"])),
                    "domain4_f1": format_float(float(domain4["f1"])),
                }
            )
            if k == 4:
                for row, cluster_id in zip(rows, clusters, strict=True):
                    out = {
                        "sample": sample,
                        "stage": sample_stage(sample),
                        "barcode": row["barcode"],
                        "annotated": row["annotated"],
                        "array_row": row["array_row"],
                        "array_col": row["array_col"],
                        "pxl_row_in_fullres": row["pxl_row_in_fullres"],
                        "pxl_col_in_fullres": row["pxl_col_in_fullres"],
                        "score_only_cluster_k4": int(cluster_id),
                        "score_only_state_label": state_labels[int(cluster_id)],
                    }
                    for score_name, column in SCORE_COLUMNS:
                        out[score_name] = row[column]
                    cluster_rows.append(out)

    write_lf_tsv(
        OUT_STATE_RECOVERY,
        recovery_rows,
        [
            "sample",
            "stage",
            "n_spots",
            "k_score_states",
            "ari_vs_author_domains",
            "nmi_vs_author_domains",
            "inertia",
            "domain3_best_cluster",
            "domain3_best_cluster_state",
            "domain3_precision",
            "domain3_recall",
            "domain3_f1",
            "domain4_best_cluster",
            "domain4_best_cluster_state",
            "domain4_precision",
            "domain4_recall",
            "domain4_f1",
        ],
    )
    write_lf_tsv(
        OUT_CLUSTER_SPOTS,
        cluster_rows,
        [
            "sample",
            "stage",
            "barcode",
            "annotated",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
            "score_only_cluster_k4",
            "score_only_state_label",
            "mechanical_border",
            "immune_fibrotic_activation",
            "fibroblast_scar_repair",
        ],
    )


def find_matrix_h5(visium_dir: Path) -> Path:
    candidates = [
        visium_dir / "filtered_feature_bc_matrix.h5",
        visium_dir / "outs" / "filtered_feature_bc_matrix.h5",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    matches = list(visium_dir.glob("**/filtered_feature_bc_matrix.h5"))
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise FileNotFoundError(f"No filtered_feature_bc_matrix.h5 found under {visium_dir}")
    raise RuntimeError(f"Multiple matrix H5 files found under {visium_dir}: {matches}")


def score_marker_modules(matrix_h5: Path, modules: dict[str, list[str]]) -> dict[str, dict[str, float]]:
    with h5py.File(matrix_h5, "r") as handle:
        matrix = handle["matrix"]
        barcodes = [decode(value) for value in matrix["barcodes"][:]]
        gene_names = [decode(value) for value in matrix["features"]["name"][:]]
        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]

    gene_to_index = {gene.casefold(): idx for idx, gene in enumerate(gene_names)}
    module_indices = {
        module: {gene_to_index[gene.casefold()] for gene in genes if gene.casefold() in gene_to_index}
        for module, genes in modules.items()
    }
    gene_to_modules: dict[int, list[str]] = defaultdict(list)
    for module, gene_indices in module_indices.items():
        for gene_index in gene_indices:
            gene_to_modules[gene_index].append(module)

    scores_by_barcode: dict[str, dict[str, float]] = {}
    for barcode_idx, barcode in enumerate(barcodes):
        start = int(indptr[barcode_idx])
        end = int(indptr[barcode_idx + 1])
        counts = data[start:end].astype(float, copy=False)
        genes = indices[start:end]
        library_size = float(counts.sum())
        sums = {module: 0.0 for module in modules}
        if library_size > 0:
            normalized = np.log1p(counts / library_size * 10000.0)
            for gene_index, value in zip(genes, normalized, strict=True):
                for module in gene_to_modules.get(int(gene_index), []):
                    sums[module] += float(value)
        scores_by_barcode[barcode] = {
            module: sums[module] / len(module_indices[module]) if module_indices[module] else 0.0
            for module in modules
        }
    return scores_by_barcode


def write_marker_scores(spot_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows_by_key = {(row["sample"], row["barcode"]): row for row in spot_rows}
    marker_rows: list[dict[str, str]] = []
    for sample, config in SAMPLES.items():
        visium_dir = locate_visium_dir(config.unpack_dir)
        scores = score_marker_modules(find_matrix_h5(visium_dir), MARKER_MODULES)
        sample_rows = [row for row in spot_rows if row["sample"] == sample]
        for row in sample_rows:
            marker_score = scores.get(row["barcode"])
            if marker_score is None:
                continue
            out = {
                "sample": sample,
                "stage": sample_stage(sample),
                "barcode": row["barcode"],
                "annotated": row["annotated"],
            }
            for score_name, column in SCORE_COLUMNS:
                out[score_name] = row[column]
            for module in MARKER_MODULES:
                out[f"celltype_{module}"] = format_float(marker_score[module])
            marker_rows.append(out)
            rows_by_key[(sample, row["barcode"])].update(out)

    for sample, rows in rows_by_sample(marker_rows).items():
        _ = sample
        for module in MARKER_MODULES:
            values = np.array([safe_float(row[f"celltype_{module}"]) for row in rows], dtype=float)
            scaled = zscore(values)
            for row, value in zip(rows, scaled, strict=True):
                row[f"z_celltype_{module}"] = format_float(float(value))

    fields = [
        "sample",
        "stage",
        "barcode",
        "annotated",
        "mechanical_border",
        "immune_fibrotic_activation",
        "fibroblast_scar_repair",
        *[f"celltype_{module}" for module in MARKER_MODULES],
        *[f"z_celltype_{module}" for module in MARKER_MODULES],
    ]
    write_lf_tsv(OUT_MARKER_SPOTS, marker_rows, fields)
    return marker_rows


def residualize(y: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    valid = np.isfinite(y) & np.all(np.isfinite(covariates), axis=1)
    residuals = np.full_like(y, np.nan, dtype=float)
    if valid.sum() <= covariates.shape[1] + 2:
        return residuals
    design = np.column_stack([np.ones(valid.sum()), covariates[valid]])
    beta = np.linalg.pinv(design) @ y[valid]
    residuals[valid] = y[valid] - design @ beta
    return residuals


def mean_or_nan(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return float("nan")
    return float(sum(finite) / len(finite))


def write_celltype_adjusted_contrasts(marker_rows: list[dict[str, str]]) -> None:
    sample_rows = []
    for mode, covariate_columns in ADJUSTMENT_MODES.items():
        for sample, rows in sorted(rows_by_sample(marker_rows).items()):
            filtered = [row for row in rows if row["annotated"] in {"3", "4"}]
            if not filtered:
                continue
            domains = np.array([row["annotated"] for row in filtered])
            covariates = np.array([[safe_float(row[column]) for column in covariate_columns] for row in filtered], dtype=float)
            for score_name, _ in SCORE_COLUMNS:
                y = np.array([safe_float(row[score_name]) for row in filtered], dtype=float)
                residuals = residualize(y, covariates)
                d3_raw = y[domains == "3"]
                d4_raw = y[domains == "4"]
                d3_adj = residuals[domains == "3"]
                d4_adj = residuals[domains == "4"]
                raw_margin = float(np.nanmean(d3_raw) - np.nanmean(d4_raw))
                adjusted_margin = float(np.nanmean(d3_adj) - np.nanmean(d4_adj))
                sample_rows.append(
                    {
                        "adjustment_mode": mode,
                        "sample": sample,
                        "stage": sample_stage(sample),
                        "score": score_name,
                        "contrast": "domain3_minus_domain4",
                        "n_domain3_spots": int((domains == "3").sum()),
                        "n_domain4_spots": int((domains == "4").sum()),
                        "raw_margin": format_float(raw_margin),
                        "marker_adjusted_residual_margin": format_float(adjusted_margin),
                    }
                )

    summary_rows = []
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in sample_rows:
        grouped[(str(row["adjustment_mode"]), str(row["stage"]), str(row["score"]))].append(row)
    for (mode, stage, score), rows in sorted(grouped.items()):
        raw_margins = [safe_float(row["raw_margin"]) for row in rows]
        adjusted_margins = [safe_float(row["marker_adjusted_residual_margin"]) for row in rows]
        raw_mean, raw_low, raw_high = bootstrap_ci(raw_margins, seed=RANDOM_SEED)
        adjusted_mean, adjusted_low, adjusted_high = bootstrap_ci(adjusted_margins, seed=RANDOM_SEED + 17)
        direction_expected = "domain3_positive" if score == "mechanical_border" else "domain4_positive"
        if score == "mechanical_border":
            preserved = sum(value > 0 for value in adjusted_margins if math.isfinite(value))
        else:
            preserved = sum(value < 0 for value in adjusted_margins if math.isfinite(value))
        summary_rows.append(
            {
                "adjustment_mode": mode,
                "stage": stage,
                "score": score,
                "contrast": "domain3_minus_domain4",
                "n_sample_pairs": len(rows),
                "expected_direction": direction_expected,
                "raw_mean_margin": format_float(raw_mean),
                "raw_ci95_low": format_float(raw_low),
                "raw_ci95_high": format_float(raw_high),
                "marker_adjusted_mean_margin": format_float(adjusted_mean),
                "marker_adjusted_ci95_low": format_float(adjusted_low),
                "marker_adjusted_ci95_high": format_float(adjusted_high),
                "n_samples_preserving_expected_direction_after_adjustment": preserved,
            }
        )

    write_lf_tsv(
        OUT_ADJUSTED_SAMPLE,
        sample_rows,
        [
            "adjustment_mode",
            "sample",
            "stage",
            "score",
            "contrast",
            "n_domain3_spots",
            "n_domain4_spots",
            "raw_margin",
            "marker_adjusted_residual_margin",
        ],
    )
    write_lf_tsv(
        OUT_ADJUSTED_SUMMARY,
        summary_rows,
        [
            "adjustment_mode",
            "stage",
            "score",
            "contrast",
            "n_sample_pairs",
            "expected_direction",
            "raw_mean_margin",
            "raw_ci95_low",
            "raw_ci95_high",
            "marker_adjusted_mean_margin",
            "marker_adjusted_ci95_low",
            "marker_adjusted_ci95_high",
            "n_samples_preserving_expected_direction_after_adjustment",
        ],
    )


def main() -> None:
    spot_rows = read_tsv(SPOT_TABLE)
    for row in spot_rows:
        add_review_score(row)
    write_score_state_recovery(spot_rows)
    marker_rows = write_marker_scores(spot_rows)
    write_celltype_adjusted_contrasts(marker_rows)
    print("Wrote reviewer-requested score-state and cell-type marker audit tables")


if __name__ == "__main__":
    main()
