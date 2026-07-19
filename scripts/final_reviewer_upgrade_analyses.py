#!/usr/bin/env python3
"""Final tutor-requested rigor upgrades for the spatial MI manuscript.

The analyses in this script are designed to address reviewer-style concerns
that remained after the first CVR-oriented revision:

1. Data provenance and spot-level QC.
2. Alternative score construction, including rank-based and
   expression-matched-control scoring.
3. Expression-matched random-signature falsification.
4. A Visium array-neighbour graph boundary analysis.
5. A conservative H&E image-intensity audit.
6. A master spot-level source table for reproducibility.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
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
from PIL import Image  # noqa: E402

from batch_map_gse214611_mi_spatial_risk import SAMPLES, SampleConfig, locate_visium_dir, sample_stage  # noqa: E402
from map_gse214611_d7_1_spatial_risk import sample_slug  # noqa: E402
from score_gse214611_visium_signatures import CORE_SCORE_COLUMNS, read_signatures  # noqa: E402
from spatial_model_utils import bootstrap_ci, read_tsv, safe_float  # noqa: E402


TABLE_DIR = ROOT / "results" / "tables"
SIGNATURES_TSV = Path(__file__).resolve().parents[1] / "config" / "spatial_cardiac_border_zone_signatures.tsv"
SPOT_TABLE = TABLE_DIR / "gse214611_d3_d7_signature_scores_by_spot.tsv"
BOUNDARY_TABLE = TABLE_DIR / "gse214611_d3_d7_domain34_spatial_features_by_spot.tsv"
CLUSTER_TABLE = TABLE_DIR / "gse214611_score_only_state_clusters_by_spot.tsv"
MARKER_TABLE = TABLE_DIR / "gse214611_celltype_marker_scores_by_spot.tsv"
HUMAN_SPOT_TABLE = TABLE_DIR / "gse214611_human_stemi_signature_scores_by_spot.tsv"
RAW_VISIUM = ROOT / "data" / "raw" / "gse214611" / "visium"

OUT_PROVENANCE_QC = TABLE_DIR / "gse214611_data_provenance_qc.tsv"
OUT_MODULE_DETECTION = TABLE_DIR / "gse214611_module_detection_by_sample.tsv"
OUT_ALTERNATIVE_SCORE_SPOTS = TABLE_DIR / "gse214611_alternative_score_sensitivity_by_spot.tsv"
OUT_ALTERNATIVE_SCORE_SUMMARY = TABLE_DIR / "gse214611_alternative_score_sensitivity_summary.tsv"
OUT_EXPR_MATCHED_RANDOM = TABLE_DIR / "gse214611_expression_matched_random_signature_controls.tsv"
OUT_GRAPH_BOUNDARY = TABLE_DIR / "gse214611_visium_graph_boundary_analysis.tsv"
OUT_GRAPH_EDGE_GRADIENTS = TABLE_DIR / "gse214611_visium_graph_boundary_edge_gradients.tsv"
OUT_HISTOLOGY_SPOTS = TABLE_DIR / "gse214611_he_image_intensity_by_spot.tsv"
OUT_HISTOLOGY_SUMMARY = TABLE_DIR / "gse214611_he_image_score_alignment.tsv"
OUT_MASTER = TABLE_DIR / "gse214611_master_spot_level_source_table.tsv"

RANDOM_SEED = 20260707
N_EXPRESSION_MATCHED_RANDOM = 500
CONTROL_GENES_PER_SIGNATURE_GENE = 10

HUMAN_SAMPLE = SampleConfig("Human_STEMI", "GSM6613090", "GSM6613090_V_Human_STEMI.zip")

MOUSE_SAMPLE_TITLES = {
    "D3_1": "Visium, Day3 MI biol rep 1",
    "D3_2": "Visium, Day3 MI biol rep 2",
    "D3_3": "Visium, Day3 MI biol rep 3",
    "D7_1": "Visium, Day7 MI biol rep 1",
    "D7_2": "Visium, Day7 MI biol rep 2",
    "D7_3": "Visium, Day7 MI biol rep 3",
}

SIGNATURE_COMPONENTS = {
    "mechanical_border": ["CM_BZ1_TRANSITION", "CM_BZ2_MECHANICAL_EDGE"],
    "immune_fibrotic_activation": ["CCR2_IL1B_MYLOID", "TGFB_SIGNALING", "FAP_POSTN_PATHO_FIBROBLAST"],
    "immune_fibrotic_overlap_reduced": ["CCR2_IL1B_MYLOID", "TGFB_SIGNALING"],
    "fibroblast_scar_repair": [
        "FAP_POSTN_PATHO_FIBROBLAST",
        "ECM_REMODELING",
        "CTHRC1_REPARATIVE_CF",
        "MYOFIBROBLAST_CONTRACTILE",
    ],
    "fibroblast_scar_overlap_reduced": [
        "ECM_REMODELING",
        "CTHRC1_REPARATIVE_CF",
        "MYOFIBROBLAST_CONTRACTILE",
    ],
}

PRIMARY_CHECKS = [
    ("day3_mi", "mechanical_border", "domain3_minus_domain4", 1.0),
    ("day7_mi", "mechanical_border", "domain3_minus_domain4", 1.0),
    ("day7_mi", "fibroblast_scar_repair", "domain4_minus_domain3", 1.0),
    ("day7_mi", "immune_fibrotic_activation", "domain4_minus_domain3", 1.0),
]


@dataclass(frozen=True)
class MatrixData:
    barcodes: list[str]
    gene_names: list[str]
    data: np.ndarray
    indices: np.ndarray
    indptr: np.ndarray


def decode(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode()
    return str(value)


def format_float(value: float, precision: int = 8) -> str:
    if not math.isfinite(value):
        return "nan"
    return f"{value:.{precision}g}"


def write_lf_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
        raise FileNotFoundError(f"No filtered_feature_bc_matrix.h5 under {visium_dir}")
    raise RuntimeError(f"Multiple filtered_feature_bc_matrix.h5 files under {visium_dir}: {matches}")


def read_10x_h5(matrix_h5: Path) -> MatrixData:
    with h5py.File(matrix_h5, "r") as handle:
        matrix = handle["matrix"]
        return MatrixData(
            barcodes=[decode(value) for value in matrix["barcodes"][:]],
            gene_names=[decode(value) for value in matrix["features"]["name"][:]],
            data=matrix["data"][:],
            indices=matrix["indices"][:],
            indptr=matrix["indptr"][:],
        )


def load_dense_lognorm(matrix_data: MatrixData) -> np.ndarray:
    dense = np.zeros((len(matrix_data.barcodes), len(matrix_data.gene_names)), dtype=np.float32)
    for barcode_idx in range(len(matrix_data.barcodes)):
        start = int(matrix_data.indptr[barcode_idx])
        end = int(matrix_data.indptr[barcode_idx + 1])
        counts = matrix_data.data[start:end].astype(np.float32, copy=False)
        genes = matrix_data.indices[start:end]
        library_size = float(np.sum(counts))
        if library_size <= 0:
            continue
        dense[barcode_idx, genes] = np.log1p(counts / library_size * 10000.0)
    return dense


def zscore(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    mean = float(np.nanmean(arr))
    sd = float(np.nanstd(arr))
    if not math.isfinite(sd) or sd == 0:
        sd = 1.0
    return (arr - mean) / sd


def matrix_qc(matrix_data: MatrixData, barcode_filter: set[str] | None = None) -> dict[str, float]:
    mito = np.array([gene.lower().startswith("mt-") for gene in matrix_data.gene_names], dtype=bool)
    umi_counts = []
    gene_counts = []
    mito_fractions = []
    for barcode_idx, barcode in enumerate(matrix_data.barcodes):
        if barcode_filter is not None and barcode not in barcode_filter:
            continue
        start = int(matrix_data.indptr[barcode_idx])
        end = int(matrix_data.indptr[barcode_idx + 1])
        counts = matrix_data.data[start:end].astype(float, copy=False)
        genes = matrix_data.indices[start:end]
        total = float(counts.sum())
        umi_counts.append(total)
        gene_counts.append(float(len(genes)))
        mito_count = float(counts[mito[genes]].sum()) if len(genes) else 0.0
        mito_fractions.append(mito_count / total if total > 0 else float("nan"))
    return {
        "n_barcodes": float(len(umi_counts)),
        "median_umis": float(np.nanmedian(umi_counts)) if umi_counts else float("nan"),
        "median_genes": float(np.nanmedian(gene_counts)) if gene_counts else float("nan"),
        "median_mito_fraction": float(np.nanmedian(mito_fractions)) if mito_fractions else float("nan"),
    }


def gene_index_map(matrix_data: MatrixData) -> dict[str, int]:
    return {gene.casefold(): idx for idx, gene in enumerate(matrix_data.gene_names)}


def signature_indices(signatures: dict[str, list[str]], matrix_data: MatrixData) -> dict[str, list[int]]:
    lookup = gene_index_map(matrix_data)
    return {
        signature_id: [lookup[gene.casefold()] for gene in genes if gene.casefold() in lookup]
        for signature_id, genes in signatures.items()
    }


def detected_gene_sets(matrix_data: MatrixData) -> set[int]:
    detected: set[int] = set()
    for barcode_idx in range(len(matrix_data.barcodes)):
        start = int(matrix_data.indptr[barcode_idx])
        end = int(matrix_data.indptr[barcode_idx + 1])
        detected.update(int(gene) for gene in matrix_data.indices[start:end])
    return detected


def domain_count_text(rows: list[dict[str, str]]) -> str:
    counts = Counter(row["annotated"] for row in rows)
    return ";".join(f"{domain}:{counts[domain]}" for domain in sorted(counts, key=str))


def build_provenance_qc(signatures: dict[str, list[str]], mouse_rows: list[dict[str, str]]) -> None:
    rows_by_sample = defaultdict(list)
    for row in mouse_rows:
        rows_by_sample[row["sample"]].append(row)

    provenance_rows = []
    module_rows = []

    all_configs: list[tuple[str, SampleConfig, str, str, str]] = [
        (sample, config, "Mus musculus", sample_stage(sample), MOUSE_SAMPLE_TITLES[sample])
        for sample, config in SAMPLES.items()
    ]
    all_configs.append(("Human_STEMI", HUMAN_SAMPLE, "Homo sapiens", "human_stemi", "Visium, STEMI"))

    human_rows = read_tsv(HUMAN_SPOT_TABLE) if HUMAN_SPOT_TABLE.exists() else []
    human_barcodes = {row["barcode"] for row in human_rows}

    for sample, config, species, stage, title in all_configs:
        visium_dir = locate_visium_dir(config.unpack_dir)
        matrix_h5 = find_matrix_h5(visium_dir)
        matrix_data = read_10x_h5(matrix_h5)
        if sample == "Human_STEMI":
            analysis_rows = human_rows
            barcode_filter = human_barcodes
            domain_counts = "not_applicable"
        else:
            analysis_rows = rows_by_sample[sample]
            barcode_filter = {row["barcode"] for row in analysis_rows}
            domain_counts = domain_count_text(analysis_rows)
        qc = matrix_qc(matrix_data, barcode_filter=barcode_filter)
        archive_sha = file_sha256(config.archive_path) if config.archive_path.exists() else "NA"
        replacement_note = (
            "GEO states processed files were replaced on 2023-05-18 because prior files were incomplete"
            if config.gsm in {"GSM6613087", "GSM6613090"}
            else "not_flagged_by_geo_replacement_note"
        )
        provenance_rows.append(
            {
                "analysis_sample": sample,
                "geo_accession": config.gsm,
                "geo_sample_title": title,
                "species": species,
                "manuscript_group": stage,
                "public_file": config.archive_name,
                "public_file_name_note": "file name contains d1 but GEO title identifies Day3 MI" if sample.startswith("D3") else "",
                "n_spots_in_analysis": len(analysis_rows),
                "n_barcodes_in_filtered_matrix": int(qc["n_barcodes"]),
                "median_umis_in_analysis_spots": format_float(qc["median_umis"]),
                "median_genes_in_analysis_spots": format_float(qc["median_genes"]),
                "median_mito_fraction_in_analysis_spots": format_float(qc["median_mito_fraction"]),
                "domain_spot_counts": domain_counts,
                "local_matrix_h5": str(matrix_h5.relative_to(ROOT)),
                "matrix_h5_sha256": file_sha256(matrix_h5),
                "local_archive": str(config.archive_path.relative_to(ROOT)) if config.archive_path.exists() else "NA",
                "archive_sha256": archive_sha,
                "geo_accessed_date": "2026-07-07",
                "geo_last_update_date_observed": "2024-08-26",
                "geo_processed_file_replacement_note": replacement_note,
            }
        )

        sig_idx = signature_indices(signatures, matrix_data)
        detected = detected_gene_sets(matrix_data)
        for signature_id, genes in signatures.items():
            present_indices = sig_idx[signature_id]
            detected_indices = [idx for idx in present_indices if idx in detected]
            module_rows.append(
                {
                    "analysis_sample": sample,
                    "geo_accession": config.gsm,
                    "species": species,
                    "manuscript_group": stage,
                    "signature_id": signature_id,
                    "n_signature_genes": len(genes),
                    "n_present_in_matrix": len(present_indices),
                    "n_detected_nonzero_in_filtered_matrix": len(detected_indices),
                    "present_genes": ";".join(matrix_data.gene_names[idx] for idx in present_indices),
                    "detected_genes": ";".join(matrix_data.gene_names[idx] for idx in detected_indices),
                }
            )

    write_lf_tsv(
        OUT_PROVENANCE_QC,
        provenance_rows,
        [
            "analysis_sample",
            "geo_accession",
            "geo_sample_title",
            "species",
            "manuscript_group",
            "public_file",
            "public_file_name_note",
            "n_spots_in_analysis",
            "n_barcodes_in_filtered_matrix",
            "median_umis_in_analysis_spots",
            "median_genes_in_analysis_spots",
            "median_mito_fraction_in_analysis_spots",
            "domain_spot_counts",
            "local_matrix_h5",
            "matrix_h5_sha256",
            "local_archive",
            "archive_sha256",
            "geo_accessed_date",
            "geo_last_update_date_observed",
            "geo_processed_file_replacement_note",
        ],
    )
    write_lf_tsv(
        OUT_MODULE_DETECTION,
        module_rows,
        [
            "analysis_sample",
            "geo_accession",
            "species",
            "manuscript_group",
            "signature_id",
            "n_signature_genes",
            "n_present_in_matrix",
            "n_detected_nonzero_in_filtered_matrix",
            "present_genes",
            "detected_genes",
        ],
    )


def gene_stats(dense: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean_expression = np.asarray(dense.mean(axis=0), dtype=float)
    detection = np.asarray((dense > 0).mean(axis=0), dtype=float)
    return mean_expression, detection


def matched_control_pool(
    signature_gene_indices: list[int],
    all_signature_indices: set[int],
    mean_expression: np.ndarray,
    detection: np.ndarray,
    rng: np.random.Generator,
    n_per_gene: int = CONTROL_GENES_PER_SIGNATURE_GENE,
) -> list[int]:
    if not signature_gene_indices:
        return []
    gene_count = len(mean_expression)
    mean_bins = np.digitize(mean_expression, np.quantile(mean_expression, np.linspace(0, 1, 11)[1:-1]))
    detect_bins = np.digitize(detection, np.quantile(detection, np.linspace(0, 1, 11)[1:-1]))
    excluded = set(all_signature_indices)
    selected: list[int] = []
    all_indices = np.arange(gene_count)
    for gene_idx in signature_gene_indices:
        target_mean_bin = mean_bins[gene_idx]
        target_detect_bin = detect_bins[gene_idx]
        mask = (mean_bins == target_mean_bin) & (detect_bins == target_detect_bin)
        candidates = [int(idx) for idx in all_indices[mask] if int(idx) not in excluded]
        if len(candidates) < n_per_gene:
            mask = mean_bins == target_mean_bin
            candidates = [int(idx) for idx in all_indices[mask] if int(idx) not in excluded]
        if len(candidates) < n_per_gene:
            candidates = [int(idx) for idx in all_indices if int(idx) not in excluded]
        replace = len(candidates) < n_per_gene
        selected.extend(rng.choice(candidates, size=n_per_gene, replace=replace).astype(int).tolist())
    return selected


def module_mean(dense: np.ndarray, indices: list[int]) -> np.ndarray:
    if not indices:
        return np.zeros(dense.shape[0], dtype=float)
    return np.asarray(dense[:, indices].mean(axis=1), dtype=float)


def module_rank_score(dense: np.ndarray, indices: list[int]) -> np.ndarray:
    if not indices:
        return np.zeros(dense.shape[0], dtype=float)
    out = np.zeros(dense.shape[0], dtype=float)
    index_set = set(indices)
    for row_idx in range(dense.shape[0]):
        values = dense[row_idx]
        nonzero = np.flatnonzero(values > 0)
        if len(nonzero) == 0:
            continue
        order = nonzero[np.argsort(values[nonzero])]
        percentiles: dict[int, float] = {}
        denom = max(len(order) - 1, 1)
        for rank, gene_idx in enumerate(order):
            if int(gene_idx) in index_set:
                percentiles[int(gene_idx)] = rank / denom
        out[row_idx] = sum(percentiles.get(gene_idx, 0.0) for gene_idx in indices) / len(indices)
    return out


def composite_from_modules(module_scores: dict[str, np.ndarray], output_name: str) -> np.ndarray:
    components = SIGNATURE_COMPONENTS[output_name]
    if not components:
        return np.zeros(next(iter(module_scores.values())).shape[0], dtype=float)
    scaled = [zscore(module_scores[component]) for component in components]
    return np.sum(np.vstack(scaled), axis=0)


def summarize_primary_contrasts(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_key: dict[tuple[str, str, str, str], dict[tuple[str, str], float]] = defaultdict(dict)
    for row in rows:
        if str(row["annotated"]) not in {"3", "4"}:
            continue
        key = (str(row["scoring_method"]), str(row["stage"]), str(row["score"]))
        by_key[key][(str(row["sample"]), str(row["annotated"]))] = safe_float(row["mean_score"])

    out = []
    for method in sorted({str(row["scoring_method"]) for row in rows}):
        for stage, score, contrast, expected_sign in PRIMARY_CHECKS:
            margins = []
            stage_samples = sorted({sample for sample, sample_stage_value in {(str(r["sample"]), str(r["stage"])) for r in rows} if sample_stage_value == stage})
            for sample in stage_samples:
                d3 = by_key.get((method, stage, score), {}).get((sample, "3"))
                d4 = by_key.get((method, stage, score), {}).get((sample, "4"))
                if d3 is None or d4 is None:
                    continue
                margin = d3 - d4 if contrast == "domain3_minus_domain4" else d4 - d3
                margins.append(margin)
            mean_margin, low, high = bootstrap_ci(margins, seed=RANDOM_SEED + len(method) + len(score))
            out.append(
                {
                    "scoring_method": method,
                    "stage": stage,
                    "score": score,
                    "contrast": contrast,
                    "n_sample_pairs": len(margins),
                    "mean_margin": format_float(mean_margin),
                    "ci95_low": format_float(low),
                    "ci95_high": format_float(high),
                    "expected_direction_preserved_n": sum(value * expected_sign > 0 for value in margins),
                    "sample_margins": ";".join(format_float(value, precision=6) for value in margins),
                }
            )
    return out


def build_alternative_scoring(signatures: dict[str, list[str]], mouse_rows: list[dict[str, str]]) -> None:
    rows_by_sample = defaultdict(list)
    for row in mouse_rows:
        rows_by_sample[row["sample"]].append(row)

    spot_out = []
    summary_input = []
    rng = np.random.default_rng(RANDOM_SEED)

    for sample, config in SAMPLES.items():
        sample_rows = rows_by_sample[sample]
        visium_dir = locate_visium_dir(config.unpack_dir)
        matrix_data = read_10x_h5(find_matrix_h5(visium_dir))
        dense = load_dense_lognorm(matrix_data)
        barcode_to_row = {barcode: idx for idx, barcode in enumerate(matrix_data.barcodes)}
        sig_idx = signature_indices(signatures, matrix_data)
        all_sig = {idx for indices in sig_idx.values() for idx in indices}
        mean_expression, detection = gene_stats(dense)

        module_mean_scores = {signature_id: module_mean(dense, indices) for signature_id, indices in sig_idx.items()}
        module_rank_scores = {signature_id: module_rank_score(dense, indices) for signature_id, indices in sig_idx.items()}
        module_control_scores = {}
        for signature_id, indices in sig_idx.items():
            controls = matched_control_pool(indices, all_sig, mean_expression, detection, rng)
            module_control_scores[signature_id] = module_mean(dense, indices) - module_mean(dense, controls)

        method_outputs = {
            "mean_expression_zsum": {
                "mechanical_border": composite_from_modules(module_mean_scores, "mechanical_border"),
                "immune_fibrotic_activation": composite_from_modules(module_mean_scores, "immune_fibrotic_activation"),
                "fibroblast_scar_repair": composite_from_modules(module_mean_scores, "fibroblast_scar_repair"),
            },
            "rank_auc_zsum": {
                "mechanical_border": composite_from_modules(module_rank_scores, "mechanical_border"),
                "immune_fibrotic_activation": composite_from_modules(module_rank_scores, "immune_fibrotic_activation"),
                "fibroblast_scar_repair": composite_from_modules(module_rank_scores, "fibroblast_scar_repair"),
            },
            "expression_matched_control_zsum": {
                "mechanical_border": composite_from_modules(module_control_scores, "mechanical_border"),
                "immune_fibrotic_activation": composite_from_modules(module_control_scores, "immune_fibrotic_activation"),
                "fibroblast_scar_repair": composite_from_modules(module_control_scores, "fibroblast_scar_repair"),
            },
            "overlap_reduced_mean_zsum": {
                "mechanical_border": composite_from_modules(module_mean_scores, "mechanical_border"),
                "immune_fibrotic_activation": composite_from_modules(module_mean_scores, "immune_fibrotic_overlap_reduced"),
                "fibroblast_scar_repair": composite_from_modules(module_mean_scores, "fibroblast_scar_overlap_reduced"),
            },
        }

        for row in sample_rows:
            matrix_idx = barcode_to_row.get(row["barcode"])
            if matrix_idx is None:
                continue
            for method, outputs in method_outputs.items():
                for score, values in outputs.items():
                    spot_out.append(
                        {
                            "sample": sample,
                            "stage": sample_stage(sample),
                            "barcode": row["barcode"],
                            "annotated": row["annotated"],
                            "scoring_method": method,
                            "score": score,
                            "score_value": format_float(float(values[matrix_idx])),
                        }
                    )

        grouped: dict[tuple[str, str, str], list[float]] = defaultdict(list)
        for row in sample_rows:
            matrix_idx = barcode_to_row.get(row["barcode"])
            if matrix_idx is None:
                continue
            for method, outputs in method_outputs.items():
                for score, values in outputs.items():
                    grouped[(method, row["annotated"], score)].append(float(values[matrix_idx]))
        for (method, domain, score), values in grouped.items():
            summary_input.append(
                {
                    "sample": sample,
                    "stage": sample_stage(sample),
                    "annotated": domain,
                    "scoring_method": method,
                    "score": score,
                    "mean_score": float(np.mean(values)),
                }
            )

    summary_rows = summarize_primary_contrasts(summary_input)
    write_lf_tsv(
        OUT_ALTERNATIVE_SCORE_SPOTS,
        spot_out,
        ["sample", "stage", "barcode", "annotated", "scoring_method", "score", "score_value"],
    )
    write_lf_tsv(
        OUT_ALTERNATIVE_SCORE_SUMMARY,
        summary_rows,
        [
            "scoring_method",
            "stage",
            "score",
            "contrast",
            "n_sample_pairs",
            "mean_margin",
            "ci95_low",
            "ci95_high",
            "expected_direction_preserved_n",
            "sample_margins",
        ],
    )


def expression_matched_random_controls(signatures: dict[str, list[str]], mouse_rows: list[dict[str, str]]) -> None:
    rows_by_sample = defaultdict(list)
    for row in mouse_rows:
        rows_by_sample[row["sample"]].append(row)

    out = []
    for stage, score, contrast, _ in PRIMARY_CHECKS:
        stage_samples = [sample for sample in SAMPLES if sample_stage(sample) == stage]
        observed_margins = []
        null_margins = []
        for sample in stage_samples:
            config = SAMPLES[sample]
            sample_rows = rows_by_sample[sample]
            visium_dir = locate_visium_dir(config.unpack_dir)
            matrix_data = read_10x_h5(find_matrix_h5(visium_dir))
            dense = load_dense_lognorm(matrix_data)
            barcode_to_row = {barcode: idx for idx, barcode in enumerate(matrix_data.barcodes)}
            sig_idx = signature_indices(signatures, matrix_data)
            all_sig = {idx for indices in sig_idx.values() for idx in indices}
            mean_expression, detection = gene_stats(dense)
            rng = np.random.default_rng(RANDOM_SEED + sum(ord(char) for char in f"{stage}:{score}:{sample}"))

            observed_modules = {signature_id: module_mean(dense, indices) for signature_id, indices in sig_idx.items()}
            observed_score = composite_from_modules(observed_modules, score)
            d3_indices = [barcode_to_row[row["barcode"]] for row in sample_rows if row["annotated"] == "3" and row["barcode"] in barcode_to_row]
            d4_indices = [barcode_to_row[row["barcode"]] for row in sample_rows if row["annotated"] == "4" and row["barcode"] in barcode_to_row]
            if not d3_indices or not d4_indices:
                continue
            observed = float(np.mean(observed_score[d3_indices]) - np.mean(observed_score[d4_indices]))
            if contrast == "domain4_minus_domain3":
                observed = -observed
            observed_margins.append(observed)

            sample_null = []
            for _ in range(N_EXPRESSION_MATCHED_RANDOM):
                random_modules = {}
                for component in SIGNATURE_COMPONENTS[score]:
                    controls = matched_control_pool(sig_idx[component], all_sig, mean_expression, detection, rng, n_per_gene=1)
                    random_modules[component] = module_mean(dense, controls)
                random_score = composite_from_modules(random_modules, score)
                margin = float(np.mean(random_score[d3_indices]) - np.mean(random_score[d4_indices]))
                if contrast == "domain4_minus_domain3":
                    margin = -margin
                sample_null.append(margin)
            null_margins.append(sample_null)

        if not observed_margins:
            continue
        observed_stage_margin = float(np.mean(observed_margins))
        null_stage = np.mean(np.asarray(null_margins, dtype=float), axis=0) if null_margins else np.array([])
        p_upper = float((np.sum(null_stage >= observed_stage_margin) + 1) / (len(null_stage) + 1)) if len(null_stage) else float("nan")
        out.append(
            {
                "stage": stage,
                "score": score,
                "contrast": contrast,
                "n_expression_matched_random_signatures": len(null_stage),
                "observed_mean_margin": format_float(observed_stage_margin),
                "observed_sample_margins": ";".join(format_float(value, precision=6) for value in observed_margins),
                "null_p50": format_float(float(np.percentile(null_stage, 50))) if len(null_stage) else "nan",
                "null_p95": format_float(float(np.percentile(null_stage, 95))) if len(null_stage) else "nan",
                "null_p99": format_float(float(np.percentile(null_stage, 99))) if len(null_stage) else "nan",
                "empirical_p_upper": format_float(p_upper),
                "exceeds_null_p95": bool(len(null_stage) and observed_stage_margin > float(np.percentile(null_stage, 95))),
                "exceeds_null_p99": bool(len(null_stage) and observed_stage_margin > float(np.percentile(null_stage, 99))),
            }
        )

    write_lf_tsv(
        OUT_EXPR_MATCHED_RANDOM,
        out,
        [
            "stage",
            "score",
            "contrast",
            "n_expression_matched_random_signatures",
            "observed_mean_margin",
            "observed_sample_margins",
            "null_p50",
            "null_p95",
            "null_p99",
            "empirical_p_upper",
            "exceeds_null_p95",
            "exceeds_null_p99",
        ],
    )


def six_nearest_edges(coords: np.ndarray) -> set[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    for idx in range(coords.shape[0]):
        delta = coords - coords[idx]
        distances = np.sqrt(np.sum(delta * delta, axis=1))
        order = np.argsort(distances)
        for neighbour in order[1:7]:
            edge = (idx, int(neighbour)) if idx < int(neighbour) else (int(neighbour), idx)
            edges.add(edge)
    return edges


def build_graph_boundary(mouse_rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, object]]:
    rows_by_sample = defaultdict(list)
    for row in mouse_rows:
        rows_by_sample[row["sample"]].append(row)

    graph_status: dict[tuple[str, str], dict[str, object]] = {}
    sample_out = []
    edge_out = []
    for sample, rows in sorted(rows_by_sample.items()):
        coords = np.array([[safe_float(row["array_row"]), safe_float(row["array_col"])] for row in rows], dtype=float)
        edges = six_nearest_edges(coords)
        domains = [row["annotated"] for row in rows]
        edge_domain_counts = Counter()
        d3_touch = set()
        d4_touch = set()
        same_domain = 0
        cross_domain = 0
        d3d4_edges = []
        score_cols = {
            "mechanical_border": "signature_mechanical_border_score",
            "immune_fibrotic_activation": "signature_fibrotic_risk",
            "fibroblast_scar_repair": "signature_fibroblast_scar_score",
        }
        gradient_values: dict[str, list[float]] = defaultdict(list)
        for i, j in edges:
            left, right = domains[i], domains[j]
            pair = "-".join(sorted([left, right], key=str))
            edge_domain_counts[pair] += 1
            if left == right:
                same_domain += 1
            else:
                cross_domain += 1
            if {left, right} == {"3", "4"}:
                d3d4_edges.append((i, j))
                if left == "3":
                    d3_idx, d4_idx = i, j
                else:
                    d3_idx, d4_idx = j, i
                d3_touch.add(d3_idx)
                d4_touch.add(d4_idx)
                for score, column in score_cols.items():
                    gradient_values[score].append(safe_float(rows[d4_idx][column]) - safe_float(rows[d3_idx][column]))

        n_d3 = sum(domain == "3" for domain in domains)
        n_d4 = sum(domain == "4" for domain in domains)
        interface_text = ";".join(f"{pair}:{edge_domain_counts[pair]}" for pair in sorted(edge_domain_counts))
        sample_out.append(
            {
                "sample": sample,
                "stage": sample_stage(sample),
                "n_spots": len(rows),
                "n_graph_edges": len(edges),
                "same_domain_edge_fraction": format_float(same_domain / len(edges) if edges else float("nan")),
                "cross_domain_edge_fraction": format_float(cross_domain / len(edges) if edges else float("nan")),
                "domain3_domain4_edge_count": len(d3d4_edges),
                "domain3_domain4_fraction_of_all_edges": format_float(len(d3d4_edges) / len(edges) if edges else float("nan")),
                "domain3_domain4_fraction_of_cross_domain_edges": format_float(len(d3d4_edges) / cross_domain if cross_domain else float("nan")),
                "domain3_graph_boundary_fraction": format_float(len(d3_touch) / n_d3 if n_d3 else float("nan")),
                "domain4_graph_boundary_fraction": format_float(len(d4_touch) / n_d4 if n_d4 else float("nan")),
                "all_interface_edge_counts": interface_text,
                "graph_definition": "symmetrized_six_nearest_neighbours_in_visium_array_coordinates",
            }
        )
        for row_idx, row in enumerate(rows):
            graph_status[(sample, row["barcode"])] = {
                "is_domain34_graph_boundary": row_idx in d3_touch or row_idx in d4_touch,
                "graph_degree_domain34": sum(1 for i, j in d3d4_edges if i == row_idx or j == row_idx),
            }
        for score, values in gradient_values.items():
            values_arr = np.asarray(values, dtype=float)
            edge_out.append(
                {
                    "sample": sample,
                    "stage": sample_stage(sample),
                    "score": score,
                    "contrast": "domain4_minus_domain3_across_graph_edges",
                    "n_domain34_graph_edges": len(values_arr),
                    "mean_edge_delta": format_float(float(np.nanmean(values_arr)) if len(values_arr) else float("nan")),
                    "median_edge_delta": format_float(float(np.nanmedian(values_arr)) if len(values_arr) else float("nan")),
                }
            )

    write_lf_tsv(
        OUT_GRAPH_BOUNDARY,
        sample_out,
        [
            "sample",
            "stage",
            "n_spots",
            "n_graph_edges",
            "same_domain_edge_fraction",
            "cross_domain_edge_fraction",
            "domain3_domain4_edge_count",
            "domain3_domain4_fraction_of_all_edges",
            "domain3_domain4_fraction_of_cross_domain_edges",
            "domain3_graph_boundary_fraction",
            "domain4_graph_boundary_fraction",
            "all_interface_edge_counts",
            "graph_definition",
        ],
    )
    write_lf_tsv(
        OUT_GRAPH_EDGE_GRADIENTS,
        edge_out,
        ["sample", "stage", "score", "contrast", "n_domain34_graph_edges", "mean_edge_delta", "median_edge_delta"],
    )
    return graph_status


def find_hires_image(visium_dir: Path) -> Path:
    candidates = [
        visium_dir / "spatial" / "tissue_hires_image.png",
        visium_dir / "tissue_hires_image.png",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    matches = list(visium_dir.glob("**/tissue_hires_image.png"))
    if len(matches) == 1:
        return matches[0]
    raise FileNotFoundError(f"No unique tissue_hires_image.png under {visium_dir}")


def find_scalefactors(visium_dir: Path) -> Path:
    candidates = [
        visium_dir / "spatial" / "scalefactors_json.json",
        visium_dir / "scalefactors_json.json",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    matches = list(visium_dir.glob("**/scalefactors_json.json"))
    if len(matches) == 1:
        return matches[0]
    raise FileNotFoundError(f"No unique scalefactors_json.json under {visium_dir}")


def build_histology_audit(mouse_rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, object]]:
    rows_by_sample = defaultdict(list)
    for row in mouse_rows:
        rows_by_sample[row["sample"]].append(row)

    histology_status: dict[tuple[str, str], dict[str, object]] = {}
    spot_out = []
    summary_out = []
    for sample, rows in sorted(rows_by_sample.items()):
        visium_dir = locate_visium_dir(SAMPLES[sample].unpack_dir)
        with find_scalefactors(visium_dir).open() as handle:
            scales = json.load(handle)
        hires_scale = float(scales.get("tissue_hires_scalef", 1.0))
        image = Image.open(find_hires_image(visium_dir)).convert("RGB")
        arr = np.asarray(image, dtype=float)
        patch_radius = max(int(0.45 * float(scales.get("spot_diameter_fullres", 90.0)) * hires_scale), 6)
        for row in rows:
            y = int(round(safe_float(row["pxl_row_in_fullres"]) * hires_scale))
            x = int(round(safe_float(row["pxl_col_in_fullres"]) * hires_scale))
            y0, y1 = max(y - patch_radius, 0), min(y + patch_radius + 1, arr.shape[0])
            x0, x1 = max(x - patch_radius, 0), min(x + patch_radius + 1, arr.shape[1])
            patch = arr[y0:y1, x0:x1]
            if patch.size == 0:
                mean_r = mean_g = mean_b = brightness = saturation = float("nan")
            else:
                mean_r, mean_g, mean_b = patch.reshape(-1, 3).mean(axis=0).tolist()
                rgb = patch.reshape(-1, 3) / 255.0
                brightness = float(rgb.mean())
                saturation = float((rgb.max(axis=1) - rgb.min(axis=1)).mean())
            blue_red_od = math.log((mean_r + 1.0) / (mean_b + 1.0)) if math.isfinite(mean_r) and math.isfinite(mean_b) else float("nan")
            out = {
                "sample": sample,
                "stage": sample_stage(sample),
                "barcode": row["barcode"],
                "annotated": row["annotated"],
                "he_mean_r": format_float(mean_r),
                "he_mean_g": format_float(mean_g),
                "he_mean_b": format_float(mean_b),
                "he_brightness": format_float(brightness),
                "he_saturation": format_float(saturation),
                "he_blue_red_optical_density_proxy": format_float(blue_red_od),
                "histology_patch_radius_hires_px": patch_radius,
            }
            spot_out.append(out)
            histology_status[(sample, row["barcode"])] = out

        sample_spots = [row for row in spot_out if row["sample"] == sample]
        scores = {
            "mechanical_border": [safe_float(row["signature_mechanical_border_score"]) for row in rows],
            "fibroblast_scar_repair": [safe_float(row["signature_fibroblast_scar_score"]) for row in rows],
            "immune_fibrotic_activation": [safe_float(row["signature_fibrotic_risk"]) for row in rows],
        }
        features = {
            "he_brightness": [safe_float(row["he_brightness"]) for row in sample_spots],
            "he_saturation": [safe_float(row["he_saturation"]) for row in sample_spots],
            "he_blue_red_optical_density_proxy": [
                safe_float(row["he_blue_red_optical_density_proxy"]) for row in sample_spots
            ],
        }
        for score_name, score_values in scores.items():
            for feature, feature_values in features.items():
                x = np.asarray(feature_values, dtype=float)
                y = np.asarray(score_values, dtype=float)
                mask = np.isfinite(x) & np.isfinite(y)
                corr = float(np.corrcoef(x[mask], y[mask])[0, 1]) if mask.sum() > 2 else float("nan")
                summary_out.append(
                    {
                        "sample": sample,
                        "stage": sample_stage(sample),
                        "score": score_name,
                        "histology_feature": feature,
                        "n_spots": int(mask.sum()),
                        "pearson_r": format_float(corr),
                        "interpretation": "image-intensity audit only; not pathologist-grade histological segmentation",
                    }
                )

    write_lf_tsv(
        OUT_HISTOLOGY_SPOTS,
        spot_out,
        [
            "sample",
            "stage",
            "barcode",
            "annotated",
            "he_mean_r",
            "he_mean_g",
            "he_mean_b",
            "he_brightness",
            "he_saturation",
            "he_blue_red_optical_density_proxy",
            "histology_patch_radius_hires_px",
        ],
    )
    write_lf_tsv(
        OUT_HISTOLOGY_SUMMARY,
        summary_out,
        ["sample", "stage", "score", "histology_feature", "n_spots", "pearson_r", "interpretation"],
    )
    return histology_status


def read_by_sample_barcode(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    if not path.exists():
        return {}
    return {(row["sample"], row["barcode"]): row for row in read_tsv(path)}


def build_master_table(mouse_rows: list[dict[str, str]], graph_status: dict[tuple[str, str], dict[str, object]]) -> None:
    boundary = read_by_sample_barcode(BOUNDARY_TABLE)
    clusters = read_by_sample_barcode(CLUSTER_TABLE)
    markers = read_by_sample_barcode(MARKER_TABLE)
    histology = read_by_sample_barcode(OUT_HISTOLOGY_SPOTS)

    rows = []
    for row in mouse_rows:
        key = (row["sample"], row["barcode"])
        out = {
            "sample": row["sample"],
            "stage": row["stage"],
            "species": "Mus musculus",
            "barcode": row["barcode"],
            "tissue_flag": "1",
            "author_domain": row["annotated"],
            "array_row": row["array_row"],
            "array_col": row["array_col"],
            "pxl_row_in_fullres": row["pxl_row_in_fullres"],
            "pxl_col_in_fullres": row["pxl_col_in_fullres"],
            "mechanical_border_score": row["signature_mechanical_border_score"],
            "immune_fibrotic_activation_score": row["signature_fibrotic_risk"],
            "fibroblast_scar_repair_score": row["signature_fibroblast_scar_score"],
            "domain34_signed_distance_fullres_px": boundary.get(key, {}).get("domain34_signed_distance_fullres_px", ""),
            "is_domain34_pixel_boundary": boundary.get(key, {}).get("is_domain34_contact_boundary", ""),
            "nearest_opposite_domain34_distance_fullres_px": boundary.get(key, {}).get(
                "nearest_opposite_domain34_distance_fullres_px", ""
            ),
            "score_only_state_label": clusters.get(key, {}).get("score_only_state_label", ""),
            "is_domain34_graph_boundary": graph_status.get(key, {}).get("is_domain34_graph_boundary", ""),
            "graph_degree_domain34": graph_status.get(key, {}).get("graph_degree_domain34", ""),
            "he_brightness": histology.get(key, {}).get("he_brightness", ""),
            "he_saturation": histology.get(key, {}).get("he_saturation", ""),
            "he_blue_red_optical_density_proxy": histology.get(key, {}).get("he_blue_red_optical_density_proxy", ""),
        }
        for module in CORE_SCORE_COLUMNS:
            out[module] = row[module]
        marker_row = markers.get(key, {})
        for marker_key, value in marker_row.items():
            if marker_key.startswith("celltype_") or marker_key.startswith("z_celltype_"):
                out[marker_key] = value
        rows.append(out)

    base_fields = [
        "sample",
        "stage",
        "species",
        "barcode",
        "tissue_flag",
        "author_domain",
        "array_row",
        "array_col",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres",
        "mechanical_border_score",
        "immune_fibrotic_activation_score",
        "fibroblast_scar_repair_score",
        "domain34_signed_distance_fullres_px",
        "is_domain34_pixel_boundary",
        "nearest_opposite_domain34_distance_fullres_px",
        "score_only_state_label",
        "is_domain34_graph_boundary",
        "graph_degree_domain34",
        "he_brightness",
        "he_saturation",
        "he_blue_red_optical_density_proxy",
        *CORE_SCORE_COLUMNS,
    ]
    extra_fields = sorted({field for row in rows for field in row if field not in base_fields})
    write_lf_tsv(OUT_MASTER, rows, [*base_fields, *extra_fields])


def main() -> None:
    signatures = read_signatures(SIGNATURES_TSV)
    mouse_rows = read_tsv(SPOT_TABLE)
    build_provenance_qc(signatures, mouse_rows)
    build_alternative_scoring(signatures, mouse_rows)
    expression_matched_random_controls(signatures, mouse_rows)
    graph_status = build_graph_boundary(mouse_rows)
    build_histology_audit(mouse_rows)
    build_master_table(mouse_rows, graph_status)
    print("Wrote final tutor-requested rigor upgrade tables")


if __name__ == "__main__":
    main()
