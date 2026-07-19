#!/usr/bin/env python3
"""Score GSE214611 D3/D7 Visium spots with project gene signatures."""

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

import h5py  # noqa: E402
import numpy as np  # noqa: E402

from batch_map_gse214611_mi_spatial_risk import SAMPLES, locate_visium_dir, sample_stage  # noqa: E402
from map_gse214611_d7_1_spatial_risk import sample_slug, safe_float  # noqa: E402


SIGNATURES_TSV = Path(__file__).resolve().parents[1] / "config" / "spatial_cardiac_border_zone_signatures.tsv"
RAW_VISIUM_DIR = ROOT / "data" / "raw" / "gse214611" / "visium"
OUT_TABLE_DIR = ROOT / "results" / "tables"

CORE_SCORE_COLUMNS = [
    "CM_BZ1_TRANSITION",
    "CM_BZ2_MECHANICAL_EDGE",
    "FAP_POSTN_PATHO_FIBROBLAST",
    "CTHRC1_REPARATIVE_CF",
    "MYOFIBROBLAST_CONTRACTILE",
    "CCR2_IL1B_MYLOID",
    "TGFB_SIGNALING",
    "ECM_REMODELING",
    "ENDOTHELIAL_INFLAMMATORY_FIBROSIS",
    "HYPOXIA_ISCHEMIA",
    "CM_ELECTRICAL_CALCIUM_FUNCTION",
    "INJURY_STRESS_IMMEDIATE_EARLY",
]


def decode(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode()
    return str(value)


def read_signatures(path: Path) -> dict[str, list[str]]:
    signatures: dict[str, list[str]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            genes = [gene.strip() for gene in row["genes_mouse"].split(";") if gene.strip()]
            signatures[row["signature_id"]] = genes
    return signatures


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


def score_10x_h5(
    matrix_h5: Path,
    signatures: dict[str, list[str]],
) -> dict[str, dict[str, float]]:
    with h5py.File(matrix_h5, "r") as handle:
        matrix = handle["matrix"]
        barcodes = [decode(value) for value in matrix["barcodes"][:]]
        gene_names = [decode(value) for value in matrix["features"]["name"][:]]
        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]

    gene_to_index = {gene.casefold(): idx for idx, gene in enumerate(gene_names)}
    signature_gene_indices: dict[str, set[int]] = {}
    for signature_id, genes in signatures.items():
        signature_gene_indices[signature_id] = {
            gene_to_index[gene.casefold()] for gene in genes if gene.casefold() in gene_to_index
        }

    gene_to_signatures: dict[int, list[str]] = defaultdict(list)
    for signature_id, gene_indices in signature_gene_indices.items():
        for gene_index in gene_indices:
            gene_to_signatures[gene_index].append(signature_id)

    scores_by_barcode: dict[str, dict[str, float]] = {}
    for barcode_idx, barcode in enumerate(barcodes):
        start = int(indptr[barcode_idx])
        end = int(indptr[barcode_idx + 1])
        counts = data[start:end].astype(float, copy=False)
        genes = indices[start:end]
        library_size = float(counts.sum())
        signature_sums = {signature_id: 0.0 for signature_id in signatures}

        if library_size > 0:
            normalized = np.log1p(counts / library_size * 10000.0)
            for gene_index, value in zip(genes, normalized, strict=True):
                for signature_id in gene_to_signatures.get(int(gene_index), []):
                    signature_sums[signature_id] += float(value)

        scores_by_barcode[barcode] = {
            signature_id: (
                signature_sums[signature_id] / len(signature_gene_indices[signature_id])
                if signature_gene_indices[signature_id]
                else 0.0
            )
            for signature_id in signatures
        }

    return scores_by_barcode


def read_spatial_risk_rows(sample: str) -> list[dict[str, str]]:
    path = OUT_TABLE_DIR / f"gse214611_{sample_slug(sample)}_spatial_risk_map.tsv"
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def add_zscore_columns(rows: list[dict[str, str]], source_columns: list[str]) -> None:
    for column in source_columns:
        values = [safe_float(row[column]) for row in rows]
        mean = sum(values) / len(values)
        sd = math.sqrt(sum((value - mean) ** 2 for value in values) / len(values)) or 1.0
        for row, value in zip(rows, values, strict=True):
            row[f"z_{column}"] = f"{(value - mean) / sd:.8g}"


def add_composite_scores(rows: list[dict[str, str]]) -> None:
    add_zscore_columns(rows, CORE_SCORE_COLUMNS)
    for row in rows:
        signature_fibrotic_risk = (
            safe_float(row["z_FAP_POSTN_PATHO_FIBROBLAST"])
            + safe_float(row["z_ECM_REMODELING"])
            + safe_float(row["z_CCR2_IL1B_MYLOID"])
            + safe_float(row["z_CM_BZ2_MECHANICAL_EDGE"])
            + safe_float(row["z_TGFB_SIGNALING"])
            - safe_float(row["z_CTHRC1_REPARATIVE_CF"])
        )
        mechanical_border_score = (
            safe_float(row["z_CM_BZ1_TRANSITION"])
            + safe_float(row["z_CM_BZ2_MECHANICAL_EDGE"])
        )
        fibroblast_scar_score = (
            safe_float(row["z_FAP_POSTN_PATHO_FIBROBLAST"])
            + safe_float(row["z_ECM_REMODELING"])
            + safe_float(row["z_CTHRC1_REPARATIVE_CF"])
            + safe_float(row["z_MYOFIBROBLAST_CONTRACTILE"])
        )
        row["signature_fibrotic_risk"] = f"{signature_fibrotic_risk:.8g}"
        row["signature_mechanical_border_score"] = f"{mechanical_border_score:.8g}"
        row["signature_fibroblast_scar_score"] = f"{fibroblast_scar_score:.8g}"


def write_rows(rows: list[dict[str, str]], out_path: Path) -> None:
    if not rows:
        raise ValueError("No rows to write")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: list[dict[str, str]], group_columns: list[str], out_path: Path) -> None:
    score_columns = [
        *CORE_SCORE_COLUMNS,
        "signature_fibrotic_risk",
        "signature_mechanical_border_score",
        "signature_fibroblast_scar_score",
        "prototype_fibrotic_risk",
        "prototype_repair_proxy",
        "prototype_immune_proxy",
        "prototype_border_activation",
    ]
    sums: dict[tuple[str, ...], dict[str, float]] = defaultdict(lambda: defaultdict(float))
    counts: dict[tuple[str, ...], int] = defaultdict(int)

    for row in rows:
        key = tuple(row[col] for col in group_columns)
        counts[key] += 1
        for col in score_columns:
            sums[key][col] += safe_float(row[col])

    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([*group_columns, "n_spots", *[f"mean_{col}" for col in score_columns]])
        for key in sorted(counts):
            values = [f"{sums[key][col] / counts[key]:.6g}" for col in score_columns]
            writer.writerow([*key, counts[key], *values])


def main() -> None:
    OUT_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    signatures = read_signatures(SIGNATURES_TSV)
    all_rows: list[dict[str, str]] = []

    for sample, config in SAMPLES.items():
        visium_dir = locate_visium_dir(config.unpack_dir)
        matrix_h5 = find_matrix_h5(visium_dir)
        print(f"[score] {sample}: {matrix_h5}")
        scores_by_barcode = score_10x_h5(matrix_h5, signatures)
        spatial_rows = read_spatial_risk_rows(sample)

        matched = 0
        for row in spatial_rows:
            barcode = row["barcode"]
            if barcode not in scores_by_barcode:
                continue
            matched += 1
            out = {
                "sample": sample,
                "stage": sample_stage(sample),
                "barcode": barcode,
                "annotated": row["annotated"],
                "array_row": row["array_row"],
                "array_col": row["array_col"],
                "pxl_row_in_fullres": row["pxl_row_in_fullres"],
                "pxl_col_in_fullres": row["pxl_col_in_fullres"],
                "prototype_fibrotic_risk": row["prototype_fibrotic_risk"],
                "prototype_repair_proxy": row["prototype_repair_proxy"],
                "prototype_immune_proxy": row["prototype_immune_proxy"],
                "prototype_border_activation": row["prototype_border_activation"],
            }
            for signature_id in CORE_SCORE_COLUMNS:
                out[signature_id] = f"{scores_by_barcode[barcode][signature_id]:.8g}"
            all_rows.append(out)

        print(f"[done] {sample}: matched {matched}/{len(spatial_rows)} spots")

    add_composite_scores(all_rows)
    write_rows(all_rows, OUT_TABLE_DIR / "gse214611_d3_d7_signature_scores_by_spot.tsv")
    summarize(
        all_rows,
        ["sample", "annotated"],
        OUT_TABLE_DIR / "gse214611_d3_d7_signature_scores_by_sample_domain.tsv",
    )
    summarize(
        all_rows,
        ["stage", "annotated"],
        OUT_TABLE_DIR / "gse214611_d3_d7_signature_scores_by_stage_domain.tsv",
    )
    print(f"[summary] scored spots: {len(all_rows)}")


if __name__ == "__main__":
    main()
