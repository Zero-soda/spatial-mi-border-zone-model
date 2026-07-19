#!/usr/bin/env python3
"""Score the GSE214611 human STEMI Visium sample with project signatures."""

from __future__ import annotations

import csv
import json
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
from PIL import Image, ImageDraw  # noqa: E402

from map_gse214611_d7_1_spatial_risk import dim_background, magma_like, percentile, safe_float  # noqa: E402
from score_gse214611_visium_signatures import CORE_SCORE_COLUMNS, decode  # noqa: E402


SIGNATURES_TSV = Path(__file__).resolve().parents[1] / "config" / "spatial_cardiac_border_zone_signatures.tsv"
HUMAN_VISIUM_DIR = ROOT / "data" / "raw" / "gse214611" / "visium" / "GSM6613090_V_Human_STEMI"
OUT_TABLE_DIR = ROOT / "results" / "tables"
OUT_FIGURE_DIR = ROOT / "results" / "figures"

SCORE_COLUMNS = [
    ("mechanical_border_score", "mechanical border"),
    ("immune_fibrotic_activation_score", "immune-fibrotic activation"),
    ("fibroblast_scar_repair_score", "fibroblast-scar repair"),
]


def read_human_signatures(path: Path) -> dict[str, list[str]]:
    signatures: dict[str, list[str]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            genes = [gene.strip() for gene in row["genes_human"].split(";") if gene.strip()]
            signatures[row["signature_id"]] = genes
    return signatures


def score_10x_h5(matrix_h5: Path, signatures: dict[str, list[str]]) -> tuple[list[str], dict[str, dict[str, float]], dict[str, int]]:
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

    detected_counts = {key: len(value) for key, value in signature_gene_indices.items()}
    return barcodes, scores_by_barcode, detected_counts


def read_positions(path: Path) -> dict[str, dict[str, str]]:
    positions: dict[str, dict[str, str]] = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if not row or row[0] == "barcode":
                continue
            barcode, in_tissue, array_row, array_col, pxl_row, pxl_col = row[:6]
            if in_tissue != "1":
                continue
            positions[barcode] = {
                "array_row": array_row,
                "array_col": array_col,
                "pxl_row_in_fullres": pxl_row,
                "pxl_col_in_fullres": pxl_col,
            }
    return positions


def add_zscores(rows: list[dict[str, str]], source_columns: list[str]) -> None:
    for column in source_columns:
        values = [safe_float(row[column]) for row in rows]
        mean = sum(values) / len(values)
        sd = math.sqrt(sum((value - mean) ** 2 for value in values) / len(values)) or 1.0
        for row, value in zip(rows, values, strict=True):
            row[f"z_{column}"] = f"{(value - mean) / sd:.8g}"


def add_composite_scores(rows: list[dict[str, str]]) -> None:
    add_zscores(rows, CORE_SCORE_COLUMNS)
    for row in rows:
        row["mechanical_border_score"] = f"{safe_float(row['z_CM_BZ1_TRANSITION']) + safe_float(row['z_CM_BZ2_MECHANICAL_EDGE']):.8g}"
        row["immune_fibrotic_activation_score"] = f"{safe_float(row['z_CCR2_IL1B_MYLOID']) + safe_float(row['z_TGFB_SIGNALING']) + safe_float(row['z_FAP_POSTN_PATHO_FIBROBLAST']):.8g}"
        row["fibroblast_scar_repair_score"] = f"{safe_float(row['z_FAP_POSTN_PATHO_FIBROBLAST']) + safe_float(row['z_ECM_REMODELING']) + safe_float(row['z_CTHRC1_REPARATIVE_CF']) + safe_float(row['z_MYOFIBROBLAST_CONTRACTILE']):.8g}"


def write_rows(rows: list[dict[str, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: list[dict[str, str]], detected_counts: dict[str, int]) -> None:
    summary_path = OUT_TABLE_DIR / "gse214611_human_stemi_signature_summary.tsv"
    score_fields = [column for column, _ in SCORE_COLUMNS]
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["sample", "GSM6613090_V_Human_STEMI"])
        writer.writerow(["n_spots", len(rows)])
        for signature_id in CORE_SCORE_COLUMNS:
            writer.writerow([f"detected_genes_{signature_id}", detected_counts.get(signature_id, 0)])
        for column in score_fields:
            values = sorted(safe_float(row[column]) for row in rows)
            writer.writerow([f"mean_{column}", f"{sum(values) / len(values):.6g}"])
            writer.writerow([f"p90_{column}", f"{percentile(values, 0.90):.6g}"])
            writer.writerow([f"p95_{column}", f"{percentile(values, 0.95):.6g}"])


def draw_score_panel(rows: list[dict[str, str]], column: str, label: str, score_min: float, score_max: float) -> Image.Image:
    image = Image.open(HUMAN_VISIUM_DIR / "spatial" / "tissue_lowres_image.png").convert("RGB")
    with (HUMAN_VISIUM_DIR / "spatial" / "scalefactors_json.json").open() as handle:
        scale = json.load(handle)["tissue_lowres_scalef"]

    panel = dim_background(image).convert("RGBA")
    draw = ImageDraw.Draw(panel, "RGBA")
    radius = max(2, round(min(image.size) / 165))
    points = [
        (
            safe_float(row["pxl_col_in_fullres"]) * scale,
            safe_float(row["pxl_row_in_fullres"]) * scale,
            safe_float(row[column]),
        )
        for row in rows
    ]
    for x, y, score in sorted(points, key=lambda item: item[2]):
        color = magma_like(score, score_min, score_max)
        draw.ellipse((x - radius, y - radius, x + radius, y + radius), fill=(*color, 218))

    out = panel.convert("RGB")
    header = Image.new("RGB", (out.width, out.height + 62), "white")
    header.paste(out, (0, 30))
    header_draw = ImageDraw.Draw(header)
    header_draw.text((10, 8), label, fill=(25, 25, 25))
    bar_x, bar_y, bar_w, bar_h = 24, out.height + 42, out.width - 48, 10
    for idx in range(bar_w):
        value = score_min + (score_max - score_min) * idx / max(1, bar_w - 1)
        header_draw.line((bar_x + idx, bar_y, bar_x + idx, bar_y + bar_h), fill=magma_like(value, score_min, score_max))
    header_draw.rectangle((bar_x, bar_y, bar_x + bar_w, bar_y + bar_h), outline=(90, 90, 90))
    header_draw.text((bar_x, bar_y + 15), f"{score_min:.2f}", fill=(35, 35, 35))
    max_label = f"{score_max:.2f}"
    bbox = header_draw.textbbox((0, 0), max_label)
    header_draw.text((bar_x + bar_w - (bbox[2] - bbox[0]), bar_y + 15), max_label, fill=(35, 35, 35))
    return header


def write_figure(rows: list[dict[str, str]]) -> Path:
    limits = {}
    for column, _ in SCORE_COLUMNS:
        values = sorted(safe_float(row[column]) for row in rows)
        limits[column] = (percentile(values, 0.02), percentile(values, 0.98))

    panels = [
        draw_score_panel(rows, column, label, *limits[column])
        for column, label in SCORE_COLUMNS
    ]
    gap = 24
    margin = 24
    title_height = 34
    panel_width = max(panel.width for panel in panels)
    panel_height = max(panel.height for panel in panels)
    canvas = Image.new(
        "RGB",
        (len(panels) * panel_width + (len(panels) - 1) * gap + 2 * margin, panel_height + title_height + 2 * margin),
        "white",
    )
    draw = ImageDraw.Draw(canvas)
    draw.text((margin, margin), "GSM6613090 human STEMI Visium signature transfer", fill=(20, 20, 20))
    for idx, panel in enumerate(panels):
        x = margin + idx * (panel_width + gap)
        canvas.paste(panel, (x, margin + title_height))

    out_path = OUT_FIGURE_DIR / "gse214611_human_stemi_signature_transfer.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out_path)
    return out_path


def main() -> None:
    signatures = read_human_signatures(SIGNATURES_TSV)
    _, scores_by_barcode, detected_counts = score_10x_h5(HUMAN_VISIUM_DIR / "filtered_feature_bc_matrix.h5", signatures)
    positions = read_positions(HUMAN_VISIUM_DIR / "spatial" / "tissue_positions_list.csv")

    rows: list[dict[str, str]] = []
    for barcode, position in positions.items():
        if barcode not in scores_by_barcode:
            continue
        row = {
            "sample": "GSM6613090_V_Human_STEMI",
            "barcode": barcode,
            **position,
        }
        for signature_id in CORE_SCORE_COLUMNS:
            row[signature_id] = f"{scores_by_barcode[barcode][signature_id]:.8g}"
        rows.append(row)

    if not rows:
        raise RuntimeError("No human STEMI spots matched between expression matrix and spatial positions")

    add_composite_scores(rows)
    write_rows(rows, OUT_TABLE_DIR / "gse214611_human_stemi_signature_scores_by_spot.tsv")
    summarize(rows, detected_counts)
    out_path = write_figure(rows)
    print(f"[summary] human STEMI spots: {len(rows)}")
    print(f"[figure] {out_path}")


if __name__ == "__main__":
    main()
