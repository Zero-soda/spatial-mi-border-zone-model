#!/usr/bin/env python3
"""Plot GSE214611 D3/D7 Visium signature score maps."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw

from batch_map_gse214611_mi_spatial_risk import SAMPLES, locate_visium_dir
from map_gse214611_d7_1_spatial_risk import (
    ROOT,
    dim_background,
    hex_to_rgb,
    magma_like,
    percentile,
    safe_float,
    spatial_file,
)


SCORES_BY_SPOT = ROOT / "results" / "tables" / "gse214611_d3_d7_signature_scores_by_spot.tsv"
OUT_FIGURE_DIR = ROOT / "results" / "figures"
RESAMPLE = getattr(Image, "Resampling", Image).LANCZOS

SCORE_COLUMNS = [
    ("signature_mechanical_border_score", "mechanical border"),
    ("signature_fibroblast_scar_score", "fibroblast scar"),
    ("signature_fibrotic_risk", "composite fibrotic risk"),
]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def rows_by_sample(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["sample"]].append(row)
    return grouped


def global_limits(rows: list[dict[str, str]]) -> dict[str, tuple[float, float]]:
    limits = {}
    for column, _ in SCORE_COLUMNS:
        values = sorted(safe_float(row[column]) for row in rows)
        limits[column] = (percentile(values, 0.02), percentile(values, 0.98))
    return limits


def load_tissue_image_and_scale(sample: str) -> tuple[Image.Image, float]:
    visium_dir = locate_visium_dir(SAMPLES[sample].unpack_dir)
    image = Image.open(spatial_file(visium_dir, "tissue_lowres_image.png")).convert("RGB")
    with spatial_file(visium_dir, "scalefactors_json.json").open() as handle:
        scale = json.load(handle)["tissue_lowres_scalef"]
    return image, scale


def score_points(rows: list[dict[str, str]], column: str, scale: float) -> list[dict[str, float]]:
    return [
        {
            "x": safe_float(row["pxl_col_in_fullres"]) * scale,
            "y": safe_float(row["pxl_row_in_fullres"]) * scale,
            "score": safe_float(row[column]),
            "domain": int(row["annotated"]),
        }
        for row in rows
    ]


def draw_score_panel(
    tissue_image: Image.Image,
    points: list[dict[str, float]],
    score_min: float,
    score_max: float,
) -> Image.Image:
    panel = dim_background(tissue_image).convert("RGBA")
    draw = ImageDraw.Draw(panel, "RGBA")
    radius = max(2, round(min(tissue_image.size) / 170))

    for point in sorted(points, key=lambda item: item["score"]):
        color = magma_like(point["score"], score_min, score_max)
        x = point["x"]
        y = point["y"]
        draw.ellipse((x - radius, y - radius, x + radius, y + radius), fill=(*color, 215))
    return panel.convert("RGB")


def centered_text(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str) -> None:
    bbox = draw.textbbox((0, 0), text)
    draw.text((xy[0] - (bbox[2] - bbox[0]) / 2, xy[1]), text, fill=(35, 35, 35))


def draw_colorbar(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    width: int,
    height: int,
    score_min: float,
    score_max: float,
) -> None:
    for idx in range(width):
        value = score_min + (score_max - score_min) * idx / max(1, width - 1)
        draw.line((x + idx, y, x + idx, y + height), fill=magma_like(value, score_min, score_max))
    draw.rectangle((x, y, x + width, y + height), outline=(90, 90, 90))
    draw.text((x, y + height + 5), f"{score_min:.2f}", fill=(35, 35, 35))
    max_label = f"{score_max:.2f}"
    bbox = draw.textbbox((0, 0), max_label)
    draw.text((x + width - (bbox[2] - bbox[0]), y + height + 5), max_label, fill=(35, 35, 35))


def write_sample_triptych(
    sample: str,
    rows: list[dict[str, str]],
    limits: dict[str, tuple[float, float]],
) -> Path:
    tissue_image, scale = load_tissue_image_and_scale(sample)
    panels = []
    for column, label in SCORE_COLUMNS:
        score_min, score_max = limits[column]
        panel = draw_score_panel(tissue_image, score_points(rows, column, scale), score_min, score_max)
        panels.append((column, label, panel, score_min, score_max))

    margin = 24
    gap = 28
    title_height = 34
    legend_height = 48
    panel_width, panel_height = panels[0][2].size
    canvas = Image.new(
        "RGB",
        (
            panel_width * len(panels) + gap * (len(panels) - 1) + margin * 2,
            title_height + panel_height + legend_height + margin * 2,
        ),
        "white",
    )
    draw = ImageDraw.Draw(canvas)

    for idx, (_, label, panel, score_min, score_max) in enumerate(panels):
        x = margin + idx * (panel_width + gap)
        y = margin + title_height
        centered_text(draw, (x + panel_width // 2, margin + 5), f"{sample} {label}")
        canvas.paste(panel, (x, y))
        draw_colorbar(draw, x + 42, y + panel_height + 16, panel_width - 84, 11, score_min, score_max)

    out_path = OUT_FIGURE_DIR / f"gse214611_{sample.lower()}_signature_score_maps.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out_path, dpi=(600, 600))
    return out_path


def write_score_contact_sheet(
    rows_grouped: dict[str, list[dict[str, str]]],
    column: str,
    label: str,
    limits: dict[str, tuple[float, float]],
) -> Path:
    score_min, score_max = limits[column]
    panels = []
    for sample in SAMPLES:
        tissue_image, scale = load_tissue_image_and_scale(sample)
        panel = draw_score_panel(
            tissue_image,
            score_points(rows_grouped[sample], column, scale),
            score_min,
            score_max,
        )
        target_width = 850
        target_height = round(panel.height * target_width / panel.width)
        panels.append((sample, panel.resize((target_width, target_height), RESAMPLE)))

    columns = 3
    rows_count = 2
    margin = 36
    gap = 30
    title_height = 38
    panel_width = max(panel.width for _, panel in panels)
    panel_height = max(panel.height for _, panel in panels)
    legend_height = 62
    canvas = Image.new(
        "RGB",
        (
            columns * panel_width + (columns - 1) * gap + 2 * margin,
            rows_count * (panel_height + title_height) + (rows_count - 1) * gap + legend_height + 2 * margin,
        ),
        "white",
    )
    draw = ImageDraw.Draw(canvas)

    for idx, (sample, panel) in enumerate(panels):
        row_idx = idx // columns
        col_idx = idx % columns
        x = margin + col_idx * (panel_width + gap)
        y = margin + row_idx * (panel_height + title_height + gap)
        draw.text((x, y), sample, fill=(25, 25, 25))
        canvas.paste(panel, (x, y + title_height))

    legend_y = canvas.height - margin - legend_height + 8
    draw.text((margin, legend_y - 20), label, fill=(25, 25, 25))
    draw_colorbar(draw, margin, legend_y, canvas.width - 2 * margin, 12, score_min, score_max)

    slug = column.replace("signature_", "").replace("_score", "")
    out_path = OUT_FIGURE_DIR / f"gse214611_d3_d7_{slug}_contact_sheet.png"
    canvas.save(out_path, dpi=(600, 600))
    canvas.save(out_path.with_suffix(".tiff"), dpi=(600, 600), compression="tiff_lzw")
    return out_path


def main() -> None:
    rows = read_rows(SCORES_BY_SPOT)
    grouped = rows_by_sample(rows)
    limits = global_limits(rows)

    for sample in SAMPLES:
        out_path = write_sample_triptych(sample, grouped[sample], limits)
        print(f"[sample-map] {sample}: {out_path}")

    for column, label in SCORE_COLUMNS:
        out_path = write_score_contact_sheet(grouped, column, label, limits)
        print(f"[contact-sheet] {column}: {out_path}")


if __name__ == "__main__":
    main()
