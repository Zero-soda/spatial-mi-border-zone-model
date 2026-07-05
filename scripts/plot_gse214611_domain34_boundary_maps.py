#!/usr/bin/env python3
"""Plot GSE214611 domain 3/4 contact-boundary maps."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw

from batch_map_gse214611_mi_spatial_risk import SAMPLES, locate_visium_dir
from map_gse214611_d7_1_spatial_risk import ROOT, dim_background, hex_to_rgb, safe_float, spatial_file


SPOT_FEATURES = ROOT / "results" / "tables" / "gse214611_d3_d7_domain34_spatial_features_by_spot.tsv"
OUT_FIGURE_DIR = ROOT / "results" / "figures"
RESAMPLE = getattr(Image, "Resampling", Image).LANCZOS


COLORS = {
    "domain3": "#f2c14e",
    "domain4": "#d95d39",
    "domain3_boundary": "#2f80ed",
    "domain4_boundary": "#00a676",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def rows_by_sample(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["sample"]].append(row)
    return grouped


def load_tissue_image_and_scale(sample: str) -> tuple[Image.Image, float]:
    visium_dir = locate_visium_dir(SAMPLES[sample].unpack_dir)
    image = Image.open(spatial_file(visium_dir, "tissue_lowres_image.png")).convert("RGB")
    with spatial_file(visium_dir, "scalefactors_json.json").open() as handle:
        scale = json.load(handle)["tissue_lowres_scalef"]
    return image, scale


def row_color(row: dict[str, str]) -> tuple[int, int, int] | None:
    is_boundary = row["is_domain34_contact_boundary"] == "True"
    domain = row["annotated"]
    if domain == "3":
        return hex_to_rgb(COLORS["domain3_boundary" if is_boundary else "domain3"])
    if domain == "4":
        return hex_to_rgb(COLORS["domain4_boundary" if is_boundary else "domain4"])
    return None


def draw_boundary_panel(sample: str, rows: list[dict[str, str]]) -> Image.Image:
    tissue_image, scale = load_tissue_image_and_scale(sample)
    panel = dim_background(tissue_image).convert("RGBA")
    draw = ImageDraw.Draw(panel, "RGBA")
    radius = max(2, round(min(tissue_image.size) / 150))

    for row in rows:
        color = row_color(row)
        if color is None:
            continue
        alpha = 235 if row["is_domain34_contact_boundary"] == "True" else 120
        x = safe_float(row["pxl_col_in_fullres"]) * scale
        y = safe_float(row["pxl_row_in_fullres"]) * scale
        draw.ellipse((x - radius, y - radius, x + radius, y + radius), fill=(*color, alpha))
    return panel.convert("RGB")


def write_sample_maps(grouped: dict[str, list[dict[str, str]]]) -> list[tuple[str, Path]]:
    outputs = []
    for sample in SAMPLES:
        panel = draw_boundary_panel(sample, grouped[sample])
        out_path = OUT_FIGURE_DIR / f"gse214611_{sample.lower()}_domain34_boundary_map.png"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        panel.save(out_path, dpi=(600, 600))
        outputs.append((sample, out_path))
        print(f"[boundary-map] {sample}: {out_path}")
    return outputs


def write_contact_sheet(grouped: dict[str, list[dict[str, str]]]) -> Path:
    panels = []
    for sample in SAMPLES:
        panel = draw_boundary_panel(sample, grouped[sample])
        target_width = 850
        target_height = round(panel.height * target_width / panel.width)
        panels.append((sample, panel.resize((target_width, target_height), RESAMPLE)))

    columns = 3
    rows_count = 2
    margin = 36
    gap = 30
    title_height = 38
    legend_height = 66
    panel_width = max(panel.width for _, panel in panels)
    panel_height = max(panel.height for _, panel in panels)
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

    legend_items = [
        ("domain 3", COLORS["domain3"]),
        ("domain 4", COLORS["domain4"]),
        ("domain 3 boundary", COLORS["domain3_boundary"]),
        ("domain 4 boundary", COLORS["domain4_boundary"]),
    ]
    legend_y = canvas.height - margin - 36
    for idx, (label, color) in enumerate(legend_items):
        x = margin + idx * 190
        draw.rectangle((x, legend_y, x + 14, legend_y + 14), fill=hex_to_rgb(color))
        draw.text((x + 20, legend_y - 1), label, fill=(35, 35, 35))

    out_path = OUT_FIGURE_DIR / "gse214611_d3_d7_domain34_boundary_contact_sheet.png"
    canvas.save(out_path, dpi=(600, 600))
    canvas.save(out_path.with_suffix(".tiff"), dpi=(600, 600), compression="tiff_lzw")
    print(f"[contact-sheet] {out_path}")
    return out_path


def main() -> None:
    rows = read_rows(SPOT_FEATURES)
    grouped = rows_by_sample(rows)
    write_sample_maps(grouped)
    write_contact_sheet(grouped)


if __name__ == "__main__":
    main()
