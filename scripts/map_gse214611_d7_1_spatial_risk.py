#!/usr/bin/env python3
"""Map the prototype GSE214611 fibrotic-risk score back onto Visium space."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from project_paths import project_root


ROOT = project_root(__file__)
SPATIAL_METADATA = (
    ROOT
    / "data"
    / "metadata"
    / "gse214611_author"
    / "spatial_object_integrated_metadata.csv"
)
DEFAULT_VISIUM_DIR = ROOT / "data" / "raw" / "gse214611" / "visium" / "GSM6613087_V_d7_1"
OUT_TABLE_DIR = ROOT / "results" / "tables"
OUT_FIGURE_DIR = ROOT / "results" / "figures"
DEFAULT_SAMPLE = "D7_1"

MODEL_COLUMNS = [
    "BZ2_genes",
    "BZ1_genes",
    "RZ_genes",
    "X33_Fib_Postn",
    "X36_Fib_Rep",
    "X11_Monos",
    "X14_Mac_C1qa",
]


def safe_float(value: str) -> float:
    try:
        return float(value) if value != "" else 0.0
    except ValueError:
        return 0.0


def read_metadata(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def zscore_columns(rows: list[dict[str, str]], columns: list[str]) -> dict[str, list[float]]:
    raw: dict[str, list[float]] = {col: [] for col in columns}
    for row in rows:
        for col in columns:
            raw[col].append(safe_float(row.get(col, "")))

    zscores: dict[str, list[float]] = {}
    for col, values in raw.items():
        mean = sum(values) / len(values)
        sd = math.sqrt(sum((value - mean) ** 2 for value in values) / len(values)) or 1.0
        zscores[col] = [(value - mean) / sd for value in values]
    return zscores


def add_scores(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    z = zscore_columns(rows, MODEL_COLUMNS)
    scored = []

    for idx, row in enumerate(rows):
        out = dict(row)
        immune_proxy = 0.5 * z["X11_Monos"][idx] + 0.5 * z["X14_Mac_C1qa"][idx]
        fibrotic_risk = (
            z["BZ2_genes"][idx]
            + z["X33_Fib_Postn"][idx]
            + immune_proxy
            - 0.5 * z["RZ_genes"][idx]
        )
        border_activation = 0.5 * z["BZ1_genes"][idx] + 0.5 * z["BZ2_genes"][idx]
        out["barcode"] = row[""].rsplit("_", 1)[0]
        out["prototype_immune_proxy"] = f"{immune_proxy:.8g}"
        out["prototype_fibrotic_risk"] = f"{fibrotic_risk:.8g}"
        out["prototype_border_activation"] = f"{border_activation:.8g}"
        out["prototype_repair_proxy"] = f"{z['X36_Fib_Rep'][idx]:.8g}"
        scored.append(out)
    return scored


def read_positions(path: Path) -> dict[str, dict[str, str]]:
    positions: dict[str, dict[str, str]] = {}
    columns = [
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres",
    ]
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        for values in reader:
            row = dict(zip(columns, values, strict=True))
            positions[row["barcode"]] = row
    return positions


def spatial_file(visium_dir: Path, filename: str) -> Path:
    candidates = [visium_dir / "spatial" / filename, visium_dir / filename]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Could not find {filename} in {visium_dir} or {visium_dir / 'spatial'}")


def merge_sample_rows(
    scored_rows: list[dict[str, str]],
    positions: dict[str, dict[str, str]],
    sample: str,
) -> list[dict[str, str]]:
    merged = []
    for row in scored_rows:
        if row["orig.ident"] != sample:
            continue
        barcode = row["barcode"]
        if barcode not in positions:
            continue
        out = dict(row)
        out.update(positions[barcode])
        merged.append(out)
    return merged


def write_rows(rows: list[dict[str, str]], path: Path) -> None:
    if not rows:
        raise ValueError("No rows to write")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_domain_summary(rows: list[dict[str, str]], path: Path) -> None:
    score_columns = [
        "prototype_fibrotic_risk",
        "prototype_repair_proxy",
        "prototype_immune_proxy",
        "prototype_border_activation",
        "BZ1_genes",
        "BZ2_genes",
        "RZ_genes",
        "X33_Fib_Postn",
        "X36_Fib_Rep",
        "CM_Score",
    ]
    counts: Counter[str] = Counter()
    sums: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))

    for row in rows:
        domain = row["annotated"]
        counts[domain] += 1
        for col in score_columns:
            sums[domain][col] += safe_float(row[col])

    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["annotated", "n_spots", *[f"mean_{col}" for col in score_columns]])
        for domain, n_spots in sorted(counts.items(), key=lambda item: item[0]):
            means = [f"{sums[domain][col] / n_spots:.6g}" for col in score_columns]
            writer.writerow([domain, n_spots, *means])


def write_top_risk_rows(rows: list[dict[str, str]], path: Path, limit: int = 50) -> None:
    keep_columns = [
        "barcode",
        "annotated",
        "prototype_fibrotic_risk",
        "prototype_repair_proxy",
        "prototype_immune_proxy",
        "prototype_border_activation",
        "BZ1_genes",
        "BZ2_genes",
        "RZ_genes",
        "X33_Fib_Postn",
        "X36_Fib_Rep",
        "CM_Score",
        "array_row",
        "array_col",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres",
    ]
    top_rows = sorted(rows, key=lambda row: safe_float(row["prototype_fibrotic_risk"]), reverse=True)[
        :limit
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keep_columns, delimiter="\t")
        writer.writeheader()
        writer.writerows({col: row[col] for col in keep_columns} for row in top_rows)


def plot_maps(rows: list[dict[str, str]], visium_dir: Path, sample: str, figure_path: Path) -> None:
    from PIL import Image, ImageDraw

    image_path = spatial_file(visium_dir, "tissue_lowres_image.png")
    scale_path = spatial_file(visium_dir, "scalefactors_json.json")
    tissue_image = Image.open(image_path).convert("RGB")
    with scale_path.open() as handle:
        scale = json.load(handle)["tissue_lowres_scalef"]

    points = [
        {
            "x": safe_float(row["pxl_col_in_fullres"]) * scale,
            "y": safe_float(row["pxl_row_in_fullres"]) * scale,
            "risk": safe_float(row["prototype_fibrotic_risk"]),
            "domain": int(row["annotated"]),
        }
        for row in rows
    ]
    risk_values = sorted(point["risk"] for point in points)
    risk_min = percentile(risk_values, 0.02)
    risk_max = percentile(risk_values, 0.98)

    risk_panel = draw_risk_panel(tissue_image, points, risk_min, risk_max)
    domain_panel = draw_domain_panel(tissue_image, points)
    write_combined_panel(
        risk_panel,
        domain_panel,
        figure_path,
        sample,
        risk_min,
        risk_max,
    )


def percentile(values: list[float], fraction: float) -> float:
    if not values:
        return 0.0
    index = max(0, min(len(values) - 1, round((len(values) - 1) * fraction)))
    return values[index]


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[idx : idx + 2], 16) for idx in (0, 2, 4))


def interpolate_rgb(
    left: tuple[int, int, int],
    right: tuple[int, int, int],
    weight: float,
) -> tuple[int, int, int]:
    return tuple(round(left[idx] + (right[idx] - left[idx]) * weight) for idx in range(3))


def magma_like(value: float, min_value: float, max_value: float) -> tuple[int, int, int]:
    palette = [
        hex_to_rgb("#000004"),
        hex_to_rgb("#3b0f70"),
        hex_to_rgb("#8c2981"),
        hex_to_rgb("#de4968"),
        hex_to_rgb("#fe9f6d"),
        hex_to_rgb("#fcfdbf"),
    ]
    if max_value <= min_value:
        return palette[-1]
    norm = (value - min_value) / (max_value - min_value)
    norm = max(0.0, min(1.0, norm))
    scaled = norm * (len(palette) - 1)
    left_idx = min(len(palette) - 2, int(scaled))
    return interpolate_rgb(palette[left_idx], palette[left_idx + 1], scaled - left_idx)


def dim_background(image: "Image.Image") -> "Image.Image":
    from PIL import Image

    white = Image.new("RGB", image.size, "white")
    return Image.blend(image, white, 0.24)


def draw_risk_panel(
    tissue_image: "Image.Image",
    points: list[dict[str, float]],
    risk_min: float,
    risk_max: float,
) -> "Image.Image":
    from PIL import ImageDraw

    panel = dim_background(tissue_image).convert("RGBA")
    draw = ImageDraw.Draw(panel, "RGBA")
    radius = max(2, round(min(tissue_image.size) / 170))

    for point in sorted(points, key=lambda item: item["risk"]):
        x = point["x"]
        y = point["y"]
        color = magma_like(point["risk"], risk_min, risk_max)
        draw.ellipse(
            (x - radius, y - radius, x + radius, y + radius),
            fill=(*color, 215),
        )
    return panel.convert("RGB")


def draw_domain_panel(tissue_image: "Image.Image", points: list[dict[str, float]]) -> "Image.Image":
    from PIL import ImageDraw

    panel = dim_background(tissue_image).convert("RGBA")
    draw = ImageDraw.Draw(panel, "RGBA")
    radius = max(2, round(min(tissue_image.size) / 170))
    colors = {
        1: hex_to_rgb("#3f7cac"),
        2: hex_to_rgb("#77b28c"),
        3: hex_to_rgb("#f2c14e"),
        4: hex_to_rgb("#d95d39"),
    }

    for point in points:
        x = point["x"]
        y = point["y"]
        color = colors.get(int(point["domain"]), hex_to_rgb("#666666"))
        draw.ellipse(
            (x - radius, y - radius, x + radius, y + radius),
            fill=(*color, 215),
        )
    return panel.convert("RGB")


def centered_text(draw: "ImageDraw.ImageDraw", xy: tuple[int, int], text: str) -> None:
    bbox = draw.textbbox((0, 0), text)
    width = bbox[2] - bbox[0]
    draw.text((xy[0] - width / 2, xy[1]), text, fill=(35, 35, 35))


def write_combined_panel(
    risk_panel: "Image.Image",
    domain_panel: "Image.Image",
    figure_path: Path,
    sample: str,
    risk_min: float,
    risk_max: float,
) -> None:
    from PIL import Image, ImageDraw

    margin = 24
    gap = 28
    title_height = 34
    legend_height = 60
    panel_width, panel_height = risk_panel.size
    canvas_width = panel_width * 2 + gap + margin * 2
    canvas_height = title_height + panel_height + legend_height + margin * 2
    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")
    draw = ImageDraw.Draw(canvas)

    left_x = margin
    right_x = margin + panel_width + gap
    panel_y = margin + title_height
    canvas.paste(risk_panel, (left_x, panel_y))
    canvas.paste(domain_panel, (right_x, panel_y))

    centered_text(draw, (left_x + panel_width // 2, margin + 5), f"{sample} prototype fibrotic-risk")
    centered_text(draw, (right_x + panel_width // 2, margin + 5), f"{sample} author spatial domains")

    legend_y = panel_y + panel_height + 18
    colorbar_x = left_x + 42
    colorbar_width = panel_width - 84
    colorbar_height = 12
    for idx in range(colorbar_width):
        value = risk_min + (risk_max - risk_min) * idx / max(1, colorbar_width - 1)
        color = magma_like(value, risk_min, risk_max)
        draw.line(
            (colorbar_x + idx, legend_y, colorbar_x + idx, legend_y + colorbar_height),
            fill=color,
        )
    draw.rectangle(
        (colorbar_x, legend_y, colorbar_x + colorbar_width, legend_y + colorbar_height),
        outline=(90, 90, 90),
    )
    draw.text((colorbar_x, legend_y + 18), f"{risk_min:.2f}", fill=(35, 35, 35))
    max_label = f"{risk_max:.2f}"
    max_bbox = draw.textbbox((0, 0), max_label)
    draw.text(
        (colorbar_x + colorbar_width - (max_bbox[2] - max_bbox[0]), legend_y + 18),
        max_label,
        fill=(35, 35, 35),
    )

    domain_colors = {
        "1 remote-like": "#3f7cac",
        "2 transition": "#77b28c",
        "3 border-zone": "#f2c14e",
        "4 fibrotic hotspot": "#d95d39",
    }
    domain_x = right_x + 34
    for idx, (label, color) in enumerate(domain_colors.items()):
        x = domain_x + (idx % 2) * 155
        y = legend_y + (idx // 2) * 24
        draw.rectangle((x, y, x + 13, y + 13), fill=hex_to_rgb(color))
        draw.text((x + 19, y - 1), label, fill=(35, 35, 35))

    figure_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(figure_path)


def sample_slug(sample: str) -> str:
    return sample.lower().replace("/", "_").replace(" ", "_")


def run_mapping(sample: str, visium_dir: Path, min_spots: int = 2000) -> dict[str, str | int]:
    OUT_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    OUT_FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    scored_rows = add_scores(read_metadata(SPATIAL_METADATA))
    positions = read_positions(spatial_file(visium_dir, "tissue_positions_list.csv"))
    merged = merge_sample_rows(scored_rows, positions, sample)

    if len(merged) < min_spots:
        raise RuntimeError(f"Expected at least {min_spots} matched {sample} spots; matched {len(merged)}")

    slug = sample_slug(sample)
    spatial_table = OUT_TABLE_DIR / f"gse214611_{slug}_spatial_risk_map.tsv"
    domain_summary = OUT_TABLE_DIR / f"gse214611_{slug}_domain_summary.tsv"
    top_risk = OUT_TABLE_DIR / f"gse214611_{slug}_top_risk_spots.tsv"
    figure = OUT_FIGURE_DIR / f"gse214611_{slug}_spatial_risk_domain_map.png"

    write_rows(merged, spatial_table)
    write_domain_summary(merged, domain_summary)
    write_top_risk_rows(merged, top_risk)
    plot_maps(merged, visium_dir, sample, figure)

    domains = ",".join(sorted({row["annotated"] for row in merged}))
    return {
        "sample": sample,
        "visium_dir": str(visium_dir),
        "n_spots": len(merged),
        "domains": domains,
        "spatial_table": str(spatial_table),
        "domain_summary": str(domain_summary),
        "top_risk": str(top_risk),
        "figure": str(figure),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", default=DEFAULT_SAMPLE, help="orig.ident value in the metadata")
    parser.add_argument(
        "--visium-dir",
        type=Path,
        default=DEFAULT_VISIUM_DIR,
        help="Unpacked 10x Visium directory containing spatial/ and filtered_feature_bc_matrix/",
    )
    parser.add_argument("--min-spots", type=int, default=2000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = run_mapping(args.sample, args.visium_dir, args.min_spots)
    print(f"Matched {result['sample']} spots: {result['n_spots']}")
    print(f"Unique domains: {result['domains']}")
    print(f"Figure: {result['figure']}")


if __name__ == "__main__":
    main()
