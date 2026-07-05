#!/usr/bin/env python3
"""Generate the Cardiovascular Research-style human STEMI transfer figure."""

from __future__ import annotations

import json
import sys
import textwrap
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    # Append rather than prepend so the bundled runtime keeps its binary wheels.
    sys.path.append(str(LOCAL_DEPS))

import matplotlib  # noqa: E402

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402


HUMAN_VISIUM_DIR = ROOT / "data" / "raw" / "gse214611" / "visium" / "GSM6613090_V_Human_STEMI"
TABLE_DIR = ROOT / "results" / "tables"
OUT_DIR = ROOT / "results" / "figures"

SPOT_TABLE = TABLE_DIR / "gse214611_human_stemi_signature_scores_by_spot.tsv"
SUMMARY_TABLE = TABLE_DIR / "gse214611_human_stemi_signature_summary.tsv"
SEPARATION_TABLE = TABLE_DIR / "gse214611_human_stemi_score_separation.tsv"

SCORES = [
    ("mechanical_border_score", "Mechanical-border", "#2F6FB5"),
    ("immune_fibrotic_activation_score", "Immune-fibrotic", "#D99032"),
    ("fibroblast_scar_repair_score", "Fibroblast-scar repair", "#3A9B65"),
]

OUTPUT_MODULES = [
    (
        "Mechanical-border",
        "#2F6FB5",
        [("CM-BZ1", "CM_BZ1_TRANSITION"), ("mechanical edge", "CM_BZ2_MECHANICAL_EDGE")],
    ),
    (
        "Immune-fibrotic",
        "#D99032",
        [
            ("CCR2/IL1B", "CCR2_IL1B_MYLOID"),
            ("TGF-beta", "TGFB_SIGNALING"),
            ("FAP/POSTN", "FAP_POSTN_PATHO_FIBROBLAST"),
        ],
    ),
    (
        "Fibroblast-scar repair",
        "#3A9B65",
        [
            ("FAP/POSTN", "FAP_POSTN_PATHO_FIBROBLAST"),
            ("ECM", "ECM_REMODELING"),
            ("CTHRC1", "CTHRC1_REPARATIVE_CF"),
            ("myofibroblast", "MYOFIBROBLAST_CONTRACTILE"),
        ],
    ),
]

PAIR_LABELS = {
    frozenset(("mechanical_border_score", "immune_fibrotic_activation_score")): "Mechanical vs immune",
    frozenset(("mechanical_border_score", "fibroblast_scar_repair_score")): "Mechanical vs scar",
    frozenset(("immune_fibrotic_activation_score", "fibroblast_scar_repair_score")): "Immune vs scar",
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 6.4,
            "axes.titlesize": 7.4,
            "axes.labelsize": 6.8,
            "xtick.labelsize": 6.2,
            "ytick.labelsize": 6.2,
            "axes.linewidth": 0.55,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def read_summary() -> dict[str, str]:
    summary = pd.read_csv(SUMMARY_TABLE, sep="\t")
    return dict(zip(summary["metric"], summary["value"], strict=True))


def label_panel(ax: plt.Axes, label: str, x: float = -0.08, y: float = 1.08) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        fontweight="bold",
        color="#111827",
    )


def add_box(ax: plt.Axes, xy: tuple[float, float], width: float, height: float, text: str, color: str) -> None:
    box = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.018,rounding_size=0.035",
        linewidth=0.85,
        edgecolor=color,
        facecolor="#FFFFFF",
    )
    ax.add_patch(box)
    ax.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=6.4,
        color="#1F2937",
        linespacing=1.15,
    )


def plot_transfer_design(ax: plt.Axes, n_spots: str) -> None:
    ax.set_axis_off()
    label_panel(ax, "a", x=-0.03, y=1.02)
    ax.text(
        0.02,
        0.98,
        "Fixed-signature transfer design",
        ha="left",
        va="top",
        fontsize=7.8,
        fontweight="bold",
        color="#111827",
    )

    xs = [0.025, 0.375, 0.725]
    colors = ["#2F6FB5", "#6B5FB5", "#3A9B65"]
    texts = [
        "Mouse MI spatial\nstate model\nD3/D7 Visium",
        "Human gene-symbol\nscoring\n(no retraining)",
        "Human STEMI\nGSM6613090\n1,551 spots",
    ]
    for x, color, text in zip(xs, colors, texts, strict=True):
        add_box(ax, (x, 0.48), 0.235, 0.30, text, color)

    for start, end in [(0.260, 0.375), (0.610, 0.725)]:
        arrow = FancyArrowPatch(
            (start, 0.63),
            (end, 0.63),
            arrowstyle="-|>",
            mutation_scale=8,
            linewidth=1.0,
            color="#6B7280",
        )
        ax.add_patch(arrow)

    ax.text(
        0.03,
        0.335,
        f"Primary human readout: {n_spots} tissue spots scored in original tissue coordinates",
        ha="left",
        va="center",
        fontsize=6.3,
        color="#4B5563",
    )
    legend_items = [
        ("Mechanical-border", "#2F6FB5", 0.04),
        ("Immune-fibrotic", "#D99032", 0.36),
        ("Fibroblast-scar repair", "#3A9B65", 0.64),
    ]
    for label, color, x in legend_items:
        ax.scatter([x], [0.185], s=20, color=color, edgecolors="#FFFFFF", linewidths=0.4, zorder=3)
        ax.text(x + 0.028, 0.185, label, ha="left", va="center", fontsize=6.2, color="#1F2937")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def plot_gene_coverage(ax: plt.Axes, summary: dict[str, str]) -> None:
    label_panel(ax, "b", x=-0.17, y=1.08)
    ax.set_axis_off()
    ax.text(
        0.0,
        1.0,
        "Detected human modules",
        ha="left",
        va="top",
        fontsize=7.8,
        fontweight="bold",
        color="#111827",
    )
    y = 0.79
    for output, color, modules in OUTPUT_MODULES:
        ax.scatter([0.025], [y + 0.006], s=28, color=color, edgecolors="#FFFFFF", linewidths=0.4, zorder=3)
        ax.text(0.06, y + 0.01, output, ha="left", va="center", fontsize=6.6, fontweight="bold", color="#1F2937")
        module_text = "; ".join(
            f"{module_label} {int(float(summary[f'detected_genes_{module_id}']))}"
            for module_label, module_id in modules
        )
        module_text = textwrap.fill(module_text, width=31, break_long_words=False)
        ax.text(
            0.06,
            y - 0.10,
            module_text,
            ha="left",
            va="top",
            fontsize=5.95,
            color="#4B5563",
            linespacing=1.16,
        )
        y -= 0.295
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def score_cmap(color: str) -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list("score_cmap", ["#F4F6F8", color])


def plot_spatial_map(ax: plt.Axes, df: pd.DataFrame, column: str, title: str, color: str, label: str | None = None) -> None:
    if label:
        label_panel(ax, label, x=-0.12, y=1.08)
    image = plt.imread(HUMAN_VISIUM_DIR / "spatial" / "tissue_lowres_image.png")
    with (HUMAN_VISIUM_DIR / "spatial" / "scalefactors_json.json").open() as handle:
        scale = float(json.load(handle)["tissue_lowres_scalef"])

    ax.imshow(image, zorder=0)
    ax.imshow(np.ones_like(image), alpha=0.48, zorder=1)
    x = df["pxl_col_in_fullres"].to_numpy(float) * scale
    y = df["pxl_row_in_fullres"].to_numpy(float) * scale
    values = df[column].to_numpy(float)
    vmin, vmax = np.percentile(values, [2, 98])
    order = np.argsort(values)
    scatter = ax.scatter(
        x[order],
        y[order],
        c=values[order],
        s=4.5,
        cmap=score_cmap(color),
        vmin=vmin,
        vmax=vmax,
        linewidths=0,
        alpha=0.94,
        zorder=2,
    )
    ax.set_title(title, loc="left", fontweight="bold", pad=2)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    colorbar = plt.colorbar(scatter, ax=ax, orientation="horizontal", fraction=0.045, pad=0.018)
    colorbar.outline.set_linewidth(0.35)
    colorbar.set_ticks([vmin, vmax])
    colorbar.set_ticklabels([f"{vmin:.1f}", f"{vmax:.1f}"])
    colorbar.ax.tick_params(length=1.5, pad=1)


def plot_score_distributions(ax: plt.Axes, df: pd.DataFrame, summary: dict[str, str]) -> None:
    label_panel(ax, "d", x=-0.070, y=1.14)
    data = [df[column].to_numpy(float) for column, _, _ in SCORES]
    parts = ax.violinplot(data, positions=np.arange(1, 4), widths=0.74, showextrema=False)
    for body, (_, _, color) in zip(parts["bodies"], SCORES, strict=True):
        body.set_facecolor(color)
        body.set_edgecolor("#FFFFFF")
        body.set_alpha(0.55)
        body.set_linewidth(0.4)

    bp = ax.boxplot(
        data,
        positions=np.arange(1, 4),
        widths=0.28,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#111827", "linewidth": 0.85},
        boxprops={"facecolor": "#FFFFFF", "edgecolor": "#374151", "linewidth": 0.65},
        whiskerprops={"color": "#374151", "linewidth": 0.6},
        capprops={"color": "#374151", "linewidth": 0.6},
    )
    _ = bp
    for idx, (column, _, color) in enumerate(SCORES, start=1):
        p90 = float(summary[f"p90_{column}"])
        p95 = float(summary[f"p95_{column}"])
        ax.scatter([idx - 0.11], [p90], marker="^", s=18, color=color, edgecolors="#FFFFFF", linewidths=0.35, zorder=4)
        ax.scatter([idx + 0.11], [p95], marker="D", s=15, color=color, edgecolors="#FFFFFF", linewidths=0.35, zorder=4)
        ax.text(idx, p95 + 0.42, f"p90 {p90:.2f}\np95 {p95:.2f}", ha="center", va="bottom", fontsize=5.8, color="#374151")

    ax.set_xticks(np.arange(1, 4))
    ax.set_xticklabels(["Mechanical", "Immune", "Scar"])
    ax.set_ylabel("Composite score")
    ax.set_title("Score dispersion across human STEMI spots", loc="left", fontweight="bold")
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.45)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.text(
        0.985,
        0.05,
        "^ p90   D p95",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.9,
        color="#4B5563",
    )


def separation_lookup(separation: pd.DataFrame, field: str) -> dict[frozenset[str], float]:
    return {
        frozenset((row["left_score"], row["right_score"])): float(row[field])
        for _, row in separation.iterrows()
    }


def plot_spearman_matrix(ax: plt.Axes, separation: pd.DataFrame) -> None:
    label_panel(ax, "e", x=-0.18, y=1.14)
    lookup = separation_lookup(separation, "spearman_r")
    columns = [column for column, _, _ in SCORES]
    labels = ["Mech.", "Immune", "Scar"]
    matrix = np.eye(3)
    for i, left in enumerate(columns):
        for j, right in enumerate(columns):
            if i == j:
                continue
            matrix[i, j] = lookup[frozenset((left, right))]
    image = ax.imshow(matrix, cmap="coolwarm", vmin=-1, vmax=1)
    for i in range(3):
        for j in range(3):
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=7.0, color="#111827")
    ax.set_xticks(np.arange(3))
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.set_yticks(np.arange(3))
    ax.set_yticklabels(labels)
    ax.set_title("Pairwise Spearman r", loc="left", fontweight="bold")
    for spine in ax.spines.values():
        spine.set_visible(False)
    colorbar = plt.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.outline.set_linewidth(0.35)
    colorbar.ax.tick_params(length=1.5, pad=1)


def plot_hotspot_overlap(ax: plt.Axes, separation: pd.DataFrame) -> None:
    label_panel(ax, "f", x=-0.045, y=1.06)
    labels = []
    values = []
    annotations = []
    colors = []
    pair_colors = {
        "Mechanical vs immune": "#8C7A5B",
        "Mechanical vs scar": "#5A877E",
        "Immune vs scar": "#67945C",
    }
    for _, row in separation.iterrows():
        label = PAIR_LABELS[frozenset((row["left_score"], row["right_score"]))]
        labels.append(label)
        values.append(float(row["top_decile_jaccard"]))
        annotations.append(f"{int(row['top_decile_overlap_spots'])}/{int(row['top_decile_union_spots'])}")
        colors.append(pair_colors[label])

    x = np.arange(len(labels))
    ax.bar(x, values, color=colors, width=0.62, edgecolor="#FFFFFF", linewidth=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=12, ha="right")
    ax.set_ylabel("Top-decile Jaccard")
    ax.set_ylim(0, max(values) * 1.45 if values else 1)
    ax.set_title("Hotspot overlap remains partial", loc="left", fontweight="bold")
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.45)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for xpos, value, annotation in zip(x, values, annotations, strict=True):
        ax.text(xpos, value + 0.014, annotation, ha="center", va="bottom", fontsize=6.0, color="#374151")
    ax.text(
        0.01,
        0.94,
        "labels = overlapping top-decile spots / union",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.9,
        color="#4B5563",
    )


def main() -> None:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(SPOT_TABLE, sep="\t")
    summary = read_summary()
    separation = pd.read_csv(SEPARATION_TABLE, sep="\t")
    n_spots = f"{int(float(summary['n_spots'])):,}"

    fig = plt.figure(figsize=(7.15, 9.10), constrained_layout=False)
    grid = fig.add_gridspec(
        nrows=4,
        ncols=3,
        height_ratios=[1.48, 2.25, 1.75, 1.18],
        width_ratios=[1, 1, 1],
        hspace=0.68,
        wspace=0.46,
    )

    ax_a = fig.add_subplot(grid[0, 0:2])
    plot_transfer_design(ax_a, n_spots)
    ax_b = fig.add_subplot(grid[0, 2])
    plot_gene_coverage(ax_b, summary)

    for col, (score_col, label, color) in enumerate(SCORES):
        ax = fig.add_subplot(grid[1, col])
        plot_spatial_map(ax, df, score_col, label, color, "c" if col == 0 else None)

    ax_d = fig.add_subplot(grid[2, 0:2])
    plot_score_distributions(ax_d, df, summary)
    ax_e = fig.add_subplot(grid[2, 2])
    plot_spearman_matrix(ax_e, separation)

    ax_f = fig.add_subplot(grid[3, 0:3])
    plot_hotspot_overlap(ax_f, separation)

    fig.suptitle(
        "Human STEMI transfer validates separable spatial state outputs",
        x=0.02,
        y=0.993,
        ha="left",
        fontsize=10.0,
        fontweight="bold",
        color="#111827",
    )
    fig.text(
        0.02,
        0.970,
        "GSM6613090 Visium; fixed human signatures; no model retraining",
        ha="left",
        va="top",
        fontsize=7.2,
        color="#4B5563",
    )
    fig.subplots_adjust(left=0.070, right=0.970, bottom=0.065, top=0.895)

    out_stem = OUT_DIR / "gse214611_human_stemi_signature_transfer"
    for suffix in [".png", ".pdf", ".svg", ".tiff"]:
        save_kwargs = {"bbox_inches": "tight", "pad_inches": 0.05}
        if suffix in {".png", ".tiff"}:
            save_kwargs["dpi"] = 600
        if suffix == ".tiff":
            save_kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(out_stem.with_suffix(suffix), **save_kwargs)
    plt.close(fig)
    print(f"[figure] wrote {out_stem.with_suffix('.png')}")
    print(f"[figure] wrote {out_stem.with_suffix('.pdf')}")
    print(f"[figure] wrote {out_stem.with_suffix('.svg')}")
    print(f"[figure] wrote {out_stem.with_suffix('.tiff')}")


if __name__ == "__main__":
    main()
