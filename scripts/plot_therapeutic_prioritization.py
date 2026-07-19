#!/usr/bin/env python3
"""Generate Figure 6: spatial-state therapeutic hypothesis prioritization."""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

from project_paths import project_root
from typing import Any


ROOT = project_root(__file__)
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.append(str(LOCAL_DEPS))

import matplotlib  # noqa: E402

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

from therapeutic_prioritization_utils import (  # noqa: E402
    extract_gene_rows,
    load_10x_h5,
    normalize_csc_log1p,
    subset_spots,
)


TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
VISIUM = ROOT / "data/raw/gse214611/visium/GSM6613090_V_Human_STEMI"
RANKING = TABLES / "gse214611_therapeutic_candidate_ranking.tsv"
PAIRS = TABLES / "gse214611_therapeutic_drug_target_pairs.tsv"
HUMAN_SPOTS = TABLES / "gse214611_human_stemi_signature_scores_by_spot.tsv"
HUMAN_H5 = VISIUM / "filtered_feature_bc_matrix.h5"
OUT_BASE = FIGURES / "Figure_6_spatial_therapeutic_prioritization"
SOURCE_DATA = TABLES / "Source_Data_Figure_6_therapeutic_prioritization.tsv"
SPATIAL_SOURCE_DATA = TABLES / "Source_Data_Figure_6_representative_human_expression.tsv"

BLUE = "#2F6FB5"
ORANGE = "#D99032"
GREEN = "#3A9B65"
RED = "#C64F48"
PURPLE = "#6B5FB5"
INK = "#18212F"
MUTED = "#667085"
LIGHT = "#E5EAF0"
STATE_COLORS = {
    "mechanical_adverse": BLUE,
    "immune_fibrotic_adverse": ORANGE,
    "scar_repair_associated": GREEN,
}
TIER_COLORS = {
    "tier_a": "#007C5A",
    "tier_b": BLUE,
    "tier_c": "#98A2B3",
    "caution_protective": RED,
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 600,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 6.2,
            "axes.titlesize": 7.2,
            "axes.labelsize": 6.5,
            "xtick.labelsize": 5.8,
            "ytick.labelsize": 5.8,
            "axes.linewidth": 0.55,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def panel_label(ax: plt.Axes, label: str, x: float = -0.08, y: float = 1.05) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        fontweight="bold",
        color=INK,
    )


def draw_box(
    ax: plt.Axes,
    x: float,
    title: str,
    body: str,
    color: str,
    width: float = 0.205,
) -> None:
    patch = FancyBboxPatch(
        (x, 0.28),
        width,
        0.47,
        boxstyle="round,pad=0.012,rounding_size=0.025",
        linewidth=0.9,
        edgecolor=color,
        facecolor="#FFFFFF",
    )
    ax.add_patch(patch)
    ax.text(x + 0.014, 0.66, title, ha="left", va="top", fontsize=6.8, fontweight="bold", color=INK)
    ax.text(x + 0.014, 0.55, body, ha="left", va="top", fontsize=5.7, color=MUTED, linespacing=1.18)


def plot_workflow(ax: plt.Axes) -> None:
    ax.set_axis_off()
    panel_label(ax, "a", x=-0.015, y=1.02)
    ax.text(
        0.02,
        0.98,
        "State- and direction-aware therapeutic hypothesis prioritization",
        ha="left",
        va="top",
        fontsize=7.8,
        fontweight="bold",
        color=INK,
    )
    boxes = [
        (0.02, "Mouse spatial evidence", "6 sections\nreplicate consistency\nboundary gradients", BLUE),
        (0.275, "Human feasibility", "1 STEMI section\ndetection + spatial\nlocalization", GREEN),
        (0.53, "Pharmacology evidence", "Open Targets\nChEMBL + DGIdb\naction direction", PURPLE),
        (0.785, "Guarded prioritization", "repair/safety penalty\nTier A/B/C\nprotective caution", RED),
    ]
    for x, title, body, color in boxes:
        draw_box(ax, x, title, body, color)
    for left, right in [(0.225, 0.275), (0.48, 0.53), (0.735, 0.785)]:
        ax.add_patch(
            FancyArrowPatch(
                (left, 0.515),
                (right, 0.515),
                arrowstyle="-|>",
                mutation_scale=8,
                linewidth=1.0,
                color="#7B8493",
            )
        )
    ax.text(
        0.02,
        0.12,
        "Output: interpretable, time-aware hypotheses; not efficacy prediction or clinical guidance",
        ha="left",
        va="center",
        fontsize=5.7,
        color=MUTED,
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def select_heatmap_rows(ranking: pd.DataFrame) -> pd.DataFrame:
    symbols = [
        "CD74",
        "COL6A1",
        "ENG",
        "LOXL2",
        "IL1B",
        "TGFB1",
        "NPPA",
        "MYH7",
        "POSTN",
        "CTHRC1",
        "LOX",
    ]
    indexed = ranking.set_index("human_gene_name")
    rows = indexed.loc[[symbol for symbol in symbols if symbol in indexed.index]].copy()
    return rows


def plot_evidence_heatmap(ax: plt.Axes, ranking: pd.DataFrame) -> None:
    panel_label(ax, "b", x=-0.15, y=1.06)
    rows = select_heatmap_rows(ranking)
    raw = np.column_stack(
        [
            rows["mouse_spatial_score"].to_numpy(float),
            rows["human_support_score"].to_numpy(float),
            rows["pharmacology_score"].to_numpy(float),
            rows["repair_safety_penalty"].to_numpy(float),
        ]
    )
    normalized = raw / np.asarray([40.0, 20.0, 30.0, 20.0])
    normalized[:, 3] *= -1
    cmap = LinearSegmentedColormap.from_list(
        "opportunity_risk", [RED, "#F6E6E4", "#FFFFFF", "#DCEAF6", BLUE]
    )
    ax.imshow(normalized, cmap=cmap, norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1), aspect="auto")
    ax.set_title("Evidence components and preservation penalty", loc="left", fontweight="bold", pad=6)
    ax.set_xticks(range(4), ["Mouse\n/40", "Human\n/20", "Pharm.\n/30", "Penalty\n/20"])
    ax.set_yticks(range(len(rows)), rows.index)
    ax.tick_params(length=0, pad=2)
    for idx, symbol in enumerate(rows.index):
        ax.get_yticklabels()[idx].set_color(TIER_COLORS[rows.loc[symbol, "priority_tier"]])
        ax.get_yticklabels()[idx].set_fontweight("bold" if rows.loc[symbol, "priority_tier"] == "tier_a" else "normal")
    for i in range(raw.shape[0]):
        for j in range(raw.shape[1]):
            color = "white" if abs(normalized[i, j]) > 0.62 else INK
            ax.text(j, i, f"{raw[i, j]:.0f}", ha="center", va="center", fontsize=5.2, color=color)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(
        0.0,
        -0.12,
        "Blue = opportunity evidence; red = repair/structural/safety deduction",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.3,
        color=MUTED,
    )


def load_representative_expression(symbols: list[str]) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    spots = pd.read_csv(HUMAN_SPOTS, sep="\t")
    raw = load_10x_h5(HUMAN_H5)
    barcode_index = {barcode: idx for idx, barcode in enumerate(raw.barcodes)}
    spots = spots.loc[spots["barcode"].isin(barcode_index)].copy()
    indices = [barcode_index[barcode] for barcode in spots["barcode"]]
    matrix = normalize_csc_log1p(subset_spots(raw, indices))
    gene_index = {name.upper(): idx for idx, name in enumerate(matrix.gene_names)}
    selected = [gene_index[symbol] for symbol in symbols]
    dense = extract_gene_rows(matrix, selected)
    return spots, {symbol: dense[idx] for idx, symbol in enumerate(symbols)}


def plot_spatial_target(
    ax: plt.Axes,
    spots: pd.DataFrame,
    values: np.ndarray,
    symbol: str,
    color: str,
    image: np.ndarray,
    scale: float,
) -> None:
    cmap = LinearSegmentedColormap.from_list(f"{symbol}_map", ["#F5F7F9", color])
    ax.imshow(image, zorder=0)
    ax.imshow(np.ones_like(image), alpha=0.56, zorder=1)
    x = spots["pxl_col_in_fullres"].to_numpy(float) * scale
    y = spots["pxl_row_in_fullres"].to_numpy(float) * scale
    positive = values[values > 0]
    vmax = float(np.percentile(positive, 98)) if len(positive) else 1.0
    order = np.argsort(values)
    scatter = ax.scatter(
        x[order],
        y[order],
        c=values[order],
        s=3.8,
        cmap=cmap,
        vmin=0,
        vmax=max(vmax, 1e-6),
        linewidths=0,
        alpha=0.95,
        zorder=2,
    )
    ax.set_title(symbol, loc="left", fontweight="bold", pad=2, color=color)
    ax.set_axis_off()
    cbar = ax.figure.colorbar(scatter, ax=ax, orientation="horizontal", fraction=0.055, pad=0.015, aspect=16)
    cbar.set_ticks([0, max(vmax, 1e-6)])
    cbar.set_ticklabels(["0", f"{vmax:.1f}"])
    cbar.ax.tick_params(labelsize=4.5, length=1, pad=1)
    cbar.outline.set_linewidth(0.35)


def plot_spatial_maps(container: plt.Axes, fig: plt.Figure) -> pd.DataFrame:
    container.set_axis_off()
    panel_label(container, "c", x=-0.05, y=1.06)
    container.set_title("Representative human STEMI spatial localization", loc="left", fontweight="bold", pad=6)
    symbols = ["NPPA", "CD74", "CTHRC1"]
    spots, expression = load_representative_expression(symbols)
    image = plt.imread(VISIUM / "spatial" / "tissue_lowres_image.png")
    with (VISIUM / "spatial" / "scalefactors_json.json").open() as handle:
        scale = float(json.load(handle)["tissue_lowres_scalef"])
    subgrid = container.get_subplotspec().subgridspec(1, 3, wspace=0.07)
    axes = [fig.add_subplot(subgrid[0, idx]) for idx in range(3)]
    colors = [BLUE, ORANGE, GREEN]
    for ax, symbol, color in zip(axes, symbols, colors, strict=True):
        plot_spatial_target(ax, spots, expression[symbol], symbol, color, image, scale)
    source = spots[["barcode", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]].copy()
    for symbol in symbols:
        source[f"{symbol}_log_normalized_expression"] = expression[symbol]
    return source


def draw_network_node(
    ax: plt.Axes,
    x: float,
    y: float,
    text: str,
    color: str,
    width: float,
    face: str = "#FFFFFF",
) -> None:
    patch = FancyBboxPatch(
        (x - width / 2, y - 0.052),
        width,
        0.104,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=0.75,
        edgecolor=color,
        facecolor=face,
    )
    ax.add_patch(patch)
    ax.text(x, y, text, ha="center", va="center", fontsize=5.3, color=INK)


def network_arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str,
    dashed: bool = False,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=6,
            linewidth=0.8,
            linestyle="--" if dashed else "-",
            color=color,
            alpha=0.9,
        )
    )


def plot_network(ax: plt.Axes, pairs: pd.DataFrame) -> None:
    ax.set_axis_off()
    panel_label(ax, "d", x=-0.07, y=1.04)
    ax.set_title("Direction-resolved target-drug hypotheses", loc="left", fontweight="bold", pad=6)
    state_y = {"Mechanical-border": 0.78, "Immune-fibrotic": 0.50, "Scar-repair": 0.20}
    for label, y, color in [
        ("Mechanical-border", state_y["Mechanical-border"], BLUE),
        ("Immune-fibrotic", state_y["Immune-fibrotic"], ORANGE),
        ("Scar-repair", state_y["Scar-repair"], GREEN),
    ]:
        draw_network_node(ax, 0.13, y, label, color, 0.22, face="#F8FAFC")
    links = [
        ("Mechanical-border", "MYH7", "MAVACAMTEN", RED, True),
        ("Immune-fibrotic", "CD74", "MILATUZUMAB", "#007C5A", False),
        ("Immune-fibrotic", "IL1B", "CANAKINUMAB", BLUE, False),
        ("Immune-fibrotic", "ENG", "CAROTUXIMAB", "#B77A1F", True),
        ("Immune-fibrotic", "TGFB1", "FRESOLIMUMAB", RED, True),
        ("Scar-repair", "LOX", "CCT365623", RED, True),
    ]
    display_drugs = {
        "MAVACAMTEN": "Mavacamten",
        "MILATUZUMAB": "Milatuzumab",
        "CANAKINUMAB": "Canakinumab",
        "CAROTUXIMAB": "Carotuximab",
        "FRESOLIMUMAB": "Fresolimumab",
        "CCT365623": "CCT365623*",
    }
    target_positions = {
        "MYH7": 0.80,
        "CD74": 0.63,
        "IL1B": 0.50,
        "ENG": 0.37,
        "TGFB1": 0.24,
        "LOX": 0.11,
    }
    for state, target, drug, color, caution in links:
        y = target_positions[target]
        draw_network_node(ax, 0.50, y, target, color, 0.15, face="#FFFFFF")
        draw_network_node(ax, 0.84, y, display_drugs[drug], color, 0.25, face="#FFFFFF")
        network_arrow(ax, (0.24, state_y[state]), (0.42, y), STATE_COLORS.get({"Mechanical-border":"mechanical_adverse","Immune-fibrotic":"immune_fibrotic_adverse","Scar-repair":"scar_repair_associated"}[state], MUTED))
        network_arrow(ax, (0.58, y), (0.705, y), color, caution)
    ax.text(0.46, 0.95, "Target", ha="center", va="center", fontsize=5.2, color=MUTED)
    ax.text(0.82, 0.95, "Drug or ligand", ha="center", va="center", fontsize=5.2, color=MUTED)
    ax.plot([0.04, 0.12], [0.01, 0.01], color="#007C5A", lw=1.2)
    ax.text(0.13, 0.01, "direction matched", va="center", fontsize=4.8, color=MUTED)
    ax.plot([0.43, 0.51], [0.01, 0.01], color=RED, lw=1.2, ls="--")
    ax.text(0.52, 0.01, "timing/repair or direction caution", va="center", fontsize=4.8, color=MUTED)
    ax.text(0.52, -0.028, "* one-source exploratory interaction", va="center", fontsize=4.35, color=MUTED)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.04, 1)


def plot_opportunity_risk(ax: plt.Axes, ranking: pd.DataFrame) -> None:
    panel_label(ax, "e", x=-0.11, y=1.05)
    opportunity = (
        ranking["mouse_spatial_score"].to_numpy(float)
        + ranking["human_support_score"].to_numpy(float)
        + ranking["pharmacology_score"].to_numpy(float)
    )
    penalty = ranking["repair_safety_penalty"].to_numpy(float)
    for tier in ["tier_c", "tier_b", "caution_protective", "tier_a"]:
        mask = ranking["priority_tier"].eq(tier).to_numpy()
        label = {
            "tier_a": "Tier A",
            "tier_b": "Tier B",
            "tier_c": "Tier C",
            "caution_protective": "Repair/structural caution",
        }[tier]
        ax.scatter(
            opportunity[mask],
            penalty[mask],
            s=19 if tier != "tier_a" else 36,
            color=TIER_COLORS[tier],
            edgecolors="#FFFFFF",
            linewidths=0.45,
            alpha=0.83,
            label=label,
            zorder=3,
        )
    xline = np.linspace(45, 90, 100)
    ax.plot(xline, xline - 65, color="#7B8493", lw=0.75, ls="--", zorder=1)
    ax.fill_between(xline, 0, np.clip(xline - 65, 0, 20), color="#DCEFE8", alpha=0.45, zorder=0)
    ax.text(84, 2.0, "Tier A rule\n(score >=65 after penalty)", ha="right", va="bottom", fontsize=4.8, color="#42705F")
    label_offsets = {
        "CD74": (1.2, 0.7),
        "COL6A1": (1.0, 0.7),
        "IL1B": (1.0, 0.8),
        "TGFB1": (1.0, -1.0),
        "ENG": (0.8, 0.6),
        "MYH7": (0.8, -1.0),
        "NPPA": (0.8, 0.6),
        "POSTN": (0.8, 0.6),
        "CTHRC1": (-3.0, -2.2),
        "LOX": (0.8, -2.0),
    }
    for _, row in ranking.iterrows():
        symbol = row["human_gene_name"]
        if symbol not in label_offsets:
            continue
        x = float(row["mouse_spatial_score"] + row["human_support_score"] + row["pharmacology_score"])
        y = float(row["repair_safety_penalty"])
        dx, dy = label_offsets[symbol]
        ax.annotate(
            symbol,
            (x, y),
            xytext=(x + dx, y + dy),
            textcoords="data",
            fontsize=4.8,
            color=TIER_COLORS[row["priority_tier"]],
            arrowprops={"arrowstyle": "-", "lw": 0.35, "color": "#98A2B3"},
        )
    ax.set_title("Therapeutic opportunity versus preservation risk", loc="left", fontweight="bold", pad=6)
    ax.set_xlabel("Pre-penalty evidence (mouse + human + pharmacology)")
    ax.set_ylabel("Repair/safety penalty")
    ax.set_xlim(35, 92)
    ax.set_ylim(-0.5, 20.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="both", color="#EDF0F3", linewidth=0.5, zorder=0)
    ax.legend(loc="upper left", frameon=False, fontsize=4.8, ncols=2, handletextpad=0.4, columnspacing=0.8)


def save_figure(fig: plt.Figure) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_BASE.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_BASE.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_BASE.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(
        OUT_BASE.with_suffix(".tiff"),
        dpi=600,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )


def main() -> None:
    setup_style()
    ranking = pd.read_csv(RANKING, sep="\t")
    pairs = pd.read_csv(PAIRS, sep="\t")
    fig = plt.figure(figsize=(7.205, 7.0), facecolor="white")
    grid = fig.add_gridspec(
        3,
        2,
        height_ratios=[0.78, 1.55, 1.65],
        width_ratios=[0.93, 1.25],
        hspace=0.37,
        wspace=0.34,
        left=0.075,
        right=0.985,
        top=0.975,
        bottom=0.065,
    )
    ax_a = fig.add_subplot(grid[0, :])
    ax_b = fig.add_subplot(grid[1, 0])
    ax_c = fig.add_subplot(grid[1, 1])
    ax_d = fig.add_subplot(grid[2, 0])
    ax_e = fig.add_subplot(grid[2, 1])
    plot_workflow(ax_a)
    plot_evidence_heatmap(ax_b, ranking)
    spatial_source = plot_spatial_maps(ax_c, fig)
    plot_network(ax_d, pairs)
    plot_opportunity_risk(ax_e, ranking)
    save_figure(fig)
    plt.close(fig)

    source_columns = [
        "human_gene_name",
        "primary_state",
        "mouse_spatial_score",
        "human_support_score",
        "pharmacology_score",
        "repair_safety_penalty",
        "final_priority_score",
        "priority_tier",
        "best_drug_name",
        "best_drug_action",
        "direction_match",
        "best_pair_source_count",
    ]
    ranking[source_columns].to_csv(SOURCE_DATA, sep="\t", index=False)
    spatial_source.to_csv(SPATIAL_SOURCE_DATA, sep="\t", index=False)
    print(f"Wrote {OUT_BASE}.png/.pdf/.svg/.tiff")
    print(f"Wrote {SOURCE_DATA} and {SPATIAL_SOURCE_DATA}")


if __name__ == "__main__":
    main()
