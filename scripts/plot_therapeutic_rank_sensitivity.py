#!/usr/bin/env python3
"""Generate Supplementary Figure S8 for therapeutic-rank robustness."""

from __future__ import annotations

import sys
from pathlib import Path

from project_paths import project_root


ROOT = project_root(__file__)
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.append(str(LOCAL_DEPS))

import matplotlib  # noqa: E402

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402


TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
RANKING = TABLES / "gse214611_therapeutic_candidate_ranking.tsv"
SENSITIVITY = TABLES / "gse214611_therapeutic_rank_sensitivity.tsv"
NULL = TABLES / "gse214611_therapeutic_spatial_score_null.tsv"
OUT_BASE = FIGURES / "Supplementary_Figure_S8_therapeutic_sensitivity"
SOURCE_DATA = TABLES / "Source_Data_Supplementary_Figure_S8_therapeutic_sensitivity.tsv"

BLUE = "#2F6FB5"
ORANGE = "#D99032"
GREEN = "#3A9B65"
RED = "#C64F48"
INK = "#18212F"
MUTED = "#667085"
STATE_COLORS = {
    "mechanical_adverse": BLUE,
    "immune_fibrotic_adverse": ORANGE,
    "scar_repair_associated": GREEN,
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 600,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 6.3,
            "axes.titlesize": 7.3,
            "axes.labelsize": 6.5,
            "xtick.labelsize": 5.7,
            "ytick.labelsize": 5.7,
            "axes.linewidth": 0.55,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.07,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        fontweight="bold",
        color=INK,
    )


def plot_loso_heatmap(ax: plt.Axes, ranking: pd.DataFrame, sensitivity: pd.DataFrame) -> None:
    panel_label(ax, "a")
    top = ranking.sort_values("overall_rank").head(12)["human_gene_name"].tolist()
    subset = sensitivity.loc[
        sensitivity["scenario"].str.startswith("leave_out_")
        & sensitivity["human_gene_name"].isin(top)
    ].copy()
    scenarios = sorted(subset["scenario"].unique())
    matrix = subset.pivot(index="human_gene_name", columns="scenario", values="scenario_rank").reindex(top)[scenarios]
    cmap = LinearSegmentedColormap.from_list("rank", ["#E7F2EE", "#F6E8C8", "#D3655A"])
    image = ax.imshow(matrix.to_numpy(float), aspect="auto", cmap=cmap, vmin=1, vmax=40)
    ax.set_title("Leave-one-section-out rank stability", loc="left", fontweight="bold", pad=6)
    ax.set_yticks(range(len(top)), top)
    labels = [scenario.replace("leave_out_", "") for scenario in scenarios]
    ax.set_xticks(range(len(labels)), labels, rotation=35, ha="right")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix.iat[i, j]
            ax.text(j, i, f"{value:.0f}", ha="center", va="center", fontsize=4.6, color=INK)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = ax.figure.colorbar(image, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("Rank", fontsize=5.4)
    cbar.ax.tick_params(labelsize=4.8, length=2)


def parse_range(text: str) -> tuple[int, int]:
    left, right = str(text).split("-", 1)
    return int(left), int(right)


def plot_weight_sensitivity(ax: plt.Axes, ranking: pd.DataFrame, sensitivity: pd.DataFrame) -> None:
    panel_label(ax, "b")
    primary = sensitivity.loc[sensitivity["scenario"].eq("primary")].set_index("human_gene_name")
    top = ranking.sort_values("overall_rank").head(15)["human_gene_name"].tolist()[::-1]
    y = np.arange(len(top))
    for idx, symbol in enumerate(top):
        row = primary.loc[symbol]
        low, high = parse_range(row["rank_range_all_scenarios"])
        median = float(row["median_rank_all_scenarios"])
        freq = float(row["top10_frequency_all_scenarios"])
        ax.plot([low, high], [idx, idx], color="#B9C1CB", lw=1.2, zorder=1)
        ax.scatter([median], [idx], s=18, color=BLUE, edgecolors="white", linewidths=0.35, zorder=2)
        ax.text(high + 1, idx, f"{freq:.0%}", va="center", fontsize=4.7, color=MUTED)
    ax.set_yticks(y, top)
    ax.set_xlabel("Rank across component omission, ±25% weights and LOSO")
    ax.set_title("Rank range and top-10 retention", loc="left", fontweight="bold", pad=6)
    ax.text(1.01, 1.01, "Top-10\nfrequency", transform=ax.transAxes, ha="left", va="bottom", fontsize=4.8, color=MUTED)
    ax.set_xlim(0, 65)
    ax.grid(axis="x", color="#EDF0F3", lw=0.5)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)


def plot_matched_null(ax: plt.Axes, null: pd.DataFrame) -> None:
    panel_label(ax, "c")
    symbols = ["CD74", "IL1B", "LOXL2", "TGFB1", "NPPA", "CTHRC1"]
    symbols = [symbol for symbol in symbols if symbol in set(null["candidate_gene"])]
    data = [
        null.loc[null["candidate_gene"].eq(symbol), "matched_mouse_spatial_score"].to_numpy(float)
        for symbol in symbols
    ]
    violins = ax.violinplot(data, positions=np.arange(len(symbols)), widths=0.75, showextrema=False)
    for body in violins["bodies"]:
        body.set_facecolor("#DCE6F0")
        body.set_edgecolor("#8DA5BE")
        body.set_alpha(0.9)
        body.set_linewidth(0.5)
    observed = [
        float(null.loc[null["candidate_gene"].eq(symbol), "candidate_mouse_spatial_score"].iloc[0])
        for symbol in symbols
    ]
    ax.scatter(np.arange(len(symbols)), observed, marker="D", s=24, color=RED, edgecolors="white", linewidths=0.4, label="Observed candidate", zorder=3)
    ax.set_xticks(np.arange(len(symbols)), symbols, rotation=30, ha="right")
    ax.set_ylabel("Mouse spatial evidence score")
    ax.set_title("Expression/detection-matched spatial-score nulls", loc="left", fontweight="bold", pad=6)
    ax.legend(loc="lower left", frameon=False, fontsize=5.0)
    ax.set_ylim(27, 40)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#EDF0F3", lw=0.5)


def plot_source_concordance(ax: plt.Axes, ranking: pd.DataFrame) -> None:
    panel_label(ax, "d")
    counts = (
        ranking.groupby(["primary_state", "pharmacology_source_count"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=[0, 1, 2, 3], fill_value=0)
    )
    x = np.arange(4)
    bottom = np.zeros(4)
    state_labels = {
        "mechanical_adverse": "Mechanical-border",
        "immune_fibrotic_adverse": "Immune-fibrotic",
        "scar_repair_associated": "Scar-repair",
    }
    for state in ["mechanical_adverse", "immune_fibrotic_adverse", "scar_repair_associated"]:
        values = counts.loc[state].to_numpy(float) if state in counts.index else np.zeros(4)
        ax.bar(
            x,
            values,
            bottom=bottom,
            width=0.68,
            color=STATE_COLORS[state],
            edgecolor="white",
            linewidth=0.4,
            label=state_labels[state],
        )
        bottom += values
    for idx, total in enumerate(bottom):
        ax.text(idx, total + 0.8, f"{int(total)}", ha="center", va="bottom", fontsize=5.2, color=INK)
    ax.set_xticks(x, ["0", "1", "2", "3"])
    ax.set_xlabel("Independent pharmacology databases per target")
    ax.set_ylabel("Candidate targets")
    ax.set_title("Pharmacology-source concordance", loc="left", fontweight="bold", pad=6)
    ax.legend(loc="upper right", frameon=False, fontsize=5.0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#EDF0F3", lw=0.5)


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
    sensitivity = pd.read_csv(SENSITIVITY, sep="\t")
    null = pd.read_csv(NULL, sep="\t")
    fig, axes = plt.subplots(2, 2, figsize=(7.205, 5.04), facecolor="white")
    plt.subplots_adjust(left=0.09, right=0.985, top=0.95, bottom=0.11, hspace=0.44, wspace=0.34)
    plot_loso_heatmap(axes[0, 0], ranking, sensitivity)
    plot_weight_sensitivity(axes[0, 1], ranking, sensitivity)
    plot_matched_null(axes[1, 0], null)
    plot_source_concordance(axes[1, 1], ranking)
    save_figure(fig)
    plt.close(fig)

    primary = sensitivity.loc[sensitivity["scenario"].eq("primary")].copy()
    source = ranking[
        [
            "human_gene_name",
            "primary_state",
            "priority_tier",
            "pharmacology_source_count",
            "overall_rank",
        ]
    ].merge(
        primary[
            [
                "human_gene_name",
                "median_rank_all_scenarios",
                "rank_range_all_scenarios",
                "top10_frequency_all_scenarios",
            ]
        ],
        on="human_gene_name",
        how="left",
    )
    source.to_csv(SOURCE_DATA, sep="\t", index=False)
    print(f"Wrote {OUT_BASE}.png/.pdf/.svg/.tiff")
    print(f"Wrote {SOURCE_DATA}")


if __name__ == "__main__":
    main()
