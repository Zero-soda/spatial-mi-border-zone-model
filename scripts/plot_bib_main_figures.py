#!/usr/bin/env python3
"""Build the five BIB-facing SSTBA figures from released source data."""

from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
OUTPUT = PROJECT / "02_main_figures"
EDITABLE = PROJECT / "06_qa" / "editable_sources"
SOURCE = ROOT / "source_data"
LEGACY_FIGURES = ROOT / "figures"

BLUE = "#2878B5"
ORANGE = "#DF8F27"
GREEN = "#2F9964"
PURPLE = "#5A50A5"
RED = "#C7473A"
INK = "#17212B"
MID = "#5A6878"
LIGHT = "#D9E0E8"
PALE = "#F5F7FA"

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.7,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
    }
)


def save_figure(fig: plt.Figure, stem: str) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    EDITABLE.mkdir(parents=True, exist_ok=True)
    fig.savefig(EDITABLE / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(EDITABLE / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(
        OUTPUT / f"{stem}.tiff",
        dpi=400,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    fig.savefig(EDITABLE / f"{stem}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.06,
        1.04,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        color=INK,
        va="bottom",
    )


def draw_box(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    title: str,
    body: str,
    colour: str,
) -> None:
    x, y = xy
    ax.add_patch(
        Rectangle(
            (x, y),
            width,
            height,
            facecolor="white",
            edgecolor=colour,
            linewidth=1.2,
        )
    )
    ax.text(
        x + 0.025,
        y + height - 0.06,
        title,
        fontsize=8.3,
        fontweight="bold",
        color=INK,
        va="top",
    )
    ax.text(
        x + 0.025,
        y + height - 0.14,
        body,
        fontsize=6.8,
        color=MID,
        va="top",
        linespacing=1.35,
    )


def figure_1() -> None:
    fig = plt.figure(figsize=(7.2, 5.35))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(
        0.035,
        0.955,
        "Figure 1 | SSTBA turns study-specific spatial atlases into auditable state and boundary outputs",
        fontsize=13,
        fontweight="bold",
        color=INK,
        va="top",
    )
    ax.text(
        0.035,
        0.905,
        "Frozen signatures, explicit graph choices, falsification diagnostics and transfer without retraining",
        fontsize=8.2,
        color=MID,
        va="top",
    )

    boxes = [
        (
            "1  Configure",
            "Expression or score table\ncoordinates + labels\nversioned signatures",
            BLUE,
        ),
        (
            "2  Score",
            "Per-gene z scores\nfrozen module means\ncoverage diagnostics",
            PURPLE,
        ),
        (
            "3  Build boundary",
            "Radius or k-nearest graph\nprespecified label pair\nedge-local gradients",
            ORANGE,
        ),
        (
            "4  Transfer",
            "Independent mouse + human\nno parameter fitting\nsame output contract",
            GREEN,
        ),
    ]
    xs = [0.04, 0.285, 0.53, 0.775]
    for x, (title, body, colour) in zip(xs, boxes, strict=True):
        draw_box(ax, (x, 0.62), 0.185, 0.205, title, body, colour)
    for left, right in zip(xs[:-1], xs[1:], strict=True):
        ax.add_patch(
            FancyArrowPatch(
                (left + 0.188, 0.722),
                (right - 0.008, 0.722),
                arrowstyle="-|>",
                mutation_scale=11,
                color=MID,
                linewidth=1.1,
            )
        )

    ax.text(
        0.04,
        0.555,
        "Three interpretable outputs",
        fontsize=9,
        fontweight="bold",
        color=INK,
    )
    output_rows = [
        ("Mechanical-border", "cardiomyocyte edge and border-zone direction", BLUE),
        ("Immune-fibrotic", "myeloid, TGF-β and fibro-inflammatory activity", ORANGE),
        ("Fibroblast-scar repair", "POSTN/CTHRC1-associated reparative scar state", GREEN),
    ]
    for index, (title, body, colour) in enumerate(output_rows):
        y = 0.49 - index * 0.085
        ax.add_patch(Rectangle((0.045, y), 0.012, 0.036, color=colour))
        ax.text(0.07, y + 0.026, title, fontsize=7.7, fontweight="bold", color=INK)
        ax.text(0.255, y + 0.026, body, fontsize=6.8, color=MID)

    ax.text(
        0.57,
        0.555,
        "Required diagnostics",
        fontsize=9,
        fontweight="bold",
        color=INK,
    )
    diagnostics = [
        "missing-gene and module coverage",
        "sample omission and module dropout",
        "matched random signatures and spatial nulls",
        "simple-score baselines and negative results",
    ]
    for index, item in enumerate(diagnostics):
        y = 0.505 - index * 0.062
        ax.text(0.575, y, f"{index + 1}", fontsize=7, fontweight="bold", color="white",
                ha="center", va="center",
                bbox={"boxstyle": "circle,pad=0.25", "facecolor": MID, "edgecolor": "none"})
        ax.text(0.602, y, item, fontsize=7.1, color=INK, va="center")

    ax.add_patch(Rectangle((0.04, 0.08), 0.92, 0.12, facecolor=PALE, edgecolor=LIGHT))
    ax.text(0.06, 0.162, "Reproducibility contract", fontsize=8.5, fontweight="bold", color=INK)
    ax.text(
        0.06,
        0.13,
        "CLI validation • deterministic test run • versions and parameters\n"
        "input/output SHA-256 hashes",
        fontsize=6.9,
        color=MID,
        va="center",
    )
    ax.text(
        0.92,
        0.158,
        "14,147 mouse spots\n+ independent mouse and human sections",
        fontsize=6.8,
        color=INK,
        ha="right",
        va="top",
    )
    save_figure(fig, "Figure_1")


def wrap_existing(source_name: str, stem: str, title: str) -> None:
    image = Image.open(LEGACY_FIGURES / source_name).convert("RGB")
    aspect = image.height / image.width
    fig = plt.figure(figsize=(7.2, 0.58 + 7.0 * aspect))
    ax = fig.add_axes([0.03, 0.03, 0.94, 0.89])
    ax.imshow(image)
    ax.axis("off")
    fig.text(0.03, 0.975, f"{stem.replace('_', ' ')} | {title}", fontsize=11.5,
             fontweight="bold", color=INK, va="top")
    save_figure(fig, stem)


def figure_2() -> None:
    spots = pd.read_csv(ROOT / "demo" / "mouse_master_spot_scores.tsv", sep="\t")
    score_columns = [
        ("mechanical_border_score", "Mechanical-border", "Blues"),
        ("immune_fibrotic_activation_score", "Immune-fibrotic", "Oranges"),
        ("fibroblast_scar_repair_score", "Fibroblast-scar repair", "Greens"),
    ]
    score_limits = {
        column: tuple(np.nanpercentile(spots[column], [2, 98]))
        for column, _, _ in score_columns
    }
    domain_colours = {
        "1": "#8FA8BD",
        "2": "#8ABF9C",
        "3": "#E3AA3B",
        "4": "#D45B43",
    }

    fig = plt.figure(figsize=(7.2, 7.15))
    grid = fig.add_gridspec(
        3,
        4,
        height_ratios=[1.0, 1.0, 0.78],
        hspace=0.38,
        wspace=0.12,
        left=0.055,
        right=0.985,
        top=0.91,
        bottom=0.09,
    )
    fig.suptitle(
        "Figure 2 | Discovery sections resolve three spatial-state outputs",
        x=0.04,
        y=0.98,
        ha="left",
        fontsize=11.5,
        fontweight="bold",
        color=INK,
    )

    representatives = [("D3_1", "Day 3 representative"), ("D7_1", "Day 7 representative")]
    for row, (sample, row_label) in enumerate(representatives):
        sample_spots = spots[spots["sample"] == sample]
        for column_index in range(4):
            ax = fig.add_subplot(grid[row, column_index])
            x = sample_spots["array_col"].to_numpy()
            y = sample_spots["array_row"].to_numpy()
            if column_index == 0:
                colours = sample_spots["author_domain"].astype(str).map(domain_colours)
                ax.scatter(x, y, c=colours, s=5.2, linewidths=0, rasterized=True)
                ax.set_title("Author domains", fontweight="bold")
                if row == 0:
                    handles = [
                        Line2D(
                            [0],
                            [0],
                            marker="o",
                            color="none",
                            markerfacecolor=domain_colours[label],
                            markeredgecolor="none",
                            label=f"Domain {label}",
                            markersize=4,
                        )
                        for label in ("1", "2", "3", "4")
                    ]
                    ax.legend(
                        handles=handles,
                        loc="lower center",
                        bbox_to_anchor=(0.5, -0.19),
                        ncol=2,
                        frameon=False,
                        fontsize=5.5,
                        handletextpad=0.2,
                        columnspacing=0.7,
                    )
            else:
                column, title, colour_map = score_columns[column_index - 1]
                low, high = score_limits[column]
                scatter = ax.scatter(
                    x,
                    y,
                    c=sample_spots[column],
                    cmap=colour_map,
                    vmin=low,
                    vmax=high,
                    s=5.2,
                    linewidths=0,
                    rasterized=True,
                )
                ax.set_title(title, fontweight="bold")
                colour_bar = fig.colorbar(
                    scatter,
                    ax=ax,
                    orientation="horizontal",
                    fraction=0.042,
                    pad=0.045,
                    aspect=24,
                )
                colour_bar.ax.tick_params(labelsize=5, length=1.5, pad=1)
                colour_bar.outline.set_linewidth(0.4)
            ax.set_aspect("equal")
            ax.invert_yaxis()
            ax.axis("off")
            if column_index == 0:
                ax.text(
                    -0.08,
                    0.5,
                    f"{row_label}\n{sample}",
                    transform=ax.transAxes,
                    ha="right",
                    va="center",
                    fontsize=7,
                    fontweight="bold",
                    color=INK,
                )

    summary = (
        spots[spots["author_domain"].astype(str).isin(["3", "4"])]
        .assign(author_domain=lambda frame: frame["author_domain"].astype(str))
        .groupby(["sample", "stage", "author_domain"], as_index=False)[
            [column for column, _, _ in score_columns]
        ]
        .mean()
    )
    stage_specs = [("day3_mi", [0, 1], "Day 3"), ("day7_mi", [3, 4], "Day 7")]
    for axis_index, (column, title, colour_map) in enumerate(
        (score_columns[0], score_columns[2])
    ):
        ax = fig.add_subplot(grid[2, axis_index * 2 : axis_index * 2 + 2])
        colour = BLUE if column == "mechanical_border_score" else GREEN
        for stage, positions, stage_label in stage_specs:
            stage_rows = summary[summary["stage"] == stage]
            for _, sample_rows in stage_rows.groupby("sample", sort=True):
                values = (
                    sample_rows.set_index("author_domain")
                    .reindex(["3", "4"])[column]
                    .to_numpy(float)
                )
                ax.plot(
                    positions,
                    values,
                    color=colour,
                    alpha=0.32,
                    linewidth=0.9,
                    marker="o",
                    markersize=2.7,
                )
            means = (
                stage_rows.groupby("author_domain")[column]
                .mean()
                .reindex(["3", "4"])
                .to_numpy(float)
            )
            ax.plot(
                positions,
                means,
                color=colour,
                linewidth=2.2,
                marker="o",
                markersize=4.5,
            )
            ax.text(
                np.mean(positions),
                ax.get_ylim()[1] if ax.get_ylim()[1] else 0,
                stage_label,
                ha="center",
                va="bottom",
                fontsize=6,
                color=MID,
            )
        ax.axhline(0, color=LIGHT, linewidth=0.7)
        ax.set_xticks(
            [0, 1, 3, 4],
            ["Domain 3", "Domain 4", "Domain 3", "Domain 4"],
        )
        ax.set_ylabel("Mean score per section")
        ax.set_title(
            f"All six sections: {title}",
            loc="left",
            fontweight="bold",
        )
        ax.text(
            0.5,
            -0.24,
            "thin: section  |  thick: stage mean",
            transform=ax.transAxes,
            ha="center",
            fontsize=5.8,
            color=MID,
        )
        panel_label(ax, "b" if axis_index == 0 else "c")
    first_map = fig.axes[0]
    panel_label(first_map, "a")
    save_figure(fig, "Figure_2")


def figure_3() -> None:
    graph = pd.read_csv(
        SOURCE / "Source_Data_Supplementary_Figure_S6_visium_graph_boundary_analysis.tsv",
        sep="\t",
    )
    gradients = pd.read_csv(
        SOURCE / "Source_Data_Supplementary_Figure_S6_visium_graph_boundary_edge_gradients.tsv",
        sep="\t",
    )
    boundary_image = Image.open(
        LEGACY_FIGURES / "Figure_4_domain34_boundary.png"
    ).convert("RGB")

    fig = plt.figure(figsize=(7.2, 7.0))
    grid = fig.add_gridspec(
        3,
        2,
        height_ratios=[2.7, 1.15, 1.15],
        width_ratios=[1, 1],
        hspace=0.55,
        wspace=0.42,
        left=0.08,
        right=0.98,
        top=0.92,
        bottom=0.08,
    )
    fig.suptitle(
        "Figure 3 | Six-nearest-neighbour analysis quantifies the domain 3-domain 4 interface",
        x=0.04,
        y=0.98,
        ha="left",
        fontsize=11.5,
        fontweight="bold",
        color=INK,
    )

    ax_a = fig.add_subplot(grid[0, :])
    ax_a.imshow(boundary_image)
    ax_a.axis("off")
    panel_label(ax_a, "a")

    ax_b = fig.add_subplot(grid[1, 0])
    x = np.arange(len(graph))
    ax_b.plot(x, graph["domain3_graph_boundary_fraction"], "o-", color=BLUE,
              label="Domain 3", linewidth=1.2, markersize=4)
    ax_b.plot(x, graph["domain4_graph_boundary_fraction"], "s-", color=GREEN,
              label="Domain 4", linewidth=1.2, markersize=4)
    ax_b.axvline(2.5, color=LIGHT, linewidth=1)
    ax_b.set_xticks(x, graph["sample"])
    ax_b.set_ylabel("Boundary fraction")
    ax_b.set_ylim(0, 0.85)
    ax_b.legend(loc="upper right", ncol=2)
    ax_b.set_title("Per-section interface fractions", loc="left", fontweight="bold")
    panel_label(ax_b, "b")

    ax_c = fig.add_subplot(grid[1, 1])
    edge_counts = graph["domain3_domain4_edge_count"]
    bars = ax_c.bar(x, edge_counts, color=[BLUE] * 3 + [GREEN] * 3, width=0.7)
    ax_c.set_xticks(x, graph["sample"])
    ax_c.set_ylabel("Domain 3-domain 4 edges")
    ax_c.set_title("Contact edges in the fixed graph", loc="left", fontweight="bold")
    panel_label(ax_c, "c")
    for bar, value in zip(bars, edge_counts, strict=True):
        ax_c.text(bar.get_x() + bar.get_width() / 2, value + 18, str(value),
                  ha="center", va="bottom", fontsize=6)

    ax_d = fig.add_subplot(grid[2, :])
    colour_map = {
        "mechanical_border": BLUE,
        "immune_fibrotic_activation": ORANGE,
        "fibroblast_scar_repair": GREEN,
    }
    label_map = {
        "mechanical_border": "Mechanical-border",
        "immune_fibrotic_activation": "Immune-fibrotic",
        "fibroblast_scar_repair": "Fibroblast-scar repair",
    }
    positions = np.arange(3)
    for offset, stage in zip((-0.13, 0.13), ("day3_mi", "day7_mi"), strict=True):
        stage_values = (
            gradients[gradients["stage"] == stage]
            .groupby("score")["mean_edge_delta"]
            .agg(["mean", "sem"])
            .reindex(colour_map)
        )
        for idx, score in enumerate(colour_map):
            ax_d.errorbar(
                idx + offset,
                stage_values.loc[score, "mean"],
                yerr=stage_values.loc[score, "sem"],
                fmt="o" if stage == "day3_mi" else "s",
                color=colour_map[score],
                markerfacecolor="white" if stage == "day3_mi" else colour_map[score],
                capsize=2,
                markersize=5,
                linewidth=1,
                label=("Day 3" if stage == "day3_mi" else "Day 7") if idx == 0 else None,
            )
    ax_d.axhline(0, color=MID, linewidth=0.7)
    ax_d.set_xticks(positions, [label_map[key] for key in colour_map])
    ax_d.set_ylabel("Mean edge delta: domain 4 − domain 3")
    ax_d.set_title(
        "Edge-local gradients retain mechanical direction and show a day 7 scar-repair direction",
        loc="left",
        fontweight="bold",
    )
    ax_d.legend(loc="upper left", ncol=2)
    panel_label(ax_d, "d")
    save_figure(fig, "Figure_3")


def short_label(value: str, width: int = 27) -> str:
    return "\n".join(textwrap.wrap(value.replace("_", " "), width=width))


def figure_4() -> None:
    loso = pd.read_csv(
        SOURCE / "Source_Data_Supplementary_Figure_S3_loso_robustness.tsv",
        sep="\t",
    )
    dropout = pd.read_csv(
        SOURCE / "Source_Data_Supplementary_Figure_S3_module_dropout_robustness.tsv",
        sep="\t",
    )
    random = pd.read_csv(
        SOURCE / "Source_Data_Supplementary_Figure_S3_random_signature_negative_control.tsv",
        sep="\t",
    )
    benchmark = pd.read_csv(
        SOURCE / "Source_Data_Figure_7_model_benchmark.tsv",
        sep="\t",
    )

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.4))
    fig.subplots_adjust(left=0.1, right=0.98, top=0.9, bottom=0.12, hspace=0.62, wspace=0.48)
    fig.suptitle(
        "Figure 4 | Protocol diagnostics expose stable directions and baseline-equivalent results",
        x=0.04,
        y=0.975,
        ha="left",
        fontsize=11.5,
        fontweight="bold",
        color=INK,
    )

    ax = axes[0, 0]
    omitted = loso[loso["omitted_sample"] != "none"]
    summary = omitted.groupby("check")["direction_preserved"].mean().sort_values()
    colours = [GREEN if value == 1 else ORANGE for value in summary]
    ax.barh(np.arange(len(summary)), summary, color=colours, height=0.65)
    ax.set_yticks(np.arange(len(summary)), [short_label(x, 24) for x in summary.index])
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Fraction of leave-one-section-out runs")
    ax.set_title("Direction preservation", loc="left", fontweight="bold")
    panel_label(ax, "a")

    ax = axes[0, 1]
    groups = []
    for (axis_name, stage), rows in dropout.groupby(["axis", "stage"], sort=False):
        baseline = float(rows.loc[rows["dropped_component"] == "none", "mean_margin"].iloc[0])
        removed = rows[rows["dropped_component"] != "none"]["mean_margin"]
        groups.append(
            {
                "axis": axis_name,
                "stage": stage,
                "baseline": baseline,
                "dropout_min": float(removed.min()),
                "dropout_max": float(removed.max()),
            }
        )
    plot_dropout = pd.DataFrame(groups)
    order = np.arange(len(plot_dropout))
    colours = plot_dropout["axis"].map(
        {
            "mechanical_border": BLUE,
            "immune_fibrotic_activation": ORANGE,
            "fibroblast_scar_repair": GREEN,
        }
    ).fillna(MID)
    ax.hlines(
        order,
        plot_dropout["dropout_min"],
        plot_dropout["dropout_max"],
        color=colours,
        linewidth=3,
        label="Range after one module is removed",
    )
    ax.scatter(
        plot_dropout["baseline"],
        order,
        c=colours,
        marker="D",
        edgecolor="white",
        linewidth=0.5,
        s=30,
        zorder=3,
        label="Complete score",
    )
    ax.axvline(0, color=MID, linewidth=0.7)
    axis_labels = {
        "mechanical_border": "Mechanical-border",
        "immune_fibrotic_activation": "Immune-fibrotic",
        "fibroblast_scar_repair": "Fibroblast-scar repair",
    }
    ax.set_yticks(
        order,
        [
            f"{axis_labels[row.axis]} | {row.stage.replace('_mi', '').replace('day', 'D')}"
            for row in plot_dropout.itertuples()
        ],
    )
    ax.invert_yaxis()
    ax.set_xlabel("Directional margin")
    ax.set_title("Module-dropout sensitivity", loc="left", fontweight="bold")
    ax.legend(loc="lower right", fontsize=5.7)
    panel_label(ax, "b")

    ax = axes[1, 0]
    y = np.arange(len(random))
    random_labels = [
        "Mechanical D3",
        "Mechanical D7",
        "Scar repair D7",
        "Immune-fibrotic D7",
    ]
    ax.hlines(y, random["random_p95_margin"], random["observed_mean_margin"],
              color=LIGHT, linewidth=2)
    ax.scatter(random["random_p95_margin"], y, marker="|", color=MID, s=80,
               label="Matched-random p95")
    status_colours = [GREEN if flag else RED for flag in random["observed_exceeds_random_p95"]]
    ax.scatter(random["observed_mean_margin"], y, color=status_colours, s=25,
               label="Observed")
    ax.set_yticks(y, random_labels)
    ax.invert_yaxis()
    ax.set_xlabel("Mean directional margin")
    ax.set_title("Matched-random signature controls", loc="left", fontweight="bold")
    ax.legend(loc="lower right")
    panel_label(ax, "c")

    ax = axes[1, 1]
    comparisons = list(dict.fromkeys(benchmark["comparison"]))
    y = np.arange(len(comparisons))
    for idx, comparison in enumerate(comparisons):
        rows = benchmark[benchmark["comparison"] == comparison]
        framework = float(
            rows.loc[rows["model"] == "three_output_framework", "rank_biserial_effect"].iloc[0]
        )
        baseline = float(
            rows.loc[
                rows["model"] == "simple_published_or_generic_baseline",
                "rank_biserial_effect",
            ].iloc[0]
        )
        ax.plot([baseline, framework], [idx, idx], color=LIGHT, linewidth=1.5)
        ax.scatter(baseline, idx, marker="s", color=MID, s=20)
        ax.scatter(framework, idx, marker="o", color=PURPLE, s=24)
    ax.axvline(0, color=MID, linewidth=0.7)
    ax.set_yticks(y, [short_label(x, 27) for x in comparisons])
    ax.invert_yaxis()
    ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel("Rank-biserial effect")
    ax.set_title("Transparent comparison with simple baselines", loc="left", fontweight="bold")
    ax.scatter([], [], marker="o", color=PURPLE, label="SSTBA")
    ax.scatter([], [], marker="s", color=MID, label="Simple baseline")
    ax.legend(loc="lower right")
    panel_label(ax, "d")
    save_figure(fig, "Figure_4")


def graphical_abstract() -> None:
    fig = plt.figure(figsize=(10.5, 5.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.05, 0.92, "SSTBA", fontsize=24, fontweight="bold", color=INK)
    ax.text(
        0.05,
        0.84,
        "Frozen spatial states + explicit boundaries + cross-study transfer",
        fontsize=12,
        color=MID,
    )
    flow = [
        ("Public spatial data", "mouse discovery\nindependent mouse + human", BLUE),
        ("Frozen state scores", "mechanical • immune-fibrotic\nfibroblast-scar repair", PURPLE),
        ("Boundary diagnostics", "six-nearest graph\nedge-local gradients + nulls", ORANGE),
        ("Auditable transfer", "no retraining\nmanifest + hashes + limitations", GREEN),
    ]
    xs = [0.05, 0.285, 0.52, 0.755]
    for x, (title, body, colour) in zip(xs, flow, strict=True):
        draw_box(ax, (x, 0.43), 0.19, 0.25, title, body, colour)
    for left, right in zip(xs[:-1], xs[1:], strict=True):
        ax.add_patch(
            FancyArrowPatch(
                (left + 0.193, 0.555),
                (right - 0.01, 0.555),
                arrowstyle="-|>",
                mutation_scale=13,
                color=MID,
                linewidth=1.3,
            )
        )
    ax.add_patch(Rectangle((0.05, 0.12), 0.895, 0.17, facecolor=PALE, edgecolor=LIGHT))
    ax.text(0.075, 0.235, "Real-data result", fontsize=10, fontweight="bold", color=INK)
    ax.text(
        0.075,
        0.17,
        "Mechanical-border and fibroblast-scar directions transferred across public MI datasets; "
        "simple fibrosis baselines remained equivalent for several scar comparisons.",
        fontsize=9,
        color=MID,
    )
    save_figure(fig, "graphical_abstract")


def main() -> None:
    figure_1()
    figure_2()
    figure_3()
    figure_4()
    wrap_existing(
        "Figure_7_cross_study_external_validation.png",
        "Figure_5",
        "Frozen outputs transfer to independent mouse and region-labelled human sections",
    )
    graphical_abstract()


if __name__ == "__main__":
    main()
