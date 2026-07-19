#!/usr/bin/env python3
"""Plot Extended Data Figure 7 for independence and human null checks.

Core conclusion: score-derived spatial states become more separable by day 7
without author-domain input, while human STEMI score maps show spatial
autocorrelation beyond expression-matched random signatures.
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

from project_paths import project_root


ROOT = project_root(__file__)
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


TABLE_DIR = ROOT / "results" / "tables"
FIGURE_DIR = ROOT / "results" / "figures"
RECOVERY_TABLE = TABLE_DIR / "gse214611_graph_smoothed_score_state_recovery.tsv"
STATE_TABLE = TABLE_DIR / "gse214611_graph_smoothed_score_states_by_spot.tsv"
HUMAN_SUMMARY = TABLE_DIR / "gse214611_human_expression_matched_spatial_null.tsv"
HUMAN_NULL = TABLE_DIR / "gse214611_human_expression_matched_spatial_null_distribution.tsv"
OUT_BASE = FIGURE_DIR / "cvr_additional_independence_checks"


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.linewidth"] = 0.7
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["legend.frameon"] = False

MECHANICAL = "#2C73B9"
IMMUNE = "#DF8B1D"
SCAR = "#3E9B68"
OVERLAP = "#7667A8"
STROMAL = "#879B53"
MIXED = "#C8CDD4"
D3_COLOR = "#4D7298"
D7_COLOR = "#C65D3B"
TEXT = "#20242A"
GRID = "#D8DDE5"

STATE_COLORS = {
    "mechanical-high": MECHANICAL,
    "scar-high": SCAR,
    "immune-high": IMMUNE,
    "mechanical/scar-high": OVERLAP,
    "mechanical/immune-high": "#8B78A8",
    "scar/immune-high": STROMAL,
    "mechanical/scar/immune-high": "#6D6A76",
    "mixed/low": MIXED,
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.08, 1.04, label, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")


def plot_state_map(ax: plt.Axes, rows: list[dict[str, str]], sample: str, title: str) -> None:
    sample_rows = [row for row in rows if row["sample"] == sample]
    ordered_states = [state for state in STATE_COLORS if any(row["graph_smoothed_state"] == state for row in sample_rows)]
    for state in ordered_states:
        subset = [row for row in sample_rows if row["graph_smoothed_state"] == state]
        ax.scatter(
            [float(row["pxl_col_in_fullres"]) for row in subset],
            [-float(row["pxl_row_in_fullres"]) for row in subset],
            s=4.0,
            c=STATE_COLORS[state],
            linewidths=0,
            alpha=0.94,
            rasterized=True,
            label=state,
        )
    ax.set_title(title, loc="left", fontsize=7.5, fontweight="bold", color=TEXT, pad=3)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def stage_means(rows: list[dict[str, str]], metric: str) -> dict[tuple[str, float], float]:
    grouped: dict[tuple[str, float], list[float]] = defaultdict(list)
    for row in rows:
        stage = "D3" if row["stage"] == "day3_mi" else "D7"
        grouped[(stage, float(row["graph_smoothing_alpha"]))].append(float(row[metric]))
    return {key: float(np.mean(values)) for key, values in grouped.items()}


def plot_metric_sensitivity(ax: plt.Axes, rows: list[dict[str, str]], metric: str, ylabel: str) -> None:
    by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_sample[row["sample"]].append(row)
    for sample, sample_rows in sorted(by_sample.items()):
        sample_rows.sort(key=lambda row: float(row["graph_smoothing_alpha"]))
        color = D3_COLOR if sample.startswith("D3") else D7_COLOR
        ax.plot(
            [float(row["graph_smoothing_alpha"]) for row in sample_rows],
            [float(row[metric]) for row in sample_rows],
            color=color,
            alpha=0.28,
            lw=0.8,
            marker="o",
            ms=2.2,
        )
    means = stage_means(rows, metric)
    alphas = sorted({float(row["graph_smoothing_alpha"]) for row in rows})
    for stage, color, marker in [("D3", D3_COLOR, "o"), ("D7", D7_COLOR, "s")]:
        mean_values = [means[(stage, alpha)] for alpha in alphas]
        ax.plot(
            alphas,
            mean_values,
            color=color,
            lw=1.8,
            marker=marker,
            ms=4,
        )
        ax.text(alphas[-1] + 0.018, mean_values[-1], stage, color=color, fontsize=6.2, va="center", fontweight="bold")
    ax.set_xlabel("Graph smoothing weight")
    ax.set_ylabel(ylabel)
    ax.set_xticks(alphas)
    ax.set_xlim(min(alphas) - 0.03, max(alphas) + 0.09)
    ax.grid(axis="y", color=GRID, linewidth=0.45)
    ax.set_axisbelow(True)


def plot_separability(ax: plt.Axes, rows: list[dict[str, str]]) -> None:
    primary = [row for row in rows if abs(float(row["graph_smoothing_alpha"]) - 0.35) < 1e-9]
    for x, stage, color in [(0, "day3_mi", D3_COLOR), (1, "day7_mi", D7_COLOR)]:
        values = [float(row["mechanical_scar_state_separable"]) for row in primary if row["stage"] == stage]
        jitter = np.linspace(-0.08, 0.08, len(values))
        ax.scatter(np.full(len(values), x) + jitter, values, s=28, color=color, edgecolor="white", linewidth=0.6, zorder=3)
        ax.text(x, np.mean(values) + (0.10 if np.mean(values) < 0.5 else -0.13), f"{int(sum(values))}/{len(values)}", ha="center", fontsize=7, fontweight="bold", color=color)
    ax.set_xlim(-0.45, 1.45)
    ax.set_ylim(-0.18, 1.18)
    ax.set_xticks([0, 1], ["D3", "D7"])
    ax.set_yticks([0, 1], ["overlap", "separate"])
    ax.set_ylabel("Mechanical vs scar state")
    ax.grid(axis="y", color=GRID, linewidth=0.45)
    ax.set_axisbelow(True)


def plot_human_null(ax: plt.Axes, summary: list[dict[str, str]], null_rows: list[dict[str, str]]) -> None:
    order = ["mechanical_border", "immune_fibrotic_activation", "fibroblast_scar_repair"]
    labels = ["Mechanical", "Immune", "Scar"]
    colors = [MECHANICAL, IMMUNE, SCAR]
    null_by_output: dict[str, list[float]] = defaultdict(list)
    for row in null_rows:
        null_by_output[row["output"]].append(float(row["graph_morans_i"]))
    summary_by_output = {row["output"]: row for row in summary}
    positions = np.arange(3)
    box = ax.boxplot(
        [null_by_output[name] for name in order],
        positions=positions,
        widths=0.56,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#4A4A4A", "linewidth": 0.9},
        whiskerprops={"color": "#8F98A3", "linewidth": 0.8},
        capprops={"color": "#8F98A3", "linewidth": 0.8},
    )
    for patch in box["boxes"]:
        patch.set_facecolor("#E6E9ED")
        patch.set_edgecolor("#9AA3AE")
        patch.set_linewidth(0.7)
    for x, name, color in zip(positions, order, colors, strict=True):
        row = summary_by_output[name]
        observed = float(row["observed_graph_morans_i"])
        p_value = float(row["moran_empirical_upper_p"])
        ax.scatter(x, observed, s=42, marker="D", color=color, edgecolor="white", linewidth=0.7, zorder=4)
        ax.text(x, observed + 0.035, f"P={p_value:.3f}", ha="center", va="bottom", fontsize=5.8, color=TEXT)
    ax.set_xticks(positions, labels, rotation=20, ha="right")
    ax.set_ylabel("Graph Moran's I")
    ax.set_ylim(-0.03, 0.75)
    ax.grid(axis="y", color=GRID, linewidth=0.45)
    ax.set_axisbelow(True)
    ax.text(0.02, 0.02, "Grey: 500 expression-matched random signatures\nDiamonds: observed fixed signatures", transform=ax.transAxes, fontsize=5.6, color="#5F6874", va="bottom")


def main() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    recovery = read_tsv(RECOVERY_TABLE)
    states = read_tsv(STATE_TABLE)
    human_summary = read_tsv(HUMAN_SUMMARY)
    human_null = read_tsv(HUMAN_NULL)

    fig = plt.figure(figsize=(7.20, 5.30), facecolor="white")
    grid = fig.add_gridspec(2, 6, height_ratios=[1.18, 0.90], hspace=0.43, wspace=0.92)
    ax_a = fig.add_subplot(grid[0, 0:3])
    ax_b = fig.add_subplot(grid[0, 3:6])
    ax_c = fig.add_subplot(grid[1, 0:2])
    ax_d = fig.add_subplot(grid[1, 2:4])
    ax_e = fig.add_subplot(grid[1, 4:6])

    plot_state_map(ax_a, states, "D3_3", "D3: score states remain partly intermingled")
    plot_state_map(ax_b, states, "D7_2", "D7: mechanical and stromal states separate")
    plot_separability(ax_c, recovery)
    plot_metric_sensitivity(ax_d, recovery, "same_state_edge_fraction", "Same-state graph edges")
    plot_human_null(ax_e, human_summary, human_null)

    for label, ax in zip("abcde", [ax_a, ax_b, ax_c, ax_d, ax_e], strict=True):
        add_panel_label(ax, label)

    handles = [
        mpl.lines.Line2D([], [], marker="o", linestyle="", color=color, markersize=4.5, label=label)
        for label, color in [
            ("Mechanical-high", MECHANICAL),
            ("Scar/immune-high", STROMAL),
            ("Mechanical/scar overlap", OVERLAP),
            ("Immune-high", IMMUNE),
            ("Mixed/low", MIXED),
        ]
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.995), ncol=5, fontsize=6.1, handletextpad=0.35, columnspacing=1.0)
    fig.suptitle(
        "Domain-independent graph states and human spatial-null checks",
        x=0.015,
        y=1.035,
        ha="left",
        fontsize=9.2,
        fontweight="bold",
        color=TEXT,
    )
    fig.text(
        0.015,
        1.005,
        "Author domain labels were excluded during state construction; human tests remain single-section internal falsification.",
        ha="left",
        fontsize=6.4,
        color="#5F6874",
    )

    fig.savefig(f"{OUT_BASE}.svg", bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.pdf", bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.png", dpi=350, bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"Wrote {OUT_BASE}.svg/.pdf/.png/.tiff")


if __name__ == "__main__":
    main()
