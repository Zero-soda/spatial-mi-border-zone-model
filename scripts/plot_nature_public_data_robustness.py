#!/usr/bin/env python3
"""Create a publication-style robustness figure for the MI border-zone model."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


TABLE_DIR = ROOT / "results" / "tables"
FIGURE_DIR = ROOT / "results" / "figures"

LOSO = TABLE_DIR / "gse214611_loso_stage_domain_robustness.tsv"
DROPOUT = TABLE_DIR / "gse214611_module_dropout_robustness.tsv"
BOUNDARY = TABLE_DIR / "gse214611_boundary_direction_robustness.tsv"
THRESHOLD = TABLE_DIR / "gse214611_boundary_threshold_sensitivity.tsv"
RANDOM_CONTROL = TABLE_DIR / "gse214611_random_signature_negative_control.tsv"
HUMAN = TABLE_DIR / "gse214611_human_stemi_score_separation.tsv"


PALETTE = {
    "blue": "#4C78A8",
    "teal": "#72B7B2",
    "orange": "#F58518",
    "purple": "#8E6C8A",
    "gray": "#6F7378",
    "light_gray": "#E6E8EB",
    "ink": "#202124",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def safe_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 7,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "legend.frameon": False,
        }
    )


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.16,
        1.14,
        label,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        va="top",
        ha="left",
    )


def plot_check_counts(
    ax: plt.Axes,
    loso: list[dict[str, str]],
    dropout: list[dict[str, str]],
    boundary: list[dict[str, str]],
    threshold: list[dict[str, str]],
    random_control: list[dict[str, str]],
) -> None:
    threshold_preserved = sum(
        row["mechanical_direction_preserved_all_samples"] == "True"
        and row["scar_direction_preserved_all_samples"] == "True"
        for row in threshold
    )
    random_preserved = sum(row["observed_exceeds_random_p95"] == "True" for row in random_control)
    groups = [
        ("LOSO", sum(row["direction_preserved"] == "True" for row in loso), len(loso), PALETTE["blue"]),
        ("Module dropout", sum(row["direction_preserved"] == "True" for row in dropout), len(dropout), PALETTE["teal"]),
        ("Boundary direction", sum(row["all_samples_preserved"] == "True" for row in boundary), len(boundary), PALETTE["orange"]),
        ("Threshold", threshold_preserved, len(threshold), PALETTE["purple"]),
        ("Random control", random_preserved, len(random_control), PALETTE["gray"]),
    ]
    values = [item[1] for item in groups]
    totals = [item[2] for item in groups]
    colors = [item[3] for item in groups]
    y = np.arange(len(groups))
    ax.barh(y, values, color=colors, height=0.50)
    for idx, (value, total) in enumerate(zip(values, totals, strict=True)):
        ax.text(value + max(totals) * 0.03, idx, f"{value}/{total}", ha="left", va="center", fontsize=6.5)
    ax.set_yticks(y, [item[0] for item in groups])
    ax.invert_yaxis()
    ax.set_xlabel("preserved checks")
    ax.set_xlim(0, max(totals) * 1.22)
    ax.set_title("Public-data robustness")
    panel_label(ax, "a")


def plot_module_dropout(ax: plt.Axes, rows: list[dict[str, str]]) -> None:
    filtered = [row for row in rows if row["dropped_component"] != "none"]
    label_map = {
        "mechanical_border": "mech",
        "fibroblast_scar_repair": "scar",
        "immune_fibrotic_activation": "immune",
    }
    colors = {
        "mechanical_border": PALETTE["blue"],
        "fibroblast_scar_repair": PALETTE["orange"],
        "immune_fibrotic_activation": PALETTE["teal"],
    }
    y = np.arange(len(filtered))
    margins = [safe_float(row["min_margin"]) for row in filtered]
    ax.axvline(0, color=PALETTE["gray"], linewidth=0.7)
    ax.scatter(margins, y, s=24, color=[colors[row["axis"]] for row in filtered], zorder=3)
    for xval, ypos, row in zip(margins, y, filtered, strict=True):
        ax.plot([0, xval], [ypos, ypos], color=colors[row["axis"]], linewidth=1.4, alpha=0.55)
    short = {
        "z_CM_BZ1_TRANSITION": "BZ1",
        "z_CM_BZ2_MECHANICAL_EDGE": "BZ2",
        "z_FAP_POSTN_PATHO_FIBROBLAST": "FAP/POSTN",
        "z_ECM_REMODELING": "ECM",
        "z_CTHRC1_REPARATIVE_CF": "CTHRC1",
        "z_MYOFIBROBLAST_CONTRACTILE": "MyoFB",
        "z_CCR2_IL1B_MYLOID": "CCR2/IL1B",
        "z_TGFB_SIGNALING": "TGFb",
    }
    stage_short = {"day3_mi": "D3", "day7_mi": "D7"}
    labels = [
        f"{stage_short.get(row['stage'], row['stage'])} {label_map[row['axis']]} drop {short.get(row['dropped_component'], row['dropped_component'])}"
        for row in filtered
    ]
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("minimum preserved margin")
    ax.set_title("No single module drives the axis call")
    panel_label(ax, "b")


def plot_boundary(ax: plt.Axes, rows: list[dict[str, str]]) -> None:
    items = []
    for row in rows:
        metric = row["metric"]
        if "mechanical" in metric:
            score = "mechanical"
            color = PALETTE["blue"]
        else:
            score = "scar repair"
            color = PALETTE["orange"]
        stage = "D3" if row["stage"] == "day3_mi" else "D7"
        items.append((f"{stage} {score}", safe_float(row["mean_value"]), safe_float(row["min_value"]), safe_float(row["max_value"]), color))
    y = np.arange(len(items))
    means = [item[1] for item in items]
    lows = [item[1] - item[2] for item in items]
    highs = [item[3] - item[1] for item in items]
    ax.axvline(0, color=PALETTE["gray"], linewidth=0.7)
    ax.errorbar(
        means,
        y,
        xerr=[lows, highs],
        fmt="o",
        color=PALETTE["ink"],
        ecolor=PALETTE["gray"],
        elinewidth=0.8,
        capsize=2,
    )
    for idx, item in enumerate(items):
        ax.scatter(item[1], idx, s=34, color=item[4], zorder=4)
    ax.set_yticks(y, [item[0] for item in items])
    ax.invert_yaxis()
    ax.set_xlabel("D4 - D3 contact gradient")
    ax.set_title("Boundary gradients\nchange with repair stage")
    panel_label(ax, "c")


def plot_human_heatmaps(ax_corr: plt.Axes, ax_jaccard: plt.Axes, rows: list[dict[str, str]]) -> None:
    labels = ["mechanical", "immune-fibrotic", "scar-repair"]
    key = {
        "mechanical_border_score": "mechanical",
        "immune_fibrotic_activation_score": "immune-fibrotic",
        "fibroblast_scar_repair_score": "scar-repair",
    }
    corr = np.eye(3)
    jac = np.eye(3)
    index = {label: idx for idx, label in enumerate(labels)}
    for row in rows:
        i = index[key[row["left_score"]]]
        j = index[key[row["right_score"]]]
        corr[i, j] = corr[j, i] = safe_float(row["spearman_r"])
        jac[i, j] = jac[j, i] = safe_float(row["top_decile_jaccard"])

    for ax, matrix, title, vmin, vmax, cmap in [
        (ax_corr, corr, "Spearman r", -1, 1, "coolwarm"),
        (ax_jaccard, jac, "Top-decile Jaccard", 0, 1, "YlGnBu"),
    ]:
        ax.imshow(matrix, vmin=vmin, vmax=vmax, cmap=cmap)
        ax.set_xticks(np.arange(3), ["mech", "immune", "scar"], rotation=0)
        ax.set_yticks(np.arange(3), labels)
        ax.set_title(title)
        for i in range(3):
            for j in range(3):
                ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=6)
        ax.tick_params(length=0)
        ax.text(1.02, -0.18, f"range {vmin} to {vmax}", transform=ax.transAxes, fontsize=5, ha="right")
    panel_label(ax_corr, "d")


def main() -> None:
    configure_matplotlib()
    loso = read_tsv(LOSO)
    dropout = read_tsv(DROPOUT)
    boundary = read_tsv(BOUNDARY)
    threshold = read_tsv(THRESHOLD)
    random_control = read_tsv(RANDOM_CONTROL)
    human = read_tsv(HUMAN)

    fig = plt.figure(figsize=(183 / 25.4, 155 / 25.4), constrained_layout=False)
    ax_a = fig.add_axes([0.09, 0.65, 0.26, 0.27])
    ax_c = fig.add_axes([0.09, 0.15, 0.26, 0.34])
    ax_b = fig.add_axes([0.48, 0.15, 0.24, 0.76])
    ax_d1 = fig.add_axes([0.80, 0.67, 0.16, 0.20])
    ax_d2 = fig.add_axes([0.80, 0.40, 0.16, 0.20])
    ax_note = fig.add_axes([0.77, 0.15, 0.21, 0.16])
    ax_note.axis("off")

    plot_check_counts(ax_a, loso, dropout, boundary, threshold, random_control)
    plot_module_dropout(ax_b, dropout)
    plot_boundary(ax_c, boundary)
    plot_human_heatmaps(ax_d1, ax_d2, human)
    ax_note.text(
        0,
        0.95,
        "Human STEMI transfer summary",
        fontsize=7,
        fontweight="bold",
        color=PALETTE["ink"],
        va="top",
    )
    ax_note.text(
        0,
        0.76,
        "Mechanical-border and fibroblast-scar axes exceeded random controls.\n"
        "The immune-fibrotic axis was directionally positive but treated as a\n"
        "contextual coupling signal. Human STEMI maps showed distinct mechanical\n"
        "hotspots and partial immune-scar overlap, supporting three outputs.",
        fontsize=6.6,
        color=PALETTE["ink"],
        va="top",
    )

    out_base = FIGURE_DIR / "nature_ed_public_data_robustness"
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out_base}.svg", bbox_inches="tight")
    fig.savefig(f"{out_base}.pdf", bbox_inches="tight")
    fig.savefig(f"{out_base}.tiff", dpi=600, bbox_inches="tight")
    fig.savefig(f"{out_base}.png", dpi=300, bbox_inches="tight")
    print(out_base)


if __name__ == "__main__":
    main()
