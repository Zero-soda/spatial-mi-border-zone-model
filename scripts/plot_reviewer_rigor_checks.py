#!/usr/bin/env python3
"""Plot reviewer-requested score-state and marker-adjustment audits."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

from project_paths import project_root

import sys


ROOT = project_root(__file__)
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.append(str(LOCAL_DEPS))

import matplotlib  # noqa: E402

matplotlib.use("Agg")

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402


TABLE_DIR = ROOT / "results" / "tables"
FIGURE_DIR = ROOT / "results" / "figures"

RECOVERY = TABLE_DIR / "gse214611_score_only_state_domain_recovery.tsv"
ADJUSTED = TABLE_DIR / "gse214611_celltype_adjusted_domain_contrasts.tsv"

OUT_STEM = FIGURE_DIR / "cvr_reviewer_rigor_checks"
FIXED_METADATA_DATE = datetime(2026, 7, 6, tzinfo=timezone.utc)

SAMPLE_ORDER = ["D3_1", "D3_2", "D3_3", "D7_1", "D7_2", "D7_3"]
SAMPLE_LABELS = {"D3_1": "D3_1", "D3_2": "D3_2", "D3_3": "D3_3", "D7_1": "D7_1", "D7_2": "D7_2", "D7_3": "D7_3"}
STAGE_COLORS = {"day3_mi": "#6B7280", "day7_mi": "#111827"}

PALETTE = {
    "mechanical_border": "#2F6FB5",
    "immune_fibrotic_activation": "#D99032",
    "fibroblast_scar_repair": "#3A9B65",
    "raw": "#111827",
    "lineage_marker_covariates": "#5B8DEF",
    "expanded_marker_covariates": "#B76E79",
    "ink": "#111827",
    "muted": "#4B5563",
    "line": "#D1D5DB",
}

SCORE_LABELS = {
    "mechanical_border": "Mechanical",
    "immune_fibrotic_activation": "Immune-fibrotic",
    "fibroblast_scar_repair": "Scar repair",
}


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 6.7,
            "axes.titlesize": 7.6,
            "axes.labelsize": 6.7,
            "xtick.labelsize": 6.0,
            "ytick.labelsize": 6.0,
            "legend.fontsize": 5.9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.55,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "svg.hashsalt": "spatial-mi-reviewer-rigor",
        }
    )


def read_table(path: Path, required: set[str]) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    data = pd.read_csv(path, sep="\t")
    missing = sorted(required.difference(data.columns))
    if missing:
        raise ValueError(f"{path.name} missing columns: {', '.join(missing)}")
    return data


def panel_label(ax: plt.Axes, label: str, x: float = -0.11, y: float = 1.10) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.8,
        fontweight="bold",
        color=PALETTE["ink"],
    )


def plot_recovery_metrics(ax: plt.Axes, recovery: pd.DataFrame) -> None:
    data = recovery[recovery["k_score_states"].astype(int) == 4].copy()
    data["sample"] = pd.Categorical(data["sample"], categories=SAMPLE_ORDER, ordered=True)
    data = data.sort_values("sample")
    x = np.arange(len(data))
    ax.plot(x, data["ari_vs_author_domains"].astype(float), marker="o", color="#2F6FB5", linewidth=1.2, markersize=4.6, label="ARI")
    ax.plot(x, data["nmi_vs_author_domains"].astype(float), marker="s", color="#3A9B65", linewidth=1.2, markersize=4.3, label="NMI")
    for idx, (_, row) in enumerate(data.iterrows()):
        ax.axvspan(idx - 0.5, idx + 0.5, color="#F3F4F6" if row["stage"] == "day3_mi" else "#FFFFFF", zorder=-2)
    ax.set_xticks(x, [SAMPLE_LABELS[sample] for sample in data["sample"].astype(str)])
    ax.set_ylim(0, 0.82)
    ax.set_ylabel("Agreement with author domains")
    ax.set_title("Score-only k = 4 states recover domain structure", loc="left", fontweight="bold")
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.45)
    ax.legend(loc="lower right", frameon=False, handlelength=1.5)
    panel_label(ax, "a")


def plot_domain_f1(ax: plt.Axes, recovery: pd.DataFrame) -> None:
    data = recovery[recovery["k_score_states"].astype(int) == 4].copy()
    data["sample"] = pd.Categorical(data["sample"], categories=SAMPLE_ORDER, ordered=True)
    data = data.sort_values("sample")
    x = np.arange(len(data))
    width = 0.34
    ax.bar(x - width / 2, data["domain3_f1"].astype(float), width=width, color="#2F6FB5", label="Domain 3")
    ax.bar(x + width / 2, data["domain4_f1"].astype(float), width=width, color="#3A9B65", label="Domain 4")
    ax.set_xticks(x, [SAMPLE_LABELS[sample] for sample in data["sample"].astype(str)])
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Best-cluster F1")
    ax.set_title("Domain 3/4 recovery without using domain labels", loc="left", fontweight="bold")
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.45)
    ax.legend(loc="lower right", frameon=False, handlelength=1.2)
    panel_label(ax, "b")


def expected_margin(row: pd.Series) -> float:
    value = float(row["marker_adjusted_mean_margin"])
    if row["score"] == "mechanical_border":
        return value
    return -value


def raw_expected_margin(row: pd.Series) -> float:
    value = float(row["raw_mean_margin"])
    if row["score"] == "mechanical_border":
        return value
    return -value


def plot_marker_adjustment(ax: plt.Axes, adjusted: pd.DataFrame) -> None:
    contrasts = [
        ("day3_mi", "mechanical_border", "D3 mech."),
        ("day7_mi", "mechanical_border", "D7 mech."),
        ("day7_mi", "fibroblast_scar_repair", "D7 scar"),
        ("day7_mi", "immune_fibrotic_activation", "D7 immune"),
    ]
    rows = []
    for stage, score, label in contrasts:
        subset = adjusted[(adjusted["stage"] == stage) & (adjusted["score"] == score)]
        raw_row = subset.iloc[0]
        rows.append({"label": label, "mode": "Raw", "value": raw_expected_margin(raw_row), "score": score})
        for mode, mode_label in [
            ("lineage_marker_covariates", "Lineage adjusted"),
            ("expanded_marker_covariates", "Expanded adjusted"),
        ]:
            row = subset[subset["adjustment_mode"] == mode].iloc[0]
            rows.append({"label": label, "mode": mode_label, "value": expected_margin(row), "score": score})

    labels = [label for _, _, label in contrasts]
    x = np.arange(len(labels))
    offsets = {"Raw": -0.24, "Lineage adjusted": 0.0, "Expanded adjusted": 0.24}
    colors = {"Raw": PALETTE["raw"], "Lineage adjusted": PALETTE["lineage_marker_covariates"], "Expanded adjusted": PALETTE["expanded_marker_covariates"]}
    for mode in offsets:
        values = [row["value"] for row in rows if row["mode"] == mode]
        ax.bar(x + offsets[mode], values, width=0.21, color=colors[mode], label=mode)
    ax.axhline(0, color=PALETTE["muted"], linewidth=0.65)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Expected-direction margin")
    ax.set_title("Marker adjustment attenuates stromal axes", loc="left", fontweight="bold")
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.45)
    ax.legend(loc="upper right", frameon=False, ncol=1, handlelength=1.1)
    panel_label(ax, "c")


def plot_interpretive_matrix(ax: plt.Axes, adjusted: pd.DataFrame) -> None:
    ax.set_axis_off()
    panel_label(ax, "d", x=-0.10, y=1.10)
    ax.text(
        0.0,
        1.0,
        "Evidence status after reviewer-requested audits",
        ha="left",
        va="top",
        fontsize=7.6,
        fontweight="bold",
        color=PALETTE["ink"],
    )
    rows = [
        ("Mechanical-border", "Strong", "Score-only recovery and adjusted margins remain positive."),
        ("Fibroblast-scar repair", "Strong tissue-state;\ncomposition-coupled", "Domain 4 enrichment is robust but marker-sensitive."),
        ("Immune-fibrotic", "Moderate contextual", "Related to scar repair and sensitive to adjustment."),
        ("Human transfer", "Feasibility only", "Single human STEMI section; no outcome labels."),
    ]
    y = 0.82
    for label, status, note in rows:
        ax.text(0.02, y, label, ha="left", va="top", fontsize=6.2, fontweight="bold", color=PALETTE["ink"])
        ax.text(0.39, y, status, ha="left", va="top", fontsize=5.95, color=PALETTE["muted"], linespacing=1.05)
        ax.text(0.39, y - 0.092, note, ha="left", va="top", fontsize=5.65, color=PALETTE["muted"], linespacing=1.05)
        ax.plot([0.02, 0.98], [y - 0.155, y - 0.155], color="#E5E7EB", linewidth=0.5)
        y -= 0.205
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def save_outputs(fig: plt.Figure) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    for suffix in [".png", ".tiff", ".pdf", ".svg"]:
        kwargs: dict[str, object] = {"bbox_inches": "tight", "pad_inches": 0.045}
        if suffix in {".png", ".tiff"}:
            kwargs["dpi"] = 600
        if suffix == ".tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        if suffix == ".pdf":
            kwargs["metadata"] = {
                "Title": "Reviewer-requested rigor checks",
                "Creator": "scripts/plot_reviewer_rigor_checks.py",
                "CreationDate": FIXED_METADATA_DATE,
                "ModDate": FIXED_METADATA_DATE,
            }
        if suffix == ".svg":
            kwargs["metadata"] = {
                "Title": "Reviewer-requested rigor checks",
                "Creator": "scripts/plot_reviewer_rigor_checks.py",
                "Date": "2026-07-06T00:00:00Z",
            }
        fig.savefig(OUT_STEM.with_suffix(suffix), **kwargs)
    svg = OUT_STEM.with_suffix(".svg")
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines()) + "\n", encoding="utf-8")


def main() -> None:
    setup_style()
    recovery = read_table(
        RECOVERY,
        {
            "sample",
            "stage",
            "k_score_states",
            "ari_vs_author_domains",
            "nmi_vs_author_domains",
            "domain3_f1",
            "domain4_f1",
        },
    )
    adjusted = read_table(
        ADJUSTED,
        {
            "adjustment_mode",
            "stage",
            "score",
            "raw_mean_margin",
            "marker_adjusted_mean_margin",
        },
    )
    fig = plt.figure(figsize=(7.1, 6.55))
    grid = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.38)
    plot_recovery_metrics(fig.add_subplot(grid[0, 0]), recovery)
    plot_domain_f1(fig.add_subplot(grid[0, 1]), recovery)
    plot_marker_adjustment(fig.add_subplot(grid[1, 0]), adjusted)
    plot_interpretive_matrix(fig.add_subplot(grid[1, 1]), adjusted)
    fig.suptitle(
        "Reviewer-requested rigor checks for domain dependence and cell composition",
        x=0.02,
        y=0.995,
        ha="left",
        fontsize=10.0,
        fontweight="bold",
        color=PALETTE["ink"],
    )
    fig.text(
        0.02,
        0.962,
        "Score-only clustering ignores author domain labels; marker adjustment is an audit, not deconvolution.",
        ha="left",
        va="top",
        fontsize=7.1,
        color=PALETTE["muted"],
    )
    fig.subplots_adjust(left=0.075, right=0.980, bottom=0.075, top=0.890)
    save_outputs(fig)
    plt.close(fig)
    print(f"[figure] wrote {OUT_STEM.with_suffix('.png')}")


if __name__ == "__main__":
    main()
