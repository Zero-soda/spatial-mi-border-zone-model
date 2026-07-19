#!/usr/bin/env python3
"""Plot CVR-style spatial statistics for the MI domain-state revision."""

from __future__ import annotations

from datetime import datetime, timezone
import sys
from pathlib import Path

from project_paths import project_root


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

EFFECT = TABLE_DIR / "gse214611_stage_domain_effect_sizes.tsv"
AUTOCORR = TABLE_DIR / "gse214611_spatial_autocorrelation.tsv"
PERMUTATION = TABLE_DIR / "gse214611_domain_label_permutation.tsv"
GRADIENT = TABLE_DIR / "gse214611_boundary_distance_gradients.tsv"
KEY_GENE = TABLE_DIR / "gse214611_key_gene_removal_sensitivity.tsv"

OUT_STEM = FIGURE_DIR / "cvr_spatial_statistics_upgrade"
FIGURE_TITLE = "CVR spatial statistics upgrade"
FIGURE_CREATOR = "scripts/plot_cvr_spatial_statistics.py"
FIXED_METADATA_DATE = datetime(2026, 7, 6, tzinfo=timezone.utc)

STAGE_ORDER = ["day3_mi", "day7_mi"]
SCORE_ORDER = [
    "mechanical_border",
    "fibroblast_scar_repair",
    "immune_fibrotic_activation_proxy",
]

STAGE_LABELS = {
    "day3_mi": "D3 MI",
    "day7_mi": "D7 MI",
}

SCORE_LABELS = {
    "mechanical_border": "Mechanical",
    "fibroblast_scar_repair": "Scar repair",
    "immune_fibrotic_activation_proxy": "Immune-fibrotic",
    "immune_fibrotic_activation": "Immune-fibrotic",
}

OUTPUT_LABELS = {
    "mechanical_border": "Mechanical",
    "fibroblast_scar_repair": "Scar repair",
    "immune_fibrotic_activation": "Immune-fibrotic",
}

PALETTE = {
    "mechanical_border": "#2F6FB5",
    "fibroblast_scar_repair": "#3A9B65",
    "immune_fibrotic_activation_proxy": "#D99032",
    "immune_fibrotic_activation": "#D99032",
    "ink": "#111827",
    "muted": "#4B5563",
    "gray": "#6B7280",
    "line": "#D1D5DB",
    "light": "#EEF2F6",
}


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 6.5,
            "axes.titlesize": 7.3,
            "axes.labelsize": 6.6,
            "xtick.labelsize": 5.9,
            "ytick.labelsize": 5.9,
            "legend.fontsize": 5.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.55,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "svg.hashsalt": "spatial-mi-cvr-statistics",
        }
    )


def read_table(path: Path, required_columns: set[str]) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    data = pd.read_csv(path, sep="\t")
    missing = sorted(required_columns.difference(data.columns))
    if missing:
        raise ValueError(f"{path.name} is missing required columns: {', '.join(missing)}")
    if data.empty:
        raise ValueError(f"{path.name} is empty")
    return data


def load_tables() -> dict[str, pd.DataFrame]:
    return {
        "effect": read_table(
            EFFECT,
            {
                "stage",
                "score",
                "n_sample_pairs",
                "mean_margin",
                "ci95_low",
                "ci95_high",
            },
        ),
        "autocorr": read_table(AUTOCORR, {"sample", "stage", "score", "morans_i"}),
        "permutation": read_table(
            PERMUTATION,
            {
                "sample",
                "stage",
                "score",
                "observed_margin",
                "null_p99",
                "block_null_p99",
            },
        ),
        "gradient": read_table(
            GRADIENT,
            {
                "sample",
                "stage",
                "score",
                "signed_distance_definition",
                "slope_per_1000_fullres_px",
            },
        ),
        "key_gene": read_table(
            KEY_GENE,
            {
                "key_gene",
                "output",
                "stage",
                "min_expected_direction_margin",
                "direction_preserved",
            },
        ),
    }


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.11) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.6,
        fontweight="bold",
        color=PALETTE["ink"],
    )


def group_label(stage: str, score: str) -> str:
    return f"{STAGE_LABELS[stage]} {SCORE_LABELS[score]}"


def ordered_groups(data: pd.DataFrame) -> list[tuple[str, str, pd.DataFrame]]:
    groups = []
    for stage in STAGE_ORDER:
        for score in SCORE_ORDER:
            subset = data[(data["stage"] == stage) & (data["score"] == score)].copy()
            if not subset.empty:
                groups.append((stage, score, subset))
    return groups


def add_y_grid(ax: plt.Axes) -> None:
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.45)
    ax.set_axisbelow(True)


def plot_effects(ax: plt.Axes, data: pd.DataFrame) -> None:
    rows = []
    for stage, score, subset in ordered_groups(data):
        row = subset.iloc[0]
        rows.append(
            {
                "label": group_label(stage, score),
                "score": score,
                "mean": float(row["mean_margin"]),
                "low": float(row["ci95_low"]),
                "high": float(row["ci95_high"]),
                "n": int(row["n_sample_pairs"]),
            }
        )

    y = np.arange(len(rows))
    means = np.array([row["mean"] for row in rows])
    lows = np.array([row["low"] for row in rows])
    highs = np.array([row["high"] for row in rows])
    colors = [PALETTE[row["score"]] for row in rows]

    ax.axvline(0, color=PALETTE["gray"], linewidth=0.65, zorder=1)
    ax.errorbar(
        means,
        y,
        xerr=[means - lows, highs - means],
        fmt="none",
        ecolor=PALETTE["gray"],
        elinewidth=0.75,
        capsize=2.2,
        capthick=0.75,
        zorder=2,
    )
    ax.scatter(means, y, s=30, color=colors, edgecolors="white", linewidths=0.45, zorder=3)

    ax.set_yticks(y, [row["label"] for row in rows])
    ax.invert_yaxis()
    ax.set_xlabel("Domain 3 - domain 4 margin")
    ax.set_title("Sample-level effects, bootstrap 95% CI", loc="left", fontweight="bold")
    ax.text(
        0.98,
        0.06,
        "n = 3 paired samples per stage",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.7,
        color=PALETTE["muted"],
    )
    ax.set_xlim(min(lows.min() - 0.45, -0.6), max(highs.max() + 0.45, 0.8))
    add_y_grid(ax)
    panel_label(ax, "a")


def plot_autocorrelation(ax: plt.Axes, data: pd.DataFrame) -> None:
    violin_data = [
        data.loc[data["score"] == score, "morans_i"].astype(float).to_numpy()
        for score in SCORE_ORDER
    ]
    positions = np.arange(len(SCORE_ORDER))
    parts = ax.violinplot(violin_data, positions=positions, widths=0.64, showextrema=False)
    for body, score in zip(parts["bodies"], SCORE_ORDER, strict=True):
        body.set_facecolor(PALETTE[score])
        body.set_edgecolor(PALETTE[score])
        body.set_alpha(0.22)
        body.set_linewidth(0.45)

    stage_offsets = {"day3_mi": -0.085, "day7_mi": 0.085}
    stage_markers = {"day3_mi": "o", "day7_mi": "s"}
    for idx, score in enumerate(SCORE_ORDER):
        subset = data[data["score"] == score].copy()
        for stage in STAGE_ORDER:
            values = subset.loc[subset["stage"] == stage, "morans_i"].astype(float).to_numpy()
            x = np.full(len(values), idx + stage_offsets[stage])
            ax.scatter(
                x,
                values,
                s=23,
                marker=stage_markers[stage],
                color=PALETTE[score],
                edgecolors="white",
                linewidths=0.42,
                zorder=3,
            )
        mean = float(subset["morans_i"].mean())
        ax.plot([idx - 0.19, idx + 0.19], [mean, mean], color=PALETTE["ink"], linewidth=0.8, zorder=4)

    ax.set_xticks(positions, [SCORE_LABELS[score] for score in SCORE_ORDER])
    ax.set_ylabel("Moran's I")
    ax.set_ylim(0.52, 0.98)
    ax.set_title("Spatial autocorrelation across samples", loc="left", fontweight="bold")
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=PALETTE["gray"], markeredgecolor="white", markersize=5, label="D3 MI"),
            Line2D([0], [0], marker="s", color="none", markerfacecolor=PALETTE["gray"], markeredgecolor="white", markersize=5, label="D7 MI"),
            Line2D([0], [0], color=PALETTE["ink"], linewidth=0.8, label="mean"),
        ],
        loc="lower right",
        handlelength=1.3,
        borderpad=0.2,
        labelspacing=0.35,
    )
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.45)
    ax.set_axisbelow(True)
    panel_label(ax, "b")


def plot_permutation(ax: plt.Axes, data: pd.DataFrame) -> None:
    groups = ordered_groups(data)
    y_positions = np.arange(len(groups))

    for y, (stage, score, subset) in zip(y_positions, groups, strict=True):
        observed = subset["observed_margin"].astype(float).to_numpy()
        label_p99 = float(subset["null_p99"].astype(float).abs().max())
        block_p99 = float(subset["block_null_p99"].astype(float).abs().max())

        ax.plot(
            [-block_p99, block_p99],
            [y, y],
            color="#C9D0DA",
            linewidth=7.0,
            solid_capstyle="round",
            zorder=1,
        )
        ax.plot(
            [-label_p99, label_p99],
            [y, y],
            color="#788392",
            linewidth=4.2,
            solid_capstyle="round",
            zorder=2,
        )
        offsets = np.linspace(-0.13, 0.13, len(observed)) if len(observed) > 1 else np.array([0.0])
        ax.scatter(
            observed,
            y + offsets,
            s=25,
            color=PALETTE[score],
            edgecolors="white",
            linewidths=0.42,
            zorder=3,
        )

    max_extent = max(
        data["observed_margin"].astype(float).abs().max(),
        data["block_null_p99"].astype(float).abs().max(),
    )
    ax.axvline(0, color=PALETTE["gray"], linewidth=0.65, zorder=0)
    ax.set_yticks(y_positions, [group_label(stage, score) for stage, score, _ in groups])
    ax.invert_yaxis()
    ax.set_xlim(-max_extent * 1.18, max_extent * 1.18)
    ax.set_xlabel("Signed D3 - D4 margin (bars = group max +/- p99 null)")
    ax.set_title("Domain-label and spatial-block permutation checks", loc="left", fontweight="bold")
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=PALETTE["gray"], markeredgecolor="white", markersize=4.8, label="observed samples"),
            Line2D([0], [0], color="#788392", linewidth=4.2, label="label null p99"),
            Line2D([0], [0], color="#C9D0DA", linewidth=7.0, label="block null p99"),
        ],
        loc="lower right",
        bbox_to_anchor=(1.0, 1.02),
        ncol=3,
        columnspacing=1.0,
        handlelength=1.8,
        borderpad=0.2,
        frameon=False,
    )
    add_y_grid(ax)
    panel_label(ax, "c", x=-0.055, y=1.08)


def plot_gradients(ax: plt.Axes, data: pd.DataFrame) -> None:
    definitions = sorted(data["signed_distance_definition"].dropna().unique())
    expected = ["domain3_side_positive_domain4_side_negative"]
    if definitions != expected:
        raise ValueError(f"Unexpected signed-distance definition: {definitions}")

    groups = ordered_groups(data)
    y_positions = np.arange(len(groups))
    for y, (stage, score, subset) in zip(y_positions, groups, strict=True):
        slopes = subset["slope_per_1000_fullres_px"].astype(float).to_numpy()
        offsets = np.linspace(-0.13, 0.13, len(slopes)) if len(slopes) > 1 else np.array([0.0])
        ax.scatter(
            slopes,
            y + offsets,
            s=25,
            color=PALETTE[score],
            edgecolors="white",
            linewidths=0.42,
            zorder=3,
        )
        ax.scatter(
            [float(np.mean(slopes))],
            [y],
            marker="D",
            s=20,
            color=PALETTE["ink"],
            edgecolors="white",
            linewidths=0.35,
            zorder=4,
        )

    max_extent = data["slope_per_1000_fullres_px"].astype(float).abs().max()
    ax.axvline(0, color=PALETTE["gray"], linewidth=0.65, zorder=1)
    ax.set_xlim(-max_extent * 1.20, max_extent * 1.20)
    ax.set_yticks(y_positions, [group_label(stage, score) for stage, score, _ in groups])
    ax.invert_yaxis()
    ax.set_xlabel("Slope per 1,000 full-res px\npositive increases toward domain 3; domain 4 side is negative")
    ax.set_title("Signed boundary-distance gradients", loc="left", fontweight="bold")
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=PALETTE["gray"], markeredgecolor="white", markersize=4.8, label="sample"),
            Line2D([0], [0], marker="D", color="none", markerfacecolor=PALETTE["ink"], markeredgecolor="white", markersize=4.8, label="mean"),
        ],
        loc="lower right",
        handlelength=1.0,
        borderpad=0.2,
        labelspacing=0.35,
    )
    add_y_grid(ax)
    panel_label(ax, "d")


def gene_order(data: pd.DataFrame) -> list[str]:
    seen: list[str] = []
    for gene in data["key_gene"]:
        if gene not in seen:
            seen.append(gene)
    return seen


def plot_key_gene(ax: plt.Axes, data: pd.DataFrame) -> None:
    genes = gene_order(data)
    x_positions = {gene: idx for idx, gene in enumerate(genes)}
    for gene in genes:
        subset = data[data["key_gene"] == gene].copy()
        offsets = np.linspace(-0.20, 0.20, len(subset)) if len(subset) > 1 else np.array([0.0])
        for offset, (_, row) in zip(offsets, subset.iterrows(), strict=True):
            output = row["output"]
            ax.scatter(
                x_positions[gene] + offset,
                float(row["min_expected_direction_margin"]),
                s=28,
                color=PALETTE[output],
                edgecolors="white",
                linewidths=0.42,
                zorder=3,
            )

    preserved = int((data["direction_preserved"].astype(str) == "True").sum())
    total = len(data)
    counts = data.groupby("key_gene")["direction_preserved"].agg(lambda values: int((values.astype(str) == "True").sum()))
    totals = data.groupby("key_gene")["direction_preserved"].size()
    y_max = float(data["min_expected_direction_margin"].max()) + 0.42

    for gene in genes:
        ax.text(
            x_positions[gene],
            y_max - 0.05,
            f"{counts[gene]}/{totals[gene]}",
            ha="center",
            va="top",
            fontsize=5.4,
            color=PALETTE["muted"],
        )

    ax.axhline(0, color=PALETTE["gray"], linewidth=0.65, zorder=1)
    ax.set_xticks(range(len(genes)), genes, rotation=36, ha="right")
    ax.set_ylim(0, y_max)
    ax.set_ylabel("Minimum expected-direction margin")
    ax.set_title("Key-gene removal sensitivity", loc="left", fontweight="bold")
    ax.text(
        0.02,
        0.92,
        f"{preserved}/{total} removals preserve direction",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.8,
        color=PALETTE["muted"],
    )
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=PALETTE[output], markeredgecolor="white", markersize=5, label=label)
            for output, label in OUTPUT_LABELS.items()
        ],
        loc="lower right",
        borderpad=0.2,
        labelspacing=0.3,
        handlelength=1.0,
    )
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.45)
    ax.set_axisbelow(True)
    panel_label(ax, "e")


def save_outputs(fig: plt.Figure) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    for suffix in [".png", ".tiff", ".pdf", ".svg"]:
        kwargs: dict[str, object] = {"bbox_inches": "tight", "pad_inches": 0.045}
        if suffix in {".png", ".tiff"}:
            kwargs["dpi"] = 600
        if suffix == ".pdf":
            kwargs["metadata"] = {
                "Title": FIGURE_TITLE,
                "Creator": FIGURE_CREATOR,
                "CreationDate": FIXED_METADATA_DATE,
                "ModDate": FIXED_METADATA_DATE,
            }
        if suffix == ".svg":
            kwargs["metadata"] = {
                "Title": FIGURE_TITLE,
                "Creator": FIGURE_CREATOR,
                "Date": "2026-07-06T00:00:00Z",
            }
        if suffix == ".tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(OUT_STEM.with_suffix(suffix), **kwargs)
    strip_trailing_whitespace(OUT_STEM.with_suffix(".svg"))


def strip_trailing_whitespace(path: Path) -> None:
    cleaned = "\n".join(line.rstrip() for line in path.read_text(encoding="utf-8").splitlines()) + "\n"
    path.write_text(cleaned, encoding="utf-8")


def main() -> None:
    setup_style()
    tables = load_tables()

    fig = plt.figure(figsize=(183 / 25.4, 218 / 25.4), constrained_layout=False)
    grid = fig.add_gridspec(
        nrows=3,
        ncols=2,
        height_ratios=[1.02, 1.22, 1.08],
        width_ratios=[1.04, 1.10],
        hspace=0.68,
        wspace=0.46,
    )

    plot_effects(fig.add_subplot(grid[0, 0]), tables["effect"])
    plot_autocorrelation(fig.add_subplot(grid[0, 1]), tables["autocorr"])
    plot_permutation(fig.add_subplot(grid[1, :]), tables["permutation"])
    plot_gradients(fig.add_subplot(grid[2, 0]), tables["gradient"])
    plot_key_gene(fig.add_subplot(grid[2, 1]), tables["key_gene"])

    fig.suptitle(
        "Spatial robustness of the domain-state framework",
        x=0.055,
        y=0.983,
        ha="left",
        fontsize=9.8,
        fontweight="bold",
        color=PALETTE["ink"],
    )
    fig.text(
        0.055,
        0.962,
        "Extended-data quantitative grid: sample effects, spatial autocorrelation, permutation checks, boundary gradients and key-gene removal.",
        ha="left",
        va="top",
        fontsize=6.5,
        color=PALETTE["muted"],
    )
    fig.subplots_adjust(left=0.112, right=0.983, bottom=0.085, top=0.908)

    save_outputs(fig)
    plt.close(fig)
    print(f"Wrote {OUT_STEM.with_suffix('.png')}")
    print(f"Wrote {OUT_STEM.with_suffix('.tiff')}")
    print(f"Wrote {OUT_STEM.with_suffix('.pdf')}")
    print(f"Wrote {OUT_STEM.with_suffix('.svg')}")


if __name__ == "__main__":
    main()
