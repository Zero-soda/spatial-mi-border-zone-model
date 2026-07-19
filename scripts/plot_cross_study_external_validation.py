#!/usr/bin/env python3
"""Generate Figure 7 and Supplementary Figure S9 for cross-study validation."""

from __future__ import annotations

import shutil
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
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402


TABLES = ROOT / "results" / "tables" / "external_validation"
FIGURES = ROOT / "results" / "figures"

MOUSE_SECTION = TABLES / "external_mouse_section_summary.tsv"
AUTHOR_AREA = TABLES / "external_mouse_author_area_summary.tsv"
HUMAN_SECTION = TABLES / "kuppe_human_section_summary.tsv"
HUMAN_SPOT = TABLES / "kuppe_human_spatial_scores_by_spot.tsv"
BENCHMARK = TABLES / "cross_study_model_benchmark.tsv"
BENCHMARK_GAIN = TABLES / "cross_study_model_benchmark_gain.tsv"
TRANSITION = TABLES / "cross_study_boundary_transition_index.tsv"

FIG7 = FIGURES / "Figure_7_cross_study_external_validation"
FIGS9 = FIGURES / "Supplementary_Figure_S9_external_validation_robustness"

BLUE = "#2F6FB5"
ORANGE = "#D99032"
GREEN = "#3A9B65"
PURPLE = "#6B5FB5"
RED = "#C64F48"
INK = "#17212F"
MUTED = "#667085"
LIGHT = "#E4E9F0"
PALE = "#F6F8FA"
GROUP_COLORS = {"control": "#8A94A3", "border": BLUE, "fibrotic": GREEN}
AXIS_COLORS = {"mechanical": BLUE, "transition": PURPLE, "scar": GREEN}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 600,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 6.2,
            "axes.titlesize": 7.2,
            "axes.labelsize": 6.4,
            "xtick.labelsize": 5.8,
            "ytick.labelsize": 5.8,
            "axes.linewidth": 0.55,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def panel_label(ax: plt.Axes, label: str, x: float = -0.08, y: float = 1.05) -> None:
    ax.text(x, y, label, transform=ax.transAxes, ha="left", va="top", fontsize=9, fontweight="bold", color=INK)


def save_figure(fig: plt.Figure, base: Path) -> None:
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(base.with_suffix(".pdf"), facecolor="white")
    fig.savefig(base.with_suffix(".svg"), facecolor="white")
    fig.savefig(
        base.with_suffix(".tiff"),
        dpi=600,
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )


def copy_source_data() -> None:
    mapping = {
        MOUSE_SECTION: TABLES.parent / "Source_Data_Figure_7_external_mouse_sections.tsv",
        AUTHOR_AREA: TABLES.parent / "Source_Data_Figure_7_author_area_gradients.tsv",
        HUMAN_SECTION: TABLES.parent / "Source_Data_Figure_7_kuppe_patient_sections.tsv",
        BENCHMARK: TABLES.parent / "Source_Data_Figure_7_model_benchmark.tsv",
        TRANSITION: TABLES.parent / "Source_Data_Figure_7_boundary_transition_index.tsv",
        HUMAN_SPOT: TABLES.parent / "Source_Data_Supplementary_Figure_S9_kuppe_spatial_scores.tsv",
        BENCHMARK_GAIN: TABLES.parent / "Source_Data_Supplementary_Figure_S9_benchmark_sensitivity.tsv",
    }
    for source, destination in mapping.items():
        shutil.copy2(source, destination)


def draw_transfer_design(ax: plt.Axes, human: pd.DataFrame) -> None:
    panel_label(ax, "a", x=-0.015, y=1.03)
    ax.set_axis_off()
    boxes = [
        (0.01, 0.20, 0.20, BLUE, "Frozen discovery", "GSE214611\n3 outputs + BTI"),
        (0.28, 0.20, 0.19, PURPLE, "Mouse time course", "GSE176092\n9 MI sections"),
        (0.53, 0.20, 0.19, ORANGE, "Author anatomy", "GSE265828\nRZ-BZ1-BZ2-IZ"),
        (
            0.78,
            0.20,
            0.20,
            GREEN,
            "Human transfer",
            f"Kuppe atlas\n9 patients; {int(human['n_spots'].sum()):,} spots",
        ),
    ]
    for x, y, w, color, title, body in boxes:
        patch = FancyBboxPatch(
            (x, y),
            w,
            0.60,
            boxstyle="round,pad=0.012,rounding_size=0.025",
            facecolor="white",
            edgecolor=color,
            linewidth=1.0,
        )
        ax.add_patch(patch)
        ax.text(x + 0.016, 0.70, title, ha="left", va="top", fontsize=6.8, fontweight="bold", color=INK)
        ax.text(x + 0.016, 0.54, body, ha="left", va="top", fontsize=5.5, color=MUTED, linespacing=1.12)
    for start, end in ((0.21, 0.28), (0.47, 0.53), (0.72, 0.78)):
        ax.add_patch(FancyArrowPatch((start, 0.50), (end, 0.50), arrowstyle="-|>", mutation_scale=8, color=MUTED, linewidth=0.9))
    ax.text(
        0.01,
        0.035,
        "No retraining | biological section/patient is the inferential unit | spot-level data support maps and within-section spatial statistics",
        ha="left",
        va="center",
        color=MUTED,
        fontsize=5.7,
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def plot_yamada_transition(ax: plt.Axes, mouse: pd.DataFrame) -> None:
    panel_label(ax, "b", x=-0.14, y=1.08)
    data = mouse.loc[mouse["dataset"] == "GSE176092"].copy()
    order = ["day1", "day7", "day14"]
    data["stage_x"] = data["stage_region"].map({stage: idx for idx, stage in enumerate(order)})
    data["framework"] = data["graph_morans_i_mechanical_border_score"] - data["graph_morans_i_fibroblast_scar_repair_score"]
    data["baseline"] = data["graph_morans_i_source_bz_baseline"] - data["graph_morans_i_generic_fibrosis_baseline"]
    jitter = np.asarray([-0.055, 0.0, 0.055])
    for column, color, marker, label, linestyle in (
        ("framework", BLUE, "o", "Three-output framework", "-"),
        ("baseline", "#9099A7", "s", "Two simple baselines", "--"),
    ):
        means = []
        for stage_idx, stage in enumerate(order):
            values = data.loc[data["stage_region"] == stage, column].to_numpy(float)
            ax.scatter(stage_idx + jitter[: len(values)], values, s=15, color=color, marker=marker, edgecolor="white", linewidth=0.35, zorder=3)
            means.append(float(np.mean(values)))
        ax.plot(range(3), means, color=color, linestyle=linestyle, linewidth=1.25, marker=marker, markersize=3.5, label=label)
    ax.axhline(0, color="#AAB2BE", linewidth=0.7)
    ax.set_xticks(range(3), ["D1", "D7", "D14"])
    ax.set_ylabel("Moran's I: mechanical - scar")
    ax.set_title("Independent mouse time course", loc="left", fontweight="bold", pad=5)
    ax.legend(loc="lower right", fontsize=5.2, handlelength=2.1)
    ax.text(0.02, 0.97, "n=3 sections/stage", transform=ax.transAxes, ha="left", va="top", color=MUTED, fontsize=5.4)


def plot_author_area_profiles(container: plt.Axes, area: pd.DataFrame) -> None:
    panel_label(container, "c", x=-0.06, y=1.20)
    container.set_axis_off()
    container.text(0.0, 1.15, "Independent author-labelled anatomy", transform=container.transAxes, ha="left", va="bottom", fontsize=7.2, fontweight="bold", color=INK)
    sub = container.get_subplotspec().subgridspec(1, 3, wspace=0.35)
    area_order = ["RZ", "BZ1", "BZ2", "IZ"]
    panels = (
        ("mean_mechanical_border_score", "Mechanical", BLUE),
        ("mean_fibroblast_scar_repair_score", "Scar repair", GREEN),
        ("mean_boundary_transition_index", "BTI", PURPLE),
    )
    for idx, (column, title, color) in enumerate(panels):
        ax = container.figure.add_subplot(sub[0, idx])
        for sample, label, linestyle, marker in (
            ("GSM8229743_dpi3", "D3", "-", "o"),
            ("GSM8229745_dpi5_male", "D5", "--", "s"),
        ):
            rows = area.loc[area["sample"] == sample].set_index("author_area").reindex(area_order)
            ax.plot(range(4), rows[column], color=color, linestyle=linestyle, marker=marker, linewidth=1.15, markersize=3.4, label=label)
        ax.axhline(0, color=LIGHT, linewidth=0.6)
        ax.set_xticks(range(4), area_order)
        ax.set_title(title, fontsize=6.4, fontweight="bold", color=color, pad=3)
        if idx != 0:
            ax.set_yticklabels([])
        if idx == 2:
            ax.legend(loc="upper left", fontsize=5.0, handlelength=1.7)


def plot_patient_transfer(container: plt.Axes, human: pd.DataFrame, transition: pd.DataFrame) -> None:
    panel_label(container, "d", x=-0.035, y=1.18)
    container.set_axis_off()
    container.text(0.0, 1.13, "Nine-patient human transfer", transform=container.transAxes, ha="left", va="bottom", fontsize=7.2, fontweight="bold", color=INK)
    joined = human.merge(
        transition.loc[transition["dataset"] == "KUPPE2022", ["sample", "section_boundary_transition_index"]],
        on="sample",
        how="left",
    )
    sub = container.get_subplotspec().subgridspec(1, 3, wspace=0.33)
    panels = (
        ("mean_raw_mechanical_border_score", "Mechanical", BLUE),
        ("mean_raw_fibroblast_scar_repair_score", "Scar repair", GREEN),
        ("section_boundary_transition_index", "BTI", PURPLE),
    )
    groups = ["control", "border", "fibrotic"]
    for idx, (column, title, color) in enumerate(panels):
        ax = container.figure.add_subplot(sub[0, idx])
        for group_idx, group in enumerate(groups):
            values = joined.loc[joined["coarse_region"] == group, column].to_numpy(float)
            jitter = np.asarray([-0.08, 0.0, 0.08])[: len(values)]
            ax.scatter(group_idx + jitter, values, color=GROUP_COLORS[group], s=25, edgecolor="white", linewidth=0.45, zorder=3)
            mean = float(np.mean(values))
            ax.plot([group_idx - 0.16, group_idx + 0.16], [mean, mean], color=INK, linewidth=1.25, zorder=4)
        ax.axhline(0, color=LIGHT, linewidth=0.6)
        ax.set_xticks(range(3), ["Control", "Border", "Fibrotic"], rotation=25, ha="right")
        ax.set_title(title, fontsize=6.5, fontweight="bold", color=color, pad=3)
        if idx == 0:
            ax.set_ylabel("Section-level score")
        ax.text(0.02, 0.97, "n=3 patients/group", transform=ax.transAxes, ha="left", va="top", fontsize=5.1, color=MUTED)


def compact_comparison(value: str) -> str:
    mapping = {
        "Yamada prespecified tissue-mean mechanical": "Yamada D1>D14\nwhole-section mech.",
        "Yamada spatial mechanical-to-scar transition": "Yamada D1>D7\nspatial transition",
        "Yamada day7 scar": "Yamada D7>D1\nscar repair",
        "Hernandez scar maturation": "Hernandez D5>early\nscar repair",
        "Kuppe border-ischaemic mechanical": "Kuppe border>control\nmechanical",
        "Kuppe fibrotic scar": "Kuppe fibrotic>control\nscar repair",
    }
    return mapping.get(value, value)


def plot_benchmark(ax: plt.Axes, benchmark: pd.DataFrame) -> None:
    panel_label(ax, "e", x=-0.055, y=1.14)
    order = list(dict.fromkeys(benchmark["comparison"]))
    y = np.arange(len(order))[::-1]
    full = benchmark.loc[benchmark["model"] == "three_output_framework"].set_index("comparison")
    base = benchmark.loc[benchmark["model"] != "three_output_framework"].set_index("comparison")
    for pos, comparison in zip(y, order, strict=True):
        f = float(full.loc[comparison, "rank_biserial_effect"])
        b = float(base.loc[comparison, "rank_biserial_effect"])
        axis = str(full.loc[comparison, "axis"])
        color = AXIS_COLORS.get(axis, BLUE)
        ax.plot([b, f], [pos, pos], color="#C8CED7", linewidth=1.0, zorder=1)
        ax.scatter(b, pos, color="#8B95A3", marker="s", s=20, edgecolor="white", linewidth=0.35, zorder=2)
        ax.scatter(f, pos, color=color, marker="o", s=24, edgecolor="white", linewidth=0.35, zorder=3)
    ax.axvline(0, color="#AAB2BE", linewidth=0.7)
    ax.set_xlim(-1.08, 1.08)
    ax.set_yticks(y, [compact_comparison(value) for value in order])
    ax.set_xlabel("Rank-biserial effect (prespecified positive direction)")
    ax.set_title("Transparent comparison with simple baselines", loc="left", fontweight="bold", pad=5)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=BLUE, markeredgecolor="white", label="Three-output framework"),
            Line2D([0], [0], marker="s", color="none", markerfacecolor="#8B95A3", markeredgecolor="white", label="Simple baseline"),
        ],
        loc="lower right",
        fontsize=5.2,
        frameon=False,
    )
    ax.text(0.01, 0.02, "Negative first row retained as prespecified falsification", transform=ax.transAxes, ha="left", va="bottom", fontsize=5.2, color=RED)


def build_main_figure(mouse: pd.DataFrame, area: pd.DataFrame, human: pd.DataFrame, benchmark: pd.DataFrame, transition: pd.DataFrame) -> plt.Figure:
    fig = plt.figure(figsize=(7.2, 7.05), facecolor="white")
    grid = fig.add_gridspec(4, 12, height_ratios=[1.0, 2.0, 2.0, 2.35], hspace=0.86, wspace=0.72)
    draw_transfer_design(fig.add_subplot(grid[0, :]), human)
    plot_yamada_transition(fig.add_subplot(grid[1, :5]), mouse)
    plot_author_area_profiles(fig.add_subplot(grid[1, 6:]), area)
    plot_patient_transfer(fig.add_subplot(grid[2, :]), human, transition)
    plot_benchmark(fig.add_subplot(grid[3, :]), benchmark)
    fig.subplots_adjust(left=0.17, right=0.98, top=0.98, bottom=0.09)
    return fig


def plot_spatial_map(ax: plt.Axes, rows: pd.DataFrame, norm: TwoSlopeNorm, cmap) -> None:
    order = np.argsort(np.abs(rows["boundary_transition_index"].to_numpy(float)))
    ax.scatter(
        rows["spatial_x"].to_numpy(float)[order],
        -rows["spatial_y"].to_numpy(float)[order],
        c=rows["boundary_transition_index"].to_numpy(float)[order],
        cmap=cmap,
        norm=norm,
        s=2.2,
        linewidths=0,
        rasterized=True,
    )
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def build_supplementary_figure(human: pd.DataFrame, spots: pd.DataFrame, benchmark_gain: pd.DataFrame) -> plt.Figure:
    fig = plt.figure(figsize=(7.2, 9.0), facecolor="white")
    grid = fig.add_gridspec(5, 3, height_ratios=[1.6, 1.6, 1.6, 1.05, 1.25], hspace=0.62, wspace=0.20)
    values = spots["boundary_transition_index"].to_numpy(float)
    limit = float(np.nanpercentile(np.abs(values), 98))
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    cmap = LinearSegmentedColormap.from_list("bti", [BLUE, "#F7F8FA", GREEN])
    group_order = pd.Categorical(human["coarse_region"], categories=["control", "border", "fibrotic"], ordered=True)
    sample_order = list(human.assign(_group_order=group_order).sort_values(["_group_order", "patient"])["sample"])
    for idx, sample in enumerate(sample_order):
        ax = fig.add_subplot(grid[idx // 3, idx % 3])
        rows = spots.loc[spots["sample"] == sample]
        plot_spatial_map(ax, rows, norm, cmap)
        meta = human.loc[human["sample"] == sample].iloc[0]
        title = f"{str(meta['coarse_region']).capitalize()} | P{int(meta['patient'])} | n={int(meta['n_spots']):,}"
        ax.set_title(title, fontsize=6.3, fontweight="bold", color=GROUP_COLORS[str(meta["coarse_region"])], pad=2)
        if idx == 0:
            panel_label(ax, "a", x=-0.10, y=1.13)
    cax = fig.add_axes([0.30, 0.405, 0.40, 0.012])
    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    colorbar = fig.colorbar(scalar, cax=cax, orientation="horizontal")
    colorbar.set_label("Spot-level boundary transition index", fontsize=5.8)
    colorbar.ax.tick_params(labelsize=5.3, length=2)

    ax_qc = fig.add_subplot(grid[3, :])
    panel_label(ax_qc, "b", x=-0.035, y=1.08)
    module_cols = [
        column
        for column in human.columns
        if column.startswith("detected_genes_")
        and not any(token in column for token in ("baseline", "surrogate"))
    ]
    matrix = human.set_index("sample").loc[sample_order, module_cols].to_numpy(float)
    image = ax_qc.imshow(matrix, cmap=LinearSegmentedColormap.from_list("detect", ["#F3F5F8", PURPLE]), aspect="auto")
    ax_qc.set_yticks(range(len(sample_order)), [f"P{int(human.set_index('sample').loc[s, 'patient'])}" for s in sample_order])
    short_names = {
        "CM_BZ1_TRANSITION": "CM BZ1",
        "CM_BZ2_MECHANICAL_EDGE": "CM BZ2",
        "FAP_POSTN_PATHO_FIBROBLAST": "FAP/POSTN",
        "CTHRC1_REPARATIVE_CF": "CTHRC1",
        "MYOFIBROBLAST_CONTRACTILE": "Myofib.",
        "CCR2_IL1B_MYLOID": "CCR2/IL1B",
        "TGFB_SIGNALING": "TGF-beta",
        "ECM_REMODELING": "ECM",
        "ENDOTHELIAL_INFLAMMATORY_FIBROSIS": "Endothelial",
        "HYPOXIA_ISCHEMIA": "Hypoxia",
        "CM_ELECTRICAL_CALCIUM_FUNCTION": "CM function",
        "INJURY_STRESS_IMMEDIATE_EARLY": "Immediate early",
    }
    short = [short_names.get(column.replace("detected_genes_", ""), column.replace("detected_genes_", "")) for column in module_cols]
    ax_qc.set_xticks(range(len(short)), short, rotation=25, ha="right")
    ax_qc.tick_params(axis="x", labelsize=4.9)
    ax_qc.set_title("Transferred signature-gene detection by patient", loc="left", fontweight="bold", pad=4)
    for spine in ax_qc.spines.values():
        spine.set_visible(False)
    fig.colorbar(image, ax=ax_qc, fraction=0.018, pad=0.012, label="Detected genes")

    ax_loo = fig.add_subplot(grid[4, :2])
    panel_label(ax_loo, "c", x=-0.07, y=1.16)
    y = np.arange(len(benchmark_gain))[::-1]
    sensitivity_labels = {
        "Yamada prespecified tissue-mean mechanical": "Yamada whole-section mech.",
        "Yamada spatial mechanical-to-scar transition": "Yamada spatial transition",
        "Yamada day7 scar": "Yamada D7 scar repair",
        "Hernandez scar maturation": "Hernandez D5 scar repair",
        "Kuppe border-ischaemic mechanical": "Kuppe border mechanical",
        "Kuppe fibrotic scar": "Kuppe fibrotic scar repair",
    }
    labels = [sensitivity_labels.get(value, value) for value in benchmark_gain["comparison"]]
    ax_loo.scatter(benchmark_gain["baseline_loo_consistency"], y - 0.08, marker="s", color="#8B95A3", s=18, label="Simple baseline")
    ax_loo.scatter(benchmark_gain["framework_loo_consistency"], y + 0.08, marker="o", color=BLUE, s=20, label="Framework")
    ax_loo.set_yticks(y, labels)
    ax_loo.set_xlim(-0.03, 1.03)
    ax_loo.set_xlabel("Leave-one-unit-out direction consistency")
    ax_loo.set_title("Sensitivity to individual sections/patients", loc="left", fontweight="bold", pad=4)
    ax_loo.legend(loc="lower right", fontsize=5.0)

    ax_counts = fig.add_subplot(grid[4, 2])
    panel_label(ax_counts, "d", x=-0.18, y=1.16)
    grouped = human.groupby("coarse_region").agg(n_patients=("patient", "nunique"), tissue_spots=("n_spots", "sum")).reindex(["control", "border", "fibrotic"])
    bars = ax_counts.bar(range(3), grouped["tissue_spots"] / 1000.0, color=[GROUP_COLORS[group] for group in grouped.index], width=0.62)
    ax_counts.set_xticks(range(3), ["Control", "Border", "Fibrotic"], rotation=25, ha="right")
    ax_counts.set_ylabel("Tissue spots (thousands)")
    ax_counts.set_title("Human transfer coverage", loc="left", fontweight="bold", pad=4)
    for bar, patients in zip(bars, grouped["n_patients"], strict=True):
        ax_counts.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.25, f"n={int(patients)}", ha="center", va="bottom", fontsize=5.4, color=INK)
    fig.subplots_adjust(left=0.17, right=0.97, top=0.98, bottom=0.08)
    return fig


def main() -> None:
    setup_style()
    mouse = pd.read_csv(MOUSE_SECTION, sep="\t")
    area = pd.read_csv(AUTHOR_AREA, sep="\t")
    human = pd.read_csv(HUMAN_SECTION, sep="\t")
    spots = pd.read_csv(HUMAN_SPOT, sep="\t")
    benchmark = pd.read_csv(BENCHMARK, sep="\t")
    benchmark_gain = pd.read_csv(BENCHMARK_GAIN, sep="\t")
    transition = pd.read_csv(TRANSITION, sep="\t")

    main_figure = build_main_figure(mouse, area, human, benchmark, transition)
    save_figure(main_figure, FIG7)
    plt.close(main_figure)

    supplementary = build_supplementary_figure(human, spots, benchmark_gain)
    save_figure(supplementary, FIGS9)
    plt.close(supplementary)
    copy_source_data()
    print(f"wrote {FIG7}.[png|pdf|svg|tiff]")
    print(f"wrote {FIGS9}.[png|pdf|svg|tiff]")


if __name__ == "__main__":
    main()
