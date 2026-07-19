#!/usr/bin/env python3
"""Plot final tutor-requested rigor upgrades for the spatial MI manuscript."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

from project_paths import project_root


ROOT = project_root(__file__)
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


TABLE_DIR = ROOT / "results" / "tables"
FIGURE_DIR = ROOT / "results" / "figures"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def safe_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def short_check(row: dict[str, str]) -> str:
    stage = "D3" if row["stage"] == "day3_mi" else "D7"
    score = {
        "mechanical_border": "Mech.",
        "fibroblast_scar_repair": "Scar",
        "immune_fibrotic_activation": "Immune",
    }[row["score"]]
    return f"{stage} {score}"


def method_label(name: str) -> str:
    return {
        "mean_expression_zsum": "Mean",
        "rank_auc_zsum": "Rank",
        "expression_matched_control_zsum": "Expr-matched",
        "overlap_reduced_mean_zsum": "No-overlap",
    }[name]


def main() -> None:
    qc = read_tsv(TABLE_DIR / "gse214611_data_provenance_qc.tsv")
    alt = read_tsv(TABLE_DIR / "gse214611_alternative_score_sensitivity_summary.tsv")
    random = read_tsv(TABLE_DIR / "gse214611_expression_matched_random_signature_controls.tsv")
    graph = read_tsv(TABLE_DIR / "gse214611_visium_graph_boundary_analysis.tsv")
    edges = read_tsv(TABLE_DIR / "gse214611_visium_graph_boundary_edge_gradients.tsv")
    hist = read_tsv(TABLE_DIR / "gse214611_he_image_score_alignment.tsv")

    mouse_qc = [row for row in qc if row["species"] == "Mus musculus"]
    alt_checks = ["D3 Mech.", "D7 Mech.", "D7 Scar", "D7 Immune"]
    methods = ["mean_expression_zsum", "rank_auc_zsum", "expression_matched_control_zsum", "overlap_reduced_mean_zsum"]
    palette = {
        "Mean": "#334155",
        "Rank": "#2563eb",
        "Expr-matched": "#0f766e",
        "No-overlap": "#d97706",
    }

    fig = plt.figure(figsize=(14.2, 10.2), dpi=220)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.05, 1.0], hspace=0.42, wspace=0.28)

    ax1 = fig.add_subplot(gs[0, 0])
    samples = [row["analysis_sample"] for row in mouse_qc]
    x = np.arange(len(samples))
    spots = np.array([safe_float(row["n_spots_in_analysis"]) for row in mouse_qc])
    median_genes = np.array([safe_float(row["median_genes_in_analysis_spots"]) for row in mouse_qc])
    ax1.bar(x - 0.18, spots / 1000.0, width=0.36, color="#64748b", label="spots (x1000)")
    ax1.bar(x + 0.18, median_genes / 1000.0, width=0.36, color="#94a3b8", label="median genes (x1000)")
    ax1.set_xticks(x, samples, rotation=0)
    ax1.set_ylabel("Count (thousands)")
    ax1.set_title("a  Provenance/QC by biological replicate", loc="left", fontweight="bold")
    ax1.legend(frameon=False, fontsize=8, ncols=2)
    ax1.spines[["top", "right"]].set_visible(False)

    ax2 = fig.add_subplot(gs[0, 1])
    width = 0.18
    for method_idx, method in enumerate(methods):
        values = []
        for check in alt_checks:
            match = [
                row
                for row in alt
                if method_label(row["scoring_method"]) == method_label(method) and short_check(row) == check
            ]
            values.append(safe_float(match[0]["mean_margin"]) if match else np.nan)
        offsets = (method_idx - 1.5) * width
        label = method_label(method)
        ax2.bar(np.arange(len(alt_checks)) + offsets, values, width=width, color=palette[label], label=label)
    ax2.axhline(0, color="#111827", linewidth=0.8)
    ax2.set_xticks(np.arange(len(alt_checks)), alt_checks)
    ax2.set_ylabel("Expected-direction margin")
    ax2.set_title("b  Alternative scoring preserves primary directions", loc="left", fontweight="bold")
    ax2.legend(frameon=False, fontsize=7, ncols=2)
    ax2.spines[["top", "right"]].set_visible(False)

    ax3 = fig.add_subplot(gs[1, 0])
    random_labels = [short_check(row) for row in random]
    obs = np.array([safe_float(row["observed_mean_margin"]) for row in random])
    p95 = np.array([safe_float(row["null_p95"]) for row in random])
    p99 = np.array([safe_float(row["null_p99"]) for row in random])
    rr = np.arange(len(random_labels))
    ax3.bar(rr - 0.22, obs, width=0.22, color="#0f766e", label="observed")
    ax3.bar(rr, p95, width=0.22, color="#bae6fd", label="matched null p95")
    ax3.bar(rr + 0.22, p99, width=0.22, color="#7dd3fc", label="matched null p99")
    ax3.set_xticks(rr, random_labels)
    ax3.set_ylabel("Stage-level margin")
    ax3.set_title("c  Expression-matched random-signature controls", loc="left", fontweight="bold")
    ax3.legend(frameon=False, fontsize=8)
    ax3.spines[["top", "right"]].set_visible(False)

    ax4 = fig.add_subplot(gs[1, 1])
    g_samples = [row["sample"] for row in graph]
    gx = np.arange(len(g_samples))
    d3 = np.array([safe_float(row["domain3_graph_boundary_fraction"]) for row in graph])
    d4 = np.array([safe_float(row["domain4_graph_boundary_fraction"]) for row in graph])
    ax4.plot(gx, d3, marker="o", color="#2563eb", label="domain 3")
    ax4.plot(gx, d4, marker="o", color="#059669", label="domain 4")
    ax4.set_xticks(gx, g_samples)
    ax4.set_ylim(0, 1.02)
    ax4.set_ylabel("Graph boundary fraction")
    ax4.set_title("d  Visium graph-based domain 3/4 contacts", loc="left", fontweight="bold")
    ax4.legend(frameon=False, fontsize=8)
    ax4.spines[["top", "right"]].set_visible(False)

    ax5 = fig.add_subplot(gs[2, 0])
    edge_scores = ["mechanical_border", "fibroblast_scar_repair", "immune_fibrotic_activation"]
    edge_labels = ["Mech.", "Scar", "Immune"]
    edge_colors = ["#2563eb", "#059669", "#d97706"]
    stage_positions = {"day3_mi": 0, "day7_mi": 1}
    for score, label, color in zip(edge_scores, edge_labels, edge_colors, strict=True):
        for stage in ["day3_mi", "day7_mi"]:
            vals = [
                safe_float(row["mean_edge_delta"])
                for row in edges
                if row["score"] == score and row["stage"] == stage
            ]
            xpos = stage_positions[stage] + {"mechanical_border": -0.18, "fibroblast_scar_repair": 0, "immune_fibrotic_activation": 0.18}[score]
            ax5.scatter([xpos] * len(vals), vals, color=color, s=32, alpha=0.85)
            ax5.plot([xpos - 0.06, xpos + 0.06], [np.mean(vals), np.mean(vals)], color=color, linewidth=2.5, label=label if stage == "day3_mi" else None)
    ax5.axhline(0, color="#111827", linewidth=0.8)
    ax5.set_xticks([0, 1], ["D3", "D7"])
    ax5.set_ylabel("Domain 4 minus domain 3 edge delta")
    ax5.set_title("e  Graph-edge score gradients", loc="left", fontweight="bold")
    ax5.legend(frameon=False, fontsize=8, ncols=3)
    ax5.spines[["top", "right"]].set_visible(False)

    ax6 = fig.add_subplot(gs[2, 1])
    selected_hist = [
        row
        for row in hist
        if row["histology_feature"] in {"he_brightness", "he_saturation"}
        and row["score"] in {"mechanical_border", "fibroblast_scar_repair"}
    ]
    labels = []
    values = []
    colors = []
    for feature in ["he_brightness", "he_saturation"]:
        for score in ["mechanical_border", "fibroblast_scar_repair"]:
            vals = [safe_float(row["pearson_r"]) for row in selected_hist if row["histology_feature"] == feature and row["score"] == score]
            labels.append(f"{feature.replace('he_', '')}\n{score.split('_')[0]}")
            values.append(float(np.nanmean(vals)))
            colors.append("#2563eb" if score == "mechanical_border" else "#059669")
    ax6.bar(np.arange(len(values)), values, color=colors, alpha=0.85)
    ax6.axhline(0, color="#111827", linewidth=0.8)
    ax6.set_xticks(np.arange(len(values)), labels)
    ax6.set_ylabel("Mean Pearson r across samples")
    ax6.set_title("f  H&E image-intensity audit (non-pathologist)", loc="left", fontweight="bold")
    ax6.spines[["top", "right"]].set_visible(False)

    fig.suptitle(
        "Final reviewer-requested rigor checks: provenance, scoring, graph boundary and image audit",
        x=0.02,
        y=0.985,
        ha="left",
        fontsize=15,
        fontweight="bold",
    )
    fig.text(
        0.02,
        0.015,
        "All analyses use public GSE214611 Visium data. Spot-level statistics are used for visualization and sensitivity analysis only.",
        fontsize=8,
        color="#475569",
    )

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    for suffix in ["png", "pdf", "svg", "tiff"]:
        save_kwargs = {
            "bbox_inches": "tight",
            "dpi": 450 if suffix in {"png", "tiff"} else None,
        }
        if suffix == "tiff":
            save_kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(FIGURE_DIR / f"cvr_final_reviewer_upgrade_checks.{suffix}", **save_kwargs)
    plt.close(fig)
    print(f"Wrote {FIGURE_DIR / 'cvr_final_reviewer_upgrade_checks.png'}")


if __name__ == "__main__":
    main()
