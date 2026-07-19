#!/usr/bin/env python3
"""Assemble upload-ready Supplementary Figures S1 and S2 from sample maps."""

from __future__ import annotations

import sys
from pathlib import Path

from project_paths import project_root


ROOT = project_root(__file__)
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))

import matplotlib.pyplot as plt  # noqa: E402


FIGURE_DIR = ROOT / "results" / "figures"
SAMPLES = ["d3_1", "d3_2", "d3_3", "d7_1", "d7_2", "d7_3"]
DISPLAY = {sample: sample.upper() for sample in SAMPLES}


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42


def save_outputs(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")


def assemble_signature_maps() -> None:
    fig, axes = plt.subplots(3, 2, figsize=(7.20, 5.05), facecolor="white")
    fig.subplots_adjust(left=0.025, right=0.99, top=0.91, bottom=0.025, hspace=0.16, wspace=0.04)
    for idx, (ax, sample) in enumerate(zip(axes.ravel(), SAMPLES, strict=True)):
        image = plt.imread(FIGURE_DIR / f"gse214611_{sample}_signature_score_maps.png")
        ax.imshow(image)
        ax.set_axis_off()
        ax.text(0.005, 1.02, chr(ord("a") + idx), transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")
        ax.text(0.055, 1.02, DISPLAY[sample], transform=ax.transAxes, fontsize=7.3, fontweight="bold", va="bottom")
    fig.suptitle("Supplementary Figure S1 | Sample-level signature maps", x=0.025, y=0.975, ha="left", fontsize=9.5, fontweight="bold")
    save_outputs(fig, FIGURE_DIR / "Supplementary_Figure_S1_sample_signature_maps")
    plt.close(fig)


def assemble_boundary_maps() -> None:
    fig, axes = plt.subplots(2, 3, figsize=(7.20, 5.20), facecolor="white")
    fig.subplots_adjust(left=0.025, right=0.99, top=0.91, bottom=0.025, hspace=0.14, wspace=0.08)
    for idx, (ax, sample) in enumerate(zip(axes.ravel(), SAMPLES, strict=True)):
        image = plt.imread(FIGURE_DIR / f"gse214611_{sample}_domain34_boundary_map.png")
        ax.imshow(image)
        ax.set_axis_off()
        ax.text(0.005, 1.02, chr(ord("a") + idx), transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")
        ax.text(0.085, 1.02, DISPLAY[sample], transform=ax.transAxes, fontsize=7.3, fontweight="bold", va="bottom")
    fig.suptitle("Supplementary Figure S2 | Domain 3/4 boundary maps by sample", x=0.025, y=0.975, ha="left", fontsize=9.5, fontweight="bold")
    save_outputs(fig, FIGURE_DIR / "Supplementary_Figure_S2_domain34_boundary_maps")
    plt.close(fig)


def main() -> None:
    assemble_signature_maps()
    assemble_boundary_maps()
    print("Wrote Supplementary Figures S1 and S2")


if __name__ == "__main__":
    main()
