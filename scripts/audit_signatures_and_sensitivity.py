#!/usr/bin/env python3
"""Signature transparency and sensitivity audits for the spatial MI model."""

from __future__ import annotations

import csv
import itertools
import math
from collections import defaultdict
from pathlib import Path

from project_paths import project_root

ROOT = project_root(__file__)
TABLE_DIR = ROOT / "results" / "tables"
CONFIG_DIR = Path(__file__).resolve().parents[1] / "config"
SPOT = TABLE_DIR / "gse214611_d3_d7_signature_scores_by_spot.tsv"
SIGNATURES = CONFIG_DIR / "spatial_cardiac_border_zone_signatures.tsv"

OUTPUT_COMPONENTS = {
    "mechanical_border": ["z_CM_BZ1_TRANSITION", "z_CM_BZ2_MECHANICAL_EDGE"],
    "immune_fibrotic_activation": ["z_CCR2_IL1B_MYLOID", "z_TGFB_SIGNALING", "z_FAP_POSTN_PATHO_FIBROBLAST"],
    "fibroblast_scar_repair": [
        "z_FAP_POSTN_PATHO_FIBROBLAST",
        "z_ECM_REMODELING",
        "z_CTHRC1_REPARATIVE_CF",
        "z_MYOFIBROBLAST_CONTRACTILE",
    ],
}

KEY_GENES = ["POSTN", "COL1A1", "COL1A2", "CTHRC1", "NPPA", "NPPB", "XIRP2", "SPP1", "LGALS3"]

OVERLAP_FIELDS = [
    "left_signature",
    "right_signature",
    "n_left_genes",
    "n_right_genes",
    "n_overlap_genes",
    "jaccard",
    "overlap_genes",
]
KEY_GENE_FIELDS = [
    "key_gene",
    "dropped_component",
    "output",
    "stage",
    "n_samples",
    "mean_expected_direction_margin",
    "min_expected_direction_margin",
    "direction_preserved",
]
COLLINEARITY_FIELDS = [
    "left_output",
    "right_output",
    "n_spots",
    "pearson_r",
    "spearman_r",
    "interpretation",
]


def mouse_gene_upper(gene: str) -> str:
    return gene.strip().upper()


def parse_gene_set(value: str) -> set[str]:
    return {mouse_gene_upper(gene) for gene in value.split(";") if gene.strip()}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def safe_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def write_tsv_lf(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    if not rows:
        raise ValueError(f"No rows to write for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def finite_values(values: list[float]) -> list[float]:
    return [value for value in values if math.isfinite(value)]


def mean_finite(values: list[float]) -> float:
    filtered = finite_values(values)
    if not filtered:
        return float("nan")
    return sum(filtered) / len(filtered)


def min_finite(values: list[float]) -> float:
    filtered = finite_values(values)
    if not filtered:
        return float("nan")
    return min(filtered)


def pearson_r(left: list[float], right: list[float]) -> float:
    paired = [(x, y) for x, y in zip(left, right) if math.isfinite(x) and math.isfinite(y)]
    if len(paired) < 2:
        return float("nan")
    xs = [x for x, _ in paired]
    ys = [y for _, y in paired]
    x_mean = sum(xs) / len(xs)
    y_mean = sum(ys) / len(ys)
    x_deltas = [x - x_mean for x in xs]
    y_deltas = [y - y_mean for y in ys]
    x_ss = sum(delta * delta for delta in x_deltas)
    y_ss = sum(delta * delta for delta in y_deltas)
    if x_ss == 0 or y_ss == 0:
        return float("nan")
    covariance = sum(x_delta * y_delta for x_delta, y_delta in zip(x_deltas, y_deltas))
    return covariance / math.sqrt(x_ss * y_ss)


def average_ranks(values: list[float]) -> list[float]:
    ranked = [0.0] * len(values)
    ordered = sorted((value, idx) for idx, value in enumerate(values))
    position = 0
    while position < len(ordered):
        tie_end = position + 1
        while tie_end < len(ordered) and ordered[tie_end][0] == ordered[position][0]:
            tie_end += 1
        rank = (position + 1 + tie_end) / 2
        for _, original_idx in ordered[position:tie_end]:
            ranked[original_idx] = rank
        position = tie_end
    return ranked


def spearman_r(left: list[float], right: list[float]) -> float:
    paired = [(x, y) for x, y in zip(left, right) if math.isfinite(x) and math.isfinite(y)]
    if len(paired) < 2:
        return float("nan")
    xs = [x for x, _ in paired]
    ys = [y for _, y in paired]
    return pearson_r(average_ranks(xs), average_ranks(ys))


def write_overlap(signatures: list[dict[str, str]]) -> None:
    module_genes = {
        row["signature_id"]: parse_gene_set(row["genes_human"]) | parse_gene_set(row["genes_mouse"]) for row in signatures
    }
    rows = []
    for left, right in itertools.combinations(sorted(module_genes), 2):
        a = module_genes[left]
        b = module_genes[right]
        overlap = a & b
        union = a | b
        rows.append(
            {
                "left_signature": left,
                "right_signature": right,
                "n_left_genes": len(a),
                "n_right_genes": len(b),
                "n_overlap_genes": len(overlap),
                "jaccard": f"{len(overlap) / len(union):.8g}" if union else "NA",
                "overlap_genes": ";".join(sorted(overlap)) if overlap else "NA",
            }
        )
    write_tsv_lf(TABLE_DIR / "gse214611_signature_overlap_audit.tsv", rows, OVERLAP_FIELDS)


def output_score(row: dict[str, str], components: list[str], excluded_component: str | None = None) -> float:
    active = [component for component in components if component != excluded_component]
    return float(sum(safe_float(row[component]) for component in active))


def write_key_gene_sensitivity(signatures: list[dict[str, str]], spot_rows: list[dict[str, str]]) -> None:
    gene_to_modules: dict[str, set[str]] = defaultdict(set)
    for sig in signatures:
        genes = parse_gene_set(sig["genes_human"]) | parse_gene_set(sig["genes_mouse"])
        for gene in genes:
            gene_to_modules[gene].add("z_" + sig["signature_id"])

    rows = []
    filtered = [r for r in spot_rows if r["stage"] == "day7_mi" and r["annotated"] in {"3", "4"}]
    samples = sorted({r["sample"] for r in filtered})
    for gene in KEY_GENES:
        dropped_components = sorted(gene_to_modules.get(gene, set()))
        for output_name, components in OUTPUT_COMPONENTS.items():
            active_drops = [component for component in dropped_components if component in components]
            if not active_drops:
                continue
            for drop_component in active_drops:
                by_sample_domain: dict[tuple[str, str], list[float]] = defaultdict(list)
                for row in filtered:
                    by_sample_domain[(row["sample"], row["annotated"])].append(
                        output_score(row, components, excluded_component=drop_component)
                    )
                margins = []
                for sample in samples:
                    if not by_sample_domain[(sample, "3")] or not by_sample_domain[(sample, "4")]:
                        continue
                    d3 = mean_finite(by_sample_domain[(sample, "3")])
                    d4 = mean_finite(by_sample_domain[(sample, "4")])
                    margin = d3 - d4 if output_name == "mechanical_border" else d4 - d3
                    if math.isfinite(margin):
                        margins.append(float(margin))
                if not margins:
                    continue
                min_margin = min_finite(margins)
                rows.append(
                    {
                        "key_gene": gene,
                        "dropped_component": drop_component,
                        "output": output_name,
                        "stage": "day7_mi",
                        "n_samples": len(margins),
                        "mean_expected_direction_margin": f"{mean_finite(margins):.8g}",
                        "min_expected_direction_margin": f"{min_margin:.8g}",
                        "direction_preserved": str(bool(min_margin > 0)),
                    }
                )
    write_tsv_lf(TABLE_DIR / "gse214611_key_gene_removal_sensitivity.tsv", rows, KEY_GENE_FIELDS)


def write_collinearity(spot_rows: list[dict[str, str]]) -> None:
    score_values = {}
    for output_name, components in OUTPUT_COMPONENTS.items():
        score_values[output_name] = [output_score(row, components) for row in spot_rows]
    rows = []
    for left, right in itertools.combinations(score_values, 2):
        x = score_values[left]
        y = score_values[right]
        pearson = pearson_r(x, y)
        spearman = spearman_r(x, y)
        rows.append(
            {
                "left_output": left,
                "right_output": right,
                "n_spots": len(x),
                "pearson_r": f"{pearson:.8g}",
                "spearman_r": f"{spearman:.8g}",
                "interpretation": "correlated_not_identical" if abs(spearman) < 0.85 else "highly_collinear",
            }
        )
    write_tsv_lf(TABLE_DIR / "gse214611_score_collinearity_audit.tsv", rows, COLLINEARITY_FIELDS)


def main() -> None:
    signatures = read_tsv(SIGNATURES)
    spot_rows = read_tsv(SPOT)
    write_overlap(signatures)
    write_key_gene_sensitivity(signatures, spot_rows)
    write_collinearity(spot_rows)
    print("Wrote signature audit and sensitivity tables")


if __name__ == "__main__":
    main()
