#!/usr/bin/env python3
"""Interpret GSE214611 spatial domains from author-provided metadata."""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SPATIAL_METADATA = (
    ROOT
    / "data"
    / "metadata"
    / "gse214611_author"
    / "spatial_object_integrated_metadata.csv"
)
OUT_DIR = ROOT / "results" / "tables"

RAW_SCORE_COLUMNS = [
    "BZ1_genes",
    "BZ2_genes",
    "RZ_genes",
    "X33_Fib_Postn",
    "X36_Fib_Rep",
    "X11_Monos",
    "X14_Mac_C1qa",
    "CM_Score",
]


def safe_float(value: str) -> float:
    try:
        return float(value) if value != "" else 0.0
    except ValueError:
        return 0.0


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def stage_from_sample(sample: str) -> str:
    if sample.startswith("sham"):
        return "sham"
    if sample in {"1hr", "4hr"}:
        return "acute_hours"
    if sample.startswith("D3"):
        return "day3_mi"
    if sample.startswith("D7"):
        return "day7_mi"
    if sample == "I/R":
        return "ischemia_reperfusion"
    if sample == "TAC":
        return "pressure_overload"
    if sample == "ISO":
        return "isoproterenol"
    if sample.startswith("NP"):
        return "needle_pass"
    return "other"


def zscore_columns(rows: list[dict[str, str]], columns: list[str]) -> dict[str, list[float]]:
    values_by_col: dict[str, list[float]] = {col: [] for col in columns}
    for row in rows:
        for col in columns:
            values_by_col[col].append(safe_float(row.get(col, "")))

    zscores: dict[str, list[float]] = {}
    for col, values in values_by_col.items():
        mean = sum(values) / len(values)
        sd = math.sqrt(sum((value - mean) ** 2 for value in values) / len(values)) or 1.0
        zscores[col] = [(value - mean) / sd for value in values]
    return zscores


def add_prototype_scores(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    z = zscore_columns(
        rows,
        [
            "BZ2_genes",
            "BZ1_genes",
            "RZ_genes",
            "X33_Fib_Postn",
            "X36_Fib_Rep",
            "X11_Monos",
            "X14_Mac_C1qa",
        ],
    )

    scored = []
    for idx, row in enumerate(rows):
        out = dict(row)
        immune_proxy = 0.5 * z["X11_Monos"][idx] + 0.5 * z["X14_Mac_C1qa"][idx]
        fibrotic_risk = (
            z["BZ2_genes"][idx]
            + z["X33_Fib_Postn"][idx]
            + immune_proxy
            - 0.5 * z["RZ_genes"][idx]
        )
        border_activation = 0.5 * z["BZ1_genes"][idx] + 0.5 * z["BZ2_genes"][idx]
        out["stage"] = stage_from_sample(row["orig.ident"])
        out["prototype_immune_proxy"] = f"{immune_proxy:.8g}"
        out["prototype_fibrotic_risk"] = f"{fibrotic_risk:.8g}"
        out["prototype_border_activation"] = f"{border_activation:.8g}"
        out["prototype_repair_proxy"] = f"{z['X36_Fib_Rep'][idx]:.8g}"
        scored.append(out)
    return scored


def summarize(rows: list[dict[str, str]], group_columns: list[str], out_path: Path) -> None:
    score_columns = [
        *RAW_SCORE_COLUMNS,
        "prototype_fibrotic_risk",
        "prototype_repair_proxy",
        "prototype_immune_proxy",
        "prototype_border_activation",
    ]
    sums: dict[tuple[str, ...], dict[str, float]] = defaultdict(lambda: defaultdict(float))
    counts: dict[tuple[str, ...], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    group_n: Counter[tuple[str, ...]] = Counter()

    for row in rows:
        key = tuple(row[col] for col in group_columns)
        group_n[key] += 1
        for col in score_columns:
            sums[key][col] += safe_float(row[col])
            counts[key][col] += 1

    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([*group_columns, "n_spots", *[f"mean_{col}" for col in score_columns]])
        for key in sorted(group_n):
            values = []
            for col in score_columns:
                values.append(f"{sums[key][col] / counts[key][col]:.6g}")
            writer.writerow([*key, group_n[key], *values])


def summarize_domain_rank(rows: list[dict[str, str]], out_path: Path) -> None:
    summarize_by_domain: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    counts: Counter[str] = Counter()
    score_columns = [
        "RZ_genes",
        "BZ1_genes",
        "BZ2_genes",
        "X33_Fib_Postn",
        "X36_Fib_Rep",
        "prototype_fibrotic_risk",
        "prototype_border_activation",
    ]

    for row in rows:
        domain = row["annotated"]
        counts[domain] += 1
        for col in score_columns:
            summarize_by_domain[domain][col] += safe_float(row[col])

    domain_means = {
        domain: {col: summarize_by_domain[domain][col] / counts[domain] for col in score_columns}
        for domain in counts
    }

    def rank_domains(col: str, reverse: bool = True) -> dict[str, int]:
        ordered = sorted(domain_means, key=lambda domain: domain_means[domain][col], reverse=reverse)
        return {domain: idx + 1 for idx, domain in enumerate(ordered)}

    ranks = {
        "rz_rank": rank_domains("RZ_genes"),
        "bz1_rank": rank_domains("BZ1_genes"),
        "bz2_rank": rank_domains("BZ2_genes"),
        "postn_rank": rank_domains("X33_Fib_Postn"),
        "repair_rank": rank_domains("X36_Fib_Rep"),
        "risk_rank": rank_domains("prototype_fibrotic_risk"),
        "border_rank": rank_domains("prototype_border_activation"),
    }

    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "annotated",
                "n_spots",
                *[f"mean_{col}" for col in score_columns],
                *ranks.keys(),
                "working_interpretation",
            ]
        )
        for domain in sorted(domain_means):
            means = [f"{domain_means[domain][col]:.6g}" for col in score_columns]
            domain_ranks = [ranks[name][domain] for name in ranks]
            writer.writerow(
                [
                    domain,
                    counts[domain],
                    *means,
                    *domain_ranks,
                    working_interpretation(domain, ranks),
                ]
            )


def working_interpretation(domain: str, ranks: dict[str, dict[str, int]]) -> str:
    if ranks["risk_rank"][domain] == 1 and ranks["repair_rank"][domain] == 1:
        return "fibrotic_reparative_hotspot"
    if ranks["border_rank"][domain] == 1 or ranks["bz2_rank"][domain] == 1:
        return "mechanical_border_zone_like"
    if ranks["rz_rank"][domain] == 1:
        return "remote_cardiomyocyte_like"
    if ranks["postn_rank"][domain] <= 2:
        return "activated_fibroblast_enriched"
    return "mixed_or_transition"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = add_prototype_scores(read_rows(SPATIAL_METADATA))
    summarize(rows, ["annotated"], OUT_DIR / "gse214611_domain_summary.tsv")
    summarize(rows, ["stage", "annotated"], OUT_DIR / "gse214611_domain_summary_by_stage.tsv")
    summarize(rows, ["orig.ident", "annotated"], OUT_DIR / "gse214611_domain_summary_by_sample.tsv")
    summarize_domain_rank(rows, OUT_DIR / "gse214611_domain_interpretation.tsv")


if __name__ == "__main__":
    main()
