#!/usr/bin/env python3
"""Build a first-pass GSE214611 border-zone fibrotic risk proxy.

This prototype intentionally uses only author-provided spatial metadata columns.
It is a quick gate before downloading expression matrices and running full
gene-set scoring.
"""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

from project_paths import project_root


ROOT = project_root(__file__)
SPATIAL_METADATA = (
    ROOT
    / "data"
    / "metadata"
    / "gse214611_author"
    / "spatial_object_integrated_metadata.csv"
)
OUT_DIR = ROOT / "results" / "tables"

MODEL_COLUMNS = [
    "BZ2_genes",
    "BZ1_genes",
    "RZ_genes",
    "X33_Fib_Postn",
    "X36_Fib_Rep",
    "X11_Monos",
    "X14_Mac_C1qa",
]


def safe_float(value: str) -> float | None:
    try:
        if value == "":
            return None
        return float(value)
    except ValueError:
        return None


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def zscore_columns(rows: list[dict[str, str]], columns: list[str]) -> dict[str, list[float]]:
    raw: dict[str, list[float]] = {col: [] for col in columns}
    for row in rows:
        for col in columns:
            raw[col].append(safe_float(row.get(col, "")) or 0.0)

    z: dict[str, list[float]] = {}
    for col, values in raw.items():
        mean = sum(values) / len(values)
        variance = sum((value - mean) ** 2 for value in values) / len(values)
        sd = math.sqrt(variance) or 1.0
        z[col] = [(value - mean) / sd for value in values]
    return z


def add_scores(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    z = zscore_columns(rows, MODEL_COLUMNS)
    scored_rows: list[dict[str, str]] = []

    for idx, row in enumerate(rows):
        bz2 = z["BZ2_genes"][idx]
        bz1 = z["BZ1_genes"][idx]
        rz = z["RZ_genes"][idx]
        postn = z["X33_Fib_Postn"][idx]
        repair = z["X36_Fib_Rep"][idx]
        monos = z["X11_Monos"][idx]
        mac_c1qa = z["X14_Mac_C1qa"][idx]

        immune_proxy = 0.5 * monos + 0.5 * mac_c1qa
        fibrotic_risk = bz2 + postn + immune_proxy - 0.5 * rz
        repair_balance = repair - postn
        border_activation = 0.5 * bz1 + 0.5 * bz2

        out = dict(row)
        out["prototype_immune_proxy"] = f"{immune_proxy:.6g}"
        out["prototype_fibrotic_risk"] = f"{fibrotic_risk:.6g}"
        out["prototype_repair_proxy"] = f"{repair:.6g}"
        out["prototype_repair_minus_postn"] = f"{repair_balance:.6g}"
        out["prototype_border_activation"] = f"{border_activation:.6g}"
        scored_rows.append(out)

    return scored_rows


def summarize(
    rows: list[dict[str, str]],
    group_columns: list[str],
    score_columns: list[str],
    out_path: Path,
) -> None:
    sums: dict[tuple[str, ...], dict[str, float]] = defaultdict(lambda: defaultdict(float))
    counts: dict[tuple[str, ...], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    group_n: Counter[tuple[str, ...]] = Counter()

    for row in rows:
        key = tuple(row.get(col, "") for col in group_columns)
        group_n[key] += 1
        for col in score_columns:
            value = safe_float(row[col])
            if value is None:
                continue
            sums[key][col] += value
            counts[key][col] += 1

    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([*group_columns, "n_spots", *[f"mean_{col}" for col in score_columns]])
        for key, n in group_n.most_common():
            means = []
            for col in score_columns:
                denom = counts[key][col]
                means.append("" if denom == 0 else f"{sums[key][col] / denom:.6g}")
            writer.writerow([*key, n, *means])


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = read_rows(SPATIAL_METADATA)
    scored_rows = add_scores(rows)

    score_columns = [
        "prototype_fibrotic_risk",
        "prototype_repair_proxy",
        "prototype_repair_minus_postn",
        "prototype_immune_proxy",
        "prototype_border_activation",
    ]

    summarize(
        scored_rows,
        ["orig.ident"],
        score_columns,
        OUT_DIR / "gse214611_prototype_risk_by_sample.tsv",
    )
    summarize(
        scored_rows,
        ["orig.ident", "annotated"],
        score_columns,
        OUT_DIR / "gse214611_prototype_risk_by_sample_annotated.tsv",
    )


if __name__ == "__main__":
    main()
