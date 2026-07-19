#!/usr/bin/env python3
"""Transfer frozen spatial-state signatures to independent mouse MI studies."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from project_paths import project_root

import numpy as np

from external_spatial_validation_utils import (
    boundary_transition_index,
    composite_framework_scores,
    read_10x_mtx_bundle,
    read_manifest,
    read_tissue_positions,
    raw_framework_scores,
    score_fixed_signatures,
)
from score_gse214611_visium_signatures import (
    read_signatures,
)
from therapeutic_prioritization_utils import (
    build_knn_edges,
    graph_morans_i,
    subset_spots,
)


ROOT = project_root(__file__)
MANIFEST = (
    Path(__file__).resolve().parents[1]
    / "config"
    / "external_validation_manifest.tsv"
)
SIGNATURES = Path(__file__).resolve().parents[1] / "config" / "spatial_cardiac_border_zone_signatures.tsv"
DATA_ROOT = ROOT / "data" / "raw" / "external_validation"
OUT_ROOT = ROOT / "results" / "tables" / "external_validation"


BASELINES = {
    "source_bz_baseline": ["Nppa", "Xirp2", "Flnc"],
    "generic_fibrosis_baseline": ["Col1a1", "Col1a2", "Col3a1", "Fn1", "Postn"],
    "cardiomyocyte_composition_surrogate": ["Tnnt2", "Tnni3", "Myh6", "Actc1"],
    "fibroblast_composition_surrogate": ["Pdgfra", "Dcn", "Lum", "Col1a1"],
}


def local_path(row: dict[str, str]) -> Path:
    if row["file_role"] == "h5ad":
        name = f"{row['sample']}.h5ad"
    else:
        name = Path(row["url"]).name
    return DATA_ROOT / row["dataset"] / row["sample"] / name


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"no rows to write: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_author_areas(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    result: dict[str, str] = {}
    for row in rows:
        barcode = row.get("") or row.get("barcode") or row.get("Barcode")
        area = row.get("new_area", "")
        if barcode and area:
            result[barcode] = area
    return result


def summarize_section(
    metadata: dict[str, str],
    scores: dict[str, np.ndarray],
    detected: dict[str, int],
    coords: np.ndarray,
) -> dict[str, object]:
    edges = build_knn_edges(coords, k=6)
    row: dict[str, object] = {
        "dataset": metadata["dataset"],
        "study_role": metadata["study_role"],
        "sample": metadata["sample"],
        "stage_region": metadata["stage_region"],
        "n_spots": len(coords),
        "n_graph_edges": len(edges),
    }
    for name, values in scores.items():
        finite = values[np.isfinite(values)]
        row[f"mean_{name}"] = float(np.mean(finite)) if len(finite) else np.nan
        row[f"median_{name}"] = float(np.median(finite)) if len(finite) else np.nan
        row[f"p90_{name}"] = float(np.percentile(finite, 90)) if len(finite) else np.nan
        if name in {
            "mechanical_border_score",
            "immune_fibrotic_activation_score",
            "fibroblast_scar_repair_score",
            "boundary_transition_index",
            "source_bz_baseline",
            "generic_fibrosis_baseline",
        }:
            row[f"graph_morans_i_{name}"] = graph_morans_i(values, edges)
    for module, count in detected.items():
        row[f"detected_genes_{module}"] = count
    return row


def main() -> None:
    manifest = read_manifest(MANIFEST)
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in manifest:
        if row["species"] == "mouse":
            grouped[(row["dataset"], row["sample"])].append(row)

    signatures = read_signatures(SIGNATURES)
    all_signatures = {**signatures, **BASELINES}
    spot_rows: list[dict[str, object]] = []
    section_rows: list[dict[str, object]] = []
    area_rows: list[dict[str, object]] = []

    for (_, _), rows in sorted(grouped.items()):
        metadata = rows[0]
        by_role = {row["file_role"]: local_path(row) for row in rows}
        missing = [path for path in by_role.values() if not path.exists()]
        if missing:
            raise FileNotFoundError(f"external files missing for {metadata['sample']}: {missing}")
        matrix = read_10x_mtx_bundle(by_role["matrix"], by_role["features"], by_role["barcodes"])
        positions = read_tissue_positions(by_role["positions"])
        author_areas = read_author_areas(by_role.get("author_metadata"))
        selected = [idx for idx, barcode in enumerate(matrix.barcodes) if positions.get(barcode, {}).get("in_tissue") == 1]
        matrix = subset_spots(matrix, selected)
        coords = np.asarray(
            [
                [positions[barcode]["array_col"], positions[barcode]["array_row"]]
                for barcode in matrix.barcodes
            ],
            dtype=float,
        )
        modules, detected = score_fixed_signatures(matrix, all_signatures, normalize=True)
        framework = composite_framework_scores(modules)
        raw = raw_framework_scores(modules)
        bti = boundary_transition_index(
            framework["mechanical_border_score"],
            framework["immune_fibrotic_activation_score"],
            framework["fibroblast_scar_repair_score"],
        )
        combined = {
            **modules,
            **framework,
            **raw,
            "boundary_transition_index": bti,
        }
        section_rows.append(summarize_section(metadata, combined, detected, coords))
        if author_areas:
            barcode_to_index = {barcode: idx for idx, barcode in enumerate(matrix.barcodes)}
            for area in ("RZ", "BZ1", "BZ2", "IZ"):
                area_indices = [
                    barcode_to_index[barcode]
                    for barcode, label in author_areas.items()
                    if label == area and barcode in barcode_to_index
                ]
                if not area_indices:
                    continue
                area_row: dict[str, object] = {
                    "dataset": metadata["dataset"],
                    "sample": metadata["sample"],
                    "stage_region": metadata["stage_region"],
                    "author_area": area,
                    "n_spots": len(area_indices),
                }
                for name, values in combined.items():
                    selected_values = values[np.asarray(area_indices, dtype=int)]
                    area_row[f"mean_{name}"] = float(np.nanmean(selected_values))
                    area_row[f"median_{name}"] = float(np.nanmedian(selected_values))
                area_rows.append(area_row)
        for spot_idx, barcode in enumerate(matrix.barcodes):
            position = positions[barcode]
            row: dict[str, object] = {
                "dataset": metadata["dataset"],
                "study_role": metadata["study_role"],
                "sample": metadata["sample"],
                "stage_region": metadata["stage_region"],
                "barcode": barcode,
                "array_row": position["array_row"],
                "array_col": position["array_col"],
                "pxl_row_in_fullres": position["pxl_row_in_fullres"],
                "pxl_col_in_fullres": position["pxl_col_in_fullres"],
                "author_area": author_areas.get(barcode, ""),
            }
            for name, values in combined.items():
                row[name] = float(values[spot_idx])
            spot_rows.append(row)
        print(
            f"[scored] {metadata['dataset']} {metadata['sample']} "
            f"spots={len(matrix.barcodes)} mechanical_genes="
            f"{detected['CM_BZ2_MECHANICAL_EDGE']} scar_genes={detected['CTHRC1_REPARATIVE_CF']}"
        )

    write_tsv(OUT_ROOT / "external_mouse_spatial_scores_by_spot.tsv", spot_rows)
    write_tsv(OUT_ROOT / "external_mouse_section_summary.tsv", section_rows)
    write_tsv(OUT_ROOT / "external_mouse_author_area_summary.tsv", area_rows)
    print(f"sections={len(section_rows)} author_areas={len(area_rows)} spots={len(spot_rows)}")


if __name__ == "__main__":
    main()
