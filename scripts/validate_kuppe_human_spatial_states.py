#!/usr/bin/env python3
"""Transfer frozen human spatial-state signatures to the Kuppe MI atlas."""

from __future__ import annotations

import csv
from pathlib import Path

from project_paths import project_root

import numpy as np

from external_spatial_validation_utils import (
    boundary_transition_index,
    composite_framework_scores,
    load_h5ad_spatial,
    raw_framework_scores,
    read_manifest,
    score_fixed_signatures,
    subset_h5ad_spatial,
)
from therapeutic_prioritization_utils import (
    build_knn_edges,
    graph_morans_i,
)


ROOT = project_root(__file__)
PRIMARY_MANIFEST = (
    Path(__file__).resolve().parents[1]
    / "config"
    / "kuppe_primary_cellxgene_manifest.tsv"
)
SIGNATURES = Path(__file__).resolve().parents[1] / "config" / "spatial_cardiac_border_zone_signatures.tsv"
DATA_ROOT = ROOT / "data" / "raw" / "external_validation" / "KUPPE2022"
OUT_ROOT = ROOT / "results" / "tables" / "external_validation"


PRIMARY_SAMPLES = (
    "Visium_control_P1",
    "Visium_control_P8",
    "Visium_control_P17",
    "Visium_RZ_BZ_P3",
    "Visium_RZ_BZ_P12",
    "Visium_IZ_BZ_P2",
    "Visium_FZ_P14",
    "Visium_FZ_P18",
    "Visium_FZ_P20",
)


BASELINES = {
    "source_bz_baseline": ["NPPA", "XIRP2", "FLNC"],
    "generic_fibrosis_baseline": ["COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN"],
    "cardiomyocyte_composition_surrogate": ["TNNT2", "TNNI3", "MYH6", "ACTC1"],
    "fibroblast_composition_surrogate": ["PDGFRA", "DCN", "LUM", "COL1A1"],
}


def coarse_region(stage_region: str) -> str:
    value = stage_region.upper()
    if "CONTROL" in value:
        return "control"
    if "BZ" in value:
        return "border"
    if "FZ" in value:
        return "fibrotic"
    if "IZ" in value:
        return "ischaemic"
    if "RZ" in value:
        return "remote"
    return "other"


def author_label_columns(obs_columns: dict[str, list[str]]) -> list[str]:
    candidates = []
    for name in obs_columns:
        lower = name.lower()
        if any(token in lower for token in ("region", "cluster", "niche", "annotation", "tissue", "spot")):
            candidates.append(name)
    return sorted(candidates)[:8]


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"no rows to write: {path}")
    fields = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    manifest = read_manifest(PRIMARY_MANIFEST)
    metadata = {
        row["sample"]: row
        for row in manifest
        if row["dataset"] == "KUPPE2022" and row["file_role"] == "h5ad"
    }
    missing_manifest = [sample for sample in PRIMARY_SAMPLES if sample not in metadata]
    if missing_manifest:
        raise ValueError(f"primary Kuppe samples absent from manifest: {missing_manifest}")

    human_signatures: dict[str, list[str]] = {}
    with SIGNATURES.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            human_signatures[row["signature_id"]] = [
                gene.strip() for gene in row["genes_human"].split(";") if gene.strip()
            ]
    all_signatures = {**human_signatures, **BASELINES}

    spot_rows: list[dict[str, object]] = []
    section_rows: list[dict[str, object]] = []
    for sample in PRIMARY_SAMPLES:
        meta = metadata[sample]
        path = DATA_ROOT / sample / f"{sample}.h5ad"
        if not path.exists():
            raise FileNotFoundError(path)
        loaded = load_h5ad_spatial(path)
        n_array_spots = loaded.matrix.shape[1]
        in_tissue = loaded.obs_columns.get("in_tissue")
        if in_tissue is not None:
            tissue_indices = [idx for idx, value in enumerate(in_tissue) if value in {"1", "True", "true"}]
            if not tissue_indices:
                raise ValueError(f"{sample} contains an in_tissue field but no tissue spots")
            loaded = subset_h5ad_spatial(loaded, tissue_indices)
        normalize = loaded.matrix_scale == "counts"
        modules, detected = score_fixed_signatures(loaded.matrix, all_signatures, normalize=normalize)
        framework = composite_framework_scores(modules)
        raw = raw_framework_scores(modules)
        bti = boundary_transition_index(
            framework["mechanical_border_score"],
            framework["immune_fibrotic_activation_score"],
            framework["fibroblast_scar_repair_score"],
        )
        combined = {**modules, **framework, **raw, "boundary_transition_index": bti}
        edges = build_knn_edges(loaded.coordinates, k=6)
        region = coarse_region(meta["stage_region"])
        section: dict[str, object] = {
            "dataset": "KUPPE2022",
            "sample": sample,
            "patient": meta["patient"],
            "stage_region": meta["stage_region"],
            "coarse_region": region,
            "n_array_spots": n_array_spots,
            "n_spots": loaded.matrix.shape[1],
            "matrix_scale": loaded.matrix_scale,
            "n_graph_edges": len(edges),
            "author_label_columns": ";".join(author_label_columns(loaded.obs_columns)),
        }
        for name, values in combined.items():
            finite = values[np.isfinite(values)]
            section[f"mean_{name}"] = float(np.mean(finite)) if len(finite) else np.nan
            section[f"median_{name}"] = float(np.median(finite)) if len(finite) else np.nan
            section[f"p90_{name}"] = float(np.percentile(finite, 90)) if len(finite) else np.nan
            if name in {
                "mechanical_border_score",
                "immune_fibrotic_activation_score",
                "fibroblast_scar_repair_score",
                "boundary_transition_index",
                "source_bz_baseline",
                "generic_fibrosis_baseline",
            }:
                section[f"graph_morans_i_{name}"] = graph_morans_i(values, edges)
        for module, count in detected.items():
            section[f"detected_genes_{module}"] = count
        section_rows.append(section)

        label_columns = author_label_columns(loaded.obs_columns)
        for spot_idx, barcode in enumerate(loaded.matrix.barcodes):
            spot: dict[str, object] = {
                "dataset": "KUPPE2022",
                "sample": sample,
                "patient": meta["patient"],
                "stage_region": meta["stage_region"],
                "coarse_region": region,
                "barcode": barcode,
                "spatial_x": float(loaded.coordinates[spot_idx, 0]),
                "spatial_y": float(loaded.coordinates[spot_idx, 1]),
                "author_labels": ";".join(
                    f"{label}={loaded.obs_columns[label][spot_idx]}" for label in label_columns
                ),
            }
            for name, values in combined.items():
                spot[name] = float(values[spot_idx])
            spot_rows.append(spot)
        print(
            f"[scored] {sample} region={region} spots={loaded.matrix.shape[1]} "
            f"scale={loaded.matrix_scale} mechanical={detected['CM_BZ2_MECHANICAL_EDGE']} "
            f"scar={detected['CTHRC1_REPARATIVE_CF']}"
        )

    write_tsv(OUT_ROOT / "kuppe_human_spatial_scores_by_spot.tsv", spot_rows)
    write_tsv(OUT_ROOT / "kuppe_human_section_summary.tsv", section_rows)
    print(f"patients={len({row['patient'] for row in section_rows})} sections={len(section_rows)} spots={len(spot_rows)}")


if __name__ == "__main__":
    main()
