#!/usr/bin/env python3
"""Create conservative external public-evidence inventory for spatial MI revision."""

from __future__ import annotations

import csv
from pathlib import Path

from project_paths import project_root

ROOT = project_root(__file__)
TABLE_DIR = ROOT / "results" / "tables"
DOCS_DIR = ROOT / "documents"

ROWS = [
    {
        "source_short_name": "Kuppe_human_MI_spatial_multiomic",
        "disease_context": "human myocardial infarction",
        "evidence_type": "human spatial multiomic atlas",
        "access_status_for_revision": "check_processed_access_before_projection",
        "use_if_accessible": "external human spatial projection of fixed signatures",
        "fallback_if_not_accessible": "literature-supported orthogonal evidence for human MI tissue-state architecture",
        "claim_strength": "external_projection_if_processed_spatial_data_available; otherwise orthogonal_support_only",
        "url": "https://www.nature.com/articles/s41586-022-05060-x",
    },
    {
        "source_short_name": "acute_mouse_MI_spatialomics_imaging",
        "disease_context": "mouse acute myocardial infarction",
        "evidence_type": "spatial transcriptomics plus immunofluorescence or imaging-associated resource",
        "access_status_for_revision": "check_repository_and_license_before_projection",
        "use_if_accessible": "orthogonal spatial/imaging support for immune and scar-repair axes",
        "fallback_if_not_accessible": "literature-supported orthogonal evidence, not validation",
        "claim_strength": "supportive_public_evidence",
        "url": "https://www.nature.com/articles/s44161-025-00717-y",
    },
    {
        "source_short_name": "public_MI_scRNA_or_snRNA",
        "disease_context": "mouse or human myocardial infarction / failing ischemic heart",
        "evidence_type": "single-cell or single-nucleus RNA-seq",
        "access_status_for_revision": "search_and_select_if_processed_matrices_available",
        "use_if_accessible": "cell-state interpretability of three signatures",
        "fallback_if_not_accessible": "omit from quantitative analysis and discuss as future direction",
        "claim_strength": "cell_state_interpretability_not_spatial_validation",
        "url": "GEO/CELLxGENE/HuBMAP search required during implementation",
    },
    {
        "source_short_name": "public_MI_bulk_timecourse",
        "disease_context": "post-MI temporal repair",
        "evidence_type": "bulk RNA-seq time course",
        "access_status_for_revision": "optional_low_priority",
        "use_if_accessible": "temporal plausibility of immune and scar-repair score trajectories",
        "fallback_if_not_accessible": "omit because bulk data cannot validate border-zone spatial localization",
        "claim_strength": "temporal_support_only",
        "url": "GEO search required during implementation",
    },
]


def write_tsv() -> None:
    path = TABLE_DIR / "spatial_mi_external_public_evidence_inventory.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=list(ROWS[0].keys()),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(ROWS)


def write_md() -> None:
    path = DOCS_DIR / "spatial_mi_external_public_evidence_inventory.md"
    lines = [
        "# External Public Evidence Inventory",
        "",
        "This inventory separates direct external projection from literature-supported orthogonal evidence. Sources listed here do not constitute validation unless processed data are accessible and the fixed signatures can be projected without parameter fitting.",
        "",
        "| Source | Evidence type | Potential use | Claim strength |",
        "|---|---|---|---|",
    ]
    for row in ROWS:
        lines.append(
            f"| {row['source_short_name']} | {row['evidence_type']} | {row['use_if_accessible']} | {row['claim_strength']} |"
        )
    lines.append("")
    lines.append("Implementation rule: if a source cannot be downloaded or processed, describe it only as orthogonal literature support and do not call it validation.")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    write_tsv()
    write_md()
    print("Wrote external public evidence inventory")


if __name__ == "__main__":
    main()
