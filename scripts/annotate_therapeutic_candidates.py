#!/usr/bin/env python3
"""Annotate spatial MI candidates with versioned public pharmacology evidence."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

from project_paths import project_root
from typing import Any, Iterable, Mapping


ROOT = project_root(__file__)
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))

from pharmacology_evidence_clients import (  # noqa: E402
    action_from_dgidb_types,
    cache_manifest,
    cardiovascular_disease_summary,
    chembl_mechanisms,
    chembl_target_search,
    dgidb_gene_records,
    dgidb_genes,
    open_targets_target,
    open_targets_target_record,
    select_human_chembl_target,
    unwrap_response,
)
from therapeutic_prioritization_utils import (  # noqa: E402
    independent_source_count,
    normalize_action,
    write_tsv,
)


RESULTS = ROOT / "results" / "tables"
CACHE = Path(__file__).resolve().parents[1] / "cache" / "pharmacology"
HUMAN_SUPPORT = RESULTS / "gse214611_therapeutic_human_support.tsv"
EVIDENCE_OUT = RESULTS / "gse214611_therapeutic_pharmacology_evidence.tsv"
TARGET_OUT = RESULTS / "gse214611_therapeutic_target_annotations.tsv"
MANIFEST_OUT = RESULTS / "gse214611_therapeutic_pharmacology_cache_manifest.tsv"

CONTROL_SYMBOLS = {
    "IL1B",
    "CCR2",
    "CCL2",
    "TGFBR1",
    "TGFB1",
    "LOX",
    "POSTN",
    "CTHRC1",
    "NPPA",
    "NPPB",
    "XIRP2",
    "FLNC",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def select_candidates(rows: list[dict[str, str]], per_state: int) -> list[dict[str, str]]:
    selected: dict[str, dict[str, str]] = {}
    by_state: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_state[row["primary_state"]].append(row)
    for state_rows in by_state.values():
        state_rows.sort(
            key=lambda row: float(row["mouse_spatial_score"])
            + float(row["human_support_score"]),
            reverse=True,
        )
        for row in state_rows[:per_state]:
            selected[row["human_gene_name"]] = row
    for row in rows:
        if row["human_gene_name"] in CONTROL_SYMBOLS:
            selected[row["human_gene_name"]] = row
    return sorted(selected.values(), key=lambda row: row["human_gene_name"])


def _phase(value: Any) -> float:
    try:
        return float(value or 0)
    except (TypeError, ValueError):
        return 0.0


def _join(values: Iterable[Any]) -> str:
    return "; ".join(sorted({str(value).strip() for value in values if str(value).strip()}))


def _open_targets_rows(
    candidate: Mapping[str, str], target: Mapping[str, Any]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    clinical = target.get("drugAndClinicalCandidates") or {}
    for item in clinical.get("rows") or []:
        drug = item.get("drug") or {}
        warnings = drug.get("drugWarnings") or []
        diseases = item.get("diseases") or []
        rows.append(
            {
                "human_gene_id": candidate["human_gene_id"],
                "human_gene_name": candidate["human_gene_name"],
                "primary_state": candidate["primary_state"],
                "source": "Open Targets",
                "source_record_id": item.get("id", ""),
                "drug_id": drug.get("id", ""),
                "drug_name": drug.get("name", ""),
                "raw_action": "",
                "normalized_action": "unknown",
                "mechanism_text": "Open Targets drug-clinical candidate association",
                "target_type": drug.get("drugType", ""),
                "max_clinical_phase": max(
                    _phase(item.get("maxClinicalStage")),
                    _phase(drug.get("maximumClinicalStage")),
                ),
                "approval_or_warning": _join(
                    f"{warning.get('warningType', '')}: {warning.get('description', '')}".strip(": ")
                    for warning in warnings
                ),
                "disease_context": _join(
                    (entry.get("disease") or {}).get("name")
                    or entry.get("diseaseFromSource", "")
                    for entry in diseases
                ),
                "upstream_evidence_sources": "Open Targets Platform",
            }
        )
    return rows


def _chembl_rows(
    candidate: Mapping[str, str], target: Mapping[str, Any], payload: Mapping[str, Any]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for item in payload.get("mechanisms") or []:
        raw_action = str(item.get("action_type", ""))
        rows.append(
            {
                "human_gene_id": candidate["human_gene_id"],
                "human_gene_name": candidate["human_gene_name"],
                "primary_state": candidate["primary_state"],
                "source": "ChEMBL",
                "source_record_id": item.get("record_id", ""),
                "drug_id": item.get("molecule_chembl_id", ""),
                "drug_name": item.get("molecule_chembl_id", ""),
                "raw_action": raw_action,
                "normalized_action": normalize_action(raw_action),
                "mechanism_text": item.get("mechanism_of_action", ""),
                "target_type": target.get("target_type", ""),
                "max_clinical_phase": _phase(item.get("max_phase")),
                "approval_or_warning": "",
                "disease_context": "",
                "upstream_evidence_sources": "ChEMBL mechanism",
            }
        )
    return rows


def _dgidb_rows(
    candidate_by_symbol: Mapping[str, Mapping[str, str]], nodes: Iterable[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for node in nodes:
        symbol = str(node.get("name", "")).upper()
        candidate = candidate_by_symbol.get(symbol)
        if not candidate:
            continue
        for item in node.get("interactions") or []:
            drug = item.get("drug") or {}
            action, raw_action = action_from_dgidb_types(item.get("interactionTypes") or [])
            rows.append(
                {
                    "human_gene_id": candidate["human_gene_id"],
                    "human_gene_name": candidate["human_gene_name"],
                    "primary_state": candidate["primary_state"],
                    "source": "DGIdb",
                    "source_record_id": item.get("interactionScore", ""),
                    "drug_id": drug.get("conceptId", ""),
                    "drug_name": drug.get("name", ""),
                    "raw_action": raw_action,
                    "normalized_action": action,
                    "mechanism_text": "DGIdb interaction type",
                    "target_type": "",
                    "max_clinical_phase": 0,
                    "approval_or_warning": "",
                    "disease_context": "",
                    "upstream_evidence_sources": _join(
                        source.get("sourceDbName", "") for source in item.get("sources") or []
                    ),
                }
            )
    return rows


def deduplicate_evidence(rows: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str, str]] = set()
    for row in rows:
        key = (
            str(row.get("source", "")),
            str(row.get("human_gene_name", "")),
            str(row.get("drug_id", "") or row.get("drug_name", "")),
            str(row.get("normalized_action", "")),
        )
        if key in seen:
            continue
        seen.add(key)
        result.append(dict(row))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=30, help="Candidates per spatial state")
    parser.add_argument("--refresh", action="store_true", help="Refresh cached API responses")
    args = parser.parse_args()

    candidates = select_candidates(read_tsv(HUMAN_SUPPORT), max(1, args.limit))
    candidate_by_symbol = {row["human_gene_name"].upper(): row for row in candidates}
    evidence: list[dict[str, Any]] = []
    annotations: list[dict[str, Any]] = []
    chembl_by_target_drug: dict[tuple[str, str], dict[str, Any]] = {}
    resource_data: dict[str, dict[str, Any]] = {}

    for index, candidate in enumerate(candidates, start=1):
        symbol = candidate["human_gene_name"]
        print(f"[{index}/{len(candidates)}] {symbol}", flush=True)
        ot_wrapper = open_targets_target(candidate["human_gene_id"], CACHE, args.refresh)
        ot_target = open_targets_target_record(ot_wrapper)
        ot_rows = _open_targets_rows(candidate, ot_target)
        evidence.extend(ot_rows)

        search_wrapper = chembl_target_search(symbol, CACHE, args.refresh)
        search_payload = unwrap_response(search_wrapper)
        chembl_target = select_human_chembl_target(search_payload, symbol)
        mechanism_payload: Mapping[str, Any] = {}
        chembl_rows: list[dict[str, Any]] = []
        if chembl_target:
            mechanism_wrapper = chembl_mechanisms(
                str(chembl_target.get("target_chembl_id", "")), CACHE, args.refresh
            )
            mechanism_payload = unwrap_response(mechanism_wrapper)
            chembl_rows = _chembl_rows(candidate, chembl_target, mechanism_payload)
            evidence.extend(chembl_rows)
            for row in chembl_rows:
                chembl_by_target_drug[(symbol, str(row["drug_id"]))] = row

        cv_score, cv_diseases = cardiovascular_disease_summary(ot_target)
        tractability = ot_target.get("tractability") or []
        safety = ot_target.get("safetyLiabilities") or []
        resource_data[symbol] = {
            "candidate": candidate,
            "ot_target": ot_target,
            "ot_rows": ot_rows,
            "chembl_target": chembl_target or {},
            "chembl_rows": chembl_rows,
            "cv_score": cv_score,
            "cv_diseases": cv_diseases,
            "tractability": tractability,
            "safety": safety,
        }

    symbols = sorted(candidate_by_symbol)
    dgidb_rows: list[dict[str, Any]] = []
    for start in range(0, len(symbols), 40):
        wrapper = dgidb_genes(symbols[start : start + 40], CACHE, args.refresh)
        dgidb_rows.extend(_dgidb_rows(candidate_by_symbol, dgidb_gene_records(wrapper)))
    evidence.extend(dgidb_rows)

    # Use target-matched ChEMBL action calls to resolve corresponding Open Targets molecules.
    for row in evidence:
        if row["source"] != "Open Targets":
            continue
        match = chembl_by_target_drug.get((row["human_gene_name"], str(row["drug_id"])))
        if match:
            row["raw_action"] = match["raw_action"]
            row["normalized_action"] = match["normalized_action"]
            row["mechanism_text"] = match["mechanism_text"]

    evidence = deduplicate_evidence(evidence)
    evidence_by_symbol: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in evidence:
        evidence_by_symbol[row["human_gene_name"]].append(row)

    dgidb_by_symbol: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in dgidb_rows:
        dgidb_by_symbol[row["human_gene_name"]].append(row)

    for symbol in sorted(resource_data):
        data = resource_data[symbol]
        candidate = data["candidate"]
        target_rows = evidence_by_symbol.get(symbol, [])
        phases = [_phase(row.get("max_clinical_phase")) for row in target_rows]
        annotations.append(
            {
                **candidate,
                "open_targets_found": bool(data["ot_target"]),
                "open_targets_tractability_labels": _join(
                    item.get("label", "")
                    for item in data["tractability"]
                    if item.get("value") is True
                ),
                "open_targets_tractability_modalities": _join(
                    item.get("modality", "")
                    for item in data["tractability"]
                    if item.get("value") is True
                ),
                "open_targets_cv_disease_score": data["cv_score"],
                "open_targets_cv_diseases": data["cv_diseases"],
                "open_targets_is_essential": data["ot_target"].get("isEssential", ""),
                "open_targets_safety_liability_count": len(data["safety"]),
                "open_targets_safety_events": _join(
                    item.get("event", "") for item in data["safety"]
                ),
                "open_targets_drug_count": len(data["ot_rows"]),
                "chembl_target_id": data["chembl_target"].get("target_chembl_id", ""),
                "chembl_target_type": data["chembl_target"].get("target_type", ""),
                "chembl_mechanism_count": len(data["chembl_rows"]),
                "dgidb_interaction_count": len(dgidb_by_symbol.get(symbol, [])),
                "pharmacology_source_count": independent_source_count(target_rows),
                "known_drug_or_interaction_count": len(target_rows),
                "highest_clinical_phase": max(phases, default=0),
                "pharmacology_status": "available" if target_rows else "unavailable",
            }
        )

    evidence_fields = [
        "human_gene_id",
        "human_gene_name",
        "primary_state",
        "source",
        "source_record_id",
        "drug_id",
        "drug_name",
        "raw_action",
        "normalized_action",
        "mechanism_text",
        "target_type",
        "max_clinical_phase",
        "approval_or_warning",
        "disease_context",
        "upstream_evidence_sources",
    ]
    write_tsv(EVIDENCE_OUT, evidence, evidence_fields)
    write_tsv(TARGET_OUT, annotations, list(annotations[0]))
    manifest = cache_manifest(CACHE)
    write_tsv(MANIFEST_OUT, manifest, list(manifest[0]))
    print(f"Wrote {len(annotations)} targets and {len(evidence)} evidence rows")


if __name__ == "__main__":
    main()
