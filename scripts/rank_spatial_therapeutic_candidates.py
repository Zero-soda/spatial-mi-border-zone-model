#!/usr/bin/env python3
"""Rank direction-aware spatial therapeutic hypotheses with repair guardrails."""

from __future__ import annotations

import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

from project_paths import project_root
from typing import Any, Iterable, Mapping


ROOT = project_root(__file__)
LOCAL_DEPS = ROOT / ".deps" / "python"
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))

import numpy as np  # noqa: E402

from compute_mouse_therapeutic_spatial_evidence import (  # noqa: E402
    summarize_gene,
)
from therapeutic_prioritization_utils import (  # noqa: E402
    assign_tier,
    write_tsv,
)


TABLES = ROOT / "results" / "tables"
ANNOTATIONS = TABLES / "gse214611_therapeutic_target_annotations.tsv"
PHARMACOLOGY = TABLES / "gse214611_therapeutic_pharmacology_evidence.tsv"
MOUSE_SUMMARY = TABLES / "gse214611_therapeutic_mouse_gene_evidence_summary.tsv"
MOUSE_BY_SAMPLE = TABLES / "gse214611_therapeutic_mouse_gene_evidence_by_sample.tsv"
RANKING_OUT = TABLES / "gse214611_therapeutic_candidate_ranking.tsv"
PAIRS_OUT = TABLES / "gse214611_therapeutic_drug_target_pairs.tsv"
SENSITIVITY_OUT = TABLES / "gse214611_therapeutic_rank_sensitivity.tsv"
NULL_OUT = TABLES / "gse214611_therapeutic_spatial_score_null.tsv"


STRUCTURAL_CAUTION = {
    "ACTA1",
    "ACTC1",
    "ACTN2",
    "ANKRD1",
    "CASQ2",
    "CKB",
    "CRIP2",
    "CRYAB",
    "CSRP3",
    "DES",
    "FLNC",
    "HSPB1",
    "HSPB7",
    "LDB3",
    "MB",
    "MYBPC3",
    "MYH7",
    "MYL6",
    "MYOZ2",
    "NPPA",
    "NPPB",
    "SLC6A6",
    "SMPX",
    "SORBS2",
    "TNNC1",
    "TNNT2",
    "TPM1",
    "TTN",
    "XIRP2",
}
REPAIR_CAUTION = {
    "AEBP1",
    "BGN",
    "CILP",
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "COL5A1",
    "COL6A1",
    "COL6A2",
    "CTHRC1",
    "DCN",
    "FBLN1",
    "FBLN2",
    "FN1",
    "FSTL1",
    "LOX",
    "POSTN",
    "SPARC",
    "TGFBR1",
    "THBS2",
    "THBS4",
}
TEMPORAL_REPAIR_CAUTION = {"TGFB1", "TGFBR1", "TGFBR2"}
VASCULAR_REPAIR_CAUTION = {"CXCL12", "ENG"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(value: Any, default: float = 0.0) -> float:
    try:
        result = float(value)
        return result if math.isfinite(result) else default
    except (TypeError, ValueError):
        return default


def as_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def canonical_drug_key(row: Mapping[str, Any]) -> str:
    drug_id = str(row.get("drug_id", "")).strip().upper()
    match = re.search(r"CHEMBL\d+", drug_id)
    if match:
        return match.group(0)
    name = str(row.get("drug_name", "")).strip().upper()
    return drug_id or name


def aggregate_drug_target_pairs(
    evidence: Iterable[Mapping[str, Any]], candidate_state: Mapping[str, str]
) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str], list[Mapping[str, Any]]] = defaultdict(list)
    for row in evidence:
        key = canonical_drug_key(row)
        if key:
            groups[(str(row.get("human_gene_name", "")), key)].append(row)

    result: list[dict[str, Any]] = []
    for (symbol, drug_key), rows in sorted(groups.items()):
        sources = sorted({str(row.get("source", "")) for row in rows})
        actions = {
            str(row.get("normalized_action", "unknown"))
            for row in rows
            if str(row.get("normalized_action", "unknown")) != "unknown"
        }
        action = next(iter(actions)) if len(actions) == 1 else "unknown"
        state = candidate_state.get(symbol, str(rows[0].get("primary_state", "unresolved")))
        desired = "inhibit" if state in {"mechanical_adverse", "immune_fibrotic_adverse"} else "context_dependent"
        caution = symbol in STRUCTURAL_CAUTION or symbol in REPAIR_CAUTION or state == "scar_repair_associated"
        if action == "unknown":
            direction = "unknown"
        elif caution and action == "inhibit":
            direction = "protective_conflict"
        elif desired == "context_dependent":
            direction = "unresolved"
        elif action == desired:
            direction = "matched"
        else:
            direction = "mismatched"
        names = [
            str(row.get("drug_name", "")).strip()
            for row in rows
            if str(row.get("drug_name", "")).strip()
            and str(row.get("drug_name", "")).strip().upper() != drug_key
        ]
        result.append(
            {
                "human_gene_name": symbol,
                "human_gene_id": rows[0].get("human_gene_id", ""),
                "primary_state": state,
                "drug_key": drug_key,
                "drug_name": sorted(names, key=len)[0] if names else drug_key,
                "sources": "; ".join(sources),
                "pair_source_count": len(sources),
                "raw_actions": "; ".join(
                    sorted({str(row.get("raw_action", "")) for row in rows if row.get("raw_action")})
                ),
                "normalized_action": action,
                "desired_action": desired,
                "direction_match": direction,
                "mechanism_text": "; ".join(
                    sorted(
                        {str(row.get("mechanism_text", "")) for row in rows if row.get("mechanism_text")}
                    )
                ),
                "max_clinical_phase": max(as_float(row.get("max_clinical_phase")) for row in rows),
                "approval_or_warning": "; ".join(
                    sorted(
                        {
                            str(row.get("approval_or_warning", ""))
                            for row in rows
                            if row.get("approval_or_warning")
                        }
                    )
                ),
                "disease_context": "; ".join(
                    sorted(
                        {str(row.get("disease_context", "")) for row in rows if row.get("disease_context")}
                    )
                ),
                "upstream_evidence_sources": "; ".join(
                    sorted(
                        {
                            str(row.get("upstream_evidence_sources", ""))
                            for row in rows
                            if row.get("upstream_evidence_sources")
                        }
                    )
                ),
            }
        )
    return result


def score_pharmacology(annotation: Mapping[str, Any]) -> dict[str, float]:
    labels = str(annotation.get("open_targets_tractability_labels", "")).lower()
    if "approved drug" in labels:
        tractability = 8.0
    elif "advanced clinical" in labels:
        tractability = 6.0
    elif "high-quality ligand" in labels:
        tractability = 4.0
    elif "small molecule binder" in labels or "structure with ligand" in labels:
        tractability = 3.0
    elif "druggable family" in labels:
        tractability = 2.0
    else:
        tractability = 0.0

    cv_score = as_float(annotation.get("open_targets_cv_disease_score"))
    disease = min(6.0, cv_score * 12.0)
    phase = as_float(annotation.get("highest_clinical_phase"))
    direct_count = int(as_float(annotation.get("known_drug_or_interaction_count")))
    if phase >= 4:
        known_drug = 6.0
    elif phase >= 3:
        known_drug = 5.0
    elif phase >= 2:
        known_drug = 4.0
    elif phase >= 1:
        known_drug = 3.0
    elif direct_count:
        known_drug = 1.5
    else:
        known_drug = 0.0

    chembl_count = int(as_float(annotation.get("chembl_mechanism_count")))
    if not chembl_count:
        chembl = 0.0
    elif phase >= 4:
        chembl = 6.0
    elif phase >= 3:
        chembl = 5.0
    elif phase >= 2:
        chembl = 4.0
    elif phase >= 1:
        chembl = 3.0
    else:
        chembl = 2.0

    source_count = int(as_float(annotation.get("pharmacology_source_count")))
    source_confirmation = 4.0 if source_count >= 3 else 2.0 if source_count >= 2 else 0.0
    total = tractability + disease + known_drug + chembl + source_confirmation
    return {
        "tractability_score": tractability,
        "cv_disease_evidence_score": disease,
        "known_drug_score": known_drug,
        "chembl_mechanism_score": chembl,
        "source_confirmation_score": source_confirmation,
        "pharmacology_score": min(30.0, total),
    }


def repair_safety_penalty(
    symbol: str,
    state: str,
    mouse: Mapping[str, Any],
    annotation: Mapping[str, Any],
    pairs: Iterable[Mapping[str, Any]],
) -> tuple[float, str]:
    penalty = 0.0
    reasons: list[str] = []
    if state == "scar_repair_associated":
        penalty += 10.0
        reasons.append("day-7 scar-repair state")
        scar_strength = max(0.0, as_float(mouse.get("scar_state_strength")))
        if scar_strength >= 0.75:
            penalty += 2.0
            reasons.append("strong scar-score association")
    if symbol in REPAIR_CAUTION:
        penalty += 4.0
        reasons.append("matrix/scar-repair caution")
    if symbol in TEMPORAL_REPAIR_CAUTION:
        penalty += 12.0
        reasons.append("TGF-beta timing and wound-healing caution")
    if symbol in VASCULAR_REPAIR_CAUTION:
        penalty += 8.0
        reasons.append("vascular repair/angiogenesis caution")
    if symbol in STRUCTURAL_CAUTION or as_bool(mouse.get("structural_or_repair_caution")):
        penalty += 12.0
        reasons.append("cardiomyocyte structural/protective caution")
    safety_count = int(as_float(annotation.get("open_targets_safety_liability_count")))
    if safety_count:
        penalty += min(4.0, float(safety_count))
        reasons.append("Open Targets safety liability")
    if as_bool(annotation.get("open_targets_is_essential")):
        penalty += 2.0
        reasons.append("essential-gene annotation")
    if any(str(pair.get("approval_or_warning", "")).strip() for pair in pairs):
        penalty += 3.0
        reasons.append("drug warning/withdrawal annotation")
    return min(20.0, penalty), "; ".join(dict.fromkeys(reasons)) or "none identified"


def choose_best_pair(pairs: list[Mapping[str, Any]]) -> Mapping[str, Any] | None:
    if not pairs:
        return None
    direction_order = {"matched": 4, "unknown": 3, "unresolved": 2, "mismatched": 1, "protective_conflict": 0}
    return max(
        pairs,
        key=lambda row: (
            direction_order.get(str(row.get("direction_match", "unknown")), 0),
            int(as_float(row.get("pair_source_count"))),
            as_float(row.get("max_clinical_phase")),
        ),
    )


def build_sensitivity(
    ranking: list[dict[str, Any]], mouse_by_sample: list[dict[str, str]]
) -> list[dict[str, Any]]:
    scenarios: dict[str, dict[str, tuple[float, float, float, float]]] = {}
    primary: dict[str, tuple[float, float, float, float]] = {}
    for row in ranking:
        primary[row["human_gene_name"]] = (
            as_float(row["mouse_spatial_score"]),
            as_float(row["human_support_score"]),
            as_float(row["pharmacology_score"]),
            as_float(row["repair_safety_penalty"]),
        )
    scenarios["primary"] = primary

    component_factors = {
        "omit_mouse": (0.0, 1.0, 1.0, 1.0),
        "omit_human": (1.0, 0.0, 1.0, 1.0),
        "omit_pharmacology": (1.0, 1.0, 0.0, 1.0),
        "omit_penalty": (1.0, 1.0, 1.0, 0.0),
        "mouse_weight_75pct": (0.75, 1.0, 1.0, 1.0),
        "mouse_weight_125pct": (1.25, 1.0, 1.0, 1.0),
        "human_weight_75pct": (1.0, 0.75, 1.0, 1.0),
        "human_weight_125pct": (1.0, 1.25, 1.0, 1.0),
        "pharmacology_weight_75pct": (1.0, 1.0, 0.75, 1.0),
        "pharmacology_weight_125pct": (1.0, 1.0, 1.25, 1.0),
        "penalty_weight_75pct": (1.0, 1.0, 1.0, 0.75),
        "penalty_weight_125pct": (1.0, 1.0, 1.0, 1.25),
    }
    for name, factors in component_factors.items():
        scenarios[name] = {
            symbol: tuple(value * factor for value, factor in zip(values, factors, strict=True))
            for symbol, values in primary.items()
        }

    rows_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in mouse_by_sample:
        rows_by_gene[row["gene_name"]].append(row)
    sample_names = sorted({row["sample"] for row in mouse_by_sample})
    mouse_gene_by_human = {row["human_gene_name"]: row["mouse_gene_name"] for row in ranking}
    mouse_id_by_human = {row["human_gene_name"]: row["mouse_gene_id"] for row in ranking}
    for sample in sample_names:
        scenario: dict[str, tuple[float, float, float, float]] = {}
        for symbol, values in primary.items():
            mouse_gene = mouse_gene_by_human[symbol]
            subset = [row for row in rows_by_gene.get(mouse_gene, []) if row["sample"] != sample]
            if subset:
                summary = summarize_gene(mouse_gene, mouse_id_by_human[symbol], subset)
                mouse_score = as_float(summary.get("mouse_spatial_score"), values[0])
            else:
                mouse_score = values[0]
            scenario[symbol] = (mouse_score, values[1], values[2], values[3])
        scenarios[f"leave_out_{sample}"] = scenario

    scenario_rows: list[dict[str, Any]] = []
    rank_history: dict[str, list[int]] = defaultdict(list)
    for scenario_name, values_by_symbol in scenarios.items():
        scores = {
            symbol: values[0] + values[1] + values[2] - values[3]
            for symbol, values in values_by_symbol.items()
        }
        order = sorted(scores, key=lambda symbol: (-scores[symbol], symbol))
        ranks = {symbol: idx for idx, symbol in enumerate(order, start=1)}
        for symbol in order:
            rank_history[symbol].append(ranks[symbol])
            values = values_by_symbol[symbol]
            scenario_rows.append(
                {
                    "scenario": scenario_name,
                    "human_gene_name": symbol,
                    "mouse_component": values[0],
                    "human_component": values[1],
                    "pharmacology_component": values[2],
                    "penalty_component": values[3],
                    "scenario_score": scores[symbol],
                    "scenario_rank": ranks[symbol],
                    "top10": ranks[symbol] <= 10,
                }
            )
    for row in scenario_rows:
        history = rank_history[row["human_gene_name"]]
        row["median_rank_all_scenarios"] = float(np.median(history))
        row["rank_range_all_scenarios"] = f"{min(history)}-{max(history)}"
        row["top10_frequency_all_scenarios"] = sum(rank <= 10 for rank in history) / len(history)
    return scenario_rows


def build_spatial_null(
    ranking: list[dict[str, Any]], mouse_summary: list[dict[str, str]], seed: int = 20260711
) -> list[dict[str, Any]]:
    selected = {row["mouse_gene_name"] for row in ranking}
    pool = [row for row in mouse_summary if row["gene_name"] not in selected]
    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []
    for candidate in ranking:
        target_detection = as_float(candidate.get("detection_fraction_median"))
        target_effect = as_float(candidate.get("effect_magnitude_percentile"))
        distances = np.asarray(
            [
                abs(as_float(row.get("detection_fraction_median")) - target_detection)
                + abs(as_float(row.get("effect_magnitude_percentile")) - target_effect)
                for row in pool
            ]
        )
        nearest = np.argsort(distances)[: min(100, len(pool))]
        if len(nearest) > 30:
            nearest = rng.choice(nearest, size=30, replace=False)
        for null_index, index in enumerate(nearest, start=1):
            control = pool[int(index)]
            rows.append(
                {
                    "candidate_gene": candidate["human_gene_name"],
                    "candidate_mouse_spatial_score": candidate["mouse_spatial_score"],
                    "null_index": null_index,
                    "matched_mouse_gene": control["gene_name"],
                    "matched_mouse_spatial_score": control["mouse_spatial_score"],
                    "detection_distance": abs(
                        as_float(control.get("detection_fraction_median")) - target_detection
                    ),
                    "effect_percentile_distance": abs(
                        as_float(control.get("effect_magnitude_percentile")) - target_effect
                    ),
                    "seed": seed,
                }
            )
    return rows


def main() -> None:
    annotations = read_tsv(ANNOTATIONS)
    pharmacology = read_tsv(PHARMACOLOGY)
    mouse_summary = read_tsv(MOUSE_SUMMARY)
    mouse_by_sample = read_tsv(MOUSE_BY_SAMPLE)
    mouse_by_id = {row["gene_id"]: row for row in mouse_summary}
    state_by_symbol = {row["human_gene_name"]: row["primary_state"] for row in annotations}
    pairs = aggregate_drug_target_pairs(pharmacology, state_by_symbol)
    pairs_by_symbol: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for pair in pairs:
        pairs_by_symbol[pair["human_gene_name"]].append(pair)

    ranking: list[dict[str, Any]] = []
    for annotation in annotations:
        symbol = annotation["human_gene_name"]
        mouse = mouse_by_id.get(annotation["mouse_gene_id"], {})
        target_pairs = pairs_by_symbol.get(symbol, [])
        best_pair = choose_best_pair(target_pairs)
        scores = score_pharmacology(annotation)
        penalty, penalty_reasons = repair_safety_penalty(
            symbol, annotation["primary_state"], mouse, annotation, target_pairs
        )
        direction = str(best_pair.get("direction_match", "unknown")) if best_pair else "unknown"
        pair_sources = int(as_float(best_pair.get("pair_source_count"))) if best_pair else 0
        row: dict[str, Any] = {
            **annotation,
            **{
                key: mouse.get(key, "")
                for key in (
                    "structural_or_repair_caution",
                    "n_relevant_samples",
                    "relevant_effect_median",
                    "relevant_direction_consistency",
                    "boundary_slope_consistency",
                    "graph_edge_consistency",
                    "leave_one_out_direction_consistency",
                    "intended_score_correlation_median",
                    "detection_fraction_median",
                    "replicate_rule_pass",
                    "mechanical_state_strength",
                    "immune_state_strength",
                    "scar_state_strength",
                    "effect_magnitude_percentile",
                )
            },
            **scores,
            "repair_safety_penalty": penalty,
            "repair_safety_reasons": penalty_reasons,
            "best_drug_key": best_pair.get("drug_key", "") if best_pair else "",
            "best_drug_name": best_pair.get("drug_name", "") if best_pair else "",
            "best_drug_action": best_pair.get("normalized_action", "unknown") if best_pair else "unknown",
            "best_pair_sources": best_pair.get("sources", "") if best_pair else "",
            "best_pair_source_count": pair_sources,
            "direction_match": direction,
        }
        row["final_priority_score"] = (
            as_float(row["mouse_spatial_score"])
            + as_float(row["human_support_score"])
            + scores["pharmacology_score"]
            - penalty
        )
        tier_input = dict(row)
        tier_input["pharmacology_source_count"] = pair_sources
        tier_input["human_detected"] = as_bool(annotation.get("human_detected"))
        tier_input["replicate_rule_pass"] = as_bool(mouse.get("replicate_rule_pass"))
        row["priority_tier"] = assign_tier(tier_input)
        row["interpretive_label"] = {
            "tier_a": "spatially anchored actionable hypothesis",
            "tier_b": "promising but context-dependent hypothesis",
            "tier_c": "exploratory hypothesis",
            "caution_protective": "repair/structural preservation caution",
        }[row["priority_tier"]]
        ranking.append(row)

    ranking.sort(key=lambda row: (-as_float(row["final_priority_score"]), row["human_gene_name"]))
    for rank, row in enumerate(ranking, start=1):
        row["overall_rank"] = rank

    rank_by_symbol = {row["human_gene_name"]: row for row in ranking}
    for pair in pairs:
        target = rank_by_symbol[pair["human_gene_name"]]
        pair["mouse_spatial_score"] = target["mouse_spatial_score"]
        pair["human_support_score"] = target["human_support_score"]
        pair["pharmacology_score"] = target["pharmacology_score"]
        pair["repair_safety_penalty"] = target["repair_safety_penalty"]
        pair["target_final_priority_score"] = target["final_priority_score"]
        pair["target_priority_tier"] = target["priority_tier"]
        pair["pair_qualifies_for_tier_a"] = (
            target["priority_tier"] == "tier_a"
            and pair["direction_match"] == "matched"
            and int(pair["pair_source_count"]) >= 2
        )
        pair["claim_boundary"] = "therapeutic hypothesis; not efficacy or treatment recommendation"

    sensitivity = build_sensitivity(ranking, mouse_by_sample)
    spatial_null = build_spatial_null(ranking, mouse_summary)
    write_tsv(RANKING_OUT, ranking, list(ranking[0]))
    write_tsv(PAIRS_OUT, pairs, list(pairs[0]))
    write_tsv(SENSITIVITY_OUT, sensitivity, list(sensitivity[0]))
    write_tsv(NULL_OUT, spatial_null, list(spatial_null[0]))
    print(
        f"Wrote {len(ranking)} ranked targets, {len(pairs)} drug-target pairs, "
        f"{len(sensitivity)} sensitivity rows and {len(spatial_null)} matched null rows"
    )


if __name__ == "__main__":
    main()
