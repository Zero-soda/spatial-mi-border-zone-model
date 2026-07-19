# Code Reproducibility README

## Analysis Scope

This package supports the Problem Solving Protocol `SSTBA: a reproducible cross-study spatial-state and boundary-analysis workflow for myocardial infarction transcriptomics`.

## Main Scripts

- `scripts/spatial_model_utils.py`: provides shared sample-level, permutation and spatial-statistical helper functions.
- `scripts/strengthen_spatial_statistics.py`: generates sample-level effect-size, spatial autocorrelation, domain-label permutation, spatial-block permutation and signed boundary-gradient tables.
- `scripts/audit_signatures_and_sensitivity.py`: audits signature overlap, high-impact gene removal sensitivity and score collinearity.
- `scripts/reviewer_rigor_score_state_celltype_audit.py`: runs score-only state clustering and cell-type marker-adjusted domain contrast audits.
- `scripts/final_reviewer_upgrade_analyses.py`: generates data-provenance/QC, module-detection, alternative-scoring, expression-matched random-control, graph-boundary, H&E image-intensity and master spot-level source tables.
- `scripts/additional_independence_analyses.py`: runs graph-smoothed domain-independent state clustering and human expression-matched spatial-null checks.
- `scripts/external_public_evidence_inventory.py`: records independent public evidence sources and claim-strength boundaries.
- `scripts/plot_bib_main_figures.py`: generates BIB Figures 1-5 and the graphical abstract from released processed data and technical component figures.
- `scripts/plot_cvr_spatial_statistics.py`: generates the spatial-statistics strengthening component retained in the supplementary materials.
- `scripts/plot_reviewer_rigor_checks.py`: generates the reviewer-requested rigor-check figure for domain dependence and cell-type composition.
- `scripts/plot_final_reviewer_upgrade_checks.py`: generates the final reviewer-requested provenance, scoring, graph-boundary and image-audit figure.
- `scripts/plot_additional_independence_checks.py`: generates Supplementary Figure S7.
- `scripts/batch_map_gse214611_mi_spatial_risk.py`: maps mouse Visium samples to spatial domains and prototype risk summaries.
- `scripts/score_gse214611_visium_signatures.py`: computes mouse spot-level signature and composite scores.
- `scripts/analyze_gse214611_domain34_spatial_features.py`: computes domain 3/4 boundary distances, boundary fractions and contact gradients.
- `scripts/score_human_stemi_spatial_signatures.py`: transfers human signatures to the GSM6613090 human STEMI Visium sample.
- `scripts/plot_cvr_human_stemi_transfer.py`: generates the single-section human STEMI feasibility figure with module coverage, spatial maps, score distributions, correlation matrix and hotspot-overlap analysis.
- `scripts/robustness_public_data_only.py`: runs leave-one-sample-out, module-dropout, boundary-direction, boundary-threshold, random-signature negative-control and human score-separation robustness checks.
- `scripts/plot_nature_public_data_robustness.py`: generates the polished public-data-only robustness figure.

- `scripts/compute_mouse_therapeutic_spatial_evidence.py`: computes all-gene mouse section-level spatial evidence.
- `scripts/map_mouse_human_orthologues.py`: creates high-confidence one-to-one mouse-human orthologue mappings.
- `scripts/compute_human_therapeutic_support.py`: computes within-section human target support and matched spatial nulls.
- `scripts/annotate_therapeutic_candidates.py`: queries and caches Open Targets, ChEMBL and DGIdb evidence.
- `scripts/rank_spatial_therapeutic_candidates.py`: applies direction, repair/safety and rank-sensitivity rules.
- `scripts/plot_therapeutic_prioritization.py`: generates the optional therapeutic-prioritization component summarized in Supplementary Figure S8.
- `scripts/plot_therapeutic_rank_sensitivity.py`: generates Supplementary Figure S8.

- `scripts/download_spatial_mi_external_validation.py`: downloads external public data with range-resume, checksum and retrieval logging.
- `scripts/validate_external_mouse_spatial_states.py`: applies frozen scores to GSE176092 and GSE265828.
- `scripts/validate_kuppe_human_spatial_states.py`: applies frozen human scores to the prespecified nine-patient Kuppe panel.
- `scripts/benchmark_cross_study_spatial_model.py`: benchmarks transfer directions against simple prespecified baselines.
- `scripts/plot_cross_study_external_validation.py`: generates the external-transfer component incorporated into BIB Figure 5 and Supplementary Figure S9.

## Input Data

Raw public data are from GSE214611, GSE176092, GSE265828 and the Kuppe human MI atlas. Exact external URLs, version identifiers, expected sizes and checksums are recorded in the retrieval manifests; derived tables are under `results/tables`.

## Reproduction Order

Run commands from the repository root. The following order separates data acquisition, frozen scoring, robustness checks, external transfer and plotting:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python scripts/batch_map_gse214611_mi_spatial_risk.py
python scripts/score_gse214611_visium_signatures.py
python scripts/analyze_gse214611_domain34_spatial_features.py
python scripts/score_human_stemi_spatial_signatures.py

python scripts/robustness_public_data_only.py
python scripts/audit_signatures_and_sensitivity.py
python scripts/strengthen_spatial_statistics.py
python scripts/reviewer_rigor_score_state_celltype_audit.py
python scripts/final_reviewer_upgrade_analyses.py
python scripts/additional_independence_analyses.py

python scripts/download_spatial_mi_external_validation.py --dry-run
python scripts/download_spatial_mi_external_validation.py
python scripts/validate_external_mouse_spatial_states.py
python scripts/validate_kuppe_human_spatial_states.py
python scripts/benchmark_cross_study_spatial_model.py
python scripts/plot_cross_study_external_validation.py
python scripts/plot_bib_main_figures.py
```

Therapeutic-prioritization scripts are run only after the spatial-state tables above have been generated:

```bash
python scripts/compute_mouse_therapeutic_spatial_evidence.py
python scripts/map_mouse_human_orthologues.py
python scripts/compute_human_therapeutic_support.py
python scripts/annotate_therapeutic_candidates.py
python scripts/rank_spatial_therapeutic_candidates.py
```

The `config/` directory contains the frozen signature table and exact external manifests. `scripts/project_paths.py` makes the same scripts usable from this standalone release and from the larger manuscript-development workspace.

## Versioned Archiving

The public repository is `https://github.com/Zero-soda/spatial-mi-border-zone-model`. Use the immutable version tag and the version-specific Zenodo DOI cited by the manuscript when reproducing the submitted analysis.
