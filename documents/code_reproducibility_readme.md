# Code Reproducibility README

## Analysis Scope

This package supports the manuscript `Spatial boundary modelling reveals a mechanical-to-fibroblast scar transition in the myocardial infarction border zone`.

## Main Scripts

- `scripts/batch_map_gse214611_mi_spatial_risk.py`: maps mouse Visium samples to spatial domains and prototype risk summaries.
- `scripts/score_gse214611_visium_signatures.py`: computes mouse spot-level signature and composite scores.
- `scripts/analyze_gse214611_domain34_spatial_features.py`: computes domain 3/4 boundary distances, boundary fractions and contact gradients.
- `scripts/score_human_stemi_spatial_signatures.py`: transfers human signatures to the GSM6613090 human STEMI Visium sample.
- `scripts/plot_cvr_human_stemi_transfer.py`: generates the Cardiovascular Research-style human STEMI transfer figure with module coverage, spatial maps, score distributions, correlation matrix and hotspot-overlap analysis.
- `scripts/robustness_public_data_only.py`: runs leave-one-sample-out, module-dropout, boundary-direction, boundary-threshold, random-signature negative-control and human score-separation robustness checks.
- `scripts/plot_nature_public_data_robustness.py`: generates the polished public-data-only robustness figure.

## Input Data

Raw public data are from GEO accession GSE214611. Processed figure-level source-data tables are provided under `source_data/`, and supplementary analysis tables are provided under `supplementary_tables/`.

## Recommended Pre-Submission Archiving

This release is deposited at https://github.com/Zero-soda/spatial-mi-border-zone-model and archived on Zenodo at https://doi.org/10.5281/zenodo.21203189. Use these identifiers in the manuscript Code Availability section.
