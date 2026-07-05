# Spatial MI Border-Zone Prediction Model

This repository-ready release supports the manuscript `Spatial boundary modelling reveals a mechanical-to-fibroblast scar transition in the myocardial infarction border zone`.

## Scope

The study is a secondary computational analysis of public spatial transcriptomic data. It re-analyses mouse MI Visium samples and a human STEMI Visium sample from GEO accession GSE214611 to build an interpretable three-output spatial boundary model.

## Main Outputs

- `mechanical_border_score`: mechanical cardiomyocyte edge and border-zone transition.
- `immune_fibrotic_activation_score`: myeloid, interferon, fibro-inflammatory and matrix activation.
- `fibroblast_scar_repair_score`: POSTN/CTHRC1-associated reparative scar maturation.

## Directory Map

- `scripts/`: Python scripts used for preprocessing, scoring, robustness analyses and figure generation.
- `source_data/`: source-data files underlying main and extended-data figures.
- `supplementary_tables/`: processed supplementary tables.
- `figures/`: high-resolution figure files for manuscript review and upload.
- `statements/`: data, code, ethics and reproducibility statements.
- `documents/`: manuscript markdown, figure legends and submission strategy drafts.

## Data Provenance

Raw third-party data should be downloaded from GEO accession GSE214611. This release does not redistribute controlled, private or newly generated hospital data. Processed tables are derived from public de-identified datasets and are provided for figure inspection and reuse.

## Reproducibility

Install the packages listed in `requirements.txt`, download GSE214611 raw files from GEO, and run the scripts in `scripts/` following the order described in `documents/code_reproducibility_readme.md`.

## Ethics

This release contains only public, de-identified secondary data and generated processed outputs. It includes no new human recruitment, patient specimens, animal experiments, identifiable clinical data or wet-lab procedures.

## Citation

Please cite the final published article and the archived Zenodo DOI: https://doi.org/10.5281/zenodo.21203189.

## Repository

GitHub: https://github.com/Zero-soda/spatial-mi-border-zone-model
Zenodo: https://doi.org/10.5281/zenodo.21203189
