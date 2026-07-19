# SSTBA: Spatial-State Transfer and Boundary Analysis

This BIB-ready local release supports the Problem Solving Protocol `SSTBA: a
reproducible cross-study spatial-state and boundary-analysis workflow for
myocardial infarction transcriptomics`.

Version `0.5.0` adds a tested, configuration-driven entry point for frozen
module scoring, coordinate or nearest-neighbour boundary analysis and
machine-readable run manifests.

## Scope

The study is a secondary computational analysis of public spatial transcriptomic data. It derives an interpretable three-output spatial-state framework in GSE214611 and transfers frozen scores to GSE176092, GSE265828 and a nine-patient Kuppe human MI panel.

## Main Outputs

- `mechanical_border_score`: mechanical cardiomyocyte edge and border-zone transition.
- `immune_fibrotic_activation_score`: myeloid, interferon, fibro-inflammatory and matrix activation.
- `fibroblast_scar_repair_score`: POSTN/CTHRC1-associated reparative scar maturation.

An optional therapeutic extension illustrates how the state outputs can feed a repair-aware hypothesis-prioritization analysis. Its tiers are experimental hypotheses, not clinical guidance or evidence of treatment efficacy.

## Directory Map

- `scripts/`: Python scripts used for preprocessing, scoring, robustness analyses and figure generation.
- `config/`: frozen signature definitions and versioned external-download manifests.
- `source_data/`: source-data files underlying main and supplementary figures.
- `supplementary_tables/`: processed supplementary tables.
- `figures/`: high-resolution figure files for manuscript review and upload.
- `statements/`: data, code, ethics and reproducibility statements.
- `documents/`: ordered code and data reproducibility instructions.

## Data Provenance

Raw third-party data should be downloaded from GSE214611, GSE176092, GSE265828 and the versioned Kuppe CELLxGENE/Zenodo resources listed in the retrieval manifests. This release does not redistribute controlled, private or newly generated hospital data. Processed tables are derived from public de-identified datasets and are provided for figure inspection and reuse.

## Reproducibility

Install the packages listed in `requirements.txt`, obtain the public inputs listed in the manifests, and run the scripts in `scripts/` following the order described in `documents/code_reproducibility_readme.md`.

Scripts resolve paths in both the manuscript-development workspace and a standalone checkout of this repository. Run commands from the repository root; downloaded raw files are expected under `data/raw`, and generated tables and figures are written under `results`.

## Quick Start

```bash
python3 -m pip install -e ".[test]"
python3 -m pytest tests -v
python3 -m sstba.cli validate --config config/demo_boundary.json
python3 -m sstba.cli run --config config/demo_boundary.json --output results/demo
python3 -m compileall -q sstba scripts
python3 scripts/generate_file_manifest.py
```

The packaged demonstration uses all six discovery sections and does not
require a raw-data download. It is intended to verify graph construction,
domain 3-domain 4 interface detection, score-gradient orientation and manifest
generation. Its 18 automated tests include a regression against the released
six-section graph summaries. Raw count-matrix scoring remains available
through `sstba score` after the public input matrices are downloaded.

## SSTBA Input and Output Contract

- `sstba score` accepts a spots-by-genes TSV matrix and a two-column
  `module`/`gene` signature table. Gene symbols are matched
  case-insensitively, zero-variance genes contribute zero and missing-module
  failures identify the affected module.
- `sstba boundary` and `sstba run` accept a JSON configuration and a spot table
  containing sample, spot, label, two coordinate and one or more score
  columns.
- Graphs may use a fixed radius or symmetrized k-nearest neighbours. The
  `exact` implementation explicitly uses NumPy quicksort ordering to reproduce
  the discovery analysis when Visium neighbours are exactly equidistant. The
  run manifest records the NumPy version. The `kdtree` option is faster for
  large irregular datasets but is not used for the frozen six-section demo.
- Boundary outputs are `boundary_summary.tsv`, `boundary_spots.tsv` and
  `run_manifest.json`. Scoring outputs are `scores.tsv`,
  `signature_coverage.json` and `run_manifest.json`.
- Each manifest stores the command, configuration or input hashes, parameters,
  signature coverage where applicable, output hashes, dependency versions and
  warnings. Output paths are relative to the manifest directory so that
  identical runs remain comparable across machines.
- `config/protocol_schema.json` is the machine-readable contract. The CLI
  performs the corresponding checks without adding a JSON-schema runtime
  dependency, including required fields, graph settings, distinct labels,
  unique score columns and input-table column availability.

These scores and boundaries are descriptive computational summaries. They are
not causal mechanisms, clinical predictions, prospective validation or
treatment recommendations.

## Ethics

This release contains only public, de-identified secondary data and generated processed outputs. It includes no new human recruitment, patient specimens, animal experiments, identifiable clinical data or wet-lab procedures.

## Citation

Please cite the final published article and the versioned archive. The stable
Zenodo concept DOI is `10.5281/zenodo.21203188`. The immutable DOI for release
`v0.5.0` is `10.5281/zenodo.21444703`.
