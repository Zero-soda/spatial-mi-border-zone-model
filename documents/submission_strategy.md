# Submission Strategy

## Recommended Primary Route

The manuscript should now be prepared primarily for Cardiovascular Research as a public-data-only spatial cardiovascular modelling Original Article. This is the most rational first submission because the study is mechanistic, spatial and translational, but it does not include local clinical samples or wet-lab perturbation validation. The strengthened version includes a Cardiovascular Research-style structured abstract, Translational Perspective, human STEMI transfer figure, leave-one-sample-out analysis, module-dropout testing, boundary-direction testing, boundary-threshold sensitivity, random-signature negative controls and human score-separation robustness analyses.

## Journal Ladder

| Priority | Journal | 2024 JIF | Fit | Probability |
|---|---:|---:|---|---|
| Recommended primary | Cardiovascular Research | 13.3 | Best balance for public spatial data, mechanistic interpretation, falsification controls and translational cardiovascular modelling | 25-35% |
| High-risk stretch | Nature Cardiovascular Research | 10.8 | Excellent fit for spatial cardiovascular biology, but risky without independent wet-lab or multi-human validation | 8-15% |
| Very high risk | Circulation Research | 16.2 | Strong cardiovascular prestige, but likely to request experimental mechanism | 5-10% |
| General biology backup | Communications Biology | 5.1 | Good for rigorous spatial biology with bounded claims and public-data robustness | 45-60% |
| Method-oriented backup | Briefings in Bioinformatics | 7.7 | Good if code, benchmark and reusable modelling workflow are emphasized | 35-50% |
| Digital/model backup | European Heart Journal - Digital Health | 4.4 | Reasonable if reframed as spatial digital biomarker or prediction workflow | 35-50% |

## Current Strengths

- Uses public spatial transcriptomics with full tissue-coordinate recovery.
- Includes six mouse MI Visium replicates across day 3 and day 7.
- Scores 14,147 mouse tissue spots and 1,551 human STEMI spots.
- Separates mechanical-border, immune-fibrotic activation and fibroblast-scar repair outputs.
- Quantifies the domain 3/4 spatial interface rather than relying only on cluster annotation.
- Adds public-data-only robustness: leave-one-sample-out 28/28 checks preserved, module-dropout 15/15 checks preserved, boundary-direction 4/4 checks preserved.
- Adds falsification checks: boundary-threshold sensitivity preserved 10/10 stage-threshold directions, and the primary mechanical-border and fibroblast-scar repair axes exceeded the 99th percentile of equal-size random-signature controls.
- Shows human STEMI score separation: mechanical-border output is largely distinct from immune-fibrotic and fibroblast-scar repair outputs, while immune-fibrotic and scar-repair outputs are related but not identical.
- Provides a direct path toward perturbation modelling in mini-border-zone systems.

## Current Weaknesses Before High-Impact Submission

- Human STEMI transfer currently uses a single spatial sample.
- No independent histology, immunostaining or clinical remodelling outcome validation is included.
- No prospective patient-derived perturbation experiment is included.
- The current model is interpretable but not yet benchmarked against graph-based alternatives.
- Current statistical inference should be conservative because spot-level observations are spatially correlated.

## Recommended Pre-Submission Additions

1. Completed: package the manuscript documents, figure files, source data, supplementary tables and reproducibility statements under `submission/`.
2. Completed: revise the manuscript to Cardiovascular Research structure with a structured abstract and Translational Perspective.
3. Completed: upgrade Figure 5 into a human STEMI transfer-validation figure with module coverage, spatial maps, score distributions, correlations and hotspot overlap.
4. Completed: add boundary-threshold sensitivity and random-signature negative-control analyses as falsification support.
5. Completed as a local release skeleton: prepare a GitHub/Zenodo-ready archive under `submission/repository_release`; insert the final public repository URL and Zenodo DOI before formal upload.
6. If another public human MI spatial dataset can be accessed without controlled data restrictions, add it as an external validation set.
7. If no extra public human dataset is practically accessible, keep the claims bounded and emphasize reproducibility, robustness and transfer feasibility.
8. Confirm remaining co-author and corresponding-author ORCID identifiers if required by the submission system, and confirm whether the institution wants a formal public-data secondary-analysis ethics exemption letter.

## Practical Recommendation

With no local hospital data and no wet-lab validation, the most rational order is Cardiovascular Research, then Communications Biology, then Briefings in Bioinformatics or iScience. Nature Cardiovascular Research can remain a high-risk stretch option only if the authors want to test the top-end novelty ceiling before moving to the more realistic primary route.
