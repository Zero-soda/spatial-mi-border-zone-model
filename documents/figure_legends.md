# Figure Legends

## Figure 1. Spatial boundary modelling framework for myocardial infarction border-zone states.

Schematic workflow showing the study design. Public myocardial infarction spatial transcriptomic data were analysed using curated cardiac injury, mechanical stress, immune-fibrotic and scar-repair signatures. Spot-level module scores were combined into three interpretable outputs: mechanical-border score, immune-fibrotic activation score and fibroblast-scar repair score. Domain 3/4 spatial contact analysis was then used to quantify boundary structure, and the same signature logic was transferred to human STEMI Visium tissue.

Source data: not applicable; conceptual workflow. Dataset and analysis provenance are provided in `submission/source_data/Source_Data_Figure_1_readme.md`, `submission/supplementary_tables/Supplementary_Table_1_dataset_sample_inventory.tsv`, `submission/supplementary_tables/Supplementary_Table_2_signature_definitions.tsv` and `submission/code_reproducibility_readme.md`.

Alt text: Workflow diagram showing public mouse MI Visium data entering signature scoring, spatial state modelling and human STEMI transfer validation.

## Figure 2. Mouse MI Visium maps identify domain 3 and domain 4 as distinct spatial programs.

Contact sheet of GSE214611 day 3 and day 7 mouse MI Visium samples showing author-derived spatial domains and prototype spatial risk structure. Six MI biological replicates were analysed: D3_1, D3_2, D3_3, D7_1, D7_2 and D7_3. Across these samples, 14,147 tissue spots were mapped and scored. Domain 3 and domain 4 show stage-dependent spatial organization, motivating explicit boundary modelling rather than one-dimensional risk ranking.

Source image: `results/figures/gse214611_d3_d7_spatial_risk_contact_sheet.png`.
Source data: `submission/source_data/Source_Data_Figure_2_mouse_spatial_risk_spot_level.tsv`; `submission/source_data/Source_Data_Figure_2_mouse_stage_domain_summary.tsv`.

Alt text: Contact sheet of six mouse MI spatial transcriptomic samples showing day 3 and day 7 domain maps and spatial risk patterns.

## Figure 3. Signature score maps distinguish mechanical-border and fibroblast-scar repair axes.

Spatial maps of mechanical-border, fibroblast-scar repair and composite fibrotic-risk scores across day 3 and day 7 mouse MI Visium samples. Scores are shown at the spot level after expression-matrix gene-set scoring and within-analysis z-score normalization. Domain 3 preferentially shows mechanical-border activation, while domain 4 acquires a strong fibroblast-scar repair program by day 7.

Source files:
`results/figures/gse214611_d3_d7_mechanical_border_contact_sheet.png`;
`results/figures/gse214611_d3_d7_fibroblast_scar_contact_sheet.png`;
`results/figures/gse214611_d3_d7_fibrotic_risk_contact_sheet.png`.
Source data: `submission/source_data/Source_Data_Figure_3_mouse_signature_scores_by_spot.tsv`; `submission/source_data/Source_Data_Figure_3_mouse_signature_scores_by_stage_domain.tsv`.

Alt text: Spatial signature score maps comparing mechanical-border, fibroblast-scar repair and fibrotic-risk signals across mouse MI samples.

## Figure 4. Domain 3/4 boundary analysis reveals temporal remodelling of the infarct interface.

Spatial boundary maps and contact-gradient analysis of the domain 3/4 interface. Boundary spots were defined using a threshold equal to 1.35 times the median nearest-neighbour distance. Quantitative summaries are calculated at the biological-replicate level across three day 3 MI and three day 7 MI Visium samples. Day 3 MI showed broader domain 3/4 intermixing, whereas day 7 MI showed increased domain 3-to-domain 4 distance and a positive fibroblast-scar repair gradient toward domain 4.

Source image: `results/figures/gse214611_d3_d7_domain34_boundary_contact_sheet.png`.
Source data: `submission/source_data/Source_Data_Figure_4_domain34_boundary_by_spot.tsv`; `submission/source_data/Source_Data_Figure_4_domain34_boundary_by_sample.tsv`; `submission/source_data/Source_Data_Figure_4_domain34_boundary_by_stage.tsv`.

Alt text: Domain 3 and domain 4 boundary maps with quantitative summaries showing temporal changes in interface mixing and scar-repair gradients.

## Figure 5. Human STEMI transfer of the spatial state model.

Transfer of human gene signatures to the GSM6613090 human STEMI Visium sample. a, Fixed-signature transfer design showing mouse MI spatial model derivation, direct human gene-symbol scoring and no model retraining. b, Detected human signature modules supporting the three model outputs. c, Spatial maps of mechanical-border, immune-fibrotic activation and fibroblast-scar repair scores across the human STEMI tissue section. d, Spot-level score distributions with 90th and 95th percentile markers. e, Pairwise Spearman correlation matrix among the three model outputs. f, Top-decile hotspot Jaccard overlap; bar labels show overlapping hotspot spots divided by the union of hotspot spots. The analysis included 1,551 tissue spots. Spearman correlations and Jaccard overlaps were calculated across spatial spots as descriptive transfer metrics; spot-level P values are not reported because Visium spots are spatially correlated. Because this analysis uses a single human spatial sample and no clinical outcome labels, it should be interpreted as migration feasibility evidence rather than definitive clinical validation.

Source image: `results/figures/gse214611_human_stemi_signature_transfer.png`; submission-ready PDF, SVG and TIFF exports are generated from the same plotting script.
Source data: `submission/source_data/Source_Data_Figure_5_human_stemi_signature_scores_by_spot.tsv`; `submission/source_data/Source_Data_Figure_5_human_stemi_summary.tsv`; `submission/source_data/Source_Data_Figure_5_human_stemi_score_separation.tsv`.

Alt text: Multi-panel human STEMI validation figure showing fixed-signature transfer, module detectability, spatial score maps, score distributions, correlation matrix and hotspot overlap.

## Extended Data Figure 1. Sample-level signature maps.

Single-sample triptych maps for each mouse MI Visium replicate showing mechanical-border, fibroblast-scar repair and composite fibrotic-risk scores.

Source files: `results/figures/gse214611_*_signature_score_maps.png`.

Alt text: Single-sample mouse MI maps showing spatial patterns of mechanical, scar-repair and fibrotic-risk scores.

## Extended Data Figure 2. Domain 3/4 boundary maps by sample.

Single-sample boundary maps showing domain 3 and domain 4 boundary spots and non-boundary spots.

Source files: `results/figures/gse214611_*_domain34_boundary_map.png`.

Alt text: Single-sample boundary maps highlighting domain 3 and domain 4 interface spots and non-boundary spots.

## Extended Data Figure 3. Public-data-only robustness and falsification summary.

Summary of robustness and falsification analyses that do not require local clinical samples or wet-lab validation. Leave-one-sample-out analysis preserved all 28 prespecified directional checks. Module-dropout analysis preserved all 15 directional checks after removing one component module at a time. Boundary-direction testing preserved all four expected stage-specific contact-gradient directions across the six mouse MI Visium biological replicates. Boundary-threshold sensitivity preserved all 10 stage-threshold directional checks across thresholds from 1.10 to 1.75 times the median nearest-neighbour distance. Random-signature negative controls showed that the two mechanical-border checks and the day 7 fibroblast-scar repair check exceeded the 99th percentile of 1,000 equal-size random signatures; the immune-fibrotic check remained directionally positive but was treated as contextual rather than domain-defining. Human STEMI score-separation analysis was performed across 1,551 human tissue spots and showed that mechanical-border hotspots were largely distinct from immune-fibrotic and fibroblast-scar repair hotspots, whereas immune-fibrotic and scar-repair outputs were correlated but only partially overlapping.

Source image: `results/figures/nature_ed_public_data_robustness.png`.
Source data: `submission/source_data/Source_Data_Extended_Data_Figure_3_loso_robustness.tsv`; `submission/source_data/Source_Data_Extended_Data_Figure_3_module_dropout_robustness.tsv`; `submission/source_data/Source_Data_Extended_Data_Figure_3_boundary_direction_robustness.tsv`; `submission/source_data/Source_Data_Extended_Data_Figure_3_boundary_threshold_sensitivity.tsv`; `submission/source_data/Source_Data_Extended_Data_Figure_3_random_signature_negative_control.tsv`; `submission/source_data/Source_Data_Extended_Data_Figure_3_human_score_separation.tsv`.

Alt text: Robustness summary figure and supplementary source data showing leave-one-sample-out, module-dropout, boundary-direction, boundary-threshold, random-signature and human score-separation checks.
