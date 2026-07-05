# Spatial boundary modelling reveals a mechanical-to-fibroblast scar transition in the myocardial infarction border zone

## Article Type

Original Article

## Target Journals

Primary target: Cardiovascular Research (Original Article; basic and translational cardiovascular research).

Formatting basis: structured abstract, Translational Perspective, conservative public-data validation claims and up to seven main figures.

## Authors

Lizhe Zhong1†, Jingwen Shi1,2†, Hongtao Liao1, Hai Deng1, Yuanhong Liang1, Jun Huang3, Yingjie Huang1, Yumei Xue1, Zili Liao1, Fangzhou Liu1, Weidong Lin1, Huiqiang Wei1, Shulin Wu1, Xianhong Fang1*, Wei Wei1*

1 Department of Cardiology, Guangdong Cardiovascular Institute, Guangdong Provincial People's Hospital (Guangdong Academy of Medical Sciences), Southern Medical University, Guangzhou 510080, China

2 The Second School of Clinical Medicine, Southern Medical University, Guangzhou 510515, China

3 Department of Geriatrics, Institute of Geriatrics, Guangdong Provincial People's Hospital (Guangdong Academy of Medical Sciences), Southern Medical University, Guangzhou 510080, China

## Corresponding Author

Equal contribution: Lizhe Zhong and Jingwen Shi contributed equally as co-first authors.

Corresponding authors: Xianhong Fang, PhD, and Wei Wei, PhD, Department of Cardiology, Guangdong Cardiovascular Institute, Guangdong Provincial People's Hospital, Guangdong Academy of Medical Sciences, 106 Zhongshan 2nd Road, Yuexiu District, Guangzhou 510080, China. Emails: drfangxh@163.com; vvsquare@163.com.

Telephone: +86-13724893386.

## Abstract

### Aims

The myocardial infarction border zone contains spatially interleaved cardiomyocyte stress, inflammation and fibroblast repair programs. We aimed to convert public spatial transcriptomic maps into an interpretable boundary-state model that separates mechanical-border injury from immune-fibrotic activation and fibroblast-scar maturation.

### Methods and Results

We re-analysed public GSE214611 Visium data, including six mouse myocardial infarction samples at day 3 and day 7 after injury and one human STEMI sample. Curated gene modules were summarized into three prespecified outputs: mechanical-border score, immune-fibrotic activation score and fibroblast-scar repair score. Across 14,147 mouse tissue spots, annotated domain 3 showed the strongest mechanical-border signal, whereas domain 4 evolved into a fibroblast-scar repair hotspot by day 7. At the stage-domain level, domain 3 mechanical-border score increased from 1.78 at day 3 to 2.84 at day 7, while domain 4 fibroblast-scar repair score increased from 0.90 to 6.96. Domain 3/4 contact analysis showed early intermixing at day 3 and a day 7 fibroblast-scar gradient toward domain 4. Leave-one-sample-out, module-dropout, boundary-direction and threshold-sensitivity analyses preserved the primary directional conclusions. Equal-size random-signature controls placed the main mechanical and scar-repair margins above the random 99th percentile. Without retraining, transfer to 1,551 human STEMI tissue spots detected all core modules and showed separable mechanical, immune-fibrotic and scar-repair outputs.

### Conclusion(s)

The infarct border zone can be modelled as a dynamic spatial state system linking mechanical injury, inflammatory activation and reparative scar maturation. Human STEMI transfer supports model portability, while remaining feasibility evidence rather than definitive clinical validation.

## Keywords

myocardial infarction; border zone; spatial transcriptomics; cardiac fibrosis; scar repair; translational modelling

## Translational Perspective

Post-infarction repair is often interpreted as a linear inflammatory-to-fibrotic process, yet clinical remodelling emerges from spatially organized tissue states. This public-data model separates mechanical-border injury, immune-fibrotic activation and fibroblast-scar repair in mouse MI and transfers the same signatures to human STEMI tissue without retraining. The framework is not yet a clinical prognostic test, but it provides a transparent spatial readout that can guide future validation in human MI cohorts, histology-aligned tissue sections and perturbable cardiac repair models.

## Introduction

Adverse left ventricular remodelling after myocardial infarction is shaped by inflammatory and reparative processes that are intrinsically spatial. The infarct core, border zone and remote myocardium differ in tissue survival, mechanical load, inflammatory infiltration, extracellular matrix deposition and electrical function. Classical models of post-infarction repair describe a transition from inflammatory clearance to fibroblast activation, matrix deposition and scar maturation, but these phases do not occur uniformly across the ventricular wall [1-3]. In clinical practice, however, these compartments are often treated as broad anatomical categories, even though the border zone itself contains multiple local microenvironments with different biological trajectories.

Recent single-cell and spatial studies have begun to resolve this complexity. Spatial transcriptomics enables gene-expression measurements to be retained in tissue coordinates [4], and human myocardial infarction atlases have mapped immune, stromal, vascular and cardiomyocyte states across injured cardiac tissue [5]. Mouse spatial transcriptomic studies have refined the infarct border zone into cardiomyocyte and stromal neighbourhoods, including a mechanical edge marked by Nppa and Xirp2 and fibroblast-rich regions enriched for matrix and matricellular programs [6,7]. Recent work has further shown that the border zone can initiate spatially clustered inflammatory programs, including type I interferon responses in mechanically stressed cardiomyocyte regions [8]. Together, these studies support the concept that the infarct border zone is a structured tissue interface rather than a simple anatomical rim.

Single-cell studies of cardiac injury and heart failure have also clarified the cellular axes that should be represented in such a model. Cardiac stromal, immune and vascular populations undergo dynamic post-injury changes [9]. In human heart failure and acute infarction, CCR2-positive myeloid cells and IL1B signalling have been linked to FAP/POSTN fibroblast trajectories [10,11]. Separately, CTHRC1-positive reparative cardiac fibroblasts have been implicated in scar formation and early post-MI repair [12,13], while periostin has long been associated with cardiac hypertrophy, fibrosis and ventricular remodelling [14]. These findings suggest that mechanical stress, inflammatory fibroblast activation and reparative scar maturation should be separated rather than collapsed into a single fibrosis score.

Despite this progress, most analyses remain descriptive or atlas-oriented. Emerging spatial pathology-score methods indicate that spatial transcriptomic data can be transformed into tissue-level disease scores [15], but myocardial infarction border-zone biology still lacks an interpretable, cardiovascular-specific model that distinguishes mechanical border stress from immune-fibrotic activation and reparative scar formation. The unresolved problem is not simply whether injured tissue can be assigned a high score; it is whether different biological routes to remodelling can be separated in tissue coordinates and tested against spatial falsification controls. A practical translational model should answer a more specific question: given a tissue location and its neighbourhood context, is that region dominated by mechanical border-zone stress, immune-fibrotic activation or fibroblast-scar repair? Such a model would provide a bridge between spatial transcriptomic atlases, in vitro perturbation systems and clinical strategies that aim to modulate infarct healing.

Here we establish an interpretable spatial boundary model for myocardial infarction border-zone states. We re-analysed GSE214611 mouse MI Visium data at day 3 and day 7 after injury [6,16] using curated gene signatures for cardiomyocyte transition, mechanical edge activation, pathogenic fibroblast activation, reparative fibroblast states, myofibroblast activation, myeloid inflammation, TGF-beta signalling, extracellular matrix remodelling, endothelial inflammation, hypoxia, cardiomyocyte electrical function and acute stress. We then combined these modules into three model outputs: mechanical-border score, immune-fibrotic activation score and fibroblast-scar repair score. The framework was designed to be transparent rather than black-box: each output is biologically interpretable, each spatial boundary can be re-quantified and each main axis can be challenged by module dropout, threshold sensitivity and random-signature controls. Finally, we quantified domain 3/4 spatial contact features and transferred the same model logic to a human STEMI Visium sample.

## Results

### A spatial modelling framework for myocardial infarction border-zone states

To model infarct border-zone heterogeneity, we curated cardiac injury and remodelling signatures from prior myocardial infarction, heart failure and fibrosis studies [1-14]. The signatures were grouped into biological axes representing cardiomyocyte border-zone transition, cardiomyocyte mechanical edge response, FAP/POSTN fibroblast activation, CTHRC1 reparative fibroblast state, contractile myofibroblast activation, CCR2/IL1B myeloid inflammation, TGF-beta signalling, extracellular matrix remodelling, endothelial inflammatory fibrosis, hypoxia-ischemia, cardiomyocyte electrical/calcium-handling function and immediate early injury response.

These signatures were converted into three interpretable spatial outputs. The mechanical-border score combined cardiomyocyte transition and mechanical edge modules. The immune-fibrotic activation score combined myeloid, TGF-beta and pathogenic fibroblast activation modules. The fibroblast-scar repair score combined FAP/POSTN fibroblast, extracellular matrix, CTHRC1 reparative fibroblast and contractile myofibroblast modules. This formulation was chosen to avoid reducing infarct healing to a single risk axis. In particular, CTHRC1 and POSTN-associated scar repair may indicate organized reparative scar maturation rather than uniformly maladaptive fibrosis.

### Mouse MI spatial transcriptomics identifies distinct domain 3 and domain 4 border-zone programs

We analysed six mouse MI Visium samples from GSE214611: three day 3 MI samples and three day 7 MI samples [6,16]. Spatial barcodes were mapped to tissue coordinates and author-provided domain annotations. Across the six samples, 14,147 tissue spots were scored using expression matrix-derived gene-set scores.

The spatial maps showed that annotated domains were not interchangeable risk labels. Domain 3 and domain 4 occupied distinct spatial regions and displayed different temporal behaviour. In day 3 MI, domain 3 showed stronger mechanical-border activation than domain 4, consistent with an early infarct-edge cardiomyocyte stress program. By day 7, domain 4 showed a marked increase in fibroblast-scar repair signal, consistent with a more mature reparative scar or activated fibroblast-rich hotspot.

### Domain 3 marks a mechanical-border axis, whereas domain 4 evolves into a fibroblast-scar repair axis

Stage-domain signature scoring supported a dual-axis interpretation. In day 3 MI, domain 3 showed a mean mechanical-border score of 1.78, whereas domain 4 showed a lower mechanical-border score. Domain 4 at day 3 showed only modest fibroblast-scar repair enrichment, with a mean fibroblast-scar repair score of 0.90. In day 7 MI, domain 3 retained a high mechanical-border score of 2.84, but domain 4 became the dominant fibroblast-scar repair compartment, with a mean fibroblast-scar repair score of 6.96.

These results indicate that domain 3 is best interpreted as a mechanical-border axis, while domain 4 is a temporally evolving fibroblast-scar axis. The distinction is important because a one-dimensional fibrotic risk score would obscure the fact that mechanical stress and reparative scar maturation have different spatial distributions and may require different therapeutic interpretations.

### Domain 3/4 boundary analysis reveals a transition from early intermixing to scar maturation

To test whether domain 3 and domain 4 form a spatially meaningful interface, we quantified nearest-neighbour contact between the two domains. Boundary spots were defined using a contact threshold equal to 1.35 times the median nearest-neighbour distance in full-resolution pixel coordinates. In the analysed data, this threshold was approximately 371 full-resolution pixels.

At day 3, the mean fraction of domain 3 spots located at the domain 3/4 boundary was 0.576, and the mean fraction of domain 4 boundary spots was 0.417. At day 7, these fractions were 0.446 and 0.398, respectively. The median domain 3-to-domain 4 distance increased from 275.7 pixels at day 3 to 410.2 pixels at day 7, suggesting that the early interface became more spatially separated as the tissue state matured.

The boundary score gradients further supported this interpretation. At the domain 3/4 contact interface, domain 4 remained lower than domain 3 for mechanical-border score at both stages. In contrast, the contact-pair fibroblast-scar repair gradient shifted toward domain 4 by day 7, with a day 7 domain 4-minus-domain 3 contact delta of 0.75. Thus, day 3 MI is characterized by a broad and intermingled mechanical-fibrotic border, whereas day 7 MI shows emergence of a fibroblast-scar repair hotspot across the domain 3/4 interface.

### Falsification and public-data-only robustness analyses support the domain-state interpretation

Because this study is designed as a public-data-only computational analysis, we performed additional robustness tests that did not require local clinical samples or wet-lab validation. First, leave-one-sample-out analysis preserved all 28 prespecified directional checks. The day 3 domain 3 mechanical-border advantage over domain 4 remained positive after each sample omission, as did the day 7 domain 3 mechanical-border advantage and the day 7 domain 4 fibroblast-scar repair advantage. The day 3 domain 4 fibroblast-scar repair signal also remained non-dominant relative to domain 3.

Second, module-dropout testing showed that no single biological module drove the main interpretation. Removing either cardiomyocyte transition or mechanical-edge module preserved the domain 3 mechanical-border advantage at day 3 and day 7. Removing FAP/POSTN fibroblast, extracellular matrix, CTHRC1 reparative fibroblast or contractile myofibroblast components preserved the day 7 domain 4 fibroblast-scar repair advantage. Removing CCR2/IL1B myeloid, TGF-beta or FAP/POSTN components also preserved the day 7 domain 4 immune-fibrotic activation advantage. Overall, all 15 module-dropout checks preserved the expected direction.

Third, boundary-direction checks were consistent across all analysed biological replicates. Mechanical-border contact gradients were negative for domain 4 minus domain 3 at both day 3 and day 7, whereas fibroblast-scar repair contact gradients were negative at day 3 but positive at day 7. Together, these analyses support the robustness of the spatial-state interpretation within the available public data.

Fourth, boundary-threshold sensitivity analysis preserved the main contact-gradient directions across threshold multipliers from 1.10 to 1.75 times the median nearest-neighbour distance. Across all 10 stage-threshold combinations, mechanical-border contact gradients remained negative for domain 4 minus domain 3, and fibroblast-scar repair gradients remained negative at day 3 but positive at day 7. Thus, the boundary interpretation was not dependent on a single arbitrary contact threshold.

Fifth, random-signature negative controls tested whether the observed domain-state margins could be reproduced by arbitrary equal-size gene modules. For the two domain 3 mechanical-border checks and the day 7 domain 4 fibroblast-scar repair check, the observed margins exceeded the 99th percentile of 1,000 random signatures of the same component number, with no random control matching or exceeding the observed margin. The day 7 immune-fibrotic domain 4 advantage remained directionally positive but did not exceed the random 95th percentile, supporting its interpretation as a contextual inflammatory-fibrotic coupling signal rather than the primary domain-defining axis.

### Human STEMI migration supports conserved spatial deployment of model outputs

To assess translational portability, we applied the same signature logic to the human STEMI Visium sample GSM6613090 from GSE214611. Human signature genes were used directly, and no model retraining was performed. The analysis included 1,551 tissue spots. All core gene signatures had detectable genes in the human matrix, including 4 genes for cardiomyocyte transition, 6 for mechanical edge, 8 for FAP/POSTN fibroblast activation, 7 for CTHRC1 reparative fibroblasts, 7 for CCR2/IL1B myeloid activation and 8 for extracellular matrix remodelling.

To make the transfer assessment explicit, we quantified signature detectability, score dispersion, pairwise score correlations and top-decile hotspot overlap in the same human tissue section. The transferred maps revealed non-uniform spatial deployment of the three model outputs. The 90th percentile values were 2.04 for mechanical-border score, 2.75 for immune-fibrotic activation score and 3.82 for fibroblast-scar repair score. Pairwise score comparisons showed that the mechanical-border output was spatially distinct from both immune-fibrotic activation and fibroblast-scar repair in the human sample, with Spearman correlations of -0.098 and -0.050, respectively. In contrast, immune-fibrotic activation and fibroblast-scar repair were correlated but not identical, with a Spearman correlation of 0.738 and a top-decile hotspot Jaccard overlap of 0.268. These results support the decision to separate inflammatory fibrosis from reparative scar maturation. Because this human analysis uses a single spatial sample and lacks independent clinical outcome labels, it should be interpreted as migration feasibility evidence rather than definitive clinical validation.

## Discussion

This study proposes an interpretable spatial model of myocardial infarction border-zone states. By re-analysing mouse MI Visium data and transferring the same signature logic to human STEMI tissue, we show that the infarct border zone can be represented as a spatial state system rather than a single high-risk compartment. The main finding is a transition from a domain 3 mechanical-border axis to a domain 4 fibroblast-scar repair axis, with the domain 3/4 interface serving as a quantifiable spatial boundary.

The results refine how infarct border-zone risk should be conceptualized. A simple fibrotic risk score can be misleading because organized reparative scar formation and maladaptive fibrosis share extracellular matrix and matricellular genes. This distinction is consistent with the dual nature of post-infarction matrix remodelling: extracellular matrix deposition is required for scar stability, but excessive or persistent fibrosis contributes to adverse remodelling [2,3]. The observed day 7 domain 4 fibroblast-scar repair enrichment, together with the human STEMI transfer maps, supports a three-output representation: mechanical border stress, immune-fibrotic activation and fibroblast-scar repair. This structure is more biologically faithful and more useful for perturbation modelling than a single aggregate risk score.

The boundary analysis provides a second conceptual advance. Rather than treating spatial domains as isolated clusters, the model quantifies their contact interface. This is important because prior spatial studies have shown that infarct border-zone biology is shaped by local mechanical destabilization, mechanosensing programs and spatially organized inflammatory responses [6-8]. Day 3 MI showed broader domain 3/4 intermixing, whereas day 7 MI showed increased distance and a fibroblast-scar repair gradient favouring domain 4. These features can serve as graph-model inputs and as readouts for engineered mini-border-zone systems in which mechanical load, hypoxia/reoxygenation, IL1B, TGF-beta or patient-derived serum are imposed as perturbations.

The public-data-only design also defines the appropriate level of claim. The strengthened analysis does not substitute for prospective human outcome data or wet-lab perturbation experiments, but it does address a common weakness of public-data reanalysis: dependence on a single marker set, a single biological replicate or a single spatial threshold. The preservation of the main conclusions in leave-one-sample-out, module-dropout, threshold-sensitivity and boundary-direction analyses suggests that the model captures a reproducible spatial signal within GSE214611 rather than an idiosyncratic feature of one sample, one marker module or one contact definition. Random-signature controls further indicate that the primary mechanical-border and fibroblast-scar repair axes are not readily reproduced by arbitrary equal-size modules. The human STEMI transfer analysis indicates that the three-output framework can be applied to human tissue coordinates, while remaining clearly bounded as migration evidence rather than clinical validation.

This study has limitations. The discovery analysis is based on public spatial transcriptomics and curated signatures rather than new experimental tissue profiling. The human STEMI migration analysis uses one Visium sample and does not yet include independent histological annotation, patient-level outcome data or wet-lab perturbation validation. The model therefore supports a spatial-state hypothesis and a draft predictive framework, but it does not yet establish clinical prognostic performance. Future work should validate the three outputs in additional human MI samples, align them with histology and imaging-derived remodelling outcomes, and test whether patient-derived inflammatory or mechanical perturbations shift mini-border-zone models along the same spatial axes.

In summary, the infarct border zone can be modelled as a dynamic spatial boundary system linking mechanical cardiomyocyte injury, immune-fibrotic activation and fibroblast-scar repair. This framework provides a tractable route from spatial transcriptomic atlases to perturbable cardiovascular prediction models.

## Methods

### Data sources

Spatial transcriptomic data were obtained from GSE214611 [16], generated in the study by Calcagno et al. [6]. Mouse MI Visium samples included three day 3 MI biological replicates and three day 7 MI biological replicates. The human STEMI migration analysis used GSM6613090, labelled as Visium, STEMI. Processed 10x Genomics Visium expression matrices and spatial coordinates were used.

### Signature curation

Gene signatures were curated for cardiomyocyte transition, cardiomyocyte mechanical edge response, FAP/POSTN fibroblast activation, CTHRC1 reparative fibroblasts, contractile myofibroblasts, CCR2/IL1B myeloid activation, TGF-beta signalling, extracellular matrix remodelling, endothelial inflammatory fibrosis, hypoxia-ischemia, cardiomyocyte electrical/calcium-handling function and injury stress. Signature selection was guided by prior work on post-MI inflammation and fibrosis [1-3], human and mouse MI spatial atlases [5-8], immune-fibroblast communication [10,11] and reparative fibroblast biology [12,13]. Mouse analyses used mouse gene symbols, and human STEMI transfer used human gene symbols from the same signature table.

### Spot-level scoring

For each Visium sample, filtered feature-barcode matrices were read from 10x HDF5 files. Counts were library-size normalized to counts per 10,000 and log-transformed using log1p. For each spot and each signature, a module score was calculated as the mean normalized expression of detected signature genes. Signature scores were z-scored within the analysed dataset before composite scores were calculated.

### Composite model outputs

The mechanical-border score was calculated as the sum of z-scored cardiomyocyte border-zone transition and cardiomyocyte mechanical-edge modules. The immune-fibrotic activation score was calculated from z-scored CCR2/IL1B myeloid, TGF-beta signalling and FAP/POSTN fibroblast modules. The fibroblast-scar repair score was calculated from z-scored FAP/POSTN fibroblast, extracellular matrix remodelling, CTHRC1 reparative fibroblast and contractile myofibroblast modules.

### Spatial boundary analysis

Domain 3/4 boundary analysis was performed using full-resolution pixel coordinates. The contact threshold was defined as 1.35 times the median nearest-neighbour distance. A domain 3 spot was considered a boundary spot if it had at least one domain 4 spot within this threshold, and vice versa. Boundary gradients were calculated for contact pairs as the score in the domain 4 spot minus the score in the paired domain 3 spot.

### Public-data-only robustness analyses

Robustness analyses were designed to avoid dependence on local clinical samples or wet-lab validation. Leave-one-sample-out analysis recalculated stage-domain directional checks after omitting each mouse MI Visium sample in turn. Module-dropout analysis recalculated composite outputs after removing one component module at a time from the mechanical-border, immune-fibrotic activation or fibroblast-scar repair score. Boundary-direction robustness tested whether each biological replicate preserved the expected sign of domain 4-minus-domain 3 contact gradients. Boundary-threshold sensitivity recalculated contact-gradient directions using thresholds from 1.10 to 1.75 times the median nearest-neighbour distance. Random-signature negative controls generated 1,000 equal-size random modules for each tested axis and compared the observed domain-state margin with the empirical random-margin distribution. Human STEMI score separation was assessed by pairwise Pearson correlation, Spearman correlation and top-decile hotspot Jaccard overlap among the three model outputs.

### Human STEMI transfer analysis

The GSM6613090 human STEMI Visium matrix and spatial coordinates were analysed using the same scoring workflow, substituting human gene symbols. The model was not retrained on the human sample. Human transfer was evaluated by the number of detected signature genes, the number of scored tissue spots and spatial deployment of the three composite model outputs.

## Statistics and Reproducibility

No statistical method was used to predetermine sample size. The mouse discovery and temporal validation analyses included all six selected public GSE214611 MI Visium samples available for the prespecified day 3 and day 7 comparison. The human transfer analysis used the public GSM6613090 human STEMI Visium sample. Spot-level values were used for spatial visualization, local-neighbourhood analysis and sensitivity analyses; because Visium spots are spatially correlated, spot-level observations were not treated as independent biological replicates for primary inference.

The primary quantitative unit was the biological sample or the stage-domain summary derived from sample-level values. Directional conclusions were assessed using leave-one-sample-out analysis, single-module dropout analysis, replicate-level domain 3/4 boundary-direction checks, boundary-threshold sensitivity and random-signature negative controls. Public-data-only robustness preserved 28 of 28 leave-one-sample-out checks, 15 of 15 module-dropout checks, 4 of 4 boundary-direction checks and 10 of 10 stage-threshold sensitivity checks. The two mechanical-border checks and the day 7 fibroblast-scar repair check exceeded the 99th percentile of equal-size random-signature controls. Randomization and blinding were not applicable to this secondary computational analysis of public transcriptomic data.

## Data Availability

This study re-analysed publicly available spatial transcriptomic data from GEO accession GSE214611 [16]. The mouse MI Visium samples analysed here correspond to GSM6613084, GSM6613085, GSM6613086, GSM6613087, GSM6613088 and GSM6613089, and the human STEMI Visium transfer analysis used GSM6613090. No new primary human, animal or clinical data were generated in this study. Source data underlying the main and extended-data figures, processed supplementary tables and repository metadata are available at GitHub (https://github.com/Zero-soda/spatial-mi-border-zone-model) and archived on Zenodo (https://doi.org/10.5281/zenodo.21203189).

## Code Availability

Analysis scripts used to generate the processed tables and figures are available in the public GitHub repository (https://github.com/Zero-soda/spatial-mi-border-zone-model) and archived on Zenodo (https://doi.org/10.5281/zenodo.21203189). The repository includes the figure-generation workflow, source-data tables, supplementary analysis tables, figure files and reproducibility notes.

## Ethics Statement

This manuscript is a secondary computational analysis of publicly available, de-identified transcriptomic datasets. The original studies obtained the relevant approvals and consent as reported in their source publications and public repositories. No new human participants, patient specimens, animal experiments, identifiable clinical data or wet-lab procedures were generated by the authors for this analysis; therefore, no additional institutional ethics approval or informed consent was required for the present study.

## AI Assistance Disclosure

AI-assisted software was used to support manuscript drafting, code organization and figure-layout refinement under author supervision. The authors reviewed and verified the analysis code, source data, figure outputs, citations, interpretation and final manuscript text, and take full responsibility for the content of the work.

## Funding

This work was supported by the Natural Science Foundation of Guangdong Province (No. 2015A030313657) awarded to W.W. The funder had no role in study design, public-data analysis, interpretation, manuscript preparation or the decision to submit the work for publication.

## Acknowledgements

The authors thank the investigators who generated and deposited the GSE214611 spatial transcriptomic dataset in the Gene Expression Omnibus.

## Competing Interests

The authors declare no competing interests.

## References

1. Frangogiannis, N. G. The inflammatory response in myocardial injury, repair, and remodelling. Nat. Rev. Cardiol. 11, 255-265 (2014). https://doi.org/10.1038/nrcardio.2014.28
2. Prabhu, S. D. & Frangogiannis, N. G. The biological basis for cardiac repair after myocardial infarction: from inflammation to fibrosis. Circ. Res. 119, 91-112 (2016). https://doi.org/10.1161/CIRCRESAHA.116.303577
3. Frangogiannis, N. G. Cardiac fibrosis. Cardiovasc. Res. 117, 1450-1488 (2021). https://doi.org/10.1093/cvr/cvaa324
4. Ståhl, P. L. et al. Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. Science 353, 78-82 (2016). https://doi.org/10.1126/science.aaf2403
5. Kuppe, C. et al. Spatial multi-omic map of human myocardial infarction. Nature 608, 766-777 (2022). https://doi.org/10.1038/s41586-022-05060-x
6. Calcagno, D. M. et al. Single-cell and spatial transcriptomics of the infarcted heart define the dynamic onset of the border zone in response to mechanical destabilization. Nat. Cardiovasc. Res. 1, 1039-1055 (2022). https://doi.org/10.1038/s44161-022-00160-3
7. Yamada, S. et al. Spatiotemporal transcriptome analysis reveals critical roles for mechano-sensing genes at the border zone in remodeling after myocardial infarction. Nat. Cardiovasc. Res. 1, 1072-1083 (2022). https://doi.org/10.1038/s44161-022-00140-7
8. Ninh, V. K. et al. Spatially clustered type I interferon responses at injury borderzones. Nature 633, 174-181 (2024). https://doi.org/10.1038/s41586-024-07806-1
9. Farbehi, N. et al. Single-cell expression profiling reveals dynamic flux of cardiac stromal, vascular and immune cells in health and injury. eLife 8, e43882 (2019). https://doi.org/10.7554/eLife.43882
10. Amrute, J. M. et al. Targeting immune-fibroblast cell communication in heart failure. Nature 635, 423-433 (2024). https://doi.org/10.1038/s41586-024-08008-5
11. Alexanian, M. et al. Chromatin remodelling drives immune cell-fibroblast communication in heart failure. Nature 635, 434-443 (2024). https://doi.org/10.1038/s41586-024-08085-6
12. Ruiz-Villalba, A. et al. Single-cell RNA sequencing analysis reveals a crucial role for CTHRC1 (collagen triple helix repeat containing 1) cardiac fibroblasts after myocardial infarction. Circulation 142, 1831-1847 (2020). https://doi.org/10.1161/CIRCULATIONAHA.119.044557
13. Hernández, S. C. et al. Single-cell and spatial transcriptomic profiling of cardiac fibroblasts following myocardial infarction. Sci. Data 13, 216 (2026). https://doi.org/10.1038/s41597-025-06533-0
14. Oka, T. et al. Genetic manipulation of periostin expression reveals a role in cardiac hypertrophy and ventricular remodeling. Circ. Res. 101, 313-321 (2007). https://doi.org/10.1161/CIRCRESAHA.107.149047
15. Rahman, M. N. et al. SPaSE: spatially resolved pathology scores using optimal transport on spatial transcriptomics data. Cell Syst. 16, 101301 (2025). https://doi.org/10.1016/j.cels.2025.101301
16. Gene Expression Omnibus. GSE214611: Single-cell and spatial transcriptomics of the infarcted heart define the dynamic onset of the border zone in response to mechanical destabilization. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214611
