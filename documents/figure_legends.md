# BIB Figure Legends and Alt Text


## Figure 1. SSTBA workflow, input-output contract and diagnostic gates.

The protocol accepts a spatial expression matrix, frozen module definitions, coordinates, optional labels and biological-unit identifiers. It produces module or composite scores, coordinate- or graph-based boundary summaries, spot-level boundary indicators and a run manifest containing parameters, versions, hashes and warnings. Diagnostic gates cover signature coverage, identifier alignment, graph definition, biological replication, null models and baseline comparison. The compact real-data demonstration begins with released processed scores; raw-matrix reconstruction uses the documented analysis scripts.

Alt text: Workflow diagram from spatial expression, coordinates, labels and frozen signatures through scoring, boundary graphs and diagnostics to score tables, boundary tables and a checksum-bearing run manifest.

## Figure 2. Discovery sections resolve three spatial-state outputs.

a, Representative day 3 (D3_1) and day 7 (D7_1) sections show author domains and mechanical-border, immune-fibrotic activation and fibroblast-scar repair scores. Each score uses one colour scale fixed to its 2nd-98th percentile across all 14,147 discovery spots. b,c, Mean mechanical-border and fibroblast-scar repair scores in author domains 3 and 4 for all six GSE214611 sections. Thin lines are biological sections and thick lines are stage means. Domain 3 is mechanically enriched at both stages, whereas domain 4 shows stronger fibroblast-scar repair at day 7. Spots are displayed for spatial localization and are not treated as independent replicates.

Alt text: Representative day 3 and day 7 mouse-heart maps show author domains and three score fields; section-level paired plots below show higher mechanical-border in domain 3 and strong day 7 scar repair in domain 4.

## Figure 3. Coordinate maps and a fixed six-nearest graph quantify the domain 3-domain 4 interface.

a, Coordinate-defined boundary maps for the six discovery sections; contact uses 1.35 times the within-section median nearest-neighbour distance. b,c, Domain-specific boundary fractions and numbers of domain 3-domain 4 contact edges in the independently specified symmetrized six-nearest-neighbour graph. d, Mean domain 4-minus-domain 3 edge-local gradients summarized across three sections per stage; error bars are standard errors across sections. Mechanical gradients remain directed toward domain 3, while scar-repair gradients are directed toward domain 4 at day 7. These are sampled-stage descriptions, not longitudinal measurements.

Alt text: Six boundary maps and graph summaries show section-level domain 3-domain 4 contacts, persistent mechanical direction toward domain 3 and a day 7 scar-repair direction toward domain 4.

## Figure 4. Falsification diagnostics and simple baselines define the performance boundary.

a, Fraction of leave-one-section-out runs retaining each prespecified discovery direction. b, Complete-score directional margins and ranges after removing one component module. c, Observed margins compared with the 95th percentile of matched random signatures; the immune-fibrotic day 7 margin does not exceed this control. d, Rank-biserial effects for SSTBA and prespecified simple baselines across external comparison tasks. The failed GSE176092 whole-section mechanical comparison and baseline-equivalent scar tasks remain visible rather than being removed after analysis.

Alt text: Four diagnostic panels show leave-one-section-out and module-dropout stability, a negative immune-fibrotic matched-random result, and external tasks where SSTBA and simple baselines are sometimes equivalent.

## Figure 5. Frozen transfer across independent mouse studies and a nine-patient human panel.

a, Frozen-transfer design and datasets. b, GSE176092 spatial mechanical-minus-scar dominance across day 1, day 7 and day 14; the prespecified whole-section day 1-versus-day 14 mechanical comparison is not supported. c, GSE265828 author-labelled remote (RZ), border-zone (BZ1 and BZ2) and infarct-zone (IZ) summaries for mechanical-border, scar-repair and Boundary Transition Index outputs. d, Patient-section summaries across three control, three border-containing and three fibrotic Kuppe patients. Horizontal bars are group means. e, Prespecified rank-biserial comparisons with simple baselines; negative and baseline-equivalent results are retained. No target dataset was used to refit signatures or weights.

Alt text: Cross-study panels show mouse temporal and anatomical summaries plus patient-level human scores, including one failed mouse comparison, improved human mechanical consistency and baseline-equivalent fibrosis separation.

## Supplementary Figure S1. Sample-level signature maps.

Single-section triptych maps for each GSE214611 mouse myocardial infarction Visium replicate show mechanical-border, fibroblast-scar repair and composite fibrotic-risk scores. Colour scales are stated within the figure.

Alt text: Six mouse-heart sections each contain three spatial maps showing mechanical-border, fibroblast-scar repair and composite fibrotic-risk score distributions.

## Supplementary Figure S2. Domain 3-domain 4 boundary maps by sample.

Single-section maps show domain 3 and domain 4 boundary spots and non-boundary spots under the coordinate-defined contact rule.

Alt text: Six mouse-heart sections show domain 3 and domain 4 spots at their mutual spatial interface, with non-boundary tissue shown separately.

## Supplementary Figure S3. Public-data robustness and falsification summary.

Leave-one-section-out, module-dropout, boundary-direction and boundary-threshold analyses preserved 28/28, 15/15, 4/4 and 10/10 prespecified directions. The two mechanical-border checks and the day 7 scar-repair check exceeded the 99th percentile of 1,000 equal-size random signatures; immune-fibrotic activation remained contextual. Human score-separation analysis used 1,551 spots and showed distinct mechanical hotspots but correlated, partially overlapping stromal outputs.

Alt text: Robustness panels summarize sample omission, module dropout, boundary sensitivity, random-signature controls and score correlation in the single human feasibility section.

## Supplementary Figure S4. Spatial-statistical strengthening of the domain-state workflow.

a, Biological-section domain 3-minus-domain 4 margins with bootstrap confidence intervals. b, Moran's I across mouse sections; Geary's C is supplied in the source data. c, Domain-label and spatial-block permutations with fixed score fields. d, Signed distance-to-boundary gradients, positive toward domain 3. e, Removal sensitivity for POSTN, COL1A1, COL1A2, CTHRC1, NPPA, NPPB, XIRP2, SPP1 and LGALS3. These are spatial falsification checks, not experimental validation.

Alt text: Five panels show section-level score margins, spatial autocorrelation, permutation nulls, boundary-distance gradients and single-gene removal sensitivity.

## Supplementary Figure S5. Audits for domain dependence and cell-type composition.

a, Adjusted Rand index and normalized mutual information between k = 4 score-only states and author domains. b, Best-cluster F1 recovery of domains 3 and 4. c, Expected-direction margins before and after cell-type marker adjustment; mechanical contrasts remained positive, whereas stromal contrasts attenuated. d, Evidence-status calibration. These analyses are not single-cell deconvolution or experimental validation.

Alt text: Domain-recovery and marker-adjustment panels show partial label-free recovery, relatively stable mechanical margins and composition-sensitive stromal margins.

## Supplementary Figure S6. Provenance, scoring and graph-boundary checks.

a, Replicate quality control; complete file checksums, domain counts and module detectability are supplied in source data. b, Expected-direction margins across four scoring methods. c, Expression-matched random signatures: mechanical checks exceeded the 99th percentile and stromal checks the 95th percentile. d,e, Six-neighbour Visium graph boundary fractions and domain 4-minus-domain 3 edge gradients. f, Public-image H&E intensity audit, interpreted only as coarse image alignment rather than histological validation.

Alt text: Six panels summarize replicate quality, alternative scoring, matched-signature nulls, graph contacts, edge gradients and conservative public-image alignment checks.

## Supplementary Figure S7. Domain-independent graph states and human spatial-null checks.

a,b, Representative day 3 and day 7 k = 4 states after one-step graph smoothing without author labels. c, Mechanical and scar maxima shared a state in all day 3 sections but separated in all day 7 sections. d, Same-state edge fractions across four smoothing weights. e, Human graph Moran's I compared with 500 expression- and detection-matched signatures. These are internal falsification checks, not independent validation.

Alt text: Label-free graph-smoothed states remain partly intermixed at day 3 but separate at day 7, and human score fields exceed matched spatial-autocorrelation nulls.

## Supplementary Figure S8. Optional repair-aware hypothesis-prioritization example.

An optional downstream workflow combines section-level spatial evidence, single-section human feasibility support, public pharmacology annotation and an explicit repair, structural and safety deduction. Sensitivity panels show candidate ranks under section omission and weight changes and concordance among pharmacology sources. The output is an experimental hypothesis list and is not evidence of target causality, drug efficacy or treatment suitability.

Alt text: Supplementary plots show how spatial evidence and public pharmacology annotations can rank experimental hypotheses while deducting repair and structural concerns; no treatment effect is displayed.

## Supplementary Figure S9. Patient-level human transfer maps and cross-study sensitivity.

a, Boundary Transition Index maps for all nine Kuppe sections on a common colour scale after restriction to deposited in-tissue spots. b, Numbers of detected genes per transferred module in each patient. c, Leave-one-biological-unit-out direction consistency for SSTBA and simple baselines across six prespecified external comparisons. d, Numbers of in-tissue human spots by deposited region group; each group contains three independent patients.

Alt text: Nine patient maps and three summary panels show transferred human score orientation, module coverage, leave-one-patient-out consistency and deposited spot counts.
