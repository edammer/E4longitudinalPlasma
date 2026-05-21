# E4 Longitudinal Plasma Proteomics Pipeline

Analysis pipeline for longitudinal plasma SomaScan 7k proteomic trajectories comparing APOE e4/e4 and APOE e3/e3 individuals across estimated years from onset (EYO). The scripts in this repository are ordered by numeric prefix and represent the manuscript analysis workflow from GNPC Harmonized Dataset extraction through proteomic harmonization, APOE genotype imputation, WGCNA network analysis, Bayesian longitudinal modeling, ontology-level trajectory summarization, CMAP-style drug connectivity analysis, and MAGMA/module enrichment testing.

## Pipeline overview by numeric prefix

| Prefix | Major purpose | Major outputs |
|---|---|---|
| `1` | Load GNPC HDS v1.3 clinical and SomaScan tables, curate traits, construct plasma abundance matrices, remove non-sample controls/calibrators, harmonize APOE and diagnosis traits, and perform two-stage proteomic regression/QC for pre-analytical and site effects. | Curated trait tables, SomaScan assay metadata with unique assay IDs, cleaned plasma abundance matrices, first- and second-pass regressed matrices, t-SNE/QC plots, variance-partition/QC summaries, volcano plots, and saved workspace objects for downstream modeling. |
| `2` | Train APOE genotype imputation models from proteomic features using binary one-vs-all learners and a final ordered six-genotype ensemble. | Ranked genotype-predictive protein features, cross-validation metrics, binary model predictions, six-class APOE genotype predictions, and saved genotype-imputed workspace objects. |
| `3` / `4` | Use known and predicted APOE information to construct final regressed plasma matrices, perform WGCNA network analysis, define plasma protein modules, and generate final network/QC summaries. | Final protein abundance matrix and metadata, WGCNA module color assignments, module eigengenes/eigenproteins, module-trait/QC plots, variance-partition outputs, differential abundance volcanoes, and final saved network/genotype workspace. |
| `5` | Fit and visualize Bayesian longitudinal restricted cubic spline models for APOE e4/e4 vs e3/e3 trajectories across EYO, including sensitivity cohorts, sex-adjusted and sex-interaction models, and post-hoc posterior summaries. | Cleaned longitudinal STAN input matrices, per-assay `stan_glm` RDS model fits, posterior trajectory plots, median and 99% credible interval summaries, posterior p-value/effect-size matrices across EYO, category labels for significant genotype × sex × EYO interaction intervals, and single-panel interaction plots. |
| `6` | Extract EYO-window-specific significant assays from posterior trajectory results and organize direction-specific signatures for ontology trajectory/waterfall and GSVA-style visualization. | Lists of significant up/down assays by EYO window, assay-count summaries, ontology/waterfall inputs, and heatmap-ready tables for downstream pathway and CMAP analyses. |
| `7` | Compute CMAP-style connectivity between e4/e4 trajectory signatures and semaglutide perturbation profiles within ontology-defined assay universes, including permutation-based significance testing and figure-ready heatmaps. | Observed C-score/WTCS/NCS matrices, permutation p-values and null summaries, ontology × EYO heatmaps, and figure-ready CMAP permutation workspaces/plots. |
| `8` | Run MAGMA/Seyfried Pipeline Adaptation enrichment tests using gene/protein list input vectors against module- or pathway-level gene sets. | MAGMA enrichment summary tables, permutation-adjusted enrichment outputs, and bar plots of enrichment across supplied input gene lists. |

## Numerically prefixed R scripts

### `1a.HDS-load_merge_clean-v1.3ms_accessed_03-27-2025.R`
Loads GNPC HDS v1.3 tables from the ADDI PostgreSQL environment, exports clinical and SomaLogic analyte tables, builds stable assay identifiers from gene symbol/UniProt/SomaID/aptamer information, handles non-protein controls and technical replicate assay IDs, and prepares raw objects used by trait and expression cleanup.

Required input is read from the GNPC Harmonized Data Set v1.3 via Postgre SQL, and cannot be executed outside the ADDI Azure-enabled Virtual Machine.
- Authorized ADDI/GNPC PostgreSQL ODBC connection named `PostgreSQL`.
- GNPC HDS tables including `ClinicalV1_3ms` and `SomalogicAnalyteInfoV1_3ms`.
- R packages/helpers: `DBI`, `odbc`, `xaputils`, `data.table`.

Non-RData/RDS files or helper code to place in `input/`: none identified; the primary source is the ADDI database connection.

### `1b.TraitClean-v1.3ms_03-27-2025.R`
Cleans and harmonizes clinical/trait metadata, repairs APOE coding, removes duplicate join-key columns, converts MoCA/MMSE values using a published crosswalk, merges cohort-specific mapping data, constructs diagnosis/control fields, and prepares curated metadata for plasma and downstream analysis.

Required input is read from the GNPC Harmonized Data Set v1.3 via Postgre SQL, and cannot be executed outside the ADDI Azure-enabled Virtual Machine.
- `loadedV1_3ms_03-27-25.RData` or equivalent workspace generated by prefix `1a`.
- Cohort-specific mapping RDS files referenced in the script, including `BH.map.RDS`, `RM.map.RDS`, `UDS.map.RDS`, and ANML matrices.

Other input files or helper code to place in `./input/`:
- `9b.4cohort.csv`

### `1c.4p13.2xRegr_1stRegr.2protPAVwithinSite_2ndReg(Site-3ways_orMedianZero)+QC_VP+tSNE+Volc_b345regrFIXED.R`
Filters to plasma SomaScan samples, removes calibrators, adds pre-analytical proxy assays, maps SomaSignal pre-analytical traits, performs t-SNE QC, segments site F into sub-batches, performs first-pass within-site regression of HNRNPA2B1/HBZ pre-analytical effects, performs second-pass site/subsite regression while protecting age/sex/APOE e4 effects, generates QC/variance-partition and volcano outputs, and prepares APOE prediction and final harmonized matrices.

Required input is read from the GNPC Harmonized Data Set v1.3 via Postgre SQL, and cannot be executed outside the ADDI Azure-enabled Virtual Machine.
- `2.saved.image_trait+human_cleanup_nm+em0_V1_3ms_03-27-25.RData` metadata/proteomics workspace generated by code in prefix `1` exists only on the ADDI VM shared drive.
- This is also the case for optional intermediate RDS/RData files referenced for checkpointing or re-entry, such as `4p13.cleanDat.22sites.RDS`, `4p13.numericMeta.22sites.RDS`, and APOE prediction RDS files.

Other input files or helper code (`./input/`):
- `BH-SomaSIgnals_forR_andHDS_BH.map.csv`
- `ROSMAP-SomaSignalsForR_andHDS_RM.map.csv`
- `UDS-SomaSignalsForR_andHDS_UDS.map.csv`
- `../samePage.aheatmap.below.R`
- `buildIgraphs.R`
- `GOparallel-FET.R`
- `geneListFET.R`
- `geneListFET_customLabels.R`
- `geneListFET_customLabels-fixedScale+15thPlotNominalP.R`
- `parANOVA.dex.R`
- `parANOVA.dex.fallback7.25.R`
- `plasma_21cohort_17moduleOntologies_unregr+color.csv`
- `RNAbindingProtein.FET.lists.csv`

### `2a.APOE_prediction-stage_1.R`
Trains one-vs-all proteomic APOE binary classifiers for each APOE genotype, using cross-validation and ensemble learners to rank genotype-predictive proteins and estimate binary classifier performance.

Required inputs:
- `4p13b4b.SecondRegressionsComplete.19sites_Fsplit+APOEgenoIMPUTATION.RData` or equivalent second-regressed matrix/metadata workspace.
- Optional prior outputs if re-entering the script: `rankedProteins.list.RDS`, `binaryPredictionMetrics.list.RDS`.

### `2b.APOE_prediction-stage_2.R`
Assembles the stage-2 APOE six-genotype ensemble using ranked proteins from stage 1, trains genotype-specific binary learners, applies ordered genotype decision logic, and generates final APOE predicted genotype calls and probabilities.

Required inputs:
- `4p13b4b.SecondRegressionsComplete.19sites_Fsplit+APOEgenoIMPUTATION.RData`.
- Stage-1 objects in the active workspace, especially `rankedProteins.prior` or equivalent ranked genotype feature lists.

### `3+4.(2ndRegr.2)4p13c1.PLASMA_known+predictedAPOE_e4carrierRegressedSite(secondRegr)-WGCNA.R`
Builds the final plasma proteomic dataset using known and imputed APOE genotype information, performs second-pass regression with APOE e4 carrier status protected, checks residual APOE signal, performs WGCNA network construction, refines modules, calculates module eigengenes/eigenproteins, and generates network QC and enrichment plots.

Required inputs:
- `saved.image-genotype_prediction_finalized.RData` or equivalent output from prefixes `1` and `2`.
- Optional intermediate matrix/metadata RDS files used for re-entry or comparison.

Non-RData/RDS files or helper code to place in `input/`:
- `../samePage.aheatmap.below.R`
- `buildIgraphs.R`
- `GOparallel-FET.R`
- `geneListFET.R`
- `geneListFET_customLabels.R`
- `geneListFET_customLabels-fixedScale+15thPlotNominalP.R`
- `parANOVA.dex.R`
- `parANOVA.dex.fallback7.25.R`
- `plasma_21cohort_17moduleOntologies_unregr+color.csv`
- `RNAbindingProtein.FET.lists.csv`

### `3b.redo VariancePartition and final QC volcanoes.R`
Reloads the finalized genotype-prediction workspace, recomputes variance-partition models on the final regressed matrix, and regenerates final QC volcano plots for age, sex, APOE e4 carrier status, and AD/control contrasts.

Required inputs:
- `saved.image-genotype_prediction_finalized.RData`.
- Volcano/differential-expression CSVs when replotting from saved summaries.

Non-RData/RDS files or helper code to place in `input/`:
- `4p13c1.QC..ANOVA_diffEx-ALL-AD.CT_volc-ANOVA_tTest.CTimputed.csv`
- `4p13c1.QC..ANOVA_diffEx-ALL-Sex_volc-ANOVA_tTest.csv`
- `4p13c1.QC..ANOVA_diffEx-ALL-e4carrier_volc-ANOVA_tTest.csv`
- `GOparallel-FET.R`
- `parANOVA.dex.R`
- `parANOVA.dex.fallback7.25.R`

### `5a.model_1_STAN_parallel.+clean(3177 sample homozygote contrast).R`
Creates the primary 3,177-sample APOE homozygote analysis dataset, applies ±5 SD outlier masking and group-wise missingness filtering, adds MMSE/CDR outcomes and module eigengenes, defines EYO and e4/e4 indicator variables, sanitizes outcome names, and fits per-outcome Bayesian restricted cubic spline models for e4/e4 vs e3/e3 trajectories in parallel.

Required inputs:
- `APOE_homozygote_2ndFinal_cleanDats+traits(3177.CT33only_4199.33inADcontinuum).RData` or equivalent containing `cleanDat.3177`, `MEs.3177`, and `numericMeta.3177`.
- Optional preprocessed trait/protein RDS files if using the commented alternate paths.

### `5a2.Plot_STAN_model_1_20250718-direct2PDF+parallel(3177) _5xSDcleaned-0727.R`
Reads Model 1 STAN fits for the 3,177-sample analysis and generates posterior trajectory plots for e4/e4 and e3/e3 groups, including difference curves across EYO and model-derived summaries.

Required inputs:
- `_numericMeta_3177_trait.RDS`
- `_full_3177_protein_dft.RDS`
- `simple.3177/name_match_table.rds`
- Per-assay Model 1 STAN RDS files in `simple.3177/`.

Non-RData/RDS files or helper code to place in `input/`:
- `plot_functions_20250718.R`
- `scatterplot_label_20250718.csv`

### `5a3.STAN_Plot_medianExtract_e44-e33.R`
Extracts posterior draws from Model 1 STAN fits, computes posterior medians and credible intervals for e3/e3 and e4/e4 trajectories and their differences across EYO, and writes matrix outputs used by downstream significance-window and ontology analyses.

Required inputs:
- `_numericMeta_3177_trait.RDS`
- `_full_3177_protein_dft.RDS`
- `simple.3177/name_match_table.rds`
- Per-assay Model 1 STAN RDS files in `simple.3177/`.

Non-RData/RDS files or helper code to place in `input/`:
- `plot_functions_20250718.R`
- `scatterplot_label_20250718.csv`

### `5b.model_1_STAN_parallel.+clean(4199 samples-with non-control e3 homozygotes).R`
Repeats Model 1 longitudinal spline fitting in the broader 4,199-sample sensitivity cohort that includes non-control e3/e3 individuals along the AD continuum, with the same outlier and missingness rules used for the primary cohort.

Required inputs:
- `APOE_homozygote_2ndFinal_cleanDats+traits(3177.CT33only_4199.33inADcontinuum).RData` or equivalent containing `cleanDat.4199`, `MEs.4199`, and `numericMeta.4199`.

### `5c.model_2_withSex_STAN_parallel+clean(3177).R`
Fits Model 2 in the primary 3,177-sample cohort by adding sex-related terms to the APOE/EYO spline model, after applying the same protein/outcome cleaning, EYO derivation, and e4/e4 indicator setup.

Required inputs:
- `APOE_homozygote_2ndFinal_cleanDats+traits(3177.CT33only_4199.33inADcontinuum).RData` or equivalent containing `cleanDat.3177`, `MEs.3177`, and `numericMeta.3177`.

### `5c2.Plot_STAN_model_2_20250718-direct2PDF+parallel(3177)_5xSDcleaned-0727+plotFixes.R`
Plots sex-adjusted Model 2 trajectories and contrasts across EYO, verifies expected coefficient order, and produces posterior summaries and visualization outputs with sex coding handled explicitly.

Required inputs:
- `_numericMeta_3177_trait.RDS`
- `_full_3177_protein_dft.RDS`
- `sexInt.3177/name_match_table.rds`
- Per-assay Model 2 STAN RDS files in `sexInt.3177/`.

Non-RData/RDS files or helper code to place in `input/`:
- `plot_functions_20250718.R`
- `scatterplot_label_20250718.csv`

### `5d.model_2_withSex_STAN_parallel+clean(4199 samples including non-control e3 homozygotes).R`
Repeats sex-adjusted Model 2 fitting in the broader 4,199-sample sensitivity cohort including non-control e3/e3 individuals.

Required inputs:
- `APOE_homozygote_2ndFinal_cleanDats+traits(3177.CT33only_4199.33inADcontinuum).RData` or equivalent containing `cleanDat.4199`, `MEs.4199`, and `numericMeta.4199`.

### `5e.model_1_STAN_without_Imputed_Genotype_Samples_5xSDcleaned(RDS)-0727.R`
Runs a sensitivity analysis excluding samples without mapped APOE genotype, then refits Model 1 using only known genotype samples after alignment checks, outlier masking, and EYO/e4/e4 indicator construction.

Required inputs:
- `_numericMeta_3177_trait.RDS`
- `_full_3177_protein_dft.RDS`
- `noImpute.simple.rerun/name_match_table.RDS`

### `5e2.Plot_STAN_model_1_withoutImputedGenotypeSamples-direct2PDF+parallel(3177) _5xSDcleaned-0727.R`
Plots the no-imputed-genotype Model 1 sensitivity analysis and extracts trajectory/difference summaries for known-genotype-only e4/e4 vs e3/e3 comparisons.

Required inputs:
- `assays+traits_forPlots.RData`
- `noImpute.simple.rerun/name_match_table.rds`
- Per-assay no-imputation Model 1 STAN RDS files in `noImpute.simple.rerun/`.

Non-RData/RDS files or helper code to place in `input/`:
- `plot_functions_20250718.R`
- `scatterplot_label_20250718.csv`

### `5f.model_3_STAN_withSex+Genotype+EYO_3way_int_terms_5xSDcleaned(RDS)-0727.R`
Fits Model 3 / Model 2.3 in the primary cohort, adding APOE genotype × sex × EYO spline interaction terms so that sex-specific APOE e4/e4 and e3/e3 longitudinal differences can be evaluated directly.

Required inputs:
- `_numericMeta_3177_trait.RDS`
- `_full_3177_protein_dft.RDS`

### `5f3.model_3_STAN_OutputSummaryStats99CI-5xSDcleaned(RDS)-0727.R`
Performs post-hoc posterior testing and summary extraction for Model 3, computing EYO-grid posterior means, credible intervals, and two-tailed Bayesian posterior probabilities for genotype × sex × EYO interaction terms, with corrected sex coding.

Required inputs:
- `_numericMeta_3177_trait.RDS`
- `name_match_table.RDS`
- Per-assay Model 3 STAN RDS files in the configured Model 3 directory.

### `5f4.Plot_STAN_model_3_99CI_withSex+Genotype+EYO_3way_int_terms-5xSDcleaned(RDS)-0727.R`
Generates single-panel Model 3 plots showing genotype × sex interaction posterior ribbons across EYO and a colored significance streak indicating the dominant sex/genotype category at significant EYO bins.

Required inputs:
- `_numericMeta_3177_trait.RDS`
- `name_match_table.RDS`
- Per-assay Model 3 STAN RDS files and post-hoc output summaries from `5f3`.

### `6.ExtractSigAssays_2directions(ALL-up+down)_z=1.645.R`
Uses posterior p-value and effect-size matrices from Model 1 to collect significant assays in sliding 5-year EYO windows, split signatures by direction, create assay-count summaries, organize ontology trajectory inputs, and support heatmap/waterfall visualization.

Required inputs:
- `../simple.3177/scatter/_99_par_diff_all_peptide_p_value.rds`
- `../simple.3177/scatter/_99_par_diff_all_peptide_up_down_notation.csv`
- `_numericMeta_3177_trait.RDS`
- `_full_3177_protein_dft.RDS`
- `name_match_table.rds`
- GO/ontology inputs such as `genelists.GO.forGSVA.RDS` and `go.obo` if running the GSVA/ontology portions.

Non-RData/RDS files or helper code to place in `input/`:
- `plot_functions_20250718.R`
- `scatterplot_label_20250718.csv`
- `ALL141terms_18categories.tsv`
- `../simple.3177/scatter/_99_par_diff_all_peptide_up_down_notation.csv`
- `go.obo`

### `7a.CMAP.permutation100000.R`
Computes ontology-restricted CMAP-style connectivity between rare/APOE e4/e4 longitudinal up/down signatures and semaglutide perturbation statistics, then performs permutation testing for each ontology × EYO epoch cell while holding the drug ranking fixed.

Required inputs:
- `../simple.3177/scatter/_99_par_diff_all_peptide_p_value.rds`
- `saved.image-SemaS6_minESsweep.RData` or equivalent upstream CMAP sweep workspace.

Non-RData/RDS files or helper code to place in `input/`:
- `SemaglutideStudy_TableS2stats.csv`
- `SemaglutideStudy_TableS6stats.csv`
- `ALL_heatmap18categoryROWorder-dataFrame.tsv`
- `Fig6A_ontology_order_110_18categories.tsv`
- `Fig6A_ontology_order_141.tsv`

### `7b.CMAP.permutation_plots.R`
Loads CMAP permutation results and renders category-ordered ComplexHeatmap figures using a fixed 18-category ontology order and figure-ready ontology/category annotation tables.

Required inputs:
- `saved.image-CMAP.perm.RData`

Non-RData/RDS files or helper code to place in `input/`:
- `Fig6A_ontology_order_110_18categories_final.tsv`
- `Fig6A_ontology_order_141.tsv`

### `8.MAGMA_wrapper.R`
Defines parameters and calls `MAGMA.SPA()` for gene/protein set enrichment analysis using semaglutide or other MAGMA-style input files and module/gene-list vectors.

Required inputs:
- Input MAGMA CSV files in `MAGMAinputDir`, including the configured examples `Sema_S6_NominalP.csv` and `Sema_S6_Qvalue.csv`.
- Required global objects expected by `MAGMA.SPA()`, especially `moduleGeneList` or module/protein list equivalents.

Non-RData/RDS files or helper code to place in `input/`:
- `MAGMA.SPA_listVectorInput.R`
- `Sema_S6_NominalP.csv`
- `Sema_S6_Qvalue.csv`

## Suggested `input/` folder contents

The scripts intentionally retain original absolute paths for provenance. For repo portability, collect helper scripts and small tabular inputs into an `input/` subfolder and update paths accordingly. Large RData/RDS workspaces are not listed here unless they are small enough and appropriate for controlled sharing.

Recommended local files to collect:

```text
input/9b.4cohort.csv
input/BH-SomaSIgnals_forR_andHDS_BH.map.csv
input/ROSMAP-SomaSignalsForR_andHDS_RM.map.csv
input/UDS-SomaSignalsForR_andHDS_UDS.map.csv
input/2.Ensemble-13cellTypes-bulkRNA_Raj59-proportionEstimates.csv
input/samePage.aheatmap.below.R
input/buildIgraphs.R
input/GOparallel-FET.R
input/geneListFET.R
input/geneListFET_customLabels.R
input/geneListFET_customLabels-fixedScale+15thPlotNominalP.R
input/parANOVA.dex.R
input/parANOVA.dex.fallback7.25.R
input/plasma_21cohort_17moduleOntologies_unregr+color.csv
input/RNAbindingProtein.FET.lists.csv
input/4p13c1.QC..ANOVA_diffEx-ALL-AD.CT_volc-ANOVA_tTest.CTimputed.csv
input/4p13c1.QC..ANOVA_diffEx-ALL-Sex_volc-ANOVA_tTest.csv
input/4p13c1.QC..ANOVA_diffEx-ALL-e4carrier_volc-ANOVA_tTest.csv
input/plot_functions_20250718.R
input/scatterplot_label_20250718.csv
input/_99_par_diff_all_peptide_up_down_notation.csv
input/ALL141terms_18categories.tsv
input/go.obo
input/SemaglutideStudy_TableS2stats.csv
input/SemaglutideStudy_TableS6stats.csv
input/ALL_heatmap18categoryROWorder-dataFrame.tsv
input/Fig6A_ontology_order_110_18categories.tsv
input/Fig6A_ontology_order_141.tsv
input/Fig6A_ontology_order_110_18categories_final.tsv
input/MAGMA.SPA_listVectorInput.R
input/Sema_S6_NominalP.csv
input/Sema_S6_Qvalue.csv
```

## Notes for use
- GNPC harmonized data set (HDS) v1.3 data is available on the ADDI Azure-implemented VM platform for authorized users.
- These scripts remain research pipeline scripts rather than package-style functions; paths and large workspace objects are intentionally left visible for provenance.
- Before rerunning on a different system, update `rootdir`/`setwd()` paths and confirm sourced helper files are available.

## Repository status

The repository currently contains numerically prefixed R scripts from prefixes `1` through `8`, including the prefix-5 STAN modeling and plotting scripts.
