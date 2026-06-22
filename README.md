# E4 Longitudinal Plasma Proteomics Pipeline

Analysis pipeline for longitudinal plasma SomaScan 7k proteomic trajectories comparing APOE e4/e4 and APOE e3/e3 individuals across estimated years from onset (EYO). The scripts in this repository are ordered by numeric prefix and represent the manuscript analysis workflow from GNPC Harmonized Dataset extraction through proteomic harmonization, APOE genotype imputation, WGCNA network analysis, Bayesian longitudinal modeling, ontology-level trajectory summarization, endophenotype/pathology overlap testing, CMAP-style drug connectivity analysis, and MAGMA/module enrichment testing.

## Pipeline overview by numeric prefix

| Prefix | Major purpose | Major outputs |
|---|---|---|
| `1` | Load GNPC HDS v1.3 clinical and SomaScan tables, curate traits, construct plasma abundance matrices, remove non-sample controls/calibrators, harmonize APOE and diagnosis traits, and perform two-stage proteomic regression/QC for pre-analytical and site effects. | Curated trait tables, SomaScan assay metadata with unique assay IDs, cleaned plasma abundance matrices, first- and second-pass regressed matrices, t-SNE/QC plots, variance-partition/QC summaries, volcano plots, and saved workspace objects for downstream modeling. |
| `2` | Train APOE genotype imputation models from proteomic features using binary one-vs-all learners and a final ordered six-genotype ensemble. | Ranked genotype-predictive protein features, cross-validation metrics, binary model predictions, six-class APOE genotype predictions, and saved genotype-imputed workspace objects. |
| `3` / `4` | Use known and predicted APOE information to construct final regressed plasma matrices, perform WGCNA network analysis, define plasma protein modules, and generate final network/QC summaries. | Final protein abundance matrix and metadata, WGCNA module color assignments, module eigengenes/eigenproteins, module-trait/QC plots, variance-partition outputs, differential abundance volcanoes, and final saved network/genotype workspace. |
| `5` | Fit and visualize Bayesian longitudinal restricted cubic spline models for APOE e4/e4 vs e3/e3 trajectories across EYO, including sensitivity cohorts, sex-adjusted and sex-interaction models, and post-hoc posterior summaries. | Cleaned longitudinal STAN input matrices, per-assay `stan_glm` RDS model fits, posterior trajectory plots, median and 99% credible interval summaries, posterior p-value/effect-size matrices across EYO, category labels for significant genotype × sex × EYO interaction intervals, and single-panel interaction plots. |
| `6` | Test whether sliding-window and before/after-age-50 significant-assay lists are enriched for proteins previously reported as pathology- or cognition-associated endophenotype markers, using one-sided Fisher's Exact Test (FET) overlap enrichment. | Cleaned (non-protein-assay-free) sliding-window hit list, and per-comparison FET heatmap PDFs + statistics tables/CSVs for overlap with Nat Aging pathology and cognitive-function marker lists. |
| `7` | Extract EYO-window-specific significant assays from posterior trajectory results (combined ALL = up+down, UP only, and DOWN only), organize direction-specific signatures, run ontology (GO/Reactome/WikiPathways/MSigDB) hypergeometric enrichment per window, and summarize/plot per-ontology Z-score trajectories and curated 18-category publication heatmaps. | Lists of significant up/down assays by EYO window, assay-count summary plots with WGCNA module overlay, per-ontology Z-score trajectory and heatmap PDFs, curated 18-category/141-term publication heatmaps (ALL/UP/DOWN), and 141×67 Z-score/genes-hit supplementary tables. |
| `8` | Run CMAP/L1000-style permutation testing of connectivity between e4/e4 trajectory signatures (within 141 curated ontology terms) and semaglutide perturbation profiles across 3 EYO epochs, including a 21-point minimum-effect-size sensitivity sweep and figure-ready heatmaps. | Observed Cscore/NCS matrices, two-tailed/directional permutation p-values and null summaries, minES sensitivity-sweep summary tables/plots, the final 18-category × 3-epoch connectivity heatmap with significance stars, and a permutation-statistics Excel workbook. |
| `9` | Run MAGMA/Seyfried Pipeline Adaptation enrichment tests using gene/protein list input vectors against module- or pathway-level gene sets. | MAGMA enrichment summary tables, permutation-adjusted enrichment outputs, and bar plots of enrichment across supplied input gene lists. |

## Numerically prefixed R scripts

### `1a.HDS-load_merge_clean-v1.3ms_accessed_03-27-2025.R`
Loads GNPC HDS v1.3 tables from the ADDI PostgreSQL environment, exports clinical and SomaLogic analyte tables, builds stable assay identifiers from gene symbol/UniProt/SomaID/aptamer information, handles non-protein controls and technical replicate assay IDs, and prepares raw objects used by trait and expression cleanup.

Required input is read from the GNPC Harmonized Data Set v1.3 via Postgre SQL, and cannot be executed outside the ADDI Azure-enabled Virtual Machine.
- Authorized ADDI/GNPC PostgreSQL ODBC connection named `PostgreSQL`.
- GNPC HDS tables including `ClinicalV1_3ms` and `SomalogicAnalyteInfoV1_3ms`.
- R packages/helpers: `DBI`, `odbc`, `xaputils`, `data.table`.

Non-RData/RDS files or helper code to place in `input/`: none identified; the primary source is the ADDI database connection.

Major outputs: raw merged clinical + SomaLogic analyte/expression objects and a saved workspace (e.g. `loadedV1_3ms_03-27-25.RData`) consumed by prefix `1b`/`1c`.

### `1b.TraitClean-v1.3ms_03-27-2025.R`
Cleans and harmonizes clinical/trait metadata, repairs APOE coding, removes duplicate join-key columns, converts MoCA/MMSE values using a published crosswalk, merges cohort-specific mapping data, constructs diagnosis/control fields, and prepares curated metadata for plasma and downstream analysis.

Required input is read from the GNPC Harmonized Data Set v1.3 via Postgre SQL, and cannot be executed outside the ADDI Azure-enabled Virtual Machine.
- `loadedV1_3ms_03-27-25.RData` or equivalent workspace generated by prefix `1a`.
- Cohort-specific mapping RDS files referenced in the script, including `BH.map.RDS`, `RM.map.RDS`, `UDS.map.RDS`, and ANML matrices.

Other input files or helper code to place in `./input/`:
- `9b.4cohort.csv`

Major outputs: curated/harmonized clinical trait table (APOE, diagnosis/control, MoCA/MMSE-crosswalked cognition fields) and saved workspace consumed by prefix `1c`.

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

Major outputs: pre-analytical-proxy-regressed and site/subsite-regressed plasma abundance matrices (first- and second-pass), t-SNE QC plots, variance-partition summaries, and volcano plots used to QC the harmonized matrix before APOE imputation and WGCNA.

### `2a.APOE_prediction-stage_1.R`
Trains one-vs-all proteomic APOE binary classifiers for each APOE genotype, using cross-validation and ensemble learners to rank genotype-predictive proteins and estimate binary classifier performance.

Required inputs:
- `4p13b4b.SecondRegressionsComplete.19sites_Fsplit+APOEgenoIMPUTATION.RData` or equivalent second-regressed matrix/metadata workspace.
- Optional prior outputs if re-entering the script: `rankedProteins.list.RDS`, `binaryPredictionMetrics.list.RDS`.

Major outputs: ranked genotype-predictive protein feature lists (`rankedProteins.list.RDS`) and binary one-vs-all classifier predictions/performance metrics (`binaryPredictionMetrics.list.RDS`) consumed by stage 2 (`2b`).

### `2b.APOE_prediction-stage_2.R`
Assembles the stage-2 APOE six-genotype ensemble using ranked proteins from stage 1, trains genotype-specific binary learners, applies ordered genotype decision logic, and generates final APOE predicted genotype calls and probabilities.

Required inputs:
- `4p13b4b.SecondRegressionsComplete.19sites_Fsplit+APOEgenoIMPUTATION.RData`.
- Stage-1 objects in the active workspace, especially `rankedProteins.prior` or equivalent ranked genotype feature lists.

Major outputs: final six-class APOE predicted genotype calls and class probabilities per sample, and a saved genotype-imputed workspace (e.g. `saved.image-genotype_prediction_finalized.RData`) consumed by prefix `3`/`4`.

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

Major outputs: final known+imputed-APOE-regressed plasma protein abundance matrix and metadata, WGCNA module color assignments, module eigengene/eigenprotein tables, module-trait association and network QC plots, and a final saved network/genotype workspace consumed by prefix `5`.

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

Major outputs: refreshed variance-partition model summaries and final QC volcano plots (age, sex, APOE e4 carrier, AD/control) for the manuscript's QC figures.

### `5a.model_1_STAN_parallel.+clean(3177 sample homozygote contrast).R`
Creates the primary 3,177-sample APOE homozygote analysis dataset, applies ±5 SD outlier masking and group-wise missingness filtering, adds MMSE/CDR outcomes and module eigengenes, defines EYO and e4/e4 indicator variables, sanitizes outcome names, and fits per-outcome Bayesian restricted cubic spline models for e4/e4 vs e3/e3 trajectories in parallel.

Required inputs:
- `APOE_homozygote_2ndFinal_cleanDats+traits(3177.CT33only_4199.33inADcontinuum).RData` or equivalent containing `cleanDat.3177`, `MEs.3177`, and `numericMeta.3177`.
- Optional preprocessed trait/protein RDS files if using the commented alternate paths.

Major outputs: per-assay `<assay>_stan_glm.rds` Bayesian model fits for the primary 3,177-sample cohort, plus cleaned trait/protein RDS objects (`_numericMeta_3177_trait.RDS`, `_full_3177_protein_dft.RDS`) consumed by `5a2`/`5a3`/`5e`/`5e2`.

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

Major outputs: per-assay scatter+ribbon/difference trajectory PDFs, and the posterior t-statistic / p-value / up-down-direction matrices (`_99_par_diff_all_peptide*.rds`/`.csv`) consumed by the prefix-`7` ExtractSigAssays scripts.

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

Major outputs: posterior median/credible-interval trajectory matrices for e3/e3, e4/e4, and their difference, across the EYO grid, for downstream significance-window/ontology analyses.

### `5b.model_1_STAN_parallel.+clean(4199 samples-with non-control e3 homozygotes).R`
Repeats Model 1 longitudinal spline fitting in the broader 4,199-sample sensitivity cohort that includes non-control e3/e3 individuals along the AD continuum, with the same outlier and missingness rules used for the primary cohort.

Required inputs:
- `APOE_homozygote_2ndFinal_cleanDats+traits(3177.CT33only_4199.33inADcontinuum).RData` or equivalent containing `cleanDat.4199`, `MEs.4199`, and `numericMeta.4199`.

Major outputs: per-assay Model 1 `stan_glm` fits for the 4,199-sample (broader, non-control-e3/e3-included) sensitivity cohort.

### `5c.model_2_withSex_STAN_parallel+clean(3177).R`
Fits Model 2 in the primary 3,177-sample cohort by adding sex-related terms to the APOE/EYO spline model, after applying the same protein/outcome cleaning, EYO derivation, and e4/e4 indicator setup.

Required inputs:
- `APOE_homozygote_2ndFinal_cleanDats+traits(3177.CT33only_4199.33inADcontinuum).RData` or equivalent containing `cleanDat.3177`, `MEs.3177`, and `numericMeta.3177`.

Major outputs: per-assay sex-adjusted Model 2 `stan_glm` fits in `sexInt.3177/` for the primary cohort.

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

Major outputs: per-assay sex-stratified trajectory/difference PDFs and posterior summary matrices for Model 2.

### `5d.model_2_withSex_STAN_parallel+clean(4199 samples including non-control e3 homozygotes).R`
Repeats sex-adjusted Model 2 fitting in the broader 4,199-sample sensitivity cohort including non-control e3/e3 individuals.

Required inputs:
- `APOE_homozygote_2ndFinal_cleanDats+traits(3177.CT33only_4199.33inADcontinuum).RData` or equivalent containing `cleanDat.4199`, `MEs.4199`, and `numericMeta.4199`.

Major outputs: per-assay sex-adjusted Model 2 `stan_glm` fits for the 4,199-sample sensitivity cohort.

### `5e.model_1_STAN_without_Imputed_Genotype_Samples_5xSDcleaned(RDS)-0727.R`  *(patched)*
Runs a sensitivity analysis excluding samples without mapped (known) APOE genotype, then refits Model 1 using only known-genotype samples. **Patched in this update**: restricted cubic spline (RCS) knot x-positions (EYO) are no longer computed from the known-genotype-only subset; instead, for each assay, knots are computed once from that assay's complete cases in the *full* 3,177-sample cohort (including imputed-genotype samples), and those fixed knot positions are then used to fit the model on the known-genotype-only data. This keeps the spline basis identical to the original full-cohort Model 1 run so the two model sets are directly comparable; the full-cohort data are used only to pick knot locations and never enter the regression itself.

Required inputs:
- `_numericMeta_3177_trait.RDS` (full 3,177-sample trait table; also used unfiltered for knot estimation)
- `_full_3177_protein_dft.RDS` (full 3,177-sample protein/outcome matrix; also used unfiltered for knot estimation)
- `./name_match_table.RDS` (prior OriginalName/CleanedName lookup, refreshed and re-saved)

Major outputs: per-assay `<assay>_stan_glm.rds` Bayesian model fits (each carrying the full-cohort spline-knot positions used to fit it) for the known-genotype-only cohort, a refreshed `./name_match_table.RDS`, and `./assays+traits_forPlots.RData` (filtered + full-cohort trait/protein objects) consumed by `5e2`.

### `5e2.Plot_STAN_model_1_withoutImputedGenotypeSamples-direct2PDF+parallel(3177) _5xSDcleaned-0727.R`  *(patched)*
Plots the no-imputed-genotype Model 1 sensitivity analysis and extracts trajectory/difference summaries for known-genotype-only e4/e4 vs e3/e3 comparisons. **Patched in this update** to match `5e`: posterior trajectories are evaluated, and ribbons drawn, using the same full-3,177-sample (imputed-genotype-included) spline knot positions that were used to fit each model in `5e`, rather than knots recomputed from the known-genotype-only plotted points.

Required inputs:
- `./assays+traits_forPlots.RData` (filtered + full-cohort trait/protein objects saved by `5e`; falls back to `_numericMeta_3177_trait.RDS` / `_full_3177_protein_dft.RDS` if the full-cohort objects are absent)
- `./name_match_table.RDS`
- Per-assay no-imputation Model 1 STAN RDS files saved by `5e`.

Non-RData/RDS files or helper code to place in `input/`:
- `plot_functions_20250718.R`
- `scatterplot_label_20250718.csv`

Major outputs: per-assay scatter+ribbon/difference trajectory PDFs (`./scatter/99par_2on1_<assay>.pdf`), the t-statistic/p-value/up-down-direction matrices (`./scatter/_99_par_diff_all_peptide*.rds`/`.csv`) consumed by the prefix-`7` ExtractSigAssays scripts, a waterfall summary PDF, and a saved workspace image.

### `5f.model_3_STAN_withSex+Genotype+EYO_3way_int_terms_5xSDcleaned(RDS)-0727.R`
Fits Model 3 / Model 2.3 in the primary cohort, adding APOE genotype × sex × EYO spline interaction terms so that sex-specific APOE e4/e4 and e3/e3 longitudinal differences can be evaluated directly.

Required inputs:
- `_numericMeta_3177_trait.RDS`
- `_full_3177_protein_dft.RDS`

Major outputs: per-assay Model 3 (genotype × sex × EYO 3-way interaction) `stan_glm` fits.

### `5f3.model_3_STAN_OutputSummaryStats99CI-5xSDcleaned(RDS)-0727.R`
Performs post-hoc posterior testing and summary extraction for Model 3, computing EYO-grid posterior means, credible intervals, and two-tailed Bayesian posterior probabilities for genotype × sex × EYO interaction terms, with corrected sex coding.

Required inputs:
- `_numericMeta_3177_trait.RDS`
- `name_match_table.RDS`
- Per-assay Model 3 STAN RDS files in the configured Model 3 directory.

Major outputs: EYO-grid posterior mean/99% credible-interval summary tables and two-tailed posterior-probability tables per genotype × sex × EYO interaction term, consumed by `5f4`.

### `5f4.Plot_STAN_model_3_99CI_withSex+Genotype+EYO_3way_int_terms-5xSDcleaned(RDS)-0727.R`
Generates single-panel Model 3 plots showing genotype × sex interaction posterior ribbons across EYO and a colored significance streak indicating the dominant sex/genotype category at significant EYO bins.

Required inputs:
- `_numericMeta_3177_trait.RDS`
- `name_match_table.RDS`
- Per-assay Model 3 STAN RDS files and post-hoc output summaries from `5f3`.

Major outputs: per-assay single-panel genotype × sex interaction PDF plots with significant-EYO-bin streak annotation.

### `6_Endophenotypes_FETSplitLists_7289.R`
Tests whether the proteins identified as significant in each 5-year EYO sliding window (and in a simple before-/after-age-50 split) are enriched for proteins previously reported as associated with neuropathology or cognitive-function endophenotypes (Nat Aging cohort marker lists), using one-sided Fisher's Exact Test (FET) overlap enrichment via the lab's general-purpose `geneListFET()` wrapper. Strips non-protein `seq.####` control-assay entries from the sliding-window hit list before testing.

Required inputs:
- `./geneListFET.R` (source()'d FET/heatmap function)
- `67_5year_sliding_windowSigHits.csv` (per-window significant-assay lists; output of the prefix-`7` ExtractSigAssays scripts)
- `SOMAbkgr_7289.csv` (full SomaScan assay background gene list for FET)
- `SplitList_BEFORE_AFTER_50yo.csv` (before/after age-50 assay split category list)
- `ROSMAP_Path_Hits(NatAgingS8).csv` (Nat Aging pathology marker reference list)
- `META_CogFnS7_AmyloidosisNoE4.S2-pLT0.05.csv` (Nat Aging cognition marker reference list, p<0.05; loaded for reuse but not run in the current script)
- `META_CogFnS7_AmyloidosisNoE4.S2-FDR_LT0.05.csv` (Nat Aging cognition marker reference list, FDR<0.05)

Major outputs: `67_5year_sliding_windowSigHits_noNonProteinAssays.csv` (cleaned hit list), plus one heatmap PDF + companion FET-statistics CSV/XLSX per `geneListFET()` call, summarizing overlap enrichment between the sliding-window hit list and each Nat Aging marker reference list.

### `7a_ExtractSigAssays_2directions_ALL-up_down___7289.R`, `7b_ExtractSigAssays_2directions_UP_-7289.R`, `7c_ExtractSigAssays_2directions_DOWN_-7289.R`
Three companion scripts sharing one pipeline, run for the combined (ALL = up+down), UP-only, and DOWN-only significance-direction subsets respectively. For 67 overlapping 5-year EYO sliding windows, each script: (1) collects the assays significant at p≤0.005 anywhere in that window for its direction subset; (2) tabulates and plots assay counts per window with a WGCNA module-color overlay (this combined-direction figure is generated once, by `7a`); (3) runs GO/Reactome/WikiPathways/MSigDB hypergeometric over-/under-representation enrichment of each window's assay list against the full assay background, in parallel, using an embedded GOparallel-style enrichment engine; (4) summarizes per-window, per-ontology enrichment Z-scores as trajectory line plots and row-clustered heatmaps; and (5) subsets to a curated, manually defined set of 141 representative ontology terms spanning 18 categories and draws the final category-ordered publication heatmap. `7a` independently curates the 141-term/18-category selection and final row order from its own (ALL) results; `7b`/`7c` reuse that same curated term list/order and simply recompute UP-only/DOWN-only Z-scores for it, so the three publication heatmaps are term-for-term comparable.

Required inputs:
- `../simple.3177/scatter/_99_par_diff_all_peptide_p_value.rds` (Model 1 posterior p-value matrix; output of `5a2`/`5e2`)
- `../simple.3177/scatter/_99_par_diff_all_peptide_up_down_notation.csv` (Model 1 up/down direction-color matrix; output of `5a2`/`5e2`)
- `Human_GO_AllPathways_noPFOCR_with_GO_iea_June_01_2025_symbol.gmt` (GO/Reactome/WikiPathways/MSigDB gene-set database; large file, downloaded/updated periodically)
- `go.obo` (Gene Ontology OBO file, for GO-term redundancy removal; downloaded automatically if missing and running interactively)
- `ALL141terms_18categories.tsv` (curated 18-category/141-term ontology selection; read by all three scripts)
- `ALL141terms_18categories-FinalOrder.tsv` (final row order for the `7a` ALL publication heatmap)
- `ALL141terms_18categories-FinalOrder2.tsv` (final row order for the `7b`/`7c` UP/DOWN publication heatmaps)
- `ALL_heatmap18categoryROWorder.txt` (overrides the row order loaded from `7a`'s saved RData, used by `7b`/`7c`)
- `ALL.df_all-141ontologies.categories18.RData` (curated 141-term/18-category Z-score selection saved by `7a`; required input to `7b` and `7c`)

Major outputs:
- `3.plot.hitCounts_plusModulesHitOverlaid-07-27_5xSDoutliersNotInModels.pdf` (assay-count-per-window plot with WGCNA module overlay; written once, by `7a`)
- `ALL_assay-Zscore_trajectories.pdf` / `UP_…` / `DOWN_…` and their `_29_SELECTED_1pp.pdf` companions (per-ontology Z-score trajectory plots)
- `ALL_assaysONLY-Z_GO-heatmaps_rowsOrdered_top100.pdf` / `UP_…` / `DOWN_…` (row-clustered ontology heatmaps)
- `3.SEPAwindows(67)ALL(top100)_18categoriesSeparated(141assays)-Z_GO-heatmaps_rowsOrdered_greyBG.pdf` and UP/DOWN equivalents (final 18-category, 141-term publication heatmaps)
- `ALL(DOWN+UP)-141Heatmap-ordered_ontologies+Zscores.csv` / `UP-141Heatmap-…csv` / `DOWN-141Heatmap-…csv` (141×67 Z-score supplementary tables)
- `ALL(DOWN+UP)-Genes_Hit(141selectedOntologies_67_5yrWindows.csv` and UP/DOWN equivalents (per-term "genes hit" annotation)
- `ALL.df_all-141ontologies.categories18.RData` (from `7a`; required input to `7b`/`7c`)
- Saved workspace images (`Ssaved.image-…ALLhits…RData`, `saved.image-…UPhitsONLY…RData`, `saved.image-…DOWNhitsONLY…RData`) consumed by prefix `8` (CMAP) and prefix `9` (MAGMA)

### `8_CMAP_perm__21x3x100000_7289.R`
Single consolidated script (combining the former permutation-computation and plotting scripts) that tests whether the up/down direction of e4/e4 protein abundance trajectories within each of the 141 curated ontology terms (from prefix `7`) is connected to the semaglutide drug signature from a published reference study, more than expected by chance, using a permutation-based weighted connectivity score (Cscore/NCS, CMAP/L1000-style). Runs a 21-point minimum-effect-size (minES) sensitivity sweep across 3 EYO epochs (pre-clinical, peri-onset, post-onset) with 100,000 permutations per (ontology × epoch) cell, then draws the final category-ordered connectivity heatmap with significance-asterisk overlays.

Required inputs:
- `SemaglutideStudy_TableS6stats.csv` (published semaglutide drug-signature T-statistics; reference CMAP "drug" ranking)
- `0727_medians_all_assays(5xSD_outliers_excluded).rds` (e4/e4 effect-size matrix; output of the prefix-`5` Model 1 analysis)
- `ALL_heatmap18categoryROWorder.FinalOrder2-dataFrame.tsv` (curated 18-category/141-ontology order table; output of the prefix-`7` ExtractSigAssays scripts)
- `Human_GO_AllPathways_noPFOCR_with_GO_iea_June_01_2025_symbol.gmt` (GO/pathway gene-set database, used to build each ontology's background gene set)
- `Fig6A_ontology_order_110_18categories_final.tsv` and `Fig6A_ontology_order_141.tsv` (ontology/category order tables for the final heatmap figure)

Major outputs:
- `minES_sensitivity_sweep_summary.csv` / `minES_sensitivity_sweep_full_results.rds`
- `minES_sensitivity_sweep_traces.pdf` / `minES_sensitivity_sweep_traces_ggplot.pdf` / `minES_sensitivity_sweep_per_epoch_traces.pdf`
- `Ontology_NCS_Fig6A_style_sweep_results_list5_finalOrder.pdf` (final 18-category, 141-ontology, 3-epoch connectivity heatmap with significance stars)
- `CMAP_permutation_stats_S6_100000perm_minES.e4_0.10.xlsx` (permutation statistics workbook: Cscore, NCS, p-values, CIs, FDR, etc.)
- `saved.image-SemaS6_minESsweep.RData` / `saved.image-CMAP.perm.RData`

### `9.MAGMA_wrapper.R`
Defines parameters and calls `MAGMA.SPA()` for gene/protein set enrichment analysis using semaglutide or other MAGMA-style input files and module/gene-list vectors.

Required inputs:
- Input MAGMA CSV files in `MAGMAinputDir`, including the configured examples `Sema_S6_NominalP.csv` and `Sema_S6_Qvalue.csv`.
- Required global objects expected by `MAGMA.SPA()`, especially `moduleGeneList` or module/protein list equivalents.

Non-RData/RDS files or helper code to place in `input/`:
- `MAGMA.SPA_listVectorInput.R`
- `Sema_S6_NominalP.csv`
- `Sema_S6_Qvalue.csv`

Major outputs: MAGMA enrichment summary tables, permutation-adjusted enrichment statistics, and bar plots of enrichment across the supplied input gene lists.

## Suggested `input/` folder contents

The scripts intentionally retain original absolute paths for provenance. For repo portability, collect helper scripts and small tabular inputs into an `input/` subfolder and update paths accordingly. Large RData/RDS workspaces are not listed here unless they are small enough and appropriate for controlled sharing.

Recommended local files to collect:

```text
input/9b.4cohort.csv *
input/BH-SomaSIgnals_forR_andHDS_BH.map.csv *
input/ROSMAP-SomaSignalsForR_andHDS_RM.map.csv *
input/UDS-SomaSignalsForR_andHDS_UDS.map.csv *
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
input/go.obo
input/Human_GO_AllPathways_noPFOCR_with_GO_iea_June_01_2025_symbol.gmt
input/67_5year_sliding_windowSigHits.csv
input/SOMAbkgr_7289.csv
input/SplitList_BEFORE_AFTER_50yo.csv
input/ROSMAP_Path_Hits(NatAgingS8).csv
input/META_CogFnS7_AmyloidosisNoE4.S2-pLT0.05.csv
input/META_CogFnS7_AmyloidosisNoE4.S2-FDR_LT0.05.csv
input/geneListFET.R
input/ALL141terms_18categories.tsv
input/ALL141terms_18categories-FinalOrder.tsv
input/ALL141terms_18categories-FinalOrder2.tsv
input/ALL_heatmap18categoryROWorder.txt
input/ALL_heatmap18categoryROWorder.FinalOrder2-dataFrame.tsv
input/SemaglutideStudy_TableS2stats.csv
input/SemaglutideStudy_TableS6stats.csv
input/Fig6A_ontology_order_110_18categories.tsv
input/Fig6A_ontology_order_141.tsv
input/Fig6A_ontology_order_110_18categories_final.tsv
input/MAGMA.SPA_listVectorInput.R
input/Sema_S6_NominalP.csv
input/Sema_S6_Qvalue.csv
```
* Protected individual information in contents; files not publicly shared -- reach out for more information.

## Notes for use
- GNPC harmonized data set (HDS) v1.3 data is available on the ADDI Azure-implemented VM platform for authorized users.
- These scripts remain research pipeline scripts rather than package-style functions; paths and large workspace objects are intentionally left visible for provenance.
- Before rerunning on a different system, update `rootdir`/`setwd()` paths and confirm sourced helper files are available.
- The large GMT gene-set database file (`Human_GO_AllPathways_noPFOCR_with_GO_iea_June_01_2025_symbol.gmt`) and `go.obo` are used by both prefix `7` and prefix `8`; prefix `7`'s scripts will download `go.obo` automatically if it is missing and the session is interactive.

## Repository status (7289 assays)

The repository currently contains numerically prefixed R scripts from prefixes `1` through `9`, including the prefix-5 STAN modeling and plotting scripts (with `5e`/`5e2` patched to fix spline knot positions to those computed from the full, imputed-genotype-included 3,177-sample cohort). Other code updated to effect removal of non-human/non-protein assays in SomaScan 7k from cleaned data are found in .R files with prefix-`6` endophenotype/pathology FET enrichment script, the three-file prefix-`7` ExtractSigAssays (ALL/UP/DOWN) ontology trajectory pipeline, and the consolidated single-file prefix-`8` CMAP permutation pipeline.
