## =============================================================================
## 5e2. Plot Model 1 - APOE-genotype-known-only sensitivity cohort
##      (ribbons drawn with spline knots fixed to full 3177-sample positions)
## =============================================================================
##
## PURPOSE
##   Companion plotting script to 5e. For every assay's known-genotype-only
##   STAN model fit (saved by 5e as <CleanedAssayName>_stan_glm.rds), this
##   script: (1) evaluates the posterior at a fine EYO grid using the SAME
##   restricted cubic spline (RCS) knot positions used during model fitting
##   (i.e. positions estimated from the full 3,177-sample/imputed-genotype
##   cohort, not from the plotted, known-genotype-only points); (2) draws a
##   scatter+ribbon trajectory plot and a difference-curve plot per assay;
##   (3) collects t-statistics, two-sided posterior p-values, and a
##   cornflowerblue/indianred3 up/down significance flag per EYO bin for
##   every assay; and (4) renders a waterfall/heatmap summary of which
##   assays first reach significance at each EYO window.
##
## STEP-BY-STEP PIPELINE
##   1. Load the filtered (known-genotype-only) and full-cohort (knot
##      estimation) trait/protein objects saved by 5e
##      (assays+traits_forPlots.RData), with a fallback to the raw RDS files
##      for older RData snapshots that pre-date the full-cohort save.
##   2. Verify sample-ID alignment between trait and protein/outcome tables,
##      for both the filtered and full-cohort objects.
##   3. Pre-compute, once, a fine 0.5-year EYO grid spanning the full-cohort
##      EYO range and a single full-cohort RCS knot set used as the fallback
##      spline basis for every assay's prediction grid.
##   4. Load the OriginalName/CleanedName lookup table and module/y-axis
##      label table used to find each model's RDS file and label its plots.
##   5. For each assay, in process_one_peptide() (run in parallel):
##        a. Read the assay's saved STAN fit and skip if it failed to fit.
##        b. Recover the full-cohort spline knots used to fit that model.
##        c. Extract posterior alpha/beta draws and build e3/e3 vs e4/e4
##           contrast matrices across the EYO grid using those knots.
##        d. Compute posterior median and 99% credible-interval trajectories
##           for each genotype group and for their difference at every EYO
##           grid point, plus a two-sided posterior p-value per EYO bin.
##        e. Draw a scatter plot of the known-genotype-only data points with
##           overlaid posterior ribbons, and a separate genotype-difference
##           plot with significance-interval shading; combine both into one
##           PDF page per assay.
##        f. Return the per-EYO t-statistics, p-values, and up/down flags.
##   6. Recombine all per-assay results into 143(EYO)x7345(assay) matrices/
##      data frames of t-statistics, p-values, and up/down direction flags,
##      and save them (RDS/CSV) for downstream sliding-window summarization.
##   7. Compute, per assay, the EYO midpoint of the first contiguous run of
##      significant EYO bins (or "constitutive" if significant everywhere),
##      as a QC/descriptive summary of the no-imputation sensitivity result.
##   8. Build a waterfall-style heatmap of -log2(p) (signed by direction) for
##      a representative subset of assays/terms, ordered by earliest onset of
##      significance, and save it as a PDF.
##   9. Save the full workspace image for provenance/downstream reuse.
##
## REQUIRED INPUTS
##   - ./assays+traits_forPlots.RData   (filtered + full-cohort trait/protein
##                                       objects saved by 5e; falls back to
##                                       _numericMeta_3177_trait.RDS and
##                                       _full_3177_protein_dft.RDS if the
##                                       full-cohort objects are absent)
##   - ./name_match_table.RDS           (OriginalName/CleanedName lookup)
##   - ./scatterplot_label_20250718.csv (per-assay y-axis/title labels)
##   - ./plot_functions_20250718.R      (source()'d scatter_plot()/diff_plot()
##                                       ggplot helper functions)
##   - Per-assay <CleanedAssayName>_stan_glm.rds STAN model fits saved by 5e,
##     in the same working directory
##
## MAJOR OUTPUTS
##   - Per-assay 2-panel scatter+ribbon / difference PDF
##     (./scatter/99par_2on1_<assay>.pdf)
##   - ./scatter/_99_par_diff_all_peptide.rds / .csv (t-statistics matrix)
##   - ./scatter/_99_par_diff_all_peptide_p_value.rds / .csv (p-value matrix;
##     consumed by the prefix-7 ExtractSigAssays scripts)
##   - ./scatter/_99_par_diff_all_peptide_up_down_notation.csv (direction flags)
##   - ./scatter/_Waterfall-NoImputedGenotypes_FixedSplineKnotsToOriginalPositions_withImputedSamples.pdf
##   - ./scatter/_saved.image.RData (full workspace snapshot)
## =============================================================================

#' All "BL" stands for "Last measure" (the longitudinal endpoint variable naming convention)

#rm(list=ls()) 

library(tidyverse)
library(rstanarm)
library(Hmisc)
library(openxlsx)
library(rstan)
library(gridExtra)
library(ggpubr)

# Comment and un-comment on HPC
#setwd("./ShijiaBian/PlasmaProteomic/Result/20250718/")
#setwd("z:/ShijiaBian/PlasmaProteomic/Result/20250727/simple.3177")
#setwd("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/")
setwd("/home/labshare/model1.noImpute/")


# ---- STEP 1. Load the filtered (known-genotype-only) and full-cohort
# (knot-estimation) trait/protein objects saved by 5e. (The original,
# pre-2026-05-27 approach of reading the raw RDS files directly and
# re-deriving EYO/filtering here is superseded by 5e's saved RData and is
# removed; a fallback for older RData snapshots is kept just below.) ----
load("assays+traits_forPlots.RData")
BL_traits_pep$ApoE_Indicator<-BL_traits$ApoE_Indicator

# CHANGED 2026-05-27: The updated 5e script saves these full-3177 objects.
# Keep a fallback for older assays+traits_forPlots.RData files by re-reading
# the source RDS files if the full-data objects are absent.  The full-data
# objects are used only for spline knot/range calculation, not for plotted
# points or model fitting.
if (!exists("BL_traits_full3177_for_knots") || !exists("BL_traits_pep_full3177_for_knots")) {
  BL_traits_full3177_for_knots <- readRDS("_numericMeta_3177_trait.RDS")
  BL_traits_pep_full3177_for_knots <- readRDS("_full_3177_protein_dft.RDS")
}

if (!"EY0" %in% names(BL_traits_full3177_for_knots)) {
  BL_traits_full3177_for_knots$EY0 <- (65.6 - BL_traits_full3177_for_knots$age_at_visit) * (-1)
}
BL_traits_full3177_for_knots$EYO <- BL_traits_full3177_for_knots$EY0
BL_traits_pep_full3177_for_knots$EYO <- BL_traits_full3177_for_knots$EYO

# ---- STEP 2. Verify sample-ID alignment between trait and protein tables ----
stopifnot(all(BL_traits_pep$sample_id == BL_traits$sample_id))
stopifnot(all(BL_traits_pep_full3177_for_knots$sample_id == BL_traits_full3177_for_knots$sample_id))




min(BL_traits_pep$EYO) # -45.6     -31.53973
max(BL_traits_pep$EYO) # 24.4       23.72 previously

# CHANGED 2026-05-27: Precompute the plotting EYO grid and a full-3177
# fallback knot set once, before parallel plotting.  New STAN model RDS files
# produced by the updated 5e script carry protein-specific full-3177 knots;
# those are used preferentially below so plotted ribbons match each model.
eyo_step_full3177 <- seq(floor(min(BL_traits_pep_full3177_for_knots$EYO, na.rm = TRUE)),
                         ceiling(max(BL_traits_pep_full3177_for_knots$EYO, na.rm = TRUE)),
                         by = 0.5)
full3177_splinefit_for_plots <- Hmisc::rcspline.eval(BL_traits_pep_full3177_for_knots$EYO,
                                                     nk = 3, norm = 2, pc = FALSE,
                                                     inclx = TRUE)
full3177_spline_knots_for_plots <- attr(full3177_splinefit_for_plots, "knots")

get_model_spline_knots <- function(stan_BL, fallback_knots) {
## ideally, would use knots for specific data modeled (minus outlier samples), but we did not do this originally when imputed samples were included, so following 3 lines are commented, to force the generic fit of knot positions for all models (which uses all samples including those with imputed genotypes)
#  model_knots <- attr(stan_BL, "full3177_spline_knots")
#  if (!is.null(model_knots)) return(model_knots)
#  if (!is.null(stan_BL$full3177_spline_knots)) return(stan_BL$full3177_spline_knots)
  fallback_knots
}

# ---- STEP 3. Load the OriginalName/CleanedName lookup and per-assay
# y-axis/title label table used to find each model's saved RDS file and
# label its plots. ----
EYO_cut = length(eyo_step_full3177) # full 3177-sample EYO grid

name_match_table <- readRDS("./name_match_table.RDS")

#*** we do not have MEs now
##name_module_label <- read.csv("../scatterplot_label_20250718.csv", header = T)
name_module_label <- read.csv("./scatterplot_label_20250718.csv", header = T)
name_module_label.names<-colnames(name_module_label)
#name_module_label.add<-cbind(name_match_table$CleanedName,name_match_table$OriginalName,name_match_table$OriginalName)
#colnames(name_module_label.add)<-colnames(name_module_label)
#name_module_label<-as.data.frame(rbind(name_module_label,name_module_label.add))
name_module_label<-cbind(name_match_table$CleanedName,name_match_table$OriginalName,name_match_table$OriginalName)
colnames(name_module_label)<-name_module_label.names
name_module_label<-as.data.frame(name_module_label)

################## -------- Fit the STAN model -------- ###########

pep_names <- names(BL_traits_pep)[c(2:7346)]  #12 modules at end, no MMSE and cdr

# ---- STEP 4. Source the ggplot2 scatter_plot()/diff_plot() helper functions
# used by process_one_peptide() below to render each assay's figure. ----
#source("z:/ShijiaBian/PlasmaProteomic/Code/CommonFunctions/plot_functions_20250718.R")   #202311235.R  changed from _Der_F to ApoE_Indicator
source("./plot_functions_20250718.R")   #202311235.R  changed from _Der_F to ApoE_Indicator


# ---- STEP 5 (per assay, run in parallel below). process_one_peptide():
# reads one assay's saved STAN fit, evaluates the e3/e3 and e4/e4 posterior
# trajectories (and their difference) across the EYO grid using the
# full-cohort spline knots, draws and saves a combined scatter+difference
# PDF, and returns that assay's per-EYO t-statistic/p-value/direction flag. ----
process_one_peptide <- function(pep_name) {
	clean_name = name_match_table$CleanedName[which(name_match_table$OriginalName == pep_name)]
	fil_stan <- file.path("./", paste(paste(clean_name, "_stan_glm",  ".rds", sep = "")))

	stan_BL <- readRDS(fil_stan)
  
	if (length(stan_BL) == 1) {
	  #next  # only valid for/while loop
	  return(list(ok = FALSE))
	}

	# CHANGED 2026-05-27: Use the full-3177 spline knots stored with each
	# updated model RDS.  If plotting an older model that lacks those attributes,
	# fall back to the full-3177 knot set computed once above.
	model_spline_knots <- get_model_spline_knots(stan_BL, full3177_spline_knots_for_plots)

	###### ----------- Plot 1: Two separate plots for carrier and non-carrier ------------- ############
	parameter_estimates <- rstan::extract(stan_BL$stanfit) # Extract the model estimates from STAN
	
	# biomarker ~ EYO_Spline_Linear + EYO_Spline_Cubic + Group + EYO_Spline_Linear * Group + EYO_Spline_Cubic 
	# Add intercept to beta weights, 4000 * 6
	factorweights <- cbind(parameter_estimates$alpha, parameter_estimates$beta) # Add intercept to beta weights 
	
	#Initialize contrast matrices
	# CHANGED 2026-05-27: Use the EYO grid derived once from the full 3177
	# samples; plotted points below remain restricted to BL_traits_pep
	# (unimputed-genotype samples).
	eyo_step = eyo_step_full3177

	#Generate output matrices, nrow = 4,000 (stan iterations); ncol = 125, 125 is the number of eyo
	spline_noncarriers = matrix(0, nrow = nrow(factorweights), ncol = length(eyo_step)) #generate blank matrices for outputs for noncarriers
	spline_carriers = matrix(0, nrow = nrow(factorweights), ncol = length(eyo_step)) #generate blank matrices for outputs for carriers
	spline_diff = matrix(0, nrow = nrow(factorweights), ncol = length(eyo_step)) #generate blank matrices for output of differences

	contrasts_non_carriers = matrix(0, nrow = length(eyo_step), ncol = dim(factorweights)[2])
	contrasts_carriers = matrix(0, nrow = length(eyo_step), ncol = dim(factorweights)[2]) #Making a blank matrix for contrasts
	# CHANGED 2026-05-27: Do not choose knots from the filtered plotting data.
	# Evaluate the prediction basis with the full-3177 knots associated with
	# the model, preserving consistency between the Stan coefficients and ribbons.
	spline_knots_for_prediction <- model_spline_knots

	for (j in 1:length(eyo_step)) {
		tempfit = Hmisc::rcspline.eval(eyo_step[j], knots = spline_knots_for_prediction, norm=2, pc=FALSE, inclx=TRUE)       #put all of the values from our plotting interval into rcscpline to get the cubic term values
		#Contrast 
#		contrasts = matrix(
#			c(1, 1, tempfit[1,1],tempfit[1,1],tempfit[1,2],tempfit[1,2], 0,1, 0, tempfit[1,1], 0,tempfit[1,2]),
#			nrow=2, 
#			ncol=6)
#		contrasts_non_carriers[j,] = contrasts[1,]
#		contrasts_carriers[j,] = contrasts[2,] #Set up matrix of contrasts for each EYO point separately for carriers and noncarriers to be used below

##              ORIGINAL MODEL
#		biomarker ~ EYO_L + EYO_C + Group + EYO_L:Group + EYO_C:Group
		# tempfit gives the spline bases at the current EYO value
		xL <- tempfit[1, 1]
		xC <- tempfit[1, 2]

		# prediction for NON-carriers (Group = 0)
		contrasts_non_carriers[j, ] <- c(1,  xL,  xC,  0,      0,        0)

		# prediction for carriers (Group = 1)
		contrasts_carriers    [j, ] <- c(1,  xL,  xC,  1,      xL,       xC)

	}

	#Begin the loop to go through every single iteration of Stan 
	for (i in 1:nrow(factorweights)) {
		weights_temp = factorweights[i,] #for every iteration select the weights

		# for every EYO in our output what would the value be using this iterations model fits
		for (j in 1:length(eyo_step)) {
			#select the appropriate contrast from the matrix made above
			contrasts = rbind(contrasts_non_carriers[j,], contrasts_carriers[j,]) 

			# Multiple beta weights by the contrast to get one value per group for this EYO and this iteration
			weights_point = contrasts %*% weights_temp 

			# Write the points out into a matrix for each group
			spline_noncarriers[i,j] = weights_point[1,]
			spline_carriers[i,j] = weights_point[2,] 

			#difference is equal to carriers - noncarriers. This distribution will be done to determine the actual deviation point
			spline_diff[i, j] = weights_point[2,] - weights_point[1,] 
		}
	} 

	## Extract median and 0.005 and 0.995 values
	#generate blank matrices we will use for lines
	noncarrier_lines = matrix(0, ncol=3, nrow=length(eyo_step))
	carrier_lines = matrix(0, ncol=3, nrow=length(eyo_step))
	diff_lines = matrix(0, ncol=3, nrow=length(eyo_step))
	noncarrier_carrier_t_stats = matrix(0, ncol=1, nrow=length(eyo_step))
	noncarrier_carrier_p_value = matrix(0, ncol=1, nrow=length(eyo_step))
	for (i in 1:ncol(spline_noncarriers)) {
		# For every EYO point (i)
		temp_non = quantile(spline_noncarriers[,i], probs = c(0.005, .5, .995))
		temp_non = unname(temp_non) #get the 0.05, median, and 99.5 percentile values
		temp_car = quantile(spline_carriers[,i], probs = c(0.005, .5, .995))
		temp_car = unname(temp_car) #get the 0.05, median, and 99.5 percentile values
		temp_diff = quantile(spline_diff[,i], probs = c(0.005, .5, .995))
		temp_diff = unname(temp_diff) #get the 0.05, median, and 99.5 percentile values

		noncarrier_lines[i,1] = temp_non[1]
		noncarrier_lines[i,2] = temp_non[2]
		noncarrier_lines[i,3] = temp_non[3] #Write out the 99% Credible intervals and median for each EYO (i) 
		carrier_lines[i,1] = temp_car[1]
		carrier_lines[i,2] = temp_car[2]
		carrier_lines[i,3] = temp_car[3] # Carriers
		diff_lines[i,1] = temp_diff[1]
		diff_lines[i,2] = temp_diff[2]
		diff_lines[i,3] = temp_diff[3] 

		temp_test_stats = t.test(spline_diff[,i])
		noncarrier_carrier_t_stats[i] = temp_test_stats$statistic[["t"]]
		noncarrier_carrier_p_value[i] = min(sum(spline_diff[,i] < 0)/length(spline_diff[,i]), 1-(sum(spline_diff[,i] < 0)/length(spline_diff[,i])))
	}

#	noncarrier_carrier_t_stats_all_pep <- cbind(noncarrier_carrier_t_stats_all_pep, noncarrier_carrier_t_stats)
#	colnames(noncarrier_carrier_t_stats_all_pep)[ncol(noncarrier_carrier_t_stats_all_pep)] = pep_name

#	noncarrier_carrier_p_value_all_pep <- cbind(noncarrier_carrier_p_value_all_pep, noncarrier_carrier_p_value)
#	colnames(noncarrier_carrier_p_value_all_pep)[ncol(noncarrier_carrier_p_value_all_pep)] = pep_name

	# Add EYO term and column names to matrices for plotting
	flipped_eyo = t(eyo_step)
	flipped_eyo = t(flipped_eyo)
	# Flip everything and make it a matrix
	combined_lines = cbind(flipped_eyo, noncarrier_lines, carrier_lines) 
	#combine carriers and noncarriers
	colnames(combined_lines) <- c("eyo","lower_non", "median_non", "upper_non", "lower_car", "median_car", "upper_car") 
	diff_lines = cbind(flipped_eyo, diff_lines)
	colnames(diff_lines) <- c("eyo","lower", "median", "upper") #separate matrix for difference, honest this could be combined with the previous one

	combined_lines = as.data.frame(combined_lines) #Make it a data frame or none of the plotting works
	diff_lines = as.data.frame(diff_lines) #Make it a data frame or none of the plotting works

	# manually set difference of medians
	diff_lines$median <- combined_lines$median_car - combined_lines$median_non

	### --------- Plot Figure 1 Starts
	### Plot Images
	temp_data_plot_before <- BL_traits_pep %>%
		select(sample_id, EYO, eval(pep_name), ApoE_Indicator) %>%  #Group_Der_F
		ungroup() %>%
		mutate(Pep = BL_traits_pep[[pep_name]] ) %>%
		dplyr::filter(complete.cases(.)) 

	before_removal_outlier <- nrow(temp_data_plot_before)

	temp_data_plot <- temp_data_plot_before	%>%
		mutate(zRT = scale(Pep)[,1]) %>%
		filter(between(zRT, -3, +3))
	post_removal_outlier <- nrow(temp_data_plot)
	number_removed <- before_removal_outlier - post_removal_outlier # record number of removed outliers

	CI_Percent <- "99%"
	
	yaxis_label_plot = name_module_label$y.Axis.Label[which(name_module_label$Raw.File.Name == clean_name)]
	title_label_plot = name_module_label$Title.Label[which(name_module_label$Raw.File.Name == clean_name)]

	ScatterPlot = scatter_plot(x = EYO, 
							   y = Pep, 
							   df = temp_data_plot, 
							   protein_name = data.frame(yaxis_label = yaxis_label_plot,
							   							 title_label = title_label_plot), 
							   number_removed = number_removed,
							   CI_Percent = CI_Percent,
							   combined_lines = combined_lines)

#        ggsave(filename = paste0("z:/ShijiaBian/PlasmaProteomic/Result/20250727/simple.3177/scatter/",clean_name,".pdf"), plot = ScatterPlot, device = cairo_pdf,
#               width = 12, height = 12, units = "in", dpi = 300)

	### --------- Plot Figure 1 Ends

	# ### --------- Plot Figure 2 Starts
	min_value = if(pep_name %in% colors()) { min(-0.001, min(diff_lines$lower))*1.3 } else { floor(min(diff_lines$lower)) }  #min range -0.001 to + 0.001 if a color (ME), else -1 to +1 (or greater range)
	max_value = if(pep_name %in% colors()) { max(0.001, max(diff_lines$upper))*1.3 } else { ceiling(max(diff_lines$upper)) }

	# Construct the diff data frame to better plot the significance difference
	diff_lines$same_sign = FALSE
	diff_lines$same_sign[sign(diff_lines$lower) == sign(diff_lines$upper)] = TRUE
	Carr_Non_EYO_diff_df = diff_lines[diff_lines$same_sign == TRUE, ]
	if (nrow(Carr_Non_EYO_diff_df) > 0) {
		Carr_Non_EYO_diff_df$Height = max_value
		Carr_Non_EYO_diff_df$Text_Height = (max_value - min_value) * 0.05
		First_Annotate <- Carr_Non_EYO_diff_df$eyo[1]
		Last_Annotate <- Carr_Non_EYO_diff_df$eyo[length(Carr_Non_EYO_diff_df$eyo)]

		# single color, regardless of # of intervals
		#Annotate_color <- "cornflowerblue"
		#if (sum(Carr_Non_EYO_diff_df$median > 0)/length(Carr_Non_EYO_diff_df$median) > 0.5) { # majortiy of the median > 0
		#	Annotate_color <- "indianred3"
		#}
                # A color for each interval
                Annotate_color <- na.omit(ifelse(Carr_Non_EYO_diff_df$upper < 0, "cornflowerblue", ifelse(Carr_Non_EYO_diff_df$lower > 0, "indianred3", NA)))

		right_limit = 24
		Last_Annotate_Number = Last_Annotate + 0.5
		if (Last_Annotate_Number > right_limit){
			right_limit = 26
		}
	} else {
		First_Annotate = Inf
		Last_Annotate = Inf
		Annotate_color = ""
		right_limit = 24
	}
	
	DiffPlot = diff_plot(x = eyo, 
							y = median, 
							df = diff_lines, 
							bar_df = Carr_Non_EYO_diff_df,
							protein_name = data.frame(yaxis_label = yaxis_label_plot,
													  title_label = title_label_plot), 
							min_value =  min_value, # -31121.83
							max_value = max_value, # 14325
							First_Annotate = First_Annotate,
							Last_Annotate = Last_Annotate,
							Annotate_color = Annotate_color,
							right_limit = right_limit,
							Last_Annotate_Number = Last_Annotate_Number,
							CI_Percent = CI_Percent)

#        ggsave(filename = paste0("z:/ShijiaBian/PlasmaProteomic/Result/20250727/simple.3177/scatter/99par_diff_",clean_name,".pdf"), plot = DiffPlot, device = cairo_pdf,
#               width = 12, height = 12, units = "in", dpi = 300)


	
	up_down_notation<-rep(NA, EYO_cut) #[count, 'Pep'] = pep_name 
	names(up_down_notation)<-as.character(eyo_step)  #as.character(Carr_Non_EYO_diff_df$eyo)
	if (dim(Carr_Non_EYO_diff_df)[1] > 0) {
		for (r in 1:dim(Carr_Non_EYO_diff_df)[1]) {
			if (Carr_Non_EYO_diff_df[r, "lower"] < 0 & Carr_Non_EYO_diff_df[r, "upper"] < 0) {
				up_down_notation[which(names(up_down_notation) == as.character(Carr_Non_EYO_diff_df[r, "eyo"]))] = "cornflowerblue"
			}
			if (Carr_Non_EYO_diff_df[r, "lower"] > 0 & Carr_Non_EYO_diff_df[r, "upper"] > 0) {
				up_down_notation[which(names(up_down_notation) == as.character(Carr_Non_EYO_diff_df[r, "eyo"]))] = "indianred3"
			}
				
		}
	}
	# ### --------- Plot Figure 2 Ends 

	ScatterPlot = scatter_plot(x = EYO, 
							   y = Pep, 
							   df = temp_data_plot, 
							   protein_name = data.frame(yaxis_label = yaxis_label_plot,
							   							 title_label = title_label_plot), 
							   number_removed = number_removed,
							   CI_Percent = CI_Percent,
							   combined_lines = combined_lines)
	DiffPlot = diff_plot(x = eyo, 
							y = median, 
							df = diff_lines, 
							bar_df = Carr_Non_EYO_diff_df,
							protein_name = data.frame(yaxis_label = yaxis_label_plot,
							   						  title_label = title_label_plot), 
							min_value =  min_value, # -31121.83
							max_value = max_value, # 14325
							First_Annotate = First_Annotate,
							Last_Annotate = Last_Annotate,
							Annotate_color = Annotate_color,
							right_limit = right_limit,
							Last_Annotate_Number = Last_Annotate_Number,
							CI_Percent = CI_Percent)

	gt = arrangeGrob(ScatterPlot, DiffPlot, ncol = 2)
	on_same_plot <- ggplotify::as.ggplot(gt, draw=FALSE)

#        ggsave(filename = paste0("z:/ShijiaBian/PlasmaProteomic/Result/20250727/simple.3177/scatter/99par_2on1_",gsub("(\\|)|(\\.)|(\\;)", "_",pep_name),".pdf"), plot = on_same_plot, device = cairo_pdf,
#        ggsave(filename = paste0("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/noImpute.simple.rerun/scatter/99par_2on1_",gsub("(\\|)|(\\.)|(\\;)", "_",pep_name),".pdf"), plot = on_same_plot, device = cairo_pdf,
        ggsave(filename = paste0("/home/labshare/model1.noImpute/scatter/99par_2on1_",gsub("(\\|)|(\\.)|(\\;)", "_",pep_name),".pdf"), plot = on_same_plot, device = cairo_pdf,               width = 22.5, height = 10, units = "in", dpi = 300)

  # return a named list with the objects you need afterwards
  list(
    t_stat   = noncarrier_carrier_t_stats,
    p_value  = noncarrier_carrier_p_value,
    up_down  = up_down_notation, #[count, ],
    ok       = TRUE                 # handy flag for foreach result binding
  )
}


library(doParallel)
ncore <- max(1, parallel::detectCores() - 1)
cl     <- makeCluster(ncore)
registerDoParallel(cl)


library(foreach)

# list of packages the workers must load
worker_pkgs <- c("tidyverse","rstanarm","Hmisc","openxlsx","rstan",
                 "gridExtra","ggpubr","ggplot2","grDevices","ggplotify")

# ---- STEP 6. Run process_one_peptide() for every assay in parallel,
# producing one combined scatter+difference PDF per assay plus per-assay
# t-statistic / p-value / direction vectors collected in `results`. ----
results <- foreach(pep_name = pep_names,
                   .packages = worker_pkgs,
                   .export   = c("process_one_peptide",      # the function
                                 "diff_plot","scatter_plot", # custom plotters
                                 "get_model_spline_knots",
                                 "BL_traits_pep","name_match_table",
                                 "name_module_label",        # big objects
                                 "EYO_cut", "eyo_step_full3177",
                                 "full3177_spline_knots_for_plots",
                                 "cairo_pdf"),     # scalars/funs
                   .errorhandling = "pass")  %dopar%  {
  process_one_peptide(pep_name)
}

stopImplicitCluster()   # or stopCluster(cl) if you built one explicitly


# ---- STEP 7. Recombine per-assay parallel results into full
# 143(EYO)x7345(assay) t-statistic, p-value, and up/down direction matrices. ----
# Filter out failed iterations
ok_idx <- vapply(results, function(x) is.list(x) && !inherits(x, "error") && x$ok, TRUE)

noncarrier_carrier_t_stats_all_pep    <- do.call(cbind,
                                                 lapply(results[ok_idx], `[[`, "t_stat"))
noncarrier_carrier_p_value_all_pep    <- do.call(cbind,
                                                 lapply(results[ok_idx], `[[`, "p_value"))
up_down_notation                      <- do.call(rbind,
                                                 lapply(results[ok_idx], `[[`, "up_down")) #[ , -1]

## Set dimnames of final data
colnames(noncarrier_carrier_t_stats_all_pep)<-colnames(noncarrier_carrier_p_value_all_pep)<-rownames(up_down_notation)<-pep_names[ok_idx]

# CHANGED 2026-05-27: Use the same full-3177 EYO grid for final matrices.
eyo_step = eyo_step_full3177
length(eyo_step)==ncol(up_down_notation)  # TRUE
rownames(noncarrier_carrier_t_stats_all_pep)<-rownames(noncarrier_carrier_p_value_all_pep)<-colnames(up_down_notation)<-as.character(eyo_step)


## ---- STEP 8. Summarize how many EYO bins/assays reach significance, and
## compute the EYO midpoint of each assay's first significant interval
## (firstSigGenoMidpoints) as a descriptive QC summary of this
## known-genotype-only / full-cohort-knots sensitivity result. ----
## Get characteristics of waterfall
apply(up_down_notation,2,function(x) table(x))

sum(apply(up_down_notation,1,function(x) length(which(!is.na(x))))>0)
#1013 using full data (with imputed samples) spline knot positions # previously, with spline knots repositioned: 1119 (on R v4.5.1, was 1124 before rerun v4.4.1) ... 881.. 886



# ============================================================
# firstSigGenoMidpoints
# For every peptide in BL_traits_pep columns 2:7346 (in that
# order), return the midpoint EYO value of the FIRST (leftmost)
# contiguous interval where the genotype Bayesian p-value
# is < 0.0051.
#   "constitutive" -> the ENTIRE EYO range is significant
#   NA             -> no EYO bin reaches significance
# The result is a character vector of length 7345 whose names
# are the peptide/column names.
# ============================================================

p_thresh        <- 0.0051
eyo_step_vals   <- as.numeric(rownames(noncarrier_carrier_p_value_all_pep))
all_pep_ordered <- names(BL_traits_pep)[2:7346]   # ordered reference

firstSigGenoMidpoints <- vapply(all_pep_ordered, function(pn) {

  # Peptide failed / was not modelled
  if (!pn %in% colnames(noncarrier_carrier_p_value_all_pep))
    return(NA_character_)

  pvec <- noncarrier_carrier_p_value_all_pep[, pn]
  sig  <- pvec < p_thresh

  if (!any(sig)) return(NA_character_)   # no bin significant ? NA
  if (all(sig))  return("constitutive") # entire range significant

  # Locate the first contiguous run of significant bins
  rl        <- rle(sig)
  ends      <- cumsum(rl$lengths)
  starts    <- ends - rl$lengths + 1L
  first_t   <- which(rl$values)[1L]     # index of first TRUE run
  run_start <- starts[first_t]
  run_end   <- ends[first_t]

  # Midpoint EYO of that run (may fall between 0.5-unit grid points)
  mid_eyo <- (eyo_step_vals[run_start] + eyo_step_vals[run_end]) / 2
  as.character(mid_eyo)

}, character(1))

# Verify length and order
stopifnot(length(firstSigGenoMidpoints) == 7345)
stopifnot(names(firstSigGenoMidpoints)  == all_pep_ordered)

# Counts
message(sprintf(
  "firstSigGenoMidpoints: %d constitutive | %d with first-sig EYO | %d NA",
  sum(firstSigGenoMidpoints == "constitutive", na.rm = TRUE),
  sum(!is.na(firstSigGenoMidpoints) & firstSigGenoMidpoints != "constitutive"),
  sum(is.na(firstSigGenoMidpoints))
))
#now: using imputed samples for spline knots (but not models)
#firstSigGenoMidpoints: 17 constitutive | 996 with first-sig EYO | 6332 NA

#previously: spline knots not using imputed samples
#firstSigGenoMidpoints: 18 constitutive | 1102 with first-sig EYO | 6225 NA


all_pep_ordered[which(firstSigGenoMidpoints=="constitutive")]
#now:
# [1] "NEFL|P07196"                         "LRRN1|Q6UXK5^SL025922@seq.11293.14"  "TBCA|O75347"                         "ARL2|P36404"                         "CTF1|Q16619^SL002783@seq.13732.79"  
# [6] "PHGDH|O43175^SL018791@seq.15548.35"  "CDA|P32320"                          "FAM50A|Q14320"                       "BCDIN3D|Q7Z5W3"                      "ST8SIA1|Q92185^SL022499@seq.21508.7"
#[11] "PPM1G|O15355"                        "FOXO1|Q12778"                        "SPC25|Q9HBM1"                        "ZW10|O43264"                         "MENT|Q9BUN1"                        
#[16] "S100A13|Q99584"                      "DCUN1D5|Q9BTE7" 

#previously
# [1] "NEFL|P07196"                         "LRRN1|Q6UXK5^SL025922@seq.11293.14"  "TBCA|O75347"                        
# [4] "ARL2|P36404"                         "CTF1|Q16619^SL002783@seq.13732.79"   "VPS29|Q9UBQ0"                       
# [7] "PHGDH|O43175^SL018791@seq.15548.35"  "CDA|P32320"                          "FAM50A|Q14320"                      
#[10] "BCDIN3D|Q7Z5W3"                      "ST8SIA1|Q92185^SL022499@seq.21508.7" "PPM1G|O15355"                       
#[13] "FOXO1|Q12778"                        "SPC25|Q9HBM1"                        "ZW10|O43264"                        
#[16] "MENT|Q9BUN1"                         "S100A13|Q99584"                      "DCUN1D5|Q9BTE7"


# ---- STEP 9. Write the final t-statistic, p-value, and direction matrices
# used downstream by the prefix-7 ExtractSigAssays sliding-window scripts. ----
# ######### ---- Write Final Outputs for Waterfall(next) step ---- ########
noncarrier_carrier_t_stats_all_pep_final <- noncarrier_carrier_t_stats_all_pep #[, -1]
#fil_pep <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files", paste(paste("99_par_diff_all_peptide", ".rds", sep = "")))
fil_pep<-"./scatter/_99_par_diff_all_peptide.rds"
saveRDS(noncarrier_carrier_t_stats_all_pep_final, fil_pep)
#fil_csv_pep <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files",paste(paste("99_par_diff_all_peptide", ".csv", sep = "")))
fil_csv_pep <- "./scatter/_99_par_diff_all_peptide.csv"
write.csv(noncarrier_carrier_t_stats_all_pep_final, fil_csv_pep, row.names = TRUE)

noncarrier_carrier_p_value_all_pep_final <- noncarrier_carrier_p_value_all_pep #[, -1]
#fil_pep <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files",paste(paste("99_par_diff_all_peptide_p_value", ".rds", sep = "")))
fil_pep<-"./scatter/_99_par_diff_all_peptide_p_value.rds"
saveRDS(noncarrier_carrier_p_value_all_pep_final, fil_pep)
#fil_csv_pep <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files",paste(paste("99_par_diff_all_peptide_p_value", ".csv", sep = "")))
fil_csv_pep <- "./scatter/_99_par_diff_all_peptide_p_value.csv"
write.csv(noncarrier_carrier_p_value_all_pep_final, fil_csv_pep, row.names = TRUE)

#fil_csv_up_down_notation <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files",paste(paste("99_par_diff_all_peptide_up_down_notation", ".csv", sep = "")))
fil_csv_up_down_notation<-"./scatter/_99_par_diff_all_peptide_up_down_notation.csv"
write.csv(up_down_notation, fil_csv_up_down_notation, row.names = TRUE)

#sink("Result/20240122/Log_File/Original_STAN_Model_Scatter_Plot.txt", append = TRUE)
#print(Sys.time())
#sink()


# ---- STEP 10. Build a waterfall-style -log2(p) heatmap (signed by
# direction) of representative assays, ordered by earliest onset of
## significance (firstSigGenoMidpoints), and save it as a PDF. ----
## Waterfall plot - in R

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)

##--- 0) Inputs expected ---------------------------------------------
## noncarrier_carrier_p_value_all_pep_final : 143 x 599 data.frame (rows = EYO -46..25 as rownames, cols = ontology terms)
## up_down_notation                         : 599 x 143 data.frame (rows = ontology terms as rownames, cols = EYO -46..25)
## ontology_clusts                          : data.frame(term, cluster)  -- only for GSVA output

##--- 1) Long-format p-values (143*599 rows) -------------------------
pvals_long <- as.data.frame(noncarrier_carrier_p_value_all_pep_final) %>%
  mutate(EYO = as.numeric(rownames(.))) %>%
  relocate(EYO) %>%
  pivot_longer(-EYO, names_to = "term", values_to = "p")

## Cap extreme p's: p < 1/4000 -> 1/8000
pvals_long <- pvals_long %>%
  mutate(p_cap = if_else(p < 1/4000, 1/8000, p))

##--- 2) Long-format "direction" (cornflowerblue / indianred3 / NA) --
dir_long <- as.data.frame(up_down_notation) %>%
  mutate(term = rownames(.)) %>%
  relocate(term) %>%
  pivot_longer(-term, names_to = "EYO", values_to = "direction") %>%
  mutate(EYO = as.numeric(EYO))

##--- 3) Merge with clusters & compute selection features ------------
dat <- pvals_long %>%
  left_join(dir_long, by = c("term","EYO")) # %>%
#  left_join(ontology_clusts, by = "term")

# Keep only positions that are "eligible significant" (p < 0.0051 AND have a direction)
sig_dat <- dat %>%
  filter(!is.na(direction), p < 0.0051)

# Per term, find earliest EYO reaching significance and the best p at that earliest EYO
term_sig_summary <- sig_dat %>%
  dplyr::group_by(term) %>%
  dplyr::slice_min(EYO, with_ties = FALSE) %>%
  dplyr::summarise(
    earliest_sig_EYO  = EYO,
    min_p_at_earliest = p,
    .groups = "drop"
  )

# Choose 1 term per cluster:
#   priority = earliest_sig_EYO (ascending), then min_p_at_earliest (ascending)
picked_terms <- term_sig_summary %>%
  arrange(earliest_sig_EYO, min_p_at_earliest) %>%
#  distinct(cluster, .keep_all = TRUE) %>%
  pull(term)

##--- 4) Prepare heatmap table for the picked terms -------------------
plot_tab <- dat %>%
  filter(term %in% picked_terms) %>%
  # signed score: -log10(p_cap) with sign from direction
  mutate(
    signed_score = case_when(
      direction == "indianred3"     ~  +(-log2(p_cap)),  # higher in e4/4 = +red
      direction == "cornflowerblue" ~  -(-log2(p_cap)),  # lower  in e4/4 = -blue
      TRUE ~ NA_real_
    )
  )

# Order Y tracks: earlier EYO first; ties by stronger (lower) p at that earliest EYO
term_order <- term_sig_summary %>%
  arrange(earliest_sig_EYO, min_p_at_earliest) %>%
  pull(term)

plot_tab <- plot_tab %>%
  mutate(term = factor(term, levels = term_order))

##--- 5) Draw heatmap -------------------------------------------------
thr=-log2(0.0051)

# Diverging scale makes white near 0, blue for negatives, red for positives
p <- ggplot(plot_tab %>% filter(!is.na(signed_score))) + scale_y_discrete(limits = rev(levels(plot_tab$term))) +
  geom_tile(aes(x = EYO, y = term, fill = signed_score), width = 0.9, height = 0.9) +
#  scale_fill_gradient2(
#    low = "cornflowerblue", mid = "white", high = "darkred", midpoint = 0,
#    name = expression(paste("signed  ", -log[2], "(p)")),
#    limits = c(-max(-log2(plot_tab$p_cap), na.rm = TRUE),
#                max(-log2(plot_tab$p_cap), na.rm = TRUE))
#  ) +
  scale_fill_gradientn(
    colors = c("cornflowerblue", "white",
               "white", "darkred"),
    values = scales::rescale(c(
      -max(-log2(plot_tab$p_cap), na.rm = TRUE),
      -thr,
      thr,
      max(-log2(plot_tab$p_cap), na.rm = TRUE)
    )),
    name = expression(paste("signed  ", -log[2], "(p)")),
    limits = c(
      -max(-log2(plot_tab$p_cap), na.rm = TRUE),
       max(-log2(plot_tab$p_cap), na.rm = TRUE)
    )
  ) +
  labs(x = "EYO", y = NULL,
       title = "Significant protein assay intervals (e4/4 vs e3/3)",
       subtitle = "white-blue = lower in e4/4, white-red = higher") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

ggsave(file="./scatter/_Waterfall-NoImputedGenotypes_FixedSplineKnotsToOriginalPositions_withImputedSamples.pdf",plot=p,width=8.5,height=11)

# ---- STEP 11. Save the complete workspace image for provenance/reuse. ----
save.image("./scatter/_saved.image.RData")
