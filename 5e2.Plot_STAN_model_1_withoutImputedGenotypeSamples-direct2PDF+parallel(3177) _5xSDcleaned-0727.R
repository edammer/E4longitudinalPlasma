################################################################################################
# Section 5 longitudinal STAN pipeline - annotated copy
# Source file: 5e2.Plot_STAN_model_1_withoutImputedGenotypeSamples-direct2PDF+parallel(3177) _5xSDcleaned-0727.R
# Pipeline process 9/12: Model 1 no-imputation sensitivity plots and EYO-bin matrices
#
# Purpose: Generate longitudinal plots and posterior EYO matrices for the Model 1 sensitivity
# run that excludes imputed genotype samples.
#
# Input:  assays+traits_forPlots.RData
# Input:  noImpute.simple.rerun/name_match_table.rds
# Input:  noImpute.simple.rerun/*_stan_glm.rds
# Input:  plot_functions_20250718.R
# Output: noImpute.simple.rerun/scatter/99par_2on1_*.pdf
# Output: noImpute.simple.rerun/scatter/_99_par_diff_all_peptide*.rds/.csv
# Output: noImpute.simple.rerun/scatter/_99_par_diff_all_peptide_up_down_notation.csv
#
# Major analysis steps in this script:
#   1. Load the filtered no-imputation traits and outcome matrices.
#   2. Read no-imputation Model 1 Stan fits.
#   3. Evaluate e3/e3 and e4/e4 posterior trajectories across EYO.
#   4. Plot and export posterior differences, intervals, p-value matrices, and direction calls.
#
# Cleanup/annotation notes:
#   - This is a cleaned, commented copy of the uploaded script; analysis logic and
#     parameter values were not intentionally changed.
#   - Files were decoded from the uploaded Windows/CP1252 text and written as UTF-8.
#   - No explicit "not run below here" block was detected in this prefix-5 file set.
################################################################################################


# ----------------------------------------------------------------------------------------
# Load packages for no-imputation Model 1 plotting.
# ----------------------------------------------------------------------------------------
#' All "BL" stands for "Last measure"

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
setwd("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/")


# ----------------------------------------------------------------------------------------
# Retained commented block documents but does not rerun preprocessing here.
# ----------------------------------------------------------------------------------------
##################### ------------ Read the Trait Data ------------ ##################
# Load the master traits data and the pep2pro data
##BL_traits <- readRDS("/home/workspace/files/EBD/Shijia_B_Derived_Data/20250709/numericMeta_3177_trait.RDS")
#BL_traits <- readRDS("_numericMeta_3177_trait.RDS")
#BL_traits$EYO<- (65.6 - BL_traits$age_at_visit)*(-1)
##BL_traits_pep<- readRDS("/home/workspace/files/EBD/Shijia_B_Derived_data/20250709/full_3177_protein_dft.RDS")
#BL_traits_pep<- readRDS("_full_3177_protein_dft.RDS")
#BL_traits_pep$EYO<-BL_traits$EYO
#BL_traits_pep$ApoE_Indicator<-BL_traits$ApoE_Indicator
#
#
#
## 2. Remove samples with missing APOE.mapped
#na_apoe_idx <- which(is.na(BL_traits$APOE.mapped))
#message(sprintf("Removing %d sample(s) with NA APOE.mapped.", length(na_apoe_idx)))
##Removing 708 sample(s) with NA APOE.mapped.
#
#if (length(na_apoe_idx) > 0) {
#  na_sample_ids          <- rownames(BL_traits)[na_apoe_idx]
#  BL_traits <- BL_traits[-na_apoe_idx, ]
#  BL_traits_pep             <- BL_traits_pep[!BL_traits_pep$sample_id %in% na_sample_ids, ]
#}
#

load("assays+traits_forPlots.RData")
BL_traits_pep$ApoE_Indicator<-BL_traits$ApoE_Indicator


# 3. Verify alignment between metadata and protein matrix
stopifnot(all(BL_traits_pep$sample_id == BL_traits$sample_id))

# ----------------------------------------------------------------------------------------
# Load the filtered no-imputation traits and outcome objects.
# ----------------------------------------------------------------------------------------



min(BL_traits_pep$EYO) # -45.6     -31.53973
max(BL_traits_pep$EYO) # 24.4       23.72 previously
EYO_cut = length(seq(-46, 25, by = 0.5)) # 142   #prev 113 for DS

name_match_table <- readRDS("./noImpute.simple.rerun/name_match_table.rds")

#*** we do not have MEs now
##name_module_label <- read.csv("../scatterplot_label_20250718.csv", header = T)
name_module_label <- read.csv("./scatterplot_label_20250718.csv", header = T)
name_module_label.names<-colnames(name_module_label)
#name_module_label.add<-cbind(name_match_table$CleanedName,name_match_table$OriginalName,name_match_table$OriginalName)
#colnames(name_module_label.add)<-colnames(name_module_label)

# ----------------------------------------------------------------------------------------
# Load no-imputation name-match and plot-label tables.
# ----------------------------------------------------------------------------------------
#name_module_label<-as.data.frame(rbind(name_module_label,name_module_label.add))
name_module_label<-cbind(name_match_table$CleanedName,name_match_table$OriginalName,name_match_table$OriginalName)
colnames(name_module_label)<-name_module_label.names
name_module_label<-as.data.frame(name_module_label)

################## -------- Fit the STAN model -------- ###########

pep_names <- names(BL_traits_pep)[c(2:7346)]  #12 modules at end, no MMSE and cdr

#count = 1

# Used for storing the final t-statistics for the difference
#noncarrier_carrier_t_stats_all_pep <- matrix(, nrow = EYO_cut)
#colnames(noncarrier_carrier_t_stats_all_pep) = "Empty"
#
#noncarrier_carrier_p_value_all_pep <- matrix(, nrow = EYO_cut)
#colnames(noncarrier_carrier_p_value_all_pep) = "Empty"
#
#up_down_notation <- matrix(, nrow = length(pep_names), ncol = EYO_cut + 1)
#colnames(up_down_notation) = c("Pep", seq(floor(min(BL_traits_pep$EYO, na.rm = T)),

# ----------------------------------------------------------------------------------------
# Source custom plotting functions.
# ----------------------------------------------------------------------------------------
#				   ceiling(max(BL_traits_pep$EYO, na.rm = T)),
#				   by=0.5) )


# ----------------------------------------------------------------------------------------
# Worker function: read one no-imputation Model 1 fit and compute trajectory summaries.
# ----------------------------------------------------------------------------------------
#source("z:/ShijiaBian/PlasmaProteomic/Code/CommonFunctions/plot_functions_20250718.R")   #202311235.R  changed from _Der_F to ApoE_Indicator
source("./plot_functions_20250718.R")   #202311235.R  changed from _Der_F to ApoE_Indicator
#sink("Original_STAN_Model_Scatter_Plot.txt", append = TRUE)
#print(Sys.time())
#sink()


process_one_peptide <- function(pep_name) {
	clean_name = name_match_table$CleanedName[which(name_match_table$OriginalName == pep_name)]
	fil_stan <- file.path("./noImpute.simple.rerun/", paste(paste(clean_name, "_stan_glm",  ".rds", sep = "")))


# ----------------------------------------------------------------------------------------
# Extract posterior draws and build model coefficient matrices.
# ----------------------------------------------------------------------------------------
	stan_BL <- readRDS(fil_stan)

	if (length(stan_BL) == 1) {
	  #next  # only valid for/while loop
	  return(list(ok = FALSE))
	}

	###### ----------- Plot 1: Two separate plots for carrier and non-carrier ------------- ############
	parameter_estimates <- rstan::extract(stan_BL$stanfit) # Extract the model estimates from STAN

# ----------------------------------------------------------------------------------------
# Define EYO grid.
# ----------------------------------------------------------------------------------------

	# biomarker ~ EYO_Spline_Linear + EYO_Spline_Cubic + Group + EYO_Spline_Linear * Group + EYO_Spline_Cubic
	# Add intercept to beta weights, 4000 * 6
	factorweights <- cbind(parameter_estimates$alpha, parameter_estimates$beta) # Add intercept to beta weights

	#Initialize contrast matrices
	eyo_step = seq(floor(min(BL_traits_pep$EYO, na.rm = T)),
				   ceiling(max(BL_traits_pep$EYO, na.rm = T)),
				   by=.5) #create a vector from x to y by the specified interval. This represents the range that we will use to plot our lines

	#Generate output matrices, nrow = 4,000 (stan iterations); ncol = 125, 125 is the number of eyo
	spline_noncarriers = matrix(0, nrow = nrow(factorweights), ncol = length(eyo_step)) #generate blank matrices for outputs for noncarriers
	spline_carriers = matrix(0, nrow = nrow(factorweights), ncol = length(eyo_step)) #generate blank matrices for outputs for carriers

# ----------------------------------------------------------------------------------------
# Evaluate spline basis across the EYO grid.
# ----------------------------------------------------------------------------------------
	spline_diff = matrix(0, nrow = nrow(factorweights), ncol = length(eyo_step)) #generate blank matrices for output of differences

	contrasts_non_carriers = matrix(0, nrow = length(eyo_step), ncol = dim(factorweights)[2])
	contrasts_carriers = matrix(0, nrow = length(eyo_step), ncol = dim(factorweights)[2]) #Making a blank matrix for contrasts
	splinefit = rcspline.eval(BL_traits_pep$EYO, nk=3, norm=2, pc=FALSE, inclx=TRUE) # Redo the spline fit (if wanting to skip earlier portions, currently omitted)

	for (j in 1:length(eyo_step)) {
		tempfit = rcspline.eval(eyo_step[j], attr(splinefit,'knots'), norm=2, pc=FALSE, inclx=TRUE)       #put all of the values from our plotting interval into rcscpline to get the cubic term values
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

# ----------------------------------------------------------------------------------------
# Compute e3/e3, e4/e4, and difference posterior trajectories.
# ----------------------------------------------------------------------------------------

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


# ----------------------------------------------------------------------------------------
# Summarize intervals and posterior tail probabilities.
# ----------------------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------------------
# Generate combined observed-data and posterior-difference plots.
# ----------------------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------------------
# Start parallel plotting/extraction.
# ----------------------------------------------------------------------------------------
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
        ggsave(filename = paste0("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/noImpute.simple.rerun/scatter/99par_2on1_",gsub("(\\|)|(\\.)|(\\;)", "_",pep_name),".pdf"), plot = on_same_plot, device = cairo_pdf,
               width = 22.5, height = 10, units = "in", dpi = 300)

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

# ----------------------------------------------------------------------------------------
# Assemble EYO-by-assay summary matrices.
# ----------------------------------------------------------------------------------------
registerDoParallel(cl)


library(foreach)

# list of packages the workers must load
worker_pkgs <- c("tidyverse","rstanarm","Hmisc","openxlsx","rstan",
                 "gridExtra","ggpubr","ggplot2","grDevices","ggplotify")

results <- foreach(pep_name = pep_names,
                   .packages = worker_pkgs,
                   .export   = c("process_one_peptide",      # the function
                                 "diff_plot","scatter_plot", # custom plotters
                                 "BL_traits_pep","name_match_table",
                                 "name_module_label",        # big objects
                                 "EYO_cut","cairo_pdf"),     # scalars/funs
                   .errorhandling = "pass")  %dopar%  {
  process_one_peptide(pep_name)
}

stopImplicitCluster()   # or stopCluster(cl) if you built one explicitly


# Filter out failed iterations
ok_idx <- vapply(results, function(x) is.list(x) && !inherits(x, "error") && x$ok, TRUE)

# ----------------------------------------------------------------------------------------
# Save p-value, effect-size, and direction matrices.
# ----------------------------------------------------------------------------------------

noncarrier_carrier_t_stats_all_pep    <- do.call(cbind,
                                                 lapply(results[ok_idx], `[[`, "t_stat"))
noncarrier_carrier_p_value_all_pep    <- do.call(cbind,
                                                 lapply(results[ok_idx], `[[`, "p_value"))
up_down_notation                      <- do.call(rbind,
                                                 lapply(results[ok_idx], `[[`, "up_down")) #[ , -1]

## Set dimnames of final data
colnames(noncarrier_carrier_t_stats_all_pep)<-colnames(noncarrier_carrier_p_value_all_pep)<-rownames(up_down_notation)<-pep_names[ok_idx]

eyo_step = seq(floor(min(BL_traits_pep$EYO, na.rm = T)),
				   ceiling(max(BL_traits_pep$EYO, na.rm = T)),
				   by=.5)
length(eyo_step)==ncol(up_down_notation)  # TRUE
rownames(noncarrier_carrier_t_stats_all_pep)<-rownames(noncarrier_carrier_p_value_all_pep)<-colnames(up_down_notation)<-as.character(eyo_step)


## Get characteristics of waterfall
apply(up_down_notation,2,function(x) table(x))

sum(apply(up_down_notation,1,function(x) length(which(!is.na(x))))>0)
#1119 (on R v4.5.1, was 1124 before rerun v4.4.1) ... 881.. 886


# ============================================================
# firstSigGenoMidpoints
# For every peptide in BL_traits_pep columns 2:7346 (in that
# order), return the midpoint EYO value of the FIRST (leftmost)
# contiguous interval where the genotype Bayesian p-value
# is < 0.0051.
#   "constitutive" ? the ENTIRE EYO range is significant
#   NA             ? no EYO bin reaches significance
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

#firstSigGenoMidpoints: 18 constitutive | 1102 with first-sig EYO | 6225 NA


all_pep_ordered[which(firstSigGenoMidpoints=="constitutive")]
# [1] "NEFL|P07196"                         "LRRN1|Q6UXK5^SL025922@seq.11293.14"  "TBCA|O75347"
# [4] "ARL2|P36404"                         "CTF1|Q16619^SL002783@seq.13732.79"   "VPS29|Q9UBQ0"
# [7] "PHGDH|O43175^SL018791@seq.15548.35"  "CDA|P32320"                          "FAM50A|Q14320"
#[10] "BCDIN3D|Q7Z5W3"                      "ST8SIA1|Q92185^SL022499@seq.21508.7" "PPM1G|O15355"
#[13] "FOXO1|Q12778"                        "SPC25|Q9HBM1"                        "ZW10|O43264"
#[16] "MENT|Q9BUN1"                         "S100A13|Q99584"                      "DCUN1D5|Q9BTE7"


# ######### ---- Write Final Outputs for Waterfall(next) step ---- ########
noncarrier_carrier_t_stats_all_pep_final <- noncarrier_carrier_t_stats_all_pep #[, -1]
#fil_pep <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files", paste(paste("99_par_diff_all_peptide", ".rds", sep = "")))
fil_pep<-"./noImpute.simple.rerun/scatter/_99_par_diff_all_peptide.rds"
saveRDS(noncarrier_carrier_t_stats_all_pep_final, fil_pep)
#fil_csv_pep <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files",paste(paste("99_par_diff_all_peptide", ".csv", sep = "")))
fil_csv_pep <- "./noImpute.simple.rerun/scatter/_99_par_diff_all_peptide.csv"
write.csv(noncarrier_carrier_t_stats_all_pep_final, fil_csv_pep, row.names = TRUE)

noncarrier_carrier_p_value_all_pep_final <- noncarrier_carrier_p_value_all_pep #[, -1]
#fil_pep <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files",paste(paste("99_par_diff_all_peptide_p_value", ".rds", sep = "")))
fil_pep<-"./noImpute.simple.rerun/scatter/_99_par_diff_all_peptide_p_value.rds"
saveRDS(noncarrier_carrier_p_value_all_pep_final, fil_pep)
#fil_csv_pep <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files",paste(paste("99_par_diff_all_peptide_p_value", ".csv", sep = "")))
fil_csv_pep <- "./noImpute.simple.rerun/scatter/_99_par_diff_all_peptide_p_value.csv"
write.csv(noncarrier_carrier_p_value_all_pep_final, fil_csv_pep, row.names = TRUE)

#fil_csv_up_down_notation <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files",paste(paste("99_par_diff_all_peptide_up_down_notation", ".csv", sep = "")))
fil_csv_up_down_notation<-"./noImpute.simple.rerun/scatter/_99_par_diff_all_peptide_up_down_notation.csv"
write.csv(up_down_notation, fil_csv_up_down_notation, row.names = TRUE)

#sink("Result/20240122/Log_File/Original_STAN_Model_Scatter_Plot.txt", append = TRUE)
#print(Sys.time())
#sink()


## Reprocess updated function -- only 12 MEs
# for (pep_name in pep_names[7334:7345]) process_one_peptide(pep_name)
