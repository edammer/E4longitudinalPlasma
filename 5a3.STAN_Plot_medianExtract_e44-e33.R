################################################################################################
# Section 5 longitudinal STAN pipeline - annotated copy
# Source file: 5a3.STAN_Plot_medianExtract_e44-e33.R
# Pipeline process 3/12: Model 1 median effect extraction for all protein assays
#
# Purpose: Extract the posterior median e4/e4 minus e3/e3 effect-size trajectory for protein
# assays from Model 1 fits.
#
# Input:  _numericMeta_3177_trait.RDS
# Input:  _full_3177_protein_dft.RDS
# Input:  simple.3177/name_match_table.rds
# Input:  simple.3177/*_stan_glm.rds
# Input:  plot_functions_20250718.R
# Output: 0727_medians_all_assays*.rds
# Output: 0727_medians_all_assays*.csv
#
# Major analysis steps in this script:
#   1. Restrict the outcome list to SomaScan protein assays.
#   2. Evaluate model-predicted e3/e3 and e4/e4 curves across the half-year EYO grid.
#   3. Compute the posterior median difference at each EYO point.
#   4. Parallelize extraction across assays and save a time-by-assay matrix for connectivity
#      analyses.
#
# Cleanup/annotation notes:
#   - This is a cleaned, commented copy of the uploaded script; analysis logic and
#     parameter values were not intentionally changed.
#   - Files were decoded from the uploaded Windows/CP1252 text and written as UTF-8.
#   - No explicit "not run below here" block was detected in this prefix-5 file set.
################################################################################################


# ----------------------------------------------------------------------------------------
# Load packages for Model 1 posterior median extraction.
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
#setwd("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/1_Ori_Model/")
setwd("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/")

##################### ------------ Read the Trait Data ------------ ##################
# Load the master traits data and the pep2pro data

# ----------------------------------------------------------------------------------------
# Load the primary-cohort trait and outcome matrices.
# ----------------------------------------------------------------------------------------
#BL_traits <- readRDS("/home/workspace/files/EBD/Shijia_B_Derived_Data/20250709/numericMeta_3177_trait.RDS")
BL_traits <- readRDS("./_numericMeta_3177_trait.RDS")
BL_traits$EYO<- (65.6 - BL_traits$age_at_visit)*(-1)
#BL_traits_pep<- readRDS("/home/workspace/files/EBD/Shijia_B_Derived_data/20250709/full_3177_protein_dft.RDS")
BL_traits_pep<- readRDS("./_full_3177_protein_dft.RDS")
BL_traits_pep$EYO<-BL_traits$EYO
BL_traits_pep$ApoE_Indicator<-BL_traits$ApoE_Indicator

min(BL_traits_pep$EYO) # -45.6     -31.53973
max(BL_traits_pep$EYO) # 24.4       23.72 previously
EYO_cut = length(seq(-46, 25, by = 0.5)) # 142   #prev 113 for DS


# ----------------------------------------------------------------------------------------
# Load the model-name lookup table.
# ----------------------------------------------------------------------------------------
name_match_table <- readRDS("./simple.3177/name_match_table.rds")
#name_module_label <- read.csv("../scatterplot_label_20250718.csv", header = T)
name_module_label.add<-cbind(name_match_table$CleanedName,name_match_table$OriginalName,name_match_table$OriginalName)
#colnames(name_module_label.add)<-colnames(name_module_label)
#name_module_label<-rbind(name_module_label,name_module_label.add)
name_module_label<-name_module_label.add


################## -------- Fit the STAN model -------- ###########

pep_names <- names(BL_traits_pep)[c(2:7334)]  #7348: 12 modules at end, then MMSE and cdr

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
#				   ceiling(max(BL_traits_pep$EYO, na.rm = T)),
#				   by=0.5) )


#source("z:/ShijiaBian/PlasmaProteomic/Code/CommonFunctions/plot_functions_20250718.R")   #202311235.R  changed from _Der_F to ApoE_Indicator

# ----------------------------------------------------------------------------------------
# Source the custom plotting/extraction helpers.
# ----------------------------------------------------------------------------------------
#source("F:/OneDrive - Emory/GV_HDS/stan/3.scatterplots/plot_functions_20250718.R")   #202311235.R  changed from _Der_F to ApoE_Indicator
source("./plot_functions_20250718.R")

process_one_peptide <- function(pep_name) {

# ----------------------------------------------------------------------------------------
# Worker function: read one saved Model 1 fit and compute posterior difference medians.
# ----------------------------------------------------------------------------------------
	clean_name = name_match_table$CleanedName[which(name_match_table$OriginalName == pep_name)]
	fil_stan <- file.path("./simple.3177/", paste(paste(clean_name, "_stan_glm",  ".rds", sep = "")))

	stan_BL <- readRDS(fil_stan)

	if (length(stan_BL) == 1) {
	  next
	}

	###### ----------- Plot 1: Two separate plots for carrier and non-carrier ------------- ############
	parameter_estimates <- rstan::extract(stan_BL$stanfit) # Extract the model estimates from STAN

# ----------------------------------------------------------------------------------------
# Extract posterior draws and construct factor-weight matrices.
# ----------------------------------------------------------------------------------------

	# biomarker ~ EYO_Spline_Linear + EYO_Spline_Cubic + Group + EYO_Spline_Linear * Group + EYO_Spline_Cubic
	# Add intercept to beta weights, 4000 * 6
	factorweights <- cbind(parameter_estimates$alpha, parameter_estimates$beta) # Add intercept to beta weights

	#Initialize contrast matrices
	eyo_step = seq(floor(min(BL_traits_pep$EYO, na.rm = T)),
				   ceiling(max(BL_traits_pep$EYO, na.rm = T)),

# ----------------------------------------------------------------------------------------
# Define the half-year EYO grid.
# ----------------------------------------------------------------------------------------
				   by=.5) #create a vector from x to y by the specified interval. This represents the range that we will use to plot our lines

	#Generate output matrices, nrow = 4,000 (stan iterations); ncol = 125, 125 is the number of eyo
	spline_noncarriers = matrix(0, nrow = nrow(factorweights), ncol = length(eyo_step)) #generate blank matrices for outputs for noncarriers
	spline_carriers = matrix(0, nrow = nrow(factorweights), ncol = length(eyo_step)) #generate blank matrices for outputs for carriers
	spline_diff = matrix(0, nrow = nrow(factorweights), ncol = length(eyo_step)) #generate blank matrices for output of differences

	contrasts_non_carriers = matrix(0, nrow = length(eyo_step), ncol = dim(factorweights)[2])
	contrasts_carriers = matrix(0, nrow = length(eyo_step), ncol = dim(factorweights)[2]) #Making a blank matrix for contrasts
	splinefit = rcspline.eval(BL_traits_pep$EYO, nk=3, norm=2, pc=FALSE, inclx=TRUE) # Redo the spline fit (if wanting to skip earlier portions, currently omitted)

	for (j in 1:length(eyo_step)) {

# ----------------------------------------------------------------------------------------
# Evaluate spline basis values at each EYO grid point.
# ----------------------------------------------------------------------------------------
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

	}

	#Begin the loop to go through every single iteration of Stan
	for (i in 1:nrow(factorweights)) {
		weights_temp = factorweights[i,] #for every iteration select the weights


# ----------------------------------------------------------------------------------------
# Compute posterior e4/e4 minus e3/e3 trajectories.
# ----------------------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------------------
# Retain the median difference as the assay effect-size trajectory.
# ----------------------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------------------
# Run parallel extraction over protein assays only.
# ----------------------------------------------------------------------------------------

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


  # return a named list with the objects you need afterwards
  list(
    median   = diff_lines$median,
    ok       = TRUE                 # handy flag for foreach result binding
  )
}


library(doParallel)
ncore <- max(1, parallel::detectCores() - 1)
cl     <- makeCluster(ncore)
registerDoParallel(cl)


library(foreach)

# ----------------------------------------------------------------------------------------
# Assemble the EYO-by-assay median effect-size matrix.
# ----------------------------------------------------------------------------------------

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
# Save the matrix for connectivity and waterfall analyses.
# ----------------------------------------------------------------------------------------
medians_all_pep    <- do.call(cbind, lapply(results[ok_idx], `[[`, "median"))


## Set dimnames of final data
colnames(medians_all_pep)<-pep_names[ok_idx]

eyo_step = seq(floor(min(BL_traits_pep$EYO, na.rm = T)),
				   ceiling(max(BL_traits_pep$EYO, na.rm = T)),
				   by=.5)
length(eyo_step)==nrow(medians_all_pep)  # TRUE
rownames(medians_all_pep)<-as.character(eyo_step)


# ######### ---- Write Final Outputs for Waterfall(next) step ---- ########
medians_all_pep_final <- medians_all_pep #[, -1]
fil_pep<-"./0727_medians_all_assays(5xSD_outliers_excluded).rds"
saveRDS(medians_all_pep_final, fil_pep)
fil_csv_pep <- "./0727_medians_all_assays(5xSD_outliers_excluded).csv"
write.csv(medians_all_pep_final, fil_csv_pep, row.names = TRUE)
