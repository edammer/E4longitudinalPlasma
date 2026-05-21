################################################################################################
# Section 5 longitudinal STAN pipeline - annotated copy
# Source file: 5c2.Plot_STAN_model_2_20250718-direct2PDF+parallel(3177)_5xSDcleaned-0727+plotFixes.R
# Pipeline process 6/12: Model 2 plotting and sex-stratified posterior summaries in the primary cohort
#
# Purpose: Post-process Model 2.2 fits into sex-stratified APOE trajectories and posterior EYO-
# bin summaries.
#
# Input:  _numericMeta_3177_trait.RDS
# Input:  _full_3177_protein_dft.RDS
# Input:  sexInt.3177/name_match_table.rds
# Input:  sexInt.3177/*_with_Sex_Interaction_stan_glm.rds
# Input:  plot_functions_20250718.R
# Output: sexInt.3177/scatter/*.pdf
# Output: sexInt.3177/scatter/*p_value*.rds/.csv
# Output: sexInt.3177/scatter/*up_down_notation*.csv
#
# Major analysis steps in this script:
#   1. Load the saved sex-interaction model fits and validate coefficient order.
#   2. Evaluate posterior curves for e3/e3 and e4/e4 groups across sex-specific strata.
#   3. Compute APOE-associated differences and posterior intervals across EYO.
#   4. Generate direct-to-PDF longitudinal plots with updated plot fixes.
#   5. Write EYO-by-assay matrices used to identify sex-dependent patterns.
#
# Cleanup/annotation notes:
#   - This is a cleaned, commented copy of the uploaded script; analysis logic and
#     parameter values were not intentionally changed.
#   - Files were decoded from the uploaded Windows/CP1252 text and written as UTF-8.
#   - No explicit "not run below here" block was detected in this prefix-5 file set.
################################################################################################


# ----------------------------------------------------------------------------------------
# Load packages for Model 2.2 plotting and posterior extraction.
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
#setwd("z:/ShijiaBian/PlasmaProteomic/Result/20250727/sexInt.3177/scatter/")
setwd("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/")

##################### ------------ Re-process the Trait Data ------------ ##################
#	load("z:/EBD/APOE_homozygote_2ndFinal_cleanDats+traits(3177.CT33only_4199.33inADcontinuum).RData")
#	ls()
#	#"cleanDat.3177. "cleanDat.4199"  "MEs.3177" "MEs.4199" "numericMeta.3177" "numericMeta.4199"
#	full_3177_protein_df <- as.data.frame(rbind(cleanDat.3177, t(MEs.3177)))
#
#	## Data cleaning - 5SD from mean max within protein, then max 20% per genotype group NA; add MMSE and cdr
#	clean_mat <- as.data.frame(t(apply(full_3177_protein_df, 1, function(row) {
#	    ## winsorise at +/- 5 SD
#	  z   <- abs(row - mean(row, na.rm = TRUE))
#	  row[z > 5 * sd(row, na.rm = TRUE)] <- NA                 # outliers -> NA
#
#	  row }                                              # otherwise return the cleaned row
#	)))
#	## restore dimnames -------------------------------------------------------
#	rownames(clean_mat) <- rownames(full_3177_protein_df)
#	colnames(clean_mat) <- colnames(full_3177_protein_df)
#
#
#
#	## grouping factor – one value per column in the expression matrix
#	grp <- numericMeta.3177$APOE.mapped.predicted
#	stopifnot(length(grp) == ncol(full_3177_protein_df))   # safety
#
#	## pre-compute the column indices that belong to each group
#	idx_by_grp <- split(seq_along(grp), grp)
#
#	## rows that fail the “<=20% NA in every group” rule  ----------------------
#	bad_row_idx <- which(
#	  apply(clean_mat, 1, function(v) {
#	    any(                                           # if *any* group exceeds 20% NA
#	      vapply(idx_by_grp,
#	             function(ix) mean(is.na(v[ix])) > 0.20,
#	             logical(1))
#	    )
#	  })
#	)
#
#	bad_row_idx        # numeric vector of offending row indices                                    # -> row is discarded
#	#named integer(0)
#
#	## If any found, remove by %in% rownames(clean_mat)[bad_row_idx]
#	if(length(bad_row_idx)>0) clean_mat<-clean_mat[which(!rownames(clean_mat) %in% rownames(clean_mat)[bad_row_idx]),]
#
#
#	## Add MMSE, cdr as additional outcomes
#	full_3177_protein_df<-as.data.frame(t(clean_mat))
#	full_3177_protein_df$MMSE<-numericMeta.3177$MMSE
#	full_3177_protein_df$cdr<-numericMeta.3177$cdr
#
#	full_3177_protein_df$sample_id <- rownames(full_3177_protein_df)
#
#
#
#	# set derived traits
#	numericMeta.3177$ApoE_Indicator<-0
#	numericMeta.3177$ApoE_Indicator[numericMeta.3177$APOE.mapped.predicted=="e4/e4"]<-1
#	table(numericMeta.3177$ApoE_Indicator)
#	#   0   1
#	#2764 413
#	numericMeta.3177$EY0 <- (65.6-numericMeta.3177$age_at_visit)*(-1)
#	range(numericMeta.3177$EY0)
#	# -45.6 24.4
#	numericMeta_3177_trait <- numericMeta.3177
#	numericMeta_3177_trait$sample_id <- rownames(numericMeta_3177_trait)

## we do not need to reprocess data, instead, load the RDS for numericMeta_3177_trait and full_3177_protein_df:
numericMeta_3177_trait<-readRDS("_numericMeta_3177_trait.RDS")
full_3177_protein_df<-readRDS("_full_3177_protein_dft.RDS")

which(!full_3177_protein_df$sample_id==numericMeta_3177_trait$sample_id)
# integer(0)

########################################################################
#' @param x A data frame (or any object with names)
#' @return The same object, but with names cleaned
sanitize_names <- function(x) {
  #Capture existing names
  nm <- names(x)
  if (is.null(nm)) return(x)

  # Replace anything that's not a letter or digit with
  clean_nm <- gsub("[^[:alnum:]]", "_", nm)
  # (Optional) collapse repeated underscores:
  clean_nm <- gsub("_+", "_", clean_nm)
  # (Optional) trim leading/trailing underscores:
  clean_nm <- gsub("^_|_$", "", clean_nm)

# ----------------------------------------------------------------------------------------
# Load or reuse the preprocessed primary cohort matrices.
# ----------------------------------------------------------------------------------------

  names(x) <- clean_nm
  x
}
##########################################################################

# For the protein file, make the person_id to be the first column
protein_df <- full_3177_protein_df %>%
  select(sample_id, everything())  #  Moves 'AAA' to the first column


# Load the master traits data and the pep2pro data
#BL_traits <- readRDS("/home/workspace/files/EBD/Shijia_B_Derived_Data/20250709/numericMeta_3177_trait.RDS")
BL_traits <- numericMeta_3177_trait  #readRDS("z:/EBD/Shijia_B_Derived_Data/20250709/numericMeta_3177_trait.RDS")
BL_traits$EYO<- (65.6 - BL_traits$age_at_visit)*(-1)
#BL_traits_pep<- readRDS("/home/workspace/files/EBD/Shijia_B_Derived_data/20250709/full_3177_protein_dft.RDS")
BL_traits_pep<- protein_df  #readRDS("z:/EBD/Shijia_B_Derived_data/20250709/full_3177_protein_dft.RDS")
BL_traits_pep$EYO<-BL_traits$EYO
BL_traits_pep$ApoE_Indicator<-BL_traits$ApoE_Indicator

min(BL_traits_pep$EYO) # -45.6     -31.53973
max(BL_traits_pep$EYO) # 24.4       23.72 previously
EYO_cut = length(seq(-46, 25, by = 0.5)) # 143   #prev 113 for DS


# ----------------------------------------------------------------------------------------
# Load primary traits and outcomes for plotting.
# ----------------------------------------------------------------------------------------
#name_match_table <- readRDS("z:/ShijiaBian/PlasmaProteomic/Result/20250727/sexInt.3177/name_match_table.rds")
name_match_table <- readRDS("./sexInt.3177/name_match_table.rds")
#name_module_label <- read.csv("z:/ShijiaBian/PlasmaProteomic/Result/20250709/scatterplot_label_20250718.csv", header = T)
name_module_label <- read.csv("./scatterplot_label_20250718.csv", header = T)
name_module_label.add<-cbind(name_match_table$CleanedName,name_match_table$OriginalName,name_match_table$OriginalName)
colnames(name_module_label.add)<-colnames(name_module_label)
name_module_label<-rbind(name_module_label,name_module_label.add)


################## -------- Fit the STAN model -------- ###########

pep_names <- names(BL_traits_pep)[c(2:7348)]  #12 modules at end

#count = 1

# Used for storing the final t-statistics for the difference
#noncarrier_carrier_t_stats_all_pep <- matrix(, nrow = EYO_cut)
#colnames(noncarrier_carrier_t_stats_all_pep) = "Empty"
#
#noncarrier_carrier_p_value_all_pep <- matrix(, nrow = EYO_cut)

# ----------------------------------------------------------------------------------------
# Load Model 2 name-match and label tables.
# ----------------------------------------------------------------------------------------
#colnames(noncarrier_carrier_p_value_all_pep) = "Empty"
#
#up_down_notation <- matrix(, nrow = length(pep_names), ncol = EYO_cut + 1)
#colnames(up_down_notation) = c("Pep", seq(floor(min(BL_traits_pep$EYO, na.rm = T)),
#				   ceiling(max(BL_traits_pep$EYO, na.rm = T)),
#				   by=0.5) )


#source("z:/ShijiaBian/PlasmaProteomic/Code/CommonFunctions/plot_functions_20250718.R")   #202311235.R  changed from _Der_F to ApoE_Indicator
source("./plot_functions_20250718.R")
#sink("Original_STAN_Model_Scatter_Plot.txt", append = TRUE)
#print(Sys.time())
#sink()


process_one_peptide <- function(pep_name) {

  ## ---------------------------------------------------------------------
  ##  0)  Read the Stan-fit that was created with the 9-parameter model
  ## ---------------------------------------------------------------------
  clean_name <- name_match_table$CleanedName[
                  match(pep_name, name_match_table$OriginalName)]

  f_stan <- file.path("./sexInt.3177/",  # "z:", "ShijiaBian", "PlasmaProteomic", "Result", "20250727","sexInt.3177",
                      sprintf("%s_with_Sex_Interaction_stan_glm.rds", clean_name))

  stan_reg <- readRDS(f_stan)
  if (length(stan_reg) == 1L)           # fit failed - skip peptide
    return(list(ok = FALSE))

  ## quick sanity-check – make sure the coefficients are in the
  ## (Intercept, EYO_L, EYO_C, ApoE, L:ApoE, C:ApoE, Sex, L:Sex, C:Sex) order

# ----------------------------------------------------------------------------------------
# Source plotting helpers.
# ----------------------------------------------------------------------------------------
  trm <- broom.mixed::tidy(stan_reg)$term
  exp_order <- c("(Intercept)",
                 "EYO_Spline_Linear", "EYO_Spline_Cubic",
                 "ApoE_Indicator",
                 "Sex",
                 "EYO_Spline_Linear:ApoE_Indicator",

# ----------------------------------------------------------------------------------------
# Worker function: read one Model 2.2 fit and compute sex-stratified trajectories.
# ----------------------------------------------------------------------------------------
                 "EYO_Spline_Cubic:ApoE_Indicator",
                 "EYO_Spline_Linear:Sex",
                 "EYO_Spline_Cubic:Sex")
  if (!identical(trm[1:9], exp_order))
      stop("Coefficient order of saved model is not what the code expects!")

  ## ---------------------------------------------------------------------
  ##  1)  Pull out posterior draws
  ## ---------------------------------------------------------------------
  post  <- rstan::extract(stan_reg$stanfit)
  wmat  <- cbind(post$alpha, post$beta)       # 4000 × 9  (iterations × beta)

  ## ---------------------------------------------------------------------
  ##  2)  Build the design matrix  (4 scenarios × 9 beta  for every EYO grid)
  ##       – Sex 0/1   ×  ApoE (Indicator) 0/1
  ## ---------------------------------------------------------------------
  eyo_grid <- seq(floor(min(BL_traits_pep$EYO, na.rm = TRUE)),
                  ceiling(max(BL_traits_pep$EYO, na.rm = TRUE)),
                  by = .5)

  nk <- 3                                # the same number of knots used for the fit
  splinefit <- spl_all <- Hmisc::rcspline.eval(BL_traits_pep$EYO, nk = nk,
                                  norm = 2, pc = FALSE, inclx = TRUE)

# ----------------------------------------------------------------------------------------
# Validate expected coefficient order before constructing contrasts.
# ----------------------------------------------------------------------------------------
  knots   <- attr(spl_all, "knots")

  # spline basis for the plotting grid – this is where eyo_L/C come in
  eyo_spline <- rcspline.eval(eyo_grid,
                              attr(splinefit, "knots"),  # reuse the knots
                              norm  = 2,
                              pc    = FALSE,
                              inclx = TRUE)
  eyo_L <- eyo_spline[, 1]   # linear-spline part
  eyo_C <- eyo_spline[, 2]   # cubic-spline part

  make_row <- function(L, C, apo, sex) {
    c(1,                    # intercept
      L, C,                 # main EYO spline bases
      apo,                  # ApoE_Indicator (0/1)
      sex,                  # Sex (0 = F, 1 = M)
      L*apo, C*apo,         # interactions with ApoE
      L*sex, C*sex)         # interactions with Sex
  }

  # Four combinations:             (apo, sex)
  combos <- list(
    F_NonCarrier = c(ApoE_Indicator = 0, Sex = 0),  # Female, non-carrier
    F_Carrier    = c(ApoE_Indicator = 1, Sex = 0),  # Female, carrier
    M_NonCarrier = c(ApoE_Indicator = 0, Sex = 1),  # Male,   non-carrier
    M_Carrier    = c(ApoE_Indicator = 1, Sex = 1)   # Male,   carrier
  )

  Xgrid <- array(NA_real_,
                 dim = c(length(eyo_grid), length(combos), 9),
                 dimnames = list(NULL, names(combos), NULL))

  for (j in seq_along(eyo_grid)) {

# ----------------------------------------------------------------------------------------
# Create the EYO grid and posterior contrast matrices.
# ----------------------------------------------------------------------------------------
    tmp <- Hmisc::rcspline.eval(eyo_grid[j], knots,
                                norm = 2, pc = FALSE, inclx = TRUE)
    L <- tmp[1, 1];  C <- tmp[1, 2]
    for (k in seq_along(combos))
      Xgrid[j, k, ] <- make_row(L, C, combos[[k]][1], combos[[k]][2])
  }

  ## ---------------------------------------------------------------------
  ##  3)  Evaluate posterior curves |  build CredInt for every curve
  ## ---------------------------------------------------------------------
  nsamp <- nrow(wmat)
  ngr   <- length(combos)
  ny    <- length(eyo_grid)

  # 3-D array : sample × grid × group
  mu_draws <- array(NA_real_,
                  dim       = c(nsamp, ny, ngr),
                  dimnames  = list(NULL,         # posterior draw
                                   NULL,         # EYO grid point
                                   names(combos) # array names
                                   ))

  ## all posterior draws as arrays / matrices

# ----------------------------------------------------------------------------------------
# Evaluate spline basis values across the EYO grid.
# ----------------------------------------------------------------------------------------
  parameter_estimates <- rstan::extract(stan_reg$stanfit)

  # fill mu_draws
  build_Xg <- function(apo, sex, L, C) {
    cbind(L, C,              # EYO spline bases    #leaving out 1, (intercept) first column
          apo,               # ApoE_Indicator
          sex,               # Sex
          L*apo, C*apo,      # interactions with ApoE
          L*sex, C*sex)      # interactions with Sex
  }

  ngr   <- length(combos)
  ny    <- length(eyo_grid)
  nsamp <- nrow(parameter_estimates$beta)

  mu_draws <- array(NA_real_,
                    dim = c(nsamp, ny, ngr),
                    dimnames = list(NULL, NULL, names(combos)))

  ## loop over the four sex × ApoE strata
  for (g in seq_len(ngr)) {
    apo <- combos[[g]]["ApoE_Indicator"]
    sex <- combos[[g]]["Sex"]

    Xg  <- build_Xg(apo, sex, eyo_L, eyo_C)          # ny × 9 (fresh each turn)

    # posterior means:  (iter × 9)  %*%  (9 × ny)  ->  iter × ny
    mu_tmp <- parameter_estimates$beta %*% t(Xg)

# ----------------------------------------------------------------------------------------
# Build posterior predictions for sex/genotype strata.
# ----------------------------------------------------------------------------------------

    # add the intercept for every iteration
    mu_draws[ , , g] <- sweep(mu_tmp, 1, parameter_estimates$alpha, "+")
  }

  # Compute the 0.5 % / 50 % / 99.5 % quantiles for each group
  Q <- lapply(seq_len(dim(mu_draws)[3]),            # loop over groups
              function(g) {
                apply(mu_draws[ , , g],             # iter × grid
                      2,                            # 2 = apply over columns = grid-points
                      quantile,
                      probs = c(.005, .5, .995))    # becomes 3 × grid matrix
              })
  names(Q) <- c("F_NonCarrier", "F_Carrier",
                "M_NonCarrier", "M_Carrier")
  # Q is now a list with names "F_NonCarrier", "F_Carrier", ...
  # Q[[g]] is 3 × ny   (rows: lower, median, upper)

  ## ---------------------------------------------------------------------
  ##  4)  differences
  ##       – carrier – non-carrier  (within each sex **and overall**)
  ##       – male – female          (within each ApoE stratum **and overall**)
  ## ---------------------------------------------------------------------

  ## helper that averages the two sex-specific (or ApoE-specific) curves
  mean2 <- function(a, b) (a + b) / 2

  diff_draws <- list(
    ## within–sex contrasts  (what you already had)
    Car_vs_Non_F = mu_draws[ , , "F_Carrier"]    - mu_draws[ , , "F_NonCarrier"],
    Car_vs_Non_M = mu_draws[ , , "M_Carrier"]    - mu_draws[ , , "M_NonCarrier"],
    M_vs_F_Non   = mu_draws[ , , "M_NonCarrier"] - mu_draws[ , , "F_NonCarrier"],
    M_vs_F_Car   = mu_draws[ , , "M_Carrier"]    - mu_draws[ , , "F_Carrier"],

    ## *** NEW overall contrasts ***
    ## (i) carrier – non-carrier    averaged over sex
    Car_vs_Non_All = mean2(mu_draws[ , , "F_Carrier"],
                           mu_draws[ , , "M_Carrier"])     -
                     mean2(mu_draws[ , , "F_NonCarrier"],
                           mu_draws[ , , "M_NonCarrier"]),

    ## (ii) male – female           averaged over ApoE status

# ----------------------------------------------------------------------------------------
# Compute posterior intervals and p-values for contrasts.
# ----------------------------------------------------------------------------------------
    M_vs_F_All     = mean2(mu_draws[ , , "M_Carrier"],
                           mu_draws[ , , "M_NonCarrier"]) -
                     mean2(mu_draws[ , , "F_Carrier"],
                           mu_draws[ , , "F_NonCarrier"])
  )

  qfun <- function(x) apply(x, 2, stats::quantile, probs = c(.005, .5, .995), names = FALSE)

  ## convert every element of diff_draws to its 0.5 %, 50 %, 99.5 % quantiles
  Qdiff <- lapply(diff_draws, qfun)   # each is a 3 × n_grid matrix


  ## pick one (e.g. F Carrier – F NonCarrier) for the “single-ribbon” demo
  diff_main_geno <- Qdiff$Car_vs_Non_F

#  ## keep using whichever one you like for the ribbon demo
#  diff_main_geno <- Qdiff$Car_vs_Non_All   # example: overall carrier – non-carrier (averaged M and F)

  ### accompanying t-statistics / tail-probs for that same contrast
  #tstat  <- apply(diff_draws$Car_vs_Non_All, 2, function(x) t.test(x)$statistic)
  #pvalue <- apply(diff_draws$Car_vs_Non_All, 2, function(x) {
  #                 p <- mean(x < 0);  pmin(p, 1 - p) })

  ## ---------------------------------------------------------------------
  ##  5)  build data.frames that the plotting helpers expect
  ## ---------------------------------------------------------------------
  ribbon_F <- data.frame(
        eyo        = eyo_grid,
        lower_non  = Q[[ "F_NonCarrier"]][1, ],
        median_non = Q[[ "F_NonCarrier"]][2, ],
        upper_non  = Q[[ "F_NonCarrier"]][3, ],
        lower_car  = Q[[ "F_Carrier"   ]][1, ],
        median_car = Q[[ "F_Carrier"   ]][2, ],
        upper_car  = Q[[ "F_Carrier"   ]][3, ])

  ribbon_M <- data.frame(
        eyo        = eyo_grid,
        lower_non  = Q[[ "M_NonCarrier"]][1, ],
        median_non = Q[[ "M_NonCarrier"]][2, ],
        upper_non  = Q[[ "M_NonCarrier"]][3, ],
        lower_car  = Q[[ "M_Carrier"   ]][1, ],
        median_car = Q[[ "M_Carrier"   ]][2, ],
        upper_car  = Q[[ "M_Carrier"   ]][3, ])

  diff_df_geno <- data.frame(
        eyo   = eyo_grid,
        lower = diff_main_geno[1, ],
        median= diff_main_geno[2, ],
        upper = diff_main_geno[3, ])


  ribbon_nonSex <- data.frame(
    eyo        = eyo_grid,
    lower_F    = Q[["F_NonCarrier"]][1, ],  median_F = Q[["F_NonCarrier"]][2, ],

# ----------------------------------------------------------------------------------------
# Prepare observed sample data for plotting.
# ----------------------------------------------------------------------------------------
    upper_F    = Q[["F_NonCarrier"]][3, ],
    lower_M    = Q[["M_NonCarrier"]][1, ],  median_M = Q[["M_NonCarrier"]][2, ],
    upper_M    = Q[["M_NonCarrier"]][3, ])

  ribbon_carSex <- data.frame(
    eyo        = eyo_grid,
    lower_F    = Q[["F_Carrier"]][1, ],     median_F = Q[["F_Carrier"]][2, ],
    upper_F    = Q[["F_Carrier"]][3, ],
    lower_M    = Q[["M_Carrier"]][1, ],     median_M = Q[["M_Carrier"]][2, ],
    upper_M    = Q[["M_Carrier"]][3, ])

#  diff_main_sex <- Qdiff$M_vs_F_All           # produced earlier in §4  - averaged both genotypes
  diff_main_sex <- Qdiff$M_vs_F_Non            # to match prior output, only show M vs F difference in e3/e3 (better powered)
  diff_df_sex   <- data.frame(
          eyo   = eyo_grid,
          lower = diff_main_sex[1, ],
          median= diff_main_sex[2, ],
          upper = diff_main_sex[3, ])

  ## ---------------------------------------------------------------------
  ##  6)  scatter-data for the current peptide
  ## ---------------------------------------------------------------------
  plot_dat <- BL_traits_pep |>
              dplyr::select(sample_id, EYO, !!sym(pep_name), ApoE_Indicator) #|>
#              tidyr::drop_na()
  names(plot_dat)[3] <- "Pep"

  # For first DiffPlot
  bar_df.geno.pre<-as.data.frame(diff_df_geno)
  bar_df.geno.pre$same_sign = FALSE
  bar_df.geno.pre$same_sign[sign(diff_df_geno$lower) == sign(diff_df_geno$upper)] = TRUE
  bar_df.geno = bar_df.geno.pre[bar_df.geno.pre$same_sign == TRUE, ]
  min_value.geno = if(pep_name %in% colors()) { min(-0.001, min(diff_df_geno$lower))*1.3 } else { floor(min(diff_df_geno$lower)) }  #min range -0.001 to + 0.001 if a color (ME), else -1 to +1 (or greater range)
                   # floor(min(ribbon_F$lower_non, ribbon_F$lower_car, ribbon_M$lower_non, ribbon_M$lower_car))
  max_value.geno = if(pep_name %in% colors()) { max(0.001, max(diff_df_geno$upper))*1.3 } else { ceiling(max(diff_df_geno$upper)) }
                   # ceiling(max(ribbon_F$upper_non, ribbon_F$upper_car, ribbon_M$upper_non, ribbon_M$upper_car))
  right_limit_EYO.geno = 24
  if (nrow(bar_df.geno) > 0) {
    bar_df.geno$Height = max_value.geno
    bar_df.geno$Text_Height = (max_value.geno - min_value.geno) * 0.05
    First_Annotate.geno <- bar_df.geno$eyo[1]
    Last_Annotate.geno <- bar_df.geno$eyo[length(bar_df.geno$eyo)]

    Annotate_color.geno <- "cornflowerblue"
    if (sum(bar_df.geno$median > 0)/length(bar_df.geno$median) > 0.5) { # majortiy of the median > 0
      Annotate_color.geno <- "indianred3"
    }

    Last_Annotate_Number.geno = Last_Annotate.geno + 0.5
    if (Last_Annotate_Number.geno > right_limit_EYO.geno){
      right_limit_EYO.geno = 26 # give more space to annotate the last annotation
    }
  } else {

# ----------------------------------------------------------------------------------------
# Generate sex-stratified trajectory and difference plots.
# ----------------------------------------------------------------------------------------
    First_Annotate.geno = Inf
    Last_Annotate.geno = Inf
    Last_Annotate_Number.geno = Inf
    Annotate_color.geno = ""
  }

  # For second DiffPlot
  bar_df.sex.pre<-as.data.frame(diff_df_sex)
  bar_df.sex.pre$same_sign = FALSE
  bar_df.sex.pre$same_sign[sign(diff_df_sex$lower) == sign(diff_df_sex$upper)] = TRUE
  bar_df.sex = bar_df.sex.pre[bar_df.sex.pre$same_sign == TRUE, ]
  min_value.sex = if(pep_name %in% colors()) { min(-0.001, min(diff_df_sex$lower))*1.3 } else { floor(min(diff_df_sex$lower)) }  #min range -0.001 to + 0.001 if a color (ME), else -1 to +1 (or greater range)
  max_value.sex = if(pep_name %in% colors()) { max(0.001, max(diff_df_sex$upper))*1.3 } else { ceiling(max(diff_df_sex$upper)) }
  right_limit_EYO.sex = 24
  if (nrow(bar_df.sex) > 0) {
    bar_df.sex$Height = max_value.sex #- 0.02*(abs(max_value.sex-min_value.sex))  # adjusted to keep in frame
    bar_df.sex$Text_Height = (max_value.sex - min_value.sex) * 0.05
    First_Annotate.sex <- bar_df.sex$eyo[1]
    Last_Annotate.sex <- bar_df.sex$eyo[length(bar_df.sex$eyo)]

    Annotate_color.sex <- "#6b8e23"  #"cornflowerblue"
    if (sum(bar_df.sex$median > 0)/length(bar_df.sex$median) > 0.5) { # majortiy of the median > 0
      Annotate_color.sex <- "#ff8c00"  #"indianred3"
    }

    Last_Annotate_Number.sex = Last_Annotate.sex + 0.5
    if (Last_Annotate_Number.sex > right_limit_EYO.sex){
      right_limit_EYO.sex = 26 # give more space to annotate the last annotation
    }
  } else {
    First_Annotate.sex = Inf
    Last_Annotate.sex = Inf
    Last_Annotate_Number.sex = Inf
    Annotate_color.sex = ""
  }

  ## ---------------------------------------------------------------------
  ##  7)  PLOTS
  ## ---------------------------------------------------------------------
  yl  <- name_module_label$y.Axis.Label[
           match(clean_name, name_module_label$Raw.File.Name)]
  ttl <- name_module_label$Title.Label[
           match(clean_name, name_module_label$Raw.File.Name)]

  plot_dat$Sex.int<-BL_traits$Sex.int


  plot_dat.F<-plot_dat[BL_traits$Sex.int==0,]
  plot_dat.F.count_before_removal_outlier <- nrow(plot_dat.F)
  temp_dat.F.plot <- plot_dat.F %>%
    mutate(zRT = scale(Pep)[,1]) %>%
    filter(between(zRT, -3, +3))
  plot_dat.F.count_post_removal_outlier <- nrow(temp_dat.F.plot)
  number_removed.F <- plot_dat.F.count_before_removal_outlier - plot_dat.F.count_post_removal_outlier  # record number of removed outliers

  scatter_plt_F <- scatter_plot_single(
        x             = EYO,
        y             = Pep,
        df            = temp_dat.F.plot,  # plot_dat[BL_traits$Sex.int==0,],
        protein_name  = data.frame(yaxis_label = yl,
                                   title_label = ttl),
        number_removed= number_removed.F,  # 0L,   # 0L if already dropped outliers
        CI_Percent    = "99%",
        combined_lines= ribbon_F,
        within        = "Female",
        compare       = "e4/4 vs e3/3",
        label         = c("e3/3", "e4/4"),
        by_group      = temp_dat.F.plot$ApoE_Indicator)  # plot_dat$ApoE_Indicator[BL_traits$Sex.int==0])


  plot_dat.M<-plot_dat[BL_traits$Sex.int==1,]
  plot_dat.M.count_before_removal_outlier <- nrow(plot_dat.M)
  temp_dat.M.plot <- plot_dat.M %>%
    mutate(zRT = scale(Pep)[,1]) %>%
    filter(between(zRT, -3, +3))
  plot_dat.M.count_post_removal_outlier <- nrow(temp_dat.M.plot)
  number_removed.M <- plot_dat.M.count_before_removal_outlier - plot_dat.M.count_post_removal_outlier  # record number of removed outliers

  scatter_plt_M <- scatter_plot_single(
        x             = EYO,
        y             = Pep,
        df            = temp_dat.M.plot,  # plot_dat[BL_traits$Sex.int==1,],
        protein_name  = data.frame(yaxis_label = yl,
                                   title_label = ttl),
        number_removed= number_removed.M,  # 0L,  # already dropped outliers
        CI_Percent    = "99%",
        combined_lines= ribbon_M,
        within        = "Male",
        compare       = "e4/4 vs e3/3",
        label         = c("e3/3", "e4/4"),
        by_group      = temp_dat.M.plot$ApoE_Indicator)  # plot_dat$ApoE_Indicator[BL_traits$Sex.int==1])

  diff_plt_geno <- diff_plot_compare(
        x               = eyo, y = median,
        df              = diff_df_geno,
        bar_df          = bar_df.geno,  #diff_df_geno[ diff_df_geno$lower*diff_df_geno$upper > 0, ],
        protein_name    = data.frame(yaxis_label = yl, title_label = ttl),
        min_value       = min_value.geno, #floor(min(diff_df_geno$lower)),
        max_value       = max_value.geno, #ceiling(max(diff_df_geno$upper)),
        First_Annotate  = First_Annotate.geno, #diff_df_geno$eyo[which(diff_df_geno$lower*diff_df_geno$upper > 0)[1] ],
        Last_Annotate   = Last_Annotate.geno, #tail(diff_df_geno$eyo[ diff_df_geno$lower*diff_df_geno$upper > 0 ], 1),
        Annotate_color  = na.omit(ifelse(diff_df_geno$upper < 0, "cornflowerblue", ifelse(diff_df_geno$lower > 0, "indianred3", NA))),   #Annotate_color.geno,
        right_limit     = right_limit_EYO.geno,
        Last_Annotate_Number = Last_Annotate_Number.geno,
        CI_Percent      = "99%",
        within          = "sex",
        compare         = "e4/4 vs e3/3")


###### second row of 3 plots

  plot_dat_non <- plot_dat[ plot_dat$ApoE_Indicator == 0, ]          # e3/e3
  plot_dat_car <- plot_dat[ plot_dat$ApoE_Indicator == 1, ]          # e4/e4


  plot_dat.non.count_before_removal_outlier <- nrow(plot_dat_non)
  temp_dat.non.plot <- plot_dat_non %>%
    mutate(zRT = scale(Pep)[,1]) %>%
    filter(between(zRT, -3, +3))
  plot_dat.non.count_post_removal_outlier <- nrow(temp_dat.non.plot)
  number_removed.non <- plot_dat.non.count_before_removal_outlier - plot_dat.non.count_post_removal_outlier  # record number of removed outliers

  scatter_nonSex <- scatter_plot_single(
          x             = EYO,  y = Pep,
          df            = temp_dat.non.plot,  # plot_dat_non,
          protein_name  = data.frame(yaxis_label = yl, title_label = ttl),
          number_removed= number_removed.non,  # 0L,
          CI_Percent    = "99%",
          combined_lines= ribbon_nonSex,
          within        = "e3/e3 (non-carrier)",
          compare       = "Female vs Male",
          label         = c("Female", "Male"),
          by_group      = temp_dat.non.plot$Sex.int)  # plot_dat_non$Sex.int)
#          ribbon_cols   = c("green3", "darkorange2"))       # optional colours


  plot_dat.car.count_before_removal_outlier <- nrow(plot_dat_car)
  temp_dat.car.plot <- plot_dat_car %>%
    mutate(zRT = scale(Pep)[,1]) %>%
    filter(between(zRT, -3, +3))
  plot_dat.car.count_post_removal_outlier <- nrow(temp_dat.car.plot)
  number_removed.car <- plot_dat.car.count_before_removal_outlier - plot_dat.car.count_post_removal_outlier  # record number of removed outliers

  scatter_carSex <- scatter_plot_single(
          x             = EYO,  y = Pep,
          df            = plot_dat_car,
          protein_name  = data.frame(yaxis_label = yl, title_label = ttl),
          number_removed= number_removed.car,  # 0L,
          CI_Percent    = "99%",

# ----------------------------------------------------------------------------------------
# Return matrices to the parallel collector.
# ----------------------------------------------------------------------------------------
          combined_lines= ribbon_carSex,
          within        = "e4/e4 (carrier)",
          compare       = "Female vs Male",
          label         = c("Female", "Male"),
          by_group      = plot_dat_car$Sex.int)
#          ribbon_cols   = c("green3", "darkorange2"))       # optional colours

  diff_plt_sex <- diff_plot_compare(
          x               = eyo, y = median,

# ----------------------------------------------------------------------------------------
# Start parallel plotting/extraction.
# ----------------------------------------------------------------------------------------
          df              = diff_df_sex,
          bar_df          = bar_df.sex,  #diff_df_sex[ diff_df_sex$lower*diff_df_sex$upper > 0, ],
          protein_name    = data.frame(yaxis_label = yl, title_label = ttl),
          min_value       = min_value.sex,  #floor(min(diff_df_sex$lower)),
          max_value       = max_value.sex,  #ceiling(max(diff_df_sex$upper)),
          First_Annotate  = First_Annotate.sex,  #diff_df_sex$eyo[ which(diff_df_sex$lower*diff_df_sex$upper > 0)[1] ],
          Last_Annotate   = Last_Annotate.sex,  #tail(diff_df_sex$eyo[ diff_df_sex$lower*diff_df_sex$upper > 0 ], 1),
          Annotate_color  = na.omit(ifelse(diff_df_sex$upper < 0, "#6b8e23", ifelse(diff_df_sex$lower > 0, "#ff8c00", NA))),   # orange if M>F, olive if F>M  # prev. single (first interval) color:  Annotate_color.sex,
          right_limit     = right_limit_EYO.sex,
          Last_Annotate_Number = Last_Annotate_Number.sex,
          CI_Percent      = "99%",
          within          = "diagnosis group",
          compare         = "Male vs Female")


  out_dir <- file.path("./sexInt.3177", "scatter")  # "z:", "ShijiaBian", "PlasmaProteomic", "Result", "20250727", "sexInt.3177", "scatter")
#                       "serialTest")        # <- choose any folder you like
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

#  ggplot2::ggsave(file.path(out_dir,
#                  sprintf("%s_scatter_Female.pdf", clean_name)),
#                  scatter_plt_F, device = cairo_pdf,
#                  width = 12, height = 12, units = "in", dpi = 300)
#
#  ggplot2::ggsave(file.path(out_dir,
#                  sprintf("%s_scatter_Male.pdf", clean_name)),
#                  scatter_plt_M, device = cairo_pdf,
#                  width = 12, height = 12, units = "in", dpi = 300)
#
#  ggplot2::ggsave(file.path(out_dir,
#                  sprintf("%s_diff.pdf", clean_name)),
#                  diff_plt_geno, device = cairo_pdf,

# ----------------------------------------------------------------------------------------
# Assemble output matrices from successful workers.
# ----------------------------------------------------------------------------------------
#                  width = 12, height = 12, units = "in", dpi = 300)


## arrange all 6 plots into one page, 2 rows, 3 columns and output to PDF
  library(gridExtra);  library(ggplotify)

  all_plts <- arrangeGrob(
                scatter_plt_F, scatter_plt_M, diff_plt_geno,   # first row
                scatter_nonSex, scatter_carSex, diff_plt_sex,  # second row
                ncol = 3)

  ggplot2::ggsave(
    filename = file.path(out_dir,
                sprintf("%s_all6plots.pdf", clean_name)),
    plot     = as_ggplot(all_plts),
    device   = cairo_pdf,
    width    = 37.5, height = 25, units = "in", dpi = 300)

  ## ---------------------------------------------------------------------
  ##  8)  everything a foreach() collector needs
  ## ---------------------------------------------------------------------
  list(t_stat.geno  = apply(diff_draws$Car_vs_Non_F, 2, function(x) t.test(x)$statistic),  #Car_vs_Non_All if averaging diff ribbons of 2 plots; otherwise _F (left plot)
       p_value.geno = apply(diff_draws$Car_vs_Non_F, 2, function(x) {
                            p <- mean(x < 0);  pmin(p, 1 - p) }),
       up_down.geno = ifelse(diff_df_geno$upper < 0, "cornflowerblue",
                        ifelse(diff_df_geno$lower > 0, "indianred3", NA)),

# ----------------------------------------------------------------------------------------
# Write posterior p-value and direction outputs.
# ----------------------------------------------------------------------------------------
       t_stat.sex  = apply(diff_draws$M_vs_F_Non, 2, function(x) t.test(x)$statistic),     #M_vs_F_All if averaging diff ribbons of 2 plots; otherwise _Non (left plot)
       p_value.sex = apply(diff_draws$M_vs_F_Non, 2, function(x) {
                           p <- mean(x < 0);  pmin(p, 1 - p) }),
       up_down.sex = ifelse(diff_df_sex$upper < 0, "#6b8e23",  #"cornflowerblue",
                        ifelse(diff_df_sex$lower > 0, "#ff8c00", NA)),  #"indianred3"
       ok      = TRUE)  # handy flag for foreach result binding
}


library(doParallel)
ncore <- max(1, parallel::detectCores() - 1)
cl     <- makeCluster(ncore)
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

genotypeContrasts_t_stats_all_pep    <- do.call(cbind,
                                                 lapply(results[ok_idx], `[[`, "t_stat.geno"))
sexContrasts_t_stats_all_pep         <- do.call(cbind,
                                                 lapply(results[ok_idx], `[[`, "t_stat.sex"))
genotype_p_value_all_pep             <- do.call(cbind,
                                                 lapply(results[ok_idx], `[[`, "p_value.geno"))
sex_p_value_all_pep                  <- do.call(cbind,
                                                 lapply(results[ok_idx], `[[`, "p_value.sex"))
up_down.geno_notation                <- do.call(rbind,
                                                 lapply(results[ok_idx], `[[`, "up_down.geno")) #[ , -1]
up_down.sex_notation                 <- do.call(rbind,
                                                 lapply(results[ok_idx], `[[`, "up_down.sex")) #[ , -1]

## Set dimnames of final data
colnames(genotypeContrasts_t_stats_all_pep)<-colnames(genotype_p_value_all_pep)<-colnames(sexContrasts_t_stats_all_pep)<-colnames(sex_p_value_all_pep)<-rownames(up_down.geno_notation)<-rownames(up_down.sex_notation)<-pep_names[ok_idx]

eyo_step = seq(floor(min(BL_traits_pep$EYO, na.rm = T)),
				   ceiling(max(BL_traits_pep$EYO, na.rm = T)),
				   by=.5)
length(eyo_step)==ncol(up_down.geno_notation)  # TRUE
rownames(genotypeContrasts_t_stats_all_pep)<-rownames(sexContrasts_t_stats_all_pep)<-rownames(genotype_p_value_all_pep)<-rownames(sex_p_value_all_pep)<-colnames(up_down.geno_notation)<-colnames(up_down.sex_notation)<-as.character(eyo_step)


## Get characteristics of waterfall
apply(up_down.geno_notation,2,function(x) table(x))
sum(apply(up_down.geno_notation,1,function(x) length(which(!is.na(x))))>0)
#880

apply(up_down.sex_notation,2,function(x) table(x))
sum(apply(up_down.sex_notation,1,function(x) length(which(!is.na(x))))>0)
#3289


# ######### ---- Write Final Outputs for Waterfall(next) step ---- ########
genotypeContrasts_t_stats_all_pep_final <- genotypeContrasts_t_stats_all_pep #[, -1]
#fil_pep <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files", paste(paste("99_par_diff_all_peptide", ".rds", sep = "")))
fil_pep<-"./sexInt.3177/scatter/_99_par_genotypeContrasts_T_diff_all_peptide.rds"
saveRDS(genotypeContrasts_t_stats_all_pep_final, fil_pep)
#fil_csv_pep <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files",paste(paste("99_par_diff_all_peptide", ".csv", sep = "")))
fil_csv_pep <- "./sexInt.3177/scatter/_99_par_genotypeContrasts_T_diff_all_peptide.csv"
write.csv(genotypeContrasts_t_stats_all_pep_final, fil_csv_pep, row.names = TRUE)

sexContrasts_t_stats_all_pep_final <- sexContrasts_t_stats_all_pep #[, -1]
#fil_pep <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files", paste(paste("99_par_diff_all_peptide", ".rds", sep = "")))
fil_pep<-"./sexInt.3177/scatter/_99_par_sexContrasts_T_diff_all_peptide.rds"
saveRDS(sexContrasts_t_stats_all_pep_final, fil_pep)
#fil_csv_pep <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files",paste(paste("99_par_diff_all_peptide", ".csv", sep = "")))
fil_csv_pep <- "./sexInt.3177/scatter/_99_par_sexContrasts_T_diff_all_peptide.csv"
write.csv(sexContrasts_t_stats_all_pep_final, fil_csv_pep, row.names = TRUE)


genotype_p_value_all_pep_final <- genotype_p_value_all_pep #[, -1]
#fil_pep <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files",paste(paste("99_par_diff_all_peptide_p_value", ".rds", sep = "")))
fil_pep<-"./sexInt.3177/scatter/_99_par_genotype_diff_all_peptide_p_value.rds"
saveRDS(genotype_p_value_all_pep_final, fil_pep)
#fil_csv_pep <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files",paste(paste("99_par_diff_all_peptide_p_value", ".csv", sep = "")))
fil_csv_pep <- "./sexInt.3177/scatter/_99_par_genotype_diff_all_peptide_p_value.csv"
write.csv(genotype_p_value_all_pep_final, fil_csv_pep, row.names = TRUE)

sex_p_value_all_pep_final <- sex_p_value_all_pep #[, -1]
#fil_pep <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files",paste(paste("99_par_diff_all_peptide_p_value", ".rds", sep = "")))
fil_pep<-"./sexInt.3177/scatter/_99_par_sex_diff_all_peptide_p_value.rds"
saveRDS(sex_p_value_all_pep_final, fil_pep)
#fil_csv_pep <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files",paste(paste("99_par_diff_all_peptide_p_value", ".csv", sep = "")))
fil_csv_pep <- "./sexInt.3177/scatter/_99_par_sex_diff_all_peptide_p_value.csv"
write.csv(sex_p_value_all_pep_final, fil_csv_pep, row.names = TRUE)


#fil_csv_up_down_notation <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files",paste(paste("99_par_diff_all_peptide_up_down_notation", ".csv", sep = "")))
fil_csv_up_down.geno_notation<-"./sexInt.3177/scatter/_99_par_genotype_diff_all_peptide_up_down_notation.csv"
write.csv(up_down.geno_notation, fil_csv_up_down.geno_notation, row.names = TRUE)

#fil_csv_up_down_notation <- file.path("Result", "20240122", "STAN_Output", "1_Ori_Model", "generated_plot_files",paste(paste("99_par_diff_all_peptide_up_down_notation", ".csv", sep = "")))
fil_csv_up_down.sex_notation<-"./sexInt.3177/scatter/_99_par_sex_diff_all_peptide_up_down_notation.csv"
write.csv(up_down.sex_notation, fil_csv_up_down.sex_notation, row.names = TRUE)


#sink("Result/20240122/Log_File/Original_STAN_Model_Scatter_Plot.txt", append = TRUE)
#print(Sys.time())
#sink()


## Reprocess updated function -- only 12 MEs
# for (pep_name in pep_names[7334:7345]) process_one_peptide(pep_name)
