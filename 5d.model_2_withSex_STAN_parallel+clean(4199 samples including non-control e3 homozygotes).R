################################################################################################
# Section 5 longitudinal STAN pipeline - annotated copy
# Source file: 5d.model_2_withSex_STAN_parallel+clean(4199 samples including non-control e3 homozygotes).R
# Pipeline process 7/12: Models 2.1 and 2.2 fitting in the 4,199-sample sensitivity cohort
#
# Purpose: Repeat sex-adjusted and sex-by-EYO Bayesian spline modeling in the e3-inclusive
# 4,199-sample cohort.
#
# Input:  APOE_homozygote_2ndFinal_cleanDats+traits(3177.CT33only_4199.33inADcontinuum).RData
# Input:  cleanDat.4199, MEs.4199, numericMeta.4199
# Output: sexInt.4199/*_with_Sex_stan_glm.rds
# Output: sexInt.4199/*_with_Sex_Interaction_stan_glm.rds
# Output: sexInt.4199/name_match_table.RDS
#
# Major analysis steps in this script:
#   1. Build the 4,199-sample cleaned outcome matrix.
#   2. Merge metadata containing ApoE_Indicator, EYO, and Sex.int.
#   3. Fit Model 2.1 and Model 2.2 with the same sampler settings as the primary analysis.
#   4. Save model outputs for sensitivity analysis of sex-EYO effects.
#
# Cleanup/annotation notes:
#   - This is a cleaned, commented copy of the uploaded script; analysis logic and
#     parameter values were not intentionally changed.
#   - Files were decoded from the uploaded Windows/CP1252 text and written as UTF-8.
#   - No explicit "not run below here" block was detected in this prefix-5 file set.
################################################################################################


# ----------------------------------------------------------------------------------------
# Set the working directory and load packages for 4,199-sample sex-interaction modeling.
# ----------------------------------------------------------------------------------------
setwd("z:/ShijiaBian/PlasmaProteomic/Result/20250727/sexInt.4199/")
library(tidyverse)
library(forcats)
library(rstanarm)
library(Hmisc)
library(openxlsx)
library(rstan)
#numericMeta_3177 trait <- readRDS("~/files/EBD/Shijia_B_Derived_Data/20250709/numericMeta_3177_trait.RDS")
#full_3177_protein_df <- readRos(·~/files/EBD/Shijia_B_Derived_Data/20250709/full_3177_protein_dft.RDS")

# ----------------------------------------------------------------------------------------
# Load the APOE homozygote workspace.
# ----------------------------------------------------------------------------------------
load("z:/EBD/APOE_homozygote_2ndFinal_cleanDats+traits(3177.CT33only_4199.33inADcontinuum).RData")
ls()
#"cleanDat.3177. "cleanDat.4199"  "MEs.3177" "MEs.4199" "numericMeta.3177" "numericMeta.4199"

# ----------------------------------------------------------------------------------------
# Build the 4,199-sample outcome matrix.
# ----------------------------------------------------------------------------------------
full_4199_protein_df <- as.data.frame(rbind(cleanDat.4199, t(MEs.4199)))


# ----------------------------------------------------------------------------------------
# Apply 5-SD outlier removal.
# ----------------------------------------------------------------------------------------
## Data cleaning - 5SD from mean max within protein, then max 20% per genotype group NA; add MMSE and cdr
clean_mat <- as.data.frame(t(apply(full_4199_protein_df, 1, function(row) {
    ## winsorise at +/- 5 SD
  z   <- abs(row - mean(row, na.rm = TRUE))
  row[z > 5 * sd(row, na.rm = TRUE)] <- NA                 # outliers -> NA

  row }                                              # otherwise return the cleaned row
)))
## restore dimnames -------------------------------------------------------
rownames(clean_mat) <- rownames(full_4199_protein_df)
colnames(clean_mat) <- colnames(full_4199_protein_df)


# ----------------------------------------------------------------------------------------
# Index samples by APOE genotype for missingness checks.
# ----------------------------------------------------------------------------------------
## grouping factor – one value per column in the expression matrix
grp <- numericMeta.4199$APOE.mapped.predicted
stopifnot(length(grp) == ncol(full_4199_protein_df))   # safety

## pre-compute the column indices that belong to each group
idx_by_grp <- split(seq_along(grp), grp)

## rows that fail the “<=20% NA in every group” rule  ----------------------
bad_row_idx <- which(
  apply(clean_mat, 1, function(v) {
    any(                                           # if *any* group exceeds 20% NA
      vapply(idx_by_grp,
             function(ix) mean(is.na(v[ix])) > 0.20,
             logical(1))
    )
  })
)

bad_row_idx        # numeric vector of offending row indices                                    # -> row is discarded
#named integer(0)

## If any found, remove by %in% rownames(clean_mat)[bad_row_idx]
if(length(bad_row_idx)>0) clean_mat<-clean_mat[which(!rownames(clean_mat) %in% rownames(clean_mat)[bad_row_idx]),]


## Add MMSE, cdr as additional outcomes
full_4199_protein_df<-as.data.frame(t(clean_mat))
full_4199_protein_df$MMSE<-numericMeta.4199$MMSE
full_4199_protein_df$cdr<-numericMeta.4199$cdr

full_4199_protein_df$sample_id <- rownames(full_4199_protein_df)


# set derived traits

# ----------------------------------------------------------------------------------------
# Define the e4/e4 binary indicator.
# ----------------------------------------------------------------------------------------
numericMeta.4199$ApoE_Indicator<-0
numericMeta.4199$ApoE_Indicator[numericMeta.4199$APOE.mapped.predicted=="e4/e4"]<-1
table(numericMeta.4199$ApoE_Indicator)
#   0   1
#2764 413

# ----------------------------------------------------------------------------------------
# Compute EYO.
# ----------------------------------------------------------------------------------------
numericMeta.4199$EY0 <- (65.6-numericMeta.4199$age_at_visit)*(-1)
range(numericMeta.4199$EY0)
# -45.6 24.4
numericMeta_4199_trait <- numericMeta.4199
numericMeta_4199_trait$sample_id <- rownames(numericMeta_4199_trait)

which(!full_4199_protein_df$sample_id==numericMeta_4199_trait$sample_id)
# integer(0)

########################################################################
#' @param x A data frame (or any object with names)
#' @return The same object, but with names cleaned
sanitize_names <- function(x) {

# ----------------------------------------------------------------------------------------
# Utility: sanitize names.
# ----------------------------------------------------------------------------------------
  #Capture existing names
  nm <- names(x)
  if (is.null(nm)) return(x)

  # Replace anything that's not a letter or digit with
  clean_nm <- gsub("[^[:alnum:]]", "_", nm)
  # (Optional) collapse repeated underscores:
  clean_nm <- gsub("_+", "_", clean_nm)
  # (Optional) trim leading/trailing underscores:
  clean_nm <- gsub("^_|_$", "", clean_nm)

  names(x) <- clean_nm
  x
}
##########################################################################

# For the protein file, make the person_id to be the first column
protein_df <- full_4199_protein_df %>%
  select(sample_id, everything())  #  Moves 'AAA' to the first column


# ----------------------------------------------------------------------------------------
# Define outcomes to model.
# ----------------------------------------------------------------------------------------
pep_names <- names(protein_df)[2:dim(protein_df)[2] ]


# ----------------------------------------------------------------------------------------
# Worker function for fitting Models 2.1 and 2.2 in the sensitivity cohort.
# ----------------------------------------------------------------------------------------
one_pepSTAN <- function(track, traits_df, protein_df) {
  pep_name=pep_names[track]
  ##### --------------------- Format the data frame -------------------- #####

# ----------------------------------------------------------------------------------------
# Join ApoE_Indicator, EYO, and Sex.int to the outcome.
# ----------------------------------------------------------------------------------------
  dat <- traits_df %>%
    select(sample_id, ApoE_Indicator, EY0, Sex.int) %>%
    left_join(
      protein_df %>% select(sample_id, all_of(pep_name)),
      by = "sample_id") %>%
    filter(complete.cases(.)) %>%
    select(-contains("sample_id")) %>%
    rename(., Sex = Sex.int)

  dat <- sanitize_names(dat)


# ----------------------------------------------------------------------------------------
# Generate EYO spline terms.
# ----------------------------------------------------------------------------------------
  splinefit = rcspline.eval(dat$EY0, nk=3, norm = 2, pc = FALSE, inclx=TRUE)
  cubic_spline_X <- as.data.frame(splinefit)
  names(cubic_spline_X) <- c("EYO_Spline_Linear", "EYO_Spline_Cubic")
  #head(cubic_spline_X)

  dat <- cbind(dat, cubic_spline_X)


# ----------------------------------------------------------------------------------------
# Model 2.1 formula.
# ----------------------------------------------------------------------------------------
  # Model 2.1: Construct the formula, add Sex to the first model
  outcome <- names(dat)[4]
  variables1 <- c("EYO_Spline_Linear", "EYO_Spline_Cubic",
                  "ApoE_Indicator", "EYO_Spline_Linear*ApoE_Indicator", "EYO_Spline_Cubic*ApoE_Indicator",
                  "Sex")
  f1 = as.formula(paste(outcome, paste(variables1, collapse = "+"), sep = "~"), env=baseenv())
  # * environment control needed for within-function stan_glm() call, per https://discourse.mc-stan.org/t/stanfit-object-fit-inside-function-explodes-in-size-when-saved-to-rds/13656
  set.seed(track)
  env <- new.env() #parent = .GlobalEnv)
  env$dat <- dat
  env$f1 <- f1
  env$track <- track
  stan_BL_1 <- with(env, {stan_glm(f1,
                                 data= dat,
                                 family=gaussian(),
                                 chains = 8,
                                 cores = 4,
                                 iter = 10000,
                                 thin = 10,
                                 refresh = 0,       # silence sampler
                                 seed = track) })  #, warning=function(w) {print(paste("Warning: ", track))})})
  fil_stan_1 <- file.path("z:/","ShijiaBian","PlasmaProteomic","Result","20250727", "sexInt.4199", paste(paste(outcome, "_with_Sex_stan_glm", ".rds", sep = "")))


# ----------------------------------------------------------------------------------------
# Model 2.2 formula with sex-by-EYO spline interactions.
# ----------------------------------------------------------------------------------------
  # Model 2.2: Construct the formula, add Sex to the first model
  variables2 <- c("EYO_Spline_Linear", "EYO_Spline_Cubic",
                  "ApoE_Indicator", "EYO_Spline_Linear*ApoE_Indicator", "EYO_Spline_Cubic*ApoE_Indicator",
                  "Sex", "EYO_Spline_Linear*Sex", "EYO_Spline_Cubic*Sex")
  f2 = as.formula(paste(outcome, paste(variables2, collapse="+"), sep = "~"), env=baseenv())
  set.seed(track+length(pep_names))
  env$f2 <- f2
  env$pep_names <- pep_names
  stan_BL_2 <- with(env, {tryCatch(stan_glm(f2,
                                 data = dat,
                                 family=gaussian(),
                                 chains = 8,
                                 cores = 4,
                                 iter = 10000,
                                 thin = 10,
                                 refresh = 0,
                                 seed = track + length(pep_names)), warning=function(w) {print(paste("Warning: ", track))})})
  fil_stan_2 <- file.path("z:/","ShijiaBian","PlasmaProteomic","Result","20250727", "sexInt.4199", paste(paste(outcome, "_with_Sex_Interaction_stan_glm", ".rds", sep = "")))

  strip_formula_envs <- function(x) {  # * environment stripping ensures small RDS file output, even when model generated from within a function (parallelized)
    rapply(x,
           # function run on every match
           f = function(obj) {
                 attr(obj, ".Environment") <- baseenv()
                 obj                       # must return the modified object
               },
           classes = c("formula", "terms"),   # what we are looking for
           how    = "replace")                # replace in situ, keep structure
  }
  stan_BL_1 <- strip_formula_envs(stan_BL_1)
  stan_BL_2 <- strip_formula_envs(stan_BL_2)

  saveRDS(stan_BL_1, fil_stan_1)
  saveRDS(stan_BL_2, fil_stan_2)

  #name_match_table_temp <-
  data.frame(OriginalName = pep_name,
             CleanedName = colnames(dat)[4],
             stringsAsFactors=FALSE)
}


# ----------------------------------------------------------------------------------------
# Start the parallel backend.
# ----------------------------------------------------------------------------------------
library(doParallel)
ncore <- 8  #max(1, parallel::detectCores() - 1)
cl <- makeCluster(ncore)
registerDoParallel(cl)

library(foreach)
worker_pkgs <- c("tidyverse","rstanarm","Hmisc","rstan","dplyr")

# ----------------------------------------------------------------------------------------
# Run model fitting across outcomes.
# ----------------------------------------------------------------------------------------

results <- foreach(track = 1:length(pep_names),
                   .packages = worker_pkgs,
                   .export = c("one_pepSTAN","sanitize_names",    # the functions
                               "pep_names"),
                   .errorhandling = "pass")  %dopar%  {
  one_pepSTAN(track, numericMeta_4199_trait, protein_df)
}

stopImplicitCluster()

name_match_table<-do.call(rbind,results)
colnames(name_match_table)<-c("OriginalName","CleanedName")

# ----------------------------------------------------------------------------------------
# Save model-name lookup table.
# ----------------------------------------------------------------------------------------

saveRDS(name_match_table,"z:/ShijiaBian/PlasmaProteomic/Result/20250727/sexInt.4199/name_match_table.RDS")
