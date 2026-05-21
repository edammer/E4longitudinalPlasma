################################################################################################
# Section 5 longitudinal STAN pipeline - annotated copy
# Source file: 5b.model_1_STAN_parallel.+clean(4199 samples-with non-control e3 homozygotes).R
# Pipeline process 4/12: Model 1 fitting in the 4,199-sample e3-inclusive sensitivity cohort
#
# Purpose: Fit the same APOE genotype-by-EYO Bayesian spline model in a sensitivity cohort that
# retains additional non-control APOE e3/e3 samples.
#
# Input:  APOE_homozygote_2ndFinal_cleanDats+traits(3177.CT33only_4199.33inADcontinuum).RData
# Input:  cleanDat.4199, MEs.4199, numericMeta.4199
# Output: simple.4199/*_stan_glm.rds
# Output: simple.4199/name_match_table.RDS
# Output: _numericMeta_4199_trait.RDS
# Output: _full_4199_protein_dft.RDS
#
# Major analysis steps in this script:
#   1. Build the 4,199-sample assay/module/outcome matrix.
#   2. Apply the same 5-SD outlier and <=20% missingness filters.
#   3. Define ApoE_Indicator and EYO on the 4,199-sample metadata.
#   4. Fit Model 1 across outcomes with the same rstanarm settings and parallel framework.
#   5. Save compact model RDS files and the cleaned data products.
#
# Cleanup/annotation notes:
#   - This is a cleaned, commented copy of the uploaded script; analysis logic and
#     parameter values were not intentionally changed.
#   - Files were decoded from the uploaded Windows/CP1252 text and written as UTF-8.
#   - No explicit "not run below here" block was detected in this prefix-5 file set.
################################################################################################


# ----------------------------------------------------------------------------------------
# Set the working directory and load packages for the 4,199-sample Model 1 sensitivity
# run.
# ----------------------------------------------------------------------------------------
setwd("z:/ShijiaBian/PlasmaProteomic/Result/20250727/simple.4199/")
library(tidyverse)
library(forcats)
library(rstanarm)
library(Hmisc)
library(openxlsx)
library(rstan)
#numericMeta_3177 trait <- readRDS("~/files/EBD/Shijia_B_Derived_Data/20250709/numericMeta_3177_trait.RDS")
#full_3177_protein_df <- readRos(·~/files/EBD/Shijia_B_Derived_Data/20250709/full_3177_protein_dft.RDS")

# ----------------------------------------------------------------------------------------
# Load the APOE homozygote workspace containing the 4,199-sample cohort.
# ----------------------------------------------------------------------------------------
load("z:/EBD/APOE_homozygote_2ndFinal_cleanDats+traits(3177.CT33only_4199.33inADcontinuum).RData")
ls()
#"cleanDat.3177. "cleanDat.4199"  "MEs.3177" "MEs.4199" "numericMeta.3177" "numericMeta.4199"

# ----------------------------------------------------------------------------------------
# Build the 4,199-sample outcome matrix from proteins and module eigengenes.
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
# Define the e4/e4 binary indicator in the 4,199-sample cohort.
# ----------------------------------------------------------------------------------------
numericMeta.4199$ApoE_Indicator<-0
numericMeta.4199$ApoE_Indicator[numericMeta.4199$APOE.mapped.predicted=="e4/e4"]<-1
table(numericMeta.4199$ApoE_Indicator)
#   0   1
#3786 413

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
# Utility: sanitize output names.
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
# Move sample_id to the first column for joins.
# ----------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------
# Define outcomes to model.
# ----------------------------------------------------------------------------------------
pep_names <- names(protein_df)[2:dim(protein_df)[2] ]


# ----------------------------------------------------------------------------------------
# Worker function for per-outcome Model 1 fitting.
# ----------------------------------------------------------------------------------------


one_pepSTAN <- function(track, traits_df, protein_df) {
  pep_name=pep_names[track]
  ##### --------------------- Format the data frame -------------------- #####
  dat <- traits_df %>%
    select(sample_id, ApoE_Indicator, EY0) %>%
    left_join(
      protein_df %>% select(sample_id, all_of(pep_name)),
      by = "sample_id") %>%
    filter(complete.cases(.)) %>%
    select(-contains("sample_id")) # %>%
#    rename(., Sex = Sex.int)


# ----------------------------------------------------------------------------------------
# Build the EYO spline basis.
# ----------------------------------------------------------------------------------------
  dat <- sanitize_names(dat)

  splinefit = rcspline.eval(dat$EY0, nk=3, norm = 2, pc = FALSE, inclx=TRUE)
  cubic_spline_X <- as.data.frame(splinefit)
  names(cubic_spline_X) <- c("EYO_Spline_Linear", "EYO_Spline_Cubic")
  #head(cubic_spline_X)

  dat <- cbind(dat, cubic_spline_X)

  # Model 1: Construct the formula
  outcome <- names(dat)[3]

# ----------------------------------------------------------------------------------------
# Model 1 formula mirrors the primary 3,177-sample analysis.
# ----------------------------------------------------------------------------------------
  variables1 <- c("EYO_Spline_Linear", "EYO_Spline_Cubic",
                  "ApoE_Indicator", "EYO_Spline_Linear*ApoE_Indicator", "EYO_Spline_Cubic*ApoE_Indicator")
  f1 = as.formula(paste(outcome, paste(variables1, collapse = "+"), sep = "~"), env=baseenv())
  # * environment control needed for within-function stan_glm() call, per https://discourse.mc-stan.org/t/stanfit-object-fit-inside-function-explodes-in-size-when-saved-to-rds/13656
  set.seed(track)
  env <- new.env() #parent = .GlobalEnv)
  env$dat <- dat
  env$f1 <- f1
  env$track <- track

# ----------------------------------------------------------------------------------------
# Fit Model 1 using rstanarm.
# ----------------------------------------------------------------------------------------
  stan_BL_1 <- with(env, {stan_glm(f1,
                                 data= dat,
                                 family=gaussian(),
                                 chains = 8,
                                 cores = 4,
                                 iter = 10000,
                                 thin = 10,
                                 refresh = 0,       # silence sampler
                                 seed = track) })  #, warning=function(w) {print(paste("Warning: ", track))})})
  fil_stan_1 <- file.path("z:/","ShijiaBian","PlasmaProteomic","Result","20250727", "simple.4199", paste(paste(outcome, "_stan_glm", ".rds", sep = "")))

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

  saveRDS(stan_BL_1, fil_stan_1)

  #name_match_table_temp <-
  data.frame(OriginalName = pep_name,
             CleanedName = colnames(dat)[3],
             stringsAsFactors=FALSE)
}

library(doParallel)
ncore <- 8  #max(1, parallel::detectCores() - 1)

# ----------------------------------------------------------------------------------------
# Start the parallel backend.
# ----------------------------------------------------------------------------------------
cl <- makeCluster(ncore)
registerDoParallel(cl)

library(foreach)
worker_pkgs <- c("tidyverse","rstanarm","Hmisc","rstan","dplyr")


# ----------------------------------------------------------------------------------------
# Run Model 1 across outcomes.
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
# Save model-name and cleaned sensitivity-cohort data products.
# ----------------------------------------------------------------------------------------
saveRDS(name_match_table,"z:/ShijiaBian/PlasmaProteomic/Result/20250727/simple.4199/name_match_table.RDS")
saveRDS(numericMeta_4199_trait,"_numericMeta_4199_trait.RDS")
saveRDS(protein_df,"_full_4199_protein_dft.RDS")
