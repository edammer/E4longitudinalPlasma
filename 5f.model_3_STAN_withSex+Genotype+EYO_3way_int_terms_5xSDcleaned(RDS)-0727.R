################################################################################################
# Section 5 longitudinal STAN pipeline - annotated copy
# Source file: 5f.model_3_STAN_withSex+Genotype+EYO_3way_int_terms_5xSDcleaned(RDS)-0727.R
# Pipeline process 10/12: Model 2.3 fitting for genotype-by-sex-by-EYO interaction terms
#
# Purpose: Fit the final three-way interaction model that decomposes genotype-by-sex effects
# into offset, linear-EYO, and cubic-EYO components.
#
# Input:  _numericMeta_3177_trait.RDS
# Input:  _full_3177_protein_dft.RDS
# Output: /home/labshare/genoSex/*_with_APOE.Sex_GenotypeXSex_stan_glm.rds
# Output: sex.APOE.intTerm/name_match_table.RDS
#
# Major analysis steps in this script:
#   1. Load the primary cleaned data and retain imputed genotype calls.
#   2. Merge each outcome with ApoE_Indicator, EYO, and Sex.int.
#   3. Construct restricted cubic spline terms for EYO.
#   4. Fit Model 2.3 with genotype, sex, genotype-by-EYO, sex-by-EYO, genotype-by-sex, and
#      genotype-by-sex-by-EYO terms.
#   5. Save compact model fits for post-hoc joint posterior testing.
#
# Cleanup/annotation notes:
#   - This is a cleaned, commented copy of the uploaded script; analysis logic and
#     parameter values were not intentionally changed.
#   - Files were decoded from the uploaded Windows/CP1252 text and written as UTF-8.
#   - No explicit "not run below here" block was detected in this prefix-5 file set.
################################################################################################


# ----------------------------------------------------------------------------------------
# Set the working directory and load packages for final Model 2.3 fitting.
# ----------------------------------------------------------------------------------------
#setwd("z:/ShijiaBian/PlasmaProteomic/Result/20250727/simple.3177/")
setwd("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/")

library(tidyverse)
library(forcats)
library(rstanarm)
library(Hmisc)
library(openxlsx)
library(rstan)


# ----------------------------------------------------------------------------------------
# Load preprocessed traits and outcome matrix; imputed genotype samples are retained.
# ----------------------------------------------------------------------------------------
# 1. Load pre-processed RDS files
numericMeta_3177_trait <- readRDS("_numericMeta_3177_trait.RDS")
protein_df             <- readRDS("_full_3177_protein_dft.RDS")
# protein_df already contains MMSE, cdr, MEs, and sample_id as first column

## DO NOT REMOVE IMPUTED GENOTYPE SAMPLES
## 2. Remove samples with missing APOE.mapped
#na_apoe_idx <- which(is.na(numericMeta_3177_trait$APOE.mapped))
#message(sprintf("Removing %d sample(s) with NA APOE.mapped.", length(na_apoe_idx)))
##Removing 708 sample(s) with NA APOE.mapped.
#
#if (length(na_apoe_idx) > 0) {
#  na_sample_ids          <- rownames(numericMeta_3177_trait)[na_apoe_idx]
#  numericMeta_3177_trait <- numericMeta_3177_trait[-na_apoe_idx, ]
#  protein_df             <- protein_df[!protein_df$sample_id %in% na_sample_ids, ]
#}


# ----------------------------------------------------------------------------------------
# Verify metadata/outcome alignment.
# ----------------------------------------------------------------------------------------
# 3. Verify alignment between metadata and protein matrix
stopifnot(all(protein_df$sample_id == numericMeta_3177_trait$sample_id))


# ----------------------------------------------------------------------------------------
# Define ApoE_Indicator and EYO.
# ----------------------------------------------------------------------------------------
# 4. Set derived traits
numericMeta_3177_trait$ApoE_Indicator <- 0
numericMeta_3177_trait$ApoE_Indicator[numericMeta_3177_trait$APOE.mapped.predicted == "e4/e4"] <- 1
table(numericMeta_3177_trait$ApoE_Indicator)

numericMeta_3177_trait$EY0 <- (65.6 - numericMeta_3177_trait$age_at_visit) * (-1)
range(numericMeta_3177_trait$EY0)

########################################################################
#' @param x A data frame (or any object with names)
#' @return The same object, but with names cleaned

# ----------------------------------------------------------------------------------------
# Utility: sanitize outcome names.
# ----------------------------------------------------------------------------------------
sanitize_names <- function(x) {
  nm <- names(x)
  if (is.null(nm)) return(x)
  clean_nm <- gsub("[^[:alnum:]]", "_", nm)
  clean_nm <- gsub("_+", "_", clean_nm)
  clean_nm <- gsub("^_|_$", "", clean_nm)
  names(x) <- clean_nm
  x
}
########################################################################


# ----------------------------------------------------------------------------------------
# Define outcomes for Model 2.3.
# ----------------------------------------------------------------------------------------
pep_names <- names(protein_df)[2:dim(protein_df)[2]]


# ----------------------------------------------------------------------------------------
# Worker function: join traits and outcome for one Model 2.3 fit.
# ----------------------------------------------------------------------------------------
one_pepSTAN <- function(track, traits_df, protein_df) {
  pep_name=pep_names[track]
  ##### --------------------- Format the data frame -------------------- #####
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
# Generate EYO spline basis.
# ----------------------------------------------------------------------------------------
  splinefit = rcspline.eval(dat$EY0, nk=3, norm = 2, pc = FALSE, inclx=TRUE)
  cubic_spline_X <- as.data.frame(splinefit)
  names(cubic_spline_X) <- c("EYO_Spline_Linear", "EYO_Spline_Cubic")
  #head(cubic_spline_X)

  dat <- cbind(dat, cubic_spline_X)
  outcome <- names(dat)[4]


# ----------------------------------------------------------------------------------------
# Commented Model 2.1/2.2 blocks document earlier model definitions but are not executed
# in this script.
# ----------------------------------------------------------------------------------------
#  # Model 2.1: APOE spline interactions + Sex main effect
#  # Model 2.1: Construct the formula, add Sex to the first model
#  variables1 <- c("EYO_Spline_Linear", "EYO_Spline_Cubic",
#                  "ApoE_Indicator", "EYO_Spline_Linear*ApoE_Indicator", "EYO_Spline_Cubic*ApoE_Indicator",
#                  "Sex")
#  f1 = as.formula(paste(outcome, paste(variables1, collapse = "+"), sep = "~"), env=baseenv())
#  # * environment control needed for within-function stan_glm() call, per https://discourse.mc-stan.org/t/stanfit-object-fit-inside-function-explodes-in-size-when-saved-to-rds/13656
#  set.seed(track)
#  env <- new.env() #parent = .GlobalEnv)
#  env$dat <- dat
#  env$f1 <- f1
#  env$track <- track
#  stan_BL_1 <- with(env, {stan_glm(f1,
#                                 data= dat,
#                                 family=gaussian(),
#                                 chains = 8,
#                                 cores = 4,
#                                 iter = 10000,
#                                 thin = 10,
#                                 refresh = 0,       # silence sampler
#                                 seed = track) })  #, warning=function(w) {print(paste("Warning: ", track))})})
#  fil_stan_1 <- file.path("f:/", "OneDrive - Emory", "Legacy", "e4_homozygoteStudy",
#                          "DL", "sex.APOE.intTerm",
#                          paste(paste(outcome, "_with_Sex_stan_glm", ".rds", sep = "")))
#
#  # Model 2.2: + Sex spline interactions
#  # Model 2.2: Construct the formula, add Sex to the first model
#  variables2 <- c("EYO_Spline_Linear", "EYO_Spline_Cubic",
#                  "ApoE_Indicator", "EYO_Spline_Linear*ApoE_Indicator", "EYO_Spline_Cubic*ApoE_Indicator",
#                  "Sex", "EYO_Spline_Linear*Sex", "EYO_Spline_Cubic*Sex")
#  f2 = as.formula(paste(outcome, paste(variables2, collapse="+"), sep = "~"), env=baseenv())
#  set.seed(track+length(pep_names))
#  env$f2 <- f2
#  env$pep_names <- pep_names
#  stan_BL_2 <- with(env, {tryCatch(stan_glm(f2,
#                                 data = dat,
#                                 family=gaussian(),
#                                 chains = 8,
#                                 cores = 4,
#                                 iter = 10000,
#                                 thin = 10,
#                                 refresh = 0,
#                                 seed = track + length(pep_names)), warning=function(w) {print(paste("Warning: ", track))})})
#
#  fil_stan_2 <- file.path("f:/", "OneDrive - Emory", "Legacy", "e4_homozygoteStudy",
#                          "DL", "sex.APOE.intTerm",
#                          paste(paste(outcome, "_with_APOE.Sex_Interaction_stan_glm", ".rds", sep = "")))


# ----------------------------------------------------------------------------------------
# Model 2.3 formula includes genotype-by-sex and genotype-by-sex-by-EYO spline terms.
# ----------------------------------------------------------------------------------------
  # Model 2.3: + ApoE_Indicator*Sex interaction term
  variables3 <- c(
    "EYO_Spline_Linear", "EYO_Spline_Cubic",
    "ApoE_Indicator",
    "EYO_Spline_Linear*ApoE_Indicator", "EYO_Spline_Cubic*ApoE_Indicator",
    "Sex",
    "EYO_Spline_Linear*Sex", "EYO_Spline_Cubic*Sex",
    "ApoE_Indicator:Sex",                        # time-invariant genotype × sex offset
    "EYO_Spline_Linear:ApoE_Indicator:Sex",      # linear EYO component of genotype × sex
    "EYO_Spline_Cubic:ApoE_Indicator:Sex"        # cubic EYO component of genotype × sex
  )
  f3      <- as.formula(paste(outcome, paste(variables3, collapse = "+"), sep = "~"), env = baseenv())
  env <- new.env() #parent = .GlobalEnv)
  env$dat <- dat
  env$f3  <- f3
  set.seed(track + 2 * length(pep_names))

# ----------------------------------------------------------------------------------------
# Fit the Model 2.3 Bayesian Gaussian model.
# ----------------------------------------------------------------------------------------
  stan_BL_3 <- with(env, {
    tryCatch(
      stan_glm(f3, data = dat, family = gaussian(),
               chains = 8, cores = 4, iter = 10000, thin = 10, refresh = 0,
               seed = track + 2 * length(pep_names)),
      warning = function(w) { print(paste("Warning (model 2.3):", track)) })
  })
#  fil_stan_3 <- file.path("f:/", "OneDrive - Emory", "Legacy", "e4_homozygoteStudy",
#                          "DL", "sex.APOE.intTerm",
  fil_stan_3 <- file.path("/home/labshare/genoSex",  # path on telomere
                          paste0(outcome, "_with_APOE.Sex_GenotypeXSex_stan_glm.rds"))


# ----------------------------------------------------------------------------------------
# Strip formula environments.
# ----------------------------------------------------------------------------------------
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


  # Strip environments and save
#  stan_BL_1 <- strip_formula_envs(stan_BL_1)
#  stan_BL_2 <- strip_formula_envs(stan_BL_2)
  stan_BL_3 <- strip_formula_envs(stan_BL_3)

#  saveRDS(stan_BL_1, fil_stan_1)
#  saveRDS(stan_BL_2, fil_stan_2)

# ----------------------------------------------------------------------------------------
# Save the Model 2.3 fit for post-hoc tests.
# ----------------------------------------------------------------------------------------
  saveRDS(stan_BL_3, fil_stan_3)

  #name_match_table_temp <-
  data.frame(OriginalName = pep_name,
             CleanedName  = outcome,
#             CleanedName = colnames(dat)[4],
             stringsAsFactors = FALSE)
}


# ----------------------------------------------------------------------------------------
# Start the parallel backend.
# ----------------------------------------------------------------------------------------
library(doParallel)
library(foreach)

ncore <- 10
cl    <- makeCluster(ncore)
registerDoParallel(cl)

worker_pkgs <- c("tidyverse", "rstanarm", "Hmisc", "rstan", "dplyr")


# ----------------------------------------------------------------------------------------
# Run Model 2.3 fitting across outcomes.
# ----------------------------------------------------------------------------------------
results <- foreach(track          = 1:length(pep_names),
                   .packages      = worker_pkgs,
                   .export        = c("one_pepSTAN", "sanitize_names", "pep_names"),
                   .errorhandling = "pass") %dopar% {
  one_pepSTAN(track, numericMeta_3177_trait, protein_df)
}

stopImplicitCluster()


# ----------------------------------------------------------------------------------------
# Save the name-match table for Model 2.3 outputs.
# ----------------------------------------------------------------------------------------
name_match_table         <- do.call(rbind, results)
colnames(name_match_table) <- c("OriginalName", "CleanedName")

saveRDS(name_match_table,
        "F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/sex.APOE.intTerm/name_match_table.RDS")
#        "/home/labshare/genoSex/name_match_table.RDS")  # on https://telomere.biochem.emory.edu/R/
