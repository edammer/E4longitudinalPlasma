## =============================================================================
## 5e. Model 1 STAN fitting - APOE-genotype-known-only sensitivity cohort
##     (spline knots fixed to full 3177-sample/imputed-genotype positions)
## =============================================================================
##
## PURPOSE
##   Sensitivity analysis for the primary Model 1 Bayesian restricted cubic
##   spline (RCS) longitudinal model: refit every assay using only the subset
##   of the 3,177-sample cohort that has a directly mapped (non-imputed) APOE
##   genotype, while keeping the RCS knot x-positions (EYO) identical to those
##   used in the original full-cohort (imputed-genotype-included) Model 1 run.
##   This isolates the effect of genotype imputation from the effect of the
##   spline basis itself when comparing trajectories between cohorts.
##
## KEY DESIGN POINT (patched 2026-05-27, see "CHANGED" comments below)
##   rcspline.eval() chooses knot locations from the EYO values of whichever
##   data are passed to it. Because each protein has a different pattern of
##   missingness, knot positions computed only on the "known genotype" subset
##   would differ assay-by-assay from the knots used in the original
##   full-cohort model, making the two model sets harder to compare directly.
##   This version therefore computes, per assay, the 3-knot RCS positions
##   from that assay's COMPLETE-CASE EYO values in the FULL 3,177-sample data
##   (including imputed-genotype samples), then evaluates the spline basis
##   using those fixed knots on the FILTERED (known-genotype-only) data that
##   is actually used to fit the Bayesian model. The full-data objects are
##   used only to pick knot locations; they never enter the regression.
##
## STEP-BY-STEP PIPELINE
##   1. Load the full 3,177-sample trait table and protein/outcome matrix
##      (RDS), and keep an unfiltered copy of each for spline-knot estimation.
##   2. Remove samples with no mapped APOE genotype (APOE.mapped == NA) from
##      the modeling copies of the trait/protein tables.
##   3. Verify sample-ID alignment between the filtered trait and protein
##      tables.
##   4. Derive EYO ("estimated years from onset") and a binary e4/e4 indicator
##      from age at visit and predicted APOE genotype, for both the filtered
##      (modeling) and full (knot-estimation) trait tables.
##   5. Re-apply (a no-op here, retained for provenance) the +/-5 SD outlier
##      censoring step used upstream, and rebuild a parallel "for knots only"
##      protein matrix from the unfiltered full-cohort data.
##   6. Drop any assay rows that fail a 20%-missing-per-genotype-group QC rule.
##   7. For each assay (in parallel), in one_pepSTAN():
##        a. Join filtered trait + outcome data, drop incomplete rows.
##        b. Estimate RCS knot positions for THIS assay from the full-cohort
##           complete cases (get_full3177_spline_info_for_peptide()).
##        c. Evaluate the spline basis on the filtered modeling data using
##           those full-cohort knots.
##        d. Fit a Bayesian linear model (rstanarm::stan_glm) of the outcome
##           on linear/cubic EYO spline terms, e4/e4 genotype, and their
##           interactions.
##        e. Strip formula environments (keeps saved RDS objects small/
##           portable) and attach the knot metadata used as object attributes.
##        f. Save the per-assay model fit to its own RDS file.
##   8. Identify any assays whose missingness changed between this run and a
##      prior run (process_again) and, if any are found, refit them serially
##      against the unfiltered protein_df instead of the cleaned matrix.
##   9. Reconcile/refresh the OriginalName <-> CleanedName lookup table used
##      by the companion plotting script (5e2) to find each assay's RDS file.
##  10. Save filtered (modeling) and full-cohort (knot-estimation) trait and
##      protein objects together to assays+traits_forPlots.RData for 5e2.
##
## REQUIRED INPUTS
##   - _numericMeta_3177_trait.RDS   (3,177-sample trait/metadata table)
##   - _full_3177_protein_dft.RDS    (3,177-sample protein/outcome matrix;
##                                    includes MMSE, cdr, and WGCNA module
##                                    eigengene columns)
##   - ./name_match_table.RDS        (prior OriginalName/CleanedName lookup,
##                                    refreshed and re-saved at the end)
##
## MAJOR OUTPUTS
##   - One <CleanedAssayName>_stan_glm.rds Bayesian model fit per assay,
##     written to the working directory, each carrying the full-cohort
##     spline-knot attributes used to fit it.
##   - ./name_match_table.RDS (refreshed OriginalName/CleanedName lookup)
##   - ./assays+traits_forPlots.RData (filtered + full-cohort trait/protein
##     objects consumed by the companion plotting script, 5e2)
## =============================================================================

#setwd("z:/ShijiaBian/PlasmaProteomic/Result/20250727/simple.3177/")
#setwd("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/")
setwd("/home/labshare/model1.noImpute/")

library(tidyverse)
library(forcats)
library(rstanarm)
library(Hmisc)
library(openxlsx)
library(rstan)

# ---- STEP 1. Load pre-processed RDS files (full 3,177-sample trait + protein/outcome tables) ----
numericMeta_3177_trait <- readRDS("_numericMeta_3177_trait.RDS")
protein_df             <- readRDS("_full_3177_protein_dft.RDS")
# protein_df already contains MMSE, cdr, MEs, and sample_id as first column

# CHANGED 2026-05-27: Preserve the unfiltered 3177-sample traits/protein
# data before removing APOE.mapped == NA samples.  These full-data objects
# are used only to estimate restricted cubic spline knot positions so the
# knot locations reflect the full cohort, including imputed-genotype samples.
numericMeta_3177_trait_full <- numericMeta_3177_trait
protein_df_full             <- protein_df
numericMeta_3177_trait_full$EY0 <- (65.6 - numericMeta_3177_trait_full$age_at_visit) * (-1)
numericMeta_3177_trait_full$EYO <- numericMeta_3177_trait_full$EY0

# ---- STEP 2. Remove samples with no mapped (known) APOE genotype -> defines the modeling cohort for this sensitivity analysis ----
na_apoe_idx <- which(is.na(numericMeta_3177_trait$APOE.mapped))
message(sprintf("Removing %d sample(s) with NA APOE.mapped.", length(na_apoe_idx)))
#Removing 708 sample(s) with NA APOE.mapped.

if (length(na_apoe_idx) > 0) {
  na_sample_ids          <- rownames(numericMeta_3177_trait)[na_apoe_idx]
  numericMeta_3177_trait <- numericMeta_3177_trait[-na_apoe_idx, ]
  protein_df             <- protein_df[!protein_df$sample_id %in% na_sample_ids, ]
}

# ---- STEP 3. Verify sample-ID alignment between the filtered trait table and protein/outcome matrix ----
stopifnot(all(protein_df$sample_id == numericMeta_3177_trait$sample_id))

# ---- STEP 4. Derive EYO (estimated years from onset) and the binary e4/e4 genotype indicator ----
numericMeta_3177_trait$ApoE_Indicator <- 0
numericMeta_3177_trait$ApoE_Indicator[numericMeta_3177_trait$APOE.mapped.predicted == "e4/e4"] <- 1
table(numericMeta_3177_trait$ApoE_Indicator)
#   0    1 
#2106  363


numericMeta_3177_trait$EY0 <- (65.6 - numericMeta_3177_trait$age_at_visit) * (-1)
# CHANGED 2026-05-27: Keep an EYO alias in the filtered, unimputed-genotype
# modeling traits so later saved plotting objects have an explicit EYO column.
numericMeta_3177_trait$EYO <- numericMeta_3177_trait$EY0
range(numericMeta_3177_trait$EY0)
#-45.6  24.4



tail(colnames(protein_df))
#[1] "magenta" "brown"   "black"   "pink"    "MMSE"    "cdr"
# avoid 5xSD Z removal with MMSE and cdr in place


# ---- STEP 5. Re-check (no-op) +/-5 SD outlier censoring already performed
# upstream before the input RDS was saved, and build the parallel full-cohort
# matrix used only to choose spline knots. ----
## Already performed on original data ingress before saving to .RDS. Trying again here.
## Data cleaning - 5SD from mean max within protein, then max 20% per genotype group NA; add MMSE and cdr
clean_mat <- as.data.frame(t(apply(protein_df[,2:(ncol(protein_df)-2)], 1, function(row) {
#    ## commented b/c input data already checked for outliers and windsorized once. ~15 assays find additional if re-winsorise at +/- 5 SD
#  z   <- abs(row - mean(row, na.rm = TRUE))
#  row[z > 5 * sd(row, na.rm = TRUE)] <- NA                 # outliers -> NA
  
  row }                                              # otherwise return the cleaned row
)))
## restore dimnames -------------------------------------------------------
rownames(clean_mat) <- rownames(protein_df)
colnames(clean_mat) <- colnames(protein_df)[2:(ncol(protein_df)-2)]

# CHANGED 2026-05-27: Build a parallel cleaned protein matrix on the full
# 3177-sample data.  It is not used for model fitting; it supplies the
# protein-specific complete-case EYO values used to choose spline knots in
# one_pepSTAN().  The same censoring operation as clean_mat is used so
# protein-specific missingness/censoring can influence knot placement.
clean_mat_full3177_for_knots <- as.data.frame(t(apply(protein_df_full[,2:(ncol(protein_df_full)-2)], 1, function(row) {
#  z <- abs(row - mean(row, na.rm = TRUE))
#  row[z > 5 * sd(row, na.rm = TRUE)] <- NA

  row
})))
rownames(clean_mat_full3177_for_knots) <- rownames(protein_df_full)
colnames(clean_mat_full3177_for_knots) <- colnames(protein_df_full)[2:(ncol(protein_df_full)-2)]



## grouping factor - one value per column in the expression matrix
grp <- numericMeta_3177_trait$APOE.mapped.predicted
#stopifnot(length(grp) == ncol(protein_df))   # safety

## pre-compute the column indices that belong to each group
idx_by_grp <- split(seq_along(grp), grp)

# ---- STEP 6. Drop any assay (row) with >20% missing values within either
# genotype group; none are dropped in this run, but the rule is retained for
# provenance/QC consistency with upstream processing. ----
## rows that fail the "<=20% NA in every group" rule  ----------------------
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



assays.2rerun <- colnames(clean_mat)[which( !apply(clean_mat,2,function(x) length(which(is.na(x)))) == apply(protein_df[,2:(ncol(protein_df)-2)],2,function(x) length(which(is.na(x)))) )]
assays.2rerun
#character(0) 

#previously: with re-removal of 5xSD outliers ('windsorization'):
# [1] "NCF2|P19878"                       "ENO1|P06733"                       "CHGA|P10645^SL002762@seq.11184.51" "CDC25A|P30304"                    
# [5] "SNAP25|P60880"                     "MRPL33|O75394"                     "IDH1|O75874"                       "SERPINA5|P05154"                  
# [9] "CNDP1|Q96KN2^SL006694@seq.5456.59" "LEAP2|Q969E1"                      "MRPL58|Q14197"                     "CNDP1|Q96KN2^SL006694@seq.7870.8" 
#[13] "OLFM2|O95897^SL012399@seq.8295.16" "CHGA|P10645^SL002762@seq.8476.11"  "PLXDC1|Q8IUK5"



## ---- STEP 7. Helper functions used by the per-assay model loop below: ----
##   sanitize_names()                      -> makes outcome names valid R names
##   get_full3177_spline_info_for_peptide() -> per-assay full-cohort knot lookup
########################################################################
#' @param x A data frame (or any object with names)
#' @return The same object, but with names cleaned
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

# CHANGED 2026-05-27: For each protein, compute restricted cubic spline
# knots from the full 3177-sample data rather than the filtered modeling
# data.  This keeps model fitting restricted to unimputed-genotype samples
# while allowing protein-specific missing/censored values in the full cohort
# to determine the x-axis knot locations.
get_full3177_spline_info_for_peptide <- function(pep_name, traits_df_for_knots, protein_df_for_knots) {
  knot_dat <- traits_df_for_knots %>%
    dplyr::select(sample_id, EY0) %>%
    dplyr::left_join(
      protein_df_for_knots %>% dplyr::select(sample_id, dplyr::all_of(pep_name)),
      by = "sample_id") %>%
    dplyr::filter(complete.cases(.))

  if (nrow(knot_dat) < 3L || length(unique(knot_dat$EY0)) < 3L) {
    stop(sprintf("Too few full-3177 complete-case EYO values to compute nk=3 spline knots for %s", pep_name))
  }

  splinefit_for_knots <- Hmisc::rcspline.eval(knot_dat$EY0, nk = 3, norm = 2, pc = FALSE, inclx = TRUE)

  list(
    knots = attr(splinefit_for_knots, "knots"),
    n_complete_full3177 = nrow(knot_dat),
    eyo_range_full3177 = range(knot_dat$EY0, na.rm = TRUE)
  )
}
########################################################################

#pep_names <- names(protein_df)[2:dim(protein_df)[2]]
pep_names <- names(clean_mat)

clean_mat<-cbind(data.frame(sample_id=protein_df$sample_id),clean_mat)
# CHANGED 2026-05-27: Add sample_id to the full-3177 knot matrix so the same
# join code can be used for model rows and full-data knot rows.
clean_mat_full3177_for_knots <- cbind(data.frame(sample_id=protein_df_full$sample_id),
                                      clean_mat_full3177_for_knots)


# ---- STEP 8 (per assay, run in parallel below). one_pepSTAN(): builds the
# modeling data frame for one assay, evaluates the restricted-cubic-spline
# EYO basis using full-cohort knots, fits the Bayesian e4/e4 x EYO spline
# model with rstanarm::stan_glm(), and saves the fit to its own RDS file. ----
one_pepSTAN <- function(track, traits_df, protein_df,
                         traits_df_for_knots, protein_df_for_knots) {
  pep_name <- pep_names[track]

  dat <- traits_df %>%
    select(sample_id, ApoE_Indicator, EY0) %>%
    left_join(
      protein_df %>% select(sample_id, all_of(pep_name)),
      by = "sample_id") %>%
    filter(complete.cases(.)) %>%
    select(-contains("sample_id"))

  dat <- sanitize_names(dat)

  # CHANGED 2026-05-27: Estimate knots from the full 3177-sample data for
  # this protein, then evaluate the spline basis only on the filtered
  # unimputed-genotype model data.  Previously rcspline.eval(dat$EY0, nk=3)
  # chose knots from the modeling rows themselves.
  spline_info   <- get_full3177_spline_info_for_peptide(
    pep_name              = pep_name,
    traits_df_for_knots   = traits_df_for_knots,
    protein_df_for_knots  = protein_df_for_knots)
  spline_knots  <- spline_info$knots
  splinefit     <- Hmisc::rcspline.eval(dat$EY0, knots = spline_knots,
                                        norm = 2, pc = FALSE, inclx = TRUE)
  cubic_spline_X  <- as.data.frame(splinefit)
  names(cubic_spline_X) <- c("EYO_Spline_Linear", "EYO_Spline_Cubic")

  dat <- cbind(dat, cubic_spline_X)

  outcome    <- names(dat)[3]
  variables1 <- c("EYO_Spline_Linear", "EYO_Spline_Cubic",
                  "ApoE_Indicator",
                  "EYO_Spline_Linear*ApoE_Indicator",
                  "EYO_Spline_Cubic*ApoE_Indicator")
  f1 <- as.formula(paste(outcome, paste(variables1, collapse = "+"), sep = "~"),
                   env = baseenv())

  set.seed(track)
  env        <- new.env()
  env$dat    <- dat
  env$f1     <- f1
  env$track  <- track

  stan_BL_1 <- with(env, {
    stan_glm(f1,
             data    = dat,
             family  = gaussian(),
             chains  = 8,
             cores   = 4,
             iter    = 10000,
             thin    = 10,
             refresh = 0,
             seed    = track)
  })

  fil_stan_1 <- file.path(".",
                          paste0(outcome, "_stan_glm.rds"))

  strip_formula_envs <- function(x) {
    rapply(x,
           f = function(obj) {
             attr(obj, ".Environment") <- baseenv()
             obj
           },
           classes = c("formula", "terms"),
           how     = "replace")
  }

  stan_BL_1 <- strip_formula_envs(stan_BL_1)
  # CHANGED 2026-05-27: Store the full-3177 knot metadata with each model so
  # downstream plotting can evaluate predictions using the exact spline basis
  # used during model fitting.
  attr(stan_BL_1, "full3177_spline_knots") <- spline_info$knots
  attr(stan_BL_1, "full3177_spline_n_complete") <- spline_info$n_complete_full3177
  attr(stan_BL_1, "full3177_spline_eyo_range") <- spline_info$eyo_range_full3177
  saveRDS(stan_BL_1, fil_stan_1)

  data.frame(OriginalName = pep_name,
             CleanedName  = names(dat)[3],
             stringsAsFactors = FALSE)
}

library(doParallel)
library(foreach)

# ---- STEP 8 (continued). Run one_pepSTAN() for every assay in parallel
# (10 worker processes x 4 STAN threads each = 40 threads), fitting the
# Bayesian e4/e4-by-EYO-spline model and saving one RDS per assay. ----
ncore <- 10  # 40 total threads will be used with 10 core setting b/c each STAN run MC uses 4 threads
cl    <- makeCluster(ncore)
registerDoParallel(cl)

worker_pkgs <- c("tidyverse", "rstanarm", "Hmisc", "rstan", "dplyr")

results <- foreach(track          = 1:length(pep_names),
                   .packages      = worker_pkgs,
                   .export        = c("one_pepSTAN", "sanitize_names",
                                      "get_full3177_spline_info_for_peptide",
                                      "pep_names", "numericMeta_3177_trait_full",
                                      "clean_mat_full3177_for_knots"),
                   .errorhandling = "pass") %dopar% {
#  one_pepSTAN(track, numericMeta_3177_trait, protein_df,
#              numericMeta_3177_trait_full, protein_df_full)
  one_pepSTAN(track, numericMeta_3177_trait, clean_mat,
              numericMeta_3177_trait_full, clean_mat_full3177_for_knots)
}

stopImplicitCluster()


# ---- STEP 9. Check whether any assay's saved RDS is missing or whether its
# missingness pattern changed relative to a prior run (process_again); none
# are flagged in this run (table below reads all FALSE / 7345). ----
processed_names <- gsub("_stan_glm.rds", "", list.files("./"))
name_check_df<-data.frame(assay=pep_names,
                          clean_name=gsub("^_|_$", "",  gsub("_+", "_",  gsub("[^[:alnum:]]", "_", pep_names)))
                         )
name_check_df$process_again<- !name_check_df$clean_name %in% processed_names
name_check_df$process_again[which(name_check_df$assay %in% assays.2rerun)] <- TRUE

table(name_check_df$process_again)
#FALSE
# 7345

#previously:
#FALSE  TRUE
# 7322    23

name_check_df$assay[which(name_check_df$process_again)]==pep_names[which(name_check_df$process_again)]       
#logical(0)  # previously all TRUE



## <SKIP> - this block is a CONDITIONAL fallback, not dead code: it only
## actually refits assays listed in name_check_df$process_again (TRUE for
## assays whose RDS is missing or whose 5xSD-outlier missingness pattern
## changed since the prior run). In this run no assays met that condition,
## so the foreach loop below has zero iterations and effectively no-ops.
ncore <- 8
cl    <- makeCluster(ncore)
registerDoParallel(cl)

worker_pkgs <- c("tidyverse", "rstanarm", "Hmisc", "rstan", "dplyr")

results2 <- foreach(track          = which(name_check_df$process_again),
                   .packages      = worker_pkgs,
                   .export        = c("one_pepSTAN", "sanitize_names",
                                      "get_full3177_spline_info_for_peptide",
                                      "pep_names", "numericMeta_3177_trait_full",
                                      "protein_df_full"),
                   .errorhandling = "pass") %dopar% {
  one_pepSTAN(track, numericMeta_3177_trait, protein_df,
              numericMeta_3177_trait_full, protein_df_full)
  # for (track in name_check_df$assay[which(name_check_df$process_again)]) one_pepSTAN(track,numericMeta_3177_trait, protein_df) #protein_df has outliers for second pass left in these 15 assays.2rerun
}

stopImplicitCluster()


## <END SKIP>


valid_names <- gsub("_stan_glm.rds", "", list.files("./"))
name_check_df$final.files<-name_check_df$clean_name %in% valid_names
table(name_check_df$final.files)
#TRUE
#7345

# ---- STEP 10. Refresh the OriginalName <-> CleanedName lookup table used by
# the companion plotting script (5e2) to locate each assay's saved RDS file,
# restricting it to assays that have a valid saved model in this run. ----
#from prior run:
name_match_table<-readRDS("./name_match_table.RDS")
name_match_table<-name_match_table[which(name_match_table$CleanedName %in% valid_names),]
dim(name_match_table)
# 7345   2

saveRDS(name_match_table,
        "./name_match_table.RDS")


BL_traits<-numericMeta_3177_trait
BL_traits_pep<-protein_df
# CHANGED 2026-05-27: Ensure the unimputed plotting data has EYO populated
# from the filtered traits; the original code referenced BL_traits$EYO even
# though EY0 was the column created above.
BL_traits$EYO <- BL_traits$EY0
BL_traits_pep$EYO<-BL_traits$EYO

# CHANGED 2026-05-27: Save full-3177 EYO/protein objects for the plotting
# script so it can recompute full-data spline knots before drawing ribbons.
# ---- STEP 11. Save filtered (modeling) and full-cohort (knot-estimation)
# trait/protein objects together for the companion plotting script, 5e2. ----
BL_traits_full3177_for_knots <- numericMeta_3177_trait_full
BL_traits_pep_full3177_for_knots <- protein_df_full
BL_traits_full3177_for_knots$EYO <- BL_traits_full3177_for_knots$EY0
BL_traits_pep_full3177_for_knots$EYO <- BL_traits_full3177_for_knots$EYO
save(BL_traits, BL_traits_pep,
     BL_traits_full3177_for_knots, BL_traits_pep_full3177_for_knots,
     file="assays+traits_forPlots.RData")
