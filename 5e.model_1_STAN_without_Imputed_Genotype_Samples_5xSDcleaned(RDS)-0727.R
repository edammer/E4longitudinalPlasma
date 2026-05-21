################################################################################################
# Section 5 longitudinal STAN pipeline - annotated copy
# Source file: 5e.model_1_STAN_without_Imputed_Genotype_Samples_5xSDcleaned(RDS)-0727.R
# Pipeline process 8/12: Model 1 sensitivity analysis excluding imputed genotype samples
#
# Purpose: Assess robustness of the Model 1 longitudinal APOE contrast after removing samples
# lacking directly observed APOE.mapped genotype.
#
# Input:  _numericMeta_3177_trait.RDS
# Input:  _full_3177_protein_dft.RDS
# Output: noImpute.simple.rerun/*_stan_glm.rds
# Output: noImpute.simple.rerun/name_match_table.RDS
# Output: assays+traits_forPlots.RData
#
# Major analysis steps in this script:
#   1. Remove samples with NA APOE.mapped from both metadata and protein data.
#   2. Verify remaining sample order and derive ApoE_Indicator and EYO.
#   3. Reapply 5-SD outlier/missingness filtering as needed.
#   4. Fit Model 1 across outcomes without imputed genotype samples.
#   5. Save rerun model fits and the filtered traits/outcome objects for plotting.
#
# Cleanup/annotation notes:
#   - This is a cleaned, commented copy of the uploaded script; analysis logic and
#     parameter values were not intentionally changed.
#   - Files were decoded from the uploaded Windows/CP1252 text and written as UTF-8.
#   - No explicit "not run below here" block was detected in this prefix-5 file set.
################################################################################################


# ----------------------------------------------------------------------------------------
# Load packages for no-imputed-genotype Model 1 sensitivity fitting.
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
# Load the preprocessed primary cohort traits and outcomes.
# ----------------------------------------------------------------------------------------
# 1. Load pre-processed RDS files
numericMeta_3177_trait <- readRDS("_numericMeta_3177_trait.RDS")
protein_df             <- readRDS("_full_3177_protein_dft.RDS")
# protein_df already contains MMSE, cdr, MEs, and sample_id as first column


# ----------------------------------------------------------------------------------------
# Remove samples without directly observed APOE.mapped genotype.
# ----------------------------------------------------------------------------------------
# 2. Remove samples with missing APOE.mapped
na_apoe_idx <- which(is.na(numericMeta_3177_trait$APOE.mapped))
message(sprintf("Removing %d sample(s) with NA APOE.mapped.", length(na_apoe_idx)))
#Removing 708 sample(s) with NA APOE.mapped.

if (length(na_apoe_idx) > 0) {
  na_sample_ids          <- rownames(numericMeta_3177_trait)[na_apoe_idx]
  numericMeta_3177_trait <- numericMeta_3177_trait[-na_apoe_idx, ]
  protein_df             <- protein_df[!protein_df$sample_id %in% na_sample_ids, ]
}

# 3. Verify alignment between metadata and protein matrix
stopifnot(all(protein_df$sample_id == numericMeta_3177_trait$sample_id))

# ----------------------------------------------------------------------------------------
# Verify metadata and outcome sample order after filtering.
# ----------------------------------------------------------------------------------------

# 4. Set derived traits
numericMeta_3177_trait$ApoE_Indicator <- 0

# ----------------------------------------------------------------------------------------
# Define ApoE_Indicator and EYO after genotype filtering.
# ----------------------------------------------------------------------------------------
numericMeta_3177_trait$ApoE_Indicator[numericMeta_3177_trait$APOE.mapped.predicted == "e4/e4"] <- 1
table(numericMeta_3177_trait$ApoE_Indicator)
#   0    1
#2106  363


numericMeta_3177_trait$EY0 <- (65.6 - numericMeta_3177_trait$age_at_visit) * (-1)
range(numericMeta_3177_trait$EY0)
#-45.6  24.4


# ----------------------------------------------------------------------------------------
# Utility: sanitize outcome names.
# ----------------------------------------------------------------------------------------


tail(colnames(protein_df))
#[1] "magenta" "brown"   "black"   "pink"    "MMSE"    "cdr"
# avoid 5xSD Z removal with MMSE and cdr in place


## Already performed on original data ingress before saving to .RDS. Trying again here.
## Data cleaning - 5SD from mean max within protein, then max 20% per genotype group NA; add MMSE and cdr
clean_mat <- as.data.frame(t(apply(protein_df[,2:(ncol(protein_df)-2)], 1, function(row) {
    ## winsorise at +/- 5 SD
  z   <- abs(row - mean(row, na.rm = TRUE))

# ----------------------------------------------------------------------------------------
# Define outcome list.
# ----------------------------------------------------------------------------------------
  row[z > 5 * sd(row, na.rm = TRUE)] <- NA                 # outliers -> NA


# ----------------------------------------------------------------------------------------
# Worker function for fitting no-imputation Model 1.
# ----------------------------------------------------------------------------------------
  row }                                              # otherwise return the cleaned row
)))
## restore dimnames -------------------------------------------------------
rownames(clean_mat) <- rownames(protein_df)

# ----------------------------------------------------------------------------------------
# Join traits and one outcome, then keep complete cases.
# ----------------------------------------------------------------------------------------
colnames(clean_mat) <- colnames(protein_df)[2:(ncol(protein_df)-2)]


## grouping factor – one value per column in the expression matrix
grp <- numericMeta_3177_trait$APOE.mapped.predicted
stopifnot(length(grp) == ncol(protein_df))   # safety

## pre-compute the column indices that belong to each group

# ----------------------------------------------------------------------------------------
# Generate EYO spline basis.
# ----------------------------------------------------------------------------------------
idx_by_grp <- split(seq_along(grp), grp)

## rows that fail the “<=20% NA in every group” rule  ----------------------
bad_row_idx <- which(
  apply(clean_mat, 1, function(v) {
    any(                                           # if *any* group exceeds 20% NA
      vapply(idx_by_grp,
             function(ix) mean(is.na(v[ix])) > 0.20,

# ----------------------------------------------------------------------------------------
# Model 1 formula used for the no-imputation sensitivity analysis.
# ----------------------------------------------------------------------------------------
             logical(1))
    )
  })
)

bad_row_idx        # numeric vector of offending row indices                                    # -> row is discarded
#named integer(0)

## If any found, remove by %in% rownames(clean_mat)[bad_row_idx]
if(length(bad_row_idx)>0) clean_mat<-clean_mat[which(!rownames(clean_mat) %in% rownames(clean_mat)[bad_row_idx]),]


# ----------------------------------------------------------------------------------------
# Fit Bayesian Gaussian Model 1.
# ----------------------------------------------------------------------------------------


assays.2rerun <- colnames(clean_mat)[which( !apply(clean_mat,2,function(x) length(which(is.na(x)))) == apply(protein_df[,2:(ncol(protein_df)-2)],2,function(x) length(which(is.na(x)))) )]
assays.2rerun
# [1] "NCF2|P19878"                       "ENO1|P06733"                       "CHGA|P10645^SL002762@seq.11184.51" "CDC25A|P30304"
# [5] "SNAP25|P60880"                     "MRPL33|O75394"                     "IDH1|O75874"                       "SERPINA5|P05154"
# [9] "CNDP1|Q96KN2^SL006694@seq.5456.59" "LEAP2|Q969E1"                      "MRPL58|Q14197"                     "CNDP1|Q96KN2^SL006694@seq.7870.8"
#[13] "OLFM2|O95897^SL012399@seq.8295.16" "CHGA|P10645^SL002762@seq.8476.11"  "PLXDC1|Q8IUK5"


########################################################################

# ----------------------------------------------------------------------------------------
# Save compact no-imputation Stan fit.
# ----------------------------------------------------------------------------------------
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

#pep_names <- names(protein_df)[2:dim(protein_df)[2]]
pep_names <- names(clean_mat)

clean_mat<-cbind(data.frame(sample_id=protein_df$sample_id),clean_mat)

# ----------------------------------------------------------------------------------------
# Start the parallel backend.
# ----------------------------------------------------------------------------------------


one_pepSTAN <- function(track, traits_df, protein_df) {
  pep_name <- pep_names[track]

  dat <- traits_df %>%
    select(sample_id, ApoE_Indicator, EY0) %>%

# ----------------------------------------------------------------------------------------
# Run Model 1 across outcomes.
# ----------------------------------------------------------------------------------------
    left_join(
      protein_df %>% select(sample_id, all_of(pep_name)),
      by = "sample_id") %>%
    filter(complete.cases(.)) %>%
    select(-contains("sample_id"))

  dat <- sanitize_names(dat)

  splinefit       <- rcspline.eval(dat$EY0, nk = 3, norm = 2, pc = FALSE, inclx = TRUE)
  cubic_spline_X  <- as.data.frame(splinefit)
  names(cubic_spline_X) <- c("EYO_Spline_Linear", "EYO_Spline_Cubic")

  dat <- cbind(dat, cubic_spline_X)

  outcome    <- names(dat)[3]
  variables1 <- c("EYO_Spline_Linear", "EYO_Spline_Cubic",
                  "ApoE_Indicator",
                  "EYO_Spline_Linear*ApoE_Indicator",
                  "EYO_Spline_Cubic*ApoE_Indicator")
  f1 <- as.formula(paste(outcome, paste(variables1, collapse = "+"), sep = "~"),

# ----------------------------------------------------------------------------------------
# Save name matching and filtered analysis objects for plotting.
# ----------------------------------------------------------------------------------------
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

  fil_stan_1 <- file.path("f:/", "OneDrive - Emory", "Legacy", "e4_homozygoteStudy",
                          "DL", "noImpute.simple.rerun",
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
  saveRDS(stan_BL_1, fil_stan_1)

  data.frame(OriginalName = pep_name,
             CleanedName  = names(dat)[3],
             stringsAsFactors = FALSE)
}

library(doParallel)
library(foreach)

ncore <- 8
cl    <- makeCluster(ncore)
registerDoParallel(cl)

worker_pkgs <- c("tidyverse", "rstanarm", "Hmisc", "rstan", "dplyr")

results <- foreach(track          = 1:length(pep_names),
                   .packages      = worker_pkgs,
                   .export        = c("one_pepSTAN", "sanitize_names", "pep_names"),
                   .errorhandling = "pass") %dopar% {
#  one_pepSTAN(track, numericMeta_3177_trait, protein_df)
  one_pepSTAN(track, numericMeta_3177_trait, clean_mat)
}

stopImplicitCluster()


processed_names <- gsub("_stan_glm.rds", "", list.files("./noImpute.simple.rerun"))
name_check_df<-data.frame(assay=pep_names,
                          clean_name=gsub("^_|_$", "",  gsub("_+", "_",  gsub("[^[:alnum:]]", "_", pep_names)))
                         )
name_check_df$process_again<- !name_check_df$clean_name %in% processed_names
name_check_df$process_again[which(name_check_df$assay %in% assays.2rerun)] <- TRUE

table(name_check_df$process_again)
#FALSE  TRUE
# 7322    23
name_check_df$assay[which(name_check_df$process_again)]==pep_names[which(name_check_df$process_again)]
#all TRUE


ncore <- 8
cl    <- makeCluster(ncore)
registerDoParallel(cl)

worker_pkgs <- c("tidyverse", "rstanarm", "Hmisc", "rstan", "dplyr")

results2 <- foreach(track          = which(name_check_df$process_again),
                   .packages      = worker_pkgs,
                   .export        = c("one_pepSTAN", "sanitize_names", "pep_names"),
                   .errorhandling = "pass") %dopar% {
  one_pepSTAN(track, numericMeta_3177_trait, protein_df)
  # for (track in name_check_df$assay[which(name_check_df$process_again)]) one_pepSTAN(track,numericMeta_3177_trait, protein_df) #protein_df has outliers for second pass left in these 15 assays.2rerun
}

stopImplicitCluster()


valid_names <- gsub("_stan_glm.rds", "", list.files("./noImpute.simple.rerun"))
name_check_df$final.files<-name_check_df$clean_name %in% valid_names
table(name_check_df$final.files)
#TRUE
#7345

#filtered_results <- lapply(results, function(x) {
#  # Skip if x is not a data frame or has fewer than 2 columns
#  if (!is.data.frame(x) || ncol(x) < 2 || nrow(x) < 1) return(NULL)
#
#  key <- x[1, 2]
#
#  # Skip if the key is not in the valid list
#  if (! key %in% valid_names) return(NULL)
#
#  # Otherwise keep the data frame
#  x
#})

#name_match_table         <- do.call(rbind, filtered_results)
#colnames(name_match_table) <- c("OriginalName", "CleanedName")

#saveRDS(name_match_table,
#        "F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/noImpute.simple.rerun/name_match_table.RDS")

#from prior run:
name_match_table<-readRDS("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/noImpute.simple.rerun/name_match_table.RDS")
name_match_table<-name_match_table[which(name_match_table$CleanedName %in% valid_names),]
dim(name_match_table)
# 7345   2

saveRDS(name_match_table,
        "F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/noImpute.simple.rerun/name_match_table.RDS")


BL_traits<-numericMeta_3177_trait
BL_traits_pep<-protein_df
BL_traits_pep$EYO<-BL_traits$EYO
save(BL_traits,BL_traits_pep,file="assays+traits_forPlots.RData")
