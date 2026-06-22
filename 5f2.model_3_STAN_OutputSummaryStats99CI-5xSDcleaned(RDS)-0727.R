################################################################################################
# Section 5 longitudinal STAN pipeline - annotated copy
# Source file: 5f2.model_3_STAN_OutputSummaryStats99CI-5xSDcleaned(RDS)-0727.R
# Pipeline process 11/12 (sex-coding corrected summary): Post-hoc Model 2.3 joint posterior significance table generation with corrected sex coding
#
# Purpose: Use saved Model 2.3 fits, without refitting, to test the composite genotype-by-sex-
# by-EYO posterior interaction at each EYO bin. Category labels aligned to the intended Sex.int
# coding of 0 = female and 1 = male.
#
# Input:  /home/labshare/genoSex/*_with_APOE.Sex_GenotypeXSex_stan_glm.rds
# Input:  /home/labshare/genoSex/name_match_table.RDS
# Input:  _numericMeta_3177_trait.RDS
# Output: posthoc_results/genoSex_joint_posterior_results.RData
# Output: posthoc_results/1_pval_sig_by_EYO_bin.csv
# Output: posthoc_results/2_effect_size_by_EYO_bin.csv
# Output: posthoc_results/3_min_pval_per_protein.csv
# Output: posthoc_results/4_first_sig_category.csv
#
# Major analysis steps in this script:
#   1. Repeat the joint posterior test using the corrected female/male coding assumption.
#   2. Retain p-value, effect-size, minimum-p-value, and first-significant-category outputs.
#   3. Write CSV and RData outputs for downstream plots and interpretation.
#
# Cleanup/annotation notes:
#   - This is a cleaned, commented copy of the uploaded script; analysis logic and
#     parameter values were not intentionally changed.
#   - Files were decoded from the uploaded Windows/CP1252 text and written as UTF-8.
#   - No explicit "not run below here" block was detected in this prefix-5 file set.
################################################################################################


# ----------------------------------------------------------------------------------------
# Header documents the corrected post-hoc joint posterior output workflow.
# ----------------------------------------------------------------------------------------
#I have run the  parallelized STAN model generation for model 2.3 on all protein assays. As you said, because the genotype × sex interaction is now spread across three coefficients, extracting a single "interaction p-value" requires a joint posterior test rather than reading off a single coefficient's credible interval. I would like to run a full joint test (e.g., computing the posterior probability that the linear combination of the three terms is non-zero at a given EYO value) as the most rigorous approach, to be done post-hoc from the saved RDS on the model for each protein assay, in a fast parallelized framework, without re-running the models.  Provide R code to do this.
#Outputs should include: (1) a data frame of proteins (rows) x EYO bins (columns) with p values for joint posterior interaction significance when below 0.05; (2) a data frame of proteins (rows) x EYO bins (columns) with magnitude or signed effect size of the genotype x sex x EYO bin interaction at that bin for that protein; (3) a vector of minimum posterior interaction significances (p values) across all EYO bins, one per protein in the same order as rows of the prior 2 outputs; (4) finally, output a vector which provides a named string for the category of the first significant interval (earliest EYO bin) reaching significance for posterior interaction significance in each protein, e.g. "Male.e4.up", "Female.e4.up", "Male.e3.up", "Female.e3.up". If no bins for a protein reach significance for the posterior interaction, output "None" for that position in the vector. Save the 4 outputs to an RData file and each to a .csv file.

# =============================================================================
# Post-hoc joint posterior test: genotype × sex × EYO interaction (Model 2.3)
#
# Background
# ----------
# Model 2.3 distributes the genotype × sex interaction across THREE coefficients:
#
#   b1  =  ApoE_Indicator:Sex                      (time-invariant offset)
#   b2  =  EYO_Spline_Linear:ApoE_Indicator:Sex    (linear EYO modulation)
#   b3  =  EYO_Spline_Cubic:ApoE_Indicator:Sex     (cubic EYO modulation)
#
# The full interaction effect at a given EYO value t is therefore:
#   interaction(t) = b1  +  b2 * spline_linear(t)  +  b3 * spline_cubic(t)
#
# A two-tailed Bayesian p-value is computed as:
#   p(t) = 2 * min( P(interaction(t) > 0 | data),  P(interaction(t) < 0 | data) )
#
# Category assignment (which group is highest at the first significant EYO bin)
# is determined by computing posterior means for all four groups from all
# relevant model coefficients, then naming whichever group has the highest mean.
#
# ASSUMPTION: Sex is coded as  0 = Female,  1 = Male.
# If Sex.int coding is reversed, swap Male/Female in the category labels
# at the bottom of the script.
# =============================================================================

library(tidyverse)
library(Hmisc)
library(doParallel)
library(foreach)

# 0. Paths & setup

# ----------------------------------------------------------------------------------------
# Set paths for saved Model 2.3 fits and post-hoc outputs.
# ----------------------------------------------------------------------------------------
rds_dir    <- "/home/labshare/genoSex"
output_dir <- "/home/labshare/genoSex/posthoc_results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Metadata needed for EYO distribution (to reconstruct reference spline knots)

# ----------------------------------------------------------------------------------------
# Load metadata and name matching.
# ----------------------------------------------------------------------------------------
numericMeta <- readRDS("_numericMeta_3177_trait.RDS")

# Name-match table maps cleaned model-output names back to original protein names
name_match <- readRDS(file.path(rds_dir, "name_match_table.RDS"))

# 1. EYO bins and reference spline basis
# Use the full EYO distribution for knot placement, consistent with model fitting.
# Per-protein knots may differ trivially due to complete-case filtering, but
# with nk=3 knots and thousands of samples this difference is negligible.

# ----------------------------------------------------------------------------------------
# Define the 0.5-year EYO grid and spline basis.
# ----------------------------------------------------------------------------------------
eyo_vals   <- numericMeta$EY0[!is.na(numericMeta$EY0)]

#eyo_bins   <- seq(floor(min(eyo_vals)), ceiling(max(eyo_vals)), by = 1L)
#n_bins     <- length(eyo_bins)
# EYO bins may be fractional, e.g. -10.5, 1.5, 5.5.  above 2 lines only gave 72 bins.
# Set this to the bin spacing used for model interpretation/output.
eyo_bin_width <- 0.5

floor_to_step <- function(x, step) floor(x / step) * step
ceiling_to_step <- function(x, step) ceiling(x / step) * step

eyo_bins <- seq(
  from = floor_to_step(min(eyo_vals), eyo_bin_width),
  to   = ceiling_to_step(max(eyo_vals), eyo_bin_width),
  by   = eyo_bin_width
)

# Avoid floating-point artifacts such as 1.499999999
eyo_bins <- round(eyo_bins, digits = 1L)
n_bins <- length(eyo_bins)


ref_spline    <- rcspline.eval(eyo_vals, nk = 3, norm = 2, pc = FALSE, inclx = TRUE)
ref_knots     <- attr(ref_spline, "knots")   # 3 knot positions

# Evaluate the same spline basis at each integer EYO bin point.
# With inclx=TRUE and nk=3, the result has 2 columns:
#   col 1 = linear basis (normalised x)
#   col 2 = single cubic restricted spline term
spline_bins        <- rcspline.eval(eyo_bins, knots = ref_knots, norm = 2, inclx = TRUE)
colnames(spline_bins) <- c("lin", "cub")

# Human-readable column labels for output data frames
#bin_labels <- paste0("EYO_", sprintf("%+.1f", round(eyo_bins,digits=1)))   # e.g. "EYO_-46", "EYO_+0"
fmt_eyo <- function(x) {
  out <- formatC(x, format = "f", digits = 6L)
  out <- sub("0+$", "", out)
  out <- sub("\\.$", "", out)
  ifelse(x >= 0, paste0("+", out), out)
}

bin_labels <- paste0("EYO_", fmt_eyo(eyo_bins))

# 2. Discover model RDS files and recover protein names
rds_files <- list.files(

# ----------------------------------------------------------------------------------------
# Discover saved Model 2.3 fit files.
# ----------------------------------------------------------------------------------------
  rds_dir,
  pattern    = "_with_APOE\\.Sex_GenotypeXSex_stan_glm\\.rds$",
  full.names = TRUE
)
stopifnot("No model 2.3 RDS files found — check rds_dir." = length(rds_files) > 0)

cleaned_from_file <- sub(
  "_with_APOE\\.Sex_GenotypeXSex_stan_glm\\.rds$", "",
  basename(rds_files)
)
protein_names <- name_match$OriginalName[match(cleaned_from_file, name_match$CleanedName)]
n_prot <- length(rds_files)

message(sprintf("Found %d model 2.3 RDS files covering %d uniquely matched proteins.",
                n_prot, sum(!is.na(protein_names))))

# 3. Helper: find a coefficient column by exact or partial name

# ----------------------------------------------------------------------------------------
# Define coefficient lookup helpers.
# ----------------------------------------------------------------------------------------
# Placed outside foreach so it can be .export-ed to workers.
find_coef <- function(cn, target) {
  # Try exact match first (fastest, most reliable)
  i <- match(target, cn)
  if (!is.na(i)) return(i)
  # Fall back to fixed substring search
  hits <- grep(target, cn, fixed = TRUE)
  if (length(hits) == 1L) return(hits[1L])
  NA_integer_
}

# 4. Parallelised post-hoc extraction
ncore <- max(1L, parallel::detectCores() - 2L)

# ----------------------------------------------------------------------------------------
# Start parallel processing.
# ----------------------------------------------------------------------------------------
cl    <- makeCluster(ncore)
registerDoParallel(cl)


# ----------------------------------------------------------------------------------------
# Evaluate each Model 2.3 fit.
# ----------------------------------------------------------------------------------------
worker_results <- foreach(
  i              = seq_len(n_prot),
  .packages      = "rstanarm",
  .export        = c("rds_files", "spline_bins", "n_bins", "find_coef"),
  .errorhandling = "pass"
) %dopar% {

  # load model
  mod <- tryCatch(readRDS(rds_files[i]), error = function(e) NULL)
  na_out <- list(pvals   = rep(NA_real_, n_bins),
                 effects = rep(NA_real_, n_bins),
                 groups  = matrix(NA_real_, nrow = n_bins, ncol = 4L))
  if (is.null(mod)) return(na_out)

  post <- as.matrix(mod)   # posterior draws (rows) × coefficients (cols)
  cn   <- colnames(post)

  # locate the three interaction coefficients (joint posterior test)

# ----------------------------------------------------------------------------------------
# Locate composite interaction coefficients.
# ----------------------------------------------------------------------------------------
  idx_b1 <- find_coef(cn, "ApoE_Indicator:Sex")
  idx_b2 <- find_coef(cn, "EYO_Spline_Linear:ApoE_Indicator:Sex")
  idx_b3 <- find_coef(cn, "EYO_Spline_Cubic:ApoE_Indicator:Sex")

  if (any(is.na(c(idx_b1, idx_b2, idx_b3)))) return(na_out)

  b1 <- post[, idx_b1]
  b2 <- post[, idx_b2]
  b3 <- post[, idx_b3]

  # locate additional coefficients for 4-group posterior means
  # Returns a column of zeros if a coefficient is absent (safe fallback)
  get_col <- function(target) {
    idx <- find_coef(cn, target)
    if (is.na(idx)) rep(0, nrow(post)) else post[, idx]
  }

# ----------------------------------------------------------------------------------------
# Extract group-specific posterior components using corrected Sex.int coding.
# ----------------------------------------------------------------------------------------

  b_int   <- get_col("(Intercept)")
  b_eL    <- get_col("EYO_Spline_Linear")
  b_eC    <- get_col("EYO_Spline_Cubic")
  b_apoe  <- get_col("ApoE_Indicator")
  b_aL    <- get_col("EYO_Spline_Linear:ApoE_Indicator")
  b_aC    <- get_col("EYO_Spline_Cubic:ApoE_Indicator")
  b_sex   <- get_col("Sex")
  b_sL    <- get_col("EYO_Spline_Linear:Sex")
  b_sC    <- get_col("EYO_Spline_Cubic:Sex")

  # per-bin calculations
  pvals   <- numeric(n_bins)
  effects <- numeric(n_bins)

  # grp_mat: posterior group means per bin
  # Column order: [1] Female.e3  [2] Male.e3  [3] Female.e4  [4] Male.e4
  # (ApoE_Indicator: 0 = non-e4/e4 labelled "e3"; 1 = e4/e4)
  # (Sex:            0 = Female;                     1 = Male)
  grp_mat <- matrix(0, nrow = n_bins, ncol = 4L)

  for (j in seq_len(n_bins)) {
    sl <- spline_bins[j, "lin"]
    sc <- spline_bins[j, "cub"]

    # Full genotype × sex interaction posterior at this EYO
    intx <- b1  +  b2 * sl  +  b3 * sc
    p_pos       <- mean(intx > 0)
    pvals[j]    <- 2 * min(p_pos, 1 - p_pos)   # two-tailed Bayesian p-value
    effects[j]  <- mean(intx)                   # posterior mean as signed effect

    # Shared EYO trajectory (identical across all 4 groups)
    base <- b_int  +  b_eL * sl  +  b_eC * sc

    # Group predicted means (posterior mean of each group's linear predictor)
    # Female   non-e4/e4  (ApoE=0, Sex=0)
    grp_mat[j, 1L] <- mean(base)
    # Male non-e4/e4  (ApoE=0, Sex=1)
    grp_mat[j, 2L] <- mean(base  +  b_sex  +  b_sL * sl  +  b_sC * sc)
    # Female   e4/e4      (ApoE=1, Sex=0)
    grp_mat[j, 3L] <- mean(base  +  b_apoe  +  b_aL * sl  +  b_aC * sc)
    # Male e4/e4      (ApoE=1, Sex=1)
    grp_mat[j, 4L] <- mean(base  +  b_sex  +  b_sL * sl  +  b_sC * sc
                           +  b_apoe  +  b_aL * sl  +  b_aC * sc
                           +  b1      +  b2   * sl  +  b3   * sc)
  }


# ----------------------------------------------------------------------------------------
# Compute posterior effects, p-values, and category calls.
# ----------------------------------------------------------------------------------------
  list(pvals = pvals, effects = effects, groups = grp_mat)
}

stopCluster(cl)

# 5. Assemble result matrices
pval_mat   <- matrix(NA_real_, nrow = n_prot, ncol = n_bins,
                     dimnames = list(protein_names, bin_labels))
effect_mat <- matrix(NA_real_, nrow = n_prot, ncol = n_bins,
                     dimnames = list(protein_names, bin_labels))

for (i in seq_len(n_prot)) {
  r <- worker_results[[i]]
  if (!inherits(r, "error") && !is.null(r[["pvals"]])) {
    pval_mat[i, ]   <- r$pvals
    effect_mat[i, ] <- r$effects
  }
}

# Output 1: p-values (significant bins only; non-sig set to NA)
pval_sig_df              <- as.data.frame(pval_mat)
pval_sig_df[pval_sig_df >= 0.0051] <- NA

# Output 2: signed posterior mean effect size per bin
effect_df <- as.data.frame(effect_mat)

# Output 3: minimum posterior p-value per protein
min_pvals <- setNames(
  apply(pval_mat, 1L, function(r) {
    v <- r[!is.na(r)]
    if (length(v) == 0L) NA_real_ else min(v)
  }),
  protein_names
)

# Output 4: category string for the first significant EYO bin
# Group order in grp_mat: [1] Female.e3  [2] Male.e3  [3] Female.e4  [4] Male.e4
group_category_names <- c("Female.e3.up", "Male.e3.up", "Female.e4.up", "Male.e4.up")

# ----------------------------------------------------------------------------------------
# Assemble matrices and vectors.
# ----------------------------------------------------------------------------------------

# Earliest, i.e. lowest-EYO, significant bin index per protein
first_sig_idx <- vapply(seq_len(n_prot), function(i) {

  pv  <- pval_mat[i, ]
  sig <- which(!is.na(pv) & pv < 0.0051)

  if (length(sig) == 0L) NA_integer_ else sig[1L]

}, FUN.VALUE = integer(1L))

names(first_sig_idx) <- protein_names

# Numeric EYO bin value at first significance; NA if no significant bin
first_sig_EYO_bin <- rep(NA_real_, length(first_sig_idx))
has_sig <- !is.na(first_sig_idx)
first_sig_EYO_bin[has_sig] <- eyo_bins[first_sig_idx[has_sig]]

names(first_sig_EYO_bin) <- protein_names

first_sig_category <- vapply(seq_len(n_prot), function(i) {

  j <- first_sig_idx[i]
  if (is.na(j)) return("None")

  r <- worker_results[[i]]

  if (inherits(r, "error") || is.null(r[["groups"]])) return("None")

  # The group with the highest posterior mean at this bin is named as "up"
  group_category_names[which.max(r$groups[j, ])]

}, FUN.VALUE = character(1L))

names(first_sig_category) <- protein_names

# 6. Report summary to console
n_sig <- sum(!is.na(min_pvals) & min_pvals < 0.0051)
message(sprintf(
  "\nSummary: %d / %d proteins have at least one EYO bin with posterior p < 0.05.",
  n_sig, n_prot))
message("Category breakdown across those proteins:")
print(table(first_sig_category))


# ----------------------------------------------------------------------------------------
# Save corrected post-hoc outputs.
# ----------------------------------------------------------------------------------------
# 7. Save all four outputs
save(pval_sig_df, effect_df, min_pvals,
     first_sig_category, first_sig_EYO_bin,
     file = file.path(output_dir, "genoSex_joint_posterior_results.RData"))

write.csv(pval_sig_df,
          file = file.path(output_dir, "1_pval_sig_by_EYO_bin.csv"),
          na   = "")

write.csv(effect_df,
          file = file.path(output_dir, "2_effect_size_by_EYO_bin.csv"),
          na   = "")

write.csv(
  data.frame(Protein       = protein_names,
             Min_PosteriorP = min_pvals,
             row.names     = NULL),
  file      = file.path(output_dir, "3_min_pval_per_protein.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(Protein            = protein_names,
             First_Sig_Category = first_sig_category,
             First_Sig_EYO_Bin  = first_sig_EYO_bin,
             Min_PosteriorP     = min_pvals,
             row.names          = NULL),
  file      = file.path(output_dir, "4_first_sig_category.csv"),
  row.names = FALSE
)

message("\nAll outputs written to: ", output_dir)
