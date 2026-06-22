################################################################################################
# Section 5 longitudinal STAN pipeline - annotated copy
# Source file: 5f3.Plot_STAN_model_3_99CI_withSex+Genotype+EYO_3way_int_terms-5xSDcleaned(RDS)-0727.R
# Pipeline process 12/12: Model 2.3 single-panel 99% interval plots and category streak outputs
#
# Purpose: Plot per-protein EYO-dependent genotype-by-sex interaction ribbons with category-
# colored significance streaks and export matrix summaries.
#
# Input:  /home/labshare/genoSex/*_with_APOE.Sex_GenotypeXSex_stan_glm.rds
# Input:  /home/labshare/genoSex/name_match_table.RDS
# Input:  _numericMeta_3177_trait.RDS
# Output: posthoc_results/model23_3way_interaction_single_panel_plots/*.pdf
# Output: model23_3way_interaction_plot_summary.csv
# Output: model23_3way_interaction_*_by_EYO.csv
# Output: model23_3way_interaction_single_panel_plot_results.rds
#
# Major analysis steps in this script:
#   1. Reconstruct the EYO grid and spline basis used for model interpretation.
#   2. For each Model 2.3 fit, extract the composite genotype-by-sex interaction posterior curve.
#   3. Compute median, 99% posterior interval, and two-tailed posterior tail probability at each
#      EYO.
#   4. Assign the dominant posterior group using Sex.int coding 0 = female, 1 = male.
#   5. Draw a single-panel ribbon plot with a colored top streak where the 99% interval excludes
#      zero.
#
# Cleanup/annotation notes:
#   - This is a cleaned, commented copy of the uploaded script; analysis logic and
#     parameter values were not intentionally changed.
#   - Files were decoded from the uploaded Windows/CP1252 text and written as UTF-8.
#   - No explicit "not run below here" block was detected in this prefix-5 file set.
################################################################################################


# ----------------------------------------------------------------------------------------
# Define the Model 2.3 single-panel plotting workflow.
# ----------------------------------------------------------------------------------------
# =============================================================================
# Model 2.3 single-panel plots:
# EYO-dependent genotype × sex interaction ribbon + category-colored sig streak
# =============================================================================

library(tidyverse)
library(Hmisc)
library(rstanarm)
library(doParallel)
library(foreach)
library(ggplot2)

# -------------------------------------------------------------------------
# 0. Paths

# ----------------------------------------------------------------------------------------
# Set directories for input models and plot outputs.
# ----------------------------------------------------------------------------------------
# -------------------------------------------------------------------------
rds_dir    <- "/home/labshare/genoSex"
output_dir <- "/home/labshare/genoSex/posthoc_results"

plot_dir <- file.path(output_dir, "model23_3way_interaction_single_panel_plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Metadata and name matching

# ----------------------------------------------------------------------------------------
# Load metadata and name matching.
# ----------------------------------------------------------------------------------------
numericMeta <- readRDS("_numericMeta_3177_trait.RDS")
name_match  <- readRDS(file.path(rds_dir, "name_match_table.RDS"))

# -------------------------------------------------------------------------
# 1. EYO grid with non-integer bins allowed
# -------------------------------------------------------------------------
eyo_col <- intersect(c("EY0", "EYO"), names(numericMeta))[1]
stopifnot("No EYO/EY0 column found in numericMeta." = length(eyo_col) == 1L)

# ----------------------------------------------------------------------------------------
# Create the non-integer EYO grid used for plots.
# ----------------------------------------------------------------------------------------

eyo_vals <- numericMeta[[eyo_col]]
eyo_vals <- eyo_vals[!is.na(eyo_vals)]

eyo_bin_width <- 0.5

floor_to_step <- function(x, step) floor(x / step) * step
ceiling_to_step <- function(x, step) ceiling(x / step) * step

eyo_grid <- seq(
  from = floor_to_step(min(eyo_vals), eyo_bin_width),
  to   = ceiling_to_step(max(eyo_vals), eyo_bin_width),
  by   = eyo_bin_width
)

eyo_grid <- round(eyo_grid, digits = 6L)
n_bins   <- length(eyo_grid)

fmt_eyo <- function(x) {
  out <- formatC(x, format = "f", digits = 6L)
  out <- sub("0+$", "", out)
  out <- sub("\\.$", "", out)
  ifelse(x >= 0, paste0("+", out), out)

# ----------------------------------------------------------------------------------------
# Reconstruct spline basis on the plotting grid.
# ----------------------------------------------------------------------------------------
}

eyo_labels <- paste0("EYO_", fmt_eyo(eyo_grid))

# Reference spline basis, matching model-fitting convention
ref_spline <- Hmisc::rcspline.eval(
  eyo_vals,
  nk    = 3,
  norm  = 2,
  pc    = FALSE,
  inclx = TRUE
)

ref_knots <- attr(ref_spline, "knots")

spline_grid <- Hmisc::rcspline.eval(
  eyo_grid,
  knots = ref_knots,
  norm  = 2,
  pc    = FALSE,
  inclx = TRUE
)

colnames(spline_grid) <- c("lin", "cub")

# -------------------------------------------------------------------------
# 2. Discover Model 2.3 files
# -------------------------------------------------------------------------
rds_files <- list.files(

# ----------------------------------------------------------------------------------------
# Discover saved Model 2.3 RDS files.
# ----------------------------------------------------------------------------------------
  rds_dir,
  pattern    = "_with_APOE\\.Sex_GenotypeXSex_stan_glm\\.rds$",
  full.names = TRUE
)

stopifnot("No Model 2.3 RDS files found — check rds_dir." = length(rds_files) > 0)

cleaned_from_file <- sub(
  "_with_APOE\\.Sex_GenotypeXSex_stan_glm\\.rds$",
  "",
  basename(rds_files)
)

protein_names <- name_match$OriginalName[
  match(cleaned_from_file, name_match$CleanedName)
]

n_prot <- length(rds_files)

message(sprintf(
  "Found %d Model 2.3 RDS files covering %d uniquely matched proteins.",
  n_prot,

# ----------------------------------------------------------------------------------------
# Define helpers for coefficient extraction and file-safe names.
# ----------------------------------------------------------------------------------------
  sum(!is.na(protein_names))
))

# -------------------------------------------------------------------------
# 3. Helpers
# -------------------------------------------------------------------------
find_coef <- function(cn, target) {
  i <- match(target, cn)
  if (!is.na(i)) return(i)

  hits <- grep(target, cn, fixed = TRUE)
  if (length(hits) == 1L) return(hits[1L])

  NA_integer_
}

get_coef_or_zero <- function(post, cn, target) {
  idx <- find_coef(cn, target)
  if (is.na(idx)) rep(0, nrow(post)) else post[, idx]
}

safe_file_stem <- function(x) {
  x <- gsub("[^[:alnum:]_\\.\\-]+", "_", x)

# ----------------------------------------------------------------------------------------
# Define group-category colors for significance streaks.
# ----------------------------------------------------------------------------------------
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# Category colors for the top streak
# These names are independent of the numeric Sex coding, but Section 4B below
# must assign posterior group means using Sex = 0 Female and Sex = 1 Male.
category_cols <- c(
  "Female.e3.up" = "maroon",        #"#ff8c00",
  "Male.e3.up"   = "darkslateblue", #"#6b8e23",
  "Female.e4.up" = "hotpink"        #"indianred3"
  "Male.e4.up"   = "darkviolet",    #"cornflowerblue",
)

group_category_names <- names(category_cols)

# FALSE follows your request exactly:
#   b2 * spline_linear(EYO) + b3 * spline_cubic(EYO)
#

# ----------------------------------------------------------------------------------------
# Configure whether the time-invariant ApoE-by-sex offset is included in the plotted
# curve.
# ----------------------------------------------------------------------------------------
# TRUE below plots the full Model 2.3 genotype × sex interaction:
#   b1 + b2 * spline_linear(EYO) + b3 * spline_cubic(EYO)
include_ApoE_Sex_offset_in_curve <- TRUE
# With FALSE, the plotted ribbon is summing only two EYO-dependent 3-way terms

# -------------------------------------------------------------------------
# 4. Per-protein worker
# -------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Worker function: create one Model 2.3 interaction plot and return matrices.
# ----------------------------------------------------------------------------------------
process_one_model23_plot <- function(i) {

  clean_name   <- cleaned_from_file[i]
  protein_name <- protein_names[i]
  rds_file     <- rds_files[i]

  out_pdf <- file.path(
    plot_dir,
    paste0(safe_file_stem(clean_name), "_model23_3way_interaction_singlePanel.pdf")
  )

  mod <- tryCatch(readRDS(rds_file), error = function(e) e)
  if (inherits(mod, "error") || is.null(mod)) {
    return(list(
      ok      = FALSE,
      protein = protein_name,
      file    = out_pdf,
      error   = paste("readRDS failed:", conditionMessage(mod))
    ))
  }

  post <- tryCatch(as.matrix(mod), error = function(e) e)
  if (inherits(post, "error") || is.null(post)) {
    return(list(
      ok      = FALSE,
      protein = protein_name,
      file    = out_pdf,
      error   = paste("as.matrix failed:", conditionMessage(post))
    ))
  }

  cn <- colnames(post)

  # Model 2.3 genotype × sex terms

# ----------------------------------------------------------------------------------------
# Locate composite interaction coefficients.
# ----------------------------------------------------------------------------------------
  idx_b1 <- find_coef(cn, "ApoE_Indicator:Sex")
  idx_b2 <- find_coef(cn, "EYO_Spline_Linear:ApoE_Indicator:Sex")
  idx_b3 <- find_coef(cn, "EYO_Spline_Cubic:ApoE_Indicator:Sex")

  if (any(is.na(c(idx_b2, idx_b3)))) {
    return(list(
      ok      = FALSE,
      protein = protein_name,
      file    = out_pdf,
      error   = "Missing one or both EYO-dependent 3-way interaction coefficients."
    ))
  }

  b1 <- if (!is.na(idx_b1)) post[, idx_b1] else rep(0, nrow(post))
  b2 <- post[, idx_b2]
  b3 <- post[, idx_b3]

  # -----------------------------------------------------------------------
  # 4A. Composite interaction posterior draws across EYO
  # -----------------------------------------------------------------------
  interaction_draws <-
    tcrossprod(b2, spline_grid[, "lin"]) +
    tcrossprod(b3, spline_grid[, "cub"])

  if (isTRUE(include_ApoE_Sex_offset_in_curve)) {
    interaction_draws <- sweep(interaction_draws, 1L, b1, "+")
  }

  qmat <- apply(
    interaction_draws,
    2L,
    stats::quantile,
    probs = c(0.005, 0.5, 0.995),
    names = FALSE,
    na.rm = TRUE
  )

  p_pos <- colMeans(interaction_draws > 0, na.rm = TRUE)
  pval  <- 2 * pmin(p_pos, 1 - p_pos)


# ----------------------------------------------------------------------------------------
# Build the plotting data frame with posterior median, 99% interval, and p-value.
# ----------------------------------------------------------------------------------------
  plot_df <- data.frame(
    eyo    = eyo_grid,
    lower  = qmat[1L, ],
    median = qmat[2L, ],
    upper  = qmat[3L, ],
    pval   = pval
  )

  plot_df$ci_excludes_zero <- with(plot_df, lower > 0 | upper < 0)

  # -----------------------------------------------------------------------
  # 4B. Determine highest posterior-mean group per EYO bin
  #
  # Correct coding:
  #   Sex = 0: Female
  #   Sex = 1: Male
  #

# ----------------------------------------------------------------------------------------
# Compute posterior group means using Sex.int = 0 female and 1 male.
# ----------------------------------------------------------------------------------------
  # ApoE_Indicator:
  #   0: e3 / non-e4-e4 reference group
  #   1: e4/e4 carrier group
  #
  # Therefore the ApoE_Indicator:Sex interaction terms apply to Male.e4,
  # because that is the group with ApoE_Indicator = 1 and Sex = 1.
  # -----------------------------------------------------------------------
  b_int  <- get_coef_or_zero(post, cn, "(Intercept)")
  b_eL   <- get_coef_or_zero(post, cn, "EYO_Spline_Linear")
  b_eC   <- get_coef_or_zero(post, cn, "EYO_Spline_Cubic")

  b_apoe <- get_coef_or_zero(post, cn, "ApoE_Indicator")
  b_aL   <- get_coef_or_zero(post, cn, "EYO_Spline_Linear:ApoE_Indicator")
  b_aC   <- get_coef_or_zero(post, cn, "EYO_Spline_Cubic:ApoE_Indicator")

  b_sex  <- get_coef_or_zero(post, cn, "Sex")
  b_sL   <- get_coef_or_zero(post, cn, "EYO_Spline_Linear:Sex")
  b_sC   <- get_coef_or_zero(post, cn, "EYO_Spline_Cubic:Sex")

  grp_mat <- matrix(
    NA_real_,
    nrow = n_bins,
    ncol = 4L,
    dimnames = list(NULL, group_category_names)
  )

  for (j in seq_len(n_bins)) {

    sl <- spline_grid[j, "lin"]
    sc <- spline_grid[j, "cub"]

    base <- b_int + b_eL * sl + b_eC * sc

    # Female.e3: ApoE = 0, Sex = 0
    grp_mat[j, "Female.e3.up"] <- mean(
      base,
      na.rm = TRUE
    )

    # Male.e3: ApoE = 0, Sex = 1
    grp_mat[j, "Male.e3.up"] <- mean(
      base +
        b_sex + b_sL * sl + b_sC * sc,
      na.rm = TRUE
    )

    # Female.e4: ApoE = 1, Sex = 0
    grp_mat[j, "Female.e4.up"] <- mean(
      base +
        b_apoe + b_aL * sl + b_aC * sc,
      na.rm = TRUE
    )

    # Male.e4: ApoE = 1, Sex = 1
    # This is the only group receiving the ApoE × Sex and EYO × ApoE × Sex terms.
    grp_mat[j, "Male.e4.up"] <- mean(
      base +
        b_sex  + b_sL * sl + b_sC * sc +
        b_apoe + b_aL * sl + b_aC * sc +
        b1     + b2   * sl + b3   * sc,
      na.rm = TRUE
    )
  }


# ----------------------------------------------------------------------------------------
# Assign category and significance-streak colors by EYO bin.
# ----------------------------------------------------------------------------------------
  plot_df$category <- group_category_names[max.col(grp_mat, ties.method = "first")]
  plot_df$streak_category <- ifelse(
    plot_df$ci_excludes_zero,
    plot_df$category,
    NA_character_
  )

  # -----------------------------------------------------------------------
  # 4C. Plot
  # -----------------------------------------------------------------------
  y_range <- range(c(plot_df$lower, plot_df$upper, 0), finite = TRUE)

  if (!all(is.finite(y_range))) {
    return(list(
      ok      = FALSE,
      protein = protein_name,
      file    = out_pdf,
      error   = "Non-finite plotting range."
    ))
  }

  y_span <- diff(y_range)
  if (!is.finite(y_span) || y_span == 0) y_span <- 0.1

  y_min <- y_range[1L] - 0.10 * y_span
  y_max <- y_range[2L] + 0.18 * y_span

  streak_y <- y_range[2L] + 0.075 * y_span
  streak_h <- 0.035 * y_span

  streak_df <- plot_df[plot_df$ci_excludes_zero, , drop = FALSE]

  plot_title <- if (!is.na(protein_name) && nzchar(protein_name)) {
    protein_name
  } else {
    clean_name
  }

  y_lab <- if (isTRUE(include_ApoE_Sex_offset_in_curve)) {
    "Full Male.e4 interaction effect beyond additive terms"
  } else {
    "EYO-dependent Male.e4 interaction effect beyond additive terms"
  }


# ----------------------------------------------------------------------------------------
# Draw the single-panel ribbon plot.
# ----------------------------------------------------------------------------------------
  p <- ggplot(plot_df, aes(x = eyo, y = median)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill  = "grey75",
      alpha = 0.70
    ) +
    geom_line(
      color     = "black",
      linewidth = 0.45
    ) +
    geom_hline(
      yintercept = 0,
      linetype   = "dashed",
      linewidth  = 0.35,
      color      = "grey30"
    ) +
    geom_tile(
      data = streak_df,
      aes(
        x    = eyo,
        y    = streak_y,
        fill = streak_category
      ),
      inherit.aes = FALSE,
      width       = eyo_bin_width * 0.95,
      height      = streak_h
    ) +
    scale_fill_manual(
      values = category_cols,
      breaks = names(category_cols),
      drop   = FALSE,
      name   = "Highest group\nwhen 99% CI excludes 0"
    ) +
    coord_cartesian(
      xlim   = range(eyo_grid),
      ylim   = c(y_min, y_max),
      expand = FALSE
    ) +
    labs(
      title    = plot_title,
      subtitle = paste0(
        "Model 2.3 composite of EYO_Spline_Linear:ApoE_Indicator:Sex and ",
        "EYO_Spline_Cubic:ApoE_Indicator:Sex. Ribbon = 99% posterior interval."
      ),
      x = "Estimated Year of Onset (EYO)",
      y = y_lab
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      plot.subtitle   = element_text(size = 9),
      axis.title      = element_text(size = 10),
      axis.text       = element_text(size = 9),
      legend.position = "right",
      legend.title    = element_text(size = 8),
      legend.text     = element_text(size = 8)
    )


# ----------------------------------------------------------------------------------------
# Save the per-protein PDF plot.
# ----------------------------------------------------------------------------------------
  ggsave(
    filename = out_pdf,
    plot     = p,
    device   = grDevices::cairo_pdf,
    width    = 7.5,
    height   = 5.2,
    units    = "in",
    dpi      = 300
  )

  first_sig_idx <- which(plot_df$ci_excludes_zero)[1L]

  first_sig_EYO_bin <- if (length(first_sig_idx) == 0L || is.na(first_sig_idx)) {
    NA_real_
  } else {
    plot_df$eyo[first_sig_idx]
  }

  first_sig_category <- if (length(first_sig_idx) == 0L || is.na(first_sig_idx)) {
    NA_character_
  } else {
    plot_df$category[first_sig_idx]
  }

  list(
    ok                 = TRUE,
    protein            = protein_name,
    clean_name         = clean_name,
    file               = out_pdf,
    first_sig_EYO_bin  = first_sig_EYO_bin,
    first_sig_category = first_sig_category,
    min_posterior_p    = suppressWarnings(min(plot_df$pval, na.rm = TRUE)),
    pval               = plot_df$pval,
    lower              = plot_df$lower,
    median             = plot_df$median,
    upper              = plot_df$upper,
    streak_category    = plot_df$streak_category
  )
}


# -------------------------------------------------------------------------
# 5. Parallel run
# -------------------------------------------------------------------------
ncore <- max(1L, parallel::detectCores() - 2L)


# ----------------------------------------------------------------------------------------
# Start parallel plot generation.
# ----------------------------------------------------------------------------------------
cl <- parallel::makeCluster(ncore)
doParallel::registerDoParallel(cl)

worker_pkgs <- c(
  "tidyverse",
  "Hmisc",
  "rstanarm",
  "ggplot2",
  "grDevices"
)

plot_results <- foreach(
  i              = seq_len(n_prot),
  .packages      = worker_pkgs,
  .export        = c(
    "process_one_model23_plot",
    "find_coef",
    "get_coef_or_zero",
    "safe_file_stem",
    "rds_files",
    "cleaned_from_file",
    "protein_names",
    "eyo_grid",
    "eyo_labels",
    "eyo_bin_width",
    "spline_grid",
    "n_bins",
    "plot_dir",
    "category_cols",
    "group_category_names",
    "include_ApoE_Sex_offset_in_curve"
  ),
  .errorhandling = "pass"
) %dopar% {
  process_one_model23_plot(i)
}

parallel::stopCluster(cl)

# -------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Collect plot outputs and summary rows.
# ----------------------------------------------------------------------------------------
# 6. Collect outputs
# -------------------------------------------------------------------------
ok_idx <- vapply(
  plot_results,
  function(x) is.list(x) && !inherits(x, "error") && isTRUE(x$ok),
  logical(1L)
)

fail_idx <- !ok_idx

plot_summary <- do.call(
  rbind,
  lapply(plot_results[ok_idx], function(x) {
    data.frame(
      Protein            = x$protein,
      CleanName          = x$clean_name,
      First_Sig_EYO_Bin  = x$first_sig_EYO_bin,
      First_Sig_Category = x$first_sig_category,
      Min_PosteriorP     = x$min_posterior_p,
      Plot_File          = x$file,
      stringsAsFactors   = FALSE
    )
  })
)


# ----------------------------------------------------------------------------------------
# Save plot summaries and EYO-by-protein matrices.
# ----------------------------------------------------------------------------------------
write.csv(
  plot_summary,
  file      = file.path(plot_dir, "model23_3way_interaction_plot_summary.csv"),
  row.names = FALSE
)

# Optional matrix outputs for downstream waterfall-style summaries
make_mat <- function(field) {
  mat <- do.call(
    rbind,
    lapply(plot_results[ok_idx], function(x) x[[field]])
  )
  rownames(mat) <- vapply(plot_results[ok_idx], `[[`, character(1L), "protein")
  colnames(mat) <- eyo_labels
  mat
}

pval_mat_plot   <- make_mat("pval")
lower_mat_plot  <- make_mat("lower")
median_mat_plot <- make_mat("median")
upper_mat_plot  <- make_mat("upper")

streak_category_mat <- do.call(
  rbind,
  lapply(plot_results[ok_idx], function(x) x$streak_category)
)

rownames(streak_category_mat) <- vapply(
  plot_results[ok_idx],
  `[[`,
  character(1L),
  "protein"
)

colnames(streak_category_mat) <- eyo_labels

write.csv(
  pval_mat_plot,
  file = file.path(plot_dir, "model23_3way_interaction_posterior_p_by_EYO.csv")
)

write.csv(
  median_mat_plot,
  file = file.path(plot_dir, "model23_3way_interaction_median_by_EYO.csv")
)

write.csv(
  lower_mat_plot,
  file = file.path(plot_dir, "model23_3way_interaction_lower99_by_EYO.csv")
)

write.csv(
  upper_mat_plot,
  file = file.path(plot_dir, "model23_3way_interaction_upper99_by_EYO.csv")
)

write.csv(
  streak_category_mat,
  file = file.path(plot_dir, "model23_3way_interaction_sig_category_streak_by_EYO.csv")
)

saveRDS(
  plot_results,
  file = file.path(plot_dir, "model23_3way_interaction_single_panel_plot_results.rds")
)

if (any(fail_idx)) {
  fail_summary <- do.call(
    rbind,

# ----------------------------------------------------------------------------------------
# Write a failure table if any files could not be plotted.
# ----------------------------------------------------------------------------------------
    lapply(plot_results[fail_idx], function(x) {
      if (inherits(x, "error")) {
        data.frame(
          Protein = NA_character_,
          Error   = conditionMessage(x),
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(
          Protein = x$protein %||% NA_character_,
          Error   = x$error %||% NA_character_,
          stringsAsFactors = FALSE
        )
      }
    })
  )

  write.csv(
    fail_summary,
    file      = file.path(plot_dir, "model23_3way_interaction_plot_failures.csv"),
    row.names = FALSE
  )
}

message("\nFinished Model 2.3 single-panel plots.")
message("Successful plots: ", sum(ok_idx), " / ", length(plot_results))
message("Output directory: ", plot_dir)
