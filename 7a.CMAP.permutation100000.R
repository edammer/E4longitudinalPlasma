##############################################################################
# Pipeline annotation header: 7a.CMAP.permutation100000.R
# Manuscript code section(s): 7
#
# Purpose:
# Compute Connectivity Map-style ontology scores between APOE e4/e4 plasma
# signatures and semaglutide proteomic effects, including permutation-based
# p values and sensitivity sweeps.
#
# Principal inputs:
#   - ALL_heatmap18categoryROWorder-dataFrame.tsv
#   - SemaglutideStudy_TableS6stats.csv
#   - SemaglutideStudy_TableS2stats.csv
#   - 0727_medians_all_assays(5xSD_outliers_excluded).rds
#
# Principal outputs:
#   - CMAP_permutation_stats*.xlsx
#   - Ontology_NCS_withPermutationStars*.pdf
#   - TopHits_2tailed*.csv
#   - minES_sensitivity_sweep*.csv/.rds/.pdf
#   - saved.image-CMAP.perm.RData
#
# Step overview:
#   1. Define weighted running-sum enrichment and weighted connectivity score
#      helper functions.
#   2. Parse ontology gene-symbol lists and semaglutide summary statistics
#      into a SomaScan-aligned ranked vector.
#   3. For each ontology and EYO epoch, form e4-up and e4-down query sets and
#      compute observed C scores.
#   4. Permute same-sized query sets from the ontology universe to estimate
#      normalized connectivity scores, confidence intervals, and two-
#      tailed/directional p values.
#   5. Repeat scoring across minimum rare-effect thresholds and export
#      sensitivity plots and summary tables.
#
# Notes added during manuscript-code alignment:
#   - Original executable code was preserved except for whitespace/encoding cleanup
#     and the explicitly marked non-run archive described below when applicable.
#   - Comments inserted by this pass are prefixed with "ANNOTATION:".
##############################################################################

# ------------------------------------------------------------------------
# ANNOTATION: Define the CMAP-style permutation testing workflow and
# document design choices.
# ------------------------------------------------------------------------
## =============================================================================
## CMAP Connectivity Score - Permutation-Based Statistical Testing
## =============================================================================
##
## PURPOSE
##   For each (ontology x EYO-epoch) cell in the connectivity heatmap, build a
##   null distribution by randomly drawing gene sets of the same size from the
##   ontology's drug-universe assays (drug ranking held fixed), then compute:
##     - Normalized C-score (NCS = obs / mean|null|, the standard CMAP metric)
##     - Two-tailed permutation p-value   (H0: |score| ? |obs| by chance)
##     - Directional permutation p-value  (H0: score is as extreme in observed direction)
##     - 95 % CI of the null distribution
##     - Null mean and SD
##
## DESIGN CHOICE - what to permute
##   We permute the QUERY (which assays are called "up" vs "down" from the rare-
##   genotype data), holding the DRUG ranking fixed.  This is identical to the
##   approach used in the original CMAP / L1000 papers and tests whether the
##   *directionality* of the rare-genotype signature matches the drug signature
##   beyond what would be expected for a random same-sized gene set from the same
##   ontology background.
##
## INPUTS (identical objects already in the workspace)
##   onts141.symbolList   - named list of gene-symbol vectors (141 ontologies)
##   drug_ranked_stats    - named numeric vector of semaglutide T-stats
##   effectSizes          - 143 x 7333 matrix of e4/e4 effect sizes (log2 abun)
##   epoch_list           - named list of row-index vectors defining EYO epochs
##
## OUTPUTS
##   A named list with one element per statistic (each a matrix ont x epoch):
##     $Cscore, $NCS, $p_two_tailed, $p_directional,
##     $CI_lower, $CI_upper, $null_mean, $null_sd,
##     $n_up, $n_down, $n_drug_assays
##
## DEPENDENCIES
##   calc_es(), calc_wtcs() are defined just below.
##   parallel, doParallel, foreach

# ------------------------------------------------------------------------
# ANNOTATION: Weighted running-sum enrichment helper for a single up/down
# query set.
# ------------------------------------------------------------------------
## =============================================================================


## Score calculation algorithm
calc_es <- function(ranked_stats, geneset, p = 1) {
  ranked_stats <- ranked_stats[order(ranked_stats, decreasing = TRUE)]
  hits <- names(ranked_stats) %in% geneset
  Nh <- sum(hits)
  N <- length(ranked_stats)
  if (Nh == 0L || Nh == N) return(0)

  w <- abs(ranked_stats)^p
  P_hit  <- cumsum(hits * w) / sum(w[hits])
  P_miss <- cumsum(!hits) / (N - Nh)
  running <- P_hit - P_miss


# ------------------------------------------------------------------------
# ANNOTATION: Weighted connectivity score helper combining up and down
# enrichment components.
# ------------------------------------------------------------------------
  mx <- max(running); mn <- min(running)
  if (abs(mx) >= abs(mn)) mx else mn
}

calc_wtcs <- function(drug_ranked_stats, up_set, down_set, p = 1) {
  es_up   <- calc_es(drug_ranked_stats, up_set,   p = p)
  es_down <- calc_es(drug_ranked_stats, down_set, p = p)

  if (sign(es_up) == sign(es_down)) return(0)
  (es_up - es_down) / 2
}


## ---------------------------------------------------------------------------
## 0.  Helper: permutation null for ONE (ontology x epoch) cell
## ---------------------------------------------------------------------------


# ------------------------------------------------------------------------
# ANNOTATION: Permutation helper for one ontology x EYO-epoch cell.
# ------------------------------------------------------------------------
.permute_cscore_cell <- function(drug_ranked_final,   # named numeric; already 1-per-gene, ontology-filtered
                                 up_used,             # observed up assay names (subset of drug_ranked_final)
                                 down_used,           # observed down assay names
                                 obs_cscore,          # pre-computed observed C-score
                                 n_perm  = 1000,
                                 p       = 1) {

  n_up   <- length(up_used)
  n_down <- length(down_used)
  n_all  <- length(drug_ranked_final)
  all_assays <- names(drug_ranked_final)

  ## Edge cases: if both sets empty, null is trivially 0
  if ((n_up + n_down) == 0L || n_all == 0L) {
    return(list(
      null_scores   = rep(0, n_perm),
      p_two_tailed  = NA_real_,
      p_directional = NA_real_,
      ci_lower      = 0, ci_upper = 0,
      null_mean     = 0, null_sd  = 0,
      ncs           = NA_real_
    ))
  }

  ## If query is larger than universe (shouldn't happen but safe)
  draw_n <- min(n_up + n_down, n_all)
  n_up_draw   <- min(n_up,   draw_n)
  n_down_draw <- min(n_down, draw_n - n_up_draw)

  null_scores <- numeric(n_perm)

  for (i in seq_len(n_perm)) {
    samp      <- sample(all_assays, draw_n, replace = FALSE)
    perm_up   <- if (n_up_draw   > 0) samp[seq_len(n_up_draw)]   else character(0)
    perm_down <- if (n_down_draw > 0) samp[seq(n_up_draw + 1L,
                                                n_up_draw + n_down_draw)] else character(0)

    up_comp   <- calc_es(drug_ranked_final, perm_up,   p = p)
    down_comp <- calc_es(drug_ranked_final, perm_down, p = p)

    null_scores[i] <- if (sign(up_comp) == sign(down_comp)) {
      0
    } else {
      (up_comp - down_comp) / 2
    }
  }

  ## -- statistics --
  p_two_tailed  <- mean(abs(null_scores) >= abs(obs_cscore) - .Machine$double.eps)
  p_directional <- if (obs_cscore >= 0) {
    mean(null_scores >= obs_cscore - .Machine$double.eps)
  } else {
    mean(null_scores <= obs_cscore + .Machine$double.eps)
  }

  ci            <- quantile(null_scores, c(0.025, 0.975), na.rm = TRUE)
  null_mean     <- mean(null_scores, na.rm = TRUE)
  null_sd       <- sd(null_scores,   na.rm = TRUE)
  mean_abs_null <- mean(abs(null_scores), na.rm = TRUE)
  ncs           <- if (is.finite(mean_abs_null) && mean_abs_null > 0) {
    obs_cscore / mean_abs_null
  } else {
    NA_real_
  }

  list(
    null_scores   = null_scores,
    p_two_tailed  = p_two_tailed,
    p_directional = p_directional,
    ci_lower      = as.numeric(ci[1]),
    ci_upper      = as.numeric(ci[2]),
    null_mean     = null_mean,
    null_sd       = null_sd,
    ncs           = ncs
  )
}


## ---------------------------------------------------------------------------
## 1.  Core function: per-ontology worker (returns a list of epoch results)
## ---------------------------------------------------------------------------
## This is exported to each parallel worker; the outer loop is over ontologies.

.ontology_permutation_worker <- function(ont,                   # ontology name
                                         this.ontology,         # character vector of gene symbols
                                         drug_ranked_stats,     # full named numeric vector
                                         epoch_es_list,         # named list: epoch_name -> ES sub-matrix (rows=tp, cols=assays)
                                         minES.rare = 0.1,
                                         n_perm     = 1000,
                                         p          = 1) {

  assay_symbol <- function(x) sub("\\|.*$", "", x)

  ## -- build ontology-restricted drug universe (1 assay per gene) --
  drug_assays_all <- names(drug_ranked_stats)
  drug_syms_all   <- assay_symbol(drug_assays_all)

# ------------------------------------------------------------------------
# ANNOTATION: Main permutation driver over ontologies and EYO epochs.
# ------------------------------------------------------------------------

  in_ont <- (drug_syms_all %in% this.ontology) &
            is.finite(drug_ranked_stats) & !is.na(drug_ranked_stats)

  ## Empty overlap ? all NAs for every epoch
  if (!any(in_ont)) {
    empty <- list(
      Cscore = NA_real_, NCS = NA_real_,
      p_two_tailed = NA_real_, p_directional = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_,
      null_mean = NA_real_, null_sd = NA_real_,
      n_up = 0L, n_down = 0L, n_drug_assays = 0L
    )
    return(setNames(replicate(length(epoch_es_list), empty, simplify = FALSE),
                    names(epoch_es_list)))
  }

  ## Collapse to 1 assay per gene by max |T|
  df_drug <- data.frame(
    assay = drug_assays_all[in_ont],
    sym   = drug_syms_all[in_ont],
    Tstat = as.numeric(drug_ranked_stats[in_ont]),
    stringsAsFactors = FALSE
  )
  df_drug <- df_drug[order(df_drug$sym, -abs(df_drug$Tstat), df_drug$assay), ]
  keep_assays       <- df_drug$assay[!duplicated(df_drug$sym)]
  drug_ranked_final <- drug_ranked_stats[keep_assays]
  n_drug_assays     <- length(drug_ranked_final)

  ## Lookup: gene symbol -> kept drug assay name
  drug_assay_by_gene <- setNames(names(drug_ranked_final),
                                 assay_symbol(names(drug_ranked_final)))

  ## -- process each epoch --
  results <- vector("list", length(epoch_es_list))
  names(results) <- names(epoch_es_list)

  for (ep_name in names(epoch_es_list)) {

    es_sub <- epoch_es_list[[ep_name]]   # matrix: rows = time-points, cols = assays

    ## Find rare-assay columns that are in this ontology
    rare_assays_all <- colnames(es_sub)
    rare_syms_all   <- assay_symbol(rare_assays_all)
    rare_in_ont     <- rare_syms_all %in% this.ontology

    if (!any(rare_in_ont)) {
      results[[ep_name]] <- list(
        Cscore = 0, NCS = NA_real_,
        p_two_tailed = NA_real_, p_directional = NA_real_,
        ci_lower = NA_real_, ci_upper = NA_real_,
        null_mean = NA_real_, null_sd = NA_real_,
        n_up = 0L, n_down = 0L, n_drug_assays = n_drug_assays
      )
      next
    }

    rare_mat_ont    <- es_sub[, rare_in_ont, drop = FALSE]
    rare_assays_ont <- colnames(rare_mat_ont)
    rare_syms_ont   <- rare_syms_all[rare_in_ont]

    ## Peak signed ES per assay (over all timepoints in this epoch)
    rare_peakES <- apply(rare_mat_ont, 2, function(v) {
      v <- as.numeric(v)
      if (all(!is.finite(v))) return(NA_real_)
      v[which.max(abs(v))]
    })
    rare_peakAbs <- abs(rare_peakES)

    keep_by_min <- is.finite(rare_peakES) & (rare_peakAbs >= minES.rare)

    if (!any(keep_by_min)) {
      results[[ep_name]] <- list(
        Cscore = 0, NCS = NA_real_,
        p_two_tailed = NA_real_, p_directional = NA_real_,
        ci_lower = NA_real_, ci_upper = NA_real_,
        null_mean = NA_real_, null_sd = NA_real_,
        n_up = 0L, n_down = 0L, n_drug_assays = n_drug_assays
      )
      next
    }

    cand_df <- data.frame(
      assay  = rare_assays_ont[keep_by_min],
      sym    = rare_syms_ont[keep_by_min],
      peakES = as.numeric(rare_peakES[keep_by_min]),
      peakAbs = as.numeric(rare_peakAbs[keep_by_min]),
      stringsAsFactors = FALSE
    )

    ## One rare assay per gene: prefer matching drug assay name, else max |peakES|
    chosen <- lapply(split(cand_df, cand_df$sym), function(df) {
      g         <- df$sym[1]
      preferred <- drug_assay_by_gene[[g]]
      if (!is.null(preferred) && preferred %in% df$assay) {
        df[df$assay == preferred, , drop = FALSE][1, ]
      } else {
        df[order(-df$peakAbs, df$assay), , drop = FALSE][1, ]
      }
    })
    chosen_df <- do.call(rbind, chosen)

    ## Intersect with drug_ranked_final so only scored assays are used
    up_raw   <- chosen_df$assay[chosen_df$peakES > 0]
    down_raw <- chosen_df$assay[chosen_df$peakES < 0]
    up_used  <- intersect(up_raw,   names(drug_ranked_final))
    down_used <- intersect(down_raw, names(drug_ranked_final))

    ## Observed C-score
    up_comp   <- calc_es(drug_ranked_final, up_used,   p = p)
    down_comp <- calc_es(drug_ranked_final, down_used, p = p)
    obs_cscore <- if (sign(up_comp) == sign(down_comp)) {
      0
    } else {
      (up_comp - down_comp) / 2
    }

    ## Permutation null
    perm_res <- .permute_cscore_cell(
      drug_ranked_final = drug_ranked_final,
      up_used           = up_used,
      down_used         = down_used,
      obs_cscore        = obs_cscore,
      n_perm            = n_perm,
      p                 = p
    )

    results[[ep_name]] <- list(
      Cscore        = obs_cscore,
      NCS           = perm_res$ncs,
      p_two_tailed  = perm_res$p_two_tailed,
      p_directional = perm_res$p_directional,
      ci_lower      = perm_res$ci_lower,
      ci_upper      = perm_res$ci_upper,
      null_mean     = perm_res$null_mean,
      null_sd       = perm_res$null_sd,
      n_up          = length(up_used),
      n_down        = length(down_used),
      n_drug_assays = n_drug_assays
    )
  }

  results
}


## ---------------------------------------------------------------------------
## 2.  Main exported function
## ---------------------------------------------------------------------------

build_Cscore_permutation_stats <- function(
    onts141.symbolList,   # named list of character vectors (gene symbols per ontology)
    drug_ranked_stats,    # named numeric: semaglutide T-stats, names = assay names
    effectSizes,          # full 143 x 7333 ES matrix (rownames = EYO, colnames = assay)
    epoch_list,           # named list of integer row-index vectors
                          #   e.g. list("EYO_pre50" = 1:62, "peri_onset" = 63:92, "post_onset" = 93:143)
    minES.rare = 0.1,     # minimum |ES| threshold for rare-genotype assay inclusion
    n_perm     = 1000,    # permutations per cell  (use >=1000 for publication; 200 for QC)
    p          = 1,       # weighting exponent in calc_es (1 = weighted, 0 = unweighted)
    ncores     = max(1L, parallel::detectCores() - 1L),
    seed       = 42       # global RNG seed (each worker gets a deterministic offset)
) {

  if (!is.list(onts141.symbolList) || is.null(names(onts141.symbolList)))
    stop("'onts141.symbolList' must be a *named* list.")
  if (!is.list(epoch_list) || is.null(names(epoch_list)))
    stop("'epoch_list' must be a *named* list of row-index vectors.")

  ont_names <- names(onts141.symbolList)
  ep_names  <- names(epoch_list)
  nO        <- length(ont_names)
  nE        <- length(ep_names)

  ## Pre-slice ES matrix into sub-matrices for each epoch (avoids shipping the
  ## full 143-row matrix to every worker on every iteration)
  es_full <- as.matrix(effectSizes)
  storage.mode(es_full) <- "numeric"

  epoch_es_list <- lapply(epoch_list, function(idx) {
    sub <- es_full[idx, , drop = FALSE]
    storage.mode(sub) <- "numeric"
    sub
  })

  ## Stat containers
  stat_names <- c("Cscore","NCS","p_two_tailed","p_directional",
                  "ci_lower","ci_upper","null_mean","null_sd",
                  "n_up","n_down","n_drug_assays")
  mats <- lapply(stat_names, function(s) {
    m <- matrix(NA_real_, nrow = nO, ncol = nE,
                dimnames = list(ont_names, ep_names))
    m
  })
  names(mats) <- stat_names

  ## ---- Parallel loop over ontologies ----------------------------------------
  library(parallel)
  library(doParallel)
  library(foreach)

  cl <- parallel::makeCluster(ncores, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  on.exit({
    foreach::registerDoSEQ()
    parallel::stopCluster(cl)
  }, add = TRUE)

  ## Export shared objects to workers once
  parallel::clusterExport(
    cl,
    varlist = c("drug_ranked_stats", "epoch_es_list", "minES.rare",
                "n_perm", "p", "seed",
                "calc_es", "calc_wtcs",
                ".permute_cscore_cell", ".ontology_permutation_worker"),
    envir = environment()
  )

  message(sprintf(
    "[%s]  Running %d permutations x %d ontologies x %d epochs on %d cores ...",
    format(Sys.time(), "%H:%M:%S"), n_perm, nO, nE, ncores
  ))

  ## Parallelise over ontologies; each worker returns a named list (epochs ? stats)
  all_results <- foreach(
    ont_idx = seq_len(nO),
    .inorder  = TRUE,
    .packages = character(0)
  ) %dopar% {

    ## Per-worker deterministic seed so results are reproducible
    set.seed(seed + ont_idx)

    ont      <- ont_names[ont_idx]
    syms     <- unique(trimws(as.character(onts141.symbolList[[ont]])))
    syms     <- syms[!is.na(syms) & nzchar(syms)]

    .ontology_permutation_worker(
      ont               = ont,
      this.ontology     = syms,
      drug_ranked_stats = drug_ranked_stats,
      epoch_es_list     = epoch_es_list,
      minES.rare        = minES.rare,
      n_perm            = n_perm,
      p                 = p
    )
  }
  names(all_results) <- ont_names

  message(sprintf("[%s]  Permutations complete. Assembling output matrices ...",
                  format(Sys.time(), "%H:%M:%S")))

  ## ---- Assemble result matrices --------------------------------------------
  for (ont_idx in seq_len(nO)) {
    ont <- ont_names[ont_idx]
    for (ep in ep_names) {
      cell <- all_results[[ont]][[ep]]
      for (s in stat_names) {
        val <- cell[[s]]
        mats[[s]][ont_idx, ep] <- if (!is.null(val) && length(val) == 1L) {
          as.numeric(val)
        } else {
          NA_real_
        }
      }
    }
  }

  mats
}


## ---------------------------------------------------------------------------
## 3.  Convenience: FDR-adjust p-values across all cells (Benjamini-Hochberg)
## ---------------------------------------------------------------------------

add_FDR_to_permutation_stats <- function(perm_stats,
                                         p_col = "p_two_tailed") {
  p_mat <- perm_stats[[p_col]]
  fdr_mat <- matrix(
    p.adjust(as.vector(p_mat), method = "BH"),
    nrow = nrow(p_mat), ncol = ncol(p_mat),
    dimnames = dimnames(p_mat)
  )
  perm_stats[[paste0("FDR_", p_col)]] <- fdr_mat

  ## Also FDR-adjust directional p
  pd_mat <- perm_stats[["p_directional"]]
  fdr_pd <- matrix(
    p.adjust(as.vector(pd_mat), method = "BH"),
    nrow = nrow(pd_mat), ncol = ncol(pd_mat),
    dimnames = dimnames(pd_mat)
  )
  perm_stats[["FDR_p_directional"]] <- fdr_pd

  perm_stats
}


## ---------------------------------------------------------------------------
## 4.  Write all stat matrices to an Excel workbook (one sheet per stat)
## ---------------------------------------------------------------------------

write_permutation_stats_excel <- function(perm_stats,
                                          file = "CMAP_permutation_stats.xlsx",
                                          overwrite = TRUE) {
  library(openxlsx)

  wb <- createWorkbook()

  for (nm in names(perm_stats)) {
    mat <- perm_stats[[nm]]
    df  <- as.data.frame(mat, check.names = FALSE)
    df  <- cbind(data.frame(ontology = rownames(df), stringsAsFactors = FALSE), df)

    sheet_nm <- substr(nm, 1, 31)
    addWorksheet(wb, sheet_nm)
    writeData(wb, sheet = sheet_nm, x = df, withFilter = TRUE)
    freezePane(wb, sheet = sheet_nm, firstRow = TRUE)
  }

  saveWorkbook(wb, file = file, overwrite = overwrite)
  invisible(file)
}


## ---------------------------------------------------------------------------
## 5.  Heatmap overlay: significance asterisks on existing ComplexHeatmap
## ---------------------------------------------------------------------------
##
## Call this inside (or after) your draw() call to overlay *, **, *** on the
## heatmap body wherever the permutation p-value passes the threshold.
##
## Example:
##   decorate_heatmap_body("Connectivity", {
##     overlay_significance_stars(perm_stats$p_two_tailed, ont_order, ep_order)
##   })

overlay_significance_stars <- function(p_mat,
                                       row_order,    # character vector matching rownames(p_mat)
                                       col_order,    # character vector matching colnames(p_mat)
                                       thresholds = c("***" = 0.001,
                                                      "**"  = 0.01,
                                                      "*"   = 0.05),
                                       fontsize = 7) {
  p_sub <- p_mat[row_order, col_order, drop = FALSE]
  nr    <- nrow(p_sub)
  nc    <- ncol(p_sub)

  thresholds <- sort(thresholds)   # ascending so tightest is last

  for (ci in seq_len(nc)) {
    for (ri in seq_len(nr)) {
      pv  <- p_sub[ri, ci]
      if (!is.finite(pv)) next

      star <- ""
      for (thr_nm in names(thresholds)) {
        if (pv < thresholds[thr_nm]) star <- thr_nm
      }
      if (!nzchar(star)) next

      x_npc <- (ci - 0.5) / nc
      y_npc <- 1 - (ri - 0.5) / nr

      grid::grid.text(
        star,
        x  = grid::unit(x_npc, "npc"),
        y  = grid::unit(y_npc, "npc"),
        gp = grid::gpar(fontsize = fontsize, fontface = "bold", col = "black")
      )
    }
  }
}


## =============================================================================
## USAGE EXAMPLE (un-comment and run after sourcing this file)
## =============================================================================

setwd("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/#manuscript/SciTranslMed_plan/CMAP.perm/")


#####################################################################
## 141 full ontologies - Split Sema S6 (or S2) ranked 7288 assay list to keep only ontology gene-linked assays.
##                       Then extract up and down lists based solely on maximum absolute effect size within window before EYO -15 (5 year windows 1:27),
##                       and calculate a connectivity score for each, saving both the up and down component scores to a data frame, as well.
##                       Also include in the table the count of up and down assays, and the assay names for up and down used.
##
##                     - cleaned all assay names in sema tables S6 and S2 to match our assay names.


# ------------------------------------------------------------------------
# ANNOTATION: Load curated 141-ontology list/order table.
# ------------------------------------------------------------------------
onts141.df<-read.table(file="F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/ALL_heatmap18categoryROWorder-dataFrame.tsv",sep="\t",header=TRUE,quote="")

GMTdatabaseFile="f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/Human_GO_AllPathways_noPFOCR_with_GO_iea_June_01_2025_symbol.gmt"

	## Load GMT file; Clean UTF-8 characters (since Dec 2023); Write clean.GMT back out
	#GMT.df <- read.delim(GMTdatabaseFile, encoding = "utf-8",quote="", sep="\t",header=FALSE)
	GMT.df <- readLines(con <- file(GMTdatabaseFile, encoding = "utf-8"))
        close(con)
        GMT.df <- unlist(sapply(GMT.df, function(x) iconv(gsub("^(PMC\\d*__.+?)\\\t(.*)$","\\1%PMC%\\2",
                                                          gsub("\\\"","", gsub("\\x83\\x80.","-",x) )),
                                                          "utf-8","ASCII", "")))
	names(GMT.df)<-NULL
        GMT.df <- lapply(GMT.df, function(x) stringr::str_split_fixed(x, pattern="\t", n=Inf))

        # Create list object that is identical to a GSC class object, just not of this class, since not loaded by the loadGSC() function in piano package.
        GSCfromGMT<-list()
        GSCfromGMT[["addInfo"]]<-do.call(rbind, lapply(GMT.df, function(x) if(grepl("^PMC.*\\%PMC\\%",x[1])) { c(x[1],gsub("^(PMC.*)\\%PMC\\%.*$","\\1",x[1])) } else { x[c(1:2)] } ))
        GSCfromGMT[["gsc"]]<-lapply(GMT.df, function(x) if(grepl("^PMC.*\\%PMC\\%",x[1])) { x[c(2:length(x))][!x[c(2:length(x))]==""] } else { x[c(3:length(x))][!x[c(3:length(x))]==""] })
        names(GSCfromGMT$gsc)<-GSCfromGMT$addInfo[,1]

        ontology=stringr::str_to_title(gsub("\\%WIKIPATHWAYS_\\d*","", gsub("\\%WP_\\d*","", gsub("\\&(.*);","\\1",gsub("<\\sI>","",gsub("<I>","", gsub("(.*)\\%.*\\%.*","\\1",names(GSCfromGMT$gsc))))))))
	ontologyType=gsub("^WP\\d*","WikiPathways", gsub(".*\\%(.*)\\%.*","\\1",names(GSCfromGMT$gsc)))

        #force all caps for ontologyType of GObp GOmf GOcc (changed in downloaded GMT files Sept 2022 and/or different in mouse GMT compared to human)
        ontologyType=gsub("GObp","GOBP",ontologyType)
        ontologyType=gsub("GOmf","GOMF",ontologyType)
        ontologyType=gsub("GOcc","GOCC",ontologyType)


length(which(onts141.df$ontology %in% ontology))
#[1] 141
length(which(ontology %in% onts141.df$ontology))
#>141

onts141.df$ontology.ontType<-paste0(onts141.df$ontology," (",onts141.df$ontologyType,")")
ontology.ontType<-paste0(ontology," (",ontologyType,")")

length(which(ontology.ontType %in% onts141.df$ontology.ontType))
#141

onts141.symbolList<-list()
for (this.ont in onts141.df$ontology.ontType) onts141.symbolList[[this.ont]] <- GSCfromGMT$gsc[[which(ontology.ontType==as.character(this.ont))]]


## get drug ranked stats (7288 assays)
sema<-read.csv(file="SemaglutideStudy_TableS6stats.csv",header=TRUE)

# ------------------------------------------------------------------------
# ANNOTATION: Read semaglutide adjusted Table S6 statistics and create
# assay-level ranked T statistics.
# ------------------------------------------------------------------------
#sema<-read.csv(file="SemaglutideStudy_TableS2stats.csv",header=TRUE)
# sema: data.frame with columns Assay.Name, effect, se (rename as needed)
sema$stat <- with(sema, effect_size / standard_error)

# clean
sema <- sema[is.finite(sema$stat) & !is.na(sema$Assay.Name) & sema$Assay.Name != "", ]

# if duplicate assay names exist, keep the entry with max |stat|
sema <- sema[order(abs(sema$stat), decreasing = TRUE), ]
sema <- sema[!duplicated(sema$Assay.Name), ]

drug_ranked_stats <- setNames(sema$stat, sema$Assay.Name)
# no BMI or A1C effect left in


## get drug ranked stats w/o BMI or A1c adjustment (7288 assays)
semaS2<-read.csv(file="SemaglutideStudy_TableS2stats.csv",header=TRUE)
# semaS2: data.frame with columns Assay.Name, effect, se (rename as needed)
semaS2$stat <- with(semaS2, effect_size / standard_error)

# clean
semaS2 <- semaS2[is.finite(semaS2$stat) & !is.na(semaS2$Assay.Name) & semaS2$Assay.Name != "", ]

# if duplicate assay names exist, keep the entry with max |stat|
semaS2 <- semaS2[order(abs(semaS2$stat), decreasing = TRUE), ]
semaS2 <- semaS2[!duplicated(semaS2$Assay.Name), ]

drug_ranked_stats.S2 <- setNames(semaS2$stat, sema$Assay.Name)
# BMI and A1C effects left in


#"f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/5.CMAPalgo/
#pVals<-readRDS("../simple.3177/scatter/_99_par_diff_all_peptide_p_value.rds")
effectSizes<-readRDS("f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/0727_medians_all_assays(5xSD_outliers_excluded).rds")
#stdErrors<-  # have not extracted this from the models yet

dim(effectSizes)

# ------------------------------------------------------------------------
# ANNOTATION: Read APOE e4/e4 effect-size trajectories across EYO.
# ------------------------------------------------------------------------
#[1]  143 7333


#if (FALSE) {

  ## ---- Define the three EYO epochs ----------------------------------------
  epoch_list <- list(
    "EYO -46 to -15.5" = 1:62,     # pre-clinical, ~50 yr before diagnosis
    "EYO -15 to -0.5"  = 63:92,    # peri-onset
    "EYO 0 to +25"     = 93:143    # post-onset
  )

  ## ---- Run permutations (S6: BMI/A1C adjusted) ----------------------------
  ##  n_perm = 200  ? quick QC run (~minutes on 8 cores)
  ##  n_perm = 1000 ? publication quality (use 31 cores as in original code)

  perm_stats_S6 <- build_Cscore_permutation_stats(
    onts141.symbolList = onts141.symbolList,
    drug_ranked_stats  = drug_ranked_stats,   # S6: BMI/A1C adjusted
    effectSizes        = effectSizes,
    epoch_list         = epoch_list,
    minES.rare         = 0.03,
    n_perm             = 100000,
    ncores             = 31,
    seed               = 42
  )

  ## ---- Add FDR columns (BH across all ontologyxepoch cells) --------------
  perm_stats_S6 <- add_FDR_to_permutation_stats(perm_stats_S6)

  ## ---- Save all stat matrices to Excel ------------------------------------
  write_permutation_stats_excel(
    perm_stats_S6,
    file = "CMAP_permutation_stats_S6_100000perm_minES.e4_0.03.xlsx"
  )
  ## Sheets: Cscore, NCS, p_two_tailed, p_directional,
  ##         ci_lower, ci_upper, null_mean, null_sd,

# ------------------------------------------------------------------------
# ANNOTATION: Run and export 100,000-permutation statistics for adjusted
# semaglutide Table S6.
# ------------------------------------------------------------------------
  ##         n_up, n_down, n_drug_assays,
  ##         FDR_p_two_tailed, FDR_p_directional

  ## ---- Repeat for S2 (unadjusted, BMI/A1C effects retained) --------------
  perm_stats_S2 <- build_Cscore_permutation_stats(
    onts141.symbolList = onts141.symbolList,
    drug_ranked_stats  = drug_ranked_stats.S2,
    effectSizes        = effectSizes,
    epoch_list         = epoch_list,
    minES.rare         = 0.1,
    n_perm             = 1000,
    ncores             = 31,
    seed               = 42
  )
  perm_stats_S2 <- add_FDR_to_permutation_stats(perm_stats_S2)
  write_permutation_stats_excel(
    perm_stats_S2,
    file = "CMAP_permutation_stats_S2_1000perm.xlsx"
  )


  ## ---- Reproduce the existing heatmap with significance stars -------------

  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(RColorBrewer)

  ## Use the NCS matrix (normalised C-score) rather than raw C-score for display
  ## - NCS accounts for how variable the null is for each ontology (small gene
  ## sets have noisier nulls; NCS corrects for that).

  ncs_mat  <- perm_stats_S6$NCS           # 141 x 3
  fdr_mat  <- perm_stats_S6$p_two_tailed  #perm_stats_S6$FDR_p_two_tailed

  ## Reorder to match the original Fig 7A row/column order
  ## (assuming ont_order and ep_order match the input file)
  ont_order <- rownames(ncs_mat)           # already in symbolList order
  ep_order  <- colnames(ncs_mat)

  col_ncs <- colorRamp2(c(-3, 0, 3), c("#0099FF", "white", "gold"))

  ht_ncs <- Heatmap(
    ncs_mat,
    name             = "NCS",
    col              = col_ncs,
    cluster_rows     = FALSE,
    cluster_columns  = FALSE,
    show_row_names   = TRUE,
    row_names_side   = "left",
    row_names_gp     = gpar(fontsize = 7),
    row_names_max_width=unit(6.5, "in"),
    width = unit(0.10, "npc"),  #10% of page width for heatmap body
    show_column_names = TRUE,
    column_names_gp  = gpar(fontsize = 8),
    border           = TRUE,
    heatmap_legend_param = list(
      title  = "NCS",
      at     = c(-3, -1, 0, 1, 3),
      labels = c("-3", "-1", "0", "1", "3")
    )
  )

  pdf("Ontology_NCS_withPermutationStars_Fig7A-100000_perm&minES_0.03.pdf",
      width = 8.5, height = 17, useDingbats = FALSE)

  draw(ht_ncs, heatmap_legend_side = "right")

# ------------------------------------------------------------------------
# ANNOTATION: Draw ontology NCS heatmaps with permutation significance
# glyphs.
# ------------------------------------------------------------------------

  ## Overlay significance stars (BH-FDR threshold 0.05)
  decorate_heatmap_body("NCS", {
    overlay_significance_stars(
      p_mat     = fdr_mat,
      row_order = ont_order,
      col_order = ep_order,
      thresholds = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
      fontsize  = 7
    )
  })

  dev.off()


  ## ---- Quick sanity check: look at top hits ---------------------------------
  top_hits <- which(perm_stats_S6$p_two_tailed < 0.05, arr.ind = TRUE)  # was $FDR_p_two_tailed < 0.05
  top_df <- data.frame(
    ontology     = rownames(perm_stats_S6$Cscore)[top_hits[, 1]],
    epoch        = colnames(perm_stats_S6$Cscore)[top_hits[, 2]],
    Cscore       = perm_stats_S6$Cscore[top_hits],
    NCS          = perm_stats_S6$NCS[top_hits],
    p_two_tailed = perm_stats_S6$p_two_tailed[top_hits],
    FDR          = perm_stats_S6$FDR_p_two_tailed[top_hits],
    n_up         = perm_stats_S6$n_up[top_hits],
    n_down       = perm_stats_S6$n_down[top_hits],
    ci_lower     = perm_stats_S6$ci_lower[top_hits],
    ci_upper     = perm_stats_S6$ci_upper[top_hits],
    stringsAsFactors = FALSE
  )
  top_df <- top_df[order(top_df$p_two_tailed), ]  #was $FDR
  print(head(top_df, 20))
  write.csv(top_df, "TopHits_2tailedPlt0.05_100000perm_minES_0.03.csv", row.names = FALSE)

#}
## =============================================================================
## END
## =============================================================================


## ---- sweep parameters -------------------------------------------------------
minES_grid <- seq(0.01, 0.20, by = 0.02)   # 10 values: 0.01 0.03 0.05 ... 0.19
sig_alpha  <- 0.05                          # p_two_tailed threshold
n_perm     <- 100000                          # permutations per cell
ncores_use <- max(1L, parallel::detectCores() - 1L)
sweep_seed <- 42

## ---- epoch definition (edit to match your analysis) ------------------------
epoch_list <- list(
  "EYO -46 to -15.5" = 1:62,
  "EYO -15 to -0.5"  = 63:92,
  "EYO 0 to +25"     = 93:143
)

## ---- result collector -------------------------------------------------------
sweep_df <- data.frame(
  minES_rare    = minES_grid,
  n_sig_neg_NCS = NA_integer_,   # significant cells with NCS < 0  (opposing drug)
  n_sig_pos_NCS = NA_integer_,   # significant cells with NCS > 0  (concordant)
  stringsAsFactors = FALSE
)

## Optional: store each full result if RAM allows (set to FALSE to save memory)
store_full_results <- TRUE
sweep_results_list <- if (store_full_results) vector("list", length(minES_grid)) else NULL
names(sweep_results_list) <- as.character(round(minES_grid, 4))

## ---- sweep loop -------------------------------------------------------------
message(sprintf("[%s]  Starting minES.rare sensitivity sweep (%d values) ...",
                format(Sys.time(), "%H:%M:%S"), length(minES_grid)))

for (k in seq_along(minES_grid)) {

  this_minES <- minES_grid[k]
  message(sprintf("[%s]  >> minES.rare = %.2f  (%d / %d)",
                  format(Sys.time(), "%H:%M:%S"), this_minES, k, length(minES_grid)))

  perm_stats <- build_Cscore_permutation_stats(
    onts141.symbolList = onts141.symbolList,
    drug_ranked_stats  = drug_ranked_stats,
    effectSizes        = effectSizes,
    epoch_list         = epoch_list,
    minES.rare         = this_minES,
    n_perm             = n_perm,
    ncores             = ncores_use,
    seed               = sweep_seed
  )

  ## Apply BH-FDR across all cells (consistent with single-run usage)

# ------------------------------------------------------------------------
# ANNOTATION: Run minimum-effect-size sensitivity sweeps for APOE rare-
# signature thresholds.
# ------------------------------------------------------------------------
  perm_stats <- add_FDR_to_permutation_stats(perm_stats)

  ## ---- count significant cells by NCS sign --------------------------------
  ## Use raw p_two_tailed here; swap for FDR_p_two_tailed if preferred.
  ## Each "cell" is one (ontology x epoch) combination.
  p_mat   <- perm_stats$p_two_tailed   # 141 x 3
  ncs_mat <- perm_stats$NCS            # 141 x 3

  is_sig     <- is.finite(p_mat)   & (p_mat   <= sig_alpha)
  is_neg_ncs <- is.finite(ncs_mat) & (ncs_mat <  0)
  is_pos_ncs <- is.finite(ncs_mat) & (ncs_mat >  0)

  sweep_df$n_sig_neg_NCS[k] <- sum(is_sig & is_neg_ncs)
  sweep_df$n_sig_pos_NCS[k] <- sum(is_sig & is_pos_ncs)

  if (store_full_results) {
    sweep_results_list[[k]] <- perm_stats
  }

  message(sprintf("         sig & NCS<0 = %d  |  sig & NCS>0 = %d",
                  sweep_df$n_sig_neg_NCS[k],
                  sweep_df$n_sig_pos_NCS[k]))
}

message(sprintf("[%s]  Sweep complete.", format(Sys.time(), "%H:%M:%S")))

## ---- inspect summary table --------------------------------------------------
print(sweep_df)
write.csv(sweep_df, "minES_sensitivity_sweep_summary.csv", row.names = FALSE)

if (store_full_results) {
  saveRDS(sweep_results_list, "minES_sensitivity_sweep_full_results.rds")
}


## =============================================================================
## PLOT: dual-trace line plot (base R + optional ggplot2 version)
## =============================================================================

## ---- colour palette (colour-blind friendly) ---------------------------------
col_neg <- "#0072B2"   # blue  ? NCS < 0 (drug opposes rare-genotype trend)
col_pos <- "#D55E00"   # vermillion ? NCS > 0 (drug concordant with rare trend)

## ---- axis limits ------------------------------------------------------------
y_max <- max(c(sweep_df$n_sig_neg_NCS, sweep_df$n_sig_pos_NCS), na.rm = TRUE)
y_max <- ceiling(y_max * 1.15)   # 15 % headroom

## ============================================================
## Option A: base-R plot (no extra packages)
## ============================================================

pdf("minES_sensitivity_sweep_traces.pdf", width = 5, height = 4, useDingbats = FALSE)

par(mar = c(4.5, 4.5, 2.5, 1.5), mgp = c(2.8, 0.7, 0))

plot(
  sweep_df$minES_rare, sweep_df$n_sig_neg_NCS,
  type = "b",
  pch  = 16, cex = 1.1,
  col  = col_neg,
  lwd  = 2,
  xlim = range(minES_grid),
  ylim = c(0, y_max),
  xlab = expression(italic("minES.rare")~"threshold"),
  ylab = sprintf("Significant ontology x epoch cells\n(p two-tailed \u2264 %.2f)", sig_alpha),
  main = "Sensitivity of connectivity score significance\nto rare-genotype effect-size threshold",
  cex.main = 0.9,
  cex.axis = 0.85,
  cex.lab  = 0.9,
  las = 1
)

lines(
  sweep_df$minES_rare, sweep_df$n_sig_pos_NCS,
  type = "b",
  pch  = 17, cex = 1.1,
  col  = col_pos,
  lwd  = 2
)

## x-axis tick marks at every tested value
axis(1, at = minES_grid, labels = sprintf("%.2f", minES_grid),
     cex.axis = 0.7, las = 2, tcl = -0.3)

## reference gridlines
abline(h = seq(0, y_max, by = 5), col = "grey88", lty = 1)
abline(v = minES_grid,             col = "grey88", lty = 1)

## re-draw traces on top of grid
lines(sweep_df$minES_rare, sweep_df$n_sig_neg_NCS, col = col_neg, lwd = 2)
points(sweep_df$minES_rare, sweep_df$n_sig_neg_NCS, pch = 16, col = col_neg, cex = 1.1)
lines(sweep_df$minES_rare, sweep_df$n_sig_pos_NCS, col = col_pos, lwd = 2)
points(sweep_df$minES_rare, sweep_df$n_sig_pos_NCS, pch = 17, col = col_pos, cex = 1.1)

legend(
  "topright",
  legend = c("NCS < 0  (drug opposes rare trend)",
             "NCS > 0  (drug mirrors rare trend)"),
  col    = c(col_neg, col_pos),
  pch    = c(16, 17),
  lwd    = 2,
  bty    = "n",
  cex    = 0.78
)

dev.off()


## ============================================================
## Option B: ggplot2 version (richer theming, easier to extend)
## ============================================================

if (requireNamespace("ggplot2", quietly = TRUE)) {

  library(ggplot2)

  ## Reshape to long format for ggplot
  sweep_long <- rbind(
    data.frame(
      minES_rare = sweep_df$minES_rare,
      count      = sweep_df$n_sig_neg_NCS,
      direction  = "NCS < 0  (drug opposes rare trend)",
      stringsAsFactors = FALSE
    ),
    data.frame(
      minES_rare = sweep_df$minES_rare,
      count      = sweep_df$n_sig_pos_NCS,
      direction  = "NCS > 0  (drug mirrors rare trend)",
      stringsAsFactors = FALSE
    )
  )
  sweep_long$direction <- factor(
    sweep_long$direction,
    levels = c("NCS < 0  (drug opposes rare trend)",
               "NCS > 0  (drug mirrors rare trend)")
  )

  dir_colours <- c(
    "NCS < 0  (drug opposes rare trend)" = col_neg,
    "NCS > 0  (drug mirrors rare trend)" = col_pos
  )
  dir_shapes <- c(
    "NCS < 0  (drug opposes rare trend)" = 16,
    "NCS > 0  (drug mirrors rare trend)" = 17
  )

  p_gg <- ggplot(sweep_long,
                 aes(x = minES_rare, y = count,
                     colour = direction, shape = direction)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = dir_colours, name = NULL) +
    scale_shape_manual( values = dir_shapes,  name = NULL) +
    scale_x_continuous(
      breaks = minES_grid,
      labels = sprintf("%.2f", minES_grid)
    ) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    labs(
      x     = expression(italic("minES.rare")~"threshold"),
      y     = sprintf("Significant ontology \u00d7 epoch cells\n(p two-tailed \u2264 %.2f)", sig_alpha),
      title = "Sensitivity of connectivity score significance\nto rare-genotype effect-size threshold"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
      legend.position  = "inside",
      legend.position.inside = c(0.97, 0.97),
      legend.justification   = c("right", "top"),
      legend.background = element_rect(fill = "white", colour = "grey70"),
      legend.key.size   = unit(0.45, "cm"),
      panel.grid.minor  = element_blank(),
      plot.title        = element_text(size = 9, face = "bold")
    )

  ggsave("minES_sensitivity_sweep_traces_ggplot.pdf",
         plot = p_gg, width = 5.5, height = 4, useDingbats = FALSE)

  message("ggplot2 version saved.")

} else {
  message("ggplot2 not available; base-R plot only.")
}

## =============================================================================
## END
## =============================================================================


## =============================================================================
## POST-HOC: per-epoch NCS significance counts from sweep_results_list
## =============================================================================
## Takes the already-computed sweep_results_list (named by minES.rare value)
## and extracts, for each epoch separately, the count of significant
## (p_two_tailed <= sig_alpha) ontology cells with NCS < 0 and NCS > 0.
##
## OUTPUT
##   epoch_sweep_df  - long-format data frame with columns:
##                     minES_rare, epoch, direction, count
##   Two PDF plots   - base-R and ggplot2 versions, one curve per epoch
##                     per direction (up to 2 x n_epochs curves total)
##
## PREREQUISITES
##   sweep_results_list  - from the sensitivity sweep loop
##   minES_grid          - the grid of tested values
##   sig_alpha           - significance threshold used previously (default 0.05)
##   add_FDR_to_permutation_stats() defined (from CMAP_permutation_stats.R)
## =============================================================================

sig_alpha <- 0.05    # adjust if you used a different threshold in the sweep

## ---------------------------------------------------------------------------
## 1.  Extract per-epoch counts from sweep_results_list
## ---------------------------------------------------------------------------

## Recover the epoch names from the first non-NULL result
ep_names <- colnames(sweep_results_list[[1]]$p_two_tailed)
n_epochs <- length(ep_names)

## Build a long-format collector
epoch_rows <- vector("list", length(sweep_results_list) * n_epochs)
row_idx <- 1L

for (k in seq_along(sweep_results_list)) {

  this_minES   <- minES_grid[k]
  perm_stats   <- sweep_results_list[[k]]

  ## Add FDR if not already present (idempotent - safe to call twice)
  if (is.null(perm_stats$FDR_p_two_tailed)) {
    perm_stats <- add_FDR_to_permutation_stats(perm_stats)
  }

  p_mat   <- perm_stats$p_two_tailed   # swap for FDR_p_two_tailed if preferred
  ncs_mat <- perm_stats$NCS

  for (ep in ep_names) {

    p_col   <- p_mat[,   ep]
    ncs_col <- ncs_mat[, ep]

    is_sig     <- is.finite(p_col)   & (p_col   <= sig_alpha)
    is_neg_ncs <- is.finite(ncs_col) & (ncs_col <  0)
    is_pos_ncs <- is.finite(ncs_col) & (ncs_col >  0)

    epoch_rows[[row_idx]] <- data.frame(
      minES_rare    = this_minES,
      epoch         = ep,
      n_sig_neg_NCS = sum(is_sig & is_neg_ncs),
      n_sig_pos_NCS = sum(is_sig & is_pos_ncs),
      stringsAsFactors = FALSE
    )
    row_idx <- row_idx + 1L
  }
}

epoch_sweep_df <- do.call(rbind, epoch_rows)
rownames(epoch_sweep_df) <- NULL

## Wide version for quick inspection
print(epoch_sweep_df)
write.csv(epoch_sweep_df,
          "minES_sensitivity_sweep_per_epoch.csv",
          row.names = FALSE)


## ---------------------------------------------------------------------------
## 2.  Reshape to long format (one row per minES x epoch x direction)
## ---------------------------------------------------------------------------

epoch_long <- rbind(
  data.frame(
    minES_rare = epoch_sweep_df$minES_rare,
    epoch      = epoch_sweep_df$epoch,
    direction  = "NCS < 0",
    count      = epoch_sweep_df$n_sig_neg_NCS,
    stringsAsFactors = FALSE
  ),
  data.frame(
    minES_rare = epoch_sweep_df$minES_rare,
    epoch      = epoch_sweep_df$epoch,
    direction  = "NCS > 0",
    count      = epoch_sweep_df$n_sig_pos_NCS,
    stringsAsFactors = FALSE
  )
)

## Factor ordering: keeps epoch order as it appears in the data
epoch_long$epoch     <- factor(epoch_long$epoch,     levels = ep_names)
epoch_long$direction <- factor(epoch_long$direction, levels = c("NCS < 0", "NCS > 0"))


## ---------------------------------------------------------------------------
## 3.  Colour / line-type palette
## ---------------------------------------------------------------------------
## One hue per epoch (colour-blind-friendly Okabe-Ito palette);
## solid line for NCS > 0, dashed for NCS < 0.
## Adjust colours / labels below if your epoch names differ.

okabe_ito <- c("#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

epoch_colours <- setNames(
  okabe_ito[seq_len(n_epochs)],
  ep_names
)

## Point shapes: one per epoch (filled circle, triangle, square, ...)
epoch_shapes <- setNames(
  c(16, 17, 15, 18, 8, 7, 12)[seq_len(n_epochs)],
  ep_names
)

## Line types: NCS > 0 = solid (1), NCS < 0 = dashed (2)
dir_lty <- c("NCS > 0" = 1, "NCS < 0" = 2)


## ---------------------------------------------------------------------------
## 4a.  Base-R plot
## ---------------------------------------------------------------------------

y_max <- ceiling(max(epoch_long$count, na.rm = TRUE) * 1.15)

pdf("minES_sensitivity_sweep_per_epoch_traces.pdf",
    width = 6.5, height = 4.5, useDingbats = FALSE)

par(mar = c(5, 4.5, 3, 1.5), mgp = c(2.8, 0.7, 0))

## Empty frame
plot(
  NA,
  xlim = range(minES_grid),
  ylim = c(0, y_max),
  xlab = expression(italic("minES.rare")~"threshold"),
  ylab = sprintf("Significant ontologies  (p \u2264 %.2f)", sig_alpha),
  main = "Per-epoch sensitivity of connectivity score significance",
  cex.main = 0.9, cex.axis = 0.82, cex.lab = 0.9,
  las = 1, xaxt = "n"
)

axis(1, at = minES_grid, labels = sprintf("%.2f", minES_grid),
     cex.axis = 0.72, las = 2, tcl = -0.3)

## Light gridlines
abline(h = seq(0, y_max, by = 5), col = "grey90", lty = 1)
abline(v = minES_grid,             col = "grey90", lty = 1)

## Draw one trace per epoch x direction
for (ep in ep_names) {
  for (dir in c("NCS > 0", "NCS < 0")) {
    sub <- epoch_long[epoch_long$epoch == ep & epoch_long$direction == dir, ]
    sub <- sub[order(sub$minES_rare), ]

    lines(
      sub$minES_rare, sub$count,
      col = epoch_colours[ep],
      lty = dir_lty[dir],
      lwd = 1.8
    )
    points(
      sub$minES_rare, sub$count,
      col = epoch_colours[ep],
      pch = epoch_shapes[ep],
      cex = 1.0
    )
  }
}

## Legend: epoch colours (left block) + line type for direction (right block)
legend(
  "topright",
  legend = ep_names,
  col    = epoch_colours,
  pch    = epoch_shapes,
  lwd    = 1.8,
  lty    = 1,
  bty    = "n",
  cex    = 0.72,
  title  = "Epoch"
)

## Small inset legend for line type
legend(
  "right",
  legend = c("NCS > 0  (concordant)", "NCS < 0  (opposing)"),
  lty    = c(1, 2),
  col    = "grey30",
  lwd    = 1.8,
  bty    = "n",
  cex    = 0.72,
  title  = "Direction"
)

dev.off()


## ---------------------------------------------------------------------------
## 4b.  ggplot2 version
## ---------------------------------------------------------------------------

if (requireNamespace("ggplot2", quietly = TRUE)) {

  library(ggplot2)

  ## Readable facet / legend labels: strip epoch names from the list names
  ## (edit the label map below if your epoch names differ from these defaults)
  epoch_labels <- setNames(ep_names, ep_names)   # identity; override as needed
  # e.g.: epoch_labels <- c("EYO -46 to -15.5" = "Pre-clinical (EYO < -15)",
  #                          "EYO -15 to -0.5"  = "Peri-onset",
  #                          "EYO 0 to +25"     = "Post-onset")

  p_gg <- ggplot(
    epoch_long,
    aes(
      x        = minES_rare,
      y        = count,
      colour   = epoch,
      linetype = direction,
      shape    = epoch
    )
  ) +
    geom_line(linewidth = 0.85) +
    geom_point(size = 2.2) +
    scale_colour_manual(
      values = epoch_colours,
      labels = epoch_labels,
      name   = "Epoch"
    ) +
    scale_shape_manual(
      values = epoch_shapes,
      labels = epoch_labels,
      name   = "Epoch"
    ) +
    scale_linetype_manual(
      values = c("NCS > 0" = "solid", "NCS < 0" = "dashed"),
      labels = c("NCS > 0" = "Concordant  (NCS > 0)",
                 "NCS < 0" = "Opposing    (NCS < 0)"),
      name   = "Direction"
    ) +
    scale_x_continuous(
      breaks = minES_grid,
      labels = sprintf("%.2f", minES_grid)
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      expand = c(0, 0)
    ) +
    labs(
      x     = expression(italic("minES.rare")~"threshold"),
      y     = sprintf("Significant ontologies  (p \u2264 %.2f)", sig_alpha),
      title = "Per-epoch sensitivity of connectivity score significance"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "right",
      legend.key.size = unit(0.5, "cm"),
      legend.text     = element_text(size = 8),
      legend.title    = element_text(size = 8, face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 9, face = "bold")
    ) +
    ## Merge colour + shape guides into one block; keep linetype separate
    guides(
      colour   = guide_legend(order = 1, override.aes = list(linewidth = 1)),
      shape    = guide_legend(order = 1),
      linetype = guide_legend(order = 2)
    )

  ggsave(
    "minES_sensitivity_sweep_per_epoch_traces_ggplot.pdf",
    plot   = p_gg,
    width  = 7,
    height = 4.5,
    useDingbats = FALSE
  )

  message("ggplot2 per-epoch plot saved.")

} else {
  message("ggplot2 not available; base-R plot only.")
}

## =============================================================================
## END
## =============================================================================
save.image("saved.image-SemaS6_minESsweep.RData")


#load("saved.image-SemaS6_minESsweep.RData")

# S2 unadj statistics sweep:


# ------------------------------------------------------------------------
# ANNOTATION: Save the sensitivity-sweep workspace for later plotting.
# ------------------------------------------------------------------------
## ---- sweep parameters -------------------------------------------------------
minES_grid <- seq(0.01, 0.20, by = 0.02)   # 10 values: 0.01 0.03 0.05 ... 0.19
sig_alpha  <- 0.05                          # p_two_tailed threshold
n_perm     <- 100000                          # permutations per cell
ncores_use <- max(1L, parallel::detectCores() - 1L)
sweep_seed <- 42

## ---- epoch definition (edit to match your analysis) ------------------------
epoch_list <- list(
  "EYO -46 to -15.5" = 1:62,
  "EYO -15 to -0.5"  = 63:92,
  "EYO 0 to +25"     = 93:143
)

## ---- result collector -------------------------------------------------------
sweep_df.S2 <- data.frame(
  minES_rare    = minES_grid,
  n_sig_neg_NCS = NA_integer_,   # significant cells with NCS < 0  (opposing drug)
  n_sig_pos_NCS = NA_integer_,   # significant cells with NCS > 0  (concordant)
  stringsAsFactors = FALSE
)

## Optional: store each full result if RAM allows (set to FALSE to save memory)
store_full_results <- TRUE
sweep_results_list <- if (store_full_results) vector("list", length(minES_grid)) else NULL
names(sweep_results_list) <- as.character(round(minES_grid, 4))

## ---- sweep loop -------------------------------------------------------------
message(sprintf("[%s]  Starting minES.rare S2 (unadj) sensitivity sweep (%d values) ...",
                format(Sys.time(), "%H:%M:%S"), length(minES_grid)))

for (k in seq_along(minES_grid)) {

  this_minES <- minES_grid[k]
  message(sprintf("[%s]  >> minES.rare = %.2f  (%d / %d)",
                  format(Sys.time(), "%H:%M:%S"), this_minES, k, length(minES_grid)))

  perm_stats <- build_Cscore_permutation_stats(
    onts141.symbolList = onts141.symbolList,
    drug_ranked_stats  = drug_ranked_stats.S2,
    effectSizes        = effectSizes,
    epoch_list         = epoch_list,
    minES.rare         = this_minES,
    n_perm             = n_perm,
    ncores             = ncores_use,
    seed               = sweep_seed
  )

  ## Apply BH-FDR across all cells (consistent with single-run usage)
  perm_stats <- add_FDR_to_permutation_stats(perm_stats)

  ## ---- count significant cells by NCS sign --------------------------------
  ## Use raw p_two_tailed here; swap for FDR_p_two_tailed if preferred.
  ## Each "cell" is one (ontology x epoch) combination.
  p_mat   <- perm_stats$p_two_tailed   # 141 x 3
  ncs_mat <- perm_stats$NCS            # 141 x 3

  is_sig     <- is.finite(p_mat)   & (p_mat   <= sig_alpha)
  is_neg_ncs <- is.finite(ncs_mat) & (ncs_mat <  0)
  is_pos_ncs <- is.finite(ncs_mat) & (ncs_mat >  0)

  sweep_df.S2$n_sig_neg_NCS[k] <- sum(is_sig & is_neg_ncs)
  sweep_df.S2$n_sig_pos_NCS[k] <- sum(is_sig & is_pos_ncs)

  if (store_full_results) {
    sweep_results_list[[k]] <- perm_stats
  }

  message(sprintf("         sig & NCS<0 = %d  |  sig & NCS>0 = %d",
                  sweep_df.S2$n_sig_neg_NCS[k],
                  sweep_df.S2$n_sig_pos_NCS[k]))
}

message(sprintf("[%s]  S2 Sweep complete.", format(Sys.time(), "%H:%M:%S")))

## ---- inspect summary table --------------------------------------------------
print(sweep_df.S2)
write.csv(sweep_df.S2, "Unadj(S2)minES_sensitivity_sweep_summary.csv", row.names = FALSE)

if (store_full_results) {
  saveRDS(sweep_results_list, "Unadj(S2)minES_sensitivity_sweep_full_results.rds")
}


#########################
## pretty Figure ready 18 category, 141 ontology ordered heatmap plotter

library(ComplexHeatmap)
library(circlize)
library(grid)
library(data.table)


ontology_order_tbl <- fread("Fig6A_ontology_order_110_18categories.tsv")


## ============================================================
## 1. Define the exact 18 category order from panel A
## ============================================================

category_levels <- c(
  "Metabolism",
  "RNA Processing",
  "GLP-1/IGFBP",

# ------------------------------------------------------------------------
# ANNOTATION: Prepare ordered category heatmap layout for manuscript-ready
# CMAP visualization.
# ------------------------------------------------------------------------
  "Lipoproteins",
  "Proteostasis",
  "Synaptic/Neuronal",
  "Eye and Retinol",
  "Steroid Metab.",
  "ECM",
  "Structural",
  "FGF/Fibronectin",
  "Hemostasis",
  "Adaptive Immune",
  "Innate Immune",
  "Wnt Signaling",
  "Apoptosis",
  "Translation",
  "Angiogenesis"
)

## ============================================================
## 2. Required: exact 110/141-row ontology/category order table
## ============================================================
## ontology_order_tbl must contain:
##   Category: one of the 18 labels above
##   Ontology: rowname exactly matching rownames(perm_stats_S6$NCS)
##
## Example:
## ontology_order_tbl <- fread("Fig6A_ontology_order_141.tsv")
##
## Or construct from a named list:
##
## ontology_categories <- list(
##   "Metabolism" = c("Alkanesulfonate Metabolic Process (GOBP)", ...),
##   "RNA Processing" = c(...),
##   ...
## )
##
## ontology_order_tbl <- data.frame(
##   Category = rep(names(ontology_categories), lengths(ontology_categories)),
##   Ontology = unlist(ontology_categories, use.names = FALSE),
##   stringsAsFactors = FALSE
## )

if (!exists("ontology_order_tbl")) {
  stop(
    "Create ontology_order_tbl first. It must be a 141-row data.frame/data.table ",
    "with columns Category and Ontology in the exact desired row order."
  )
}

ontology_order_tbl <- as.data.frame(ontology_order_tbl, stringsAsFactors = FALSE)

stopifnot(all(c("Category", "Ontology") %in% colnames(ontology_order_tbl)))

if (nrow(ontology_order_tbl) != 110L) {
  stop("ontology_order_tbl must contain exactly 110 rows; found ", nrow(ontology_order_tbl), ".")
}

if (anyDuplicated(ontology_order_tbl$Ontology)) {
  stop("ontology_order_tbl$Ontology contains duplicate ontology names.")
}

ontology_order_tbl$Category <- as.character(ontology_order_tbl$Category)

unknown_categories <- setdiff(unique(ontology_order_tbl$Category), category_levels)
if (length(unknown_categories) > 0) {
  stop(
    "Unknown category labels in ontology_order_tbl: ",
    paste(unknown_categories, collapse = ", ")
  )
}

## Enforce that the categories are contiguous blocks and occur in panel-A order.
category_rle <- rle(ontology_order_tbl$Category)

if (!identical(category_rle$values, category_levels)) {
  stop(
    "Category blocks are not in the expected 18-category panel-A order, ",
    "or at least one category is split into non-contiguous blocks.\n\nObserved order:\n",
    paste(category_rle$values, collapse = " | "), "\n\nExpected order:\n",
    paste(category_levels, collapse = " | ")
  )
}

## ============================================================
## 3. Pull and strictly reorder matrices from sweep_results_list
## ============================================================

sweep_ids <- seq_along(sweep_results_list)

ncs_list <- lapply(sweep_ids, function(k) {
  m <- as.matrix(sweep_results_list[[k]]$NCS)
  storage.mode(m) <- "numeric"
  m[ont_order, , drop = FALSE]
})

fdr_list <- lapply(sweep_ids, function(k) {
  m <- as.matrix(sweep_results_list[[k]]$p_two_tailed)
  storage.mode(m) <- "numeric"
  m[ont_order, colnames(ncs_list[[k]]), drop = FALSE]
})

names(ncs_list) <- names(sweep_results_list) #sprintf("Sweep_%02d", sweep_ids)
names(fdr_list) <- names(ncs_list)

## Optional sanity check: each sweep contributes 3 heatmap columns
stopifnot(all(vapply(ncs_list, ncol, integer(1)) == 3L))
stopifnot(all(vapply(ncs_list, nrow, integer(1)) == length(ont_order)))

## ============================================================
## 4. Significance stars for each sweep
## ============================================================

star_for_p <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(
      p <= 0.001, "***",
      ifelse(
        p <= 0.01, "**",
        ifelse(p <= 0.05, "*", "")
      )
    )
  )
}

star_list <- lapply(fdr_list, function(fdr_mat) {
  matrix(
    star_for_p(fdr_mat),
    nrow = nrow(fdr_mat),
    ncol = ncol(fdr_mat),
    dimnames = dimnames(fdr_mat)
  )
})

## ============================================================
## 5a. Figure-A style labels / layout
## ============================================================

col_ncs <- colorRamp2(
  c(-3, 0, 3),
  c("#0099FF", "white", "gold")
)

ontology_fontsize <- 4.7
category_fontsize <- 9

ha_category <- rowAnnotation(
  Category = anno_block(
    labels = category_levels,
    labels_gp = gpar(fontsize = category_fontsize, fontface = "bold"),
    labels_rot = 0,
    gp = gpar(fill = NA, col = NA),
    width = unit(1.25, "in")
  ),
  show_annotation_name = FALSE
)

## PATCH: thick black vertical/block separator to right of category names
ha_category_rule <- rowAnnotation(
  CategoryRule = anno_block(
    labels = rep("", length(category_levels)),
    gp = gpar(fill = NA, col = "black", lwd = 1.25),
    width = unit(.01, "mm")
  ),
  show_annotation_name = FALSE
)

ha_ontology <- rowAnnotation(
  Ontology = anno_text(
    ont_order,
    gp = gpar(fontsize = ontology_fontsize),
    just = "right",
    location = unit(1, "npc"),
    width = max_text_width(
      ont_order,
      gp = gpar(fontsize = ontology_fontsize)
    ) + unit(1.5, "mm")
  ),
  show_annotation_name = FALSE
)

## ============================================================
## 5b. Build one 3-column heatmap per sweep_results_list element
## ============================================================

ht_ncs_list <- lapply(names(ncs_list), function(nm) {
  this_mat  <- ncs_list[[nm]]
  this_star <- star_list[[nm]]

  Heatmap(
    this_mat,
    name = paste0("NCS_", nm),
    col = col_ncs,

    cluster_rows = FALSE,
    cluster_columns = FALSE,

    row_split = row_split,
    row_gap = unit(1.6, "mm"),

    ## remove duplicated rotated row-split labels
    row_title = NULL,
    row_title_gp = gpar(fontsize = 0),

    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 7),
    column_names_rot = 45,

    column_title = nm,
    column_title_gp = gpar(fontsize = 8, fontface = "bold"),

    width = unit(1.5, "cm"),

    border = TRUE,
    rect_gp = gpar(col = "#CFCFCF", lwd = 0.35),

    heatmap_legend_param = list(
      title = "NCS",
      at = c(-3, -1, 0, 1, 3),
      labels = c("-3", "-1", "0", "1", "3")
    ),

    cell_fun = function(j, i, x, y, width, height, fill) {
      s <- this_star[i, j]
      if (!is.na(s) && nzchar(s)) {
        grid.text(
          s,
          x = x,
          y = y,
          gp = gpar(fontsize = 7, fontface = "bold")
        )
      }
    }
  )
})

names(ht_ncs_list) <- names(ncs_list)

ht_list <- ha_category + ha_category_rule + ha_ontology

for (nm in names(ht_ncs_list)) {
  ht_list <- ht_list + ht_ncs_list[[nm]]
}

## ============================================================
## 6. Draw PDF with sweep groups left-to-right
## ============================================================

pdf(
  "Ontology_NCS_Fig6A_style_sweep_results_list.pdf",
  width = 15,
  height = 17,
  useDingbats = FALSE
)

draw(
  ht_list,
  heatmap_legend_side = "right",
  merge_legend = TRUE,
  ht_gap = unit(3, "mm"),
  padding = unit(c(3, 3, 3, 3), "mm")
)

## Draw black rectangle around each category slice for every 3-column NCS block
for (nm in names(ht_ncs_list)) {
  ht_name <- paste0("NCS_", nm)

  for (si in seq_along(category_levels)) {
    decorate_heatmap_body(ht_name, slice = si, {
      grid.rect(gp = gpar(fill = NA, col = "black", lwd = 0.8))
    })
  }
}

dev.off()


## ---- sweep parameters -------------------------------------------------------
minES_grid2 <- seq(0, 0.20, by = 0.10)   # 10 values: 0.01 0.03 0.05 ... 0.19
sig_alpha  <- 0.05                          # p_two_tailed threshold
n_perm     <- 100000                          # permutations per cell
ncores_use <- max(1L, parallel::detectCores() - 1L)
sweep_seed <- 42

## ---- epoch definition (edit to match your analysis) ------------------------
epoch_list <- list(
  "EYO -46 to -15.5" = 1:62,
  "EYO -15 to -0.5"  = 63:92,
  "EYO 0 to +25"     = 93:143
)

## ---- result collector -------------------------------------------------------
sweep_df2 <- data.frame(
  minES_rare    = minES_grid2,
  n_sig_neg_NCS = NA_integer_,   # significant cells with NCS < 0  (opposing drug)
  n_sig_pos_NCS = NA_integer_,   # significant cells with NCS > 0  (concordant)
  stringsAsFactors = FALSE
)

## Optional: store each full result if RAM allows (set to FALSE to save memory)
store_full_results <- TRUE
sweep_results_list2 <- if (store_full_results) vector("list", length(minES_grid2)) else NULL
names(sweep_results_list2) <- as.character(round(minES_grid2, 4))

## ---- sweep loop -------------------------------------------------------------
message(sprintf("[%s]  Starting minES.rare sensitivity sweep (%d values) ...",
                format(Sys.time(), "%H:%M:%S"), length(minES_grid2)))

for (k in seq_along(minES_grid2)) {

  this_minES <- minES_grid2[k]
  message(sprintf("[%s]  >> minES.rare = %.2f  (%d / %d)",
                  format(Sys.time(), "%H:%M:%S"), this_minES, k, length(minES_grid2)))

  perm_stats2 <- build_Cscore_permutation_stats(
    onts141.symbolList = onts141.symbolList,
    drug_ranked_stats  = drug_ranked_stats,
    effectSizes        = effectSizes,
    epoch_list         = epoch_list,
    minES.rare         = this_minES,
    n_perm             = n_perm,
    ncores             = ncores_use,
    seed               = sweep_seed
  )

  ## Apply BH-FDR across all cells (consistent with single-run usage)
  perm_stats2 <- add_FDR_to_permutation_stats(perm_stats2)

  ## ---- count significant cells by NCS sign --------------------------------
  ## Use raw p_two_tailed here; swap for FDR_p_two_tailed if preferred.
  ## Each "cell" is one (ontology x epoch) combination.
  p_mat   <- perm_stats2$p_two_tailed   # 141 x 3
  ncs_mat <- perm_stats2$NCS            # 141 x 3

  is_sig     <- is.finite(p_mat)   & (p_mat   <= sig_alpha)
  is_neg_ncs <- is.finite(ncs_mat) & (ncs_mat <  0)
  is_pos_ncs <- is.finite(ncs_mat) & (ncs_mat >  0)

  sweep_df2$n_sig_neg_NCS[k] <- sum(is_sig & is_neg_ncs)
  sweep_df2$n_sig_pos_NCS[k] <- sum(is_sig & is_pos_ncs)

  if (store_full_results) {
    sweep_results_list2[[k]] <- perm_stats2
  }

  message(sprintf("         sig & NCS<0 = %d  |  sig & NCS>0 = %d",
                  sweep_df2$n_sig_neg_NCS[k],
                  sweep_df2$n_sig_pos_NCS[k]))
}

message(sprintf("[%s]  Sweep complete.", format(Sys.time(), "%H:%M:%S")))


save.image("saved.image-CMAPperm10+3_minES_in2lists(sweep_results_list).RData")


names(sweep_results_list)
names(sweep_results_list2)

sweep_results_list3<-list('0'=sweep_results_list2[["0"]], '0.05'=sweep_results_list[["0.05"]], '0.1'=sweep_results_list2[["0.1"]], '0.15'=sweep_results_list[["0.15"]], '0.2'=sweep_results_list2[["0.2"]])

save.image("saved.image-CMAP.perm.RData")
