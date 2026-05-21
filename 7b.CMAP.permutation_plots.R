##############################################################################
# Pipeline annotation header: 7b.CMAP.permutation_plots.R
# Manuscript code section(s): 7
#
# Purpose:
# Replot CMAP permutation results as manuscript-ready ordered heatmaps and
# extended sensitivity sweep visualizations.
#
# Principal inputs:
#   - saved.image-CMAP.perm.RData
#   - Fig6A_ontology_order_110_18categories_final.tsv
#
# Principal outputs:
#   - Manuscript-ready CMAP NCS heatmaps
#   - CMAP_permutation_stats_S6_100000perm_minES.e4_0.10.xlsx
#   - CMAPperm21_minES_in1list*.RDS
#   - minES_sensitivity_sweep_per_epoch_traces(21points).pdf
#
# Step overview:
#   1. Validate the 18-category ontology ordering table and enforce contiguous
#      category blocks.
#   2. Extract NCS and p-value matrices from each sensitivity-sweep result
#      object and arrange each sweep as three heatmap columns.
#   3. Add significance glyphs, category label blocks, and black slice borders
#      with ComplexHeatmap.
#   4. Export all statistic matrices to Excel and generate additional per-
#      epoch sensitivity traces.
#
# Notes added during manuscript-code alignment:
#   - Original executable code was preserved except for whitespace/encoding cleanup
#     and the explicitly marked non-run archive described below when applicable.
#   - Comments inserted by this pass are prefixed with "ANNOTATION:".
##############################################################################

# ------------------------------------------------------------------------
# ANNOTATION: Load saved CMAP permutation results and plotting packages for
# manuscript re-rendering.
# ------------------------------------------------------------------------
## =============================================================================
## CMAP Connectivity Score - Permutation-Based Statistical Testing  - replot
## =============================================================================

#########################
## pretty Figure ready 18 category, 141 ontology ordered heatmap plotter


setwd("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/#manuscript/SciTranslMed_plan/CMAP.perm/")


load("saved.image-CMAP.perm.RData")


library(ComplexHeatmap)
library(circlize)
library(grid)
library(data.table)


ontology_order_tbl <- fread("Fig6A_ontology_order_110_18categories_final.tsv")


## ============================================================

# ------------------------------------------------------------------------
# ANNOTATION: Define the 18 ontology category order used in the final
# figure.
# ------------------------------------------------------------------------
## 1. Define the exact 18 category order from panel A
## ============================================================

category_levels <- c(
  "Metabolism",
  "RNA Processing",
  "GLP-1/IGFBP",
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

# ------------------------------------------------------------------------
# ANNOTATION: Validate the ordered ontology/category table used for row
# layout.
# ------------------------------------------------------------------------
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


ont_order<-ontology_order_tbl$Ontology  #overwrites, correcting a couple of disordered terms


## ============================================================

# ------------------------------------------------------------------------
# ANNOTATION: Extract and row-order NCS matrices from each sensitivity-sweep
# result.
# ------------------------------------------------------------------------
## 3. Pull and strictly reorder matrices from sweep_results_list
## ============================================================

sweep_ids <- seq_along(sweep_results_list3)

ncs_list <- lapply(sweep_ids, function(k) {
  m <- as.matrix(sweep_results_list3[[k]]$NCS)
  storage.mode(m) <- "numeric"
  m[ont_order, , drop = FALSE]
})

fdr_list <- lapply(sweep_ids, function(k) {
  m <- as.matrix(sweep_results_list3[[k]]$p_two_tailed)
  storage.mode(m) <- "numeric"
  m[ont_order, colnames(ncs_list[[k]]), drop = FALSE]
})

names(ncs_list) <- names(sweep_results_list3) #sprintf("Sweep_%02d", sweep_ids)
names(fdr_list) <- names(ncs_list)

## Optional sanity check: each sweep contributes 3 heatmap columns
stopifnot(all(vapply(ncs_list, ncol, integer(1)) == 3L))
stopifnot(all(vapply(ncs_list, nrow, integer(1)) == length(ont_order)))

## ============================================================

# ------------------------------------------------------------------------
# ANNOTATION: Create significance glyph matrices from two-tailed permutation
# p values.
# ------------------------------------------------------------------------
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

# ------------------------------------------------------------------------
# ANNOTATION: Build figure labels, color scales, and category side
# annotations.
# ------------------------------------------------------------------------
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
## 5b. Build one 3-column heatmap per sweep_results_list3 element

# ------------------------------------------------------------------------
# ANNOTATION: Build one three-column ComplexHeatmap block per sensitivity
# sweep.
# ------------------------------------------------------------------------
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
    row_gap = unit(3, "mm"),  #1.6 prev

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

# ------------------------------------------------------------------------
# ANNOTATION: Draw the final PDF and decorate each category block with black
# borders.
# ------------------------------------------------------------------------
## 6. Draw PDF with sweep groups left-to-right
## ============================================================

pdf(
  "Ontology_NCS_Fig6A_style_sweep_results_list3_finalOrder.pdf",
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


## ---------------------------------------------------------------------------
## 4.  Write all stat matrices to an Excel workbook (one sheet per stat)
## ---------------------------------------------------------------------------

perm_stats<-sweep_results_list3$'0.1'


# ------------------------------------------------------------------------
# ANNOTATION: Export all permutation statistic matrices to an Excel
# workbook.
# ------------------------------------------------------------------------
write_permutation_stats_excel <- function(perm_stats,
                                          file = "CMAP_permutation_stats(0.10_ESthresh).xlsx",
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


  ## ---- Save all stat matrices to Excel ------------------------------------
  write_permutation_stats_excel(
    perm_stats_S6,
    file = "CMAP_permutation_stats_S6_100000perm_minES.e4_0.10.xlsx"
  )
  ## Sheets: Cscore, NCS, p_two_tailed, p_directional,
  ##         ci_lower, ci_upper, null_mean, null_sd,
  ##         n_up, n_down, n_drug_assays,
  ##         FDR_p_two_tailed, FDR_p_directional


## RUN BELOW HERE 05/19/2026 (on telomere /R/)


## ---- sweep parameters -------------------------------------------------------
minES_grid2 <- c(0.02, 0.04, 0.06, 0.08, 0.12, 0.14, 0.16, 0.18)  #seq(0, 0.20, by = 0.10)   # 10 values: 0.01 0.03 0.05 ... 0.19
sig_alpha  <- 0.05                          # p_two_tailed threshold
n_perm     <- 100000                          # permutations per cell

# ------------------------------------------------------------------------
# ANNOTATION: Optional/run-on-server block for computing a denser 21-point
# sensitivity sweep.
# ------------------------------------------------------------------------
ncores_use <- max(1L, parallel::detectCores() - 1L)
sweep_seed <- 42

## ---- epoch definition (edit to match your analysis) ------------------------
epoch_list <- list(
  "EYO -46 to -15.5" = 1:62,
  "EYO -15 to -0.5"  = 63:92,
  "EYO 0 to +25"     = 93:143
)

## ---- result collector -------------------------------------------------------
sweep_df4 <- data.frame(
  minES_rare    = minES_grid2,
  n_sig_neg_NCS = NA_integer_,   # significant cells with NCS < 0  (opposing drug)
  n_sig_pos_NCS = NA_integer_,   # significant cells with NCS > 0  (concordant)
  stringsAsFactors = FALSE
)

## Optional: store each full result if RAM allows (set to FALSE to save memory)
store_full_results <- TRUE
sweep_results_list4 <- if (store_full_results) vector("list", length(minES_grid2)) else NULL
names(sweep_results_list4) <- as.character(round(minES_grid2, 4))

## ---- sweep loop -------------------------------------------------------------
message(sprintf("[%s]  Starting minES.rare sensitivity sweep (%d values) ...",
                format(Sys.time(), "%H:%M:%S"), length(minES_grid2)))

for (k in seq_along(minES_grid2)) {

  this_minES <- minES_grid2[k]
  message(sprintf("[%s]  >> minES.rare = %.2f  (%d / %d)",
                  format(Sys.time(), "%H:%M:%S"), this_minES, k, length(minES_grid2)))

  perm_stats4 <- build_Cscore_permutation_stats(
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
  perm_stats4 <- add_FDR_to_permutation_stats(perm_stats4)

  ## ---- count significant cells by NCS sign --------------------------------
  ## Use raw p_two_tailed here; swap for FDR_p_two_tailed if preferred.
  ## Each "cell" is one (ontology x epoch) combination.
  p_mat   <- perm_stats4$p_two_tailed   # 141 x 3
  ncs_mat <- perm_stats4$NCS            # 141 x 3

  is_sig     <- is.finite(p_mat)   & (p_mat   <= sig_alpha)
  is_neg_ncs <- is.finite(ncs_mat) & (ncs_mat <  0)
  is_pos_ncs <- is.finite(ncs_mat) & (ncs_mat >  0)

  sweep_df4$n_sig_neg_NCS[k] <- sum(is_sig & is_neg_ncs)
  sweep_df4$n_sig_pos_NCS[k] <- sum(is_sig & is_pos_ncs)

  if (store_full_results) {
    sweep_results_list4[[k]] <- perm_stats4
  }

  message(sprintf("         sig & NCS<0 = %d  |  sig & NCS>0 = %d",
                  sweep_df4$n_sig_neg_NCS[k],
                  sweep_df4$n_sig_pos_NCS[k]))
}

message(sprintf("[%s]  Sweep complete.", format(Sys.time(), "%H:%M:%S")))


sweep_results_list5<-list('0'=sweep_results_list3[["0"]], '0.01'=sweep_results_list[["0.01"]], '0.02'=sweep_results_list4[["0.02"]], '0.03'=sweep_results_list[["0.03"]], '0.04'=sweep_results_list4[["0.04"]],
                       '0.05'=sweep_results_list3[["0.05"]],'0.06'=sweep_results_list4[["0.04"]],'0.07'=sweep_results_list[["0.07"]], '0.08'=sweep_results_list4[["0.08"]],
                       '0.09'=sweep_results_list[["0.09"]], '0.10'=sweep_results_list3[["0.1"]],'0.11'=sweep_results_list[["0.11"]], '0.12'=sweep_results_list4[["0.12"]],
                       '0.13'=sweep_results_list[["0.13"]],'0.14'=sweep_results_list4[["0.14"]],'0.15'=sweep_results_list3[["0.15"]], '0.16'=sweep_results_list4[["0.16"]],
                       '0.17'=sweep_results_list[["0.17"]],'0.18'=sweep_results_list4[["0.18"]],'0.19'=sweep_results_list[["0.19"]], '0.20'=sweep_results_list3[["0.2"]])
#saveRDS(sweep_results_list5,"CMAPperm21_minES_in1list(sweep_results_list5).RDS")


names(sweep_results_list)
names(sweep_results_list2)
names(sweep_results_list5)


#sweep_results_list3<-list('0'=sweep_results_list2[["0"]], '0.05'=sweep_results_list[["0.05"]], '0.1'=sweep_results_list2[["0.1"]], '0.15'=sweep_results_list[["0.15"]], '0.2'=sweep_results_list2[["0.2"]])

#save.image("saved.image-CMAP.perm.RData")


## =============================================================================
## POST-HOC: per-epoch NCS significance counts from sweep_results_list
## =============================================================================
## Takes the already-computed sweep_results_list (named by minES.rare value)
## and extracts, for each epoch separately, the count of significant
## (p_two_tailed <= sig_alpha) ontology cells with NCS < 0 and NCS > 0.
##
## OUTPUT

# ------------------------------------------------------------------------
# ANNOTATION: Post-hoc plotting of per-epoch NCS significance counts across
# sensitivity thresholds.
# ------------------------------------------------------------------------
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
minES_grid5<-seq(0, 0.20, by = 0.01)

sig_alpha <- 0.05    # adjust if you used a different threshold in the sweep

## ---------------------------------------------------------------------------
## 1.  Extract per-epoch counts from sweep_results_list
## ---------------------------------------------------------------------------

## Recover the epoch names from the first non-NULL result
ep_names <- colnames(sweep_results_list5[[1]]$p_two_tailed)
n_epochs <- length(ep_names)

## Build a long-format collector
epoch_rows <- vector("list", length(sweep_results_list5) * n_epochs)
row_idx <- 1L

for (k in seq_along(sweep_results_list5)) {

  this_minES   <- minES_grid5[k]
  perm_stats4   <- sweep_results_list5[[k]]

  ## Add FDR if not already present (idempotent - safe to call twice)
  if (is.null(perm_stats4$FDR_p_two_tailed)) {
    perm_stats4 <- add_FDR_to_permutation_stats(perm_stats4)
  }

  p_mat   <- perm_stats4$p_directional #two_tailed   # swap for FDR_p_two_tailed if preferred
  ncs_mat <- perm_stats4$NCS

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
          "minES_sensitivity_sweep_per_epoch(21_ES_thresh-oneTailedP).csv",
          row.names = FALSE)

saveRDS(sweep_results_list5,"CMAPperm21_minES_in1list(sweep_results_list5).RDS")


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

y_max <- ceiling(max(epoch_long$count, na.rm = TRUE) * 1.15)


## ---------------------------------------------------------------------------
## 4a.  Base-R plot
## ---------------------------------------------------------------------------

pdf("minES_sensitivity_sweep_per_epoch_traces(21points).pdf",
    width = 6.5, height = 4.5, useDingbats = FALSE)

par(mar = c(5, 4.5, 3, 1.5), mgp = c(2.8, 0.7, 0))

## Empty frame
plot(
  NA,
  xlim = range(minES_grid5),
  ylim = c(0, y_max),
  xlab = expression(italic("minES.rare")~"threshold"),
  ylab = sprintf("Significant ontologies  (p \u2264 %.2f)", sig_alpha),
  main = "Per-epoch sensitivity of connectivity score significance",
  cex.main = 0.9, cex.axis = 0.82, cex.lab = 0.9,
  las = 1, xaxt = "n"
)

axis(1, at = minES_grid5, labels = sprintf("%.2f", minES_grid5),
     cex.axis = 0.72, las = 2, tcl = -0.3)

## Light gridlines
abline(h = seq(0, y_max, by = 5), col = "grey90", lty = 1)
abline(v = minES_grid5,             col = "grey90", lty = 1)

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
      breaks = minES_grid5,
      labels = sprintf("%.2f", minES_grid5)
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
    "minES_sensitivity_sweep_per_epoch_traces(21points-1tailedP)_ggplot.pdf",
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
