##############################################################################
# Pipeline annotation header: 2b.APOE_prediction-stage_2.R
# Manuscript code section(s): 2
#
# Purpose:
# Stage 2 APOE genotype imputation: train the final six-genotype stacked
# ensemble using the stage 1 ranked feature sets and evaluate with nested
# hold-out folds.
#
# Principal inputs:
#   - rankedProteins.prior
#   - Known APOE genotype samples and first-pass adjusted protein matrix
#
# Principal outputs:
#   - saved.image-genotype_prediction_finalized_2026rerun.RData
#   - genotype_prediction_100outerCVfolds-stats_2026rerun.RDS
#   - confusion matrix summaries
#
# Step overview:
#   1. Force genotype call order to e2/e4, e2/e2, e4/e4, e2/e3, e3/e4, and
#      e3/e3.
#   2. Train one binary ensemble per genotype using the top 50 features
#      selected in stage 1.
#   3. Assemble six binary learners into a multiclass caller, using
#      thresholded calls where available and maximum probability otherwise.
#   4. Run nested outer cross-validation with frequency-aware hold-out sets,
#      then summarize fold-level confusion matrices, accuracy, and macro-F1.
#   5. Train the final all-data model and apply it to samples with unknown
#      APOE genotype.
#
# Notes added during manuscript-code alignment:
#   - Original executable code was preserved except for whitespace/encoding cleanup
#     and the explicitly marked non-run archive described below when applicable.
#   - Comments inserted by this pass are prefixed with "ANNOTATION:".
##############################################################################

# ------------------------------------------------------------------------
# ANNOTATION: Load ML packages and configure the final stage 2 APOE genotype
# imputation workflow.
# ------------------------------------------------------------------------
# Stage 2; continuing in loaded session on VM Windows03
setwd("Z:/EBD/grid/4p13b3forAPOEpredict+2ndRegrAgain/26predict/")
#load("z:/EBD/4p13b4b.SecondRegressionsComplete.19sites_Fsplit+APOEgenoIMPUTATION.RData") #13.7 GB

# Prompt: Given 6 lists of features (log2(protein relative abundances) which predict each of 6 APOE genotypes (e2/e2, e2/e3, e2/e4, e3/e3, e3/e4, and e4/e4) in a binary fashion with 97 percent accuracy in a large (15000 samples) test data set subjected to the attached ML function, I want to devise a strategy (adapted algorithm) and implement it it R code to assemble the most accurate 6-genotype predictor. Because the binary predictors were most accurate for e2/e4 and e2/e2, I suggest predicting e2/e4 first, then e2/e2, then e4/e4, followed by e2/e3, e3/4, and e3/e3. Given the known accuracy of individual binary predictors all 97 percent or better, suggest what the accuracy of the 6-genotype prediction algorithm could be, at best.
library(caret)
library(glmnet)
library(xgboost)
library(ranger)
library(progressr)
library(dplyr)
library(future)
library(future.apply)
library(doFuture)
library(doRNG)

## Prevent RStudio-specific error at very end of function ("error in evaluating the argument 'x' in selecting a method for function 'print': object 'spaces' not found")
#if (!exists("spaces", envir = asNamespace("cli")))
#    assign("spaces", c("", vapply(1:20, strrep, "", x = " ")),
#           envir = asNamespace("cli"))
## Just DO NOT RUN in RStudio.


############################################
## STAGE 2 - Final Ensemble Learner Assembly with 0.80 (helper function's internal floor above 0.50 setting) using top 50 features for each genotype

genotypes <- c("e2/e4","e2/e2","e4/e4","e2/e3","e3/e4","e3/e3")

# ------------------------------------------------------------------------
# ANNOTATION: Define the final binary learner used inside the six-genotype
# stage 2 ensemble.
# ------------------------------------------------------------------------
## force the same ordering as the `genotypes` vector
rankedProteins.prior <- rankedProteins.prior[genotypes]


################################################################################
## 1.  Binary learner ----------------------------------------------------------
################################################################################
fit_APOE_binary <- function(expr, APOE_gt,
                            target,                    # e.g. "e4/e4"
                            target_ppv = 0.95,         # PPV for hard-threshold
                            nfold = 5, nrep = 5,
                            ncores = parallel::detectCores() - 1,
                            seed   = 42)
{
  set.seed(seed)

  ## ---------- helpers --------------------------------------------------------
  prep_expr <- function(mat, ref = NULL) {   #Z scaling (run intra-fold) and handling of missing data (there is none in our input)
    keep <- which(colMeans(is.na(mat)) <= 0.20)
    if (is.null(ref)) {
      X <- scale(mat[, keep, drop = FALSE])

# ------------------------------------------------------------------------
# ANNOTATION: Prepare fold-specific z-scaled expression matrices and reuse
# the training-set scaling on test data.
# ------------------------------------------------------------------------
      list(x = X,
           center = attr(X, "scaled:center"),
           scale  = attr(X, "scaled:scale"),
           vars   = colnames(X))
    } else {
      X <- scale(mat[, ref$vars, drop = FALSE],
                 center = ref$center,
                 scale  = ref$scale)
      X[is.na(X)] <- 0
      list(x = X)
    }
  }

  pick_thr <- function(prob, truth, target_ppv = 0.95,
                       min_tp = 30, floor = 0.80) {
    ok <- !is.na(prob)
    prob <- prob[ok]; truth <- truth[ok]
    if (sum(truth) == 0 || all(prob == 0)) return(1)
    ord <- order(prob, decreasing = TRUE)
    tp <- cumsum(truth[ord] == 1); fp <- cumsum(truth[ord] == 0)
    ppv <- tp / (tp + fp)
    idx <- which(ppv >= target_ppv & tp >= min_tp)
    thr <- if (length(idx)) prob[ord[max(idx)]] else
                         quantile(prob[truth == 1], .90, na.rm = TRUE)
    max(thr, floor)
  }

  ## ---------- inner CV  ------------------------------------------------------
  cvIndex <- caret::createMultiFolds(APOE_gt, k = nfold, times = nrep)

  ## -- future backend (no nested progress objects) ------------

# ------------------------------------------------------------------------
# ANNOTATION: Fit genotype-specific glmnet/xgboost/ranger learners and
# probability thresholds.
# ------------------------------------------------------------------------
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(multisession, workers = ncores)
  doFuture::registerDoFuture()
  doRNG::registerDoRNG(seed)                # reproducible but no %dorng

    cv_stats <- foreach(fold = seq_along(cvIndex),
            .packages = c("glmnet", "xgboost", "ranger", "caret", "progressr", "dplyr", "cli"),
            .combine  = dplyr::bind_rows,
            .export   = c("pick_thr", "APOE_gt", "target", "seed"),
            .options.RNG = seed) %dopar% {

      expr_loc <- expr
      tr <- cvIndex[[fold]]
      te <- setdiff(seq_len(nrow(expr_loc)), tr)

      ## ---------- scaling ----------------------------------------------------
      X_tr <- prep_expr(expr[tr, ])$x
      ref  <- list(center = attr(X_tr, "scaled:center"),
                   scale  = attr(X_tr, "scaled:scale"),
                   vars   = colnames(X_tr))
      X_te <- prep_expr(expr[te, ], ref)$x

      ## ---------- binarise labels -------------------------------------------
      y_tr <- factor(ifelse(APOE_gt[tr] == target, "pos", "neg"))
      y_te <- factor(ifelse(APOE_gt[te] == target, "pos", "neg"))
      w_tr <- ifelse(y_tr == "pos", 8, 1)

      if (length(unique(y_tr)) < 2) {
#        p(message = sprintf("fold %d - skipped (single class)", fold), amount = 3)
        return(data.frame(Precision = NA, Recall = NA))
      }

      ## ---------- glmnet -----------------------------------------------------
      glm_cv <- glmnet::cv.glmnet(X_tr, y_tr, family = "binomial",
                                  weights = w_tr, type.measure = "class")
      p_glm <- drop(predict(glm_cv, X_te, s = "lambda.min",
                            type = "response"))
#      p(amount = 1)

      ## ---------- xgboost ----------------------------------------------------
      dtr <- xgboost::xgb.DMatrix(X_tr, label = as.numeric(y_tr) - 1,
                                  weight = w_tr)
      dte <- xgboost::xgb.DMatrix(X_te, label = as.numeric(y_te) - 1)
      xpar <- list(objective = "binary:logistic",
                   eta = 0.1, max_depth = 6,
                   subsample = 0.8, colsample_bytree = 0.8,
                   nthread = 1, eval_metric = "logloss")
      bst <- xgboost::xgb.train(xpar, dtr,
                                watchlist = list(train = dtr, eval = dte),
                                nrounds = 200, verbose = 0,
                                early_stopping_rounds = 20)
      p_xgb <- drop(predict(bst, dte))
#      p(amount = 1)

      ## ---------- ranger -----------------------------------------------------
      rf <- ranger::ranger(y_tr ~ ., data = data.frame(y_tr, X_tr),
                           probability   = TRUE,
                           num.trees     = 500,
                           num.threads   = 1,
                           class.weights = c(neg = 1, pos = 8))
      p_rf <- predict(rf, data.frame(X_te))$predictions[, "pos"]
#      p(amount = 1)

      ## ---------- evaluate ---------------------------------------------------
      p_avg_tr <- rowMeans(cbind(
        drop(predict(glm_cv, X_tr, s = "lambda.min", type = "response")),
        drop(predict(bst, dtr)),
        predict(rf, data.frame(X_tr))$predictions[, "pos"]))

      thr <- pick_thr(p_avg_tr, APOE_gt[tr] == target, target_ppv)

      p_avg_te <- rowMeans(cbind(p_glm, p_xgb, p_rf))
      pred <- ifelse(p_avg_te >= thr, "pos", "neg")

      tp <- sum(pred == "pos" & y_te == "pos")
      fp <- sum(pred == "pos" & y_te == "neg")
      fn <- sum(pred == "neg" & y_te == "pos")

      data.frame(
        Precision = if ((tp + fp) > 0) tp / (tp + fp) else NA_real_,
        Recall    = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_)
    } #}) no more with_progress

  message(sprintf("[%s]  CV Precision %.3f ± %.3f | Recall %.3f ± %.3f",
                  target,
                  mean(cv_stats$Precision, na.rm = TRUE),
                  sd  (cv_stats$Precision, na.rm = TRUE),
                  mean(cv_stats$Recall,    na.rm = TRUE),
                  sd  (cv_stats$Recall,    na.rm = TRUE)))

  ## ---------- fit once on all data ------------------------------------------
  prep_all <- prep_expr(expr)
  X_all <- prep_all$x
  y_all <- factor(ifelse(APOE_gt == target, "pos", "neg"))
  w_all <- ifelse(y_all == "pos",
                  if (target %in% c("e2/e2", "e2/e4", "e4/e4")) 12 else 8,
                  1)

  glm_fit <- glmnet::cv.glmnet(X_all, y_all, family = "binomial",
                               weights = w_all, type.measure = "class")

  d_all  <- xgboost::xgb.DMatrix(X_all, label = as.numeric(y_all) - 1,
                                 weight = w_all)
  xpar_b <- list(objective = "binary:logistic",
                 eta = 0.1, max_depth = 6,
                 subsample = 0.8, colsample_bytree = 0.8,
                 nthread = ncores, eval_metric = "logloss")
  xgb_fit <- xgboost::xgb.train(xpar_b, d_all,
                                watchlist = list(train = d_all),
                                nrounds   = 200,
                                verbose   = 0,
                                early_stopping_rounds = 20)

  rf_fit <- ranger::ranger(y_all ~ ., data = data.frame(y_all, X_all),
                           probability   = TRUE,
                           num.trees     = 500,
                           class.weights = c(neg = 1,
                                             pos = if (target %in%
                                                       c("e2/e2","e2/e4","e4/e4"))
                                                        12 else 8),
                           num.threads   = ncores)

  p_all <- rowMeans(cbind(
    drop(predict(glm_fit, X_all, s = "lambda.min", type = "response")),
    drop(predict(xgb_fit, d_all)),
    rf_fit$predictions[, "pos"]))

  thr <- pick_thr(p_all, APOE_gt == target, target_ppv)

  ## ---------- wrapper with probability accessor -----------------------------
  prob_fun <- function(new_expr) {
    new_expr <- prep_expr(new_expr, prep_all)$x
    rowMeans(cbind(
      drop(predict(glm_fit, new_expr, s = "lambda.min", type = "response")),
      drop(predict(xgb_fit, new_expr)),
      predict(rf_fit, data.frame(new_expr))$predictions[, "pos"]))
  }

  predict_fun <- function(new_expr) {
    p_bin <- prob_fun(new_expr)
    res   <- factor(ifelse(p_bin >= thr, target, NA_character_),
                    levels = target)
    attr(res, "prob") <- p_bin
    res
  }

  attr(predict_fun, "prob")     <- prob_fun
  attr(predict_fun, "thr")      <- thr
  attr(predict_fun, "features") <- colnames(expr)

  predict_fun
}

################################################################################
## 2.  Helper - train the 6 binary models & build a 6-class predictor ---------
################################################################################
train_APOE6_full <- function(expr, APOE_gt,
                             rankedProteins.prior,        # list of top-50 feats
                             target_ppv = 0.95,
                             ncores = parallel::detectCores() - 1,
                             seed = 42)

# ------------------------------------------------------------------------
# ANNOTATION: Train the six-class ensemble from top-ranked features and
# apply ordered genotype calling.
# ------------------------------------------------------------------------
{
  genotypes <- names(rankedProteins.prior)

  ## one binary fit per genotype - each with *its own* 50 columns ------------
  fit_fns <- lapply(genotypes, function(gt) {
               feats <- rankedProteins.prior[[gt]][1:50, "feature"]
               fit_APOE_binary(expr   = expr[, feats, drop = FALSE],
                               APOE_gt = APOE_gt,
                               target  = gt,
                               target_ppv = target_ppv,
                               ncores  = ncores,
                               seed    = seed)
             })
  names(fit_fns) <- genotypes

  ## ----------- helper to query *all* probabilities for one sample ----------
  score_all <- function(new_expr_row) {
    sapply(genotypes, function(gt) {
      fn    <- fit_fns[[gt]]
      feats <- attr(fn, "features")
      res   <- fn(new_expr_row[, feats, drop = FALSE])
      attr(res, "prob")
    })
  }

  ## ----------- final 6-class wrapper ---------------------------------------
  predict_APOE6 <- function(new_expr,
                            ncores = parallel::detectCores()) {

    if (is.data.frame(new_expr)) new_expr <- as.matrix(new_expr)

    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(multisession, workers = ncores)

    prob_mat <- t(future.apply::future_apply(
                    new_expr, 1L, future.seed = TRUE,
                    future.packages = c("glmnet", "xgboost", "ranger"),
                    function(row_vec) {
                      mat <- matrix(row_vec, nrow=1, dimnames=list(NULL, colnames(new_expr)))
                      score_all(mat)
                    }))

    res <- character(nrow(prob_mat))
    for (k in seq_len(nrow(prob_mat))) {
      p_row <- prob_mat[k, ]; names(p_row) <- genotypes
      above_thr <- mapply(function(gt, p)
                            p >= attr(fit_fns[[gt]], "thr"),
                          gt = genotypes, p = p_row)
      if (any(above_thr)) {
        margin <- p_row[above_thr] -
                  vapply(genotypes[above_thr],
                         function(gt) attr(fit_fns[[gt]], "thr"), 0)
        res[k] <- names(which.max(margin))
      } else {
        res[k] <- names(which.max(p_row))
      }
    }
    factor(res, levels = genotypes)
  }

  ## ship out: the 6-class wrapper plus access to internals ------------------
  attr(predict_APOE6, "fit_fns")   <- fit_fns
  attr(predict_APOE6, "genotypes") <- genotypes
  predict_APOE6
}

################################################################################
## 3.  20 % hold-out x 100 outer splits (“nested-CV”) -------------------------
################################################################################
APOE6_nestedCV <- function(expr, APOE_gt,
                           rankedProteins.prior,
                           n_outer = 100, hold_frac = 0.20,
                           target_ppv = 0.95,
                           ncores = parallel::detectCores() - 1,
                           seed = 42)
{
  genotypes <- names(rankedProteins.prior)
  outer_res <- vector("list", n_outer)
  set.seed(seed)

  future::plan(multisession, workers = ncores)
  handlers("progress")

  with_progress({
    p <- progressor(steps = n_outer)

    for (i in seq_len(n_outer)) {
      hold_idx  <- caret::createDataPartition(APOE_gt, p = hold_frac,
                                              list = FALSE)[, 1]
      train_idx <- setdiff(seq_len(nrow(expr)), hold_idx)

      ## --- fit on the training 80 % ----------------------------------------
      apoe6 <- train_APOE6_full(expr[train_idx, , drop = FALSE],
                                APOE_gt[train_idx],
                                rankedProteins.prior,
                                target_ppv = target_ppv,
                                ncores = ncores,
                                seed = seed + i)

      ## --- predict on the 20 % hold-out ------------------------------------
      y_hat <- apoe6(expr[hold_idx, , drop = FALSE], ncores = ncores)
      cm <- caret::confusionMatrix(
              y_hat,

# ------------------------------------------------------------------------
# ANNOTATION: Run outer cross-validation with hold-out folds and collect
# confusion matrices/performance metrics.
# ------------------------------------------------------------------------
              factor(APOE_gt[hold_idx], levels = genotypes))

      outer_res[[i]] <- list(
        Accuracy = cm$overall["Accuracy"],
        MacroF1  = mean(cm$byClass[, "F1"], na.rm = TRUE),
        byClass  = as.data.frame(cm$byClass),
        ConfMat  = cm$table)

      p(message = sprintf("outer %3d/%d done", i, n_outer), amount = 1)
    }
  })

  ## ----------- summarise ----------------------------------------------------
  acc   <- sapply(outer_res, `[[`, "Accuracy")
  macro <- sapply(outer_res, `[[`, "MacroF1")

  message(sprintf("\nRepeated 20 %% hold-out (n = %d)\n",
                  n_outer))
  message(sprintf("Accuracy  %.3f ± %.3f",
                  mean(acc), sd(acc)))
  message(sprintf("Macro-F1  %.3f ± %.3f",
                  mean(macro), sd(macro)))

  invisible(outer_res)
}

################################################################################
## 4.  Usage  -----------------------------------------------------------------
################################################################################
## (1)  train the final model on *all* data for real-world deployment on unknown samples
apoe6_full <- train_APOE6_full(t(training.cleanDat.noNA),
                               training.gt.APOE,
                               rankedProteins.prior,
                               target_ppv = 0.50,    # preferred PPV (floor is 0.80 in internal helper function)
                               ncores = 8)

  memLimit=10*1024^3
  options(future.globals.maxSize= memLimit)  #4GB Total size of all global objects that need to be exported - up from 500MB
  Sys.setenv(R_FUTURE_GLOBALS_MAXSIZE=memLimit) #inherited by workers


## (2)  predict training samples
train.APOEpred<-apoe6_full(t(training.cleanDat.noNA), ncores=32)

## Naive statistics; mean accuracy and confusion matrix:
mean(train.APOEpred == training.gt.APOE)
#0.9888196 (original run); rerun: 0.9892262 due to non-determinism
table(train.APOEpred, training.gt.APOE)  #confusion matrix -- should (nearly) match below before rebuild of above function system
#               training.gt.APOE
#train.APOEpred e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#         e2/e4     0     0   351     1     0     0
#         e2/e2    50     0     0     0     0     0
#         e4/e4     0     1     0     4     0   841
#         e2/e3     0  1368     0     1     5     1
#         e3/e4     0     1     4    39  4676     7
#         e3/e3     0    17     5  7307    67    12


## (3)  experimental nested-CV
cv_folds <- APOE6_nestedCV(t(training.cleanDat.noNA),
                           training.gt.APOE,
                           rankedProteins.prior,
                           n_outer = 100,
                           hold_frac = 0.20,
                           target_ppv = 0.50,
                           ncores = 32)

#Accuracy  0.943 ± 0.004
#Macro-F1  0.915 ± 0.013


## (4)  predict unknown samples
#5715 unknown prediction using predictor trained on full noNA data 14758
unknown.APOEpred <- apoe6_full(t(cleanDat.4p13b3[,names(gt.APOE)[which(is.na(gt.APOE))] ]), ncores=8)
table(unknown.APOEpred)
#e2/e4 e2/e2 e4/e4 e2/e3 e3/e4 e3/e3
#  109    18   163   579  1589  3257
names(unknown.APOEpred)<-names(gt.APOE)[which(is.na(gt.APOE))]


## End of genotype prediction, stage 2


######################################################
## Post-prediction integration of genotypes


numericMeta.reg.b345$APOE.mapped<-NA
numericMeta.reg.b345[names(gt.APOE),"APOE.mapped"]<-gt.APOE  #16677 !is.na/22392
library(dplyr)

filtered_data.3097 <- numericMeta.reg.b345 %>%
  filter(
    (Group.withCTimputed == "CT" & APOE.mapped == "e3/e3") |
    (Group.withCTimputed %in% c("AsymAD", "AD", "CT", "MCI", "CI.Other") & APOE.mapped == "e4/e4")
  ) %>%
  group_by(person_id) %>%
  slice_min(order_by = sequential_visit_number, n = 1, with_ties = FALSE) %>%
  ungroup()


numericMeta.reg.b345$APOE.mapped.predicted<-numericMeta.reg.b345$APOE.mapped
numericMeta.reg.b345[names(unknown.APOEpred),"APOE.mapped.predicted"]<-as.character(unknown.APOEpred)
numericMeta.reg.b345$APOE.mapped.predicted[numericMeta.reg.b345$APOE.mapped.predicted==22]<-"e2/e2"
numericMeta.reg.b345$APOE.mapped.predicted[numericMeta.reg.b345$APOE.mapped.predicted==23]<-"e2/e3"
numericMeta.reg.b345$APOE.mapped.predicted[numericMeta.reg.b345$APOE.mapped.predicted==24]<-"e2/e4"
numericMeta.reg.b345$APOE.mapped.predicted[numericMeta.reg.b345$APOE.mapped.predicted==33]<-"e3/e3"
numericMeta.reg.b345$APOE.mapped.predicted[numericMeta.reg.b345$APOE.mapped.predicted==34]<-"e3/e4"
numericMeta.reg.b345$APOE.mapped.predicted[numericMeta.reg.b345$APOE.mapped.predicted==44]<-"e4/e4"

filtered_data.plus5715 <- numericMeta.reg.b345 %>%
  filter(
    (Group.withCTimputed == "CT" & APOE.mapped.predicted == "e3/e3") |
    (Group.withCTimputed %in% c("AsymAD", "AD", "CT", "MCI", "CI.Other") & APOE.mapped.predicted == "e4/e4")
  ) %>%
  group_by(person_id) %>%
  slice_min(order_by = sequential_visit_number, n = 1, with_ties = FALSE) %>%
  ungroup()

table(filtered_data.3097$APOE.mapped, filtered_data.3097$Group.withCTimputed)

table(filtered_data.plus5715$APOE.mapped.predicted, filtered_data.plus5715$Group.withCTimputed)


filtered_data.3097.last <- numericMeta.reg.b345 %>%
  filter(
    (Group.withCTimputed == "CT" & APOE.mapped == "e3/e3") |
    (Group.withCTimputed %in% c("AsymAD", "AD", "CT", "MCI", "CI.Other") & APOE.mapped == "e4/e4")
  ) %>%
  group_by(person_id) %>%
  slice_max(order_by = sequential_visit_number, n = 1, with_ties = FALSE) %>%
  ungroup()


filtered_data.plus5715.last <- numericMeta.reg.b345 %>%
  filter(
    (Group.withCTimputed == "CT" & APOE.mapped.predicted == "e3/e3") |
    (Group.withCTimputed %in% c("AsymAD", "AD", "CT", "MCI", "CI.Other") & APOE.mapped.predicted == "e4/e4")
  ) %>%
  group_by(person_id) %>%
  slice_max(order_by = sequential_visit_number, n = 1, with_ties = FALSE) %>%
  ungroup()

table(filtered_data.3097.last$APOE.mapped, filtered_data.3097.last$Group.withCTimputed)

table(filtered_data.plus5715.last$APOE.mapped.predicted, filtered_data.plus5715.last$Group.withCTimputed)

dim(filtered_data.3097.last)
dim(filtered_data.plus5715.last)


filtered_data.plus5715.33inclAD <- numericMeta.reg.b345 %>%
  filter(
    (Group.withCTimputed %in% c("AsymAD", "AD", "CT", "MCI", "CI.Other") & APOE.mapped.predicted == "e3/e3") |
    (Group.withCTimputed %in% c("AsymAD", "AD", "CT", "MCI", "CI.Other") & APOE.mapped.predicted == "e4/e4")
  ) %>%
  group_by(person_id) %>%
  slice_min(order_by = sequential_visit_number, n = 1, with_ties = FALSE) %>%
  ungroup()


table(filtered_data.plus5715.33inclAD$APOE.mapped.predicted, filtered_data.plus5715.33inclAD$Group.withCTimputed)

dim(filtered_data.plus5715.33inclAD)

#regress 22392 for site, protecting age, sex, and APOE4 carrier status
#subset first sequential visits 3177 (e3/e3 CT only) and 4199 (same Dx's as e4/e4) for Shijia's analysis
#build network on 22392 regressed/protected cleanDat
#deliver MEs subset to same as above for Shijia

save.image("saved.image-genotype_prediction_finalized_2026rerun.RData") # does not include cv_folds var with 100 outer CV fold statistics


# Extract the "byClass" data frames into a new list
byClass_list <- lapply(cv_folds, function(x) x$byClass)


# ------------------------------------------------------------------------
# ANNOTATION: Save the finalized genotype prediction workspace after
# retraining on known genotype samples.
# ------------------------------------------------------------------------
# Optionally ensure all are data frames and of the same dimensions
# Convert each to a numeric matrix
numeric_matrices <- lapply(byClass_list, function(df) {
  as.matrix(suppressWarnings(apply(df, 2, as.numeric)))  # ensure numeric, suppress coercion warnings
})

# Check dimensions are consistent
stopifnot(all(sapply(numeric_matrices, function(x) all(dim(x) == c(6, 11)))))

# Stack into a 3D array: 6 x 11 x 100
array_data <- simplify2array(numeric_matrices)

# Calculate mean and standard deviation across the third dimension
mean_matrix <- apply(array_data, c(1, 2), mean, na.rm = TRUE)
sd_matrix <- apply(array_data, c(1, 2), sd, na.rm = TRUE)

saveRDS(cv_folds,"genotype_prediction_100outerCVfolds-stats_2026rerun.RDS") # on Windows03


## Pipeline continues with final regression pass on 22392 total samples with age, sex, and mapped+predicted APOE genotype



# ------------------------------------------------------------------------
# ANNOTATION: Persist 100-outer-fold summary statistics for manuscript
# tables and reproducibility.
# ------------------------------------------------------------------------
###################################################################################
## Summarize confusion matrix on 100x 20% outer CV folds (test data only: cv_folds)

summarize_cv_confmat <- function(cv_folds, genotypes = NULL) {

  # pull the per-fold confusion tables
  cms <- lapply(cv_folds, function(x) x$ConfMat)
  cms <- Filter(Negate(is.null), cms)
  if (!length(cms)) stop("No ConfMat tables found in cv_folds.")

  # establish genotype order
  if (is.null(genotypes)) {
    genotypes <- colnames(as.matrix(cms[[1]]))
    if (is.null(genotypes)) stop("Could not infer genotypes from ConfMat dimnames.")
  }

  K <- length(genotypes)
  N <- length(cms)

  # 3D arrays: counts and per-fold column-wise percents
  count_arr <- array(0, dim = c(K, K, N),
                     dimnames = list(pred = genotypes, truth = genotypes, fold = seq_len(N)))
  pct_arr   <- array(NA_real_, dim = c(K, K, N),
                     dimnames = list(pred = genotypes, truth = genotypes, fold = seq_len(N)))

  for (i in seq_len(N)) {
    m <- as.matrix(cms[[i]])

    # enforce full KxK layout (fill missing with 0)
    m_full <- matrix(0, nrow = K, ncol = K, dimnames = list(genotypes, genotypes))
    rr <- intersect(rownames(m), genotypes)
    cc <- intersect(colnames(m), genotypes)
    m_full[rr, cc] <- m[rr, cc, drop = FALSE]

    count_arr[, , i] <- m_full

    # column-wise percent of ground-truth totals within this fold
    col_tot <- colSums(m_full)
    pct <- sweep(m_full, 2, col_tot, FUN = "/") * 100
    pct[, col_tot == 0] <- NA_real_  # absent truth genotype in this fold
    pct_arr[, , i] <- pct
  }

  # ---- mean/sd in percent space (mean-of-folds) ----
  mean_pct <- apply(pct_arr,   c(1, 2), mean, na.rm = TRUE)
  sd_pct   <- apply(pct_arr,   c(1, 2), sd,   na.rm = TRUE)
  n_valid_pct <- apply(pct_arr, c(1, 2), function(z) sum(is.finite(z)))

  # ---- mean/sd of raw counts (mean-of-folds) ----
  mean_count <- apply(count_arr, c(1, 2), mean)
  sd_count   <- apply(count_arr, c(1, 2), sd)

  # ---- pooled aggregate confusion matrix (sum counts across folds) ----
  pooled_count <- apply(count_arr, c(1, 2), sum)

  # pooled column-wise percent (optional, but usually what you want to look at)

# ------------------------------------------------------------------------
# ANNOTATION: Summarize and plot confusion matrices across folds and full-
# data predictions.
# ------------------------------------------------------------------------
  pooled_col_tot <- colSums(pooled_count)
  pooled_pct <- sweep(pooled_count, 2, pooled_col_tot, FUN = "/") * 100
  pooled_pct[, pooled_col_tot == 0] <- NA_real_

  # return as matrices + also as data.frames (requested)
  list(
    ## per-fold percent summaries
    mean_pct = mean_pct,
    sd_pct   = sd_pct,
    n_valid_pct = n_valid_pct,
    pct_arr  = pct_arr,

    ## per-fold raw-count summaries (no percent conversion)
    mean_count = mean_count,
    sd_count   = sd_count,
    mean_count_df = as.data.frame(mean_count),
    sd_count_df   = as.data.frame(sd_count),

    ## pooled aggregates
    pooled_count = pooled_count,
    pooled_count_df = as.data.frame(pooled_count),
    pooled_pct = pooled_pct,
    pooled_pct_df = as.data.frame(pooled_pct),

    ## raw counts array if you ever want to debug / slice by fold
    count_arr = count_arr
  )
}

## ---- usage ---------------------------------------------------------------
genotypes <- c("e2/e2","e2/e3","e2/e4","e3/e3","e3/e4","e4/e4")
cv_sum <- summarize_cv_confmat(cv_folds, genotypes = genotypes)

round(cv_sum$mean_pct, 2)
round(cv_sum$sd_pct, 2)

# pooled (aggregate) versions
cv_sum$pooled_count
round(cv_sum$pooled_pct, 2)

# mean/sd of integer counts as data frames
cv_sum$mean_count_df
cv_sum$sd_count_df


#####################
macroF1_from_cv_folds <- function(cv_folds, genotypes = NULL, recompute_if_needed = TRUE) {
  # Helper: compute macro-F1 from a confusion matrix table
  macroF1_from_table <- function(tab, genotypes) {
    m <- as.matrix(tab)
    # enforce full square layout
    K <- length(genotypes)
    mf <- matrix(0, nrow = K, ncol = K, dimnames = list(genotypes, genotypes))
    rr <- intersect(rownames(m), genotypes)
    cc <- intersect(colnames(m), genotypes)
    mf[rr, cc] <- m[rr, cc, drop = FALSE]

    # per-class precision/recall/F1 (one-vs-all, by class label)
    tp <- diag(mf)
    fp <- rowSums(mf) - tp
    fn <- colSums(mf) - tp

    precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_)
    recall    <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)
    f1        <- ifelse(is.finite(precision) & is.finite(recall) & (precision + recall) > 0,
                        2 * precision * recall / (precision + recall),
                        NA_real_)
    mean(f1, na.rm = TRUE)
  }

  # infer genotypes if needed
  if (is.null(genotypes)) {
    first_cm <- NULL
    for (x in cv_folds) {
      if (!is.null(x$ConfMat)) { first_cm <- x$ConfMat; break }
      if (!is.null(x$byClass)) { # can infer from byClass rownames
        rn <- rownames(x$byClass)
        if (!is.null(rn)) { genotypes <- rn; break }
      }
    }
    if (is.null(genotypes) && !is.null(first_cm)) genotypes <- colnames(as.matrix(first_cm))
    if (is.null(genotypes)) stop("Could not infer genotypes; please supply genotypes=.")
  }

  macro_vec <- vapply(seq_along(cv_folds), function(i) {
    x <- cv_folds[[i]]

    # 1) if you stored MacroF1, use it
    if (!is.null(x$MacroF1) && is.finite(as.numeric(x$MacroF1))) {
      return(as.numeric(x$MacroF1))
    }

    # 2) else, if byClass exists with F1 column
    if (!is.null(x$byClass)) {
      bc <- x$byClass
      if (is.data.frame(bc) && "F1" %in% colnames(bc)) {
        return(mean(bc[, "F1"], na.rm = TRUE))
      }
    }

    # 3) else, recompute from ConfMat if requested
    if (recompute_if_needed && !is.null(x$ConfMat)) {
      return(macroF1_from_table(x$ConfMat, genotypes))
    }

    NA_real_
  }, numeric(1))

  list(
    per_fold = macro_vec,
    mean = mean(macro_vec, na.rm = TRUE),
    sd   = sd(macro_vec, na.rm = TRUE),
    n_valid = sum(is.finite(macro_vec)),
    genotypes = genotypes
  )
}

## ---- usage ---------------------------------------------------------------
genotypes <- c("e2/e4","e2/e2","e4/e4","e2/e3","e3/e4","e3/e3")
mf1 <- macroF1_from_cv_folds(cv_folds, genotypes = genotypes)

mf1$mean
mf1$sd
mf1$n_valid
# mf1$per_fold  # vector length 100 (with NAs possible if a fold couldn't compute)


##############
summarize_accuracy_from_cv_sum <- function(cv_sum, genotypes = NULL) {
  pct_arr   <- cv_sum$pct_arr
  count_arr <- cv_sum$count_arr

  if (length(dim(pct_arr)) != 3) stop("cv_sum$pct_arr must be a 3D array [pred, truth, fold].")
  if (length(dim(count_arr)) != 3) stop("cv_sum$count_arr must be a 3D array [pred, truth, fold].")
  if (dim(pct_arr)[3] != dim(count_arr)[3]) stop("pct_arr and count_arr must have same number of folds.")

  if (is.null(genotypes)) {
    genotypes <- dimnames(pct_arr)$pred
    if (is.null(genotypes)) stop("Could not infer genotypes from cv_sum$pct_arr dimnames; please supply genotypes=.")
  }

  K <- length(genotypes)
  N <- dim(pct_arr)[3]

  ## ---------------------------
  ## 1) Per-genotype % accuracy (diagonal of column-normalized % matrix)
  ##    This is per-genotype recall/sensitivity: P(pred=gt | truth=gt)
  ## ---------------------------
  per_gt_pct_by_fold <- vapply(seq_len(N), function(i) {
    diag(as.matrix(pct_arr[, , i]))   # <-- key fix
  }, numeric(K))

  rownames(per_gt_pct_by_fold) <- genotypes
  colnames(per_gt_pct_by_fold) <- paste0("fold", seq_len(N))

  per_genotype <- data.frame(
    genotype = genotypes,
    mean_pct_accuracy = rowMeans(per_gt_pct_by_fold, na.rm = TRUE),
    sd_pct_accuracy   = apply(per_gt_pct_by_fold, 1, sd, na.rm = TRUE),
    n_valid_folds     = rowSums(is.finite(per_gt_pct_by_fold)),
    row.names = NULL
  )

  ## ---------------------------
  ## 2) Overall accuracy per fold from raw counts
  ## ---------------------------
  overall_acc_by_fold <- vapply(seq_len(N), function(i) {
    m <- as.matrix(count_arr[, , i])
    tot <- sum(m)
    if (tot == 0) return(NA_real_)
    sum(diag(m)) / tot * 100
  }, numeric(1))

  overall <- list(
    per_fold_pct  = overall_acc_by_fold,
    mean_pct      = mean(overall_acc_by_fold, na.rm = TRUE),
    sd_pct        = sd(overall_acc_by_fold, na.rm = TRUE),
    n_valid_folds = sum(is.finite(overall_acc_by_fold))
  )

  ## (optional) pooled overall accuracy across all folds (not the same as mean-of-folds)
  pooled <- NULL
  if (!is.null(cv_sum$pooled_count)) {
    pc <- as.matrix(cv_sum$pooled_count)
    pooled <- list(overall_pct = sum(diag(pc)) / sum(pc) * 100)
  }

  list(
    per_genotype = per_genotype,
    per_genotype_pct_by_fold = per_gt_pct_by_fold,
    overall = overall,
    pooled = pooled
  )
}

## ---- usage ---------------------------------------------------------------
acc_sum <- summarize_accuracy_from_cv_sum(cv_sum, genotypes = genotypes)

acc_sum$per_genotype
acc_sum$overall$mean_pct
acc_sum$overall$sd_pct
