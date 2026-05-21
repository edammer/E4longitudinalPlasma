##############################################################################
# Pipeline annotation header: 2a.APOE_prediction-stage_1.R
# Manuscript code section(s): 2
#
# Purpose:
# Stage 1 APOE genotype machine-learning feature selection: train one-vs-all
# genotype classifiers and rank proteins by ensemble feature importance.
#
# Principal inputs:
#   - Regression-adjusted log2 protein abundance matrix
#   - Harmonized APOE genotype labels for samples with known genotype
#
# Principal outputs:
#   - rankedProteins / rankedProteins.prior objects with top-ranked features
#     per genotype
#
# Step overview:
#   1. Define fit_APOE_binary.95Acc(), an ensemble wrapper around glmnet,
#      xgboost, and ranger.
#   2. Within repeated cross-validation, z-scale features from the training
#      fold and use fold-specific scaling on held-out data.
#   3. Estimate genotype-specific decision thresholds to meet target positive
#      predictive value with a minimum probability floor.
#   4. Fit final all-data learners for each genotype and aggregate scaled
#      feature importances across algorithms.
#   5. Rank features separately for e2/e2, e2/e3, e2/e4, e3/e3, e3/e4, and
#      e4/e4.
#
# Notes added during manuscript-code alignment:
#   - Original executable code was preserved except for whitespace/encoding cleanup
#     and the explicitly marked non-run archive described below when applicable.
#   - Comments inserted by this pass are prefixed with "ANNOTATION:".
##############################################################################

# ------------------------------------------------------------------------
# ANNOTATION: Set the working directory and load ML packages used by the
# APOE stage 1 feature-selection ensemble.
# ------------------------------------------------------------------------
setwd("Z:/EBD/grid/4p13b3forAPOEpredict+2ndRegrAgain/26predict/")
#load("z:/EBD/4p13b4b.SecondRegressionsComplete.19sites_Fsplit+APOEgenoIMPUTATION.RData") #13.7 GB

# Prompt: Given 6 lists of features (log2(protein relative abundances) which predict each of 6 APOE genotypes (e2/e2, e2/e3, e2/e4, e3/e3, e3/e4, and e4/e4) in a binary fashion with 97 percent accuracy in a large (15000 samples) test data set subjected to the attached ML function, I want to devise a strategy (adapted algorithm) and implement it it R code to assemble the most accurate 6-genotype predictor. Because the binary predictors were most accurate for e2/e4 and e2/e2, I suggest predicting e2/e4 first, then e2/e2, then e4/e4, followed by e2/e3, e3/4, and e3/e3. Given the known accuracy of individual binary predictors all 97 percent or better, suggest what the accuracy of the 6-genotype prediction algorithm could be, at best.
library(caret)  #ADDI Windows03 VM, R v4.1.3
library(glmnet)
library(xgboost)
library(ranger)
library(progressr)
library(dplyr)
library(future)
library(doFuture)
library(doRNG)

## Prevent RStudio-specific error at very end of function ("error in evaluating the argument 'x' in selecting a method for function 'print': object 'spaces' not found")
#if (!exists("spaces", envir = asNamespace("cli")))
#    assign("spaces", c("", vapply(1:20, strrep, "", x = " ")),
#           envir = asNamespace("cli"))
## Just don't run in RStudio.

fit_APOE_binary.95Acc <- function(expr, APOE_gt,
                                  nfold = 5, nrep = 5,
                                  ncores = parallel::detectCores() - 1,
                                  target=c("e4/e4"),   # run one genotype vs all others
                                  target_ppv=0.95,
                                  seed   = 1) {

# ------------------------------------------------------------------------
# ANNOTATION: Define the repeated-CV binary classifier wrapper used for one-
# vs-all genotype feature ranking.
# ------------------------------------------------------------------------

  set.seed(seed)

  if (!"package:cli" %in% search()) suppressPackageStartupMessages(library(cli))

  memLimit=4*1024^3
  options(future.globals.maxSize= memLimit)  #4GB Total size of all global objects that need to be exported - up from 500MB
  Sys.setenv(R_FUTURE_GLOBALS_MAXSIZE=memLimit) #inherited by workers

  ## Helper function - threshold minimum
  pick_thr <- function(prob, truth, target_ppv = 0.95,
                       min_tp = 30, floor = 0.80) {

# ------------------------------------------------------------------------
# ANNOTATION: Select hard-call thresholds from held-out probabilities to
# target high PPV with a minimum probability floor.
# ------------------------------------------------------------------------
    ok  <- !is.na(prob)                       # drop rows with NA prob
    prob <- prob[ok]
    truth <- truth[ok]

    if (sum(truth) == 0)  return(1)           # no positives
    if (all(prob == 0))   return(1)           # degenerate

    o   <- order(prob, decreasing = TRUE)
    tp  <- cumsum(truth[o] == 1)
    fp  <- cumsum(truth[o] == 0)
    ppv <- tp / (tp + fp)
    ok1  <- which(ppv >= target_ppv & tp >= min_tp)
    thr <- if (length(ok1)) prob[o[max(ok1)]] else quantile(prob[truth==1], .90, na.rm=TRUE)
    max(thr, floor)
  }


# ------------------------------------------------------------------------
# ANNOTATION: Prepare expression matrices using fold-specific z-scaling and
# missing-value handling.
# ------------------------------------------------------------------------
  ## ------------------------------------------------------------------------
  ## 0.  create resampling indices once
  ## ------------------------------------------------------------------------
  cvIndex <- createMultiFolds(APOE_gt, k = nfold, times = nrep)
  nTasks  <- length(cvIndex)
  ## ------------------------------------------------------------------------
  ## 1.  start the cluster
  ## ------------------------------------------------------------------------
  handlers("progress")
  handlers(global=TRUE)  # set handler for progress bar before the cluster

  plan(multisession, workers = ncores)  #, globals.maxSize=memLimit)  # or multicore on linux/macOS
  registerDoFuture()

  ## ------------------------------------------------------------------------
  ## 2.  run each fold in a worker  -----------------------------------------
  ## ------------------------------------------------------------------------

  with_progress({                              # << all progress lives here
    n_sub <- 3   # Progress bar increments 3x per fold
    p <- progressor(steps=length(cvIndex) * n_sub )     #along = cvIndex)          # one step per fold

    p(sprintf("Initializing %d workers...", ncores), amount=0)


  metrics <- foreach(fold = seq_along(cvIndex),
                     .combine   = rbind,

# ------------------------------------------------------------------------
# ANNOTATION: Train glmnet, xgboost, and ranger learners inside repeated
# cross-validation folds.
# ------------------------------------------------------------------------
                     .options.future = list(  #expr not exported explicitly; captured only once.
                        globals = list(APOE_gt  = APOE_gt,
                                       pick_thr = pick_thr,
                                       target   = target,
                                       seed     = seed)
                     ),
                     .export = "p",            # let workers see 'p'
                     .packages  = c("progressr","glmnet","xgboost","ranger","dplyr","cli")) %dorng% {

    ## ----------------- announce fold start -----------------------------
    p(sprintf("fold %d/%d  -  started", fold, length(cvIndex)), amount=0)

    set.seed(seed + fold)                        # reproducible inside worker

    tr <- cvIndex[[fold]]
    te <- setdiff(seq_len(nrow(expr)), tr)

    tr_idx <- cvIndex[[fold]]
    te_idx <- setdiff(seq_len(nrow(expr)), tr_idx)

    ## --------- preprocessing  (Z scale, no missing data in input) -------
    prep <- function(m) scale(m[, colMeans(is.na(m)) <= .20, drop = FALSE])
    X_tr <- prep(expr[tr, ])
    X_te <- scale(expr[te, colnames(X_tr)],
                  center = attr(X_tr,"scaled:center"),
                  scale  = attr(X_tr,"scaled:scale"))
    X_te[is.na(X_te)] <- 0

    y_bin_tr <- factor(ifelse(APOE_gt[tr] == target, "pos", "neg"))
    y_bin_te <- factor(ifelse(APOE_gt[te] == target, "pos", "neg"))
    w_bin    <- ifelse(y_bin_tr == "pos", 8, 1)   # keep your weight idea

    # detect folds without positives OR without negatives
    if (length(unique(y_bin_tr)) < 2) {
      p(sprintf("fold %d/%d  •  skipped (only one class)", fold,
                length(cvIndex)), amount = 3)             # step the bar
      return(data.frame(Precision = NA, Recall = NA, Fold = fold))
    }

    # ----------   learner 1   glmnet  -------------------------------------
    glm_cv <- cv.glmnet(X_tr, y_bin_tr, family = "binomial",
                        weights = w_bin, type.measure = "class")
    p_glm <-    drop(predict(glm_cv, X_te, s = "lambda.min", type = "response"))
    p(sprintf("fold %d/%d  •  glmnet done", fold, length(cvIndex)), amount=1)

    # ----------   learner 2   xgboost  ------------------------------------
    dtr <- xgb.DMatrix(X_tr, label = as.numeric(y_bin_tr) - 1, weight = w_bin)
    dte <- xgb.DMatrix(X_te, label = as.numeric(y_bin_te) - 1)
    xpar <- list(objective = "binary:logistic", eta = 0.1,
                 max_depth = 6, subsample = 0.8, colsample_bytree = 0.8,
                 nthread = 1, eval_metric = "logloss")
    xgb <- xgb.train(xpar, dtr, watchlist=list(train=dtr, eval=dte), nrounds = 200,
                     verbose = 0, early_stopping_rounds = 20)
    p_xgb <- drop(predict(xgb, dte))
    p(sprintf("fold %d/%d  •  xgboost done", fold, length(cvIndex)), amount=1)

    # ----------   learner 3   ranger   ------------------------------------
    rf  <- ranger(y_bin_tr ~ ., data = data.frame(y_bin_tr, X_tr),
                  probability   = TRUE,
                  num.trees     = 500, num.threads = 1,
                  class.weights = c(neg = 1, pos = 8))
    p_rf <- predict(rf, data.frame(X_te))$predictions[,"pos"]
    p(sprintf("fold %d/%d  •  ranger done", fold, length(cvIndex)), amount=1)

    # ----------   averaged prob & hard threshold --------------------------
    # probabilities on training set
    p_glm_tr <- drop(predict(glm_cv, X_tr, s = "lambda.min", type = "response"))
    p_xgb_tr <- drop(predict(xgb, dtr))
    p_rf_tr  <-        predict(rf, data.frame(X_tr))$predictions[,"pos"]
    p_avg_tr <- (p_glm_tr + p_xgb_tr + p_rf_tr) / 3    # length == length(y_bin_tr)

    thr   <- pick_thr(p_avg_tr, y_bin_tr == "pos", target_ppv)   # use helper

    # probabilities on test set
    p_avg_te <- (p_glm + p_xgb + p_rf) / 3  # already computed on X_te
    pred  <- factor(ifelse(p_avg_te >= thr, target, NA_character_),
                    levels = c(target))

    tp <- sum(pred == target & y_bin_te == "pos", na.rm = TRUE)
    fp <- sum(pred == target & y_bin_te == "neg", na.rm = TRUE)
    fn <- sum(pred != target & y_bin_te == "pos", na.rm = TRUE)

    prec <- if ((tp+fp) > 0) tp/(tp+fp) else NA_real_
    rec  <- if ((tp+fn) > 0) tp/(tp+fn) else NA_real_

    data.frame(Precision = prec, Recall = rec, Fold = fold)
  } # foreach
  }) # with_progress

#  stopCluster(cl)                               # tidy up
#  registerDoSEQ()                               # back to sequential

  cat(sprintf("CV %s - Precision %.3f ± %.3f | Recall %.3f ± %.3f\n\n",
              target,
              mean(metrics$Precision, na.rm=TRUE), sd(metrics$Precision, na.rm=TRUE),
              mean(metrics$Recall,    na.rm=TRUE), sd(metrics$Recall,    na.rm=TRUE)))

  if (!exists("binaryPredictionMetrics", envir=.GlobalEnv)) assign("binaryPredictionMetrics",list(), envir=.GlobalEnv)
  binaryPredictionMetrics[[target]] <<- data.frame(Precision=metrics$Precision, Recall=metrics$Recall)

  ## ------------------------------------------------------------------------
  ## 3.  fit final ensemble on all data  (sequential) -----------------------
  ## ------------------------------------------------------------------------
  prep_expr <- function(mat) scale(mat[, colMeans(is.na(mat)) <= .20, drop = FALSE])
  X_all <- prep_expr(expr)
  y_all <- APOE_gt                                   # keep as character

  with_progress({
    p2 <- progressor(steps=length(target)*3)
    p2(sprintf("Starting binary fit of final ensemble on all data (serial/non-parallel)..."), amount=0)

    ## ----------  each genotype specified in fallback-targets, vs rest  -------
    ovr <- list()

    for (tg in target) {

      y_bin  <- factor(ifelse(y_all == tg, "pos", "neg"))
      w_pos  <- if (tg %in% c("e2/e2","e2/e4","e4/e4")) 12 else 8  # higher weight for rarer genotypes
      w_bin  <- ifelse(y_bin == "pos", w_pos, 1)

      ## ---- (i) glmnet --------------------------------------------------------
      glm_tg <- cv.glmnet(X_all, y_bin, family = "binomial",
                          weights = w_bin, type.measure = "class")
      p2(sprintf("GLM fit on full data  •  finished"), amount=1)

      ## ---- (ii) xgboost ------------------------------------------------------
      d_bin  <- xgb.DMatrix(X_all, label = as.numeric(y_bin) - 1, weight = w_bin)
      xpar_b <- list(objective = "binary:logistic", eta = 0.1,
                     max_depth = 6, subsample = 0.8,
                     colsample_bytree = 0.8, nthread = ncores, eval_metric = "logloss")
      xgb_tg <- xgb.train(xpar_b, d_bin, watchlist=list(train=d_bin), nrounds = 200,
                          verbose = 0, early_stopping_rounds = 20)
      p2(sprintf("XGboost fit on full data  •  finished"), amount=1)

      ## ---- (iii) ranger ------------------------------------------------------
      rf_tg  <- ranger(y_bin ~ ., data = data.frame(y_bin, X_all),
                       probability = TRUE, num.trees = 500,
                       class.weights = c("neg" = 1, "pos" = w_pos),
                       num.threads  = ncores)
      p2(sprintf("Random Forest fit on full data  •  finished"), amount=1)

      ## ---- averaged prob on the *training* set ------------------------------
      p_glm <-      drop(predict(glm_tg, X_all, s = "lambda.min", type="response"))
      p_xgb <-      drop(predict(xgb_tg, d_bin))
      p_rf  <- rf_tg$predictions[,"pos"]

      p_avg <- (p_glm + p_xgb + p_rf) / 3

      ## learn a threshold that guarantees >= target_ppv on the 15 k labelled rows
      thr   <- pick_thr(p_avg, y_bin == "pos", target_ppv)

      ovr[[tg]] <- list(glm = glm_tg, xgb = xgb_tg, rf = rf_tg, thr = thr)
    }

    p2(sprintf("Binary fit on full data  •  finished"), amount=0)


## ---------------- feature importance -------------------------------- ###
    imp_glm <- abs(as.matrix(coef(ovr[[target]]$glm, s="lambda.min"))[-1,1])
    imp_xgb <- {
         g <- xgb.importance(model = ovr[[target]]$xgb)
         setNames(g$Gain, g$Feature)
    }
    imp_rf  <- ovr[[target]]$rf$variable.importance

    ## put everything on the same scale and average
    all_feats <- union(names(imp_glm), union(names(imp_xgb), names(imp_rf)))
    top_k=length(all_feats)

    imp_mat   <- cbind(
       glm = imp_glm[all_feats], xgb = imp_xgb[all_feats], rf = imp_rf[all_feats])
    imp_mat[is.na(imp_mat)] <- 0

# ------------------------------------------------------------------------
# ANNOTATION: Aggregate held-out predictions and summarize binary
# performance for each genotype target.
# ------------------------------------------------------------------------
    imp_scaled <- scale(imp_mat)
    imp_mean   <- rowMeans(imp_scaled)
    top_feats  <- sort(imp_mean, decreasing = TRUE) #[1:top_k]

  }) # with_progress

  ## -------- export top_k predictive proteins --------------------------
    if (!exists("rankedProteins", envir=.GlobalEnv)) assign("rankedProteins",list(), envir=.GlobalEnv)
    rankedProteins[[target]] <<- data.frame(feature=names(sort(top_feats, decreasing = TRUE)), importance=sort(top_feats, decreasing = TRUE))
    names(top_feats)<-gsub("\\|","_",names(top_feats))
#    print(knitr::kable(data.frame(Protein = names(top_feats)[c(1:10,(top_k-9):top_k)],
#                                  Importance = round(top_feats,3)[c(1:10,(top_k-9):top_k)]),
#                       caption = sprintf("Top 10 and bottom 10 (of top %d) predictive proteins - (all exported to list rankedProteins)", top_k)))
  ### ------------------------------------------------------------------- ###


  ## ------------------------------------------------------------------------
  ## 4.  prediction wrapper -------------------------------------------------
  ## ------------------------------------------------------------------------
    ## helper - probability only ---------------------------------------------
     prob_fun <- function(new_expr){
       # keep exactly the columns the model was trained with ------------------
       new_expr <- as.matrix(new_expr)[ , colnames(X_all), drop = FALSE]
       # use the same centering/scaling that is stored inside X_all ----------
       new_expr <- scale(new_expr,
                         center = attr(X_all, "scaled:center"),
                         scale  = attr(X_all, "scaled:scale"))
       new_expr[is.na(new_expr)] <- 0

       p_g <-      drop(predict(ovr[[target]]$glm, new_expr,
                                s = "lambda.min", type = "response"))
      p_x <-      drop(predict(ovr[[target]]$xgb, new_expr))
       p_r <- predict(ovr[[target]]$rf , data.frame(new_expr))$predictions[,"pos"]

       (p_g + p_x + p_r) / 3                # numeric vector of probabilities
     }

    ## main wrapper - returns factor but *carries* the probability -----------
     predict_wrapper <- function(new_expr){
       p_bin <- prob_fun(new_expr)
       pred  <- factor(ifelse(p_bin >= thr, target, NA_character_),
                       levels = target)
       attr(pred, "prob") <- p_bin          # <-- hand the prob back to caller
       pred
     }

     ## expose helper & threshold on the function object itself ---------------
     attr(predict_wrapper, "prob") <- prob_fun
     attr(predict_wrapper, "thr")  <- thr
     ## ---- remember the 50 features that were used to train this model ----
     attr(predict_wrapper, "features") <- colnames(X_all) #rankedProteins.prior[[target]][1:50, "feature"]

     predict_wrapper
}

# The above function is used to generate the 6 genotype-specific feature importance lists with all 7334 assays as input, and target ppv 0.95

# ------------------------------------------------------------------------
# ANNOTATION: Fit final all-data learners and combine algorithm-specific
# importances into ensemble feature ranks.
# ------------------------------------------------------------------------
#  (see code after genotype accounting below):
#binaryPredictionMetrics<-readRDS("z:/EBD/binaryPredictionMetrics.list.RDS")
#rankedProteins<-readRDS("z:/EBD/rankedProteins.list.RDS")
## (Don't have these yet): binaryPredictionMetrics, rankedProteins


#############################################################
##  Genotype cleanup / Accounting of ground truth (gt)
table(numericMeta.reg.b345$APOE)
#22   23  24   33   34   44
#60 1492 347 7587 4657  780
table(gt.APOE)
#22   23  24   33   34   44
#60 1494 359 7671 4695  784


# BUILD gt.APOE vector from scratch (some discrepancies found vs. numericMetareg.b345$APOE (without any mapping)  -- avoid mapping our 3 cohorts yet
gt.APOE.old<-gt.APOE

gt.APOE<-numericMeta.reg.b345$APOE
names(gt.APOE)<-rownames(numericMeta.reg.b345)

gt.APOE[gt.APOE==22]<-"e2/e2"
gt.APOE[gt.APOE==23]<-"e2/e3"
gt.APOE[gt.APOE==24]<-"e2/e4"
gt.APOE[gt.APOE==33]<-"e3/e3"
gt.APOE[gt.APOE==34]<-"e3/e4"
gt.APOE[gt.APOE==44]<-"e4/e4"

gt.APOE.noMapping<-gt.APOE


### MAPPED APOE genotypes from our 3 cohorts
#UDS.map
UDS.map.APOE<-UDS.map
UDS.map.APOE$APOE.4cohort<-as.numeric(gsub("e","",traits.4cohort$APOE[match(UDS.map[,2],traits.4cohort$LoadedSampleName)]))
#UDS.map.APOE$APOE.predict<-imputed.APOE.b3.all.noNA[match(UDS.map.APOE[,1],names(imputed.APOE.b3.all.noNA))]

UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE$APOE.4cohort==22)]<-"e2/e2"
UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE$APOE.4cohort==23)]<-"e2/e3"
UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE$APOE.4cohort==33)]<-"e3/e3"
UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE$APOE.4cohort==24)]<-"e2/e4"
UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE$APOE.4cohort==34)]<-"e3/e4"
UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE$APOE.4cohort==44)]<-"e4/e4"

table(paste0(UDS.map.APOE$APOE.4cohort,".",UDS.map.APOE$APOE.predict))


BH.map.APOE<-BH.map
BH.map.APOE$APOE.4cohort<-as.numeric(gsub("e","",traits.4cohort$APOE[match(BH.map[,2],traits.4cohort$LoadedSampleName)]))
#BH.map.APOE$APOE.predict<-imputed.APOE.b3.all.noNA[match(BH.map.APOE[,1],names(imputed.APOE.b3.all.noNA))]

BH.map.APOE$APOE.4cohort[which(BH.map.APOE$APOE.4cohort==22)]<-"e2/e2"
BH.map.APOE$APOE.4cohort[which(BH.map.APOE$APOE.4cohort==23)]<-"e2/e3"
BH.map.APOE$APOE.4cohort[which(BH.map.APOE$APOE.4cohort==33)]<-"e3/e3"
BH.map.APOE$APOE.4cohort[which(BH.map.APOE$APOE.4cohort==24)]<-"e2/e4"
BH.map.APOE$APOE.4cohort[which(BH.map.APOE$APOE.4cohort==34)]<-"e3/e4"
BH.map.APOE$APOE.4cohort[which(BH.map.APOE$APOE.4cohort==44)]<-"e4/e4"

table(paste0(BH.map.APOE$APOE.4cohort,".",BH.map.APOE$APOE.predict))


RM.map.APOE<-RM.map
RM.map.APOE$APOE.4cohort<-as.numeric(gsub("e","",traits.4cohort$APOE[match(RM.map[,2],traits.4cohort$LoadedSampleName)]))
#RM.map.APOE$APOE.predict<-imputed.APOE.b3.all.noNA[match(RM.map.APOE[,1],names(imputed.APOE.b3.all.noNA))]

RM.map.APOE$APOE.4cohort[which(RM.map.APOE$APOE.4cohort==22)]<-"e2/e2"
RM.map.APOE$APOE.4cohort[which(RM.map.APOE$APOE.4cohort==23)]<-"e2/e3"
RM.map.APOE$APOE.4cohort[which(RM.map.APOE$APOE.4cohort==33)]<-"e3/e3"
RM.map.APOE$APOE.4cohort[which(RM.map.APOE$APOE.4cohort==24)]<-"e2/e4"
RM.map.APOE$APOE.4cohort[which(RM.map.APOE$APOE.4cohort==34)]<-"e3/e4"
RM.map.APOE$APOE.4cohort[which(RM.map.APOE$APOE.4cohort==44)]<-"e4/e4"

table(paste0(RM.map.APOE$APOE.4cohort,".",RM.map.APOE$APOE.predict))

RM.map.APOE$APOE.HDS<-numericMeta.reg$APOE[match(RM.map[,1],rownames(numericMeta.reg))]
table(paste0(RM.map.APOE$APOE.4cohort,".",RM.map.APOE$APOE.HDS))


##  Add known mapped APOE genotypes to ground truth vector (for training)
length(which(is.na(gt.APOE)))
# 7469
mapped.APOE.known.sampleID<-unique(c(intersect(names(gt.APOE)[which(is.na(gt.APOE))],UDS.map.APOE[which(!is.na(UDS.map.APOE$APOE.4cohort)),1]),
                                     intersect(names(gt.APOE)[which(is.na(gt.APOE))],BH.map.APOE[which(!is.na(BH.map.APOE$APOE.4cohort)),1]),
                                     intersect(names(gt.APOE)[which(is.na(gt.APOE))],RM.map.APOE[which(!is.na(RM.map.APOE$APOE.4cohort)),1]) ))
length(mapped.APOE.known.sampleID)
# 1754

#mapped.APOE.known<-rep(NA,length(mapped.APOE.known.sampleID))
#names(mapped.APOE.known)<-mapped.APOE.known.sampleID
#mapped.APOE.known[which(UDS.map.APOE[,1] %in% names(mapped.APOE.known))]<-UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE[,1] %in% names(mapped.APOE.known))]
#mapped.APOE.known[which(BH.map.APOE[,1] %in% names(mapped.APOE.known))]<-BH.map.APOE$APOE.4cohort[which(BH.map.APOE[,1] %in% names(mapped.APOE.known))]
#mapped.APOE.known[which(RM.map.APOE[,1] %in% names(mapped.APOE.known))]<-RM.map.APOE$APOE.4cohort[which(RM.map.APOE[,1] %in% names(mapped.APOE.known))]
mapped.APOE.known<-c(UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE[,1] %in% mapped.APOE.known.sampleID)],
                     BH.map.APOE$APOE.4cohort[which(BH.map.APOE[,1] %in% mapped.APOE.known.sampleID)],
                     RM.map.APOE$APOE.4cohort[which(RM.map.APOE[,1] %in% mapped.APOE.known.sampleID)] )

names(mapped.APOE.known)<-c(UDS.map.APOE[which(UDS.map.APOE[,1] %in% mapped.APOE.known.sampleID),1],
                     BH.map.APOE[which(BH.map.APOE[,1] %in% mapped.APOE.known.sampleID),1],
                     RM.map.APOE[which(RM.map.APOE[,1] %in% mapped.APOE.known.sampleID),1] )

length(which(is.na(mapped.APOE.known)))
# 0

length(which(!is.na(gt.APOE[names(mapped.APOE.known)])))
# 0
gt.APOE[names(mapped.APOE.known)]<-mapped.APOE.known

length(which(is.na(gt.APOE)))
#5715


## Align training set non-missing feature data to gt.APOE (ground truth)
training.cleanDat.noNA<-t(na.omit(t(cleanDat.4p13b3[which(!rownames(cleanDat.4p13b3) %in% c("IRF6|O14896")), which(!is.na(gt.APOE))])))
dim(training.cleanDat.noNA)
#  7334 14758 -- with Mapped APOE   # 7335 13004 with IRF6 in (site K does not have APOE genotypes) -- otherwise it would be 13410 samples for training
training.gt.APOE<-gt.APOE[colnames(training.cleanDat.noNA)]

#############################################################


## Perform Stage 1 determination of binary fit functions -- the variables  binaryPredictionMetrics AND rankedProteins  are saved to the global environment as the list of learning predictors "fit_fns" is calculated.

## -------------------------------------------------------------------
##  (1)  Fit / load the six 1-vs-rest ensembles ----------------------
## -------------------------------------------------------------------
genotypes <- c("e2/e4","e2/e2","e4/e4","e2/e3","e3/e4","e3/e3")

# one list entry per genotype; each entry is the *predict()* function
# returned by the function fit_APOE_binary.95Acc(); run time ~ 12-16 hours on 3.2 GHz multicore CPU for the 6 genotypes

fit_fns <- lapply(genotypes, function(gt)
    fit_APOE_binary.95Acc(expr   = t(training.cleanDat.noNA[                                           ,]),     # FIRST TIME - ALL FEATURES
#    fit_APOE_binary.95Acc(expr   = t(training.cleanDat.noNA[rankedProteins.prior[[gt]][1:50,"feature"],]),     # your 15 000 x p matrix
                          APOE_gt= training.gt.APOE,  # factor of length 15 000
                          target = gt,
                          ncores=8,
                          target_ppv=0.95,  # HERE WE USE DEFAULT, 0.95; in stage 2, with only 50 features, we set this to 0.50 (but the floor is 0.80 in the helper)
                          seed   = 42))

#fit_fns<-list(e2e4=predict_APOE.b3.allAndMapped.e24only.95Acc,
#              e2e2=predict_APOE.b3.allAndMapped.e22only.95Acc,
#              e4e4=predict_APOE.b3.allAndMapped.e44only.95Acc,
#              e2e3=predict_APOE.b3.allAndMapped.e23only.95Acc,
#              e3e4=predict_APOE.b3.allAndMapped.e34only.95Acc,
#              e3e3=predict_APOE.b3.allAndMapped.e34only.95Acc)
#
names(fit_fns) <- genotypes

#fit_fns_0.95ppv<-fit_fns


rankedProteins.prior<-rankedProteins
## ^ importance-ranked assays are input for stage 2


##########################################################
## Below code is not needed to proceed to stage II, now that we have the top 50 ranked important features for predicting each genotype
##
## (OK to skip to APOE_prediction-stage_2.R)



# ------------------------------------------------------------------------
# ANNOTATION: Optional scorer utilities below can apply trained binary
# learners to new samples or matrices.
# ------------------------------------------------------------------------
## -------------------------------------------------------------------
##  (2)  A helper that asks *all* six models for probabilities -------
## -------------------------------------------------------------------
score_all <- function(new_expr_row) {
  # new_expr_row is a 1 x p data.frame or matrix
  stopifnot(is.matrix(new_expr_row) || is.data.frame(new_expr_row))

  sapply(genotypes, function(gt) {

    ## --- extract the binary-model & its feature set ----
    fit_fn  <- fit_fns[[gt]]
    feats   <- attr(fit_fn, "features")

    ## --- give the model exactly the columns it expects ---
    row_sub <- as.matrix(new_expr_row)[ , feats, drop = FALSE]

    # add zero columns for missing proteins (harmless because data are z-scored)
    miss <- setdiff(feats, colnames(row_sub))
    if (length(miss))
      row_sub <- cbind(row_sub, matrix(0, nrow = 1, ncol = length(miss),
                                       dimnames = list(NULL, miss)))

    row_sub <- row_sub[ , feats, drop = FALSE]

    ## --- run the model and read its probability ----------
    res     <- fit_fn(row_sub)          # factor with attributes
    attr(res, "prob")                   # numeric scalar
  })
}

## ------------------------------------------------------------------
## ---------------------------------------------------------------
## (1)  make sure we have the sample x feature matrix
## ---------------------------------------------------------------
expr_all <- t(training.cleanDat.noNA)          # 14 758 x  7334

## ---------------------------------------------------------------
## (2)  collect probabilities from the six binary models
## ---------------------------------------------------------------
prob_train <- sapply(genotypes, function(gt) {
  ## each wrapper returns a factor; its "prob" attribute is the numeric vector
  attr(fit_fns[[gt]](expr_all), "prob")        # length == nrow(expr_all)
})

dim(prob_train)        # 14 758 x 6
colnames(prob_train) <- genotypes
rownames(prob_train) <- rownames(expr_all)     # sample names

## ---------------------------------------------------------------
## (3)  reconstruct the 6-way prediction rule
## ---------------------------------------------------------------
pick_one <- function(p_row) {
  names(p_row) <- genotypes                    # safety belt

  ## any genotype that passes its own PPV-/threshold?
  above <- mapply(function(gt, p)
                    p >= attr(fit_fns[[gt]], "thr"),
                  gt = genotypes, p = p_row)

  if (any(above)) {
    margin <- p_row[above] -
              vapply(genotypes[above],
                     function(gt) attr(fit_fns[[gt]], "thr"),
                     numeric(1))
    names(which.max(margin))                   # largest safety margin wins
  } else {
    names(which.max(p_row))                    # otherwise max-probability
  }
}

pred_train <- factor(apply(prob_train, 1, pick_one), levels = genotypes)

## ---------------------------------------------------------------
## (4)  quick sanity check
## ---------------------------------------------------------------
mean(pred_train == training.gt.APOE)            # overall training accuracy
table(pred_train, training.gt.APOE)             # confusion matrix

summary(prob_train)
hist(prob_train[ , "e4/e4"])


## -------------------------------------------------------------------
##  (3)  The 6-class prediction wrapper ------------------------------
## -------------------------------------------------------------------
predict_APOE6 <- function(new_expr,
                          ncores = parallel::detectCores()) {
  stopifnot(is.matrix(new_expr) || is.data.frame(new_expr))

  if (is.data.frame(new_expr))
      new_expr <- as.matrix(new_expr)

  features <- colnames(new_expr)        # reuse later for each row

  # preserve a possibly pre-existing future plan
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  ## ------------------------------------------------------------------------
  ##  (a)  set up a parallel pool for the *row* computations
  ## ------------------------------------------------------------------------
  future::plan(multisession, workers = ncores)

  ## ------------------------------------------------------------------------
  ##  (b)  build the N x 6 probability matrix in parallel
  ## ------------------------------------------------------------------------
  prob_mat <- t(future.apply::future_apply(
                  new_expr,                  # the full matrix
                  MARGIN= 1L,                # row by row
                  future.seed = TRUE,        # keep RNG reproducible
                  function(x) {
                    x.df<-as.data.frame(x)
                    rownames(x.df)<-features # keep column order / names
                    score_all(t(x.df))       # returns the 6-probability vector
                  }))

  ## ------------------------------------------------------------------------
  ##  (c)  Post-processing
  ## ------------------------------------------------------------------------
  res <- character(nrow(prob_mat))                      # output vector

  for (k in seq_len(nrow(prob_mat))) {
    p_row <- prob_mat[k, ]
    names(p_row)<-genotypes
    # STEP 1  - genotypes whose p >= private thr
    above_thr <- mapply(function(gt, p) p >= attr(fit_fns[[gt]], "thr"),
                        gt = genotypes, p = p_row)

    if (any(above_thr)) {
      # choose the one with the largest (p - thr) margin
      margin <- p_row[above_thr] -
                vapply(genotypes[above_thr],
                       function(gt) attr(fit_fns[[gt]], "thr"), 0)
      res[k] <- names(which.max(margin))
    } else {
      # STEP 2  - fallback: simply pick the max probability
      res[k] <- names(which.max(p_row))
    }
  }

  factor(res, levels = genotypes)
}


#attr(fit_fun, "thr")   # numeric scalar
#attr(fit_fun, "prob")  # internal helper set at first use

train.APOEpred_0.95ppv <- predict_APOE6(t(training.cleanDat.noNA))      # 15 000 predictions

table(train.APOEpred_0.95ppv, training.gt.APOE)  # binary predictors targeting 0.95 ppv (originally)
#                       training.gt.APOE
#train.APOEpred_0.95ppv e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#                 e2/e4     4    16   276     9    14     3
#                 e2/e2     6     0     0     0     0     0
#                 e4/e4     0     2     2    27    51   745
#                 e2/e3    39   931     8    37    36     1
#                 e3/e4     0   168    64   234  4354    77
#                 e3/e3     1   270    10  7045   293    35
mean(train.APOEpred_0.95ppv == training.gt.APOE)
#0.9050684
