## =============================================================================
## 7c. ExtractSigAssays_2directions_DOWN - sliding-window
##     significant-assay extraction (DOWN direction only) and
##     ontology-level Z-score trajectory/heatmap summarization
## =============================================================================
##
## PURPOSE
##   Companion to 7a (which covers ALL = up+down combined): repeats the same
##   sliding-window significant-assay extraction and ontology hypergeometric
##   enrichment pipeline, but restricts each EYO window's test-set assay list
##   to those whose Model 1 posterior trajectory is down-regulated in
##   APOE e4/e4 vs e3/e3 (color-coded "cornflowerblue" in the up/down direction
##   matrix from 5a2/5e2). See 7a's header for the full shared pipeline
##   description; only DOWN-specific differences are called out below.
##
## KEY DIFFERENCES FROM 7a (ALL)
##   - The per-window gene list fed to the enrichment engine is
##     collect_sigDN_by_window(), i.e. down-only assays.
##   - The enrichment-heatmap significance threshold is z = 1.3 (vs 1.645 in 7a).
##   - The combined assay-count + WGCNA-module-overlay PDF
##     (3.plot.hitCounts_plusModulesHitOverlaid-...pdf) is identical across
##     7a/7b/7c (it always plots all three direction series together); it is
##     written once by 7a only, and the redundant duplicate write that used
##     to exist here (wrapped in <SKIP>...<END SKIP>) has been removed.
##   - After enrichment, rather than independently re-curating a term list,
##     this script LOADS the 141-term/18-category ontology selection and
##     row order already curated from the ALL (7a) results
##     (ALL.df_all-141ontologies.categories18.RData), then recomputes and
##     plots DOWN-only Z-scores for those SAME 141 terms, so that the
##     ALL/UP/DOWN publication heatmaps are directly, term-for-term
##     comparable.
##
## REQUIRED INPUTS
##   - ../simple.3177/scatter/_99_par_diff_all_peptide_p_value.rds  (Model 1
##     posterior p-value matrix; output of 5a2/5e2)
##   - ../simple.3177/scatter/_99_par_diff_all_peptide_up_down_notation.csv
##     (Model 1 up/down direction-color matrix; output of 5a2/5e2)
##   - GMT gene-set database file(s) and go.obo (same as 7a; downloaded
##     automatically if missing and running interactively)
##   - F:/.../SEPA.all_redelivered/ALL.df_all-141ontologies.categories18.RData
##     (curated 141-term/18-category ontology selection; produced by 7a)
##   - F:/.../3.Five_yr_slidingWindow/ALL_heatmap18categoryROWorder.txt and
##     .../ALL141terms_18categories-FinalOrder2.tsv (final row/category order
##     tables for the publication heatmap)
##
## MAJOR OUTPUTS
##   - ./SEPA.all_redelivered/DOWN_assay-Zscore_trajectories.pdf /
##     ..._29_SELECTED_1pp.pdf
##   - ./SEPA.all_redelivered/DOWN_assaysONLY-Z_GO-heatmaps_rowsOrdered_top100.pdf /
##     ..._29_SELECTED_1pp.pdf
##   - .../3.SEPAwindows(67)DOWN(top100)_18categoriesSeparated-Z_GO-heatmaps_rowsOrdered.pdf
##     (final 18-category, 141-term DOWN-only publication heatmap)
##   - ./SEPA.all_redelivered/DOWN-141Heatmap-ordered_ontologies+Zscores.csv
##   - ./SEPA.all_redelivered/DOWN-Genes_Hit(141selectedOntologies_67_5yrWindows.csv
##   - ./SEPA.all_redelivered/saved.image-67_5yearInterval_slidingWindow_GOparallel.DOWNhitsONLY-complete-7289_redelivered.RData
## =============================================================================

#####################################################################################################################################################################################
## Start Ontologies spaghetti and heatmap plots (both directions in same lists)        ##   GSVA list generation and waterfall on MMS-PROTEOMICS-KB VM starts below line 700.
##
## - adding logic to split gene lists to up and down by also reading color annotation for waterfall input
#####################################################################################################################################################################################
setwd("f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/")
# ---- STEP 1. Load the Model 1 posterior p-value matrix and drop non-protein
# "seq.####" control assays. ----
pVals<-readRDS("../simple.3177/scatter/_99_par_diff_all_peptide_p_value.rds")
dim(pVals)
#[1]  143 7347
# last protein assay is 7333.


nonProteinAssays<-colnames(pVals)[which(grepl("^seq\\.",colnames(pVals)))]
length(nonProteinAssays)
# 46


pVals<-pVals[,which(!colnames(pVals) %in% nonProteinAssays)]
dim(pVals)
#  143 7301  # 7347 - 46 = 7301


# ---- STEP 2. Define the sliding-window significant-assay collector
# functions (shared with 7a; only the direction-specific one is used below). ----
# In R, a data frame of p values, pVals, has rows named by 1/2 year intervals from "-46", "-45.5", ... to "25" (+25). There are 7345 columns. Provide R code that collects a vector of all colnames that have at least one p value less than or equal to 0.005 within a 5-year interval starting at "-46" and ending at "-41" into a list element named by the midpoint of the 5 year interval ("-43.5"), then repeat the process incrementing the interval plus 1 year.  The end result of the code or function should be a list with 67 elements named "-43.5" to "22.5" (character strings), each a vector of column names from the pVals data frame.
## pVals: data.frame or matrix of p-values
collect_sig_by_window <- function(pVals,
                                  alpha  = 0.005,
                                  start  = -46,
                                  end    = 25,
                                  window = 5,
                                  shift  = 1) {

  yrs <- as.numeric(rownames(pVals))   # row names like "-46", "-45.5", ..., "25"
  starts <- seq(start, end - window, by = shift)    # -46, -45, ..., 20
  mids   <- starts + window/2                       # -43.5, -42.5, ..., 22.5

  out <- setNames(
    lapply(starts, function(s) {
      idx <- yrs >= s & yrs <= (s + window)
      if (!any(idx)) return(character(0))
      keep <- colSums(pVals[idx, , drop = FALSE] <= alpha, na.rm = TRUE) > 0
      colnames(pVals)[keep]
    }),
    as.character(mids)
  )
  out
}

## ---- run it ---------------------------------------------------------------
sig_list <- collect_sig_by_window(pVals) #[,1:7333]

## (Optional sanity checks)
length(sig_list)        # should be 67
head(names(sig_list), 1)  # "-43.5"
tail(names(sig_list), 1)  # "22.5"


assayCount<-unlist(lapply(sig_list,length))
assayCount
## Current (5xSD model data)
#-43.5 -42.5 -41.5 -40.5 -39.5 -38.5 -37.5 -36.5 -35.5 -34.5 -33.5 -32.5 -31.5 -30.5 -29.5 -28.5 -27.5 -26.5 -25.5 -24.5 -23.5 -22.5 -21.5 -20.5 -19.5 -18.5 -17.5 -16.5 
#  164   164   165   168   168   172   173   173   180   178   181   180   182   177   179   182   181   185   187   191   200   200   207   218   231   236   248   259 
#-15.5 -14.5 -13.5 -12.5 -11.5 -10.5  -9.5  -8.5  -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5   0.5   1.5   2.5   3.5   4.5   5.5   6.5   7.5   8.5   9.5  10.5  11.5 
#  271   285   298   315   332   344   345   347   359   368   372   381   389   396   405   411   424   449   480   507   544   572   588   599   600   587   556   522 
# 12.5  13.5  14.5  15.5  16.5  17.5  18.5  19.5  20.5  21.5  22.5 
#  473   425   364   325   283   260   243   227   209   195   182 


## Previously
#-43.5 -42.5 -41.5 -40.5 -39.5 -38.5 -37.5 -36.5 -35.5 -34.5 -33.5 -32.5 -31.5 -30.5 -29.5 -28.5 -27.5 -26.5 -25.5 -24.5 -23.5 -22.5 -21.5 
#  164   164   165   168   168   172   173   173   180   178   181   180   182   177   179   182   181   185   187   191   200   200   207 
#-20.5 -19.5 -18.5 -17.5 -16.5 -15.5 -14.5 -13.5 -12.5 -11.5 -10.5  -9.5  -8.5  -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5   0.5   1.5 
#  218   231   236   248   259   271   285   298   316   333   345   346   348   360   369   373   382   390   396   405   411   424   449 
#  2.5   3.5   4.5   5.5   6.5   7.5   8.5   9.5  10.5  11.5  12.5  13.5  14.5  15.5  16.5  17.5  18.5  19.5  20.5  21.5  22.5 
#  480   507   544   572   588   599   600   587   557   523   475   427   366   327   285   262   244   228   210   196   183 

assayCount<-unlist(lapply(sig_list,length))

direction_colors<-t(read.csv(file="../simple.3177/scatter/_99_par_diff_all_peptide_up_down_notation.csv",header=TRUE,row.names=1,check.names=FALSE))


direction_colors<-direction_colors[,which(!colnames(direction_colors) %in% nonProteinAssays)]

all(colnames(pVals)==colnames(direction_colors))
#TRUE  -- assays in same order

collect_sigUP_by_window <- function(directionColors,
                                  alpha  = 0.005,
                                  start  = -46,
                                  end    = 25,
                                  window = 5,
                                  shift  = 1) {

  yrs <- as.numeric(rownames(directionColors))   # row names like "-46", "-45.5", ..., "25"
  starts <- seq(start, end - window, by = shift)    # -46, -45, ..., 20
  mids   <- starts + window/2                       # -43.5, -42.5, ..., 22.5

  out <- setNames(
    lapply(starts, function(s) {
      idx <- yrs >= s & yrs <= (s + window)
      if (!any(idx)) return(character(0))
      keep <- colSums(directionColors[idx, , drop = FALSE] == "indianred3", na.rm = TRUE) > 0
      colnames(directionColors)[keep]
    }),
    as.character(mids)
  )
  out
}

collect_sigDN_by_window <- function(directionColors,
                                  alpha  = 0.005,
                                  start  = -46,
                                  end    = 25,
                                  window = 5,
                                  shift  = 1) {

  yrs <- as.numeric(rownames(directionColors))   # row names like "-46", "-45.5", ..., "25"
  starts <- seq(start, end - window, by = shift)    # -46, -45, ..., 20
  mids   <- starts + window/2                       # -43.5, -42.5, ..., 22.5

  out <- setNames(
    lapply(starts, function(s) {
      idx <- yrs >= s & yrs <= (s + window)
      if (!any(idx)) return(character(0))
      keep <- colSums(directionColors[idx, , drop = FALSE] == "cornflowerblue", na.rm = TRUE) > 0
      colnames(directionColors)[keep]
    }),
    as.character(mids)
  )
  out
}
## ---- run it again - UP, DN -----------------------------------------------
sig_listUP <- collect_sigUP_by_window(direction_colors) #[,1:7333]

## (Optional sanity checks)
length(sig_listUP)        # should be 67
head(names(sig_listUP), 1)  # "-43.5"
tail(names(sig_listUP), 1)  # "22.5"


assayCountUP<-unlist(lapply(sig_listUP,length))
assayCountUP
## Current
#-43.5 -42.5 -41.5 -40.5 -39.5 -38.5 -37.5 -36.5 -35.5 -34.5 -33.5 -32.5 -31.5 -30.5 -29.5 -28.5 -27.5 -26.5 -25.5 -24.5 -23.5 -22.5 -21.5 -20.5 -19.5 -18.5 -17.5 -16.5 
#   71    70    70    70    67    67    68    67    69    68    68    66    68    69    69    69    67    68    69    73    76    76    78    83    89    93   102   112 
#-15.5 -14.5 -13.5 -12.5 -11.5 -10.5  -9.5  -8.5  -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5   0.5   1.5   2.5   3.5   4.5   5.5   6.5   7.5   8.5   9.5  10.5  11.5 
#  117   122   130   140   149   157   159   162   171   177   182   188   199   202   211   214   221   239   260   279   304   324   330   344   349   341   322   299 
# 12.5  13.5  14.5  15.5  16.5  17.5  18.5  19.5  20.5  21.5  22.5 
#  269   245   204   183   154   142   138   129   119   112   103

## Previous
#-43.5 -42.5 -41.5 -40.5 -39.5 -38.5 -37.5 -36.5 -35.5 -34.5 -33.5 -32.5 -31.5 -30.5 -29.5 -28.5 -27.5 -26.5 -25.5 -24.5 -23.5 -22.5 -21.5 -20.5 -19.5 -18.5 
#   71    70    70    70    67    67    68    67    69    68    68    66    68    69    69    69    67    68    69    73    76    76    78    83    89    93 
#-17.5 -16.5 -15.5 -14.5 -13.5 -12.5 -11.5 -10.5  -9.5  -8.5  -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5   0.5   1.5   2.5   3.5   4.5   5.5   6.5   7.5 
#  102   112   117   122   130   141   150   158   160   163   172   178   183   189   200   202   211   214   221   239   260   279   304   324   330   344 
#  8.5   9.5  10.5  11.5  12.5  13.5  14.5  15.5  16.5  17.5  18.5  19.5  20.5  21.5  22.5 
#  349   341   323   300   271   247   206   185   156   144   139   130   120   113   104

assayCountUP<-unlist(lapply(sig_listUP,length))


sig_listDN <- collect_sigDN_by_window(direction_colors) #[,1:7333]

## (Optional sanity checks)
length(sig_listDN)        # should be 67
head(names(sig_listDN), 1)  # "-43.5"
tail(names(sig_listDN), 1)  # "22.5"


assayCountDN<-unlist(lapply(sig_listDN,length))
assayCountDN
## Current
#-43.5 -42.5 -41.5 -40.5 -39.5 -38.5 -37.5 -36.5 -35.5 -34.5 -33.5 -32.5 -31.5 -30.5 -29.5 -28.5 -27.5 -26.5 -25.5 -24.5 -23.5 -22.5 -21.5 -20.5 -19.5 -18.5 -17.5 -16.5 
#   98   100   100   103   106   108   109   108   112   115   114   116   117   111   113   115   116   122   122   121   126   129   133   136   144   146   152   154 
#-15.5 -14.5 -13.5 -12.5 -11.5 -10.5  -9.5  -8.5  -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5   0.5   1.5   2.5   3.5   4.5   5.5   6.5   7.5   8.5   9.5  10.5  11.5 
#  157   164   169   177   184   190   190   188   189   194   194   194   192   196   195   199   204   212   227   231   243   253   260   260   254   247   235   226 
# 12.5  13.5  14.5  15.5  16.5  17.5  18.5  19.5  20.5  21.5  22.5 
#  206   180   163   143   131   120   105    99    94    88    84

## Previous
#-43.5 -42.5 -41.5 -40.5 -39.5 -38.5 -37.5 -36.5 -35.5 -34.5 -33.5 -32.5 -31.5 -30.5 -29.5 -28.5 -27.5 -26.5 -25.5 -24.5 -23.5 -22.5 -21.5 -20.5 -19.5 -18.5 
#   98   100   100   103   106   108   109   108   112   115   114   116   117   111   113   115   116   122   122   121   126   129   133   136   144   146 
#-17.5 -16.5 -15.5 -14.5 -13.5 -12.5 -11.5 -10.5  -9.5  -8.5  -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5   0.5   1.5   2.5   3.5   4.5   5.5   6.5   7.5 
#  152   154   157   164   169   177   184   190   190   188   189   194   194   194   192   196   195   199   204   212   227   231   243   253   260   260 
#  8.5   9.5  10.5  11.5  12.5  13.5  14.5  15.5  16.5  17.5  18.5  19.5  20.5  21.5  22.5 
#  254   247   235   226   206   180   163   143   131   120   105    99    94    88    84 


assayCountDN<-unlist(lapply(sig_listDN,length))





# ---- STEP 3. Plot the number of significant assays per EYO window
# (overall = black, UP = red, DOWN = blue); same figure as in 7a. ----
## Plot of Assay Counts reaching significance - All (black); Up (red); and Down (blue)
plot.new()
par(mar=c(4.5,5.5,5,2))
plot(as.numeric(names(assayCount)),assayCount, ylab="Number of dysregulated assays\nin e4/4 homozygotes vs e3/3",xlab="EYO (years)",
     main="Protein Counts Significant within 5-year Intervals relative to EYO\n(5-year Sliding Windows, 1 yr resolution)", type="b", col="black", ylim=c(0,605), lwd=2)
lines(as.numeric(names(assayCount)),assayCountUP, type="b", col="indianred3", lwd=2)
lines(as.numeric(names(assayCount)),assayCountDN, type="b", col="cornflowerblue", lwd=2)


# ---- STEP 4. For each of the 12 standard WGCNA module colors, find the
# first/last EYO window containing a significant assay of that color, and
# overlay those ranges as colored blocks on the assay-count plot. ----
# For each of 12 colors defined by the function WGCNA::labels2colors(1:12), check sig_list for the first and last named list element in which each is found. Return a list modsHit with an element for each color that is found in sig_list, naming each element the color found; that element should contain a vector of two values which are the first and last sig_list list element names for list elements in which they were found, converted from strings to as.numeric().
library(WGCNA)
colors12 <- unique(labels2colors(1:12))  # the 12 standard WGCNA colors
colorBlockHeight=80
modsHit <- list()
nms <- names(sig_list)

for (clr in colors12) {
  hits <- which(vapply(sig_listDN, function(v) clr %in% v, logical(1)))
  if (length(hits) > 0) {
    first_name <- nms[min(hits)]
    last_name  <- nms[max(hits)]
    modsHit[[clr]] <- as.numeric(c(first_name, last_name))  # convert names to numeric
  }
}

# modsHit is a list with one element per color found; each value is c(first, last)
str(modsHit)
#List of 3
# $ blue   : num [1:2] -43.5 -40.5
# $ magenta: num [1:2] -43.5 22.5
# $ purple : num [1:2] -4.5 12.5

# Following generation of modsHit, suggest R code that can add to an existing base R plot a very transparent overlaid borderless rectangle of that color for the full y range and x ranging from the first vector value for that color stored in modsHit to the second vector value stored in modsHit for that color.
usr <- par("usr")               # c(xmin, xmax, ymin, ymax)
ybot <- usr[3]; #ytop <- usr[4]

for (clr in names(modsHit)) {
  xr <- modsHit[[clr]]
  rect(xleft = xr[1], ybottom = ybot,
       xright = xr[2], ytop = ybot+colorBlockHeight,
       col = adjustcolor(clr, alpha.f = 0.8),  # very transparent
       border = NA)                             # no border
  text(x=(xr[1]+xr[2])/2,y=ybot+colorBlockHeight/2,labels=clr,col="white", cex=1.4, font=2)
  ybot=ybot+colorBlockHeight
}



modsHit <- list()
nms <- names(sig_list)

for (clr in colors12) {
  hits <- which(vapply(sig_listUP, function(v) clr %in% v, logical(1)))
  if (length(hits) > 0) {
    first_name <- nms[min(hits)]
    last_name  <- nms[max(hits)]
    modsHit[[clr]] <- as.numeric(c(first_name, last_name))  # convert names to numeric
  }
}

# modsHit is a list with one element per color found; each value is c(first, last)
str(modsHit)
#List of 3
# $ blue   : num [1:2] -43.5 -40.5
# $ magenta: num [1:2] -43.5 22.5
# $ purple : num [1:2] -4.5 12.5

# Following generation of modsHit, suggest R code that can add to an existing base R plot a very transparent overlaid borderless rectangle of that color for the full y range and x ranging from the first vector value for that color stored in modsHit to the second vector value stored in modsHit for that color.
usr <- par("usr")               # c(xmin, xmax, ymin, ymax)
ytop <- usr[4] #; ybot <- usr[3]

for (clr in names(modsHit)) {
  xr <- modsHit[[clr]]
  rect(xleft = xr[1], ybottom = ytop-colorBlockHeight,
       xright = xr[2], ytop = ytop,
       col = adjustcolor(clr, alpha.f = 0.8),  # very transparent
       border = NA)                             # no border
  text(x=(xr[1]+xr[2])/2,y=ytop-colorBlockHeight/2,labels=clr,col="white", cex=1.4, font=2)
  ytop=ytop-colorBlockHeight
}

legend("topright",c("Down Assays", "Up Assays", "All sig Assays"), col=c("cornflowerblue","indianred3","black"),pch=21,bg="white", lwd=2, cex=1.45)


# NOTE: the combined assay-count + module-overlay PDF
# (3.plot.hitCounts_plusModulesHitOverlaid-07-27_5xSDoutliersNotInModels.pdf)
# is identical across the ALL/UP/DOWN variants of this script (it always
# plots all three direction series) and is already written once by 7a; the
# redundant duplicate write here (previously wrapped in <SKIP>...<END SKIP>)
# has been removed.


# ---- STEP 5. Configure and run the embedded GOparallel-style ontology
# hypergeometric enrichment engine (shared with 7a), once per EYO window,
# using each window's DOWN-only significant-assay list as the test set. ----
##############################
## Run GOparallel  (modulesData is the list usually filled by reading input .csv with columns of lists)

inputFile <- "dummyFilename.csv"
filePath <- "f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/"   #gsub("//","/",outputfigs)
            #Folder that (may) contain the input file specified above, and which will contain the outFilename project Folder.
outFilename <- "SEPA.down"  #SUBFOLDER WITH THIS NAME WILL BE CREATED, and .PDF + .csv file using the same name will be created within this folder.
outputGOeliteInputs=FALSE  #If TRUE, GO Elite background file and module or list-specific input files will be created in the outFilename subfolder.
maxBarsPerOntology=25      #Ontologies per ontology type, used for generating the PDF report; does not limit tabled output
GMTdatabaseFile="f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/Human_GO_AllPathways_noPFOCR_with_GO_iea_June_01_2025_symbol.gmt"
GO.OBOfile<-"f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/go.obo"
modulesInMemory=FALSE
ANOVAgroups=FALSE
parallelThreads=31
WGCNAinput=FALSE

#modulesData<-sig_list.SomaOnly <- collect_sigUP_by_window(direction_colors[,1:(7333-length(nonProteinAssays))])  #<-just do the UP lists here
modulesData<-sig_list.SomaOnly <- collect_sigDN_by_window(direction_colors[,1:(7333-length(nonProteinAssays))])  #<-just do the UP lists here
#modulesData<-sig_list.SomaOnly <- collect_sig_by_window(pVals[,1:(7333-length(nonProteinAssays))])  #previously: any significant hits (UP+DOWN)
modulesData$bkgr<-colnames(direction_colors)[1:(7333-length(nonProteinAssays))]

                          #changed to 3 from 5
	minHitsPerOntology=3  # ontologies hit by fewer than this number of genes in an input list will not be plotted in Z score barplots.

	if(!exists("filePath")) { cat(paste0("- filePath not set. Using current working directory: ",getwd(),"\n")); filePath=getwd(); }
	## Clean out spaces and escaped backslashes from folder paths (folder names with spaces should not be used on non-windows systems with this script)
	#filePath=paste0(paste( sapply(do.call(c,strsplit(filePath,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(x,paste0(substr(gsub(" ","",x),1,6),"~1"),x) } else { x } } ),collapse="/"),"/")
	filePath=paste0(paste( sapply(do.call(c,strsplit(filePath,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(" ","\ ",x) } else { x } } ),collapse="/"),"/")
	if(!dir.exists(filePath)) { cat(paste0("- filePath set to ",filePath," ...this path was not found. Using current working directory: ",getwd(),"\n")); filePath=getwd(); }
	
	if(!exists("GMTdatabaseFile")) { cat(paste0("- GMTdatabaseFile variable not specified. Current BaderLab .GMT database file will be downloaded to ",filePath,"\n")); GMTdatabaseFile=paste0(filePath,"nonexistent.file"); }
	GMTdatabaseFile=paste0(paste( sapply(do.call(c,strsplit(GMTdatabaseFile,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(" ","\ ",x) } else { x } } ),collapse="/"),"")
	if(!exists("GO.OBOfile")) { cat(paste0("- go.obo file not specified; if needed, current working directory will be checked and if not present, will be downloaded...\n")); GO.OBOfile=paste0(getwd(),"/go.obo"); }
	GO.OBOfile=paste0(paste( sapply(do.call(c,strsplit(GO.OBOfile,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(" ","\ ",x) } else { x } } ),collapse="/"),"")
	
	#pythonPath=paste0(paste( sapply(do.call(c,strsplit(pythonPath,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(x,paste0(substr(gsub(" ","",x),1,6),"~1"),x) } else { x } } ),collapse="/"),"/")
	#GOeliteFolder=paste0(paste( sapply(do.call(c,strsplit(GOeliteFolder,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(x,paste0(substr(gsub(" ","",x),1,6),"~1"),x) } else { x } } ),collapse="/"),"/")
	
	
	## The files we create here are input files for GO-Elite, text files with the gene list as the 1st column, a symbol identified (gene symbol, uniprot etc) as the 2nd column
	## Different accepted inputs are given in the tutorial
	## Commonly used symbols - Gene Symbol - Sy (example of input file below)
	### GeneSymbol		SystemCode (Symbol format)
	###	  GFAP		Sy
	###	  APOE		Sy
	## All input files are placed in one folder
	
	## The background file is prepared similarly and is placed in a separate folder
	## The initial part of the code prepares files for GO-Elite. This can be skipped if the files are being made manually as described above.
	## The second part of the code runs GO-ELite either from R (using the system command) or can be run using the terminal (in mac)
	## The second part requires GO-Elite to be installed and path to the GO-Elite installation site indicated following python
	## The 3rd part of the code plots the results from the GO-Elite results folder. When using the GUI the 1st 2 parts can be skipped and only the 3rd part can be used for plotting
	
	##-------------------------------##
	## Preparing files for GO-Elite ##
	## Takes in the module assignment file as input with 1st column having gene names, 2nd column having color assignments followed by kME values
	

	if (!exists("filePath")) { cat(paste0(" - filePath variable not specified. Input/Output will take place in the current working directory: ",getwd(),"/ ...\n")); filePath==paste0(getwd(),"/"); }
	if (!exists("outFilename")) { cat(paste0("- outFilename variable not specified. Output files will be saved to: ",getwd(),"/GOparallel/ ...\n")); outFilename="GOparallel"; }
	if (!dir.exists(file.path(filePath, outFilename))) dir.create(file.path(filePath, outFilename))
	
	if(!file.exists(GMTdatabaseFile)) {
		if (interactive()) {
			suppressPackageStartupMessages(require(rvest,quietly=TRUE))
			species.links <- html_attr(html_nodes(read_html("http://download.baderlab.org/EM_Genesets/current_release/"), xpath="//a"), "href")
			species.links <- species.links[grepl("^[A-z].*\\/$",species.links)]
			cat("- GMT File not found: ", GMTdatabaseFile,"\n\n")
			print(data.frame(Species=species.links))
			input.idx <- readline(paste0("[INTERACTIVE]\nChoose one of the above species from http://download.baderlab.org/EM_Genesets/current_release/ [1-",length(species.links),"]: "))
			input.idx <- as.integer(input.idx)
		
			find.symbol.in.links <- html_attr(html_nodes(read_html(paste0("http://download.baderlab.org/EM_Genesets/current_release/",species.links[input.idx])), xpath="//a"), "href")
			find.symbol.in.links <- find.symbol.in.links[grepl("[Ss][Yy][Mm][Bb][Oo][Ll]\\/",find.symbol.in.links)]
			gmt.links <- html_attr(html_nodes(read_html(paste0("http://download.baderlab.org/EM_Genesets/current_release/",species.links[as.integer(input.idx)],find.symbol.in.links)), xpath="//a"), "href")
			file.candidates1<-which(grepl("*\\_GO\\_AllPathways\\_.*\\.[Gg][Mm][Tt]",gmt.links))    # main regEx filter for file to download
			file.candidates2<-which(grepl("*\\_noPFOCR.*\\.[Gg][Mm][Tt]",gmt.links))                # files after March 2024 are a subset, excluding PMID-linked lists
			file.candidates3<-which(grepl("*\\_with\\_GO\\_iea\\_.*\\.[Gg][Mm][Tt]",gmt.links))     # take files with automated ontologies
			this.file.idx=if(length(file.candidates2)>0) { intersect(file.candidates2, intersect(file.candidates1,file.candidates3)) } else { intersect(file.candidates1,file.candidates3) }
			if(length(this.file.idx)<1) stop(paste0("Web scraping of the Bader Lab Website could not find an expected GMT filename pattern match.\nDownload and specify a GMTdatabaseFile prior to running this function."))

			full.dl.file=gmt.links[this.file.idx[1] ]
		
			GMTtargetPath=gsub("\\/\\/","/", gsub("(.*\\/).*$","\\1",GMTdatabaseFile) )
			gmt.url<-paste0("http://download.baderlab.org/EM_Genesets/current_release/",species.links[input.idx],find.symbol.in.links,full.dl.file)
			if(file.exists(file.path(GMTtargetPath,full.dl.file))) {
				cat(paste0("- Found that the full current GMT file online matches a file name you already have:\n  ",full.dl.file," [skipping download]\n"))
				GMTdatabaseFile=paste0(GMTtargetPath,full.dl.file)
			} else {
				cat("Found full current GMT file online:  ",gmt.url,"\n")
				cat("Download this file to folder:  ",GMTtargetPath,"\n")
				input.dlYN <- readline("[Y/n]?")
				if(input.dlYN == "Y" | input.dlYN == "y" | input.dlYN == "") {
					suppressPackageStartupMessages(require(curl,quietly=TRUE))
					if (!dir.exists(GMTtargetPath)) dir.create(GMTtargetPath)
					curr.dir<-getwd()
					setwd(GMTtargetPath)
					cat("Downloading .gmt file for ",species.links[input.idx],"...\n")
					curl_download(url=gmt.url, destfile=full.dl.file, quiet = TRUE, mode = "w")
					setwd(curr.dir)
					cat("Using new downloaded .gmt file: ", paste0(GMTtargetPath,full.dl.file),"\n")
					GMTdatabaseFile=paste0(GMTtargetPath,full.dl.file)
				}
			}
		} else { stop(paste0("This is not an interactive session and required GMT file not found.\n",GMTdatabaseFile," must be downloaded interactively or prior to running this function.")) }
	}

	if(!exists("removeRedundantGOterms")) { cat("- removeRedundantGOterms not specified TRUE/FALSE. Removing them as the default, using go.obo and ontologyIndex package.\n"); removeRedundantGOterms=TRUE; }
	if(removeRedundantGOterms) {
		if (!file.exists(GO.OBOfile)) {
			suppressPackageStartupMessages(require(curl,quietly=TRUE))
			OBOtargetPath=gsub("(.*\\/).*$","\\1",GO.OBOfile)
			if (!dir.exists(OBOtargetPath)) dir.create(OBOtargetPath)
			curr.dir<-getwd()
			setwd(OBOtargetPath)
			cat(paste0("- Downloading go.obo file for main GO term redundancy cleanup...\n...to location:  ",OBOtargetPath,"go.obo\n"))
			curl_download(url="http://current.geneontology.org/ontology/go.obo", destfile="go.obo", quiet = TRUE, mode = "w")
			setwd(curr.dir)
			cat("GO.OBOfile set to downloaded file: ", paste0(OBOtargetPath,"go.obo"),"\n")
			GO.OBOfile=paste0(OBOtargetPath,"go.obo")
		}
	}
	

###################################################
	    nModules <- length(names(modulesData))
	    semicolonsFound=FALSE
	    for (a in 1:nModules) {
	      modulesData[[a]] <- unique(modulesData[[a]][modulesData[[a]] != ""])
	      modulesData[[a]] <- modulesData[[a]][!is.na(modulesData[[a]])]
	      modulesData[[a]] <- suppressWarnings(do.call("rbind",strsplit(as.character(modulesData[[a]]), "[|]"))[,1])
	      if(length(which(grepl(";",modulesData[[a]])))>0) {
	        modulesData[[a]]<-suppressWarnings(do.call("rbind",strsplit(as.character(modulesData[[a]]), "[;]"))[,1])
	      }
	    }
	    if(semicolonsFound) cat("- *Found some gene symbols have semicolons! Splitting these and keeping only symbol *before* semicolon.\n")

	    ## Creating background file for GO Elite analysis
	    background <- modulesData[order(sapply(modulesData,length),decreasing=TRUE)][[1]]
	    background <- unique(background)
	    background <- cbind(background,rep("Sy",length=length(background)))
	    colnames(background) <- c("GeneSymbol","SystemCode")
	    if(outputGOeliteInputs) dir.create(file.path(paste0(filePath,outFilename),"background"))
	    if(outputGOeliteInputs) write.table(background,paste0(filePath,outFilename,"/background/background.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
	    
	    # Separate Symbol Lists into independent module txt files for analysis by GO-Elite (not performed by this script)  (CREATE INPUT FILES)
	    modulesData[[ names(modulesData[order(sapply(modulesData,length),decreasing=TRUE)])[1] ]] <- NULL
	    nModules = nModules -1 #no background
	    listNames <- uniquemodcolors <- names(modulesData)
	    for (i in listNames) {
	      listName <- i
	      listInfo <- cbind(modulesData[[listName]],rep("Sy",length=length(modulesData[[listName]])))
	      colnames(listInfo) <- c("GeneSymbol","SystemCode")
	      if(outputGOeliteInputs) write.table(unique(listInfo),file=paste(filePath,outFilename,"/",listName,".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
	    }

#####################################################
	##2. GSA FET (parallelized within R, must have parallelThreads>1 to work currently)
	####----------------------- piano package and dependencies required ------------------------------------#####
	
	suppressPackageStartupMessages(require(piano,quietly=TRUE))
	
	## Adapted version of piano::runGSAhyper() function with depletion p value also calculated (for signed Z score if we will use it to cocluster by module, e.g.)
	runGSAhyper.twoSided <- function(genes, pvalues, pcutoff, universe, gsc, gsSizeLim = c(1,Inf), adjMethod = "fdr") {
	    if (length(gsSizeLim) != 2) 
	        stop("argument gsSizeLim should be a vector of length 2")
	    if (missing(genes)) {
	        stop("argument genes is required")
	    } else {
	        genes <- as.vector(as.matrix(genes))
	        if (!is(genes, "character")) 
	            stop("argument genes should be a character vector")
	        if (length(unique(genes)) != length(genes)) 
	            stop("argument genes should contain no duplicated entries")
	    }
	    if (missing(pvalues)) {
	        pvalues <- rep(0, length(genes))
	    } else {
	        pvalues <- as.vector(as.matrix(pvalues))
	        if (!is(pvalues, "numeric")) 
	            stop("argument pvalues should be a numeric vector")
	        if (length(pvalues) != length(genes)) 
	            stop("argument pvalues should be the same length as argument genes")
	        if (max(pvalues) > 1 | min(pvalues) < 0) 
	            stop("pvalues need to lie between 0 and 1")
	    }
	    if (missing(pcutoff)) {
	        if (all(pvalues %in% c(0, 1))) {
	            pcutoff <- 0
	        } else {
	            pcutoff <- 0.05
	        }
	    } else {
	        if (length(pcutoff) != 1 & !is(pcutoff, "numeric")) 
	            stop("argument pcutoff should be a numeric of length 1")
	        if (max(pcutoff) > 1 | min(pcutoff) < 0) 
	            stop("argument pcutoff needs to lie between 0 and 1")
	    }
	    if (missing(gsc)) {
	        stop("argument gsc needs to be given")
#	    } else {
#	        if (!is(gsc, "GSC")) 
#	            stop("argument gsc should be of class GSC, as returned by the loadGSC function")  # disabled since the list we create is not of GSC class
	    }
	    if (missing(universe)) {
	        if (!all(pvalues == 0)) {
	            universe <- genes
	            message("Using all genes in argument genes as universe.")
	        } else {
	            universe <- unique(unlist(gsc$gsc))
	            message("Using all genes present in argument gsc as universe.")
	        }
	    } else {
	        if (!is(universe, "character")) 
	            stop("argument universe should be a character vector")
	        if (!all(pvalues == 0)) 
	            stop("if universe is given, genes should be only the genes of interest, i.e. pvalues should all be set to 0.")
	    }
	    if (!all(unique(unlist(gsc$gsc)) %in% universe)) 
	        warning("there are genes in gsc that are not in the universe, these will be removed before analysis")
	    if (!all(genes %in% universe)) {
	        warning("not all genes given by argument genes are present in universe, these will be added to universe")
	        universe <- c(universe, genes[!genes %in% universe])
	    }
	    if (length(unique(universe)) != length(universe)) 
	        stop("argument universe should contain no duplicated entries")
	    tmp <- try(adjMethod <- match.arg(adjMethod, c("holm", 
	        "hochberg", "hommel", "bonferroni", 
	        "BH", "BY", "fdr", "none"), several.ok = FALSE), 
	        silent = TRUE)
	    if (is(tmp, "try-error")) {
	        stop("argument adjMethod set to unknown method")
	    }
	    pvalues[pvalues == 0] <- -1e-10
	    goi <- genes[pvalues < pcutoff]
	    if (length(goi) < 1) {
	        cat("\nrunGSEAhyper: no genes selected due to too strict pcutoff. (no genes of interest made an input list)\n")
	        res<-list()
	        res$resTab <- NA
	        res$gsc <- NA
	        return(res)
	    }
	    bg <- universe[!universe %in% goi]
	    gsc <- gsc$gsc
	    delInd <- vector()
	    for (i in 1:length(gsc)) {
	        gs <- gsc[[i]]
	        gs <- gs[gs %in% universe]
	        if (length(gs) < gsSizeLim[1] | length(gs) > gsSizeLim[2]) 
	            delInd <- c(delInd, i)
	        gsc[[i]] <- gs
	    }
	    gsc <- gsc[!c(1:length(gsc)) %in% delInd]
	    message(paste("Analyzing the overrepresentation of ", 
	        length(goi), " genes of interest in ", length(gsc), 
	        " gene sets, using a background of ", length(bg), 
	        " non-interesting genes.", sep = ""))
	    p <- p.depletion <- rep(NA, length(gsc))
	    names(p) <- names(p.depletion) <- names(gsc)
	    padj <- rep(NA, length(gsc))
	    names(padj) <- names(gsc)
	    contTabList <- list()
	    resTab <- matrix(nrow = length(gsc), ncol = 8)  #added 8th column to hold "Genes.Hit"
	    colnames(resTab) <- c("Pvalue.Enrichment", "Adjusted.Enr.Pvalue", "Pvalue.Depletion",
	        "Significant (in gene set)", "Non-significant (in gene set)", 
	        "Significant (not in gene set)", "Non-significant (not in gene set)", "Genes.Hit")
	    rownames(resTab) <- names(gsc)
	    for (i in 1:length(gsc)) {
	        gs <- gsc[[i]]
	        nogs <- universe[!universe %in% gs]
	        ctab <- rbind(c(sum(goi %in% gs), sum(goi %in% nogs)), 
	            c(sum(bg %in% gs), sum(bg %in% nogs)))
	        p[i] <- fisher.test(ctab, alternative = "greater")$p.value
	        p.depletion[i] <- fisher.test(ctab, alternative = "less")$p.value
	        rownames(ctab) <- c("Significant", "Non-significant")
	        colnames(ctab) <- c("Genes in gene set", "Genes not in gene set")
	        contTabList[[i]] <- ctab
	        resTab[i, ] <- c(p[i], NA, p.depletion[i], sum(goi %in% gs), sum(bg %in% 
	            gs), sum(goi %in% nogs), sum(bg %in% nogs), paste0(goi[goi %in% gs],collapse=";"))  #*** added semicolon separated Genes.Hit to 8th column
	    }
	    padj.greater <- p.adjust(p, method = adjMethod)
	    resTab[, 2] <- padj.greater
	    res <- list()
	    res$pvalues.greater <- p
	    res$p.adj.greater <- padj.greater
	    res$pvalues.depletion <- p.depletion
	    res$resTab <- resTab   #*** includes Genes.Hit in 8th column.
	    res$contingencyTable <- contTabList
	    res$gsc <- gsc
	    return(res)
	}
	
	
	## Set up parallel backend.
	suppressPackageStartupMessages(require("doParallel",quietly=TRUE))
	clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
	registerDoParallel(clusterLocal)
	
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
        
        # Time and memory overhead are too great to write and read back in a clean.GMT.  We process the provided .GMT with UTF-8 and inconsistencies every time this script is run.
        #write.table(GMT.df,file="clean.GMT",sep='\t',quote=FALSE, col.names=FALSE, row.names=FALSE)
        #GSCfromGMT<-loadGSC(file="clean.GMT")  # loadGSC(file=GMTdatabaseFile)
	
	## Be sure cluster nodes for parallel processing inherit needed variables from both .GlobalEnv and current function environment (error seen in R 4.2.1 in RStudio on Windows).
	if(!exists("DEXlistsForGO")) DEXlistsForGO<-list()
	parallel::clusterExport(cl=clusterLocal, list("ANOVAgroups","WGCNAinput","background","DEXlistsForGO","GSCfromGMT"), envir=environment())   ## avoid error during foreach below:  Error in { : task 1 failed - "object 'ANOVAgroups' not found"

	
	## Output piano package GSA FET output tables as list assembly
	GSA.FET.outlist<-list()
	
	#  colnames(modulesData)[3:(ncol(modulesData)-1)]
	if (ANOVAgroups) uniquemodcolors=names(DEXlistsForGO)  #otherwise, already set above.
	
	# parallelized to speed up.
	cat("\nRunning FET overlap statistics in parallel for ",length(uniquemodcolors)," symbol lists using up to ", parallelThreads," threads...\n\n")

	#for (this.geneList in uniquemodcolors) {
	GSA.FET.outlist <- foreach(this.geneList=uniquemodcolors) %dopar% {
	#  this.geneList=uniquemodcolors[i]
	  zeroToKeep.idx= if (WGCNAinput) { which( background[,"GeneSymbol"] %in% unique(modulesData[which(modulesData$net.colors==this.geneList),"Unique.ID"]) ) } else {
	                                    if(ANOVAgroups) { which( background[,"GeneSymbol"] %in% DEXlistsForGO[[this.geneList]] ) } else {
	                                            which( background[,"GeneSymbol"] %in% unique(modulesData[[this.geneList]]) ) }}  #Handles file-based input modulesData
	  zeroToKeep=rep(1,nrow(background))
	  zeroToKeep[zeroToKeep.idx]<- 0
	  cat( paste0("\n",this.geneList, "... ") )  #" (n=",length(zeroToKeep.idx)," gene symbols) now processing: ") )
	
	  thislist <- runGSAhyper.twoSided(genes=background[which(zeroToKeep==0),"GeneSymbol"],universe=background[,"GeneSymbol"],gsc=GSCfromGMT,gsSizeLim=c(minHitsPerOntology,Inf),adjMethod="BH",)
	  #above line runs in time, ~30 sec/list, or 50 sec/list for .twoSided
	  #list [[this.color]][["pvalues.greater"]] is enrichment p value vector
	  #list [[this.color]][["padj.greater"]] is FDR vector
	  #list [[this.color]][["resTab"]] is same-ordered (rows) matrix of: p-value, FDR, ... with rownames equal to the ontology name%ontology type%OntologyID
	  #list [[this.color]][["pvalues.less"]] is depletion p value vector (relevant for signed Z score calculation), only in customized function
	
	  return(list(thislist[["resTab"]], thislist[["gsc"]]))
	}
	
	  # re-combine list elements from two outputs, over all uniquemodcolors
	  GSA.FET.resTab.list           = do.call(list,lapply(GSA.FET.outlist,function(x){x[[1]]}))
	  GSA.FET.genesByOntology.list  = do.call(list,lapply(GSA.FET.outlist,function(x){x[[2]]}))
	
	  names(GSA.FET.resTab.list) <- names(GSA.FET.genesByOntology.list) <- uniquemodcolors
	
	
	#Add signed Zscore, pull out ontologyType, ontology (description, title case)
	GSA.FET.outSimple <- lapply(GSA.FET.resTab.list, function(x) { 
	  if(is.na(x[1])) {    #*** occurs when "no genes selected due to too strict pcutoff" in GSEA-FET piano
	    NA
	  } else {
	    #PMC\\d*__F\\d* ontologyTypes for PMC gene sets added December 2023 -- collapse to ontologyType "PMC"
	    #rownames(x)=gsub("^(PMC\\d*__.*)\\t(.*)\\t(.*)",iconv("\\2%PMC%\\1\\t\\3", "latin1","ASCII", ""),rownames(x))
	    ontology=stringr::str_to_title(gsub("\\%WIKIPATHWAYS_\\d*","", gsub("\\%WP_\\d*","", gsub("\\&(.*);","\\1",gsub("<\\sI>","",gsub("<I>","", gsub("(.*)\\%.*\\%.*","\\1",rownames(x))))))))
	    ontologyType=gsub("^WP\\d*","WikiPathways", gsub(".*\\%(.*)\\%.*","\\1",rownames(x)))

	    #force all caps for ontologyType of GObp GOmf GOcc (changed in downloaded GMT files Sept 2022 and/or different in mouse GMT compared to human)
	    ontologyType=gsub("GObp","GOBP",ontologyType)
	    ontologyType=gsub("GOmf","GOMF",ontologyType)
	    ontologyType=gsub("GOcc","GOCC",ontologyType)

	    ZscoreSign=rep(1,nrow(x))
	    ZscoreSign[ as.numeric(x[,"Pvalue.Depletion"]) < as.numeric(x[,"Pvalue.Enrichment"]) ] <- -1
	    Zscore=apply(x, 1, function(p) qnorm(min(as.numeric(p["Pvalue.Enrichment"]), as.numeric(p["Pvalue.Depletion"]))/2, lower.tail=FALSE))
	    out=as.data.frame(x)
	    out$Zscore=Zscore*ZscoreSign
	    out$ontologyType=ontologyType
	    out$ontology=ontology
	    out
	  }
	})



# Given GSA.FET.outSimple, a list of data frames in R, with columns named "Zscore", "ontologyType", and "ontology" and over 13000 rows for each data frame, find the top 25 "ontology" values (character strings) which fall under each of 6 categories stored in ontologyType ("GOBP","GOMF","GOCC","REACTOME","WikiPathways", and "MSIGDB_C2"), with the 25 highest positive Zscores in any of the data frames in the list. If there are ties, keep all tied values.  For each of the 6 ontologyType values, create a data.frame as a list element named by ontologyType, with rownames that are the character string values for each ontology, colnames that are values converted with as.numeric from the names of list elements in GSA.FET.outSimple, and values in the data frame that correspond to the Zscore for each ontology (row) from each named data frame in GSA.FET.outSimple (column).  Call the list variable holding output Zout.list .
library(data.table)

## ---------------- inputs ----------------
## GSA.FET.outSimple: named list of data.frames/data.tables
## each with columns: "Zscore", "ontologyType", "ontology"
## names(GSA.FET.outSimple) are numeric-like (e.g., "-43.5", ...)

## bind all into one long table, tagging each row with its list name as numeric
wins_num <- suppressWarnings(as.numeric(names(GSA.FET.outSimple)))
longDT <- rbindlist(
  Map(function(df, w) {
        dt <- as.data.table(df)
        dt[, win := w]                      # numeric window label (column key)
        dt
      },
      GSA.FET.outSimple, wins_num),
  use.names = TRUE, fill = TRUE
)

## keep only finite, positive Z-scores and required columns
longDT <- longDT[is.finite(Zscore) & Zscore > 0, .(ontologyType, ontology, Zscore, win)]

## collapse duplicates within the same list element (ontology x win)
longDT <- longDT[, .(Zscore = max(Zscore)), by = .(ontologyType, ontology, win)]

## target categories
cats <- c("GOBP","GOMF","GOCC","REACTOME","WikiPathways","MSIGDB_C2")

Zout.list <- setNames(vector("list", length(cats)), cats)

for (cat in cats) {
  dtc <- longDT[ontologyType == cat]
  if (nrow(dtc) == 0L) {
    Zout.list[[cat]] <- data.frame()
    next
  }

  ## --- pick top 100 ontologies by their best (max) Z per ontology across all wins
  maxZ <- dtc[, .(maxZ = max(Zscore)), by = ontology][order(-maxZ)]
  if (nrow(maxZ) == 0L) {
    Zout.list[[cat]] <- data.frame()
    next
  }
  n_top  <- min(100L, nrow(maxZ))
  thresh <- maxZ$maxZ[n_top]                     # include ties at the 25th rank
  top_on <- maxZ[maxZ >= thresh, ontology]

  ## --- wide matrix: rows = ontology, cols = numeric list names, values = Zscore
  wide <- dcast(dtc[ontology %in% top_on],
                ontology ~ win, value.var = "Zscore")

  rn <- wide$ontology
  wide[, ontology := NULL]

  ## ensure columns are ordered by numeric win and named as character(as.numeric())
  col_nums <- as.numeric(names(wide))
  ord <- order(col_nums)
  wide <- wide[, ..ord]
  colnames(wide) <- as.character(col_nums[ord])

  ## finalize as data.frame with rownames
  df <- as.data.frame(wide, check.names = FALSE)
  rownames(df) <- rn

  Zout.list[[cat]] <- df
}

## Zout.list is the requested output:
## - one element per ontologyType (name preserved)
## - rows = ontology strings (top 25 + ties)
## - columns = as.numeric(names(GSA.FET.outSimple)) as strings
## - values = Zscore (NA if ontology absent in that window)

# Intermediate checkpoint: save the workspace right after the (slow,
# parallelized) enrichment step completes, before plotting begins.
save.image("./SEPA.all_redelivered/saved.image-67_5yearInterval_slidingWindow_GOparallel.DOWNhitsONLY-complete-7289_redelivered.RData")


# ---- STEP 7. Plot smoothed Z-score trajectory line plots and row-clustered
# enrichment heatmaps for the top-ranked DOWN-only terms in each
# ontology category (threshold z = 1.3, vs 1.645 in 7a). ----
# With a data frame from a Zout.list element as input, use ggplot2 to plot smoothed lines (curves) for Zscore trajectories across the numbered columns, representing values along the x-axis, with the y-axis representing Zscore for each ontology in one of the categories. Use a minimal theme, but add text to label each Zscore trendline with the rowname of the Zout.list data.frame from which the data for the Zscores is derived. Use different line types and colors to distinguish the 25+ lines drawn. Save the plot to a page of an open PDF with dimensions 8" high by 14" wide, and repeat for each of the 6 Zout.list elements.
## ---- packages ---------------------------------------------------------------
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(ggrepel)

## ---- helper: plot one Zout df ----------------------------------------------
plot_Zout_df <- function(df, title = NULL, span = 0.5) {
  # Ensure plain data.frame
  df <- as.data.frame(df, check.names = FALSE)

  long <- df |>
    rownames_to_column("ontology") |>
    pivot_longer(-ontology, names_to = "win", values_to = "Z") |>
    mutate(win = suppressWarnings(as.numeric(win))) |>
    filter(is.finite(win), !is.na(Z))

  # LOESS needs >= 2 points per group
  long <- long |>
    group_by(ontology) |>
    filter(n() >= 2) |>
    ungroup()

  # Labels at the right-most available x for each ontology
  lab_dat <- long |>
    group_by(ontology) |>
    slice_max(order_by = win, n = 1, with_ties = FALSE) |>
    ungroup()

  # Build enough linetypes for many curves (will recycle if > 6)
  lt_vals <- rep(c("solid","dashed","dotted","dotdash","longdash","twodash"),
                 length.out = dplyr::n_distinct(long$ontology))

  p <- ggplot(long, aes(win, Z, color = ontology, linetype = ontology, group = ontology)) +
    geom_smooth(se = FALSE, method = "loess", span = span, na.rm = TRUE, linewidth = 0.7) +
    geom_text_repel(data = lab_dat, aes(label = ontology),
                    size = 3, hjust = 0, direction = "y",
                    nudge_x = diff(range(long$win, na.rm = TRUE)) * 0.02,
                    segment.size = 0.25, show.legend = FALSE, max.overlaps = Inf) +
    scale_linetype_manual(values = lt_vals, guide = "none") +
    scale_x_continuous(name = "Window (numeric)",
                       expand = expansion(mult = c(0.02, 0.08))) +
    ylab("Z-score") +
    ggtitle(title %||% "e4/4 UP Z-score trajectories") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

  p
}

`%||%` <- function(x, y) if (is.null(x)) y else x

## ---- draw all six pages to one PDF ------------------------------------------
# Zout.list: named list of 6 data.frames (GOBP, GOMF, GOCC, REACTOME, WikiPathways, MSIGDB_C2)
pdf("./SEPA.all_redelivered/DOWN_assay-Zscore_trajectories.pdf", width = 14, height = 8, onefile = TRUE)
for (nm in names(Zout.list)) {
  df <- Zout.list[[nm]]
  if (!is.null(df) && nrow(df) > 0 && ncol(df) > 0) {
    print( plot_Zout_df(df, title = paste0("DOWN assay Z-score trajectories: ", nm), span = 0.5) )
  }
}
dev.off()



## With a data frame from a Zout.list element as input, use the complexHeatmap::pheatmap() function to plot Z scores above 1.30 as a color gradient that goes from black below 1.30, to a deep dark red at 1.30, until reaching yellow at the maximum Z score in the data.frame. The progression of color from left to right should correspond to the increments of columns in the Zout.list element data.frame for each row (ontology). Do not outline cells in the heatmap, and all missing values (NA) should be black. Do not cluster columns, but rows (ontologies) may be clustered. Show the text for ongtology labels to the left of each row, and label all columns divisible by five after subtracting 0.5. Save the plot to a page of an open PDF Z_GO_heatmaps.pdf with dimensions 8" high by 17.5" wide, and repeat for each of the 6 Zout.list elements.

library(ComplexHeatmap)

#---- helper: plot one Zout.list data frame as ComplexHeatmap::pheatmap ----------
plot_Zout_heatmap <- function(df, title = NULL, thr = 1.3) {
  stopifnot(is.data.frame(df) || is.matrix(df))
  mat <- as.matrix(df)
  rn <- rownames(mat)

  # Order columns numerically; drop non-numeric column names if any
  cx  <- suppressWarnings(as.numeric(colnames(mat)))
  keep <- !is.na(cx)
  mat  <- mat[, keep, drop = FALSE]
  cx   <- cx[keep]

  ord <- order(cx)
  mat <- mat[, ord, drop = FALSE]
  cx  <- cx[ord]
  colnames(mat) <- as.character(cx)   # keep as character labels

  # Column labels: only those where (x - 0.5) %% 5 == 0
  lab_col <- ifelse(((cx - 0.5) %% 5) == 0, as.character(cx), "")

  # --- ROW SORT: earliest first-hit where Z >= thr (left to right)
  first_hit <- apply(mat, 1, function(x) {
    wh <- which(!is.na(x) & x >= thr)
    if (length(wh)) wh[1] else Inf
  })
  row_max <- apply(mat, 1, function(x) suppressWarnings(max(x, na.rm = TRUE)))
  row_max[!is.finite(row_max)] <- -Inf

  ord_rows <- order(first_hit, -row_max, rn)
  mat <- mat[ord_rows, , drop = FALSE]

  # --- color mapping: <= 1.30 = black; 1.30 = deep dark red; up to max = yellow
  maxZ <- suppressWarnings(max(mat, na.rm = TRUE))
  minZ <- suppressWarnings(min(mat, na.rm = TRUE))
  if (!is.finite(maxZ)) maxZ <- thr
  if (!is.finite(minZ)) minZ <- thr

  # ensure we have a bin for everything below 1.30 mapping to black
  minBreak <- min(minZ, thr)
  if (maxZ <= thr) {
    # nothing above threshold: still draw with two bins so 1.30 maps to dark red
    breaks <- c(minBreak, thr, thr+1e-6)
    cols   <- c("black", "#8B0000")
  } else {
    n_grad <- 100L
    breaks <- c(minBreak, thr, seq(thr, maxZ, length.out = n_grad + 1L)[-1L])
    cols   <- c("black", colorRampPalette(c("#8B0000", "yellow"))(n_grad))
  }

  # Draw with ComplexHeatmap::pheatmap() (no cell borders, NAs black, no col clustering)
  ht <- ComplexHeatmap::pheatmap(
    mat,
    color          = cols,
    breaks         = breaks,
    na_col         = "black",
    row_names_max_width = unit(6.75, "in"),
    cluster_cols   = FALSE,
    cluster_rows   = FALSE,  # NA not tolerated
    clustering_distance_rows = "euclidean",
    clustering_method = "complete",
    show_rownames  = TRUE,
    show_colnames  = TRUE,
    labels_col     = lab_col,
    border_color   = "#FFFFFFFF",
    main           = title
  )

  # In ComplexHeatmap, pheatmap() returns a Heatmap object; draw it explicitly
  ComplexHeatmap::draw(ht)
}

#---- write all six heatmaps (one per Zout.list element) to a single PDF ---------
# Zout.list: named list with 6 elements (e.g., "GOBP","GOMF","GOCC","REACTOME","WikiPathways","MSIGDB_C2")
pdf("./SEPA.all_redelivered/DOWN_assaysONLY-Z_GO-heatmaps_rowsOrdered.pdf", width = 17.5, height = 35, onefile = TRUE)
#grid::grid.newpage()
#par(mar=c(1,1,1,15))
#par(oma=c(1,1,1,30))
for (nm in names(Zout.list)) {
  df <- Zout.list[[nm]]
  if (!is.null(df) && nrow(df) > 0 && ncol(df) > 0) {
    plot_Zout_heatmap(df, title = paste0("DOWN ASSAYS Z-score heatmap: ", nm))
  }
}
dev.off()


## ---- select separated top terms from all six pages to put on one page : first spaghetti plot to PDF ------------------------------------------
# Zout.list: named list of 6 data.frames (GOBP, GOMF, GOCC, REACTOME, WikiPathways, MSIGDB_C2)

## Eric's 11 selected terms for UP
#selectedTerms<-c("MSIGDB_C2.Pid_hif1_tfpathway","MSIGDB_C2.Pid_cmyb_pathway","GOMF.Steroid Binding","GOMF.Lipid Transporter Activity","GOBP.Postsynapse Assembly","GOBP.Synapse Assembly","GOBP.Positive Regulation Of Camp/Pka Signal Transduction","GOMF.Glycosaminoglycan Binding",
#                 "MSIGDB_C2.Biocarta_hes_pathway","REACTOME.Pre-Notch Transcription And Translation","GOCC.Extracellular Matrix")

## Erik's 29 selected terms for ALL
selectedTerms<-c("GOBP.Alkanesulfonate Metabolic Process","GOBP.Taurine Metabolic Process","GOBP.Hormone Metabolic Process","GOBP.Neuron Development","GOBP.Neurogenesis","GOBP.Negative Regulation Of Canonical Wnt Signaling Pathway","GOBP.Negative Regulation Of Coagulation","GOBP.Negative Regulation Of Wound Healing","GOBP.Kidney Vasculature Morphogenesis","GOMF.Deubiquitinase Activity","GOMF.Ubiquitin-Like Protein Peptidase Activity","GOMF.Extracellular Matrix Binding","GOMF.Insulin-Like Growth Factor I Binding","GOMF.Glycosaminoglycan Binding","GOMF.Heparin Binding","GOMF.Steroid Binding","GOMF.Integrin Binding","GOMF.Insulin-Like Growth Factor Ii Binding","GOMF.Cellular Response To Vascular Endothelial Growth Factor Stimulus","GOCC.Low-Density Lipoprotein Particle","GOCC.High-Density Lipoprotein Particle","GOCC.Triglyceride-Rich Plasma Lipoprotein Particle","GOCC.Very-Low-Density Lipoprotein Particle","GOCC.Chylomicron","GOCC.Extracellular Matrix",
                 "REACTOME.Heparan Sulfate Heparin (Hs-Gag) Metabolism","REACTOME.Complement Cascade","REACTOME.Ncam1 Interactions","REACTOME.Notch-Hlh Transcription Pathway")

Zout.selected<-do.call(rbind,Zout.list)
Zout.selected<-Zout.selected[which(rownames(Zout.selected) %in% selectedTerms),]

pdf("./SEPA.all_redelivered/DOWN_assay-Zscore_trajectories_29_SELECTED_1pp.pdf", width = 14, height = 8, onefile = TRUE)
  df <- Zout.selected
  if (!is.null(df) && nrow(df) > 0 && ncol(df) > 0) {
    print( plot_Zout_df(df, title = paste0("DOWN assay Z-score trajectories: ", "Selected Terms"), span = 0.5) )
  }
dev.off()

# Replot same 11 selected terms as heatmap
pdf("./SEPA.all_redelivered/DOWN_assaysONLY-Z_GO-heatmaps_rowsOrdered_29_SELECTED_1pp.pdf", width = 17.5, height = 8, onefile = TRUE)
  if (!is.null(df) && nrow(df) > 0 && ncol(df) > 0) {
    plot_Zout_heatmap(df, title = paste0("DOWN ASSAYS Z-score heatmap: ", "Selected Terms"))
  }
dev.off()



##############################################################################
## Redefine plot_Zout_heatmap2() with the DOWN-direction (blue-gradient)
## color scale, used by the "GO back and curate" step below.



plot_Zout_heatmap2 <- function(mat, title=NULL, thr=1.3, rn = rownames(mat), draw=TRUE, cols=NULL, breaks=NULL, reorderRows=TRUE,...) {
  # Ensure rownames exist
  if (is.null(rn)) rn <- seq_len(nrow(mat))

  stopifnot(is.data.frame(df) || is.matrix(df))
  mat <- as.matrix(df)
  rn <- rownames(mat)

  # Order columns numerically; drop non-numeric column names if any
  cx  <- suppressWarnings(as.numeric(colnames(mat)))
  keep <- !is.na(cx)
  mat  <- mat[, keep, drop = FALSE]
  cx   <- cx[keep]

  ord <- order(cx)
  mat <- mat[, ord, drop = FALSE]
  cx  <- cx[ord]
  colnames(mat) <- as.character(cx)   # keep as character labels

  # Column labels: only those where (x - 0.5) %% 5 == 0
  lab_col <- ifelse(((cx - 0.5) %% 5) == 0, as.character(cx), "")

  ## --- ROW SORT LOGIC ---
  # 1. Identify rows where ALL values are significant (>= thr) and not NA
  all_sig <- apply(mat, 1, function(x) all(!is.na(x) & x >= thr))

  # 2. For each row, find index of first significant value
  first_hit <- apply(mat, 1, function(x) {
    wh <- which(!is.na(x) & x >= thr)
    if (length(wh)) wh[1] else Inf
  })

  # 3. Z value at that first significant position
  first_hit_val <- mapply(function(row, idx) {
    if (is.finite(idx)) mat[row, idx] else -Inf
  }, row = seq_len(nrow(mat)), idx = first_hit)

  # 4. First non-significant (black/grey) after first_hit
  #    Black/grey = NA or < thr
  first_non_sig_after_hit <- mapply(function(row, idx) {
    if (!is.finite(idx) || idx == ncol(mat)) return(Inf)
    wh <- which(is.na(mat[row, (idx+1):ncol(mat)]) |
                  mat[row, (idx+1):ncol(mat)] < thr)
    if (length(wh)) (idx + wh[1]) else Inf
  }, row = seq_len(nrow(mat)), idx = first_hit)

  # 5. Order:
  #    - all_sig rows first
  #    - then by first_hit (earlier columns first)
  #    - then by higher Z at that first hit
  #    - then by later first_non_sig_after_hit (later = higher priority)
  #    - then by row name
  if (reorderRows) {
    ord_rows <- order(!all_sig, first_hit, -first_hit_val,
                      -first_non_sig_after_hit, rn)

    mat <- mat[ord_rows, , drop = FALSE]
  }

  ## --- COLOR MAPPING ---
  if (is.null(cols) | is.null(breaks)) {
    maxZ <- suppressWarnings(max(mat, na.rm = TRUE))
    minZ <- suppressWarnings(min(mat, na.rm = TRUE))
    if (!is.finite(maxZ)) maxZ <- thr
    if (!is.finite(minZ)) minZ <- thr

    minBreak <- min(minZ, thr)
    if (maxZ <= thr) {
      breaks <- c(minBreak, thr, thr + 1e-6)
      cols   <- c("white", "#CDCDCD")  # was 8B in black-dkred-gold
    } else {
      n_grad <- 100L
      breaks <- c(minBreak, 0.000001, thr-0.000001, seq(thr, maxZ, length.out = n_grad + 1L)[-1L])
      cols   <- c("white", "#CDCDCD", "#CDCDCD", colorRampPalette(c("#5555CD", "yellow"))(n_grad))
    }
  }

  ## --- DRAW HEATMAP ---
  ht <- pheatmap(
    mat,
    color          = cols,
    breaks         = breaks,
    na_col         = "#FFFFFFFF",  # "black",
    row_names_max_width = unit(6.75, "in"),
    cluster_cols   = FALSE,
    cluster_rows   = FALSE,
    show_rownames  = TRUE,
    show_colnames  = TRUE,
    labels_col     = lab_col,
    border_color   = "#EDEDED",
    main           = title,
    name           = "Z-score"  #legend scale name
  )

  if(draw) { ComplexHeatmap::draw(ht) } else { ht }
}


##############################################################################
# ---- STEP 8. Reuse the 141-term/18-category ontology selection already
# curated from the ALL (7a) results, and recompute/plot Z-scores for those
# SAME terms using only this script's DOWN-direction enrichment
# results, so the ALL/UP/DOWN publication heatmaps are term-for-term
# comparable. ----
## GO back and curate the same 141 terms from ALL (UP+DOWN) but Z scores are from DOWN hits only.

load("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/SEPA.all_redelivered/ALL.df_all-141ontologies.categories18.RData")
cats2<-cats[1:4] #GOBP, GOMF, GOCC, REACTOME
ALL.ontologyNames<-rownames(df_all)
order141<- read.delim(file="F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/ALL_heatmap18categoryROWorder.txt",sep="\t",header=FALSE)[,1]  # overwrite what is in the above loaded RData


library(data.table)

Zout.list2 <- setNames(vector("list", length(cats2)), cats2)

for (cat in cats2) {
  dtc <- longDT[which(longDT$ontologyType == cat),]
  if (nrow(dtc) == 0L) {
    Zout.list2[[cat]] <- data.frame()
    next
  }

  dtc$ontologyType.ontology<-paste0(dtc$ontologyType,".",dtc$ontology)
  ## --- keep only ontologies in ALL.ontologyNames
  keep_on <- intersect(unique(dtc$ontologyType.ontology), ALL.ontologyNames)

  if (length(keep_on) == 0L) {
    Zout.list2[[cat]] <- data.frame()
    next
  }

  ## --- wide matrix: rows = ontology, cols = numeric list names, values = Zscore
  wide <- dcast(dtc[ontologyType.ontology %in% keep_on],
                ontology ~ win, value.var = "Zscore")

  rn <- wide$ontology
  wide[, ontology := NULL]

  ## ensure columns are ordered by numeric win and named as character(as.numeric())
  col_nums <- as.numeric(names(wide))
  ord <- order(col_nums)
  wide <- wide[, ..ord]
  colnames(wide) <- as.character(col_nums[ord])

  ## finalize as data.frame with rownames
  df <- as.data.frame(wide, check.names = FALSE)
  rownames(df) <- rn

  Zout.list2[[cat]] <- df
}

## Zout.list2 is the requested output:
## - one element per ontologyType (name preserved)
## - rows = ontology strings (selected from ALL (UP+DOWN) 141, regardless of Z score max)
## - columns = as.numeric(names(GSA.FET.outSimple)) as strings
## - values = Zscore (NA if ontology absent in that window)

Zout.selected3<-do.call(rbind,Zout.list2)


## 18 categories splitting 141 of the above 144 terms. Plot each separately but on one page.
#categories18<-read.delim("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/ALL141terms_18categories.tsv",sep="\t",header=TRUE,row.names=NULL)
categs<-c('Metabolism','RNA Processing','GLP-1 and IGFBP','Lipoproteins','Proteostasis','Synaptic and Neuronal','Eye and Retinol','Steroid Metabolism','ECM','Structural','FGF and Fibronectin','Hemostasis','Adaptive Immune','Innate Immune','Wnt signaling','Apoptosis','Translation','Angiogenesis')

library(ComplexHeatmap)
library(circlize)
library(grid)

# 1. Compute global min/max across all categories
all_rows <- unlist(lapply(categs, function(categ) {
  rn <- sub("^([^\\.]+)\\.(.*)$", "\\2 (\\1)", rownames(Zout.selected3))
  idx <- rn %in% categories18[categories18[,2] == categ, 1]
  rownames(Zout.selected3)[idx]
}))

df_all <- Zout.selected3[all_rows, , drop = FALSE]
# All ontology names you want represented
all_needed <- unique(unlist(lapply(categs, function(categ) {
  categories18[categories18[,2] == categ, 1]
})))

# Which ones are missing from df_all?
missing_rows <- setdiff(all_needed, sub("^([^\\.]+)\\.(.*)$", "\\2 (\\1)",rownames(df_all)))

if (length(missing_rows) > 0) {
  # Create a zero-filled data frame with same columns
  zero_mat <- matrix(0, nrow = length(missing_rows), ncol = ncol(df_all),
                     dimnames = list(missing_rows, colnames(df_all)))
  zero_df <- as.data.frame(zero_mat, check.names = FALSE)
  
  # Append to the end
  df_all <- rbind(df_all, zero_df)
}

df_all[is.na(df_all)] <- 0  # change NA values (grey) to 0 color.
global_minZ <- min(df_all, na.rm = TRUE)
global_maxZ <- max(df_all, na.rm = TRUE)

# 2. Define one shared color function
thr=1.3
if (!is.finite(global_maxZ)) global_maxZ <- thr
if (!is.finite(global_minZ)) global_minZ <- thr

minBreak <- min(global_minZ, thr)
if (global_maxZ <= thr) {
  breaks <- c(minBreak, thr, thr + 1e-6)
  cols   <- c("white", "#5555CD")
} else {
  n_grad <- 100L
  breaks <- c(minBreak, 0.000001, thr-0.000001, seq(thr, global_maxZ, length.out = n_grad + 1L)[-1L])
  cols   <- c("#FFFFFF", "#CDCDCD", "#CDCDCD", colorRampPalette(c("#5555CD", "yellow"))(n_grad))
}

## Looped function call - 1 PDF page one ggplot2 object assembled, drawn by the ComplexHeatmap::draw() function
make_title_ht <- function(title, ncols) {
  mat <- matrix(NA, nrow = 1, ncol = ncols)
  Heatmap(
    mat,
    col = NA,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_heatmap_legend = FALSE,
    row_names_side = "left",
    row_labels = title,
    row_names_gp = gpar(fontsize = 16, fontface = "bold"),
    column_names_gp = gpar(fontsize = 0),
    rect_gp = gpar(col = NA, fill = NA)
  )
}



# ---- STEP 9. Draw the final, curated 18-category / 141-term DOWN-only
# publication heatmap using the fixed FinalOrder2 category/term table. ----
## Heatmap rendered in final row order below 141x67


##############################################################################################
# ---- STEP 10. Write the final 141x67 DOWN-only Z-score matrix and
# matching "genes hit" annotation table to CSV for the manuscript's
# supplementary tables. ----
## Write Publication tables of Zscores + Genes Hit (141x67)

#df_all.clean<-df_all
#rownames(df_all.clean) <- sub("^([^\\.]+)\\.(.*)$", "\\2 (\\1)", rownames(df_all.clean))
#write.csv(df_all.clean,"Zscores_UP_redGT1.30_forSuppTable.csv")


cats2<-cats[1:4] #GOBP, GOMF, GOCC, REACTOME


#order141<-readRDS("F:\\OneDrive - Emory\\Legacy\\e4_homozygoteStudy\\DL\\3.Five_yr_slidingWindow\\SEPA.all_redelivered\\order141.RDS")
categories18<-read.delim("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/SEPA.all_redelivered/ALL141terms_18categories-FinalOrder2.tsv",sep="\t",header=TRUE,row.names=NULL)
categs<-c('Metabolism','RNA Processing','GLP-1 and IGFBP','Lipoproteins','Proteostasis','Synaptic and Neuronal','Eye and Retinol','Steroid Metabolism','ECM','Structural','FGF and Fibronectin','Hemostasis','Adaptive Immune','Innate Immune','Wnt signaling','Apoptosis','Translation','Angiogenesis')

categories18 #[match(order141,categories18$NewAnnotation),]
#  in correct order 1-141


exactGO.141<-categories18 #[match(order141,categories18$NewAnnotation),]
rownames(exactGO.141)<-NULL



############## Final Order Heatmap Drawn here ###
ht_list <- NULL
#order141<-vector()


for (categ in categs) {
  df <- df_all #Zout.selected3
  rownames(df) <- sub("^([^\\.]+)\\.(.*)$", "\\2 (\\1)", rownames(df))
    #df <- df[rownames(df) %in% exactGO.141[exactGO.141[,2] == categ, 1], ]
  df<-df[match(exactGO.141[exactGO.141[,2] == categ, 1], rownames(df), nomatch = 0), , drop = FALSE]
      #df[order141[which(order141 %in% rownames(df))],]

  
  if (nrow(df) > 0 && ncol(df) > 0) {
    title_ht <- make_title_ht(categ, ncol(df))
    ht <- plot_Zout_heatmap2(df, draw = FALSE, breaks = breaks, cols = cols, reorderRows=FALSE)
    section <- title_ht %v% ht
    ht_list <- if (is.null(ht_list)) section else ht_list %v% section
  }
}


pdf("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/SEPA.all_redelivered/3.SEPAwindows(67)DOWN(141ALL)_18categoriesSeparated-Z_GO-heatmaps_rowsOrdered.pdf", width = 20, height = 30, onefile = TRUE)
  draw(ht_list, merge_legend = TRUE)
dev.off()
############################





library(data.table)

## ---------------- inputs ----------------
## GSA.FET.outSimple: named list of data.frames/data.tables
## each with columns: "Zscore", "ontologyType", "ontology"
## names(GSA.FET.outSimple) are numeric-like (e.g., "-43.5", ...)

## bind all into one long table, tagging each row with its list name as numeric
wins_num <- suppressWarnings(as.numeric(names(GSA.FET.outSimple)))
longDT.genesHit <- rbindlist(
  Map(function(df, w) {
        dt <- as.data.table(df)
        dt[, win := w]                      # numeric window label (column key)
        dt
      },
      GSA.FET.outSimple, wins_num),
  use.names = TRUE, fill = TRUE
)

## keep only finite, positive Z-scores and required columns
longDT.genesHit <- longDT.genesHit[is.finite(Zscore) & Zscore > 0, .(ontologyType, ontology, Zscore, Genes.Hit, win)]

## collapse duplicates within the same list element (ontology x win)
longDT.genesHit <- longDT.genesHit[, .(Zscore = max(Zscore)), by = .(ontologyType, ontology, Genes.Hit, win)]

longDT.genesHit$GOexact=paste0(longDT.genesHit$ontology," (",longDT.genesHit$ontologyType,")")


dim(longDT.genesHit)
#[1] 271973     6  (current; 7289 assays)



############### Get Z score matrix (141x67)
library(data.table)

Zout.list2 <- setNames(vector("list", length(cats2)), cats2)

for (cat in cats2) {
  dtc <- longDT.genesHit[which(longDT.genesHit$ontologyType == cat),]
  if (nrow(dtc) == 0L) {
    Zout.list2[[cat]] <- data.frame()
    next
  }

  dtc$ontologyType.ontology<-paste0(dtc$ontologyType,".",dtc$ontology)
  ## --- keep only ontologies in ALL.ontologyNames
  keep_on <- intersect(unique(dtc$GOexact), categories18$NewAnnotation)

  if (length(keep_on) == 0L) {
    Zout.list2[[cat]] <- data.frame()
    next
  }

  ## --- wide matrix: rows = ontology, cols = numeric list names, values = Zscore
  wide <- dcast(dtc[GOexact %in% keep_on],
                ontology ~ win, value.var = "Zscore")

  rn <- wide$ontology
  wide[, ontology := NULL]

  ## ensure columns are ordered by numeric win and named as character(as.numeric())
  col_nums <- as.numeric(names(wide))
  ord <- order(col_nums)
  wide <- wide[, ..ord]
  colnames(wide) <- as.character(col_nums[ord])

  ## finalize as data.frame with rownames
  df <- as.data.frame(wide, check.names = FALSE)
  rownames(df) <- rn

  Zout.list2[[cat]] <- df
}

## Zout.list2 is the requested output:
## - one element per ontologyType (name preserved)
## - rows = ontology strings (selected from ALL (UP+DOWN) 141, regardless of Z score max)
## - columns = as.numeric(names(GSA.FET.outSimple)) as strings
## - values = Zscore (NA if ontology absent in that window)

Zout.selected3<-do.call(rbind,Zout.list2)


### 18 categories splitting 141 of the above 144 terms. Plot each separately but on one page.
##categories18<-read.delim("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/ALL141terms_18categories.tsv",sep="\t",header=TRUE,row.names=NULL)
#categs<-c('Metabolism','RNA Processing','GLP-1 and IGFBP','Lipoproteins','Proteostasis','Synaptic and Neuronal','Eye and Retinol','Steroid Metabolism','ECM','Structural','FGF and Fibronectin','Hemostasis','Adaptive Immune','Innate Immune','Wnt signaling','Apoptosis','Translation','Angiogenesis')

library(ComplexHeatmap)
library(circlize)
library(grid)

# 1. Compute global min/max across all categories
all_rows <- unlist(lapply(categs, function(categ) {
  rn <- sub("^([^\\.]+)\\.(.*)$", "\\2 (\\1)", rownames(Zout.selected3))
  idx <- rn %in% categories18[categories18[,2] == categ, 1]
  rownames(Zout.selected3)[idx]
}))

df_all.z <- Zout.selected3[all_rows, , drop = FALSE]
# All ontology names you want represented
all_needed <- unique(unlist(lapply(categs, function(categ) {
  categories18[categories18[,2] == categ, 1]
})))

# Which ones are missing from df_all.z?
missing_rows <- setdiff(all_needed, sub("^([^\\.]+)\\.(.*)$", "\\2 (\\1)",rownames(df_all.z)))

if (length(missing_rows) > 0) {
  # Create a zero-filled data frame with same columns
  zero_mat <- matrix(0, nrow = length(missing_rows), ncol = ncol(df_all.z),
                     dimnames = list(missing_rows, colnames(df_all.z)))
  zero_df <- as.data.frame(zero_mat, check.names = FALSE)
  
  # Append to the end
  df_all.z <- rbind(df_all.z, zero_df)
}


old.df_all.rownames<-ALL.ontologyNames<-rownames(df_all.z)

rownames(df_all.z)<-gsub("^(.*?)\\.(.*)$", "\\2 (\\1)",old.df_all.rownames)

all(rownames(df_all.z) %in% categories18[,1])
#TRUE




# now safe to reorder -- missed ontologies in mat have NA for values.
mat <- df_all.z[categories18$NewAnnotation, as.character(wins_num), drop=FALSE]
## enforce row and column order
#mat <- mat[exactGO.141$NewAnnotation, as.character(wins_num)]

# convert to data.frame if desired
z.141 <- as.data.frame.matrix(mat)
dim(z.141)
#[1] 141  67


# Write Z scores in 141x67 matrix to CSV
write.csv(z.141[exactGO.141$NewAnnotation,],"./SEPA.all_redelivered/DOWN-141Heatmap-ordered_ontologies+Zscores.csv")

max(apply(z.141,1,function(x) length(which(is.na(x)))))
#57  out of 67 are NA at most in any of the 141 rows

all(rownames(z.141) %in% exactGO.141$NewAnnotation)
#[1] TRUE





# Write genes hit in 141x67 matrix to CSV
########################


dim(longDT.genesHit)
#[1] 271973      6  (current)




#library(data.table)

Zout.list2 <- setNames(vector("list", length(cats2)), cats2)

for (cat in cats2) {
  dtc <- longDT.genesHit[which(longDT.genesHit$ontologyType == cat),]
  if (nrow(dtc) == 0L) {
    Zout.list2[[cat]] <- data.frame()
    next
  }

  dtc$ontologyType.ontology<-paste0(dtc$ontologyType,".",dtc$ontology)
  ## --- keep only ontologies in ALL.ontologyNames
  keep_on <- intersect(unique(dtc$GOexact), categories18$NewAnnotation)

  if (length(keep_on) == 0L) {
    Zout.list2[[cat]] <- data.frame()
    next
  }

  ## --- wide matrix: rows = ontology, cols = numeric list names, values = Zscore
  wide <- dcast(dtc[GOexact %in% keep_on],
                ontology ~ win, value.var = "Genes.Hit")

  rn <- wide$ontology
  wide[, ontology := NULL]

  ## ensure columns are ordered by numeric win and named as character(as.numeric())
  col_nums <- as.numeric(names(wide))
  ord <- order(col_nums)
  wide <- wide[, ..ord]
  colnames(wide) <- as.character(col_nums[ord])

  ## finalize as data.frame with rownames
  df <- as.data.frame(wide, check.names = FALSE)
  rownames(df) <- rn

  Zout.list2[[cat]] <- df
}

## Zout.list2 is the requested output:
## - one element per ontologyType (name preserved)
## - rows = ontology strings (selected from ALL (UP+DOWN) 141, regardless of Z score max)
## - columns = as.numeric(names(GSA.FET.outSimple)) as strings
## - values = Zscore (NA if ontology absent in that window)

Zout.selected3<-do.call(rbind,Zout.list2)


### 18 categories splitting 141 of the above 144 terms. Plot each separately but on one page.
##categories18<-read.delim("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/ALL141terms_18categories.tsv",sep="\t",header=TRUE,row.names=NULL)
#categs<-c('Metabolism','RNA Processing','GLP-1 and IGFBP','Lipoproteins','Proteostasis','Synaptic and Neuronal','Eye and Retinol','Steroid Metabolism','ECM','Structural','FGF and Fibronectin','Hemostasis','Adaptive Immune','Innate Immune','Wnt signaling','Apoptosis','Translation','Angiogenesis')

library(ComplexHeatmap)
library(circlize)
library(grid)

# 1. Compute global min/max across all categories
all_rows <- unlist(lapply(categs, function(categ) {
  rn <- sub("^([^\\.]+)\\.(.*)$", "\\2 (\\1)", rownames(Zout.selected3))
  idx <- rn %in% categories18[categories18[,2] == categ, 1]
  rownames(Zout.selected3)[idx]
}))

df_all <- Zout.selected3[all_rows, , drop = FALSE]
# All ontology names you want represented
all_needed <- unique(unlist(lapply(categs, function(categ) {
  categories18[categories18[,2] == categ, 1]
})))

# Which ones are missing from df_all?
missing_rows <- setdiff(all_needed, sub("^([^\\.]+)\\.(.*)$", "\\2 (\\1)",rownames(df_all)))

if (length(missing_rows) > 0) {
  # Create a zero-filled data frame with same columns
  zero_mat <- matrix(0, nrow = length(missing_rows), ncol = ncol(df_all),
                     dimnames = list(missing_rows, colnames(df_all)))
  zero_df <- as.data.frame(zero_mat, check.names = FALSE)
  
  # Append to the end
  df_all <- rbind(df_all, zero_df)
}


rownames(df_all)<-gsub("^(.*?)\\.(.*)$", "\\2 (\\1)",old.df_all.rownames)

all(rownames(df_all) %in% categories18[,1])
#TRUE




# now safe to reorder -- missed ontologies in mat have NA for values.
mat <- df_all[exactGO.141$NewAnnotation, as.character(wins_num), drop=FALSE]
## enforce row and column order
#mat <- mat[exactGO.141$NewAnnotation, as.character(wins_num)]

# convert to data.frame if desired
genesHit.141 <- as.data.frame.matrix(mat)
dim(genesHit.141)
#[1] 141  67


# Write genes hit in 141x67 matrix to CSV
write.csv(genesHit.141,"./SEPA.all_redelivered/DOWN-Genes_Hit(141selectedOntologies_67_5yrWindows.csv")


all(rownames(genesHit.141)==exactGO.141$NewAnnotation)
#TRUE




# ---- STEP 11. Save the complete workspace image for reuse by downstream
# CMAP (prefix 8) and MAGMA (prefix 9) analyses. ----
save.image("./SEPA.all_redelivered/saved.image-67_5yearInterval_slidingWindow_GOparallel.DOWNhitsONLY-complete-7289_redelivered.RData")
#overwrites RData saved after GOparallel.


## Down ontology heat SEPA has been run as up, separately

