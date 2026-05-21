##############################################################################
# Pipeline annotation header: 6.ExtractSigAssays_2directions(ALL-up+down)_z=1.645.R
# Manuscript code section(s): 6
#
# Purpose:
# Create EYO-window significant assay lists, direction-aware ontology
# enrichment matrices, curated ontology heatmaps, GSVA term-level data, and
# waterfall visualizations.
#
# Principal inputs:
#   - simple.3177/scatter/_99_par_diff_all_peptide_p_value.rds
#   - simple.3177/scatter/_99_par_diff_all_peptide_up_down_notation.csv
#   - SEPA.all saved image
#   - Bader Lab GMT and GO OBO files
#   - _full_3177_protein_dft.RDS
#   - _numericMeta_3177_trait.RDS
#
# Principal outputs:
#   - EYO-window hit-count plots
#   - ALL_assay-Zscore_trajectories*.pdf
#   - GO/Reactome heatmaps
#   - ALL.top100.genelists.GO.forGSVA.RDS
#   - GSVA STAN outputs
#   - GSVA waterfall PDF
#
# Step overview:
#   1. Collect assays with p <= 0.005 within sliding 5-year EYO windows and
#      track directionality from up/down annotations.
#   2. Run or load GOparallel enrichment results for all, increased, and
#      decreased assay sets.
#   3. Order and visualize ontology Z-score matrices with pre-specified row
#      orders and category groupings.
#   4. Convert selected ontology gene sets into GSVA-style aggregate features.
#   5. Fit STAN/rstanarm models to ontology-level GSVA trajectories and render
#      waterfall heatmaps across EYO.
#
# Notes added during manuscript-code alignment:
#   - Original executable code was preserved except for whitespace/encoding cleanup
#     and the explicitly marked non-run archive described below when applicable.
#   - Comments inserted by this pass are prefixed with "ANNOTATION:".
##############################################################################

# ------------------------------------------------------------------------
# ANNOTATION: Load EYO x assay p-value and directionality matrices from
# longitudinal modeling outputs.
# ------------------------------------------------------------------------
#####################################################################################################################################################################################
## Start Ontologies spaghetti and heatmap plots (both directions in same lists)        ##   GSVA list generation and waterfall on MMS-PROTEOMICS-KB VM starts below line 700.
##
## - adding logic to split gene lists to up and down by also reading color annotation for waterfall input
#####################################################################################################################################################################################
setwd("f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/")
pVals<-readRDS("../simple.3177/scatter/_99_par_diff_all_peptide_p_value.rds")
dim(pVals)
#[1]  143 7347
# last protein assay is 7333.

# In R, a data frame of p values, pVals, has rows named by 1/2 year intervals from "-46", "-45.5", ... to "25" (+25). There are 7345 columns. Provide R code that collects a vector of all colnames that have at least one p value less than or equal to 0.005 within a 5-year interval starting at "-46" and ending at "-41" into a list element named by the midpoint of the 5 year interval ("-43.5"), then repeat the process incrementing the interval plus 1 year.  The end result of the code or function should be a list with 67 elements named "-43.5" to "22.5" (character strings), each a vector of column names from the pVals data frame.
## pVals: data.frame or matrix of p-values
collect_sig_by_window <- function(pVals,
                                  alpha  = 0.005,
                                  start  = -46,

# ------------------------------------------------------------------------
# ANNOTATION: Collect significant assays in sliding 5-year EYO windows.
# ------------------------------------------------------------------------
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
#-43.5 -42.5 -41.5 -40.5 -39.5 -38.5 -37.5 -36.5 -35.5 -34.5 -33.5 -32.5 -31.5 -30.5 -29.5 -28.5 -27.5 -26.5 -25.5 -24.5 -23.5 -22.5 -21.5 -20.5 -19.5 -18.5 -17.5 -16.5 -15.5
#  164   164   165   168   168   172   173   173   180   178   181   180   182   177   179   182   181   185   187   191   200   200   207   218   231   236   248   259   271
#-14.5 -13.5 -12.5 -11.5 -10.5  -9.5  -8.5  -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5   0.5   1.5   2.5   3.5   4.5   5.5   6.5   7.5   8.5   9.5  10.5  11.5  12.5  13.5
#  285   298   316   333   345   346   348   360   369   373   382   390   396   405   411   424   449   480   507   544   572   588   599   600   587   557   523   475   427
# 14.5  15.5  16.5  17.5  18.5  19.5  20.5  21.5  22.5
#  366   327   285   262   244   228   210   196   183


## Previously
#-43.5 -42.5 -41.5 -40.5 -39.5 -38.5 -37.5 -36.5 -35.5 -34.5 -33.5 -32.5 -31.5 -30.5 -29.5 -28.5 -27.5 -26.5 -25.5 -24.5 -23.5 -22.5 -21.5
#  164   164   165   168   168   172   173   173   180   178   181   180   182   177   179   182   181   185   187   191   200   200   207
#-20.5 -19.5 -18.5 -17.5 -16.5 -15.5 -14.5 -13.5 -12.5 -11.5 -10.5  -9.5  -8.5  -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5   0.5   1.5
#  218   231   236   248   259   271   285   298   316   333   345   346   348   360   369   373   382   390   396   405   411   424   449
#  2.5   3.5   4.5   5.5   6.5   7.5   8.5   9.5  10.5  11.5  12.5  13.5  14.5  15.5  16.5  17.5  18.5  19.5  20.5  21.5  22.5
#  480   507   544   572   588   599   600   587   557   523   475   427   366   327   285   262   244   228   210   196   183

assayCount<-unlist(lapply(sig_list,length))


# ------------------------------------------------------------------------
# ANNOTATION: Read direction annotations so increased and decreased assay
# signatures can be handled jointly or separately.
# ------------------------------------------------------------------------
direction_colors<-t(read.csv(file="../simple.3177/scatter/_99_par_diff_all_peptide_up_down_notation.csv",header=TRUE,row.names=1,check.names=FALSE))
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
#-43.5 -42.5 -41.5 -40.5 -39.5 -38.5 -37.5 -36.5 -35.5 -34.5 -33.5 -32.5 -31.5 -30.5 -29.5 -28.5 -27.5 -26.5 -25.5 -24.5 -23.5 -22.5 -21.5 -20.5 -19.5 -18.5 -17.5 -16.5 -15.5
#   71    70    70    70    67    67    68    67    69    68    68    66    68    69    69    69    67    68    69    73    76    76    78    83    89    93   102   112   117
#-14.5 -13.5 -12.5 -11.5 -10.5  -9.5  -8.5  -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5   0.5   1.5   2.5   3.5   4.5   5.5   6.5   7.5   8.5   9.5  10.5  11.5  12.5  13.5
#  122   130   141   150   158   160   163   172   178   183   189   200   202   211   214   221   239   260   279   304   324   330   344   349   341   323   300   271   247
# 14.5  15.5  16.5  17.5  18.5  19.5  20.5  21.5  22.5
#  206   185   156   144   139   130   120   113   104

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
#-43.5 -42.5 -41.5 -40.5 -39.5 -38.5 -37.5 -36.5 -35.5 -34.5 -33.5 -32.5 -31.5 -30.5 -29.5 -28.5 -27.5 -26.5 -25.5 -24.5 -23.5 -22.5 -21.5 -20.5 -19.5 -18.5 -17.5 -16.5 -15.5
#   98   100   100   103   106   108   109   108   112   115   114   116   117   111   113   115   116   122   122   121   126   129   133   136   144   146   152   154   157
#-14.5 -13.5 -12.5 -11.5 -10.5  -9.5  -8.5  -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5   0.5   1.5   2.5   3.5   4.5   5.5   6.5   7.5   8.5   9.5  10.5  11.5  12.5  13.5
#  164   169   177   184   190   190   188   189   194   194   194   192   196   195   199   204   212   227   231   243   253   260   260   254   247   235   226   206   180
# 14.5  15.5  16.5  17.5  18.5  19.5  20.5  21.5  22.5
#  163   143   131   120   105    99    94    88    84

## Previous
#-43.5 -42.5 -41.5 -40.5 -39.5 -38.5 -37.5 -36.5 -35.5 -34.5 -33.5 -32.5 -31.5 -30.5 -29.5 -28.5 -27.5 -26.5 -25.5 -24.5 -23.5 -22.5 -21.5 -20.5 -19.5 -18.5
#   98   100   100   103   106   108   109   108   112   115   114   116   117   111   113   115   116   122   122   121   126   129   133   136   144   146
#-17.5 -16.5 -15.5 -14.5 -13.5 -12.5 -11.5 -10.5  -9.5  -8.5  -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5   0.5   1.5   2.5   3.5   4.5   5.5   6.5   7.5
#  152   154   157   164   169   177   184   190   190   188   189   194   194   194   192   196   195   199   204   212   227   231   243   253   260   260
#  8.5   9.5  10.5  11.5  12.5  13.5  14.5  15.5  16.5  17.5  18.5  19.5  20.5  21.5  22.5
#  254   247   235   226   206   180   163   143   131   120   105    99    94    88    84


assayCountDN<-unlist(lapply(sig_listDN,length))


## Plot of Assay Counts reaching significance - All (black); Up (red); and Down (blue)
plot.new()
par(mar=c(4.5,5.5,5,2))
plot(as.numeric(names(assayCount)),assayCount, ylab="Number of dysregulated assays\nin e4/4 homozygotes vs e3/3",xlab="EYO (years)",
     main="Protein Counts Significant within 5-year Intervals relative to EYO\n(5-year Sliding Windows, 1 yr resolution)", type="b", col="black", ylim=c(0,605), lwd=2)
lines(as.numeric(names(assayCount)),assayCountUP, type="b", col="indianred3", lwd=2)
lines(as.numeric(names(assayCount)),assayCountDN, type="b", col="cornflowerblue", lwd=2)


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

setwd("SEPA.all") # empty folder for collection of all outputs - now thresholded in later below plots at min z=1.645

plot.hitCounts_plusModOverlay<-recordPlot()
pdf("3.plot.hitCounts_plusModulesHitOverlaid-07-27_5xSDoutliersNotInModels.pdf",width=22,height=8.5)
  print(plot.hitCounts_plusModOverlay)
dev.off()



# ------------------------------------------------------------------------
# ANNOTATION: Plot hit counts across EYO and overlay module/gene-set
# summaries.
# ------------------------------------------------------------------------
load("../SEPA.all(1.30zThreshold)/saved.image-67_5yearInterval_slidingWindow_GOparallel.ALLhits(up+dn_together)-complete.RData")
# SKIP TO LINE 659 (below GOparallel framework)



# ------------------------------------------------------------------------
# ANNOTATION: Load prior GOparallel enrichment results or continue with the
# embedded GOparallel framework below.
# ------------------------------------------------------------------------
##############################
## Run GOparallel  (modulesData is the list usually filled by reading input .csv with columns of lists)

inputFile <- "dummyFilename.csv"
filePath <- "f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/"   #gsub("//","/",outputfigs)
            #Folder that (may) contain the input file specified above, and which will contain the outFilename project Folder.
outFilename <- "SEPA.all"  #SUBFOLDER WITH THIS NAME WILL BE CREATED, and .PDF + .csv file using the same name will be created within this folder.
outputGOeliteInputs=FALSE  #If TRUE, GO Elite background file and module or list-specific input files will be created in the outFilename subfolder.
maxBarsPerOntology=25      #Ontologies per ontology type, used for generating the PDF report; does not limit tabled output
GMTdatabaseFile="f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/Human_GO_AllPathways_noPFOCR_with_GO_iea_June_01_2025_symbol.gmt"
GO.OBOfile<-"f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/go.obo"
modulesInMemory=FALSE
ANOVAgroups=FALSE
parallelThreads=31
WGCNAinput=FALSE

#modulesData<-sig_list.SomaOnly <- collect_sigUP_by_window(direction_colors[,1:7333])  #<-just do the UP lists here
#modulesData<-sig_list.SomaOnly <- collect_sigDN_by_window(direction_colors[,1:7333])  #<-just do the UP lists here
modulesData<-sig_list.SomaOnly <- collect_sig_by_window(pVals[,1:7333])  #previously: any significant hits (UP+DOWN)
modulesData$bkgr<-colnames(direction_colors)[1:7333]

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


# ------------------------------------------------------------------------
# ANNOTATION: Embedded GOparallel code starts here; it can generate GO-Elite
# inputs, run enrichment, and plot outputs.
# ------------------------------------------------------------------------
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

## save.image("saved.image-67_5yearInterval_slidingWindow_GOparallel.ALLhits(up+dn_together)-complete.RData")


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
pdf("ALL_assay-Zscore_trajectories.pdf", width = 14, height = 8, onefile = TRUE)
for (nm in names(Zout.list)) {
  df <- Zout.list[[nm]]
  if (!is.null(df) && nrow(df) > 0 && ncol(df) > 0) {
    print( plot_Zout_df(df, title = paste0("ALL assay Z-score trajectories: ", nm), span = 0.5) )
  }

# ------------------------------------------------------------------------
# ANNOTATION: Plot assay-level Z-score trajectories across EYO.
# ------------------------------------------------------------------------
}
dev.off()


## With a data frame from a Zout.list element as input, use the complexHeatmap::pheatmap() function to plot Z scores above 1.6450 as a color gradient that goes from black below 1.6450, to a deep dark red at 1.6450, until reaching yellow at the maximum Z score in the data.frame. The progression of color from left to right should correspond to the increments of columns in the Zout.list element data.frame for each row (ontology). Do not outline cells in the heatmap, and all missing values (NA) should be black. Do not cluster columns, but rows (ontologies) may be clustered. Show the text for ongtology labels to the left of each row, and label all columns divisible by five after subtracting 0.5. Save the plot to a page of an open PDF Z_GO_heatmaps.pdf with dimensions 8" high by 17.5" wide, and repeat for each of the 6 Zout.list elements.

library(ComplexHeatmap)

#---- helper: plot one Zout.list data frame as ComplexHeatmap::pheatmap ----------
plot_Zout_heatmap <- function(df, title = NULL, thr = 1.645) {
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

  # --- color mapping: <= 1.6450 = black; 1.6450 = deep dark red; up to max = yellow
  maxZ <- suppressWarnings(max(mat, na.rm = TRUE))
  minZ <- suppressWarnings(min(mat, na.rm = TRUE))
  if (!is.finite(maxZ)) maxZ <- thr
  if (!is.finite(minZ)) minZ <- thr

  # ensure we have a bin for everything below 1.6450 mapping to black
  minBreak <- min(minZ, thr)
  if (maxZ <= thr) {
    # nothing above threshold: still draw with two bins so 1.6450 maps to dark red
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
pdf("ALL_assaysONLY-Z_GO-heatmaps_rowsOrdered_top100.pdf", width = 17.5, height = 35, onefile = TRUE)
#grid::grid.newpage()
#par(mar=c(1,1,1,15))
#par(oma=c(1,1,1,30))
for (nm in names(Zout.list)) {
  df <- Zout.list[[nm]]
  if (!is.null(df) && nrow(df) > 0 && ncol(df) > 0) {

# ------------------------------------------------------------------------
# ANNOTATION: Render ordered ontology heatmaps for top ontology terms.
# ------------------------------------------------------------------------
    plot_Zout_heatmap(df, title = paste0("ALL ASSAYS Z-score heatmap: ", nm))
  }
}
dev.off()


## ---- select separated top terms from all six pages to put on one page : first spaghetti plot to PDF ------------------------------------------
# Zout.list: named list of 6 data.frames (GOBP, GOMF, GOCC, REACTOME, WikiPathways, MSIGDB_C2)

## Eric's 11 selected terms for UP
#selectedTerms<-c("MSIGDB_C2.Pid_hif1_tfpathway","MSIGDB_C2.Pid_cmyb_pathway","GOMF.Steroid Binding","GOMF.Lipid Transporter Activity","GOBP.Postsynapse Assembly","GOBP.Synapse Assembly","GOBP.Positive Regulation Of Camp/Pka Signal Transduction","GOMF.Glycosaminoglycan Binding",
#                 "MSIGDB_C2.Biocarta_hes_pathway","REACTOME.Pre-Notch Transcription And Translation","GOCC.Extracellular Matrix")

## Erik's 29 selected terms for ALL from top 25 lists per ontologyType
selectedTerms<-c("GOBP.Alkanesulfonate Metabolic Process","GOBP.Taurine Metabolic Process","GOBP.Hormone Metabolic Process","GOBP.Neuron Development","GOBP.Neurogenesis","GOBP.Negative Regulation Of Canonical Wnt Signaling Pathway","GOBP.Negative Regulation Of Coagulation","GOBP.Negative Regulation Of Wound Healing","GOBP.Kidney Vasculature Morphogenesis","GOMF.Deubiquitinase Activity","GOMF.Ubiquitin-Like Protein Peptidase Activity","GOMF.Extracellular Matrix Binding","GOMF.Insulin-Like Growth Factor I Binding","GOMF.Glycosaminoglycan Binding","GOMF.Heparin Binding","GOMF.Steroid Binding","GOMF.Integrin Binding","GOMF.Insulin?Like Growth Factor Ii Binding","GOMF.Cellular Response To Vascular Endothelial Growth Factor Stimulus","GOCC.Low-Density Lipoprotein Particle","GOCC.High-Density Lipoprotein Particle","GOCC.Triglyceride-Rich Plasma Lipoprotein Particle","GOCC.Very-Low-Density Lipoprotein Particle","GOCC.Chylomicron","GOCC.Extracellular Matrix",
                 "REACTOME.Heparan Sulfate Heparin (Hs-Gag) Metabolism","REACTOME.Complement Cascade","REACTOME.Ncam1 Interactions","REACTOME.Notch-Hlh Transcription Pathway")

Zout.selected<-do.call(rbind,Zout.list)
Zout.selected<-Zout.selected[which(rownames(Zout.selected) %in% selectedTerms),]

pdf("ALL_assay-Zscore_trajectories_29_SELECTED_1pp.pdf", width = 14, height = 8, onefile = TRUE)
  df <- Zout.selected
  if (!is.null(df) && nrow(df) > 0 && ncol(df) > 0) {
    print( plot_Zout_df(df, title = paste0("ALL assay Z-score trajectories: ", "Selected Terms"), span = 0.5) )
  }
dev.off()

# Replot same 20/29 selected terms as heatmap
pdf("ALL_assaysONLY-Z_GO-heatmaps_rowsOrdered_29_SELECTED_1pp.pdf", width = 17.5, height = 8, onefile = TRUE)
  if (!is.null(df) && nrow(df) > 0 && ncol(df) > 0) {
    plot_Zout_heatmap(df, title = paste0("ALL ASSAYS Z-score heatmap: ", "Selected Terms"))
  }
dev.off()


## Go back, select expanded terms from the top 100 (ECBJ, 09/06/2025), and improve plotting function ver 2
selectTerms<-c("GOBP.Steroid Catabolic Process","GOBP.Synapse Assembly","GOBP.Alkanesulfonate Metabolic Process","GOBP.Taurine Metabolic Process","GOBP.Regulation Of Synapse Organization","GOBP.Steroid Metabolic Process","GOBP.Negative Regulation Of Wnt Signaling Pathway","GOBP.Presynaptic Membrane Organization","GOBP.Energy Homeostasis","GOBP.Neuron Projection Regeneration","GOBP.Protein Neddylation","GOBP.Negative Regulation Of Mitophagy","GOBP.Protein Deubiquitination","GOBP.Postsynaptic Membrane Organization","GOBP.Rna 5'-End Processing","GOBP.Hormone Metabolic Process","GOBP.Negative Regulation Of Autophagy","GOBP.Regulation Of Canonical Wnt Signaling Pathway","GOBP.Negative Regulation Of Interferon-Alpha Production","GOBP.Cholesterol Efflux","GOBP.Plasma Lipoprotein Particle Assembly","GOBP.Protein Deneddylation","GOBP.Regulation Of Protein Neddylation","GOBP.Negative Regulation Of Neuron Projection Regeneration","GOBP.Negative Regulation Of Interleukin-12 Production","GOBP.Neuron Projection Development","GOBP.Neuron Development","GOBP.Negative Regulation Of Macroautophagy","GOBP.Positive Regulation Of Protein Monoubiquitination","GOBP.Acute-Phase Response","GOBP.Neuron Differentiation","GOBP.Positive Regulation Of Cholesterol Efflux","GOBP.Negative Regulation Of Response To Wounding","GOBP.Negative Regulation Of Canonical Wnt Signaling Pathway","GOBP.Negative Regulation Of Platelet Activation","GOBP.T Cell Proliferation","GOBP.Regulation Of Leukocyte Migration","GOBP.Negative Regulation Of Coagulation","GOBP.Regulation Of Hemostasis","GOBP.Acute Inflammatory Response","GOBP.Extracellular Matrix Organization","GOBP.Kidney Vasculature Morphogenesis","GOBP.Regulation Of Fibrinolysis","GOMF.Fibroblast Growth Factor Binding","GOMF.Deubiquitinase Activity","GOMF.Ubiquitin-Like Protein Peptidase Activity","GOMF.Heparan Sulfate Proteoglycan Binding","GOMF.Proteoglycan Binding","GOMF.Lipoprotein Particle Receptor Binding","GOMF.Oxidoreductase Activity, Acting On The Ch-Oh Group Of Donors, Nad Or Nadp As Acceptor","GOMF.Structural Constituent Of Eye Lens","GOMF.Low-Density Lipoprotein Particle Receptor Binding","GOMF.Very-Low-Density Lipoprotein Particle Receptor Binding","GOMF.Steroid Dehydrogenase Activity, Acting On The Ch-Oh Group Of Donors, Nad Or Nadp As Acceptor","GOMF.All-Trans-Retinol Dehydrogenase (Nad+) Activity","GOMF.Alcohol Dehydrogenase (Nad+) Activity","GOMF.Aminoacyl-Trna Ligase Activity","GOMF.Opsonin Binding","GOMF.Sterol Transport","GOMF.Lipid Transfer Activity","GOMF.Neurexin Family Protein Binding","GOMF.[Heparan Sulfate]-Glucosamine 3-Sulfotransferase Activity","GOMF.Heparan Sulfate Sulfotransferase Activity","GOMF.Lipid Transport","GOMF.Fibronectin Binding","GOMF.Extracellular Matrix Binding","GOMF.Insulin-Like Growth Factor I Binding","GOMF.Complement Binding","GOMF.Insulin-Like Growth Factor Binding","GOMF.Complement Component C3b Binding","GOMF.Glycosaminoglycan Binding","GOMF.Heparin Binding","GOMF.Integrin Binding","GOMF.Insulin-Like Growth Factor Ii Binding","GOMF.Amylase Activity","GOMF.Extracellular Matrix Structural Constituent","GOMF.Hyaluronic Acid Binding","GOMF.Vascular Endothelial Growth Factor Receptor Activity","GOCC.Endocytic Vesicle Lumen","GOCC.Intermediate-Density Lipoprotein Particle","GOCC.Spliceosomal Tri-Snrnp Complex","GOCC.U4/U6 X U5 Tri-Snrnp Complex","GOCC.Precatalytic Spliceosome","GOCC.U2-Type Precatalytic Spliceosome","GOCC.U4 Snrnp","GOCC.Neurofilament","GOCC.U1 Snrnp","GOCC.Mhc Class Ii Protein Complex","GOCC.Sno(S)Rna-Containing Ribonucleoprotein Complex","GOCC.Ubiquitin Ligase Complex","GOCC.Spliceosomal Snrnp Complex","GOCC.Death-Inducing Signaling Complex","GOCC.Low-Density Lipoprotein Particle","GOCC.High-Density Lipoprotein Particle","GOCC.Intermediate Filament Cytoskeleton","GOCC.Triglyceride-Rich Plasma Lipoprotein Particle","GOCC.Very-Low-Density Lipoprotein Particle","GOCC.Chylomicron","GOCC.Lysosomal Lumen","GOCC.Synaptic Cleft","GOCC.Perikaryon","GOCC.Synaptobrevin 2-Snap-25-Syntaxin-1a-Complexin Ii Complex","GOCC.Extracellular Matrix",
               "GOCC.Blood Microparticle","GOCC.Axon Initial Segment","GOCC.Membrane Attack Complex","GOCC.Vacuolar Lumen","GOCC.Azurophil Granule Lumen","GOCC.Perisynaptic Extracellular Matrix","GOCC.Microfibril","GOCC.Schaffer Collateral - Ca1 Synapse","GOCC.Late Endosome","GOCC.Lysosome","GOCC.Insulin-Like Growth Factor Binding Protein Complex","REACTOME.The Canonical Retinoid Cycle In Rods (Twilight Vision)","REACTOME.Incretin Synthesis, Secretion, And Inactivation","REACTOME.Akt-Mediated Inactivation Of Foxo1a","REACTOME.Regulation Of Tnfr1 Signaling","REACTOME.Synthesis, Secretion, And Inactivation Of Glucose-Dependent Insulinotropic Polypeptide (Gip)","REACTOME.Synthesis, Secretion, And Inactivation Of Glucagon-Like Peptide-1 (Glp-1)","REACTOME.Deubiquitination","REACTOME.Ligand-Dependent Caspase Activation","REACTOME.Wnt Ligand Biogenesis And Trafficking","REACTOME.Ra Biosynthesis Pathway","REACTOME.Cytosolic Trna Aminoacylation","REACTOME.Ldl Remodeling","REACTOME.Hs-Gag Biosynthesis","REACTOME.Activation Of Matrix Metalloproteinases","REACTOME.Heparan Sulfate Heparin (Hs-Gag) Metabolism","REACTOME.Interleukin-4 And Interleukin-13 Signaling","REACTOME.Complement Cascade","REACTOME.Regulation Of Complement Cascade","REACTOME.Ncam1 Interactions","REACTOME.Neurofascin Interactions","REACTOME.Notch-Hlh Transcription Pathway","REACTOME.Interleukin-12 Signaling","REACTOME.L1cam Interactions","REACTOME.Nrcam Interactions","REACTOME.Terminal Pathway Of Complement","REACTOME.Glycosaminoglycan Metabolism","REACTOME.Extracellular Matrix Organization","REACTOME.Intrinsic Pathway Of Fibrin Clot Formation","REACTOME.Sema3a-Plexin Repulsion Signaling By Inhibiting Integrin Adhesion","REACTOME.Neurophilin Interactions With Vegf And Vegfr")
#note: selectedTerms above is length=29
length(selectTerms)
#144
Zout.selected2<-do.call(rbind,Zout.list)
Zout.selected2<-Zout.selected2[which(rownames(Zout.selected2) %in% selectTerms),]

plot_Zout_heatmap2 <- function(mat, title=NULL, thr=1.645, rn = rownames(mat), draw=TRUE, cols=NULL, breaks=NULL,...) {
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
  ord_rows <- order(!all_sig, first_hit, -first_hit_val,
                    -first_non_sig_after_hit, rn)

  mat <- mat[ord_rows, , drop = FALSE]

  ## --- COLOR MAPPING ---
  if (is.null(cols) | is.null(breaks)) {
    maxZ <- suppressWarnings(max(mat, na.rm = TRUE))
    minZ <- suppressWarnings(min(mat, na.rm = TRUE))
    if (!is.finite(maxZ)) maxZ <- thr
    if (!is.finite(minZ)) minZ <- thr

    minBreak <- min(minZ, thr)
    if (maxZ <= thr) {
      breaks <- c(minBreak, thr, thr + 1e-6)
      cols   <- c("#FFFFFF", "#C9C9C9")  #c("black", "#8B0000")
    } else {
      n_grad <- 100L
      breaks <- c(0, 0.000001, thr-0.000001, thr, seq(thr, maxZ, length.out = n_grad + 1L)[-1L])     ## c(seq(minBreak, thr-0.000001, length.out=as.integer(n_grad*(100*thr/maxZ))+1L)[-1L], thr, seq(thr, maxZ, length.out = as.integer(n_grad*(100-(100*thr/maxZ))) + 1L)[-1L])
      cols   <- c("#FFFFFF","#C9C9C9","#C9C9C9", colorRampPalette(c("#8B338B", "yellow"))(n_grad))  ## c(colorRampPalette(c("#CCCCCC","#333333"))(as.integer(n_grad*(100*thr/maxZ))), colorRampPalette(c("#8B3333", "yellow"))(as.integer(n_grad*(100-(100*thr/maxZ)))))
    }
  }

  ## --- DRAW HEATMAP ---
  ht <- ComplexHeatmap::pheatmap(
    mat,
    color          = cols,
    breaks         = breaks,
    na_col         = "black",  # does not work
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

# plot gold-purple-grey-white heatmap
pdf("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/SEPA.all/3.SEPAwindows(67)ALL(top100)_assays-Z_GO-heatmaps_rowsOrdered_144_SELECTED_1pp_greyBG.pdf", width = 17.5, height = 32, onefile = TRUE)
  df <- Zout.selected2
  rownames(df) <- sub("^([^\\.]+)\\.(.*)$", "\\2 (\\1)", rownames(df))
  if (!is.null(df) && nrow(df) > 0 && ncol(df) > 0) {
    df[is.na(df)] <- 0
    plot_Zout_heatmap2(df, title = paste0("ALL ASSAYS Z-score heatmap: ", "144 Selected Terms"))
  }
dev.off()


## 18 categories splitting 141 of the above 144 terms. Plot each separately but on one page.
categories18<-read.delim("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/ALL141terms_18categories.tsv",sep="\t",header=TRUE,row.names=NULL)
categs<-c('Metabolism','RNA Processing','GLP-1 and IGFBP','Lipoproteins','Proteostasis','Synaptic and Neuronal','Eye and Retinol','Steroid Metabolism','ECM','Structural','FGF and Fibronectin','Hemostasis','Adaptive Immune','Innate Immune','Wnt signaling','Apoptosis','Translation','Angiogenesis')

library(ComplexHeatmap)
library(circlize)
library(grid)

# 1. Compute global min/max across all categories
all_rows <- unlist(lapply(categs, function(categ) {
  rn <- sub("^([^\\.]+)\\.(.*)$", "\\2 (\\1)", rownames(Zout.selected2))
  idx <- rn %in% categories18[categories18[,2] == categ, 1]
  rownames(Zout.selected2)[idx]
}))

df_all <- Zout.selected2[all_rows, , drop = FALSE]

df_all[is.na(df_all)] <- 0

global_minZ <- min(df_all, na.rm = TRUE)
global_maxZ <- max(df_all, na.rm = TRUE)

# 2. Define one shared color function
thr=1.645
if (!is.finite(global_maxZ)) global_maxZ <- thr
if (!is.finite(global_minZ)) global_minZ <- thr

minBreak <- min(global_minZ, thr)
if (global_maxZ <= thr) {
  breaks <- c(minBreak, thr, thr + 1e-6)
  cols   <- c("#FFFFFF", "#C9C9C9")  # c("black", "#8B0000")
} else {
  n_grad <- 100L
#  breaks <- c(seq(minBreak, thr-0.000001, length.out=as.integer(n_grad*(100*thr/global_maxZ))+1L)[-1L], thr, seq(thr, global_maxZ, length.out = as.integer(n_grad*(100-(100*thr/global_maxZ))) + 1L)[-1L])  # c(minBreak, thr, seq(thr, global_global_maxZ, length.out = n_grad + 1L)[-1L])
  breaks <- c(0, 0.000001, thr-0.000001, thr, seq(thr, global_maxZ, length.out = n_grad + 1L)[-1L])     ## c(seq(minBreak, thr-0.000001, length.out=as.integer(n_grad*(100*thr/maxZ))+1L)[-1L], thr, seq(thr, maxZ, length.out = as.integer(n_grad*(100-(100*thr/maxZ))) + 1L)[-1L])
#  cols   <- c(colorRampPalette(c("#CCCCCC","#333333"))(as.integer(n_grad*(100*thr/global_maxZ))), colorRampPalette(c("#8B3333", "yellow"))(as.integer(n_grad*(100-(100*thr/global_maxZ)))))  # c("black", colorRampPalette(c("#8B0000", "yellow"))(n_grad))
  cols   <- c("#FFFFFF","#C9C9C9","#C9C9C9", colorRampPalette(c("#8B338B", "yellow"))(n_grad))  ## c(colorRampPalette(c("#CCCCCC","#333333"))(as.integer(n_grad*(100*thr/maxZ))), colorRampPalette(c("#8B3333", "yellow"))(as.integer(n_grad*(100-(100*thr/maxZ)))))
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

ht_list <- NULL
order141<-vector()

for (categ in categs) {
  df <- Zout.selected2
  rownames(df) <- sub("^([^\\.]+)\\.(.*)$", "\\2 (\\1)", rownames(df))
  df <- df[rownames(df) %in% categories18[categories18[,2] == categ, 1], ]
  df[is.na(df)] <- 0
  order141<-c(order141,rownames(df))

  if (nrow(df) > 0 && ncol(df) > 0) {
    title_ht <- make_title_ht(categ, ncol(df))
    ht <- plot_Zout_heatmap2(df, draw = FALSE, breaks = breaks, cols = cols)
    section <- title_ht %v% ht
    ht_list <- if (is.null(ht_list)) section else ht_list %v% section
  }
}


pdf("F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/SEPA.all/3.SEPAwindows(67)ALL(top100)_18categoriesSeparated(141assays)-Z_GO-heatmaps_rowsOrdered_greyBG.pdf", width = 20, height = 30, onefile = TRUE)
  draw(ht_list, merge_legend = TRUE)
dev.off()


## Save RData with categories18 and df_all (141 up+down Z score n=67 5-year window trajectories -- for aligned use in up only and down only processing
save(categories18,df_all,order141,file="F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/SEPA.all/ALL.df_all-141ontologies.categories18.RData")


## ALL run complete; DOWN would be run as above for UP


#####################################################################################################################################################################################
## Start GSVA
#####################################################################################################################################################################################
## Extract top 100 ontologies hit (for each category) in the 67 5 year windows' genes reaching significance for e4/4 difference from e3/3 in STAN 99% CI ribbons
library(data.table)

cats <- c("GOBP","GOMF","GOCC","REACTOME","WikiPathways","MSIGDB_C2")

Zout.100.list <- setNames(vector("list", length(cats)), cats)

for (cat in cats) {
  dtc <- longDT[ontologyType == cat]
  if (nrow(dtc) == 0L) {
    Zout.100.list[[cat]] <- data.frame()
    next
  }

  ## --- pick top 100 ontologies by their best (max) Z per ontology across all wins
  maxZ <- dtc[, .(maxZ = max(Zscore)), by = ontology][order(-maxZ)]
  if (nrow(maxZ) == 0L) {
    Zout.100.list[[cat]] <- data.frame()
    next
  }
  n_top  <- min(100L, nrow(maxZ))
  thresh <- maxZ$maxZ[n_top]                     # include ties at the 100th rank
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

  Zout.100.list[[cat]] <- df
}


ontologies.forGSVA <- lapply(Zout.100.list,rownames)
# GSCfromGMT$gsc is the GOparallel list of 29000 ontologies (names), genes in vector elements
genelists.GO.forGSVA<-vector("list", length(ontologies.forGSVA))
names(genelists.GO.forGSVA) <- names(ontologies.forGSVA)

gsc_names <- names(GSCfromGMT$gsc)


for (cat in names(genelists.GO.forGSVA)) {
  # ontology names in this category (ensure a plain character vector)
  onto_vec <- ontologies.forGSVA[[cat]]

  idx<-vector()
  gene_lists<-vector("list", length(ontologies.forGSVA[[cat]]))
  names(gene_lists)<-ontologies.forGSVA[[cat]]
  for (onto in onto_vec) {
    # build exact keys like "MY_ONTOLOGY%Category"
    onto<-gsub("\\(","\\\\\\(", onto)  # replace open parentheses with escaped version "\\("...
    onto<-gsub("\\)","\\\\\\)",onto)
    onto<-gsub("\\+","\\\\\\+",onto)
    onto<-gsub("\\[","\\\\\\[",onto)
    onto<-gsub("\\]","\\\\\\]",onto)
    keys <- paste0(toupper(onto), "%", toupper(cat))

    #grepl match into the 29k-entry gsc list
    idx <- c(idx, which(grepl(paste0("^",keys),names(GSCfromGMT$gsc))) )

    gene_lists<-GSCfromGMT$gsc[idx]
  }

  # handle any not-found ontologies gracefully
  if (anyNA(idx)) {
    gene_lists[is.na(idx)] <- list(character(0))
    warning(sprintf("In category '%s', %d ontology(ies) not found in GSC; set to character(0).",
                    cat, sum(is.na(idx))))
  }

  genelists.GO.forGSVA[[cat]] <- gene_lists
}

lapply(genelists.GO.forGSVA, length)
#$GOBP
#[1] 100
#$GOMF
#[1] 101
#$GOCC
#[1] 100
#$REACTOME
#[1] 100
#$WikiPathways
#[1] 102
#$MSIGDB_C2
#[1] 101

genelists.GO.forGSVA<-lapply(genelists.GO.forGSVA, function(cat_list) {
  names(cat_list)=stringr::str_to_title(gsub("\\%WIKIPATHWAYS_\\d*","", gsub("\\%WP_\\d*","", gsub("\\&(.*);","\\1",gsub("<\\sI>","",gsub("<I>","", gsub("(.*)\\%.*\\%.*","\\1", names(cat_list) )))))))
  cat_list })

saveRDS(genelists.GO.forGSVA,"f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/ALL.top100.genelists.GO.forGSVA.RDS")


# ## On VM Windows03 (AD Workbench)
# setwd("z:/ShijiaBian/PlasmaProteomic/Result/20250727/
# genelists.GO.forGSVA<-readRDS("genelists.GO.forGSVA.RDS")

#setwd("z:/ShijiaBian/PlasmaProteomic/Result/20250727/simple.3177/")
library(tidyverse)
library(forcats)
library(rstanarm)
library(Hmisc)
library(openxlsx)
library(rstan)

# ------------------------------------------------------------------------
# ANNOTATION: Save selected ontology gene lists for downstream GSVA
# analyses.
# ------------------------------------------------------------------------
#numericMeta_3177 trait <- readRDS("~/files/EBD/Shijia_B_Derived_Data/20250709/numericMeta_3177_trait.RDS")
#full_3177_protein_df <- readRos(·~/files/EBD/Shijia_B_Derived_Data/20250709/full_3177_protein_dft.RDS")
# load("z:/EBD/APOE_homozygote_2ndFinal_cleanDats+traits(3177.CT33only_4199.33inADcontinuum).RData")
# ls()
#"cleanDat.3177. "cleanDat.4199"  "MEs.3177" "MEs.4199" "numericMeta.3177" "numericMeta.4199"
#full_3177_protein_df <- as.data.frame(rbind(cleanDat.3177, t(MEs.3177)))  #not used
cleanDat.3177<- t(readRDS("../_full_3177_protein_dft.RDS"))
numericMeta.3177<- readRDS("../_numericMeta_3177_trait.RDS")

## Collapse 3177 7333 assays to unique gene product assays (maxVar)
UniqueIDs<-rownames(cleanDat.3177)
Symbols.7333<-as.data.frame(do.call("rbind",strsplit(UniqueIDs,"[|]")))[,1]
library(WGCNA)
cleanDat.3177.collapsed<-collapseRows(cleanDat.3177,rowGroup=Symbols.7333,rowID=UniqueIDs,method="maxRowVariance")
cleanDat.3177.collapsed<-cleanDat.3177.collapsed$datETcollapsed
symbols.collapsed<-as.data.frame(do.call("rbind",strsplit(rownames(cleanDat.3177.collapsed),"[|]")))[,1]
#cleanRelAbun.Unreg.Collapsed.ENSGIDs<-EnsgIDlookup[match(symbols.collapsed,EnsgIDlookup$hgnc_symbol),"ensembl_gene_id"]
#rownames(cleanRelAbun.Unreg.Collapsed)<-cleanRelAbun.Unreg.Collapsed.ENSGIDs
dim(cleanDat.3177.collapsed)
# 6408 3177  # previously: 6393 3177


library(GSVA)
library(BiocParallel)

## 1) Expression data: genes (rows) x samples (cols)

# ------------------------------------------------------------------------
# ANNOTATION: Build GSVA-style ontology score matrices from the clean
# 3,177-sample protein dataset.
# ------------------------------------------------------------------------
expr <- as.matrix(cleanDat.3177.collapsed)

## 2) Flatten the 6 category lists into one gene-set list
##    - Keep names unique by prefixing with the category
##    - Ensure each element is a character vector of gene symbols
flatten_sets <- function(glists) {
  out <- list()
  for (cat in names(glists)) {
    sub <- glists[[cat]]
    # ensure sub-elements are named; if not, make a fallback name
    if (is.null(names(sub)) || any(names(sub) == "" | is.na(names(sub)))) {
      names(sub) <- paste0("SET_", seq_along(sub))
    }
    # prefix ontology names with category
    new_names <- paste0(cat, "_", names(sub))
    # coerce each to unique character vector
    sub <- lapply(sub, function(v) unique(as.character(v)))
    names(sub) <- new_names
    out <- c(out, sub)
  }
  out
}

gene_sets_all <- flatten_sets(genelists.GO.forGSVA)

## Merge unique symbols in 2 duplicate named gene sets (431, 432: WikiPathways_Estrogen Metabolism); (472, 473) WikiPathways_Notch Signaling)
gene_sets_all[[472]]<-unique(c(gene_sets_all[[472]],gene_sets_all[[473]]))
gene_sets_all[[431]]<-unique(c(gene_sets_all[[431]],gene_sets_all[[432]]))

gene_sets_all<-gene_sets_all[which(!duplicated(names(gene_sets_all)))]
length(gene_sets_all)
#602


## For RDS loaded cleanDat collapsed only:
max(apply(expr,1,function(x) length(which(is.na(x)))))
#3177
which(apply(expr,1,function(x) length(which(is.na(x))))==3177)
#sample_id
#     5014
expr<-expr[which(!rownames(expr) %in% c("sample_id","MMSE","cdr")),]
dim(expr)
#[1] 6405 3177
which(apply(expr,1,function(x) length(which(is.na(x))))==3177)
#named integer(0)
max(apply(expr,1,function(x) length(which(is.na(x)))))
#52  #row with 52/3177 missing 5xSD outliers set to NA

# Count rows with zero variance (all values identical, including NAs ignored)
expr_mat <- as.matrix(expr)
# Get ranks per sample (ties averaged)
ranks <- apply(expr_mat, 2, rank, ties.method = "average")
# Genes whose rank vector is identical across samples
constant_rank <- apply(ranks, 1, function(x) sd(x) == 0)

sum(constant_rank)  # This should match GSVA's reported constant count


## 3) Filter gene sets by overlap with your expression matrix
min_size <- 3  # choose as you prefer
gene_sets_filtered <- lapply(gene_sets_all, function(gs) intersect(gs, rownames(expr)))
gene_sets_filtered <- gene_sets_filtered[lengths(gene_sets_filtered) >= min_size]

length(gene_sets_all)      # total gene sets provided
length(gene_sets_filtered) # gene sets that overlap & pass size filter
#602, in both cases (if min_size=5, 550)

## 4) Parallel setup (choose one appropriate for your OS)
## On Linux/macOS:
# bpp <- MulticoreParam(workers = 8)
## On Windows:
bpp <- SnowParam(workers = 31, type = "SOCK")

## 5) Run GSVA
## For log2 proteomics, kcdf="Gaussian" is typically appropriate.
## method can be "gsva" (default) or "ssgsea" if you prefer faster, rank-based scoring.

## NEWER version on R v4.5.1
#gsva_es <- gsva(gsvaParam(expr, gene_sets_filtered, kcdf="Gaussian", maxDiff=TRUE),
#                BPPARAM  = bpp)
## On R v4.2.3
gsva_es <- gsva(impute::impute.knn(expr)$data, gene_sets_filtered, kcdf="Gaussian", mx.diff=TRUE,
                method   = "gsva",
                min.sz   = min_size,
                max.sz   = Inf,
                BPPARAM  = bpp)

## gsva_es is a matrix: gene-sets (rows) x samples (columns)
dim(gsva_es)

## 6) (Optional) Split results back by category, if helpful
#cat_of <- sub("_.*$", "", rownames(gsva_es))
#gsva_by_category <- split(gsva_es, f = cat_of)  # each element is a matrix of sets for that category


## Data cleaning - 5SD from mean max within protein, then max 20% per genotype group NA; add MMSE and cdr
clean_mat <- as.data.frame(t(apply(gsva_es, 1, function(row) {
    ## winsorise at ± 5 SD
  z   <- abs(row - mean(row, na.rm = TRUE))
  row[z > 5 * sd(row, na.rm = TRUE)] <- NA                 # outliers -> NA

  row }                                              # otherwise return the cleaned row
)))
## restore dimnames -------------------------------------------------------
rownames(clean_mat) <- rownames(gsva_es)
colnames(clean_mat) <- colnames(gsva_es)


## grouping factor - one value per column in the expression matrix
grp <- numericMeta.3177$APOE.mapped.predicted
stopifnot(length(grp) == ncol(gsva_es))   # safety

## pre-compute the column indices that belong to each group
idx_by_grp <- split(seq_along(grp), grp)

## rows that fail the “<=20% NA in every group” rule  ----------------------
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


## Add MMSE, cdr as additional outcomes
full_3177_protein_df<-as.data.frame(t(clean_mat))
#full_3177_protein_df$MMSE<-numericMeta.3177$MMSE
#full_3177_protein_df$cdr<-numericMeta.3177$cdr

full_3177_protein_df$sample_id <- rownames(full_3177_protein_df)


# set derived traits
numericMeta.3177$ApoE_Indicator<-0
numericMeta.3177$ApoE_Indicator[numericMeta.3177$APOE.mapped.predicted=="e4/e4"]<-1
table(numericMeta.3177$ApoE_Indicator)
#   0   1
#2764 413
numericMeta.3177$EY0 <- (65.6-numericMeta.3177$age_at_visit)*(-1)
range(numericMeta.3177$EY0)
# -45.6 24.4
numericMeta_3177_trait <- numericMeta.3177
numericMeta_3177_trait$sample_id <- rownames(numericMeta_3177_trait)

which(!full_3177_protein_df$sample_id==numericMeta_3177_trait$sample_id)
# integer(O)

########################################################################
#' @param x A data frame (or any object with names)
#' @return The same object, but with names cleaned
sanitize_names <- function(x) {
  #Capture existing names
  nm <- names(x)
  if (is.null(nm)) return(x)

  # Replace anything that's not a letter or digit with
  clean_nm <- gsub("[^[:alnum:]]", "_", nm)
  # (Optional) collapse repeated underscores:
  clean_nm <- gsub("_+", "_", clean_nm)
  # (Optional) trim leading/trailing underscores:
  clean_nm <- gsub("^_|_$", "", clean_nm)

  names(x) <- clean_nm
  x
}
##########################################################################

# For the protein file, make the person_id to be the first column
protein_df <- full_3177_protein_df %>%
  select(sample_id, everything())  #  Moves 'AAA' to the first column


pep_names <- names(protein_df)[2:dim(protein_df)[2] ]


one_pepSTAN <- function(track, traits_df, protein_df) {
  pep_name=pep_names[track]
  ##### --------------------- Format the data frame -------------------- #####
  dat <- traits_df %>%
    select(sample_id, ApoE_Indicator, EY0) %>%
    left_join(
      protein_df %>% select(sample_id, all_of(pep_name)),
      by = "sample_id") %>%
    filter(complete.cases(.)) %>%
    select(-contains("sample_id")) # %>%
#    rename(., Sex = Sex.int)

  dat <- sanitize_names(dat)

  splinefit = rcspline.eval(dat$EY0, nk=3, norm = 2, pc = FALSE, inclx=TRUE)
  cubic_spline_X <- as.data.frame(splinefit)
  names(cubic_spline_X) <- c("EYO_Spline_Linear", "EYO_Spline_Cubic")
  #head(cubic_spline_X)

  dat <- cbind(dat, cubic_spline_X)

  # Model 1: Construct the formula
  outcome <- names(dat)[3]
  variables1 <- c("EYO_Spline_Linear", "EYO_Spline_Cubic",
                  "ApoE_Indicator", "EYO_Spline_Linear*ApoE_Indicator", "EYO_Spline_Cubic*ApoE_Indicator")
  f1 = as.formula(paste(outcome, paste(variables1, collapse = "+"), sep = "~"), env=baseenv())
  # * environment control needed for within-function stan_glm() call, per https://discourse.mc-stan.org/t/stanfit-object-fit-inside-function-explodes-in-size-when-saved-to-rds/13656
  set.seed(track)
  env <- new.env() #parent = .GlobalEnv)
  env$dat <- dat
  env$f1 <- f1
  env$track <- track
  stan_BL_1 <- with(env, {stan_glm(f1,
                                 data= dat,
                                 family=gaussian(),
                                 chains = 8,
                                 cores = 4,
                                 iter = 10000,
                                 thin = 10,
                                 refresh = 0,       # silence sampler
                                 seed = track) })  #, warning=function(w) {print(paste("Warning: ", track))})})
## On Windows03 (Workbench)
#  fil_stan_1 <- file.path("z:/","ShijiaBian","PlasmaProteomic","Result","20250727", "simple.3177.GSVA", paste(paste(outcome, "_stan_glm", ".rds", sep = "")))
  fil_stan_1 <- file.path("f:/","OneDrive - Emory","Legacy","e4_homozygoteStudy","DL", "3.Five_yr_slidingWindow","GSVA", paste(paste(outcome, "_stan_glm", ".rds", sep = "")))

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
  stan_BL_1 <- strip_formula_envs(stan_BL_1)

  saveRDS(stan_BL_1, fil_stan_1)

  #name_match_table_temp <-
  data.frame(OriginalName = pep_name,
             CleanedName = colnames(dat)[3],
             stringsAsFactors=FALSE)
}

library(doParallel)
ncore <- 8  #max(1, parallel::detectCores() - 1)
cl <- makeCluster(ncore)
registerDoParallel(cl)

library(foreach)
worker_pkgs <- c("tidyverse","rstanarm","Hmisc","rstan","dplyr")

results <- foreach(track = 1:length(pep_names),
                   .packages = worker_pkgs,
                   .export = c("one_pepSTAN","sanitize_names",    # the functions
                               "pep_names"),
                   .errorhandling = "pass")  %dopar%  {
  one_pepSTAN(track, numericMeta_3177_trait, protein_df)
}

stopImplicitCluster()

name_match_table<-do.call(rbind,results)
colnames(name_match_table)<-c("OriginalName","CleanedName")

## on Workbench Windows03
#saveRDS(name_match_table,"z:/ShijiaBian/PlasmaProteomic/Result/20250727/simple.3177.GSVA/name_match_table.RDS")
#saveRDS(numericMeta_3177_trait,"z:/ShijiaBian/PlasmaProteomic/Result/20250727/simple.3177.GSVA/_numericMeta_3177_trait.RDS")
#saveRDS(protein_df,"z:/ShijiaBian/PlasmaProteomic/Result/20250727/simple.3177.GSVA/_full_3177_protein_dft.GSVA.RDS")
saveRDS(name_match_table,"f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/GSVA/name_match_table.RDS")
saveRDS(numericMeta_3177_trait,"f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/GSVA/_numericMeta_3177_trait.RDS")

# ------------------------------------------------------------------------
# ANNOTATION: Save rstanarm model fits for ontology-level GSVA trajectories.
# ------------------------------------------------------------------------
saveRDS(protein_df,"f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/GSVA/_full_3177_protein_dft.GSVA.RDS")

## END OF MODEL GENERATION


## START SCATTERPLOTTING
library(tidyverse)
library(rstanarm)
library(Hmisc)
library(openxlsx)
library(rstan)
library(gridExtra)
library(ggpubr)

# Comment and un-comment on HPC
## On AD Workbench Windows03
#setwd("z:/ShijiaBian/PlasmaProteomic/Result/20250727/simple.3177.GSVA/")
setwd("f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/GSVA/")

##################### ------------ Read the Trait Data ------------ ##################
# Load the master traits data and the pep2pro data
#BL_traits <- readRDS("/home/workspace/files/EBD/Shijia_B_Derived_Data/20250709/numericMeta_3177_trait.RDS")
BL_traits <- readRDS("_numericMeta_3177_trait.RDS")
BL_traits$EYO<- (65.6 - BL_traits$age_at_visit)*(-1)
#BL_traits_pep<- readRDS("/home/workspace/files/EBD/Shijia_B_Derived_data/20250709/full_3177_protein_dft.RDS")
BL_traits_pep<- readRDS("_full_3177_protein_dft.GSVA.RDS")
BL_traits_pep$EYO<-BL_traits$EYO
BL_traits_pep$ApoE_Indicator<-BL_traits$ApoE_Indicator

min(BL_traits_pep$EYO) # -45.6     -31.53973
max(BL_traits_pep$EYO) # 24.4       23.72 previously
EYO_cut = length(seq(-46, 25, by = 0.5)) # 142   #prev 113 for DS

name_match_table <- readRDS("name_match_table.rds")
name_module_label <- read.csv("../../scatterplot_label_20250718.csv", header = T)
name_module_label.add<-cbind(name_match_table$CleanedName,name_match_table$OriginalName,name_match_table$OriginalName)
colnames(name_module_label.add)<-colnames(name_module_label)
name_module_label<-rbind(name_module_label,name_module_label.add)


################## -------- Fit the STAN model -------- ###########

pep_names <- names(BL_traits_pep)[c(2:603)]  #EYO and ApoE_Indicator at end (605),

## On Workbench Windows03
#source("z:/ShijiaBian/PlasmaProteomic/Code/CommonFunctions/plot_functions_20250718.R")   #202311235.R  changed from _Der_F to ApoE_Indicator
source("f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/plot_functions_20250718.R")


process_one_peptide <- function(pep_name) {
	clean_name = name_match_table$CleanedName[which(name_match_table$OriginalName == pep_name)]
	fil_stan <- file.path(".", paste(paste(clean_name, "_stan_glm",  ".rds", sep = "")))

	stan_BL <- readRDS(fil_stan)

	if (length(stan_BL) == 1) {
	  next
	}

	###### ----------- Plot 1: Two separate plots for carrier and non-carrier ------------- ############
	parameter_estimates <- rstan::extract(stan_BL$stanfit) # Extract the model estimates from STAN

	# biomarker ~ EYO_Spline_Linear + EYO_Spline_Cubic + Group + EYO_Spline_Linear * Group + EYO_Spline_Cubic
	# Add intercept to beta weights, 4000 * 6
	factorweights <- cbind(parameter_estimates$alpha, parameter_estimates$beta) # Add intercept to beta weights

	#Initialize contrast matrices
	eyo_step = seq(floor(min(BL_traits_pep$EYO, na.rm = T)),
				   ceiling(max(BL_traits_pep$EYO, na.rm = T)),
				   by=.5) #create a vector from x to y by the specified interval. This represents the range that we will use to plot our lines

	#Generate output matrices, nrow = 4,000 (stan iterations); ncol = 125, 125 is the number of eyo
	spline_noncarriers = matrix(0, nrow = nrow(factorweights), ncol = length(eyo_step)) #generate blank matrices for outputs for noncarriers
	spline_carriers = matrix(0, nrow = nrow(factorweights), ncol = length(eyo_step)) #generate blank matrices for outputs for carriers
	spline_diff = matrix(0, nrow = nrow(factorweights), ncol = length(eyo_step)) #generate blank matrices for output of differences

	contrasts_non_carriers = matrix(0, nrow = length(eyo_step), ncol = dim(factorweights)[2])
	contrasts_carriers = matrix(0, nrow = length(eyo_step), ncol = dim(factorweights)[2]) #Making a blank matrix for contrasts
	splinefit = rcspline.eval(BL_traits_pep$EYO, nk=3, norm=2, pc=FALSE, inclx=TRUE) # Redo the spline fit (if wanting to skip earlier portions, currently omitted)

	for (j in 1:length(eyo_step)) {
		tempfit = rcspline.eval(eyo_step[j], attr(splinefit,'knots'), norm=2, pc=FALSE, inclx=TRUE)       #put all of the values from our plotting interval into rcscpline to get the cubic term values
		#Contrast
#		contrasts = matrix(
#			c(1, 1, tempfit[1,1],tempfit[1,1],tempfit[1,2],tempfit[1,2], 0,1, 0, tempfit[1,1], 0,tempfit[1,2]),
#			nrow=2,
#			ncol=6)
#		contrasts_non_carriers[j,] = contrasts[1,]
#		contrasts_carriers[j,] = contrasts[2,] #Set up matrix of contrasts for each EYO point separately for carriers and noncarriers to be used below

##              ORIGINAL MODEL
#		biomarker ~ EYO_L + EYO_C + Group + EYO_L:Group + EYO_C:Group
		# tempfit gives the spline bases at the current EYO value
		xL <- tempfit[1, 1]
		xC <- tempfit[1, 2]

		# prediction for NON-carriers (Group = 0)
		contrasts_non_carriers[j, ] <- c(1,  xL,  xC,  0,      0,        0)

		# prediction for carriers (Group = 1)
		contrasts_carriers    [j, ] <- c(1,  xL,  xC,  1,      xL,       xC)

	}

	#Begin the loop to go through every single iteration of Stan
	for (i in 1:nrow(factorweights)) {
		weights_temp = factorweights[i,] #for every iteration select the weights

		# for every EYO in our output what would the value be using this iterations model fits
		for (j in 1:length(eyo_step)) {
			#select the appropriate contrast from the matrix made above
			contrasts = rbind(contrasts_non_carriers[j,], contrasts_carriers[j,])

			# Multiple beta weights by the contrast to get one value per group for this EYO and this iteration
			weights_point = contrasts %*% weights_temp

			# Write the points out into a matrix for each group
			spline_noncarriers[i,j] = weights_point[1,]
			spline_carriers[i,j] = weights_point[2,]

			#difference is equal to carriers - noncarriers. This distribution will be done to determine the actual deviation point
			spline_diff[i, j] = weights_point[2,] - weights_point[1,]
		}
	}

	## Extract median and 0.005 and 0.995 values
	#generate blank matrices we will use for lines
	noncarrier_lines = matrix(0, ncol=3, nrow=length(eyo_step))
	carrier_lines = matrix(0, ncol=3, nrow=length(eyo_step))
	diff_lines = matrix(0, ncol=3, nrow=length(eyo_step))
	noncarrier_carrier_t_stats = matrix(0, ncol=1, nrow=length(eyo_step))
	noncarrier_carrier_p_value = matrix(0, ncol=1, nrow=length(eyo_step))
	for (i in 1:ncol(spline_noncarriers)) {
		# For every EYO point (i)
		temp_non = quantile(spline_noncarriers[,i], probs = c(0.005, .5, .995))
		temp_non = unname(temp_non) #get the 0.05, median, and 99.5 percentile values
		temp_car = quantile(spline_carriers[,i], probs = c(0.005, .5, .995))
		temp_car = unname(temp_car) #get the 0.05, median, and 99.5 percentile values
		temp_diff = quantile(spline_diff[,i], probs = c(0.005, .5, .995))
		temp_diff = unname(temp_diff) #get the 0.05, median, and 99.5 percentile values

		noncarrier_lines[i,1] = temp_non[1]
		noncarrier_lines[i,2] = temp_non[2]
		noncarrier_lines[i,3] = temp_non[3] #Write out the 99% Credible intervals and median for each EYO (i)
		carrier_lines[i,1] = temp_car[1]
		carrier_lines[i,2] = temp_car[2]
		carrier_lines[i,3] = temp_car[3] # Carriers
		diff_lines[i,1] = temp_diff[1]
		diff_lines[i,2] = temp_diff[2]
		diff_lines[i,3] = temp_diff[3]

		temp_test_stats = t.test(spline_diff[,i])
		noncarrier_carrier_t_stats[i] = temp_test_stats$statistic[["t"]]
		noncarrier_carrier_p_value[i] = min(sum(spline_diff[,i] < 0)/length(spline_diff[,i]), 1-(sum(spline_diff[,i] < 0)/length(spline_diff[,i])))
	}

#	noncarrier_carrier_t_stats_all_pep <- cbind(noncarrier_carrier_t_stats_all_pep, noncarrier_carrier_t_stats)
#	colnames(noncarrier_carrier_t_stats_all_pep)[ncol(noncarrier_carrier_t_stats_all_pep)] = pep_name

#	noncarrier_carrier_p_value_all_pep <- cbind(noncarrier_carrier_p_value_all_pep, noncarrier_carrier_p_value)
#	colnames(noncarrier_carrier_p_value_all_pep)[ncol(noncarrier_carrier_p_value_all_pep)] = pep_name

	# Add EYO term and column names to matrices for plotting
	flipped_eyo = t(eyo_step)
	flipped_eyo = t(flipped_eyo)
	# Flip everything and make it a matrix
	combined_lines = cbind(flipped_eyo, noncarrier_lines, carrier_lines)
	#combine carriers and noncarriers
	colnames(combined_lines) <- c("eyo","lower_non", "median_non", "upper_non", "lower_car", "median_car", "upper_car")
	diff_lines = cbind(flipped_eyo, diff_lines)
	colnames(diff_lines) <- c("eyo","lower", "median", "upper") #separate matrix for difference, honest this could be combined with the previous one

	combined_lines = as.data.frame(combined_lines) #Make it a data frame or none of the plotting works
	diff_lines = as.data.frame(diff_lines) #Make it a data frame or none of the plotting works

	# manually set difference of medians
	diff_lines$median <- combined_lines$median_car - combined_lines$median_non

	### --------- Plot Figure 1 Starts
	### Plot Images
	temp_data_plot_before <- BL_traits_pep %>%
		select(sample_id, EYO, eval(pep_name), ApoE_Indicator) %>%  #Group_Der_F
		ungroup() %>%
		mutate(Pep = BL_traits_pep[[pep_name]] ) %>%
		dplyr::filter(complete.cases(.))

	before_removal_outlier <- nrow(temp_data_plot_before)

	temp_data_plot <- temp_data_plot_before	%>%
		mutate(zRT = scale(Pep)[,1]) %>%
		filter(between(zRT, -3, +3))
	post_removal_outlier <- nrow(temp_data_plot)
	number_removed <- before_removal_outlier - post_removal_outlier # record number of removed outliers

	CI_Percent <- "99%"

	yaxis_label_plot = name_module_label$y.Axis.Label[which(name_module_label$Raw.File.Name == clean_name)]
	title_label_plot = name_module_label$Title.Label[which(name_module_label$Raw.File.Name == clean_name)]

	ScatterPlot = scatter_plot(x = EYO,
							   y = Pep,
							   df = temp_data_plot,
							   protein_name = data.frame(yaxis_label = yaxis_label_plot,
							   							 title_label = title_label_plot),
							   number_removed = number_removed,
							   CI_Percent = CI_Percent,
							   combined_lines = combined_lines)

#        ggsave(filename = paste0("z:/ShijiaBian/PlasmaProteomic/Result/20250727/simple.3177/scatter/",clean_name,".pdf"), plot = ScatterPlot, device = cairo_pdf,
#               width = 12, height = 12, units = "in", dpi = 300)

	### --------- Plot Figure 1 Ends

	# ### --------- Plot Figure 2 Starts
	min_value = floor(min(diff_lines$lower))
	max_value = ceiling(max(diff_lines$upper))

	# Construct the diff data frame to better plot the significance difference
	diff_lines$same_sign = FALSE
	diff_lines$same_sign[sign(diff_lines$lower) == sign(diff_lines$upper)] = TRUE
	Carr_Non_EYO_diff_df = diff_lines[diff_lines$same_sign == TRUE, ]
	if (nrow(Carr_Non_EYO_diff_df) > 0) {
		Carr_Non_EYO_diff_df$Height = max_value
		Carr_Non_EYO_diff_df$Text_Height = (max_value - min_value) * 0.05
		First_Annotate <- Carr_Non_EYO_diff_df$eyo[1]
		Last_Annotate <- Carr_Non_EYO_diff_df$eyo[length(Carr_Non_EYO_diff_df$eyo)]

		# single color, regardless of # of intervals
		#Annotate_color <- "cornflowerblue"
		#if (sum(Carr_Non_EYO_diff_df$median > 0)/length(Carr_Non_EYO_diff_df$median) > 0.5) { # majortiy of the median > 0
		#	Annotate_color <- "indianred3"
		#}
                # A color for each interval
                Annotate_color <- na.omit(ifelse(Carr_Non_EYO_diff_df$upper < 0, "cornflowerblue", ifelse(Carr_Non_EYO_diff_df$lower > 0, "indianred3", NA)))

		right_limit = 24
		Last_Annotate_Number = Last_Annotate + 0.5
		if (Last_Annotate_Number > right_limit){
			right_limit = 26
		}
	} else {
		First_Annotate = Inf
		Last_Annotate = Inf
		Annotate_color = ""
		right_limit = 24
	}

	DiffPlot = diff_plot(x = eyo,
							y = median,
							df = diff_lines,
							bar_df = Carr_Non_EYO_diff_df,
							protein_name = data.frame(yaxis_label = yaxis_label_plot,
													  title_label = title_label_plot),
							min_value =  min_value, # -31121.83
							max_value = max_value, # 14325
							First_Annotate = First_Annotate,
							Last_Annotate = Last_Annotate,
							Annotate_color = Annotate_color,
							right_limit = right_limit,
							Last_Annotate_Number = Last_Annotate_Number,
							CI_Percent = CI_Percent)

#        ggsave(filename = paste0("z:/ShijiaBian/PlasmaProteomic/Result/20250727/simple.3177/scatter/99par_diff_",clean_name,".pdf"), plot = DiffPlot, device = cairo_pdf,
#               width = 12, height = 12, units = "in", dpi = 300)


	up_down_notation<-rep(NA, EYO_cut) #[count, 'Pep'] = pep_name
	names(up_down_notation)<-as.character(eyo_step)  #as.character(Carr_Non_EYO_diff_df$eyo)
	if (dim(Carr_Non_EYO_diff_df)[1] > 0) {
		for (r in 1:dim(Carr_Non_EYO_diff_df)[1]) {
			if (Carr_Non_EYO_diff_df[r, "lower"] < 0 & Carr_Non_EYO_diff_df[r, "upper"] < 0) {
				up_down_notation[which(names(up_down_notation) == as.character(Carr_Non_EYO_diff_df[r, "eyo"]))] = "cornflowerblue"
			}
			if (Carr_Non_EYO_diff_df[r, "lower"] > 0 & Carr_Non_EYO_diff_df[r, "upper"] > 0) {
				up_down_notation[which(names(up_down_notation) == as.character(Carr_Non_EYO_diff_df[r, "eyo"]))] = "indianred3"
			}

		}
	}
	# ### --------- Plot Figure 2 Ends

	ScatterPlot = scatter_plot(x = EYO,
							   y = Pep,
							   df = temp_data_plot,
							   protein_name = data.frame(yaxis_label = yaxis_label_plot,
							   							 title_label = title_label_plot),
							   number_removed = number_removed,
							   CI_Percent = CI_Percent,
							   combined_lines = combined_lines)
	DiffPlot = diff_plot(x = eyo,
							y = median,
							df = diff_lines,
							bar_df = Carr_Non_EYO_diff_df,
							protein_name = data.frame(yaxis_label = yaxis_label_plot,
							   						  title_label = title_label_plot),
							min_value =  min_value, # -31121.83
							max_value = max_value, # 14325
							First_Annotate = First_Annotate,
							Last_Annotate = Last_Annotate,
							Annotate_color = Annotate_color,
							right_limit = right_limit,
							Last_Annotate_Number = Last_Annotate_Number,
							CI_Percent = CI_Percent)

	gt = arrangeGrob(ScatterPlot, DiffPlot, ncol = 2)
	on_same_plot <- ggplotify::as.ggplot(gt, draw=FALSE)

#        ggsave(filename = paste0("z:/ShijiaBian/PlasmaProteomic/Result/20250727/simple.3177.GSVA/scatter/99par_2on1_",gsub("(\\|)|(\\.)|(\\;)", "_",pep_name),".pdf"), plot = on_same_plot, device = cairo_pdf,
        ggsave(filename = paste0("f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/GSVA/scatter/99par_2on1_",gsub("(\\|)|(\\.)|(\\;)", "_",pep_name),".pdf"), plot = on_same_plot, device = cairo_pdf,
               width = 22.5, height = 10, units = "in", dpi = 300)

  # return a named list with the objects you need afterwards
  list(
    t_stat   = noncarrier_carrier_t_stats,
    p_value  = noncarrier_carrier_p_value,
    up_down  = up_down_notation, #[count, ],
    ok       = TRUE                 # handy flag for foreach result binding
  )
}


library(doParallel)
ncore <- max(1, parallel::detectCores() - 1)
cl     <- makeCluster(ncore)
registerDoParallel(cl)


library(foreach)

# list of packages the workers must load
worker_pkgs <- c("tidyverse","rstanarm","Hmisc","openxlsx","rstan",
                 "gridExtra","ggpubr","ggplot2","grDevices","ggplotify")

results <- foreach(pep_name = pep_names,
                   .packages = worker_pkgs,
                   .export   = c("process_one_peptide",      # the function
                                 "diff_plot","scatter_plot", # custom plotters
                                 "BL_traits_pep","name_match_table",
                                 "name_module_label",        # big objects
                                 "EYO_cut","cairo_pdf"),     # scalars/funs
                   .errorhandling = "pass")  %dopar%  {
  process_one_peptide(pep_name)
}

stopImplicitCluster()   # or stopCluster(cl) if you built one explicitly


# Filter out failed iterations
ok_idx <- vapply(results, function(x) is.list(x) && !inherits(x, "error") && x$ok, TRUE)

noncarrier_carrier_t_stats_all_pep    <- do.call(cbind,
                                                 lapply(results[ok_idx], `[[`, "t_stat"))
noncarrier_carrier_p_value_all_pep    <- do.call(cbind,
                                                 lapply(results[ok_idx], `[[`, "p_value"))
up_down_notation                      <- do.call(rbind,
                                                 lapply(results[ok_idx], `[[`, "up_down")) #[ , -1]

## Set dimnames of final data
colnames(noncarrier_carrier_t_stats_all_pep)<-colnames(noncarrier_carrier_p_value_all_pep)<-rownames(up_down_notation)<-pep_names[ok_idx]

eyo_step = seq(floor(min(BL_traits_pep$EYO, na.rm = T)),
				   ceiling(max(BL_traits_pep$EYO, na.rm = T)),
				   by=.5)
length(eyo_step)==ncol(up_down_notation)  # TRUE
rownames(noncarrier_carrier_t_stats_all_pep)<-rownames(noncarrier_carrier_p_value_all_pep)<-colnames(up_down_notation)<-as.character(eyo_step)


## Get characteristics of waterfall
apply(up_down_notation,2,function(x) table(x))

sum(apply(up_down_notation,1,function(x) length(which(!is.na(x))))>0)
#285  # of 602 ontologies reach significance


# ######### ---- Write Final Outputs for Waterfall(next) step ---- ########
noncarrier_carrier_t_stats_all_pep_final <- noncarrier_carrier_t_stats_all_pep #[, -1]
fil_pep<-"f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/GSVA/scatter/_99_par_diff_all_peptide.rds"
saveRDS(noncarrier_carrier_t_stats_all_pep_final, fil_pep)
fil_csv_pep <- "f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/GSVA/scatter/_99_par_diff_all_peptide.csv"
write.csv(noncarrier_carrier_t_stats_all_pep_final, fil_csv_pep, row.names = TRUE)

noncarrier_carrier_p_value_all_pep_final <- noncarrier_carrier_p_value_all_pep #[, -1]
fil_pep<-"f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/GSVA/scatter/_99_par_diff_all_peptide_p_value.rds"
saveRDS(noncarrier_carrier_p_value_all_pep_final, fil_pep)
fil_csv_pep <- "f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/GSVA/scatter/_99_par_diff_all_peptide_p_value.csv"
write.csv(noncarrier_carrier_p_value_all_pep_final, fil_csv_pep, row.names = TRUE)

fil_csv_up_down_notation<-"f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/GSVA/scatter/_99_par_diff_all_peptide_up_down_notation.csv"
write.csv(up_down_notation, fil_csv_up_down_notation, row.names = TRUE)


#####################
## Visualization - waterfall (R)

## Selected terms 09/06/2025 by ECBJ
selectTerms<-c("GOBP.Steroid Catabolic Process","GOBP.Synapse Assembly","GOBP.Alkanesulfonate Metabolic Process","GOBP.Taurine Metabolic Process","GOBP.Regulation Of Synapse Organization","GOBP.Steroid Metabolic Process","GOBP.Negative Regulation Of Wnt Signaling Pathway","GOBP.Presynaptic Membrane Organization","GOBP.Energy Homeostasis","GOBP.Neuron Projection Regeneration","GOBP.Protein Neddylation","GOBP.Negative Regulation Of Mitophagy","GOBP.Protein Deubiquitination","GOBP.Postsynaptic Membrane Organization","GOBP.Rna 5'-End Processing","GOBP.Hormone Metabolic Process","GOBP.Negative Regulation Of Autophagy","GOBP.Regulation Of Canonical Wnt Signaling Pathway","GOBP.Negative Regulation Of Interferon-Alpha Production","GOBP.Cholesterol Efflux","GOBP.Plasma Lipoprotein Particle Assembly","GOBP.Protein Deneddylation","GOBP.Regulation Of Protein Neddylation","GOBP.Negative Regulation Of Neuron Projection Regeneration","GOBP.Negative Regulation Of Interleukin-12 Production","GOBP.Neuron Projection Development","GOBP.Neuron Development","GOBP.Negative Regulation Of Macroautophagy","GOBP.Positive Regulation Of Protein Monoubiquitination","GOBP.Acute-Phase Response","GOBP.Neuron Differentiation","GOBP.Positive Regulation Of Cholesterol Efflux","GOBP.Negative Regulation Of Response To Wounding","GOBP.Negative Regulation Of Canonical Wnt Signaling Pathway","GOBP.Negative Regulation Of Platelet Activation","GOBP.T Cell Proliferation","GOBP.Regulation Of Leukocyte Migration","GOBP.Negative Regulation Of Coagulation","GOBP.Regulation Of Hemostasis","GOBP.Acute Inflammatory Response","GOBP.Extracellular Matrix Organization","GOBP.Kidney Vasculature Morphogenesis","GOBP.Regulation Of Fibrinolysis","GOMF.Fibroblast Growth Factor Binding","GOMF.Deubiquitinase Activity","GOMF.Ubiquitin-Like Protein Peptidase Activity","GOMF.Heparan Sulfate Proteoglycan Binding","GOMF.Proteoglycan Binding","GOMF.Lipoprotein Particle Receptor Binding","GOMF.Oxidoreductase Activity, Acting On The Ch-Oh Group Of Donors, Nad Or Nadp As Acceptor","GOMF.Structural Constituent Of Eye Lens","GOMF.Low-Density Lipoprotein Particle Receptor Binding","GOMF.Very-Low-Density Lipoprotein Particle Receptor Binding","GOMF.Steroid Dehydrogenase Activity, Acting On The Ch-Oh Group Of Donors, Nad Or Nadp As Acceptor","GOMF.All-Trans-Retinol Dehydrogenase (Nad+) Activity","GOMF.Alcohol Dehydrogenase (Nad+) Activity","GOMF.Aminoacyl-Trna Ligase Activity","GOMF.Opsonin Binding","GOMF.Sterol Transport","GOMF.Lipid Transfer Activity","GOMF.Neurexin Family Protein Binding","GOMF.[Heparan Sulfate]-Glucosamine 3-Sulfotransferase Activity","GOMF.Heparan Sulfate Sulfotransferase Activity","GOMF.Lipid Transport","GOMF.Fibronectin Binding","GOMF.Extracellular Matrix Binding","GOMF.Insulin-Like Growth Factor I Binding","GOMF.Complement Binding","GOMF.Insulin-Like Growth Factor Binding","GOMF.Complement Component C3b Binding","GOMF.Glycosaminoglycan Binding","GOMF.Heparin Binding","GOMF.Integrin Binding","GOMF.Insulin-Like Growth Factor Ii Binding","GOMF.Amylase Activity","GOMF.Extracellular Matrix Structural Constituent","GOMF.Hyaluronic Acid Binding","GOMF.Vascular Endothelial Growth Factor Receptor Activity","GOCC.Endocytic Vesicle Lumen","GOCC.Intermediate-Density Lipoprotein Particle","GOCC.Spliceosomal Tri-Snrnp Complex","GOCC.U4/U6 X U5 Tri-Snrnp Complex","GOCC.Precatalytic Spliceosome","GOCC.U2-Type Precatalytic Spliceosome","GOCC.U4 Snrnp","GOCC.Neurofilament","GOCC.U1 Snrnp","GOCC.Mhc Class Ii Protein Complex","GOCC.Sno(S)Rna-Containing Ribonucleoprotein Complex","GOCC.Ubiquitin Ligase Complex","GOCC.Spliceosomal Snrnp Complex","GOCC.Death-Inducing Signaling Complex","GOCC.Low-Density Lipoprotein Particle","GOCC.High-Density Lipoprotein Particle","GOCC.Intermediate Filament Cytoskeleton","GOCC.Triglyceride-Rich Plasma Lipoprotein Particle","GOCC.Very-Low-Density Lipoprotein Particle","GOCC.Chylomicron","GOCC.Lysosomal Lumen","GOCC.Synaptic Cleft","GOCC.Perikaryon","GOCC.Synaptobrevin 2-Snap-25-Syntaxin-1a-Complexin Ii Complex","GOCC.Extracellular Matrix",
               "GOCC.Blood Microparticle","GOCC.Axon Initial Segment","GOCC.Membrane Attack Complex","GOCC.Vacuolar Lumen","GOCC.Azurophil Granule Lumen","GOCC.Perisynaptic Extracellular Matrix","GOCC.Microfibril","GOCC.Schaffer Collateral - Ca1 Synapse","GOCC.Late Endosome","GOCC.Lysosome","GOCC.Insulin-Like Growth Factor Binding Protein Complex","REACTOME.The Canonical Retinoid Cycle In Rods (Twilight Vision)","REACTOME.Incretin Synthesis, Secretion, And Inactivation","REACTOME.Akt-Mediated Inactivation Of Foxo1a","REACTOME.Regulation Of Tnfr1 Signaling","REACTOME.Synthesis, Secretion, And Inactivation Of Glucose-Dependent Insulinotropic Polypeptide (Gip)","REACTOME.Synthesis, Secretion, And Inactivation Of Glucagon-Like Peptide-1 (Glp-1)","REACTOME.Deubiquitination","REACTOME.Ligand-Dependent Caspase Activation","REACTOME.Wnt Ligand Biogenesis And Trafficking","REACTOME.Ra Biosynthesis Pathway","REACTOME.Cytosolic Trna Aminoacylation","REACTOME.Ldl Remodeling","REACTOME.Hs-Gag Biosynthesis","REACTOME.Activation Of Matrix Metalloproteinases","REACTOME.Heparan Sulfate Heparin (Hs-Gag) Metabolism","REACTOME.Interleukin-4 And Interleukin-13 Signaling","REACTOME.Complement Cascade","REACTOME.Regulation Of Complement Cascade","REACTOME.Ncam1 Interactions","REACTOME.Neurofascin Interactions","REACTOME.Notch-Hlh Transcription Pathway","REACTOME.Interleukin-12 Signaling","REACTOME.L1cam Interactions","REACTOME.Nrcam Interactions","REACTOME.Terminal Pathway Of Complement","REACTOME.Glycosaminoglycan Metabolism","REACTOME.Extracellular Matrix Organization","REACTOME.Intrinsic Pathway Of Fibrin Clot Formation","REACTOME.Sema3a-Plexin Repulsion Signaling By Inhibiting Integrin Adhesion","REACTOME.Neurophilin Interactions With Vegf And Vegfr")
#note: selectedTerms above is length=29
length(selectTerms)
#144

gene_sets_present<- lapply(gene_sets_filtered, function(set) { set[which(set %in% symbols.collapsed)] })
## Can a kappa score for similarity of ontologies based on the gene membership (gene symbols) only be calculated similar to what is done to cluster ontology terms in Metascape? If so, implement the calculation in R given a list of gene ontologies with vectors of gene symbols stored in the list variable gene_sets_present. Cluster the ontologies with a branch cut threshold of kappa=0.30, returning the ontology names in a column of a data frame, with a second column containing the cluster number for each.

# gene_sets_present: named list of character vectors (genes per ontology)
# returns: data.frame(term, cluster) using hierarchical clustering at kappa >= 0.30
kappa_cluster_ontologies <- function(gene_sets_present, kappa_cut = 0.30,
                                     method = c("hier", "graph")) {
  method <- match.arg(method)

  # --- 1) Build universe and index encoding ---
  terms <- names(gene_sets_present)
  if (is.null(terms)) stop("gene_sets_present must be a *named* list.")
  gene_sets_present <- lapply(gene_sets_present, unique)
  U <- sort(unique(unlist(gene_sets_present)))

# ------------------------------------------------------------------------
# ANNOTATION: Write STAN-derived ontology p-value/effect-size matrices and
# direction calls.
# ------------------------------------------------------------------------
  N <- length(U)
  idx_map <- setNames(seq_along(U), U)
  sets_idx <- lapply(gene_sets_present, function(gs) sort(idx_map[gs]))
  n <- length(sets_idx)
  sizes <- vapply(sets_idx, length, integer(1))

  # --- 2) Pairwise intersections (upper triangle) ---
  # Efficiently compute |Gi (is member of set) Gj| without building a huge incidence matrix
  intersec_mat <- matrix(0L, n, n, dimnames = list(terms, terms))
  for (i in seq_len(n-1)) {
    gi <- sets_idx[[i]]
    for (j in (i+1):n) {
      intersec_mat[i, j] <- length(intersect(gi, sets_idx[[j]]))
    }
  }
  intersec_mat <- intersec_mat + t(intersec_mat) + diag(sizes) # fill symmetric + diagonal (|Gi|)

  # --- 3) Cohen’s kappa from 2x2 counts over the universe ---
  # For i,j:
  # a=|Gi?Gj|; b=|Gi\Gj|=|Gi|-a; c=|Gj\Gi|=|Gj|-a; d=N - (a+b+c)
  # p0=(a+d)/N; pe=[(|Gi|/N * |Gj|/N) + ((N-|Gi|)/N * (N-|Gj|)/N)]
  # kappa=(p0-pe)/(1-pe)
  A <- intersec_mat
  Bi <- matrix(sizes, n, n, byrow = FALSE) - A
  Cj <- matrix(sizes, n, n, byrow = TRUE)  - A
  D  <- N - (A + Bi + Cj)

  p0 <- (A + D) / N
  pe <- (outer(sizes, sizes, "*") + outer(N - sizes, N - sizes, "*")) / (N^2)

  K  <- (p0 - pe) / (1 - pe)
  diag(K) <- 1
  K[is.nan(K)] <- 0  # guard rare degenerate cases

  if (method == "hier") {
    # --- 4A) Hierarchical clustering with branch cut at kappa_cut ---
    # distance = 1 - kappa, cut height = 1 - kappa_cut
    Dmat <- as.dist(pmax(0, 1 - K))  # ensure non-negative
    hc <- hclust(Dmat, method = "average")
    cl <- cutree(hc, h = 1 - kappa_cut)
    out <- data.frame(term = terms, cluster = cl[terms], row.names = NULL)
    return(out[order(out$cluster), ])
  } else {
    # --- 4B) Graph clustering via connected components at kappa threshold ---
    # Edge if kappa >= kappa_cut
    adj <- (K >= kappa_cut) & !diag(n)
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("Install 'igraph' for graph-based clustering or use method='hier'.")
    }
    g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
    comp <- igraph::components(g)$membership
    out <- data.frame(term = terms, cluster = comp, row.names = NULL)
    return(out[order(out$cluster), ])
  }
}

## Example usage:
# Use method="hier" if you specifically want a “branch cut” (cut height = 0.70).
# Use method="graph" for Metascape-style kappa-network clustering (connected components at kappa >= 0.30).
# result_hier  <- kappa_cluster_ontologies(gene_sets_present, kappa_cut = 0.30, method = "hier")
# result_graph <- kappa_cluster_ontologies(gene_sets_present, kappa_cut = 0.30, method = "graph")

## If we select by kappa:
#ontology_clusts<-kappa_cluster_ontologies(gene_sets_present, kappa_cut = 0.30, method = "graph")
## Here, select all the 144 terms manually curated by ECBJ
selectTerms <- sub("\\.", "_", selectTerms)  #replace first "." only
ontology_clusts <- data.frame(
  term    = selectTerms,
  cluster = 1:length(selectTerms),
  stringsAsFactors = FALSE
)
# All all selected terms in our data?
which(!selectTerms %in% names(gene_sets_present))
#integer(0)


## Using ggplot2, given inputs of (1) noncarrier_carrier_p_value_all_pep_final, a data frame carrying 599 named columns for 599 ontologies significance (p value) of difference in e4/4 homozygotes vs. e3/3 homozygotes at each of 143 1/2 year intervals from (rownames) "-46" to "25"; (2) up_down_notation, a complimentary data frame of transposed dimensions 599 rows x 143 columns, indicating for each position whether the difference, when significant, is lower in e4/4 ("cornflowerblue"), or higher ("indianred3"). When not significant, values in this data frame are NA. And (3) ontology_clusts, a data frame with column "term" carrying the names of available ontologies, and a second column named "cluster", with a number assigned to all ontologies so that similar ontologies have the same cluster number, produce the following graphical heatmap representation of significant ontology intervals over the time period -46 to +25 years, on the x-axis, labeled as "EYO".
## The heatmap will represent at most one ontology from each of the clusters given in ontology_clusts, provided that there is at least 1 half-year interval reaching p=<0.005. P values less than 1/4000 should be set to 1/8000. A white-to-blue palette for the heatmap will be used in positions along EYO (x) when the up_down_notation for that position is indicated as "cornflowerblue", and white-to-darkred palette for the heatmap track will be used in positions along EYO (x) when up_down_notation for that position indicates "indianred3".  The heatmap tracks along y, from top to bottom, should be sorted in order of earlier EYO intervals reaching significance first; when there are ties for significance at the same EYO interval, the track with higher significance (lower p) in that interval should be plotted first.

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)

##--- 0) Inputs expected ---------------------------------------------
## noncarrier_carrier_p_value_all_pep_final : 143 x 599 data.frame (rows = EYO -46..25 as rownames, cols = ontology terms)
## up_down_notation                         : 599 x 143 data.frame (rows = ontology terms as rownames, cols = EYO -46..25)
## ontology_clusts                          : data.frame(term, cluster)

##--- 1) Long-format p-values (143*599 rows) -------------------------
pvals_long <- as.data.frame(noncarrier_carrier_p_value_all_pep_final) %>%
  mutate(EYO = as.numeric(rownames(.))) %>%
  relocate(EYO) %>%
  pivot_longer(-EYO, names_to = "term", values_to = "p")

## Cap extreme p’s: p < 1/4000 ? 1/8000
pvals_long <- pvals_long %>%
  mutate(p_cap = if_else(p < 1/4000, 1/8000, p))

##--- 2) Long-format “direction” (cornflowerblue / indianred3 / NA) --
dir_long <- as.data.frame(up_down_notation) %>%
  mutate(term = rownames(.)) %>%
  relocate(term) %>%
  pivot_longer(-term, names_to = "EYO", values_to = "direction") %>%
  mutate(EYO = as.numeric(EYO))

##--- 3) Merge with clusters & compute selection features ------------
dat <- pvals_long %>%
  left_join(dir_long, by = c("term","EYO")) %>%
  left_join(ontology_clusts, by = "term")

# Keep only positions that are “eligible significant” (p ? 0.005 AND have a direction)
sig_dat <- dat %>%
  filter(!is.na(direction), p <= 0.005)

# Per term, find earliest EYO reaching significance and the best p at that earliest EYO
term_sig_summary <- sig_dat %>%
  group_by(term, cluster) %>%
  dplyr::summarize(earliest_sig_EYO = min(EYO, na.rm = TRUE), .groups = "drop") %>%
  left_join(sig_dat, by = c("term","cluster")) %>%
  filter(EYO == earliest_sig_EYO) %>%
  group_by(term, cluster, earliest_sig_EYO) %>%
  dplyr::summarize(min_p_at_earliest=min(p, na.rm=TRUE), .groups="drop")

# Choose <=1 term per cluster:
#   priority = earliest_sig_EYO (ascending), then min_p_at_earliest (ascending)
picked_terms <- term_sig_summary %>%
  arrange(earliest_sig_EYO, min_p_at_earliest) %>%
  distinct(cluster, .keep_all = TRUE) %>%
  pull(term)

##--- 4) Prepare heatmap table for the picked terms -------------------
plot_tab <- dat %>%
  filter(term %in% picked_terms) %>%
  # signed score: -log10(p_cap) with sign from direction
  mutate(
    signed_score = case_when(
      direction == "indianred3"     ~  +(-log2(p_cap)),  # higher in e4/4 -> +red
      direction == "cornflowerblue" ~  -(-log2(p_cap)),  # lower  in e4/4 -> -blue
      TRUE ~ NA_real_
    )
  )

# Order Y tracks: earlier EYO first; ties by stronger (lower) p at that earliest EYO
term_order <- term_sig_summary %>%
  arrange(earliest_sig_EYO, min_p_at_earliest) %>%
  pull(term)

plot_tab <- plot_tab %>%
  mutate(term = factor(term, levels = rev(term_order)))

##--- 5) Draw heatmap -------------------------------------------------
# Diverging scale makes white near 0, blue for negatives, red for positives
waterfall.all<-ggplot(plot_tab %>% filter(!is.na(signed_score))) +
  geom_tile(aes(x = EYO, y = term, fill = signed_score), width = 0.9, height = 0.9) +
#  scale_fill_gradient2(
#    low = "darkslateblue", mid = "white", high = "darkred", midpoint = 0,
#    name = expression(paste("signed  ", -log[2], "(p)")),
##    limits = c(-max(-log2(plot_tab$p_cap), na.rm = TRUE),
##                max(-log2(plot_tab$p_cap), na.rm = TRUE))
#  limits = c(-log2(1/8192), log2(1/8192)),
#  breaks = c(-log2(0.005), 0, log2(0.005)),
#  labels = c(
#    expression(-log[2](0.005)), "0", expression(log[2](0.005)) )
#  ) +
scale_fill_gradientn(
  colors = c("darkslateblue", "white", "white", "white", "darkred"),
  values = scales::rescale(c(-13, -7.4, 0, 7.4, 13)),
  name = expression(paste("signed  ", -log[2], "(p)")),
  limits = c(-13, 13)
) +
guides(
  fill = guide_colorbar(
    barheight = unit(20, "cm"),  # Adjust height here
    barwidth = unit(0.5, "cm")   # Optional: make it narrow
  )
) +
# Thin vertical lines every 5 units from -45 to +25
geom_vline(
  xintercept = seq(-45, 25, by = 5),
  color = "gray70",
  size = 0.3
) +
# Thicker black vertical line at x = 0
geom_vline(
  xintercept = 0,
  color = "black",
  size = 0.8
) +
  labs(x = "EYO", y = NULL,
       title = "Significant ontology intervals (e4/4 vs e3/3)",
       subtitle = "One representative ontology per cluster; white->blue = lower in e4/4, white->red = higher") +
  theme_minimal(base_size = 12.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 11),
    legend.position = "right"
  )


#ggsave(filename = paste0("z:/ShijiaBian/PlasmaProteomic/Result/20250727/simple.3177.GSVA_waterfall-ALL.pdf"), plot = waterfall.all, device = cairo_pdf,
ggsave(filename = paste0("f:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/3.GSVA-waterfall-simple.3177-ALL-selected_144_terms.pdf"), plot = waterfall.all, device = cairo_pdf,
 width = 19, height = 35, units = "in", dpi = 300)
