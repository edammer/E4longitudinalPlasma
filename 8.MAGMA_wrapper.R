##############################################################################
# Pipeline annotation header: 8.MAGMA_wrapper.R
# Manuscript code section(s): 8
#
# Purpose:
# Configure and run the list-vector adaptation of MAGMA.SPA to test
# semaglutide-significant assay enrichment in user-supplied protein/gene
# modules.
#
# Principal inputs:
#   - Sema_S6_NominalP.csv
#   - Sema_S6_Qvalue.csv
#   - moduleGeneList(881+322).csv
#   - moduleGeneList(332+804+881).csv
#   - MAGMA.SPA_listVectorInput.R
#
# Principal outputs:
#   - MAGMA.SPA XLSX/PDF/permutation-statistics outputs with Sema prefix
#
# Step overview:
#   1. Set MAGMA.SPA input paths, p-value/FDR thresholds, colors, and parallel
#      thread count.
#   2. Read moduleGeneList CSVs as ordered list vectors instead of deriving
#      modules from cleanDat/net$colors.
#   3. Source MAGMA.SPA_listVectorInput.R and run MAGMA.SPA() for nominal and
#      q-value semaglutide hit lists.
#   4. Repeat enrichment over alternative module-list definitions spanning
#      selected EYO/age protein signatures.
#
# Notes added during manuscript-code alignment:
#   - Original executable code was preserved except for whitespace/encoding cleanup
#     and the explicitly marked non-run archive described below when applicable.
#   - Comments inserted by this pass are prefixed with "ANNOTATION:".
##############################################################################

# ------------------------------------------------------------------------
# ANNOTATION: Set MAGMA.SPA wrapper-level inputs and user-adjustable
# analysis parameters.
# ------------------------------------------------------------------------
#MAGMA-SPA (Seyfried Pipeline Adaptation for MAGMA)
#---------------------------------

# Required parameters, variables, and data must be set as shown above in .GlobalEnv before calling function; currently no defaults are automatic.
##################################
MAGMAinputDir= "F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/#manuscript/SciTranslMed_plan/FET_Sema_881_322/"

MAGMAinputs= c(	#"AD_GWAS_ENSEMBLE_averageMinusLogP(Plt0.05).csv",
                "Sema_S6_NominalP.csv",
                "Sema_S6_Qvalue.csv")     #These files must be in MAGMAinputDir

maxP=0.05                 #no genes with a MAGMA summarized p value greater than this will be considered even if in the MAGMA-derived input files.
FDR=0.10                  #FDR or q value (0 < FDR < 1); recommend 0.10, i.e. 10%
barcolors= c("darkslateblue","mediumorchid") #,"hotpink")  #specify one unique color for each of above MAGMAinputs
                          #common colors: "darkslateblue","mediumorchid", "seagreen3","hotpink","goldenrod","darkorange","darkmagenta", ...
relatednessOrderBar=FALSE  #Plot mean scaled enrichment bar plot in column order (relatedness) of MEs?  If FALSE, they will be plotted in size rank order M1, M2, ...


# ------------------------------------------------------------------------
# ANNOTATION: Document the objects that the original Seyfried Analysis
# Pipeline version would have required.
# ------------------------------------------------------------------------
# Data created during the Seyfried Analysis Pipeline
##################################
#NETcolors= net$colors     #module color assignments, vector of length equal to number of rows in cleanDat; should have all colors for modules from 1:minimumSizeRank as printed by WGCNA::labels2colors(1:nModules)
#MEs= MEs                 #Module Eigengenes (or Eigenproteins) with columns of MEs ordered in relatedness order
#cleanDat= cleanDat       #rownames must start with HUMAN gene symbols, separated by any other rowname information using ';' or '|' character


# ------------------------------------------------------------------------
# ANNOTATION: Set output naming, plotting, permutation, and parallelization
# controls.
# ------------------------------------------------------------------------
# Other variables
#################################
outFilePrefix="Sema"         #Filename prefix; step in the pipeline -- for file sorting by name.
outFileSuffix="SignificanceEnrichment.in.EYO881_or322before50_hits"
parallelThreads=8         #Each permutation analysis is run on a separate thread simultaneously, up to this many threads.
calculateMEs=FALSE         #Recalculate MEs and their relatedness order, even if the data already exists.
plotOnly=FALSE            #If plotOnly is TRUE, the variables created by MAGMA.SPA function holding plot data should already exist (xlabels, allBarData).
##################################

setwd(MAGMAinputDir)


# ------------------------------------------------------------------------
# ANNOTATION: Read the first moduleGeneList CSV and run the adapted
# MAGMA.SPA function.
# ------------------------------------------------------------------------
moduleGeneList<-as.list(read.csv(file="moduleGeneList(881+322).csv", header=TRUE,row.names=NULL))
moduleGeneList<-lapply(moduleGeneList,na.omit)


# Run the permutation analysis and generate all outputs
source("MAGMA.SPA_listVectorInput.R")
MAGMAoutList <- MAGMA.SPA()
# Outputs XLSX, PDF, and list of barplot y values (allBarData), barplot labels (xlabels), and all permutation statistics and gene symbol hits (all_output)



# ------------------------------------------------------------------------
# ANNOTATION: Repeat the MAGMA.SPA run with an alternative moduleGeneList
# definition.
# ------------------------------------------------------------------------
## Repeat with p value lists only, <50, >50, all
source("MAGMA.SPA_listVectorInput.R")
outFileSuffix="SignificanceEnrichment.in.EYO881_322before50_or_804after50-hits"
MAGMAinputs= c( "Sema_S6_NominalP.csv") #, "Sema_S6_Qvalue.csv")
barcolors= c("darkslateblue") #,"mediumorchid")

moduleGeneList<-as.list(read.csv(file="moduleGeneList(332+804+881).csv", header=TRUE, row.names=NULL))
moduleGeneList<-lapply(moduleGeneList,na.omit)

MAGMAoutList.3wayPonly <- MAGMA.SPA()
