##############################################################################
# Pipeline annotation header: 1c.4p13.2xRegr_1stRegr.2protPAVwithinSite_2ndReg(Site-3ways_orMedianZero)+QC_VP+tSNE+Volc_b345regrFIXED.R
# Manuscript code section(s): 1 / 2 / 3 QC and regression development
#
# Purpose:
# Generate plasma-only QC, pre-analytical-variance adjustment, site/batch
# regression candidates, site-F sub-batch labels, APOE imputation
# development runs, and final regression/QC products.
#
# Principal inputs:
#   - 2.saved.image_trait+human_cleanup_nm+em0_V1_3ms_03-27-25.RData
#   - 4p13.cleanDat.22sites.RDS
#   - 4p13.numericMeta.22sites.RDS
#   - SomaSignals cohort map CSV files
#   - parANOVA.dex.fallback7.25.R and related plotting helpers
#
# Principal outputs:
#   - Multiple 4p13a/b/c QC PDFs
#   - normExpr.reg_sites1-19*.RDS
#   - variancePartition RDS/PDF outputs
#   - APOE prediction RDS files
#   - saved.image-genotype_prediction_finalized.RData
#
# Step overview:
#   1. Restrict to SomaScan plasma samples, remove calibrators, and add pre-
#      analytical proxy assays (HNRNPA2B1 and HBZ) to metadata.
#   2. Create t-SNE QC plots across sites and traits, excluding low-coverage
#      sites and segmenting site F into F1/F2/F3 sub-batches.
#   3. Regress the two PAV proxy assays within site and compare
#      residualization options by variance partitioning and t-SNE.
#   4. Fit candidate site regressions that protect age, sex, and APOE-related
#      covariates for biological downstream analysis.
#   5. Develop and assess APOE genotype imputation candidates and finalize the
#      protected APOE e4 carrier status used in the second-pass regression.
#   6. Run final QC volcano plots and save finalized objects for stage 2 APOE
#      prediction and WGCNA.
#
# Notes added during manuscript-code alignment:
#   - Original executable code was preserved except for whitespace/encoding cleanup
#     and the explicitly marked non-run archive described below when applicable.
#   - Removed 1801 lines after the explicit "NOT RUN BELOW HERE" marker from
#   - this annotated copy; the original upload remains unchanged.
#   - Comments inserted by this pass are prefixed with "ANNOTATION:".
##############################################################################

# ------------------------------------------------------------------------
# ANNOTATION: Load curated traits/protein data and initialize the plasma
# harmonization workspace.
# ------------------------------------------------------------------------
##################################################
# 4. Plasma network exploration

rootdir="z:/EBD/"
#rootdir="c:/Users/edammer_4ceb6ff/Downloads/"
setwd(rootdir)
#load("4p6.unreg.saved.image.PLASMA_18sitesUnregressed(7335x22547)_WGCNA.mms10.ds4.pwr11.RData")
#load("4p9.2xRegr.saved.image.PLASMA_18siteRegr(7333x22547)_WGCNA.mms10.ds4.pwr7.RData")

load("2.saved.image_trait+human_cleanup_nm+em0_V1_3ms_03-27-25.RData")
   #("3.saved.image_trait+human_cleanup2_2fluidSplit+person_ids_CSFandPlasma_(both).RData")

library(purrr)
cleanDat<-exprMat0[,match(rownames(numericMeta.blood),colnames(exprMat0))]
numericMeta<-numericMeta.blood

table(numericMeta$sample_type)
#Calibrator     Sample
#      1115      27850

numericMeta$is_somalogic[which(numericMeta$sample_type=="Calibrator")]<-1


#NETWORK NOT RUN YET -- NOTE- sample_type=="Calibrator" samples are here - but we are not going to do TAMPOR (mode 1), or any TAMPOR in this cleanup iteration
#numericMeta$Batch<-numericMeta$contributor_code
dim(cleanDat)
#7335 28965   # previously 27878 w/o calibrators but some extra samples compared to here.
which(is.na(rownames(cleanDat)))
# integer(0)
#cleanDat<-cleanDat[which(!is.na(rownames(cleanDat))),]
#dim(cleanDat)
#7334 27878  previously  #CSF:3233 samples/columns


## First subset of data appropriate for TAMPOR: 3 sites with full 7k (7334) assays
table(apply(cleanDat,2,function(x) length(which(!is.na(x)))))
# 1267  3709  5031  7333  7334   - with this many values not NA:
#  135  1061  2254   406 25109   - count of samples  now

# ------------------------------------------------------------------------
# ANNOTATION: Remove non-SomaScan and calibrator rows, leaving plasma
# samples for batch/QC analyses.
# ------------------------------------------------------------------------
#   95  1040  2210   406 24127   - count of samples  previously

table(numericMeta$contributor_code[which(apply(cleanDat,2,function(x) length(which(!is.na(x))))==7335)])
# here (PLASMA)
#   A    B    C    D    E    F    G    I    J    L    M    N    P    Q    R    S    T    U
#1058 1303 2120  788  678 4076 1330 1411  814 1191  991  839 1596 1455  710  215  322 4212  now
# 983 1228 2000  743  648 3966 1330 1333  827 1101  931  789 1491 1370  670  200  302 4215  previously

# previously (CSF) before calibrants
#   N    Q    T
# 307 1370  278

table(numericMeta$contributor_code)
# here (PLASMA)
#   A    B    C    D    E    F    G    H    I    J    K    L    M    N    P    Q    R    S    T    U    V    W
#1058 1303 2120  788  678 4076 1330 2254 1411  814  406 1191  991  839 1596 1455  710  215  322 4212 1061  135  now
# 983 1228 2000  743  648 3966 1330 2210 1333  827  406 1101  931  789 1491 1370  670  200  302 4215 1040   95  previously without calibrators but some extra samples compared to now

# previously (CSF) before calibrants added
#   J    N    O    Q    T
# 122  307 1156 1370  278


# Remove ms samples (!is_somalogic==1)
length(which(numericMeta$is_somalogic==1))
#28965   # previously: 27850 without Calibrators
numericMeta<-numericMeta[which(numericMeta$is_somalogic==1),]
cleanDat<-cleanDat[,match(rownames(numericMeta),colnames(cleanDat))]


table(numericMeta$contributor_code) #--without MS
# here (PLASMA)
#   A    B    C    D    E    F    G    H    I    J    K    L    M    N    P    Q    R    S    T    U    V    W
#1058 1303 2120  788  678 4076 1330 2254 1411  814  406 1191  991  839 1596 1455  710  215  322 4212 1061  135  now - with calibrator samples
# 983 1228 2000  743  638 3966 1330 2254 1331  814  406 1101  931  789 1491 1370  670  200  302 4212  996   95  now - without calibrator samples (missed is_somalogic==1 setting for Calibrator samples)
# 983 1228 2000  743  648 3966 1330 2210 1333  827  406 1101  931  789 1491 1370  670  200  302 4215 1040   95  previously without calibrators but some extra samples compared to now


table(numericMeta$sample_type)
#Calibrator     Sample
#      1115      27850

# Remove Calibrators
numericMeta<-numericMeta[which(!numericMeta$sample_type=="Calibrator"),]

# ------------------------------------------------------------------------
# ANNOTATION: Add PAV proxy assays and derived binary clinical/amyloid
# variables to the metadata.
# ------------------------------------------------------------------------
cleanDat<-cleanDat[,match(rownames(numericMeta),colnames(cleanDat))]


table(numericMeta$contributor_code) #--without MS, and without SOMAscan calibrators
#   A    B    C    D    E    F    G    H    I    J    K    L    M    N    P    Q    R    S    T    U    V    W
# 983 1228 2000  743  638 3966 1330 2254 1331  814  406 1101  931  789 1491 1370  670  200  302 4212  996   95  now - without calibrator samples


## TRAIT FINALIZATION / CLEANUP (EARLIER, BEFORE REGRESSION OF BATCH?SITE
# Preanalytical factor proxies - as traits before their regression
HBZ.Soma=cleanDat["HBZ|P02008",]
HNRNPA2B1.Soma=cleanDat["HNRNPA2B1|P22626",]

length(which(is.na(HBZ.Soma)))
#0
length(which(is.na(HNRNPA2B1.Soma)))
#0

numericMeta$RegrBloodPreanalyticFactor.HNRNPA2B1<-HNRNPA2B1.Soma[match(colnames(cleanDat),rownames(numericMeta))]
numericMeta$RegrBloodPreanalyticFactor.HBZ<-HBZ.Soma[match(colnames(cleanDat),rownames(numericMeta))]


# Amyloid Positivity
numericMeta$AmyloidPositivity.01<-NA
numericMeta$AmyloidPositivity.01[which(numericMeta$AmyloidPositivity.withRM=="NEGATIVE")]<-0
numericMeta$AmyloidPositivity.01[which(numericMeta$AmyloidPositivity.withRM=="POSITIVE")]<-1


## MMSE values above 30 are only remaining trait not ready for Global Network Plots
length(which(numericMeta$MMSE>30))
# 24

numericMeta$MMSE[which(numericMeta$MMSE>30)]<-NA


##########################################################################################
## SomaSignals in BH

#BioHermes
SomaSignal.traits<-read.csv(file="BH-SomaSIgnals_forR_andHDS_BH.map.csv",header=TRUE,row.names=1,check.names=TRUE)
SomaSignal.traits$ShortID<-gsub("-",".",SomaSignal.traits$ShortID)

# ------------------------------------------------------------------------
# ANNOTATION: Read external SomaSignal maps and merge time-to-processing
# variables from mapped BioHermes, ROSMAP, and UDS samples.
# ------------------------------------------------------------------------
SomaSignal.traits$HDSid<-BH.map$SiteA.sample[match(SomaSignal.traits$ShortID,BH.map$BH.sample)]

SomaSignal.traits<-SomaSignal.traits[match(colnames(cleanDat),SomaSignal.traits$HDSid),]

length(which(!is.na(SomaSignal.traits[,2])))
#975 mapped B-H samples


#ROSMAP
SomaSignal.traits.RM<-read.csv(file="ROSMAP-SomaSignalsForR_andHDS_RM.map.csv",header=TRUE,row.names=1,check.names=TRUE)
SomaSignal.traits.RM$ShortID<-gsub("-",".",SomaSignal.traits.RM$ShortID)
SomaSignal.traits.RM$HDSid<-RM.map$SiteR.sample[match(SomaSignal.traits.RM$ShortID,RM.map$RM.sample)]

SomaSignal.traits.RM<-SomaSignal.traits.RM[match(colnames(cleanDat),SomaSignal.traits.RM$HDSid),]

length(which(!is.na(SomaSignal.traits.RM[,2])))
#670 mapped RM samples


#UDS
SomaSignal.traits.UDS<-read.csv(file="UDS-SomaSignalsForR_andHDS_UDS.map.csv",header=TRUE,row.names=1,check.names=TRUE)
SomaSignal.traits.UDS$ShortID<-gsub("-",".",SomaSignal.traits.UDS$ShortID)
SomaSignal.traits.UDS$HDSid<-UDS.map$SiteD.sample[match(SomaSignal.traits.UDS$ShortID,UDS.map$UDS.sample)]

SomaSignal.traits.UDS<-SomaSignal.traits.UDS[match(colnames(cleanDat),SomaSignal.traits.UDS$HDSid),]

length(which(!is.na(SomaSignal.traits.UDS[,2])))
#481 mapped UDS samples


## Merged, in order for new columns of numericMeta  ;  taking the first non-NA value at each position
SomaSignals.for.numericMeta<- data.frame(TimeToSpin=dplyr::coalesce(SomaSignal.traits$TimeToSpinPlasma_R1143.V01, SomaSignal.traits.RM$TimeToSpinPlasma_R1143.V01, SomaSignal.traits.UDS$TimeToSpinPlasma_R1143.V01),
                                         TimeToDecant=dplyr::coalesce(SomaSignal.traits$TimeToDecantPlasma_R1144.V01, SomaSignal.traits.RM$TimeToDecantPlasma_R1144.V01, SomaSignal.traits.UDS$TimeToDecantPlasma_R1144.V01),
                                         TimeToFreeze=dplyr::coalesce(SomaSignal.traits$TimeToFreezePlasma_R1145.V01, SomaSignal.traits.RM$TimeToFreezePlasma_R1145.V01, SomaSignal.traits.UDS$TimeToFreezePlasma_R1145.V01),
                                         FedFastedTime=dplyr::coalesce(SomaSignal.traits$FedFastedPlasma_R1152.V01, SomaSignal.traits.RM$FedFastedPlasma_R1152.V01, SomaSignal.traits.UDS$FedFastedPlasma_R1152.V01),
                                         FreezeThawCycles=dplyr::coalesce(SomaSignal.traits$FreezeThawPlasma_R1154.V01, SomaSignal.traits.RM$FreezeThawPlasma_R1154.V01, SomaSignal.traits.UDS$FreezeThawPlasma_R1154.V01))
length(which(!is.na(SomaSignals.for.numericMeta$TimeToSpin)))
#2126 mapped values in 3 cohorts A, D, and R


numericMeta$TimeToSpin<-SomaSignals.for.numericMeta$TimeToSpin
numericMeta$TimeToDecant<-SomaSignals.for.numericMeta$TimeToDecant
numericMeta$TimeToFreeze<-SomaSignals.for.numericMeta$TimeToFreeze
numericMeta$FedFastedTime<-SomaSignals.for.numericMeta$FedFastedTime
numericMeta$FreezeThawCycles<-SomaSignals.for.numericMeta$FreezeThawCycles

##########################################################################################

numericMeta.22sites<-numericMeta
#Backup traits - 96 columns, 27850 samples (no calibrators, somascan only):  numericMeta.22sites<-readRDS("4p13.numericMeta.22sites.RDS")
#Backup log2(RFU) matrix:  cleanDat.22sites<-readRDS("4p13.cleanDat.22sites.RDS")


##################################################
# 4p13a. Initial QC (tSNE plots)

#rootdir="z:/EBD/"
##rootdir="c:/Users/edammer_4ceb6ff/Downloads/"
#setwd(rootdir)

####################
## Examine tSNE of PLASMA - Human only assays (rows) and samples (columns);  no missing data by row.

#Group.3mappedCohortsPlusCTimputed)
Group=as.factor(numericMeta$Group.withCTimputed)


# ------------------------------------------------------------------------
# ANNOTATION: Generate starting t-SNE maps of plasma samples colored by
# site, age, sample type, and APOE e4 dose.
# ------------------------------------------------------------------------
numericMeta.plasma<-numericMeta #[which(numericMeta$sample_matrix=="CSF"),]
dim(numericMeta.plasma)
#27850    96   #previously (16 sites with Calibrators):  18739   93
exprMat.plasma<-cleanDat

dim(as.data.frame(na.omit(exprMat.plasma)))
#  1167 27850   # only 1167 nonmissing rows in site W

tSNE.list.plasma <- Rtsne::Rtsne(t(as.data.frame(na.omit(exprMat.plasma))),perplexity=20)
dim(tSNE.list.plasma$Y)
#27850     2
tSNE.plasma.xy<-as.data.frame(tSNE.list.plasma$Y)
colnames(tSNE.plasma.xy)<-c('x','y')
tSNE.4p13a.plasma.xy<-tSNE.plasma.xy
rownames(tSNE.4p13a.plasma.xy)<-rownames(numericMeta.22sites)
library(ggplot2)
#library(ggpubr) - rlang upgrade required, ggplot2 upgrade required. cannot install from older source!
library(ggrepel)

# Get the indices of the first occurrence of each unique value in contributor_code
first_occurrence_indices.4p13a<-first_occurrence_indices <- match(unique(numericMeta.22sites$contributor_code), numericMeta.22sites$contributor_code)

labels.4p13a<-labels<-numericMeta.22sites$contributor_code[first_occurrence_indices]

tSNE.plasma.samples.sites<-ggplot2::ggplot(tSNE.4p13a.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.22sites$contributor_code), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13a.plasma.xy[first_occurrence_indices.4p13a, ],
                  aes(x=x,y=y, label = labels.4p13a), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.age<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$age_at_visit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleMatrix<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_matrix), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleType<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_type), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.apoe4dose<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$APOE4.Dose), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.MMSE<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$MMSE), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.Group<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=Group), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )


pdf(file="4p13a.Starting_log2(RFU).22sites.tSNE-Plasma(7335x27850)_sitesUVW_in-samples_coloredByTraits.pdf",width=11,height=9)
  print(tSNE.plasma.samples.sites)
  print(tSNE.plasma.samples.age)
  print(tSNE.plasma.samples.sampleMatrix)
  print(tSNE.plasma.samples.sampleType)
  print(tSNE.plasma.samples.apoe4dose)
  print(tSNE.plasma.samples.MMSE)
  print(tSNE.plasma.samples.Group)
dev.off()


tSNE.4p13a.plasma.samples.sites<-tSNE.plasma.samples.sites
#retained plot for reference


##################################################
# 4p13a2. Initial QC (tSNE plots)  - no site W

#rootdir="z:/EBD/"
##rootdir="c:/Users/edammer_4ceb6ff/Downloads/"
#setwd(rootdir)

numericMeta.21sites<-numericMeta<-numericMeta.22sites[which(!numericMeta.22sites$contributor_code %in% c("W")),]
cleanDat<-cleanDat[,match(rownames(numericMeta.21sites),colnames(cleanDat))]


####################
## Examine tSNE of PLASMA - Human only assays (rows) and samples (columns);  no missing data by row.

#Group.3mappedCohortsPlusCTimputed)
Group=as.factor(numericMeta$Group.withCTimputed)

numericMeta.plasma<-numericMeta.21sites
dim(numericMeta.plasma)
#27755    96   #previously (16 sites with Calibrators):  18739   93
exprMat.plasma<-cleanDat

dim(as.data.frame(na.omit(exprMat.plasma)))
#  3415 27755   # only 3415 nonmissing rows in site V; 5031 in H; 7334 in K; 7335 in all others.

tSNE.list.plasma <- Rtsne::Rtsne(t(as.data.frame(na.omit(exprMat.plasma))),perplexity=20)
dim(tSNE.list.plasma$Y)
#27850     2
tSNE.plasma.xy<-as.data.frame(tSNE.list.plasma$Y)
colnames(tSNE.plasma.xy)<-c('x','y')
tSNE.4p13a2.plasma.xy<-tSNE.plasma.xy
rownames(tSNE.4p13a2.plasma.xy)<-rownames(numericMeta.21sites)
library(ggplot2)
#library(ggpubr) - rlang upgrade required, ggplot2 upgrade required. cannot install from older source!
library(ggrepel)

# Get the indices of the first occurrence of each unique value in contributor_code
first_occurrence_indices.4p13a2<-first_occurrence_indices <- match(unique(numericMeta.21sites$contributor_code), numericMeta.21sites$contributor_code)

labels.4p13a2<-labels<-numericMeta.21sites$contributor_code[first_occurrence_indices.4p13a2]

tSNE.plasma.samples.sites<-ggplot2::ggplot(tSNE.4p13a2.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.21sites$contributor_code), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13a2.plasma.xy[first_occurrence_indices.4p13a2, ],
                  aes(x=x,y=y, label = labels.4p13a2), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.age<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$age_at_visit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleMatrix<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_matrix), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleType<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_type), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )


# ------------------------------------------------------------------------
# ANNOTATION: Filter out sites with serum-only or low assay/sample coverage
# before regression modeling.
# ------------------------------------------------------------------------
tSNE.plasma.samples.apoe4dose<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$APOE4.Dose), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.MMSE<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$MMSE), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.Group<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=Group), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )


pdf(file="4p13a2.Starting_log2(RFU).21sites.tSNE-Plasma(7335x27755)_sitesUV_in-samples_coloredByTraits.pdf",width=11,height=9)
  print(tSNE.plasma.samples.sites)
  print(tSNE.plasma.samples.age)
  print(tSNE.plasma.samples.sampleMatrix)
  print(tSNE.plasma.samples.sampleType)
  print(tSNE.plasma.samples.apoe4dose)
  print(tSNE.plasma.samples.MMSE)
  print(tSNE.plasma.samples.Group)
dev.off()


tSNE.4p13a2.plasma.samples.sites<-tSNE.plasma.samples.sites


##################################################
# 4p13a3. Initial QC (tSNE plots)  - no site U,V,W (19 sites)

#rootdir="z:/EBD/"
##rootdir="c:/Users/edammer_4ceb6ff/Downloads/"
#setwd(rootdir)

numericMeta.19sites<-numericMeta<-numericMeta[which(!numericMeta$contributor_code %in% c("U","V")),]
cleanDat<-cleanDat[,match(rownames(numericMeta),colnames(cleanDat))]
dim(cleanDat)
#  7335 22547

####################
## Examine tSNE of PLASMA - Human only assays (rows) and samples (columns);  no missing data by row.

#Group.3mappedCohortsPlusCTimputed)
Group=as.factor(numericMeta$Group.withCTimputed)

numericMeta.plasma<-numericMeta.19sites #[which(numericMeta$sample_matrix=="CSF"),]
dim(numericMeta.plasma)
#22547    96   #previously (16 sites with Calibrators):  18739   93
exprMat.plasma<-cleanDat

dim(as.data.frame(na.omit(exprMat.plasma)))
#  5031 22547   # 5031 in H; 7334 in K; 7335 in all others.

tSNE.list.plasma <- Rtsne::Rtsne(t(as.data.frame(na.omit(exprMat.plasma))),perplexity=20)
dim(tSNE.list.plasma$Y)
#27850     2
tSNE.plasma.xy<-as.data.frame(tSNE.list.plasma$Y)
colnames(tSNE.plasma.xy)<-c('x','y')
tSNE.4p13a3.plasma.xy<-tSNE.plasma.xy
rownames(tSNE.4p13a3.plasma.xy)<-rownames(numericMeta.19sites)
library(ggplot2)
#library(ggpubr) - rlang upgrade required, ggplot2 upgrade required. cannot install from older source!
library(ggrepel)

# Get the indices of the first occurrence of each unique value in contributor_code
first_occurrence_indices.4p13a3<-first_occurrence_indices <- match(unique(numericMeta.19sites$contributor_code), numericMeta.19sites$contributor_code)

labels.4p13a3<-labels<-numericMeta.19sites$contributor_code[first_occurrence_indices.4p13a3]

tSNE.plasma.samples.sites<-ggplot2::ggplot(tSNE.4p13a3.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.19sites$contributor_code), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13a3.plasma.xy[first_occurrence_indices.4p13a3, ],
                  aes(x=x,y=y, label = labels.4p13a3), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.age<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$age_at_visit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleMatrix<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_matrix), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleType<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_type), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.apoe4dose<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$APOE4.Dose), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.MMSE<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$MMSE), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.Group<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=Group), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )


pdf(file="4p13a3.Starting_log2(RFU).19sites.tSNE-Plasma(7335x22547)-samples_coloredByTraits.pdf",width=11,height=9)
  print(tSNE.plasma.samples.sites)
  print(tSNE.plasma.samples.age)
  print(tSNE.plasma.samples.sampleMatrix)
  print(tSNE.plasma.samples.sampleType)
  print(tSNE.plasma.samples.apoe4dose)
  print(tSNE.plasma.samples.MMSE)
  print(tSNE.plasma.samples.Group)
dev.off()


tSNE.4p13a3.plasma.samples.sites<-tSNE.plasma.samples.sites


## Reduce footprint in RAM
rm(numericMeta.plasma)
rm(exprMat.plasma)


#############################################################################
## 4p13b. 2PAV Protein, intrasite, Regression for Preanalytical Factors (QC)


names(table(numericMeta.19sites$contributor_code))  # 19 sites:
#A, B, C, D, E, F, G, H, I, J, K, L, M, N, P, Q, R, S, T
save.image("4p13a4.inputForFirstRegressions.19sites.RData")
#load("4p13a4.MINIMALinputForFirstRegressions.19sites.RData")
# Reduced footprint:
cleanDat<-readRDS("4p13.cleanDat.22sites.RDS")
numericMeta.22sites<-readRDS("4p13.numericMeta.22sites.RDS")
numericMeta.19sites<-numericMeta.22sites[which(!numericMeta.22sites$contributor_code %in% c("U","V","W")),]
cleanDat<-cleanDat[,match(rownames(numericMeta.19sites),colnames(cleanDat))]
dim(cleanDat)
#  7335 22547


## Set up parallel backend
library("doParallel")
parallelThreads=4  #now Windows02  #max is number of processes that can run on your computer at one time
stopCluster(clusterLocal)
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="PSOCK")

registerDoParallel(clusterLocal)



# ------------------------------------------------------------------------
# ANNOTATION: Save the 19-site plasma dataset used as input for the first
# within-site PAV regression.
# ------------------------------------------------------------------------
normExpr.reg<-list()
for (site in names(table(numericMeta.19sites$contributor_code))[1:19]) {

  ## Subset to single site's samples
  cleanDat.unreg<-cleanDat[,which(numericMeta.19sites$contributor_code==site)]
  regvars<-data.frame(HNRNPA2B1.TimeToSpin=numericMeta.19sites$RegrBloodPreanalyticFactor.HNRNPA2B1[which(numericMeta.19sites$contributor_code==site)],
                         HBZ.TimeToDecant=numericMeta.19sites$RegrBloodPreanalyticFactor.HBZ[which(numericMeta.19sites$contributor_code==site)])

  ## Run the regression
  normExpr.reg[[site]] <- matrix(NA,nrow=nrow(cleanDat.unreg),ncol=ncol(cleanDat.unreg))
  #coefmat <- matrix(NA,nrow=nrow(cleanDat.unreg),ncol=9)  #ncol(regvars)+1) ## change this to ncol(regvars)+2 when condition has 2 levels if BOOT=TRUE, +1 if BOOT=FALSE

  #RNG seed set for reproducibility
  set.seed(1234567);
  #** coefmat.residuals.list <-
  normExpr.reg[[site]] <-  foreach (i=1:nrow(cleanDat.unreg), .combine=rbind, .packages="stats") %dopar% {  #** .combine=c, .multicombine=TRUE, .packages="stats", .export="regvars") %dopar% {
        set.seed(1234567)
        options(stringsAsFactors=FALSE)
        tryCatch({
          lmmod1 <- lm(as.numeric(cleanDat.unreg[i,])~HNRNPA2B1.TimeToSpin+HBZ.TimeToDecant,data=regvars)
          #** list(list(coef=coef(lmmod1),residuals=lmmod1$residuals))  #return a list of (1) vector of coefficients; and (2) vector of residuals length nrow(cleanDat.unreg)
          ##datpred <- predict(object=lmmod1,newdata=regvars)
          coef <- coef(lmmod1)
          #coefmat[i,] <- coef
#          coef[1] + coef[2]*regvars[,"HNRNPA2B1.TimeToSpin"] + coef[3]*regvars[,"HBZ.TimeToDecant"] + lmmod1$residuals ## The full data + with undesired covariates

# ------------------------------------------------------------------------
# ANNOTATION: Run the first-pass within-site regression of HNRNPA2B1 and HBZ
# PAV proxy assays for every protein.
# ------------------------------------------------------------------------
          coef[1] + lmmod1$residuals ## The full data - the undesired covariates
        }, error=function(e) { rep(NA_real_,ncol(cleanDat.unreg)) })
  }
  #** coefmat<-as.data.frame(do.call(rbind, lapply(coefmat.residuals.list, `[[`, "coef"))) #as.data.frame(do.call(rbind, lapply(coefmat.residuals.list,function(x){x[[1]]}))) #t(sapply(coefmat.residuals.list, `[[`, 1))
  #dim(coefmat)
  ##  7335    9
  #residualsMat<-as.data.frame(do.call(rbind,lapply(coefmat.residuals.list, `[[`, "residuals")))  #t(sapply(coefmat.residuals.list, `[[`, 2))
  #dim(residualsMat)
  ##  7335 22392
  #
  #normExpr.reg <- foreach (i=1:nrow(cleanDat.unreg), .combine=rbind) %dopar% {
  #      coefmat[i,1] + coefmat[i,2]*regvars[,"Age"] + coefmat[i,3]*as.numeric(regvars[,"Sex"]) + residualsMat[i,] ## The full data - the undesired covariates
  #      ## Also equivalent to <- thisexp - coef*var expression above
  #      #cat('Done for Protein ',i,'\n')
  #}
  rownames(normExpr.reg[[site]]) <- rownames(cleanDat.unreg)
  colnames(normExpr.reg[[site]]) <- colnames(cleanDat.unreg)

  cat(paste0("Finished 2PAV regression of site ",site,"\n"))

}
saveRDS(normExpr.reg,"normExpr.reg_sites1-19.RDS")

cleanDat.reg<-do.call(cbind, normExpr.reg)
numericMeta.reg<-numericMeta.19sites[match(colnames(cleanDat.reg),rownames(numericMeta.19sites)),]


#############################################################################
## 4p13. Variance Partition post (2PAV) Regression for Preanalytical Factors (QC)

regvars.vp<-data.frame(numericMeta.reg)
regvars.vp$Sex<-factor(abs(regvars.vp$sex -2))
regvars.vp$Age<-as.numeric(regvars.vp$age_at_visit)
regvars.vp$contributor_code<-factor(regvars.vp$contributor_code)
regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1)
regvars.vp$RegrBloodPreanalyticFactor.HBZ<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HBZ)

# too many missing values:
#form <- ~ age_at_visit+(1|Sex)+(1|raceAA)+(1|recruited_control)+(1|ad)+(1|ftd)+(1|pd)+(1|als)+(1|mci_sci)+(1|mi)+(1|C9Orf72)+(1|GRN)+(1|MAPT)+APOE4.Dose+(1|contributor_code)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+MMSE
#form <- ~ (1|contributor_code) +RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ
form <- ~ Age+(1|Sex)+(1|contributor_code)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ

library(variancePartition)

# (If regressed) REMOVE regressed proteins to avoid "Response variable 4641 has a variance of 0":
#regrProts.idx<-c(which(grepl("^HNRNPA2B1\\|",rownames(cleanDat))),which(grepl("^HBZ\\|",rownames(cleanDat))))
#regrProts.idx
#5569 6087 -- take out
#previously: 6156 4641
#integer(0)  # if already removed

#previously: removed 2 regr proteins in line:  varPart.reg <- fitExtractVarPartModel(impute::impute.knn(cleanDat[-c(4641,6156),])$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = 8, type = "SOCK"))
varPart.reg <- fitExtractVarPartModel(impute::impute.knn(as.matrix(cleanDat.reg[,]))$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp1 <- sortCols(varPart.reg,FUN=median,last= c("Residuals"))

pdf(file="4p13b.19sitesRegr1x.2PAV-VariancePartition-PLASMA-7335x22547.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp1, main="HDS 1.3ms - 4p13 - Plasma 2 PAV within Each Site Regressed" )

	SexSortOrder<-order(vp1[["Sex"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp1)) { vp1[[i]]<-vp1[[i]][SexSortOrder]; }
	rownames(vp1)<-rownames(vp1)[SexSortOrder]

	plotPercentBars( vp1[1:50,]) + ggtitle( "Top Sex-covariates" )


	AgeSortOrder<-order(vp1[["Age"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp1)) { vp1[[i]]<-vp1[[i]][AgeSortOrder]; }
	rownames(vp1)<-rownames(vp1)[AgeSortOrder]

	plotPercentBars( vp1[1:50,]) + ggtitle( "Top Age-covariates" )



# ------------------------------------------------------------------------
# ANNOTATION: Quantify variance explained after two-PAV within-site
# regression using variancePartition.
# ------------------------------------------------------------------------
        BatchSortOrder<-order(vp1[["RegrBloodPreanalyticFactor.HNRNPA2B1"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp1)) { vp1[[i]]<-vp1[[i]][BatchSortOrder]; }
        rownames(vp1)<-rownames(vp1)[BatchSortOrder]

        plotPercentBars( vp1[1:50,]) + ggtitle( "Top Time To Spin (HNRNPA2B1)-covariates" )


        BatchSortOrder<-order(vp1[["RegrBloodPreanalyticFactor.HBZ"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp1)) { vp1[[i]]<-vp1[[i]][BatchSortOrder]; }
        rownames(vp1)<-rownames(vp1)[BatchSortOrder]

        plotPercentBars( vp1[1:50,]) + ggtitle( "Top Preanalytical Factor 2 (HBZ)-covariates" )


        BatchSortOrder<-order(vp1[["contributor_code"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp1)) { vp1[[i]]<-vp1[[i]][BatchSortOrder]; }
        rownames(vp1)<-rownames(vp1)[BatchSortOrder]

        plotPercentBars( vp1[1:50,]) + ggtitle( "Top Contributor (site)-covariates" )


#	BatchSortOrder<-order(vp1[["MMSE"]],decreasing=TRUE)
#	#rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
#	for (i in ls(vp1)) { vp1[[i]]<-vp1[[i]][BatchSortOrder]; }
#	rownames(vp1)<-rownames(vp1)[BatchSortOrder]
#
#	plotPercentBars( vp1[1:50,]) + ggtitle( "MMSE Cog. Score-covariates" )

dev.off()


saveRDS(varPart.reg,"4p13b.1x19sitesRegr2PAVs.varPart.reg.RDS")


stopCluster(clusterLocal)


## Set up parallel backend
library("doParallel")
parallelThreads=8  #now Windows02  #max is number of processes that can run on your computer at one time
#stopCluster(clusterLocal)
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="PSOCK")

registerDoParallel(clusterLocal)


#############################################################################
## 4p13. Variance Partition unregressed (pre 2PAV regression) for Preanalytical Factors (QC)

regvars.vp<-data.frame(numericMeta.19sites)
regvars.vp$Sex<-factor(abs(regvars.vp$sex -2))
regvars.vp$Age<-as.numeric(regvars.vp$age_at_visit)
regvars.vp$contributor_code<-factor(regvars.vp$contributor_code)
regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1)
regvars.vp$RegrBloodPreanalyticFactor.HBZ<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HBZ)

# too many missing values:
#form <- ~ age_at_visit+(1|Sex)+(1|raceAA)+(1|recruited_control)+(1|ad)+(1|ftd)+(1|pd)+(1|als)+(1|mci_sci)+(1|mi)+(1|C9Orf72)+(1|GRN)+(1|MAPT)+APOE4.Dose+(1|contributor_code)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+MMSE
#form <- ~ (1|contributor_code) +RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ
form <- ~ Age+(1|Sex)+(1|contributor_code)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ

library(variancePartition)

# (If regressed) REMOVE regressed proteins to avoid "Response variable 4641 has a variance of 0":
#regrProts.idx<-c(which(grepl("^HNRNPA2B1\\|",rownames(cleanDat))),which(grepl("^HBZ\\|",rownames(cleanDat))))
#regrProts.idx
#5569 6087 -- take out
#previously: 6156 4641
#integer(0)  # if already removed

#previously: removed 2 regr proteins in line:  varPart.reg <- fitExtractVarPartModel(impute::impute.knn(cleanDat[-c(4641,6156),])$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = 8, type = "SOCK"))
varPart.unreg <- fitExtractVarPartModel(impute::impute.knn(as.matrix(cleanDat[,]))$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp0 <- sortCols(varPart.unreg,FUN=median,last= c("Residuals"))

pdf(file="4p13a3.Unregressed-VariancePartition-PLASMA-7335x22547.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp0, main="HDS 1.3ms - 4p13a3 - 19 sites - Unregressed" )

	SexSortOrder<-order(vp0[["Sex"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp0)) { vp0[[i]]<-vp0[[i]][SexSortOrder]; }
	rownames(vp0)<-rownames(vp0)[SexSortOrder]

	plotPercentBars( vp0[1:50,]) + ggtitle( "Top Sex-covariates" )


	AgeSortOrder<-order(vp0[["Age"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp0)) { vp0[[i]]<-vp0[[i]][AgeSortOrder]; }
	rownames(vp0)<-rownames(vp0)[AgeSortOrder]

	plotPercentBars( vp0[1:50,]) + ggtitle( "Top Age-covariates" )


        BatchSortOrder<-order(vp0[["RegrBloodPreanalyticFactor.HNRNPA2B1"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp0)) { vp0[[i]]<-vp0[[i]][BatchSortOrder]; }
        rownames(vp0)<-rownames(vp0)[BatchSortOrder]

        plotPercentBars( vp0[1:50,]) + ggtitle( "Top Time To Spin (HNRNPA2B1)-covariates" )


        BatchSortOrder<-order(vp0[["RegrBloodPreanalyticFactor.HBZ"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp0)) { vp0[[i]]<-vp0[[i]][BatchSortOrder]; }
        rownames(vp0)<-rownames(vp0)[BatchSortOrder]

        plotPercentBars( vp0[1:50,]) + ggtitle( "Top Preanalytical Factor 2 (HBZ)-covariates" )


        BatchSortOrder<-order(vp0[["contributor_code"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp0)) { vp0[[i]]<-vp0[[i]][BatchSortOrder]; }
        rownames(vp0)<-rownames(vp0)[BatchSortOrder]

        plotPercentBars( vp0[1:50,]) + ggtitle( "Top Contributor (site)-covariates" )


#	BatchSortOrder<-order(vp0[["MMSE"]],decreasing=TRUE)
#	#rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
#	for (i in ls(vp0)) { vp0[[i]]<-vp0[[i]][BatchSortOrder]; }
#	rownames(vp0)<-rownames(vp0)[BatchSortOrder]
#
#	plotPercentBars( vp0[1:50,]) + ggtitle( "MMSE Cog. Score-covariates" )

dev.off()

varPart.unreg<-vp0
saveRDS(varPart.unreg,"4p13a3.19sitesUnregressed.varPart.unreg.RDS")


##################################################
# 4p13b. 19 sites (no site U,V,W) QC (tSNE plots) - following 2PAV protein regression

#rootdir="z:/EBD/"
##rootdir="c:/Users/edammer_4ceb6ff/Downloads/"
#setwd(rootdir)


####################
## Examine tSNE of PLASMA - Human only assays (rows) and samples (columns);  no missing data by row.

#Group.3mappedCohortsPlusCTimputed)
Group=as.factor(numericMeta.reg$Group.withCTimputed)

numericMeta.plasma<-numericMeta.reg #[which(numericMeta$sample_matrix=="CSF"),]
dim(numericMeta.plasma)
#22547   96   #previously (16 sites with Calibrators):  18739   93
exprMat.plasma<-cleanDat.reg

dim(as.data.frame(na.omit(exprMat.plasma)))
#  5031 22547   # 5031 in H; 7334 in K; 7335 in all others.

tSNE.list.plasma <- Rtsne::Rtsne(t(as.data.frame(na.omit(exprMat.plasma))),perplexity=20)
dim(tSNE.list.plasma$Y)
#27850     2
tSNE.plasma.xy<-as.data.frame(tSNE.list.plasma$Y)
colnames(tSNE.plasma.xy)<-c('x','y')
tSNE.4p13b.plasma.xy<-tSNE.plasma.xy
rownames(tSNE.4p13b.plasma.xy)<-rownames(numericMeta.19sites)
library(ggplot2)
#library(ggpubr) - rlang upgrade required, ggplot2 upgrade required. cannot install from older source!
library(ggrepel)

# Get the indices of the first occurrence of each unique value in contributor_code
first_occurrence_indices.4p13b<-first_occurrence_indices <- match(unique(numericMeta.reg$contributor_code), numericMeta.reg$contributor_code)

labels.4p13b<-labels<-numericMeta.reg$contributor_code[first_occurrence_indices.4p13b]

tSNE.plasma.samples.sites<-ggplot2::ggplot(tSNE.4p13b.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.reg$contributor_code), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13b.plasma.xy[first_occurrence_indices.4p13b, ],
                  aes(x=x,y=y, label = labels.4p13b), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.age<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$age_at_visit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleMatrix<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_matrix), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleType<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_type), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.apoe4dose<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$APOE4.Dose), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.MMSE<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$MMSE), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.Group<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=Group), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )


pdf(file="4p13b.2PAVproteinIntrasiteRegressed_log2(RFU).19sites.tSNE-Plasma(7335x22547)-samples_coloredByTraits.pdf",width=11,height=9)
  print(tSNE.plasma.samples.sites)
  print(tSNE.plasma.samples.age)
  print(tSNE.plasma.samples.sampleMatrix)
  print(tSNE.plasma.samples.sampleType)
  print(tSNE.plasma.samples.apoe4dose)
  print(tSNE.plasma.samples.MMSE)
  print(tSNE.plasma.samples.Group)
dev.off()


tSNE.4p13b.plasma.samples.sites<-tSNE.plasma.samples.sites


## Reduce footprint in RAM
rm(numericMeta.plasma)
rm(exprMat.plasma)
#rm(normExpr.reg)

#save.image("4p13b.1xRegr.2PAVproteinRegr.saved.image.PLASMA_(7335x22547)_noWGCNA.RData")
#load("4p12.regr.6SVregr.saved.image.PLASMA_(7335x22392)_noWGCNA.RData")


#xy coordinates in tSNE:  tSNE.4p13b.plasma.xy
#segment1
x1=-12.5

# ------------------------------------------------------------------------
# ANNOTATION: Use t-SNE to evaluate residual site structure after two-PAV
# regression.
# ------------------------------------------------------------------------
y1=20
x2=1.5
y2=12.5

#segment2
x1.2=1.5
y1.2=12.5
x2.2=6
y2.2=19

#segment3
x1.3=1.5
y1.3=12.5
x2.3=6.7
y2.3=4

#Visualize the division lines between subsets of the F site samples
print(tSNE.4p13b.plasma.samples.sites + # last_plot() +  #if in session, actively plotted
  annotate("segment",
           x = x1,  y = y1,
           xend = x2, yend = y2,
           colour = "red", size = 1) +
  annotate("segment",
           x = x1.2,  y = y1.2,
           xend = x2.2, yend = y2.2,
           colour = "red", size = 1) +
  annotate("segment",
           x = x1.3,  y = y1.3,
           xend = x2.3, yend = y2.3,
           colour = "red", size = 1)
)

# tSNE.4p13b.plasma.xy is a data frame with columns "x" and "y" defining points in a tSNE plot and having the same number of rows as length of numericMeta.reg$contributor_code.
# We want to set values in a vector "F.subset" that are identical to values in numericMeta.reg$contributor_code, and replace contributor_code values of F, setting them to "F1"
# if the points are below the line defined by segment 1 and to the right of the line defined by segment 3; setting them to "F2" if the points are above the line defined by
# segment 1 and to the left of the line defined by segment 2; and finally, setting the remaining points that are to the right of the line defined by segment 2 and above the line
# defined by segment 3.

## ------------------------------------------------------------------------------
## 1.  Define the three boundary lines
## ------------------------------------------------------------------------------
seg <- list(
  s1 = cbind(x = c(x1,  x2), y = c(y1, y2)),   # segment-1   ---
  s2 = cbind(x = c(  x1.2,  x2.2), y = c(y1.2, y2.2)), # segment-2   --
  s3 = cbind(x = c(  x1.3,  x2.3), y = c(y1.3,  y2.3  ))  # segment-3   --
)

get_mb <- function(p) {                    # slope-intercept pair
  m <- diff(p[,"y"]) / diff(p[,"x"])
  b <- p[1,"y"] - m * p[1,"x"]
  c(m = m, b = b)
}
mb <- lapply(seg, get_mb)

## ------------------------------------------------------------------------------
## 2.  Helper predicates (vectorised over the whole t-SNE data set)
## ------------------------------------------------------------------------------
dat <- tSNE.4p13b.plasma.xy                # cols x, y  (same row-order as meta)

# below / above segment-1  (y versus m-·x+b-)
below1 <- with(dat, y < mb$s1["m"] * x + mb$s1["b.y"])
above1 <- !below1                            # mutually exclusive
x_on_1 <- with(dat, (y - mb$s1["b.y"]) / mb$s1["m"])
right1 <- dat$x > x_on_1

# left / right segment-2
x_on_2 <- with(dat, (y - mb$s2["b.y"]) / mb$s2["m"])
left2  <- dat$x < x_on_2
right2 <- !left2
above2 <- with(dat, y > mb$s2["m"] * x + mb$s2["b.y"])

# left / right + above segment-3
x_on_3 <- with(dat, (y - mb$s3["b.y"]) / mb$s3["m"])
right3 <- dat$x > x_on_3
left3 <- !right3
below3 <- with(dat, y < mb$s3["m"] * x + mb$s3["b.y"])
above3 <- !below3

below1and3<-above1and2<-right2and3<-right1andLeft2<-right2andAbove3<-rep(FALSE,length(numericMeta.reg$contributor_code))
below1and3[which(below1 & below3)]<-TRUE
above1and2[which(above1 & above2)]<-TRUE
right2and3[which(right2 & right3)]<-TRUE
right1andLeft2[which(right1 & left2)]<-TRUE
right2andAbove3[which(right2 & above3)]<-TRUE

## ------------------------------------------------------------------------------
## 3.  Re-label contributor_code == "F"
## ------------------------------------------------------------------------------
F.subset <- numericMeta.reg$contributor_code   # start with original labels

##   F1  :  below seg-1   &   right seg-3
sel <- which(F.subset == "F" & below1and3)
F.subset[sel] <- "F1"

##   F2  :  above seg-1   &   left  seg-2
sel <- which(F.subset == "F" & right1andLeft2)
F.subset[sel] <- "F2"

##   F3  :  right seg-2   &   above seg-3   (the remaining requested region)
sel <- which(F.subset == "F" & right2andAbove3)
F.subset[sel] <- "F3"

## F.subset now contains the updated contributor codes
table(F.subset)


## SHOW F samples as split by segments in plot
tempCol=rep("black",22547)
tempCol[F.subset=="F"]<-"maroon"
tempCol[F.subset=="F1"]<-"yellowgreen"
tempCol[F.subset=="F3"]<-"darkturquoise"

pdf(file="4p13b.2PAVproteinIntrasiteRegressed_log2(RFU).19sites.tSNE-Plasma(7335x22547)-samples_SiteF_SplitF1_F2_F3.pdf",width=11,height=9)

plot(tSNE.4p13b.plasma.xy[numericMeta.reg$contributor_code=="F",], col=tempCol[numericMeta.reg$contributor_code=="F"])
with(seg, {
  segments(s1[1, "x"], s1[1, "y"], s1[2, "x"], s1[2, "y"], col = "red", lwd = 2)
  segments(s2[1, "x"], s2[1, "y"], s2[2, "x"], s2[2, "y"], col = "red", lwd = 2)
  segments(s3[1, "x"], s3[1, "y"], s3[2, "x"], s3[2, "y"], col = "red", lwd = 2)
})

## Show in context of all intrasite 2PAV-regressed sites tSNE
print(
ggplot2::ggplot(tSNE.4p13b.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=F.subset), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13b.plasma.xy[first_occurrence_indices.4p13b, ],
                  aes(x=x,y=y, label = labels.4p13b), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  ) +
  annotate("segment",
           x = x1,  y = y1,
           xend = x2, yend = y2,
           colour = "red", size = 1) +

# ------------------------------------------------------------------------
# ANNOTATION: Manually segment site F into F1/F2/F3 sub-batches based on
# t-SNE coordinates.
# ------------------------------------------------------------------------
  annotate("segment",
           x = x1.2,  y = y1.2,
           xend = x2.2, yend = y2.2,
           colour = "red", size = 1) +
  annotate("segment",
           x = x1.3,  y = y1.3,
           xend = x2.3, yend = y2.3,
           colour = "red", size = 1)
)
dev.off()


numericMeta.reg$contributor_Fsplit<-F.subset


## 4p13b1. Set the medians within all the sites (batches, contributor_code s) to exactly the same; use site "A" as template IF injecting back a set of medians from one site to all others.

#templateMedians.unlog=2^apply(cleanDat.reg[,numericMeta.reg$contributor_code=="A"],1,function(x) median(x,na.rm=T))
#names(templateMedians.unlog)<-rownames(cleanDat.reg)
cleanRelAbun.sameMedianRank.list<-list()
for(site in names(table(numericMeta.reg$contributor_Fsplit))) cleanRelAbun.sameMedianRank.list[[site]]<-2^t(apply(cleanDat.reg[,which(numericMeta.reg$contributor_Fsplit==site)],1,function(x) x - median(x,na.rm=T)))
#for(site in names(table(numericMeta$contributor_code))) cleanRelAbun.sameMedianRank.list[[site]]<-sweep(cleanRelAbun.sameMedianRank.list[[site]],1,templateMedians.unlog,"*")
#here we strip without injecting medians. all medians within site are now ZERO
#
### sanity check: A should be unchanged
#all.equal(log2(cleanRelAbun.sameMedianRank.list[["A"]]),cleanDat[,numericMeta$contributor_code=="A"])
##TRUE


cleanDat.zeroSiteMedian<-log2(do.call("cbind", cleanRelAbun.sameMedianRank.list))
cleanDat.zeroSiteMedian<-cleanDat.zeroSiteMedian[,match(rownames(numericMeta.reg),colnames(cleanDat.zeroSiteMedian))]


##save RAM
rm(cleanRelAbun.sameMedianRank.list)
#rm(templateMedians.unlog)


commonRows<-which(apply(cleanDat.zeroSiteMedian,1,function(x) length(which(is.na(x))))==0)  # rows with no missing values
length(commonRows)
#5031

# How poor or good was median (rowwise) correlation across sites before the above fix?
cor(apply(cleanDat.reg[commonRows,numericMeta.reg$contributor_code=="F"],1,function(x) median(x,na.rm=T)),
    apply(cleanDat.reg[commonRows,numericMeta.reg$contributor_code=="A"],1,function(x) median(x,na.rm=T)), method="spearman")
#0.7908457 for the split tSNE cluster count cohort F vs. full depth cohort A (BioHermes).

cor(apply(cleanDat.zeroSiteMedian[commonRows,numericMeta.reg$contributor_code=="E"],1,function(x) median(x,na.rm=T)),
    apply(cleanDat.zeroSiteMedian[commonRows,numericMeta.reg$contributor_code=="A"],1,function(x) median(x,na.rm=T)), method="spearman")
#SD is zero.


## Set up parallel backend
library("doParallel")
parallelThreads=8  #now Windows02  #max is number of processes that can run on your computer at one time
stopCluster(clusterLocal)
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="PSOCK")

registerDoParallel(clusterLocal)


#############################################################################
## 4p13b1. Variance Partition regressed (2PAV regression intrasite)+Sitewise RowMedian Zeroed (QC)

regvars.vp<-data.frame(numericMeta.reg)
regvars.vp$Sex<-factor(abs(regvars.vp$sex -2))
regvars.vp$Age<-as.numeric(regvars.vp$age_at_visit)
regvars.vp$contributor_Fsplit<-factor(regvars.vp$contributor_Fsplit)
regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1)
regvars.vp$RegrBloodPreanalyticFactor.HBZ<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HBZ)

# too many missing values:
#form <- ~ age_at_visit+(1|Sex)+(1|raceAA)+(1|recruited_control)+(1|ad)+(1|ftd)+(1|pd)+(1|als)+(1|mci_sci)+(1|mi)+(1|C9Orf72)+(1|GRN)+(1|MAPT)+APOE4.Dose+(1|contributor_code)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+MMSE
#form <- ~ (1|contributor_code) +RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ
form <- ~ Age+(1|Sex)+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ

library(variancePartition)

# (If regressed) REMOVE regressed proteins to avoid "Response variable 4641 has a variance of 0":
#regrProts.idx<-c(which(grepl("^HNRNPA2B1\\|",rownames(cleanDat))),which(grepl("^HBZ\\|",rownames(cleanDat))))
#regrProts.idx
#5569 6087 -- take out
#previously: 6156 4641
#integer(0)  # if already removed

#previously: removed 2 regr proteins in line:  varPart.reg <- fitExtractVarPartModel(impute::impute.knn(cleanDat[-c(4641,6156),])$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = 8, type = "SOCK"))
varPart.b1 <- fitExtractVarPartModel(impute::impute.knn(as.matrix(cleanDat.zeroSiteMedian[,]))$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp.b1 <- sortCols(varPart.b1,FUN=median,last= c("Residuals"))

pdf(file="4p13b1.19sites1xPAVregr+SitewiseZeroRowMedians-VariancePartition-PLASMA-7335x22547.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp.b1, main="HDS 1.3ms - 4p13b1 - 19 sites - 1xRegr + SiteRowMediansZero" )

	SexSortOrder<-order(vp.b1[["Sex"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b1)) { vp.b1[[i]]<-vp.b1[[i]][SexSortOrder]; }
	rownames(vp.b1)<-rownames(vp.b1)[SexSortOrder]

	plotPercentBars( vp.b1[1:50,]) + ggtitle( "Top Sex-covariates" )


	AgeSortOrder<-order(vp.b1[["Age"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b1)) { vp.b1[[i]]<-vp.b1[[i]][AgeSortOrder]; }
	rownames(vp.b1)<-rownames(vp.b1)[AgeSortOrder]

	plotPercentBars( vp.b1[1:50,]) + ggtitle( "Top Age-covariates" )


        BatchSortOrder<-order(vp.b1[["RegrBloodPreanalyticFactor.HNRNPA2B1"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b1)) { vp.b1[[i]]<-vp.b1[[i]][BatchSortOrder]; }
        rownames(vp.b1)<-rownames(vp.b1)[BatchSortOrder]

        plotPercentBars( vp.b1[1:50,]) + ggtitle( "Top Time To Spin (HNRNPA2B1)-covariates" )



# ------------------------------------------------------------------------
# ANNOTATION: Evaluate the alternative site-wise median-zeroing strategy by
# variance partitioning and t-SNE.
# ------------------------------------------------------------------------
        BatchSortOrder<-order(vp.b1[["RegrBloodPreanalyticFactor.HBZ"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b1)) { vp.b1[[i]]<-vp.b1[[i]][BatchSortOrder]; }
        rownames(vp.b1)<-rownames(vp.b1)[BatchSortOrder]

        plotPercentBars( vp.b1[1:50,]) + ggtitle( "Top Preanalytical Factor 2 (HBZ)-covariates" )


        BatchSortOrder<-order(vp.b1[["contributor_Fsplit"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b1)) { vp.b1[[i]]<-vp.b1[[i]][BatchSortOrder]; }
        rownames(vp.b1)<-rownames(vp.b1)[BatchSortOrder]

        plotPercentBars( vp.b1[1:50,]) + ggtitle( "Top Contributor (site)-covariates" )


#	BatchSortOrder<-order(vp.b1[["MMSE"]],decreasing=TRUE)
#	#rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
#	for (i in ls(vp.b1)) { vp.b1[[i]]<-vp.b1[[i]][BatchSortOrder]; }
#	rownames(vp.b1)<-rownames(vp.b1)[BatchSortOrder]
#
#	plotPercentBars( vp.b1[1:50,]) + ggtitle( "MMSE Cog. Score-covariates" )

dev.off()

varPart.b1<-vp.b1
saveRDS(varPart.b1,"4p13b1.19sites1xPAVregr+SitewiseZeroRowMedians.varPart.unreg.RDS")
#note:  varPart.b1_noFsplit saved prior without splitting site F above
#tSNE.4p13b1.plasma.samples.sites_noFsplit<-tSNE.4p13b1.plasma.samples.sites


##################################################
# 4p13b1. 19 sites (no site U,V,W) 1xPAVregr+sitewiseZeroRowMedians QC (tSNE plots)

#rootdir="z:/EBD/"
##rootdir="c:/Users/edammer_4ceb6ff/Downloads/"
#setwd(rootdir)


####################
## Examine tSNE of PLASMA - Human only assays (rows) and samples (columns);  no missing data by row.

#Group.3mappedCohortsPlusCTimputed)
Group=as.factor(numericMeta.reg$Group.withCTimputed)

numericMeta.plasma<-numericMeta.reg #[which(numericMeta$sample_matrix=="CSF"),]
dim(numericMeta.plasma)
#22547   96   #previously (16 sites with Calibrators):  18739   93
exprMat.plasma<-cleanDat.zeroSiteMedian

dim(as.data.frame(na.omit(exprMat.plasma)))
#  5031 22547   # 5031 in H; 7334 in K; 7335 in all others.

tSNE.list.plasma <- Rtsne::Rtsne(t(as.data.frame(na.omit(exprMat.plasma))),perplexity=20)
dim(tSNE.list.plasma$Y)
#27850     2
tSNE.plasma.xy<-as.data.frame(tSNE.list.plasma$Y)
colnames(tSNE.plasma.xy)<-c('x','y')
tSNE.4p13b1.plasma.xy<-tSNE.plasma.xy
rownames(tSNE.4p13b1.plasma.xy)<-rownames(numericMeta.reg)
library(ggplot2)
#library(ggpubr) - rlang upgrade required, ggplot2 upgrade required. cannot install from older source!
library(ggrepel)

# Get the indices of the first occurrence of each unique value in contributor_Fsplit
first_occurrence_indices.4p13b1<-first_occurrence_indices <- match(unique(numericMeta.reg$contributor_Fsplit), numericMeta.reg$contributor_Fsplit)

labels.4p13b1<-labels<-numericMeta.reg$contributor_Fsplit[first_occurrence_indices.4p13b1]

tSNE.plasma.samples.sites<-ggplot2::ggplot(tSNE.4p13b1.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.reg$contributor_Fsplit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13b1.plasma.xy[first_occurrence_indices.4p13b1, ],
                  aes(x=x,y=y, label = labels.4p13b1), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.age<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$age_at_visit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleMatrix<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_matrix), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleType<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_type), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.apoe4dose<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$APOE4.Dose), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.MMSE<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$MMSE), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.Group<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=Group), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )


pdf(file="4p13b1.2PAVproteinIntrasiteRegressed+ZeroSitewiseRowMedians.19sites.tSNE-Plasma(7335x22547)-samples_coloredByTraits.pdf",width=11,height=9)
  print(tSNE.plasma.samples.sites)
  print(tSNE.plasma.samples.age)
  print(tSNE.plasma.samples.sampleMatrix)
  print(tSNE.plasma.samples.sampleType)
  print(tSNE.plasma.samples.apoe4dose)
  print(tSNE.plasma.samples.MMSE)
  print(tSNE.plasma.samples.Group)
dev.off()


tSNE.4p13b1.plasma.samples.sites<-tSNE.plasma.samples.sites


## Single site per page tSNE replots
numericMeta.plasma<-numericMeta.reg
dim(numericMeta.plasma)
library(ggplot2)
library(ggrepel)

## tSNE plot for Each site alone - colored by Age (19 pages)
tSNE.plasma.samples.age.1site<-Age<-list()
for(site in names(table(numericMeta.plasma$contributor_Fsplit))) {

  Age[[site]]=numericMeta.plasma$age_at_visit[which(numericMeta.plasma$contributor_Fsplit==site)]

  tSNE.plasma.samples.age.1site[[site]]<-ggplot2::ggplot(tSNE.plasma.xy[which(numericMeta.plasma$contributor_Fsplit==site),],label=1:length(which(numericMeta.plasma$contributor_Fsplit==site))) + geom_point(aes(x=x,y=y, color=Age[[site]]), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
    geom_text_repel(data = tSNE.plasma.xy[which(numericMeta.plasma$contributor_Fsplit==site)[1], ],
                    aes(x=x,y=y, label = site), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
    theme_minimal() +  # Minimal theme
    theme(
      panel.background = element_blank(),  # Remove plot area color
  #    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
      legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
      legend.key = element_rect(fill = "white"),  # Keep legend keys clean
      axis.title.x = element_text(size = 28),  # Double x-axis label text size
      axis.title.y = element_text(size = 28)
    )
}


pdf(file="4p13b1.AgeScale_tSNE-SingleSitePerPage.2PAVproteinIntrasiteRegressed+ZeroSitewiseRowMedians.19sites-Plasma(7335x22547).pdf",width=11,height=9)
  for(site in names(table(numericMeta.plasma$contributor_Fsplit))) print(tSNE.plasma.samples.age.1site[[site]])
dev.off()


## Redo-color by Dx group (with CT imputed and 3 Emory study sites updated Dx group.pathCog
tSNE.plasma.samples.group.1site<-GroupDx<-list()
for(site in names(table(numericMeta.plasma$contributor_Fsplit))) {

  #Age[[site]]=numericMeta.plasma$age_at_visit[which(numericMeta.plasma$contributor_Fsplit==site)]
  GroupDx[[site]]=as.factor(numericMeta$Group.withCTimputed[which(numericMeta.plasma$contributor_Fsplit==site)])

  tSNE.plasma.samples.group.1site[[site]]<-ggplot2::ggplot(tSNE.plasma.xy[which(numericMeta.plasma$contributor_Fsplit==site),],label=1:length(which(numericMeta.plasma$contributor_Fsplit==site))) + geom_point(aes(x=x,y=y, color=GroupDx[[site]]), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
    geom_text_repel(data = tSNE.plasma.xy[which(numericMeta.plasma$contributor_Fsplit==site)[1], ],
                    aes(x=x,y=y, label = site), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
    theme_minimal() +  # Minimal theme
    theme(
      panel.background = element_blank(),  # Remove plot area color
  #    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
      legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
      legend.key = element_rect(fill = "white"),  # Keep legend keys clean
      axis.title.x = element_text(size = 28),  # Double x-axis label text size
      axis.title.y = element_text(size = 28)
    )
}

pdf(file="4p13b1.DxGroupScale_tSNE-SingleSitePerPage.2PAVproteinIntrasiteRegressed+ZeroSitewiseRowMedians.19sites-Plasma(7335x22547).pdf",width=11,height=9)
  for(site in names(table(numericMeta.plasma$contributor_Fsplit))) print(tSNE.plasma.samples.group.1site[[site]])
dev.off()

#####################################################################
## Reduce footprint in RAM
rm(numericMeta.plasma)
rm(exprMat.plasma)
#rm(normExpr.reg)


cleanDat<-cleanDat.reg
save.image("4p13b1.inputForSecondRegressions.19sites.RData")

cleanDat<-cleanDat.reg


#####################################################################
## Check APOE e4 carrier best correlates in data - we will protect APOE e4 carrier status, but do not have all sample genotypes
numericMeta.reg$APOE4.carrier<-NA
numericMeta.reg$APOE4.carrier[numericMeta.reg$APOE4.Dose==0]<-0
numericMeta.reg$APOE4.carrier[numericMeta.reg$APOE4.Dose>0]<-1

table(numericMeta.reg$APOE4.carrier)
#   0    1
#9225 5838

library(WGCNA)
#e4.bicor.to.siteCorr.assays<-apply(cleanDat.zeroSiteMedian,1,function(x) bicor(x,numericMeta.reg$APOE4.carrier, use='p'))
e4.bicor.to.siteCorr.assays<-bicor(t(cleanDat.zeroSiteMedian),numericMeta.reg$APOE4.carrier, use='p')
#sort(e4.bicor.to.siteCorr.assays)
e4.bicor.to.siteCorr.assays[order(unlist(t(e4.bicor.to.siteCorr.assays)),decreasing=TRUE),][c(1:5,7331:7335)]
#SPC25|Q9HBM1   LRRN|Q6UXK5^SL025922@seq.11293.14   CTF1|Q16619^SL002783@seq.13732.79               TBCA|O75347
#   0.7466951                           0.6161494                           0.4517782      ...       -0.5981416

e4.bicor.tointrasiteCorrOnly.assays<-bicor(t(cleanDat.reg),numericMeta.reg$APOE4.carrier, use='p')
e4.bicor.to.intrasiteCorrOnly.assays[order(unlist(t(e4.bicor.to.intrasiteCOrrOnly.assays)),decreasing=TRUE),][c(1:5,7331:7335)]
#SPC25|Q9HBM1   LRRN|Q6UXK5^SL025922@seq.11293.14                      EPB41L1|Q9H4G0               TBCA|O75347             S100A13|Q99584
#   0.5754525                           0.3670285                           0.1583231      ...       -0.3557982                 -0.4231963
#####################################################################


## Start round b (second) regressions (lm):
# 4p13b2: site (factored, A as reference level) with F split to F1, F2, F3 -- no protect
# 4p13b3: site + protect Age+Sex (lose 155 samples; new numericMeta.4p13b.AgeSexNoNA)
# 4p13b4: site + protect Age+Sex+SPC25 (from cleanDat.reg, only intrasite PAV regressed)
# 4p13b5: site + protect Age+Sex+(SPC25-TBCA); APOE e4 inversely correlated pair, equiv to log2(ratio))
#############################################################################
## 4p13b2. Site Regression for intersite Preanalytical Factors (no protection)

names(table(numericMeta.reg$contributor_Fsplit))  # 22 sites (F split into F1, F2, F3):
#A, B, C, D, E, F1, F2, F3, G, H, I, J, K, L, M, N, P, Q, R, S, T
cleanDat<-cleanDat.reg
cleanDat<-cleanDat[,match(rownames(numericMeta.reg),colnames(cleanDat))]  # numericMeta.19sites - different order
dim(cleanDat)
#  7335 22547


## Set up parallel backend
library("doParallel")
parallelThreads=31  #now Windows02  #max is number of processes that can run on your computer at one time
stopCluster(clusterLocal)
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="PSOCK")

registerDoParallel(clusterLocal)

  ## regression variables for first of the second regression passes
  cleanDat.unreg<-cleanDat
  Sex=as.integer(abs(regvars.vp$sex -2))  #0=F; 1=M  -- will only be used as factor
  regvars<-data.frame(Site=factor(numericMeta.reg$contributor_Fsplit), Age=as.numeric(numericMeta.reg$age_at_visit), Sex=relevel(factor(Sex), ref="0"), Sex.int=as.integer(Sex),
                      APOEe4.LRRN1=as.numeric(cleanDat.unreg["LRRN1|Q6UXK5^SL025922@seq.11293.14",]), APOEe4.LRRN1.TBCAlog2Ratio=as.numeric(cleanDat.unreg["LRRN1|Q6UXK5^SL025922@seq.11293.14",] - cleanDat.unreg["TBCA|O75347",]))


  ## Run the regression (4p13 b2) - Site with no protection
  normExpr.reg <- matrix(NA,nrow=nrow(cleanDat.unreg),ncol=ncol(cleanDat.unreg), dimnames = dimnames(cleanDat.unreg))
  #coefmat <- matrix(NA,nrow=nrow(cleanDat.unreg),ncol=9)  #ncol(regvars)+1) ## change this to ncol(regvars)+2 when condition has 2 levels if BOOT=TRUE, +1 if BOOT=FALSE

  normExpr.reg <-  foreach (i=1:nrow(cleanDat.unreg), .combine=rbind, .packages="stats") %dopar% {  #** .combine=c, .multicombine=TRUE, .packages="stats", .export="regvars") %dopar% {
    y <- as.numeric(cleanDat.unreg[i, ])
    keep <- which(!is.na(y))             # columns that have a value

    ## initialise result for this row as all NA
    adj <- rep(NA_real_, length(y))

    if (length(keep) > 1) {              # need >=2 points to fit lm
        fit <- tryCatch(
            lm(y[keep]~Site, data = regvars[keep, , drop = FALSE]),
            error = function(e) NULL)

        if (!is.null(fit)) {
            coef <- coef(fit)
            ## coefficient[1] + residuals -> adjusted expression
            adj[keep] <- coef[1] + residuals(fit)
        }
    }
    adj                                    # returned to foreach
  }
  dimnames(normExpr.reg) <- dimnames(cleanDat.unreg)

  cat(paste0("Finished Pass 4p13 b2 regression of intersite variance.\n"))

# sanity check: some sites have missing values:
table(apply(normExpr.reg,2,function(x) length(which(is.na(x)))))
#    0     1  2303  <- # of missing values in column
#19887   406  2254  <- # of columns with that many NA

range(normExpr.reg,na.rm=T)
# -23.65559  42.18707

saveRDS(normExpr.reg,"4p13b2.normExpr.reg_sites1-19_Fsplit.RDS")

# ------------------------------------------------------------------------
# ANNOTATION: Fit second-pass site regression models with site F split and
# different protected covariates.
# ------------------------------------------------------------------------
cleanDat.4p13b2<-normExpr.reg
#numericMeta.reg still valid, ordered for this cleanDat


  ## regression variables for b3, b4, b5 second regression passes
  cleanDat.unreg<-cleanDat
  Sex=as.integer(abs(regvars.vp$sex -2))  #0=F; 1=M  -- will only be used as factor
  regvars.b345<-data.frame(Site=factor(numericMeta.reg$contributor_Fsplit), Age=as.numeric(numericMeta.reg$age_at_visit), Sex=relevel(factor(Sex), ref="0"), Sex.int=as.integer(Sex),
                           APOEe4.LRRN1=as.numeric(cleanDat.unreg["LRRN1|Q6UXK5^SL025922@seq.11293.14",]), APOEe4.LRRN1.TBCAlog2Ratio=as.numeric(cleanDat.unreg["LRRN1|Q6UXK5^SL025922@seq.11293.14",] - cleanDat.unreg["TBCA|O75347",]))
  rownames(regvars.b345)<-rownames(numericMeta.reg)
  regvars.b345<-na.omit(regvars.b345)  # removes samples with missing age, sex
  dim(regvars.b345)
  # 22392     6
  cleanDat.unreg.b345<-cleanDat.unreg[,match(rownames(regvars.b345),colnames(cleanDat.unreg))]


  ##  Run the regression (4p13 b3) - Site with Age+Sex protection

  ## covariate names you want to protect
  cov_keep <- c("Age", "Sex")   # Technically, we cleared out NAs, but one could list protection covariates that have missing values here

  normExpr.reg <- matrix(NA,nrow=nrow(cleanDat.unreg.b345),ncol=ncol(cleanDat.unreg.b345), dimnames = dimnames(cleanDat.unreg.b345))
  good_samp=which(complete.cases(regvars.b345[, cov_keep]))
  normExpr.reg <-  foreach (i=1:nrow(cleanDat.unreg.b345), .combine=rbind, .packages="stats") %dopar% {  #** .combine=c, .multicombine=TRUE, .packages="stats", .export="regvars") %dopar% {
    y <- as.numeric(cleanDat.unreg.b345[i, ])
    ## keep = samples where y and all covariates are present
    keep <- intersect(which(!is.na(y)), good_samp)             # columns that have a value, and regress with variables that have a value (we know Age and Sex are missing some)

    ## initialise result for this row as all NA
    adj <- rep(NA_real_, length(y))

    fit <- tryCatch(
        lm(y[keep]~Age+Sex +Site, data = regvars.b345[keep, , drop = FALSE]),
        error = function(e) NULL)

    if (!is.null(fit)) {
        coef <- coef(fit)
        ## coefficient[1] + residuals -> adjusted expression
        adj[keep] <- coef["(Intercept)"] + coef["Age"]*regvars.b345[keep,"Age"] + coef[which(grepl("^Sex", names(coef)))]*abs(regvars.b345[keep,"Sex.int"]) + residuals(fit)
    }

#        if (!is.null(fit)) {
#            ## build a model matrix with Site *zeroed out*
#            X     <- model.matrix(fit)
#            X[ , grep("^Site", colnames(X)) ] <- 0     # drop site contribution
#            y_hat <- as.numeric(X %*% coef(fit))       # fitted without Site
#
#            ## protected expression = y_hat + residuals  ( == y - Site effect )
#            adj[keep] <- y_hat + residuals(fit)
#        }
    adj                                    # returned to foreach
  }
  dimnames(normExpr.reg) <- dimnames(cleanDat.unreg.b345)

  cat(paste0("Finished Pass 4p13 b3 regression of intersite variance.\n"))

# sanity check: some sites have missing values:
table(apply(normExpr.reg,2,function(x) length(which(is.na(x)))))
#    0     1  2303  <- # of missing values in column
#19732   406  2254  <- # of columns with that many NA

range(normExpr.reg,na.rm=T)
# -23.65559  42.18470 (above code, using explicit sex, Males as reference level)
# -23.65559  42.18470 (above code, using y_hat with missingness in check [commented])
# -23.63702  42.18528 (prior with explicit sex component addition, using Females as non-reference level)

saveRDS(normExpr.reg,"4p13b3.normExpr.reg_sites1-19_Fsplit.RDS")
cleanDat.4p13b3<-normExpr.reg
#numericMeta.reg.b345 still valid, ordered for this cleanDat


  ##  Run the regression (4p13 b4) - Site with Age+Sex+APOE.e4 (LRRN1) protection

  ## covariate names you want to protect
  cov_keep <- c("Age", "Sex", "APOEe4.LRRN1")   # Technically, we cleared out NAs, but one could list protection covariates that have missing values here

  normExpr.reg <- matrix(NA,nrow=nrow(cleanDat.unreg.b345),ncol=ncol(cleanDat.unreg.b345), dimnames = dimnames(cleanDat.unreg.b345))
  good_samp=which(complete.cases(regvars.b345[, cov_keep]))
  normExpr.reg <-  foreach (i=1:nrow(cleanDat.unreg.b345), .combine=rbind, .packages="stats") %dopar% {  #** .combine=c, .multicombine=TRUE, .packages="stats", .export="regvars") %dopar% {
    y <- as.numeric(cleanDat.unreg.b345[i, ])
    ## keep = samples where y and all covariates are present
    keep <- intersect(which(!is.na(y)), good_samp)             # columns that have a value, and regress with variables that have a value (we know Age and Sex are missing some)

    ## initialise result for this row as all NA
    adj <- rep(NA_real_, length(y))

    fit <- tryCatch(
        lm(y[keep]~Age+Sex+APOEe4.LRRN1 +Site, data = regvars.b345[keep, , drop = FALSE]),
        error = function(e) NULL)

    if (!is.null(fit)) {
        coef <- coef(fit)
        ## coefficient[1] + residuals -> adjusted expression
        adj[keep] <- coef["(Intercept)"] + coef["Age"]*regvars.b345[keep,"Age"] + coef[which(grepl("^Sex", names(coef)))]*abs(regvars.b345[keep,"Sex.int"]) + coef["APOEe4.LRRN1"]*regvars.b345[keep,"APOEe4.LRRN1"] + residuals(fit)
    }

#        if (!is.null(fit)) {
#            ## build a model matrix with Site *zeroed out*
#            X     <- model.matrix(fit)
#            X[ , grep("^Site", colnames(X)) ] <- 0     # drop site contribution
#            y_hat <- as.numeric(X %*% coef(fit))       # fitted without Site
#
#            ## protected expression = y_hat + residuals  ( == y - Site effect )
#            adj[keep] <- y_hat + residuals(fit)
#        }
    adj                                    # returned to foreach
  }
  dimnames(normExpr.reg) <- dimnames(cleanDat.unreg.b345)

  cat(paste0("Finished Pass 4p13 b4 regression of intersite variance.\n"))

# sanity check: some sites have missing values:
table(apply(normExpr.reg,2,function(x) length(which(is.na(x)))))
#    0     1  2303  <- # of missing values in column
#19732   406  2254  <- # of columns with that many NA

saveRDS(normExpr.reg,"4p13b4.normExpr.reg_sites1-19_Fsplit.RDS")
cleanDat.4p13b4<-normExpr.reg
#numericMeta.reg.b345 still valid, ordered for this cleanDat


  ##  Run the regression (4p13 b5) - Site with Age+Sex+APOE.e4 log2(LRRN1/TBCA) protection

  ## covariate names you want to protect
  cov_keep <- c("Age", "Sex", "APOEe4.LRRN1.TBCAlog2Ratio")   # Technically, we cleared out NAs, but one could list protection covariates that have missing values here

  normExpr.reg <- matrix(NA,nrow=nrow(cleanDat.unreg.b345),ncol=ncol(cleanDat.unreg.b345), dimnames = dimnames(cleanDat.unreg.b345))
  good_samp=which(complete.cases(regvars.b345[, cov_keep]))
  normExpr.reg <-  foreach (i=1:nrow(cleanDat.unreg.b345), .combine=rbind, .packages="stats") %dopar% {  #** .combine=c, .multicombine=TRUE, .packages="stats", .export="regvars") %dopar% {
    y <- as.numeric(cleanDat.unreg.b345[i, ])
    ## keep = samples where y and all covariates are present
    keep <- intersect(which(!is.na(y)), good_samp)             # columns that have a value, and regress with variables that have a value (we know Age and Sex are missing some)

    ## initialise result for this row as all NA
    adj <- rep(NA_real_, length(y))

    fit <- tryCatch(
        lm(y[keep]~Age+Sex+APOEe4.LRRN1.TBCAlog2Ratio +Site, data = regvars.b345[keep, , drop = FALSE]),
        error = function(e) NULL)

    if (!is.null(fit)) {
        coef <- coef(fit)
        ## coefficient[1] + residuals -> adjusted expression
        adj[keep] <- coef["(Intercept)"] + coef["Age"]*regvars.b345[keep,"Age"] + coef[which(grepl("^Sex", names(coef)))]*abs(regvars.b345[keep,"Sex.int"]) + coef["APOEe4.LRRN1.TBCAlog2Ratio"]*regvars.b345[keep,"APOEe4.LRRN1.TBCAlog2Ratio"] + residuals(fit)
    }

#        if (!is.null(fit)) {
#            ## build a model matrix with Site *zeroed out*
#            X     <- model.matrix(fit)
#            X[ , grep("^Site", colnames(X)) ] <- 0     # drop site contribution
#            y_hat <- as.numeric(X %*% coef(fit))       # fitted without Site
#
#            ## protected expression = y_hat + residuals  ( == y - Site effect )
#            adj[keep] <- y_hat + residuals(fit)
#        }
    adj                                    # returned to foreach
  }
  dimnames(normExpr.reg) <- dimnames(cleanDat.unreg)

  cat(paste0("Finished Pass 4p13 b5 regression of intersite variance.\n"))

# sanity check: some sites have missing values:
table(apply(normExpr.reg,2,function(x) length(which(is.na(x)))))
#expected:
#    0     1  2303  <- # of missing values in column
#19732   406  2254  <- # of columns with that many NA


saveRDS(normExpr.reg,"4p13b5.normExpr.reg_sites1-19_Fsplit.RDS")
cleanDat.4p13b5<-normExpr.reg
#numericMeta.reg.b345 still valid, ordered for this cleanDat


##########################################
## Sanity check -- missingness after first regression intrasite matches after 2nd regression passes:
cleanDat.22siteUnreg<-readRDS("4p13.cleanDat.22sites.RDS")
#numericMeta.22sites<-readRDS("4p13.numericMeta.22sites.RDS")
numericMeta.19sites<-numericMeta.22sites[which(!numericMeta.22sites$contributor_code %in% c("U","V","W")),]
cleanDat.19siteUnreg<-cleanDat.22siteUnreg[,match(rownames(numericMeta.19sites),colnames(cleanDat.22siteUnreg))]
dim(cleanDat.19siteUnreg)
#  7335 22547
table(apply(cleanDat.19siteUnreg,2,function(x) length(which(is.na(x)))))
#    0     1  2303  <- # of missing values in column
#19887   406  2254  <- # of columns with that many NA

rm(cleanDat.22siteUnreg)
rm(cleanDat.19siteUnreg)

## Recheck

table(apply(cleanDat.zeroSiteMedian,2,function(x) length(which(is.na(x)))))
#19887   406   2254
table(apply(cleanDat.4p13b2,2,function(x) length(which(is.na(x)))))
#19887   406   2254
table(apply(cleanDat.4p13b3,2,function(x) length(which(is.na(x)))))
#19732   406   2254
table(apply(cleanDat.4p13b4,2,function(x) length(which(is.na(x)))))
#19732   406   2254
table(apply(cleanDat.4p13b5,2,function(x) length(which(is.na(x)))))
#19732   406   2254

# The 155 samples were already dropped from the last 3, and aligned traits also dropping these are paired to them.

## 155 samples already dropped from cleanDat.4p13 b3, b4, b5:
length(which(is.na(numericMeta.reg$age_at_visit) | is.na(numericMeta.reg$sex)))
#155
numericMeta.reg.b345<-numericMeta.reg[which(!is.na(numericMeta.reg$age_at_visit) & !is.na(numericMeta.reg$sex)),]
which(!rownames(regvars.b345)==rownames(numericMeta.reg.b345))
#should be integer(0)

which(!colnames(cleanDat.4p13b3)==rownames(numericMeta.reg.b345))  #<-cleanDat.4p13b3[,match(rownames(numericMeta.reg.b345),colnames(cleanDat.4p13b3))]
which(!colnames(cleanDat.4p13b4)==rownames(numericMeta.reg.b345))  #<-cleanDat.4p13b4[,match(rownames(numericMeta.reg.b345),colnames(cleanDat.4p13b4))]
which(!colnames(cleanDat.4p13b5)==rownames(numericMeta.reg.b345))  #<-cleanDat.4p13b5[,match(rownames(numericMeta.reg.b345),colnames(cleanDat.4p13b5))]
dim(cleanDat.4p13b3)
#  7335 22392
dim(cleanDat.4p13b4)
#  7335 22392
dim(cleanDat.4p13b5)
#  7335 22392


###########################################
## Check e4 carrier top assay correlations

e4.bicor.to.siteCorr.b2.assays<-bicor(t(cleanDat.4p13b2),numericMeta.reg$APOE4.carrier, use='p')
# sort by bicor
e4.bicor.to.siteCorr.b2.assays[order(unlist(t(e4.bicor.to.siteCorr.b2.assays)),decreasing=TRUE),][c(1:5,7331:7335)]
#SPC25|Q9HBM1   LRRN1|Q6UXK5^SL025922@seq.11293.14   CTF1|Q16619^SL002783@seq.13732.79            NEFL|P07196        S100A13|Q99584                TBCA|O75347
#   0.7551699                            0.6467373                           0.4497459      ...    -0.5681812            -0.5724336                 -0.6075718

e4.bicor.to.siteCorr.b3.assays<-bicor(t(cleanDat.4p13b3),numericMeta.reg.b345$APOE4.carrier, use='p')
# sort by bicor
e4.bicor.to.siteCorr.b3.assays[order(unlist(t(e4.bicor.to.siteCorr.b3.assays)),decreasing=TRUE),][c(1:5,7331:7335)]
#SPC25|Q9HBM1   LRRN1|Q6UXK5^SL025922@seq.11293.14   CTF1|Q16619^SL002783@seq.13732.79            NEFL|P07196        S100A13|Q99584                TBCA|O75347
#   0.7568582                            0.6478450                           0.4488890      ...    -0.5650291            -0.5672361                 -0.6022569

e4.bicor.to.siteCorr.b4.assays<-bicor(t(cleanDat.4p13b4),numericMeta.reg.b345$APOE4.carrier, use='p')
# sort by bicor
e4.bicor.to.siteCorr.b4.assays[order(unlist(t(e4.bicor.to.siteCorr.b4.assays)),decreasing=TRUE),][c(1:5,7331:7335)]
#SPC25|Q9HBM1   CTF1|Q16619^SL002783@seq.13732.79                       OTULIN|Q96BN8           FOXO1|Q12778        BCDIN3D|Q725W3   ST8SIA1|Q92185^SL022499@sdeq.21508.7
#   0.5759814                           0.3927323                           0.3174930      ...    -0.4425611            -0.4712035                             -0.4754549

e4.bicor.to.siteCorr.b5.assays<-bicor(t(cleanDat.4p13b5),numericMeta.reg.b345$APOE4.carrier, use='p')
# sort by bicor
e4.bicor.to.siteCorr.b5.assays[order(unlist(t(e4.bicor.to.siteCorr.b5.assays)),decreasing=TRUE),][c(1:5,7331:7335)]
#SPC25|Q9HBM1   CTF1|Q16619^SL002783@seq.13732.79                       OTULIN|Q96BN8           FOXO1|Q12778        BCDIN3D|Q725W3   ST8SIA1|Q92185^SL022499@sdeq.21508.7
#   0.5065277                           0.3925127                           0.3360795      ...    -0.4073110            -0.4247937                             -0.4418332


## Histogram of e4 carrier breakdown within distribution of best correlate: LRRN1/TBCA

LRRN1.med=median(cleanDat.4p13b3["LRRN1|Q6UXK5^SL025922@seq.11293.14",])
TBCA.med=median(cleanDat.4p13b3["TBCA|O75347",])
hist.data=((cleanDat.4p13b3["LRRN1|Q6UXK5^SL025922@seq.11293.14",] - LRRN1.med) - (cleanDat.4p13b3["TBCA|O75347",] - TBCA.med))
hist(hist.data, breaks=100, xlab="log2(abundance ratio):  LRRN1(median-centered) / TBCA(median-centered)", main="APOE e4 Carrier Status Best Nonmissing Correlate")
hist(hist.data[which(numericMeta.reg.b345$APOE4.carrier==1)],breaks=100,col="#FFBBBB99",add=T)  #red overlay
hist(hist.data[which(numericMeta.reg.b345$APOE4.carrier==0)],breaks=100,col="#BBBBFF99",add=T)  #blue
legend("topright",c("E4 Carrier","Non-E4"),fill=c("#FFBBBB99","#BBBBFF99"))
abline(v=0.24,col="maroon",lty=2)

## Originally: impute e4 carrier binary status based on above--for regression and variance partition in b2, b3, b4, and b5
## (we will re-regress b4(b) using this instead of protein proxy for e4)

LRRN1.med.b2=median(cleanDat.4p13b2["LRRN1|Q6UXK5^SL025922@seq.11293.14",])
TBCA.med.b2=median(cleanDat.4p13b2["TBCA|O75347",])
hist.data.b2=((cleanDat.4p13b2["LRRN1|Q6UXK5^SL025922@seq.11293.14",] - LRRN1.med.b2) - (cleanDat.4p13b2["TBCA|O75347",] - TBCA.med.b2))

table(numericMeta.reg$APOE4.carrier)
#   0    1
#9225 5838
numericMeta.reg$APOE4.carrier.imputed<-numericMeta.reg$APOE4.carrier
numericMeta.reg$APOE4.carrier.imputed[which(is.na(numericMeta.reg$APOE4.carrier) & hist.data.b2>=0.25)]<-1
numericMeta.reg$APOE4.carrier.imputed[which(is.na(numericMeta.reg$APOE4.carrier) & hist.data.b2<0.25)]<-0
table(numericMeta.reg$APOE4.carrier.imputed)
#    0     1
#13889  8658
numericMeta.reg.b345$APOE4.carrier.imputed<-numericMeta.reg$APOE4.carrier.imputed[match(rownames(numericMeta.reg.b345),rownames(numericMeta.reg))]
table(numericMeta.reg.b345$APOE4.carrier.imputed)
#    0     1
#13794  8598


## after computing regression pass b4b using imputed e4 carrier binary status, we use ML to impute APOE epsilon full genotypes, and will perform pass b6 regression


#######################LAST REGRESSION - b4b
  regvars.b345$APOE4.carrier<-factor(numericMeta.reg.b345$APOE4.carrier.imputed)
  regvars.b345$APOE4.carrier.int<-as.integer(numericMeta.reg.b345$APOE4.carrier.imputed)
  regvars.b345<-na.omit(regvars.b345)  # removes samples with missing age, sex
  dim(regvars.b345)
  # 22392     8

  ##  Run the regression (4p13 b4b) - Site with Age+Sex+APOE.e4 (APOE4.carrier.imputed) protection

  ## covariate names you want to protect
  cov_keep <- c("Age", "Sex", "APOE4.carrier")   # Technically, we cleared out NAs, but one could list protection covariates that have missing values here

  normExpr.reg <- matrix(NA,nrow=nrow(cleanDat.unreg.b345),ncol=ncol(cleanDat.unreg.b345), dimnames = dimnames(cleanDat.unreg.b345))
  good_samp=which(complete.cases(regvars.b345[, cov_keep]))
  normExpr.reg <-  foreach (i=1:nrow(cleanDat.unreg.b345), .combine=rbind, .packages="stats") %dopar% {  #** .combine=c, .multicombine=TRUE, .packages="stats", .export="regvars") %dopar% {
    y <- as.numeric(cleanDat.unreg.b345[i, ])
    ## keep = samples where y and all covariates are present
    keep <- intersect(which(!is.na(y)), good_samp)             # columns that have a value, and regress with variables that have a value (we know Age and Sex are missing some)

    ## initialise result for this row as all NA
    adj <- rep(NA_real_, length(y))

    fit <- tryCatch(
        lm(y[keep]~Age+Sex+APOE4.carrier +Site, data = regvars.b345[keep, , drop = FALSE]),
        error = function(e) NULL)

    if (!is.null(fit)) {
        coef <- coef(fit)
        ## coefficient[1] + residuals -> adjusted expression
        adj[keep] <- coef["(Intercept)"] + coef["Age"]*regvars.b345[keep,"Age"] + coef[which(grepl("^Sex", names(coef)))]*abs(regvars.b345[keep,"Sex.int"]) + coef["APOE4.carrier1"]*regvars.b345[keep,"APOE4.carrier.int"] + residuals(fit)
    }

#        if (!is.null(fit)) {
#            ## build a model matrix with Site *zeroed out*
#            X     <- model.matrix(fit)
#            X[ , grep("^Site", colnames(X)) ] <- 0     # drop site contribution
#            y_hat <- as.numeric(X %*% coef(fit))       # fitted without Site
#
#            ## protected expression = y_hat + residuals  ( == y - Site effect )
#            adj[keep] <- y_hat + residuals(fit)
#        }
    adj                                    # returned to foreach
  }
  dimnames(normExpr.reg) <- dimnames(cleanDat.unreg.b345)

  cat(paste0("Finished Pass 4p13 b4b regression of intersite variance.\n"))

# sanity check: some sites have missing values:
table(apply(normExpr.reg,2,function(x) length(which(is.na(x)))))
#    0     1  2303  <- # of missing values in column
#19732   406  2254  <- # of columns with that many NA

saveRDS(normExpr.reg,"4p13b4b.normExpr.reg_sites1-19_Fsplit_e4carrierBinaryIMPUTED.RDS")
cleanDat.4p13b4b<-normExpr.reg
#numericMeta.reg.b345 still valid, ordered for this cleanDat


## Sanity check
e4.bicor.to.siteCorr.b4b.assays<-bicor(t(cleanDat.4p13b4b),numericMeta.reg.b345$APOE4.carrier.imputed, use='p')
# sort by bicor
e4.bicor.to.siteCorr.b4b.assays[order(unlist(t(e4.bicor.to.siteCorr.b4b.assays)),decreasing=TRUE),][c(1:5,7331:7335)]
#b4b (current regression)
#SPC25|Q9HBM1   LRRN1|Q6UXK5^SL025922@seq.11293.14   CTF1|Q16619^SL002783@seq.13732.79            NEFL|P07196        S100A13|Q99584                TBCA|O75347
#   0.7228271                            0.6734026                           0.4010664      ...    -0.5761990            -0.6053804                 -0.6564871

#b2 (previously)
#SPC25|Q9HBM1   LRRN1|Q6UXK5^SL025922@seq.11293.14   CTF1|Q16619^SL002783@seq.13732.79            NEFL|P07196        S100A13|Q99584                TBCA|O75347
#   0.7568582                            0.6478450                           0.4488890      ...    -0.5650291            -0.5672361                 -0.6022569


##################################################
# 4p13b2. 19 sites (no site U,V,W; site F split to F1, F2, F3) 1xPAVregr+sitewise Regress, no Protection QC (tSNE plots)

####################
## Examine tSNE of PLASMA - Human only assays (rows) and samples (columns);  no missing data by row.

#Group.3mappedCohortsPlusCTimputed alternate
Group=as.factor(numericMeta.reg$Group.withCTimputed)

numericMeta.plasma<-numericMeta.reg
dim(numericMeta.plasma)
#22547   98
exprMat.plasma<-cleanDat.4p13b2

dim(as.data.frame(na.omit(exprMat.plasma)))
#  5031 22547   # 5031 in H; 7334 in K; 7335 in all others.

tSNE.list.plasma <- Rtsne::Rtsne(t(as.data.frame(na.omit(exprMat.plasma))),perplexity=20)
dim(tSNE.list.plasma$Y)
#22547     2
tSNE.plasma.xy<-as.data.frame(tSNE.list.plasma$Y)
colnames(tSNE.plasma.xy)<-c('x','y')
tSNE.4p13b2.plasma.xy<-tSNE.plasma.xy
rownames(tSNE.4p13b2.plasma.xy)<-rownames(numericMeta.reg)
library(ggplot2)
#library(ggpubr) - rlang upgrade required, ggplot2 upgrade required. cannot install from older source!
library(ggrepel)

# Get the indices of the first occurrence of each unique value in contributor_Fsplit
first_occurrence_indices.4p13b2<-first_occurrence_indices <- match(unique(numericMeta.reg$contributor_Fsplit), numericMeta.reg$contributor_Fsplit)

labels.4p13b2<-labels<-numericMeta.reg$contributor_Fsplit[first_occurrence_indices.4p13b2]

tSNE.plasma.samples.sites<-ggplot2::ggplot(tSNE.4p13b2.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.reg$contributor_Fsplit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13b2.plasma.xy[first_occurrence_indices.4p13b2, ],
                  aes(x=x,y=y, label = labels.4p13b2), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.age<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$age_at_visit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleMatrix<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_matrix), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleType<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_type), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.apoe4dose<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$APOE4.Dose), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.MMSE<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$MMSE), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.Group<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=Group), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )


pdf(file="4p13b2.2PAVproteinIntrasiteRegressed+Regress19sites_Fsplit.tSNE-Plasma(7335x22547)-samples_coloredByTraits.pdf",width=11,height=9)
  print(tSNE.plasma.samples.sites)
  print(tSNE.plasma.samples.age)
  print(tSNE.plasma.samples.sampleMatrix)
  print(tSNE.plasma.samples.sampleType)
  print(tSNE.plasma.samples.apoe4dose)
  print(tSNE.plasma.samples.MMSE)
  print(tSNE.plasma.samples.Group)
dev.off()


tSNE.4p13b2.plasma.samples.sites<-tSNE.plasma.samples.sites


##################################################
# 4p13b3. 19 sites (no site U,V,W; site F split to F1, F2, F3) 1xPAVregr+sitewise Regress, Protect Age+Sex QC (tSNE plots)

####################
## Examine tSNE of PLASMA - Human only assays (rows) and samples (columns);  no missing data by row.

#Group.3mappedCohortsPlusCTimputed)
Group=as.factor(numericMeta.reg.b345$Group.withCTimputed)

numericMeta.plasma<-numericMeta.reg.b345
dim(numericMeta.plasma)
#22392   98   #previously (16 sites with Calibrators):  18739   93
exprMat.plasma<-cleanDat.4p13b3

dim(as.data.frame(na.omit(exprMat.plasma)))
#  5031 22392   # 5031 in H; 7334 in K; 7335 in all others.

tSNE.list.plasma <- Rtsne::Rtsne(t(as.data.frame(na.omit(exprMat.plasma))),perplexity=20)
dim(tSNE.list.plasma$Y)
#22392     2
tSNE.plasma.xy<-as.data.frame(tSNE.list.plasma$Y)
colnames(tSNE.plasma.xy)<-c('x','y')
tSNE.4p13b3.plasma.xy<-tSNE.plasma.xy
rownames(tSNE.4p13b3.plasma.xy)<-rownames(numericMeta.reg.b345)
library(ggplot2)
#library(ggpubr) - rlang upgrade required, ggplot2 upgrade required. cannot install from older source!
library(ggrepel)

# Get the indices of the first occurrence of each unique value in contributor_Fsplit
first_occurrence_indices.4p13b3<-first_occurrence_indices <- match(unique(numericMeta.reg.b345$contributor_Fsplit), numericMeta.reg.b345$contributor_Fsplit)

labels.4p13b3<-labels<-numericMeta.reg.b345$contributor_Fsplit[first_occurrence_indices.4p13b3]

# ------------------------------------------------------------------------
# ANNOTATION: Visualize second-pass regression candidates by t-SNE across
# biological and technical covariates.
# ------------------------------------------------------------------------

tSNE.plasma.samples.sites<-ggplot2::ggplot(tSNE.4p13b3.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.reg.b345$contributor_Fsplit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13b3.plasma.xy[first_occurrence_indices.4p13b3, ],
                  aes(x=x,y=y, label = labels.4p13b3), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.age<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$age_at_visit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleMatrix<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_matrix), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleType<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_type), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.apoe4dose<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$APOE4.Dose), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.MMSE<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$MMSE), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.Group<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=Group), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )


pdf(file="4p13b3.2PAVproteinIntrasiteRegressed+Regress19sites_protectAgeSex_Fsplit.tSNE-Plasma(7335x22392)-samples_coloredByTraits.pdf",width=11,height=9)
  print(tSNE.plasma.samples.sites)
  print(tSNE.plasma.samples.age)
  print(tSNE.plasma.samples.sampleMatrix)
  print(tSNE.plasma.samples.sampleType)
  print(tSNE.plasma.samples.apoe4dose)
  print(tSNE.plasma.samples.MMSE)
  print(tSNE.plasma.samples.Group)
dev.off()


tSNE.4p13b3.plasma.samples.sites<-tSNE.plasma.samples.sites


##################################################
# 4p13b4. 19 sites (no site U,V,W; site F split to F1, F2, F3) 1xPAVregr+sitewise Regress, Protect Age+Sex+APOE.e4 proxy SPC25 QC (tSNE plots)

####################
## Examine tSNE of PLASMA - Human only assays (rows) and samples (columns);  no missing data by row.

numericMeta.plasma<-numericMeta.reg.b345
dim(numericMeta.plasma)
#22392   98
exprMat.plasma<-cleanDat.4p13b4

dim(as.data.frame(na.omit(exprMat.plasma)))
#  5031 22392   # 5031 in H; 7334 in K; 7335 in all others.
#exprMat.plasma<-exprMat.plasma[,which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))]
#numericMeta.plasma<-numericMeta.plasma[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma))))),]

#Group.3mappedCohortsPlusCTimputed alternate
Group=as.factor(numericMeta.reg.b345$Group.withCTimputed ) #[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))])


tSNE.list.plasma <- Rtsne::Rtsne(t(as.data.frame(na.omit(exprMat.plasma))),perplexity=20)
dim(tSNE.list.plasma$Y)
#22392     2
tSNE.plasma.xy<-as.data.frame(tSNE.list.plasma$Y)
colnames(tSNE.plasma.xy)<-c('x','y')
tSNE.4p13b4.plasma.xy<-tSNE.plasma.xy
rownames(tSNE.4p13b4.plasma.xy)<-rownames(numericMeta.reg.b345) #[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))]
library(ggplot2)
#library(ggpubr) - rlang upgrade required, ggplot2 upgrade required. cannot install from older source!
library(ggrepel)

# Get the indices of the first occurrence of each unique value in contributor_Fsplit
first_occurrence_indices.4p13b4<-first_occurrence_indices <- match(unique(numericMeta.reg.b345$contributor_Fsplit), numericMeta.reg.b345$contributor_Fsplit)

labels.4p13b4<-labels<-numericMeta.reg.b345$contributor_Fsplit[first_occurrence_indices.4p13b4]

tSNE.plasma.samples.sites<-ggplot2::ggplot(tSNE.4p13b4.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.reg.b345$contributor_Fsplit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13b4.plasma.xy[first_occurrence_indices.4p13b4, ],
                  aes(x=x,y=y, label = labels.4p13b4), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.age<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$age_at_visit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleMatrix<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_matrix), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleType<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_type), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.apoe4dose<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$APOE4.Dose), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.MMSE<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$MMSE), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.Group<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=Group), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )


pdf(file="4p13b4.2PAVproteinIntrasiteRegressed+Regress19sites_protectAgeSexAPOE.e4_LRRN1_Fsplit.tSNE-Plasma(7335x22392)-samples_coloredByTraits.pdf",width=11,height=9)
  print(tSNE.plasma.samples.sites + labs(colour="Site"))
  print(tSNE.plasma.samples.age)
  print(tSNE.plasma.samples.sampleMatrix)
  print(tSNE.plasma.samples.sampleType)
  print(tSNE.plasma.samples.apoe4dose)
  print(tSNE.plasma.samples.MMSE)
  print(tSNE.plasma.samples.Group + labs(colour="Group"))
dev.off()


tSNE.4p13b4.plasma.samples.sites<-tSNE.plasma.samples.sites


##################################################
# 4p13b5. 19 sites (no site U,V,W; site F split to F1, F2, F3) 1xPAVregr+sitewise Regress, Protect Age+Sex+APOE.e4 proxy SPC25.TBCAlog2ratio QC (tSNE plots)

####################
## Examine tSNE of PLASMA - Human only assays (rows) and samples (columns);  no missing data by row.

#Group.3mappedCohortsPlusCTimputed alternate
Group=as.factor(numericMeta.reg.b345$Group.withCTimputed)

numericMeta.plasma<-numericMeta.reg.b345
dim(numericMeta.plasma)
#22392   98
exprMat.plasma<-cleanDat.4p13b5

dim(as.data.frame(na.omit(exprMat.plasma)))
#  5031 22392   # 5031 in H; 7334 in K; 7335 in all others.
exprMat.plasma<-exprMat.plasma[,which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))]
numericMeta.plasma<-numericMeta.plasma  #[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma))))),]

tSNE.list.plasma <- Rtsne::Rtsne(t(as.data.frame(na.omit(exprMat.plasma))),perplexity=20)  #)))[!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))),]
dim(tSNE.list.plasma$Y)
#22392     2
tSNE.plasma.xy<-as.data.frame(tSNE.list.plasma$Y)
colnames(tSNE.plasma.xy)<-c('x','y')
tSNE.4p13b5.plasma.xy<-tSNE.plasma.xy
rownames(tSNE.4p13b5.plasma.xy)<-rownames(numericMeta.reg.b345) #[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))]
library(ggplot2)
#library(ggpubr) - rlang upgrade required, ggplot2 upgrade required. cannot install from older source!
library(ggrepel)

# Get the indices of the first occurrence of each unique value in contributor_Fsplit
first_occurrence_indices.4p13b5<-first_occurrence_indices <- match(unique(numericMeta.reg.b345$contributor_Fsplit), numericMeta.reg.b345$contributor_Fsplit)

labels.4p13b5<-labels<-numericMeta.reg.b345$contributor_Fsplit[first_occurrence_indices.4p13b5]

tSNE.plasma.samples.sites<-ggplot2::ggplot(tSNE.4p13b5.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.reg.b345$contributor_Fsplit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13b5.plasma.xy[first_occurrence_indices.4p13b5, ],
                  aes(x=x,y=y, label = labels.4p13b5), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.age<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$age_at_visit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleMatrix<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_matrix), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleType<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_type), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.apoe4dose<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$APOE4.Dose), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.MMSE<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$MMSE), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.Group<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=Group), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )


pdf(file="4p13b5.2PAVproteinIntrasiteRegressed+Regress19sites_protectAgeSexAPOE.e4_LRRN1.TBCAlog2ratio_Fsplit.tSNE-Plasma(7335x22392)-samples_coloredByTraits.pdf",width=11,height=9)
  print(tSNE.plasma.samples.sites + labs(colour="Site"))
  print(tSNE.plasma.samples.age)
  print(tSNE.plasma.samples.sampleMatrix)
  print(tSNE.plasma.samples.sampleType)
  print(tSNE.plasma.samples.apoe4dose)
  print(tSNE.plasma.samples.MMSE)
  print(tSNE.plasma.samples.Group)
dev.off()


tSNE.4p13b5.plasma.samples.sites<-tSNE.plasma.samples.sites


##################################################
# 4p13b4b. 19 sites (no site U,V,W; site F split to F1, F2, F3) 1xPAVregr+sitewise Regress, Protect Age+Sex+APOE.e4 carrier status (binary, imputed NA) QC (tSNE plots)

####################
## Examine tSNE of PLASMA - Human only assays (rows) and samples (columns);  no missing data by row.

numericMeta.plasma<-numericMeta.reg.b345
dim(numericMeta.plasma)
#22392   99
exprMat.plasma<-cleanDat.4p13b4b

dim(as.data.frame(na.omit(exprMat.plasma)))
#  5031 22392   # 5031 in H; 7334 in K; 7335 in all others.
#exprMat.plasma<-exprMat.plasma[,which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))]
#numericMeta.plasma<-numericMeta.plasma[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma))))),]

#Group.3mappedCohortsPlusCTimputed alternate
Group=as.factor(numericMeta.reg.b345$Group.withCTimputed ) #[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))])


tSNE.list.plasma <- Rtsne::Rtsne(t(as.data.frame(na.omit(exprMat.plasma))),perplexity=20)
dim(tSNE.list.plasma$Y)
#22392     2
tSNE.plasma.xy<-as.data.frame(tSNE.list.plasma$Y)
colnames(tSNE.plasma.xy)<-c('x','y')
tSNE.4p13b4b.plasma.xy<-tSNE.plasma.xy
rownames(tSNE.4p13b4b.plasma.xy)<-rownames(numericMeta.reg.b345) #[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))]
library(ggplot2)
#library(ggpubr) - rlang upgrade required, ggplot2 upgrade required. cannot install from older source!
library(ggrepel)

# Get the indices of the first occurrence of each unique value in contributor_Fsplit
first_occurrence_indices.4p13b4b<-first_occurrence_indices <- match(unique(numericMeta.reg.b345$contributor_Fsplit), numericMeta.reg.b345$contributor_Fsplit)

labels.4p13b4b<-labels<-numericMeta.reg.b345$contributor_Fsplit[first_occurrence_indices.4p13b4b]

tSNE.plasma.samples.sites<-ggplot2::ggplot(tSNE.4p13b4b.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.reg.b345$contributor_Fsplit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13b4b.plasma.xy[first_occurrence_indices.4p13b4b, ],
                  aes(x=x,y=y, label = labels.4p13b4b), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.age<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$age_at_visit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleMatrix<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_matrix), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleType<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_type), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.apoe4dose<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$APOE4.Dose), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.MMSE<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$MMSE), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.Group<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=Group), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )


pdf(file="4p13b4b.2PAVproteinIntrasiteRegressed+Regress19sites_protectAgeSexAPOE.e4_carrierBinaryImputed_Fsplit.tSNE-Plasma(7335x22392)-samples_coloredByTraits.pdf",width=11,height=9)
  print(tSNE.plasma.samples.sites + labs(colour="Site"))
  print(tSNE.plasma.samples.age)
  print(tSNE.plasma.samples.sampleMatrix)
  print(tSNE.plasma.samples.sampleType)
  print(tSNE.plasma.samples.apoe4dose)
  print(tSNE.plasma.samples.MMSE)
  print(tSNE.plasma.samples.Group + labs(colour="Group"))
dev.off()


tSNE.4p13b4b.plasma.samples.sites<-tSNE.plasma.samples.sites


## Set up parallel backend
library("doParallel")
parallelThreads=31  #now Windows02  #max is number of processes that can run on your computer at one time
stopCluster(clusterLocal)
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="PSOCK")

registerDoParallel(clusterLocal)


#############################################################################
## 4p13b2. Variance Partition regressed (2PAV regression intrasite)+Site regressed; no protection (QC)

regvars.vp<-data.frame(numericMeta.reg)
regvars.vp$Sex<-factor(abs(regvars.vp$sex -2))
regvars.vp$Age<-as.numeric(regvars.vp$age_at_visit)
regvars.vp$contributor_Fsplit<-factor(regvars.vp$contributor_Fsplit)
regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1)
regvars.vp$RegrBloodPreanalyticFactor.HBZ<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HBZ)
#regvars.vp$APOE.E4carrier.Proxy.LRRN1<-as.numeric(cleanDat.4p13b2["LRRN1|Q6UXK5^SL025922@seq.11293.14",])
regvars.vp$APOE.E4carrier.imputed<-factor(numericMeta.reg$APOE4.carrier.imputed)

# too many missing values:
#form <- ~ age_at_visit+(1|Sex)+(1|raceAA)+(1|recruited_control)+(1|ad)+(1|ftd)+(1|pd)+(1|als)+(1|mci_sci)+(1|mi)+(1|C9Orf72)+(1|GRN)+(1|MAPT)+APOE4.Dose+(1|contributor_code)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+MMSE
#form <- ~ (1|contributor_code) +RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ
form <- ~ Age+(1|Sex)+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+(1|APOE.E4carrier.imputed)  #+APOE.E4carrier.Proxy.LRRN1

library(variancePartition)

# (If regressed) REMOVE regressed proteins to avoid "Response variable 4641 has a variance of 0":
#regrProts.idx<-c(which(grepl("^HNRNPA2B1\\|",rownames(cleanDat))),which(grepl("^HBZ\\|",rownames(cleanDat))))
#regrProts.idx
#5569 6087 -- take out
#previously: 6156 4641
#integer(0)  # if already removed

#previously: removed 2 regr proteins in line:  varPart.reg <- fitExtractVarPartModel(impute::impute.knn(cleanDat[-c(4641,6156),])$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = 8, type = "SOCK"))
varPart.b2 <- fitExtractVarPartModel(impute::impute.knn(as.matrix(cleanDat.4p13b2[,]))$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp.b2 <- sortCols(varPart.b2,FUN=median,last= c("Residuals"))

pdf(file="4p13b2.19sites1xPAVregr+SiteRegress_noProtect-VariancePartition-PLASMA-7335x22547.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp.b2, main="HDS 1.3ms - 4p13b2 - 19 sites 2x Regr(2PAV) + Site Regr, No Protect" )

	SexSortOrder<-order(vp.b2[["Sex"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][SexSortOrder]; }
	rownames(vp.b2)<-rownames(vp.b2)[SexSortOrder]

	plotPercentBars( vp.b2[1:50,]) + ggtitle( "Top Sex-covariates" )


	AgeSortOrder<-order(vp.b2[["Age"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][AgeSortOrder]; }
	rownames(vp.b2)<-rownames(vp.b2)[AgeSortOrder]

	plotPercentBars( vp.b2[1:50,]) + ggtitle( "Top Age-covariates" )


        BatchSortOrder<-order(vp.b2[["RegrBloodPreanalyticFactor.HNRNPA2B1"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][BatchSortOrder]; }
        rownames(vp.b2)<-rownames(vp.b2)[BatchSortOrder]

        plotPercentBars( vp.b2[1:50,]) + ggtitle( "Top Time To Spin (HNRNPA2B1)-covariates" )


        BatchSortOrder<-order(vp.b2[["RegrBloodPreanalyticFactor.HBZ"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][BatchSortOrder]; }
        rownames(vp.b2)<-rownames(vp.b2)[BatchSortOrder]

        plotPercentBars( vp.b2[1:50,]) + ggtitle( "Top Preanalytical Factor 2 (HBZ)-covariates" )


        BatchSortOrder<-order(vp.b2[["contributor_Fsplit"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][BatchSortOrder]; }
        rownames(vp.b2)<-rownames(vp.b2)[BatchSortOrder]

        plotPercentBars( vp.b2[1:50,]) + ggtitle( "Top Contributor (site)-covariates" )


        BatchSortOrder<-order(vp.b2[["APOE.E4carrier.imputed"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][BatchSortOrder]; }
        rownames(vp.b2)<-rownames(vp.b2)[BatchSortOrder]

        plotPercentBars( vp.b2[1:50,]) + ggtitle( "Top APOE e4 carrier (Binary, imputed NA)-covariates" )


#	BatchSortOrder<-order(vp.b2[["MMSE"]],decreasing=TRUE)
#	#rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
#	for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][BatchSortOrder]; }
#	rownames(vp.b2)<-rownames(vp.b2)[BatchSortOrder]
#
#	plotPercentBars( vp.b2[1:50,]) + ggtitle( "MMSE Cog. Score-covariates" )

dev.off()

varPart.b2<-vp.b2
saveRDS(varPart.b2,"4p13b2.19sites1xPAVregr+SiteRegress_noProtect.varPart.b2.RDS")

# ------------------------------------------------------------------------
# ANNOTATION: Compare second-pass candidate regressions by variance
# partitioning.
# ------------------------------------------------------------------------


#############################################################################
## 4p13b3. Variance Partition regressed (2PAV regression intrasite)+Site regressed; protect Age+Sex (QC)

regvars.vp<-data.frame(numericMeta.reg.b345)
regvars.vp$Sex<-factor(abs(regvars.vp$sex -2))
regvars.vp$Age<-as.numeric(regvars.vp$age_at_visit)
regvars.vp$contributor_Fsplit<-factor(regvars.vp$contributor_Fsplit)
regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1)
regvars.vp$RegrBloodPreanalyticFactor.HBZ<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HBZ)
#regvars.vp$APOE.E4carrier.Proxy.LRRN1<-as.numeric(cleanDat.4p13b2["LRRN1|Q6UXK5^SL025922@seq.11293.14",])
regvars.vp$APOE.E4carrier.imputed<-factor(regvars.vp$APOE4.carrier.imputed)

# too many missing values:
#form <- ~ age_at_visit+(1|Sex)+(1|raceAA)+(1|recruited_control)+(1|ad)+(1|ftd)+(1|pd)+(1|als)+(1|mci_sci)+(1|mi)+(1|C9Orf72)+(1|GRN)+(1|MAPT)+APOE4.Dose+(1|contributor_code)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+MMSE
#form <- ~ (1|contributor_code) +RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ
form <- ~ Age+(1|Sex)+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+(1|APOE.E4carrier.imputed)  #+APOE.E4carrier.Proxy.LRRN1

library(variancePartition)

# (If regressed) REMOVE regressed proteins to avoid "Response variable 4641 has a variance of 0":
#regrProts.idx<-c(which(grepl("^HNRNPA2B1\\|",rownames(cleanDat))),which(grepl("^HBZ\\|",rownames(cleanDat))))
#regrProts.idx
#5569 6087 -- take out
#previously: 6156 4641
#integer(0)  # if already removed

#previously: removed 2 regr proteins in line:  varPart.reg <- fitExtractVarPartModel(impute::impute.knn(cleanDat[-c(4641,6156),])$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = 8, type = "SOCK"))
varPart.b3 <- fitExtractVarPartModel(impute::impute.knn(as.matrix(cleanDat.4p13b3[,]))$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp.b3 <- sortCols(varPart.b3,FUN=median,last= c("Residuals"))

pdf(file="4p13b3.19sites1xPAVregr+SiteRegress_Protect_age+sex-VariancePartition-PLASMA-7335x22392.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp.b3, main="HDS 1.3ms - 4p13b3 - 19 sites 2x Regr(2PAV) + Site Regr Prot. age+sex" )

	SexSortOrder<-order(vp.b3[["Sex"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b3)) { vp.b3[[i]]<-vp.b3[[i]][SexSortOrder]; }
	rownames(vp.b3)<-rownames(vp.b3)[SexSortOrder]

	plotPercentBars( vp.b3[1:50,]) + ggtitle( "Top Sex-covariates" )


	AgeSortOrder<-order(vp.b3[["Age"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b3)) { vp.b3[[i]]<-vp.b3[[i]][AgeSortOrder]; }
	rownames(vp.b3)<-rownames(vp.b3)[AgeSortOrder]

	plotPercentBars( vp.b3[1:50,]) + ggtitle( "Top Age-covariates" )


        BatchSortOrder<-order(vp.b3[["RegrBloodPreanalyticFactor.HNRNPA2B1"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b3)) { vp.b3[[i]]<-vp.b3[[i]][BatchSortOrder]; }
        rownames(vp.b3)<-rownames(vp.b3)[BatchSortOrder]

        plotPercentBars( vp.b3[1:50,]) + ggtitle( "Top Time To Spin (HNRNPA2B1)-covariates" )


        BatchSortOrder<-order(vp.b3[["RegrBloodPreanalyticFactor.HBZ"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b3)) { vp.b3[[i]]<-vp.b3[[i]][BatchSortOrder]; }
        rownames(vp.b3)<-rownames(vp.b3)[BatchSortOrder]

        plotPercentBars( vp.b3[1:50,]) + ggtitle( "Top Preanalytical Factor 2 (HBZ)-covariates" )


        BatchSortOrder<-order(vp.b3[["contributor_Fsplit"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b3)) { vp.b3[[i]]<-vp.b3[[i]][BatchSortOrder]; }
        rownames(vp.b3)<-rownames(vp.b3)[BatchSortOrder]

        plotPercentBars( vp.b3[1:50,]) + ggtitle( "Top Contributor (site)-covariates" )


        BatchSortOrder<-order(vp.b3[["APOE.E4carrier.imputed"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b3)) { vp.b3[[i]]<-vp.b3[[i]][BatchSortOrder]; }
        rownames(vp.b3)<-rownames(vp.b3)[BatchSortOrder]

        plotPercentBars( vp.b3[1:50,]) + ggtitle( "Top APOE e4 carrier (Binary, imputed NA)-covariates" )


#	BatchSortOrder<-order(vp.b3[["MMSE"]],decreasing=TRUE)
#	#rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
#	for (i in ls(vp.b3)) { vp.b3[[i]]<-vp.b3[[i]][BatchSortOrder]; }
#	rownames(vp.b3)<-rownames(vp.b3)[BatchSortOrder]
#
#	plotPercentBars( vp.b3[1:50,]) + ggtitle( "MMSE Cog. Score-covariates" )

dev.off()

varPart.b3<-vp.b3
saveRDS(varPart.b3,"4p13b3.19sites1xPAVregr+SiteRegress_Protect_age+sex.varPart.b3.RDS")


#############################################################################
## 4p13b4. Variance Partition regressed (2PAV regression intrasite)+Site regressed; protect Age+Sex+LRRN1 (e4 carrier proxy) (QC)

regvars.vp<-data.frame(numericMeta.reg.b345)
regvars.vp$Sex<-factor(abs(regvars.vp$sex -2))
regvars.vp$Age<-as.numeric(regvars.vp$age_at_visit)
regvars.vp$contributor_Fsplit<-factor(regvars.vp$contributor_Fsplit)
regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1)
regvars.vp$RegrBloodPreanalyticFactor.HBZ<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HBZ)
#regvars.vp$APOE.E4carrier.Proxy.LRRN1<-as.numeric(cleanDat.4p13b2["LRRN1|Q6UXK5^SL025922@seq.11293.14",])
regvars.vp$APOE.E4carrier.imputed<-factor(regvars.vp$APOE4.carrier.imputed)

# too many missing values:
#form <- ~ age_at_visit+(1|Sex)+(1|raceAA)+(1|recruited_control)+(1|ad)+(1|ftd)+(1|pd)+(1|als)+(1|mci_sci)+(1|mi)+(1|C9Orf72)+(1|GRN)+(1|MAPT)+APOE4.Dose+(1|contributor_code)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+MMSE
#form <- ~ (1|contributor_code) +RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ
form <- ~ Age+(1|Sex)+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+(1|APOE.E4carrier.imputed)  #+APOE.E4carrier.Proxy.LRRN1

library(variancePartition)

# (If regressed) REMOVE regressed proteins to avoid "Response variable 4641 has a variance of 0":
#regrProts.idx<-c(which(grepl("^HNRNPA2B1\\|",rownames(cleanDat))),which(grepl("^HBZ\\|",rownames(cleanDat))))
#regrProts.idx
#5569 6087 -- take out
#previously: 6156 4641
#integer(0)  # if already removed

#previously: removed 2 regr proteins in line:  varPart.reg <- fitExtractVarPartModel(impute::impute.knn(cleanDat[-c(4641,6156),])$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = 8, type = "SOCK"))
varPart.b4 <- fitExtractVarPartModel(impute::impute.knn(as.matrix(cleanDat.4p13b4[,]))$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp.b4 <- sortCols(varPart.b4,FUN=median,last= c("Residuals"))

pdf(file="4p13b4.19sites1xPAVregr+SiteRegress_Protect_age+sex+LRRN1-VariancePartition-PLASMA-7335x22392.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp.b4, main="HDS 1.3ms - 4p13b4 - 19 sites 2x Regr(2PAV) + Site Regr Prot. age+sex+LRRN1" )

	SexSortOrder<-order(vp.b4[["Sex"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b4)) { vp.b4[[i]]<-vp.b4[[i]][SexSortOrder]; }
	rownames(vp.b4)<-rownames(vp.b4)[SexSortOrder]

	plotPercentBars( vp.b4[1:50,]) + ggtitle( "Top Sex-covariates" )


	AgeSortOrder<-order(vp.b4[["Age"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b4)) { vp.b4[[i]]<-vp.b4[[i]][AgeSortOrder]; }
	rownames(vp.b4)<-rownames(vp.b4)[AgeSortOrder]

	plotPercentBars( vp.b4[1:50,]) + ggtitle( "Top Age-covariates" )


        BatchSortOrder<-order(vp.b4[["RegrBloodPreanalyticFactor.HNRNPA2B1"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4)) { vp.b4[[i]]<-vp.b4[[i]][BatchSortOrder]; }
        rownames(vp.b4)<-rownames(vp.b4)[BatchSortOrder]

        plotPercentBars( vp.b4[1:50,]) + ggtitle( "Top Time To Spin (HNRNPA2B1)-covariates" )


        BatchSortOrder<-order(vp.b4[["RegrBloodPreanalyticFactor.HBZ"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4)) { vp.b4[[i]]<-vp.b4[[i]][BatchSortOrder]; }
        rownames(vp.b4)<-rownames(vp.b4)[BatchSortOrder]

        plotPercentBars( vp.b4[1:50,]) + ggtitle( "Top Preanalytical Factor 2 (HBZ)-covariates" )


        BatchSortOrder<-order(vp.b4[["contributor_Fsplit"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4)) { vp.b4[[i]]<-vp.b4[[i]][BatchSortOrder]; }
        rownames(vp.b4)<-rownames(vp.b4)[BatchSortOrder]

        plotPercentBars( vp.b4[1:50,]) + ggtitle( "Top Contributor (site)-covariates" )


        BatchSortOrder<-order(vp.b4[["APOE.E4carrier.imputed"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4)) { vp.b4[[i]]<-vp.b4[[i]][BatchSortOrder]; }
        rownames(vp.b4)<-rownames(vp.b4)[BatchSortOrder]

        plotPercentBars( vp.b4[1:50,]) + ggtitle( "Top APOE e4 carrier (Binary, imputed NA)-covariates" )


#	BatchSortOrder<-order(vp.b4[["MMSE"]],decreasing=TRUE)
#	#rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
#	for (i in ls(vp.b4)) { vp.b4[[i]]<-vp.b4[[i]][BatchSortOrder]; }
#	rownames(vp.b4)<-rownames(vp.b4)[BatchSortOrder]
#
#	plotPercentBars( vp.b4[1:50,]) + ggtitle( "MMSE Cog. Score-covariates" )

dev.off()

varPart.b4<-vp.b4
saveRDS(varPart.b4,"4p13b4.19sites1xPAVregr+SiteRegress_Protect_age+sex+LRRN1.varPart.b4.RDS")


#############################################################################
## 4p13b5. Variance Partition regressed (2PAV regression intrasite)+Site regressed; protect Age+Sex+log2(LRRN1/TBCA) (e4 carrier proxy) (QC)

regvars.vp<-data.frame(numericMeta.reg.b345)
regvars.vp$Sex<-factor(abs(regvars.vp$sex -2))
regvars.vp$Age<-as.numeric(regvars.vp$age_at_visit)
regvars.vp$contributor_Fsplit<-factor(regvars.vp$contributor_Fsplit)
regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1)
regvars.vp$RegrBloodPreanalyticFactor.HBZ<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HBZ)
#regvars.vp$APOE.E4carrier.Proxy.LRRN1<-as.numeric(cleanDat.4p13b2["LRRN1|Q6UXK5^SL025922@seq.11293.14",])
regvars.vp$APOE.E4carrier.imputed<-factor(regvars.vp$APOE4.carrier.imputed)

# too many missing values:
#form <- ~ age_at_visit+(1|Sex)+(1|raceAA)+(1|recruited_control)+(1|ad)+(1|ftd)+(1|pd)+(1|als)+(1|mci_sci)+(1|mi)+(1|C9Orf72)+(1|GRN)+(1|MAPT)+APOE4.Dose+(1|contributor_code)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+MMSE
#form <- ~ (1|contributor_code) +RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ
form <- ~ Age+(1|Sex)+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+(1|APOE.E4carrier.imputed)  #+APOE.E4carrier.Proxy.LRRN1

library(variancePartition)

# (If regressed) REMOVE regressed proteins to avoid "Response variable 4641 has a variance of 0":
#regrProts.idx<-c(which(grepl("^HNRNPA2B1\\|",rownames(cleanDat))),which(grepl("^HBZ\\|",rownames(cleanDat))))
#regrProts.idx
#5569 6087 -- take out
#previously: 6156 4641
#integer(0)  # if already removed

#previously: removed 2 regr proteins in line:  varPart.reg <- fitExtractVarPartModel(impute::impute.knn(cleanDat[-c(4641,6156),])$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = 8, type = "SOCK"))
varPart.b5 <- fitExtractVarPartModel(impute::impute.knn(as.matrix(cleanDat.4p13b5[,]))$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp.b5 <- sortCols(varPart.b5,FUN=median,last= c("Residuals"))

pdf(file="4p13b5.19sites1xPAVregr+SiteRegress_Protect_age+sex+LRRN1minusTBCA-VariancePartition-PLASMA-7335x22392.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp.b5, main="HDS 1.3ms - 4p13b5 - 19 sites 2x Regr(2PAV) + Site Regr Prot. age+sex+(LRRN1-TBCA)" )

	SexSortOrder<-order(vp.b5[["Sex"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b5)) { vp.b5[[i]]<-vp.b5[[i]][SexSortOrder]; }
	rownames(vp.b5)<-rownames(vp.b5)[SexSortOrder]

	plotPercentBars( vp.b5[1:50,]) + ggtitle( "Top Sex-covariates" )


	AgeSortOrder<-order(vp.b5[["Age"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b5)) { vp.b5[[i]]<-vp.b5[[i]][AgeSortOrder]; }
	rownames(vp.b5)<-rownames(vp.b5)[AgeSortOrder]

	plotPercentBars( vp.b5[1:50,]) + ggtitle( "Top Age-covariates" )


        BatchSortOrder<-order(vp.b5[["RegrBloodPreanalyticFactor.HNRNPA2B1"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b5)) { vp.b5[[i]]<-vp.b5[[i]][BatchSortOrder]; }
        rownames(vp.b5)<-rownames(vp.b5)[BatchSortOrder]

        plotPercentBars( vp.b5[1:50,]) + ggtitle( "Top Time To Spin (HNRNPA2B1)-covariates" )


        BatchSortOrder<-order(vp.b5[["RegrBloodPreanalyticFactor.HBZ"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b5)) { vp.b5[[i]]<-vp.b5[[i]][BatchSortOrder]; }
        rownames(vp.b5)<-rownames(vp.b5)[BatchSortOrder]

        plotPercentBars( vp.b5[1:50,]) + ggtitle( "Top Preanalytical Factor 2 (HBZ)-covariates" )


        BatchSortOrder<-order(vp.b5[["contributor_Fsplit"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b5)) { vp.b5[[i]]<-vp.b5[[i]][BatchSortOrder]; }
        rownames(vp.b5)<-rownames(vp.b5)[BatchSortOrder]

        plotPercentBars( vp.b5[1:50,]) + ggtitle( "Top Contributor (site)-covariates" )


        BatchSortOrder<-order(vp.b5[["APOE.E4carrier.imputed"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b5)) { vp.b5[[i]]<-vp.b5[[i]][BatchSortOrder]; }
        rownames(vp.b5)<-rownames(vp.b5)[BatchSortOrder]

        plotPercentBars( vp.b5[1:50,]) + ggtitle( "Top APOE e4 carrier (Binary, imputed NA)-covariates" )


#	BatchSortOrder<-order(vp.b5[["MMSE"]],decreasing=TRUE)
#	#rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
#	for (i in ls(vp.b5)) { vp.b5[[i]]<-vp.b5[[i]][BatchSortOrder]; }
#	rownames(vp.b5)<-rownames(vp.b5)[BatchSortOrder]
#
#	plotPercentBars( vp.b5[1:50,]) + ggtitle( "MMSE Cog. Score-covariates" )

dev.off()

varPart.b5<-vp.b5
saveRDS(varPart.b5,"4p13b5.19sites1xPAVregr+SiteRegress_Protect_age+sex+LRRN1minusTBCA.varPart.b5.RDS")


#############################################################################
## 4p13b4b. Variance Partition regressed (2PAV regression intrasite)+Site regressed; protect Age+Sex+E4 carrier binary status (NA imputed) (QC)

regvars.vp<-data.frame(numericMeta.reg.b345)
regvars.vp$Sex<-factor(abs(regvars.vp$sex -2))
regvars.vp$Age<-as.numeric(regvars.vp$age_at_visit)
regvars.vp$contributor_Fsplit<-factor(regvars.vp$contributor_Fsplit)
regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1)
regvars.vp$RegrBloodPreanalyticFactor.HBZ<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HBZ)
#regvars.vp$APOE.E4carrier.Proxy.LRRN1<-as.numeric(cleanDat.4p13b2["LRRN1|Q6UXK5^SL025922@seq.11293.14",])
regvars.vp$APOE.E4carrier.imputed<-factor(regvars.vp$APOE4.carrier.imputed)

# too many missing values:
#form <- ~ age_at_visit+(1|Sex)+(1|raceAA)+(1|recruited_control)+(1|ad)+(1|ftd)+(1|pd)+(1|als)+(1|mci_sci)+(1|mi)+(1|C9Orf72)+(1|GRN)+(1|MAPT)+APOE4.Dose+(1|contributor_code)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+MMSE
#form <- ~ (1|contributor_code) +RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ
form <- ~ Age+(1|Sex)+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+(1|APOE.E4carrier.imputed)  #+APOE.E4carrier.Proxy.LRRN1

library(variancePartition)

# (If regressed) REMOVE regressed proteins to avoid "Response variable 4641 has a variance of 0":
#regrProts.idx<-c(which(grepl("^HNRNPA2B1\\|",rownames(cleanDat))),which(grepl("^HBZ\\|",rownames(cleanDat))))
#regrProts.idx
#5569 6087 -- take out
#previously: 6156 4641
#integer(0)  # if already removed

#previously: removed 2 regr proteins in line:  varPart.reg <- fitExtractVarPartModel(impute::impute.knn(cleanDat[-c(4641,6156),])$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = 8, type = "SOCK"))
varPart.b4b <- fitExtractVarPartModel(impute::impute.knn(as.matrix(cleanDat.4p13b4b[,]))$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp.b4b <- sortCols(varPart.b4b,FUN=median,last= c("Residuals"))

pdf(file="4p13b4b.19sites1xPAVregr+SiteRegress_Protect_age+sex+LRRN1-VariancePartition-PLASMA-7335x22392.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp.b4b, main="HDS 1.3ms - 4p13b4b - 19 sites 2x Regr(2PAV) + Site Regr Prot. age+sex+(e4 binary)" )

	SexSortOrder<-order(vp.b4b[["Sex"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b4b)) { vp.b4b[[i]]<-vp.b4b[[i]][SexSortOrder]; }
	rownames(vp.b4b)<-rownames(vp.b4b)[SexSortOrder]

	plotPercentBars( vp.b4b[1:50,]) + ggtitle( "Top Sex-covariates" )


	AgeSortOrder<-order(vp.b4b[["Age"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b4b)) { vp.b4b[[i]]<-vp.b4b[[i]][AgeSortOrder]; }
	rownames(vp.b4b)<-rownames(vp.b4b)[AgeSortOrder]

	plotPercentBars( vp.b4b[1:50,]) + ggtitle( "Top Age-covariates" )


        BatchSortOrder<-order(vp.b4b[["RegrBloodPreanalyticFactor.HNRNPA2B1"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4b)) { vp.b4b[[i]]<-vp.b4b[[i]][BatchSortOrder]; }
        rownames(vp.b4b)<-rownames(vp.b4b)[BatchSortOrder]

        plotPercentBars( vp.b4b[1:50,]) + ggtitle( "Top Time To Spin (HNRNPA2B1)-covariates" )


        BatchSortOrder<-order(vp.b4b[["RegrBloodPreanalyticFactor.HBZ"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4b)) { vp.b4b[[i]]<-vp.b4b[[i]][BatchSortOrder]; }
        rownames(vp.b4b)<-rownames(vp.b4b)[BatchSortOrder]

        plotPercentBars( vp.b4b[1:50,]) + ggtitle( "Top Preanalytical Factor 2 (HBZ)-covariates" )


        BatchSortOrder<-order(vp.b4b[["contributor_Fsplit"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4b)) { vp.b4b[[i]]<-vp.b4b[[i]][BatchSortOrder]; }
        rownames(vp.b4b)<-rownames(vp.b4b)[BatchSortOrder]

        plotPercentBars( vp.b4b[1:50,]) + ggtitle( "Top Contributor (site)-covariates" )


        BatchSortOrder<-order(vp.b4b[["APOE.E4carrier.imputed"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4b)) { vp.b4b[[i]]<-vp.b4b[[i]][BatchSortOrder]; }
        rownames(vp.b4b)<-rownames(vp.b4b)[BatchSortOrder]

        plotPercentBars( vp.b4b[1:50,]) + ggtitle( "Top APOE e4 carrier (Binary, imputed NA)-covariates" )


#	BatchSortOrder<-order(vp.b4b[["MMSE"]],decreasing=TRUE)
#	#rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
#	for (i in ls(vp.b4b)) { vp.b4b[[i]]<-vp.b4b[[i]][BatchSortOrder]; }
#	rownames(vp.b4b)<-rownames(vp.b4b)[BatchSortOrder]
#
#	plotPercentBars( vp.b4b[1:50,]) + ggtitle( "MMSE Cog. Score-covariates" )

dev.off()

varPart.b4b<-vp.b4b
saveRDS(varPart.b4b,"4p13b4b.19sites1xPAVregr+SiteRegress_Protect_age+sex+e4_carrierBinaryImputedStatus.varPart.b4b.RDS")


## Rerun the above 5 plots with contributor_code without Fsplit


#############################################################################
## 4p13b2. Variance Partition regressed (2PAV regression intrasite)+Site regressed; no protection (QC)

regvars.vp<-data.frame(numericMeta.reg)
regvars.vp$Sex<-factor(abs(regvars.vp$sex -2))
regvars.vp$Age<-as.numeric(regvars.vp$age_at_visit)
regvars.vp$contributor_code<-factor(regvars.vp$contributor_code)
regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1)
regvars.vp$RegrBloodPreanalyticFactor.HBZ<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HBZ)
#regvars.vp$APOE.E4carrier.Proxy.LRRN1<-as.numeric(cleanDat.4p13b2["LRRN1|Q6UXK5^SL025922@seq.11293.14",])
regvars.vp$APOE.E4carrier.imputed<-factor(numericMeta.reg$APOE4.carrier.imputed)

# too many missing values:
#form <- ~ age_at_visit+(1|Sex)+(1|raceAA)+(1|recruited_control)+(1|ad)+(1|ftd)+(1|pd)+(1|als)+(1|mci_sci)+(1|mi)+(1|C9Orf72)+(1|GRN)+(1|MAPT)+APOE4.Dose+(1|contributor_code)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+MMSE
#form <- ~ (1|contributor_code) +RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ
form <- ~ Age+(1|Sex)+(1|contributor_code)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+(1|APOE.E4carrier.imputed)  #+APOE.E4carrier.Proxy.LRRN1

library(variancePartition)

# (If regressed) REMOVE regressed proteins to avoid "Response variable 4641 has a variance of 0":
#regrProts.idx<-c(which(grepl("^HNRNPA2B1\\|",rownames(cleanDat))),which(grepl("^HBZ\\|",rownames(cleanDat))))
#regrProts.idx
#5569 6087 -- take out
#previously: 6156 4641
#integer(0)  # if already removed

#previously: removed 2 regr proteins in line:  varPart.reg <- fitExtractVarPartModel(impute::impute.knn(cleanDat[-c(4641,6156),])$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = 8, type = "SOCK"))
varPart.b2 <- fitExtractVarPartModel(impute::impute.knn(as.matrix(cleanDat.4p13b2[,]))$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp.b2 <- sortCols(varPart.b2,FUN=median,last= c("Residuals"))

pdf(file="4p13b2.contributor_code_FnotSplit_19sites1xPAVregr+SiteRegress_noProtect-VariancePartition-PLASMA-7335x22547.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp.b2, main="HDS 1.3ms - 4p13b2 - 19 sites 2x Regr(2PAV) + Site Regr, No Protect" )

	SexSortOrder<-order(vp.b2[["Sex"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][SexSortOrder]; }
	rownames(vp.b2)<-rownames(vp.b2)[SexSortOrder]

	plotPercentBars( vp.b2[1:50,]) + ggtitle( "Top Sex-covariates" )


	AgeSortOrder<-order(vp.b2[["Age"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][AgeSortOrder]; }
	rownames(vp.b2)<-rownames(vp.b2)[AgeSortOrder]

	plotPercentBars( vp.b2[1:50,]) + ggtitle( "Top Age-covariates" )


        BatchSortOrder<-order(vp.b2[["RegrBloodPreanalyticFactor.HNRNPA2B1"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][BatchSortOrder]; }
        rownames(vp.b2)<-rownames(vp.b2)[BatchSortOrder]

        plotPercentBars( vp.b2[1:50,]) + ggtitle( "Top Time To Spin (HNRNPA2B1)-covariates" )


        BatchSortOrder<-order(vp.b2[["RegrBloodPreanalyticFactor.HBZ"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][BatchSortOrder]; }
        rownames(vp.b2)<-rownames(vp.b2)[BatchSortOrder]

        plotPercentBars( vp.b2[1:50,]) + ggtitle( "Top Preanalytical Factor 2 (HBZ)-covariates" )


        BatchSortOrder<-order(vp.b2[["contributor_code"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][BatchSortOrder]; }
        rownames(vp.b2)<-rownames(vp.b2)[BatchSortOrder]

        plotPercentBars( vp.b2[1:50,]) + ggtitle( "Top Contributor (site)-covariates" )


        BatchSortOrder<-order(vp.b2[["APOE.E4carrier.imputed"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][BatchSortOrder]; }
        rownames(vp.b2)<-rownames(vp.b2)[BatchSortOrder]

        plotPercentBars( vp.b2[1:50,]) + ggtitle( "Top APOE e4 carrier (Binary, imputed NA)-covariates" )


#	BatchSortOrder<-order(vp.b2[["MMSE"]],decreasing=TRUE)
#	#rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
#	for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][BatchSortOrder]; }
#	rownames(vp.b2)<-rownames(vp.b2)[BatchSortOrder]
#
#	plotPercentBars( vp.b2[1:50,]) + ggtitle( "MMSE Cog. Score-covariates" )

dev.off()

varPart.b2<-vp.b2
saveRDS(varPart.b2,"4p13b2.contributor_code_FnotSplit_19sites1xPAVregr+SiteRegress_noProtect.varPart.b2.RDS")

#<REMOVED CODE FOR b3, b4, b5, b4b using impute.knn on respective cleanDat and contributor_code (F not Fsplit)> - artificial site effect due to imputation(?)


## Rerun with na.omit instead of impute.knn - Fsplit -- "stats" not loading in worker warning -- suspect computation with contributor_Fsplit not seen out of form/model.
## (Rerun later for accurate output)

#############################################################################
## 4p13b2. Variance Partition regressed (2PAV regression intrasite)+Site regressed; no protection (QC)

regvars.vp<-data.frame(numericMeta.reg)
regvars.vp$Sex<-factor(abs(regvars.vp$sex -2))
regvars.vp$Age<-as.numeric(regvars.vp$age_at_visit)
regvars.vp$contributor_Fsplit<-factor(regvars.vp$contributor_Fsplit)
regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1)
regvars.vp$RegrBloodPreanalyticFactor.HBZ<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HBZ)
#regvars.vp$APOE.E4carrier.Proxy.LRRN1<-as.numeric(cleanDat.4p13b2["LRRN1|Q6UXK5^SL025922@seq.11293.14",])
regvars.vp$APOE.E4carrier.imputed<-factor(numericMeta.reg$APOE4.carrier.imputed)

# too many missing values:
#form <- ~ age_at_visit+(1|Sex)+(1|raceAA)+(1|recruited_control)+(1|ad)+(1|ftd)+(1|pd)+(1|als)+(1|mci_sci)+(1|mi)+(1|C9Orf72)+(1|GRN)+(1|MAPT)+APOE4.Dose+(1|contributor_code)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+MMSE
#form <- ~ (1|contributor_code) +RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ
form <- ~ Age+(1|Sex)+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+(1|APOE.E4carrier.imputed)  #+APOE.E4carrier.Proxy.LRRN1

library(variancePartition)

# (If regressed) REMOVE regressed proteins to avoid "Response variable 4641 has a variance of 0":
#regrProts.idx<-c(which(grepl("^HNRNPA2B1\\|",rownames(cleanDat))),which(grepl("^HBZ\\|",rownames(cleanDat))))
#regrProts.idx
#5569 6087 -- take out
#previously: 6156 4641
#integer(0)  # if already removed

#previously: removed 2 regr proteins in line:  varPart.reg <- fitExtractVarPartModel(impute::impute.knn(cleanDat[-c(4641,6156),])$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = 8, type = "SOCK"))
varPart.b2 <- fitExtractVarPartModel(na.omit(as.matrix(cleanDat.4p13b2[,])), form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp.b2 <- sortCols(varPart.b2,FUN=median,last= c("Residuals"))

pdf(file="4p13b2.19sites1xPAVregr+SiteRegress_noProtect-VariancePartition(na.omit)-PLASMA-7335x22547.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp.b2, main="HDS 1.3ms - 4p13b2 - 19 sites 2x Regr(2PAV) + Site Regr, No Protect" )

	SexSortOrder<-order(vp.b2[["Sex"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][SexSortOrder]; }
	rownames(vp.b2)<-rownames(vp.b2)[SexSortOrder]

	plotPercentBars( vp.b2[1:50,]) + ggtitle( "Top Sex-covariates" )


	AgeSortOrder<-order(vp.b2[["Age"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][AgeSortOrder]; }
	rownames(vp.b2)<-rownames(vp.b2)[AgeSortOrder]

	plotPercentBars( vp.b2[1:50,]) + ggtitle( "Top Age-covariates" )


        BatchSortOrder<-order(vp.b2[["RegrBloodPreanalyticFactor.HNRNPA2B1"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][BatchSortOrder]; }
        rownames(vp.b2)<-rownames(vp.b2)[BatchSortOrder]

        plotPercentBars( vp.b2[1:50,]) + ggtitle( "Top Time To Spin (HNRNPA2B1)-covariates" )


        BatchSortOrder<-order(vp.b2[["RegrBloodPreanalyticFactor.HBZ"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][BatchSortOrder]; }
        rownames(vp.b2)<-rownames(vp.b2)[BatchSortOrder]

        plotPercentBars( vp.b2[1:50,]) + ggtitle( "Top Preanalytical Factor 2 (HBZ)-covariates" )


        BatchSortOrder<-order(vp.b2[["contributor_Fsplit"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][BatchSortOrder]; }
        rownames(vp.b2)<-rownames(vp.b2)[BatchSortOrder]

        plotPercentBars( vp.b2[1:50,]) + ggtitle( "Top Contributor (site)-covariates" )


        BatchSortOrder<-order(vp.b2[["APOE.E4carrier.imputed"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][BatchSortOrder]; }
        rownames(vp.b2)<-rownames(vp.b2)[BatchSortOrder]

        plotPercentBars( vp.b2[1:50,]) + ggtitle( "Top APOE e4 carrier (Binary, imputed NA)-covariates" )


#	BatchSortOrder<-order(vp.b2[["MMSE"]],decreasing=TRUE)
#	#rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
#	for (i in ls(vp.b2)) { vp.b2[[i]]<-vp.b2[[i]][BatchSortOrder]; }
#	rownames(vp.b2)<-rownames(vp.b2)[BatchSortOrder]
#
#	plotPercentBars( vp.b2[1:50,]) + ggtitle( "MMSE Cog. Score-covariates" )

dev.off()

varPart.b2<-vp.b2
saveRDS(varPart.b2,"4p13b2.19sites1xPAVregr+SiteRegress_noProtect.na.omit-varPart.b2.RDS")


#############################################################################
## 4p13b3. Variance Partition regressed (2PAV regression intrasite)+Site regressed; protect Age+Sex (QC)

regvars.vp<-data.frame(numericMeta.reg.b345)
regvars.vp$Sex<-factor(abs(regvars.vp$sex -2))
regvars.vp$Age<-as.numeric(regvars.vp$age_at_visit)
regvars.vp$contributor_Fsplit<-factor(regvars.vp$contributor_Fsplit)
regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1)
regvars.vp$RegrBloodPreanalyticFactor.HBZ<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HBZ)
#regvars.vp$APOE.E4carrier.Proxy.LRRN1<-as.numeric(cleanDat.4p13b2["LRRN1|Q6UXK5^SL025922@seq.11293.14",])
regvars.vp$APOE.E4carrier.imputed<-factor(regvars.vp$APOE4.carrier.imputed)

# too many missing values:
#form <- ~ age_at_visit+(1|Sex)+(1|raceAA)+(1|recruited_control)+(1|ad)+(1|ftd)+(1|pd)+(1|als)+(1|mci_sci)+(1|mi)+(1|C9Orf72)+(1|GRN)+(1|MAPT)+APOE4.Dose+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+MMSE
#form <- ~ (1|contributor_Fsplit) +RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ
form <- ~ Age+(1|Sex)+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+(1|APOE.E4carrier.imputed)  #+APOE.E4carrier.Proxy.LRRN1

library(variancePartition)

# (If regressed) REMOVE regressed proteins to avoid "Response variable 4641 has a variance of 0":
#regrProts.idx<-c(which(grepl("^HNRNPA2B1\\|",rownames(cleanDat))),which(grepl("^HBZ\\|",rownames(cleanDat))))
#regrProts.idx
#5569 6087 -- take out
#previously: 6156 4641
#integer(0)  # if already removed

#previously: removed 2 regr proteins in line:  varPart.reg <- fitExtractVarPartModel(impute::impute.knn(cleanDat[-c(4641,6156),])$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = 8, type = "SOCK"))
varPart.b3 <- fitExtractVarPartModel(na.omit(as.matrix(cleanDat.4p13b3[,])), form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp.b3 <- sortCols(varPart.b3,FUN=median,last= c("Residuals"))

pdf(file="4p13b3.contributor_Fsplit_19sites1xPAVregr+SiteRegress_Protect_age+sex-VariancePartition-PLASMA-7335x22392.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp.b3, main="HDS 1.3ms - 4p13b3 - 19 sites 2x Regr(2PAV) + Site Regr Prot. age+sex" )

	SexSortOrder<-order(vp.b3[["Sex"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b3)) { vp.b3[[i]]<-vp.b3[[i]][SexSortOrder]; }
	rownames(vp.b3)<-rownames(vp.b3)[SexSortOrder]

	plotPercentBars( vp.b3[1:50,]) + ggtitle( "Top Sex-covariates" )


	AgeSortOrder<-order(vp.b3[["Age"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b3)) { vp.b3[[i]]<-vp.b3[[i]][AgeSortOrder]; }
	rownames(vp.b3)<-rownames(vp.b3)[AgeSortOrder]

	plotPercentBars( vp.b3[1:50,]) + ggtitle( "Top Age-covariates" )


        BatchSortOrder<-order(vp.b3[["RegrBloodPreanalyticFactor.HNRNPA2B1"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b3)) { vp.b3[[i]]<-vp.b3[[i]][BatchSortOrder]; }
        rownames(vp.b3)<-rownames(vp.b3)[BatchSortOrder]

        plotPercentBars( vp.b3[1:50,]) + ggtitle( "Top Time To Spin (HNRNPA2B1)-covariates" )


        BatchSortOrder<-order(vp.b3[["RegrBloodPreanalyticFactor.HBZ"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b3)) { vp.b3[[i]]<-vp.b3[[i]][BatchSortOrder]; }
        rownames(vp.b3)<-rownames(vp.b3)[BatchSortOrder]

        plotPercentBars( vp.b3[1:50,]) + ggtitle( "Top Preanalytical Factor 2 (HBZ)-covariates" )


        BatchSortOrder<-order(vp.b3[["contributor_Fsplit"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b3)) { vp.b3[[i]]<-vp.b3[[i]][BatchSortOrder]; }
        rownames(vp.b3)<-rownames(vp.b3)[BatchSortOrder]

        plotPercentBars( vp.b3[1:50,]) + ggtitle( "Top Contributor (site)-covariates" )


        BatchSortOrder<-order(vp.b3[["APOE.E4carrier.imputed"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b3)) { vp.b3[[i]]<-vp.b3[[i]][BatchSortOrder]; }
        rownames(vp.b3)<-rownames(vp.b3)[BatchSortOrder]

        plotPercentBars( vp.b3[1:50,]) + ggtitle( "Top APOE e4 carrier (Binary, imputed NA)-covariates" )


#	BatchSortOrder<-order(vp.b3[["MMSE"]],decreasing=TRUE)
#	#rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
#	for (i in ls(vp.b3)) { vp.b3[[i]]<-vp.b3[[i]][BatchSortOrder]; }
#	rownames(vp.b3)<-rownames(vp.b3)[BatchSortOrder]
#
#	plotPercentBars( vp.b3[1:50,]) + ggtitle( "MMSE Cog. Score-covariates" )

dev.off()

varPart.b3<-vp.b3
saveRDS(varPart.b3,"4p13b3.contributor_Fsplit_19sites1xPAVregr+SiteRegress_Protect_age+sex.na.omit-varPart.b3.RDS")


#############################################################################
## 4p13b4. Variance Partition regressed (2PAV regression intrasite)+Site regressed; protect Age+Sex+LRRN1 (e4 carrier proxy) (QC)

regvars.vp<-data.frame(numericMeta.reg.b345)
regvars.vp$Sex<-factor(abs(regvars.vp$sex -2))
regvars.vp$Age<-as.numeric(regvars.vp$age_at_visit)
regvars.vp$contributor_Fsplit<-factor(regvars.vp$contributor_Fsplit)
regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1)
regvars.vp$RegrBloodPreanalyticFactor.HBZ<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HBZ)
#regvars.vp$APOE.E4carrier.Proxy.LRRN1<-as.numeric(cleanDat.4p13b2["LRRN1|Q6UXK5^SL025922@seq.11293.14",])
regvars.vp$APOE.E4carrier.imputed<-factor(regvars.vp$APOE4.carrier.imputed)

# too many missing values:
#form <- ~ age_at_visit+(1|Sex)+(1|raceAA)+(1|recruited_control)+(1|ad)+(1|ftd)+(1|pd)+(1|als)+(1|mci_sci)+(1|mi)+(1|C9Orf72)+(1|GRN)+(1|MAPT)+APOE4.Dose+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+MMSE
#form <- ~ (1|contributor_Fsplit) +RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ
form <- ~ Age+(1|Sex)+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+(1|APOE.E4carrier.imputed)  #+APOE.E4carrier.Proxy.LRRN1

library(variancePartition)

# (If regressed) REMOVE regressed proteins to avoid "Response variable 4641 has a variance of 0":
#regrProts.idx<-c(which(grepl("^HNRNPA2B1\\|",rownames(cleanDat))),which(grepl("^HBZ\\|",rownames(cleanDat))))
#regrProts.idx
#5569 6087 -- take out
#previously: 6156 4641
#integer(0)  # if already removed

#previously: removed 2 regr proteins in line:  varPart.reg <- fitExtractVarPartModel(impute::impute.knn(cleanDat[-c(4641,6156),])$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = 8, type = "SOCK"))
varPart.b4 <- fitExtractVarPartModel(na.omit(as.matrix(cleanDat.4p13b4[,])), form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp.b4 <- sortCols(varPart.b4,FUN=median,last= c("Residuals"))

pdf(file="4p13b4.contributor_Fsplit_19sites1xPAVregr+SiteRegress_Protect_age+sex+LRRN1-VariancePartition-PLASMA-7335x22392.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp.b4, main="HDS 1.3ms - 4p13b4 - 19 sites 2x Regr(2PAV) + Site Regr Prot. age+sex+LRRN1" )

	SexSortOrder<-order(vp.b4[["Sex"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b4)) { vp.b4[[i]]<-vp.b4[[i]][SexSortOrder]; }
	rownames(vp.b4)<-rownames(vp.b4)[SexSortOrder]

	plotPercentBars( vp.b4[1:50,]) + ggtitle( "Top Sex-covariates" )


	AgeSortOrder<-order(vp.b4[["Age"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b4)) { vp.b4[[i]]<-vp.b4[[i]][AgeSortOrder]; }
	rownames(vp.b4)<-rownames(vp.b4)[AgeSortOrder]

	plotPercentBars( vp.b4[1:50,]) + ggtitle( "Top Age-covariates" )


        BatchSortOrder<-order(vp.b4[["RegrBloodPreanalyticFactor.HNRNPA2B1"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4)) { vp.b4[[i]]<-vp.b4[[i]][BatchSortOrder]; }
        rownames(vp.b4)<-rownames(vp.b4)[BatchSortOrder]

        plotPercentBars( vp.b4[1:50,]) + ggtitle( "Top Time To Spin (HNRNPA2B1)-covariates" )


        BatchSortOrder<-order(vp.b4[["RegrBloodPreanalyticFactor.HBZ"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4)) { vp.b4[[i]]<-vp.b4[[i]][BatchSortOrder]; }
        rownames(vp.b4)<-rownames(vp.b4)[BatchSortOrder]

        plotPercentBars( vp.b4[1:50,]) + ggtitle( "Top Preanalytical Factor 2 (HBZ)-covariates" )


        BatchSortOrder<-order(vp.b4[["contributor_Fsplit"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4)) { vp.b4[[i]]<-vp.b4[[i]][BatchSortOrder]; }
        rownames(vp.b4)<-rownames(vp.b4)[BatchSortOrder]

        plotPercentBars( vp.b4[1:50,]) + ggtitle( "Top Contributor (site)-covariates" )


        BatchSortOrder<-order(vp.b4[["APOE.E4carrier.imputed"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4)) { vp.b4[[i]]<-vp.b4[[i]][BatchSortOrder]; }
        rownames(vp.b4)<-rownames(vp.b4)[BatchSortOrder]

        plotPercentBars( vp.b4[1:50,]) + ggtitle( "Top APOE e4 carrier (Binary, imputed NA)-covariates" )


#	BatchSortOrder<-order(vp.b4[["MMSE"]],decreasing=TRUE)
#	#rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
#	for (i in ls(vp.b4)) { vp.b4[[i]]<-vp.b4[[i]][BatchSortOrder]; }
#	rownames(vp.b4)<-rownames(vp.b4)[BatchSortOrder]
#
#	plotPercentBars( vp.b4[1:50,]) + ggtitle( "MMSE Cog. Score-covariates" )

dev.off()

varPart.b4<-vp.b4
saveRDS(varPart.b4,"4p13b4.contributor_Fsplit_19sites1xPAVregr+SiteRegress_Protect_age+sex+LRRN1.na.omit-varPart.b4.RDS")


#############################################################################
## 4p13b5. Variance Partition regressed (2PAV regression intrasite)+Site regressed; protect Age+Sex+log2(LRRN1/TBCA) (e4 carrier proxy) (QC)

regvars.vp<-data.frame(numericMeta.reg.b345)
regvars.vp$Sex<-factor(abs(regvars.vp$sex -2))
regvars.vp$Age<-as.numeric(regvars.vp$age_at_visit)
regvars.vp$contributor_Fsplit<-factor(regvars.vp$contributor_Fsplit)
regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1)
regvars.vp$RegrBloodPreanalyticFactor.HBZ<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HBZ)
#regvars.vp$APOE.E4carrier.Proxy.LRRN1<-as.numeric(cleanDat.4p13b2["LRRN1|Q6UXK5^SL025922@seq.11293.14",])
regvars.vp$APOE.E4carrier.imputed<-factor(regvars.vp$APOE4.carrier.imputed)

# too many missing values:
#form <- ~ age_at_visit+(1|Sex)+(1|raceAA)+(1|recruited_control)+(1|ad)+(1|ftd)+(1|pd)+(1|als)+(1|mci_sci)+(1|mi)+(1|C9Orf72)+(1|GRN)+(1|MAPT)+APOE4.Dose+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+MMSE
#form <- ~ (1|contributor_Fsplit) +RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ
form <- ~ Age+(1|Sex)+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+(1|APOE.E4carrier.imputed)  #+APOE.E4carrier.Proxy.LRRN1

library(variancePartition)

# (If regressed) REMOVE regressed proteins to avoid "Response variable 4641 has a variance of 0":
#regrProts.idx<-c(which(grepl("^HNRNPA2B1\\|",rownames(cleanDat))),which(grepl("^HBZ\\|",rownames(cleanDat))))
#regrProts.idx
#5569 6087 -- take out
#previously: 6156 4641
#integer(0)  # if already removed

#previously: removed 2 regr proteins in line:  varPart.reg <- fitExtractVarPartModel(impute::impute.knn(cleanDat[-c(4641,6156),])$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = 8, type = "SOCK"))
varPart.b5 <- fitExtractVarPartModel(na.omit(as.matrix(cleanDat.4p13b5[,])), form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp.b5 <- sortCols(varPart.b5,FUN=median,last= c("Residuals"))

pdf(file="4p13b5.contributor_Fsplit_19sites1xPAVregr+SiteRegress_Protect_age+sex+LRRN1minusTBCA-VariancePartition-PLASMA-7335x22392.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp.b5, main="HDS 1.3ms - 4p13b5 - 19 sites 2x Regr(2PAV) + Site Regr Prot. age+sex+(LRRN1-TBCA)" )

	SexSortOrder<-order(vp.b5[["Sex"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b5)) { vp.b5[[i]]<-vp.b5[[i]][SexSortOrder]; }
	rownames(vp.b5)<-rownames(vp.b5)[SexSortOrder]

	plotPercentBars( vp.b5[1:50,]) + ggtitle( "Top Sex-covariates" )


	AgeSortOrder<-order(vp.b5[["Age"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b5)) { vp.b5[[i]]<-vp.b5[[i]][AgeSortOrder]; }
	rownames(vp.b5)<-rownames(vp.b5)[AgeSortOrder]

	plotPercentBars( vp.b5[1:50,]) + ggtitle( "Top Age-covariates" )


        BatchSortOrder<-order(vp.b5[["RegrBloodPreanalyticFactor.HNRNPA2B1"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b5)) { vp.b5[[i]]<-vp.b5[[i]][BatchSortOrder]; }
        rownames(vp.b5)<-rownames(vp.b5)[BatchSortOrder]

        plotPercentBars( vp.b5[1:50,]) + ggtitle( "Top Time To Spin (HNRNPA2B1)-covariates" )


        BatchSortOrder<-order(vp.b5[["RegrBloodPreanalyticFactor.HBZ"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b5)) { vp.b5[[i]]<-vp.b5[[i]][BatchSortOrder]; }
        rownames(vp.b5)<-rownames(vp.b5)[BatchSortOrder]

        plotPercentBars( vp.b5[1:50,]) + ggtitle( "Top Preanalytical Factor 2 (HBZ)-covariates" )


        BatchSortOrder<-order(vp.b5[["contributor_Fsplit"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b5)) { vp.b5[[i]]<-vp.b5[[i]][BatchSortOrder]; }
        rownames(vp.b5)<-rownames(vp.b5)[BatchSortOrder]

        plotPercentBars( vp.b5[1:50,]) + ggtitle( "Top Contributor (site)-covariates" )


        BatchSortOrder<-order(vp.b5[["APOE.E4carrier.imputed"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b5)) { vp.b5[[i]]<-vp.b5[[i]][BatchSortOrder]; }
        rownames(vp.b5)<-rownames(vp.b5)[BatchSortOrder]

        plotPercentBars( vp.b5[1:50,]) + ggtitle( "Top APOE e4 carrier (Binary, imputed NA)-covariates" )


#	BatchSortOrder<-order(vp.b5[["MMSE"]],decreasing=TRUE)
#	#rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
#	for (i in ls(vp.b5)) { vp.b5[[i]]<-vp.b5[[i]][BatchSortOrder]; }
#	rownames(vp.b5)<-rownames(vp.b5)[BatchSortOrder]
#
#	plotPercentBars( vp.b5[1:50,]) + ggtitle( "MMSE Cog. Score-covariates" )

dev.off()

varPart.b5<-vp.b5
saveRDS(varPart.b5,"4p13b5.contributor_Fsplit_19sites1xPAVregr+SiteRegress_Protect_age+sex+LRRN1minusTBCA.na.omit-varPart.b5.RDS")


#############################################################################
## 4p13b4b. Variance Partition regressed (2PAV regression intrasite)+Site regressed; protect Age+Sex+E4 carrier binary status (NA imputed) (QC)

regvars.vp<-data.frame(numericMeta.reg.b345)
regvars.vp$Sex<-factor(abs(regvars.vp$sex -2))
regvars.vp$Age<-as.numeric(regvars.vp$age_at_visit)
regvars.vp$contributor_Fsplit<-factor(regvars.vp$contributor_Fsplit)
regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1)
regvars.vp$RegrBloodPreanalyticFactor.HBZ<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HBZ)
#regvars.vp$APOE.E4carrier.Proxy.LRRN1<-as.numeric(cleanDat.4p13b2["LRRN1|Q6UXK5^SL025922@seq.11293.14",])
regvars.vp$APOE.E4carrier.imputed<-factor(regvars.vp$APOE4.carrier.imputed)

# too many missing values:
#form <- ~ age_at_visit+(1|Sex)+(1|raceAA)+(1|recruited_control)+(1|ad)+(1|ftd)+(1|pd)+(1|als)+(1|mci_sci)+(1|mi)+(1|C9Orf72)+(1|GRN)+(1|MAPT)+APOE4.Dose+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+MMSE
#form <- ~ (1|contributor_Fsplit) +RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ
form <- ~ Age+(1|Sex)+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+(1|APOE.E4carrier.imputed)  #+APOE.E4carrier.Proxy.LRRN1

library(variancePartition)

# (If regressed) REMOVE regressed proteins to avoid "Response variable 4641 has a variance of 0":
#regrProts.idx<-c(which(grepl("^HNRNPA2B1\\|",rownames(cleanDat))),which(grepl("^HBZ\\|",rownames(cleanDat))))
#regrProts.idx
#5569 6087 -- take out
#previously: 6156 4641
#integer(0)  # if already removed

#previously: removed 2 regr proteins in line:  varPart.reg <- fitExtractVarPartModel(impute::impute.knn(cleanDat[-c(4641,6156),])$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = 8, type = "SOCK"))
varPart.b4b <- fitExtractVarPartModel(na.omit(as.matrix(cleanDat.4p13b4b[,])), form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp.b4b <- sortCols(varPart.b4b,FUN=median,last= c("Residuals"))

pdf(file="4p13b4b.contributor_Fsplit_19sites1xPAVregr+SiteRegress_Protect_age+sex+carrierBinaryImputedStatus-VariancePartition-PLASMA-7335x22392.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp.b4b, main="HDS 1.3ms - 4p13b4b - 19 sites 2x Regr(2PAV) + Site Regr Prot. age+sex+(e4 binary)" )

	SexSortOrder<-order(vp.b4b[["Sex"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b4b)) { vp.b4b[[i]]<-vp.b4b[[i]][SexSortOrder]; }
	rownames(vp.b4b)<-rownames(vp.b4b)[SexSortOrder]

	plotPercentBars( vp.b4b[1:50,]) + ggtitle( "Top Sex-covariates" )


	AgeSortOrder<-order(vp.b4b[["Age"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b4b)) { vp.b4b[[i]]<-vp.b4b[[i]][AgeSortOrder]; }
	rownames(vp.b4b)<-rownames(vp.b4b)[AgeSortOrder]

	plotPercentBars( vp.b4b[1:50,]) + ggtitle( "Top Age-covariates" )


        BatchSortOrder<-order(vp.b4b[["RegrBloodPreanalyticFactor.HNRNPA2B1"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4b)) { vp.b4b[[i]]<-vp.b4b[[i]][BatchSortOrder]; }
        rownames(vp.b4b)<-rownames(vp.b4b)[BatchSortOrder]

        plotPercentBars( vp.b4b[1:50,]) + ggtitle( "Top Time To Spin (HNRNPA2B1)-covariates" )


        BatchSortOrder<-order(vp.b4b[["RegrBloodPreanalyticFactor.HBZ"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4b)) { vp.b4b[[i]]<-vp.b4b[[i]][BatchSortOrder]; }
        rownames(vp.b4b)<-rownames(vp.b4b)[BatchSortOrder]

        plotPercentBars( vp.b4b[1:50,]) + ggtitle( "Top Preanalytical Factor 2 (HBZ)-covariates" )


        BatchSortOrder<-order(vp.b4b[["contributor_Fsplit"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4b)) { vp.b4b[[i]]<-vp.b4b[[i]][BatchSortOrder]; }
        rownames(vp.b4b)<-rownames(vp.b4b)[BatchSortOrder]

        plotPercentBars( vp.b4b[1:50,]) + ggtitle( "Top Contributor (site)-covariates" )


        BatchSortOrder<-order(vp.b4b[["APOE.E4carrier.imputed"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4b)) { vp.b4b[[i]]<-vp.b4b[[i]][BatchSortOrder]; }
        rownames(vp.b4b)<-rownames(vp.b4b)[BatchSortOrder]

        plotPercentBars( vp.b4b[1:50,]) + ggtitle( "Top APOE e4 carrier (Binary, imputed NA)-covariates" )


#	BatchSortOrder<-order(vp.b4b[["MMSE"]],decreasing=TRUE)
#	#rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
#	for (i in ls(vp.b4b)) { vp.b4b[[i]]<-vp.b4b[[i]][BatchSortOrder]; }
#	rownames(vp.b4b)<-rownames(vp.b4b)[BatchSortOrder]
#
#	plotPercentBars( vp.b4b[1:50,]) + ggtitle( "MMSE Cog. Score-covariates" )

dev.off()

varPart.b4b<-vp.b4b
saveRDS(varPart.b4b,"4p13b4b.contributor_Fsplit_19sites1xPAVregr+SiteRegress_Protect_age+sex+e4_carrierBinaryImputedStatus.na.omit-varPart.b4b.RDS")


## Better histogram using top 6 correlates (3 + and 3 -)
SPC.med=median(cleanDat.4p13b3["SPC25|Q9HBM1",],na.rm=T)
CTF1.med=median(cleanDat.4p13b3["CTF1|Q16619^SL002783@seq.13732.79",],na.rm=T)
NEFL.med=median(cleanDat.4p13b3["NEFL|P07196",],na.rm=T)
S100A13.med=median(cleanDat.4p13b3["S100A13|Q99584",],na.rm=T)
hist.data.SPC=((cleanDat.4p13b3["SPC25|Q9HBM1",] - SPC.med) + (cleanDat.4p13b3["LRRN1|Q6UXK5^SL025922@seq.11293.14",] - LRRN1.med) + (cleanDat.4p13b3["CTF1|Q16619^SL002783@seq.13732.79",] - CTF1.med) - (cleanDat.4p13b3["NEFL|P07196",] - NEFL.med) - (cleanDat.4p13b3["S100A13|Q99584",] - S100A13.med) - (cleanDat.4p13b3["TBCA|O75347",] - TBCA.med))
hist(hist.data.SPC, breaks=100, xlab="log2(abundance ratio):  SPC+LRRN1+CTF1(median-centered) / TBCA+S100A13+NEFL(median-centered)", main="APOE e4 Carrier Status Best Nonmissing Correlate")
hist(hist.data.SPC[which(numericMeta.reg.b345$APOE==24 | numericMeta.reg.b345$APOE==34)],breaks=60,col="#BBFFBB40",add=T)  #green overlay
hist(hist.data.SPC[which(numericMeta.reg.b345$APOE==44)],breaks=70,col="#FFBBBB99",add=T)  #red overlay
hist(hist.data.SPC[which(numericMeta.reg.b345$APOE==33)],breaks=100,col="#BBBBFF99",add=T)  #blue
legend("topright",c("e4/e4","e3/e3", "e4 het"),fill=c("#FFBBBB99","#BBBBFF99","#BBFFBB40"))
abline(v=0.625,col="darkgreen",lty=2, lwd=2.2)
abline(v=4.25,col="maroon",lty=2, lwd=2.2)
# saved plot to powerpoint (capture)

## APOE genotype epsilon allele Imputation via ML training 3 multi-class learners: 1) elastic-net multinomial GLM; 2) extreme-gradient boosting; 3) random forest within 5-fold CV framework
## Ensemble the three by soft voting, average class probability; collect accuracy, macro-F1, and confusion matrices; refit the ensemble on all samples return a predict()
## wrapper for genotype imputation on new expression log2(abundance) matrix.

#https://chatgpt.com/share/681a76bb-5604-8007-a568-60ee16f5639d

## Ground Truth APOE genotypes
gt.APOE<-rep(NA,length(numericMeta.reg.b345$APOE))
gt.APOE[which(numericMeta.reg.b345$APOE==22)]<-"e2/e2"
gt.APOE[which(numericMeta.reg.b345$APOE==23)]<-"e2/e3"
gt.APOE[which(numericMeta.reg.b345$APOE==33)]<-"e3/e3"
gt.APOE[which(numericMeta.reg.b345$APOE==24)]<-"e2/e4"
gt.APOE[which(numericMeta.reg.b345$APOE==34)]<-"e3/e4"
gt.APOE[which(numericMeta.reg.b345$APOE==44)]<-"e4/e4"
names(gt.APOE)<-rownames(numericMeta.reg.b345)


#save.image("c:/Users/workspace/Downloads/4p13b4b.SecondRegressionsComplete.19sites_Fsplit+APOEgenoIMPUTATION_Windows02.RData")  ## Full run (16 GB)
#save.image("c:/Users/workspace/Downloads/4p13b4b.SecondRegressionsComplete.19sites_Fsplit+APOEgenoIMPUTATION.RData")  ## Saved Rerun above, without tSNE, VP QC steps (14 GB)


## VIM kNN disabled (impute input to the function)
##############################################################################
##  Parallel fit_APOE_ensemble()  -------------------------------------------
##############################################################################
library(caret)
library(glmnet)
library(xgboost)
library(ranger)
library(progressr)
#library(VIM)
library(dplyr)
#library(doParallel)
library(future)
library(doFuture)
library(doRNG)

fit_APOE_ensemble_par <- function(expr, APOE_gt,
                                  nfold = 5, nrep = 5,
                                  ncores = parallel::detectCores() - 1,
                                  seed   = 1) {

  set.seed(seed)

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

  plan(multisession, workers = ncores)  # or multicore on linux/macOS
  registerDoFuture()
#  handlers("txtprogressbar")
#  handlers("cli")
#  handlers(list(
#    handler_txtprogressbar(),   # classic bar
#    handler_newline()           # prints the p() message lines
#  ))
#  cl <- makeCluster(ncores, outfile="")  # outfile="" prints worker msg
##  cl <- makeCluster(c(rep("localhost",ncores)),type="PSOCK",)
#  registerDoParallel(cl)

  ## ------------------------------------------------------------------------
  ## 2.  run each fold in a worker  -----------------------------------------
  ## ------------------------------------------------------------------------

  with_progress({                              # << all progress lives here
    n_sub <- 3   # Progress bar increments 3x per fold
    p <- progressor(steps=length(cvIndex) * n_sub )     #along = cvIndex)          # one step per fold

    p(sprintf("Initializing %d workers...", ncores), amount=0)


  metrics <- foreach(fold = seq_along(cvIndex),
                     .combine   = rbind,
                     .export = "p",            # let workers see 'p'
                     .packages  = c("progressr","glmnet","xgboost","ranger","dplyr")) %dorng% {

    ## ----------------- announce fold start -----------------------------
    p(sprintf("fold %d/%d  -  started", fold, length(cvIndex)), amount=0)

    set.seed(seed + fold)                        # reproducible inside worker
    tr_idx <- cvIndex[[fold]]
    te_idx <- setdiff(seq_len(nrow(expr)), tr_idx)

    ## ---------------- preprocessing ---------------------------------------
#    prep_expr <- function(mat, thr_miss = 0.20) {
#      keep <- which(colMeans(is.na(mat)) <= thr_miss)
#      mat  <- mat[, keep, drop = FALSE]
##      mat  <- kNN(mat, k = 5, imp_var = FALSE)   # saved time - preimputed
#      scale(mat)
#    }

#    X_tr <- prep_expr(expr[tr_idx, ])
    ## ---------- preprocessing  (no kNN, already imputed) --------------
    prep <- function(m) scale(m[, colMeans(is.na(m)) <= .20, drop = FALSE])
    X_tr <- prep(expr[tr_idx, ])

# ------------------------------------------------------------------------
# ANNOTATION: Develop APOE genotype/proxy prediction functions and
# imputation candidate calls from the first-pass regressed data.
# ------------------------------------------------------------------------
    X_te <- scale(expr[te_idx, colnames(X_tr)],
                  center = attr(X_tr, "scaled:center"),
                  scale  = attr(X_tr, "scaled:scale"))
    X_te[is.na(X_te)] <- 0

    levels_all <- sort(unique(APOE_gt))            # six possible APOE genotypes
    num_class   <- length(levels_all)              # == 6

    y_tr <- factor(APOE_gt[tr_idx], levels = levels_all)
    y_te <- factor(APOE_gt[te_idx], levels = levels_all)

    #For robustness, drop samples and labels whose labels are NA (should not have been input)
    keep_tr <- !is.na(y_tr)
    keep_te <- !is.na(y_te)

    X_tr <- X_tr[keep_tr, , drop = FALSE]
    y_tr <- y_tr[keep_tr]

    X_te <- X_te[keep_te, , drop = FALSE]
    y_te <- y_te[keep_te]

    ## ---------------- glmnet ---------------------------------------------
    cv_glm <- cv.glmnet(X_tr, y_tr, family = "multinomial",
                        type.measure = "class", parallel = FALSE)
    p(sprintf("fold %d/%d  •  glmnet done", fold, length(cvIndex)), amount=1)

    p_glm  <- predict(cv_glm, X_te, s = "lambda.min",
                      type = "response")[,,1]

    ## ---------------- xgboost --------------------------------------------
    dtr <- xgb.DMatrix(X_tr, label = as.numeric(y_tr) - 1)
    dval<- xgb.DMatrix(X_te, label = as.numeric(y_te) - 1)

    xpar <- list(eta = 0.1, max_depth = 6, subsample = 0.8,
                 colsample_bytree = 0.8,
                 objective = "multi:softprob",
                 num_class = num_class,     # length(levels(y_tr)),
                 nthread = 1)               # <- single thread per worker

    bst  <- xgb.train(params=xpar, data=dtr, watchlist=list(train=dtr, eval=dval), nrounds = 200,
                      verbose = 0, early_stopping_rounds = 20)
    p_xgb <- matrix(predict(bst, X_te),
                    ncol = length(levels(y_tr)), byrow = TRUE)

    p(sprintf("fold %d/%d  •  xgboost done", fold, length(cvIndex)), amount=1)

    colnames(p_xgb) <- levels(y_tr)

    ## ---------------- random-forest --------------------------------------
    rf  <- ranger(y_tr ~ ., data = data.frame(y_tr, X_tr),
                  probability = TRUE, num.trees = 500,
                  num.threads = 1)          # <- single thread
    p(sprintf("fold %d/%d  •  ranger done", fold, length(cvIndex)), amount=1)

    p_rf <- predict(rf, data.frame(X_te))$predictions
    colnames(p_rf) <- rf$forest$levels

    ## ---------------- ensemble vote --------------------------------------
    p_avg <- (p_glm + p_xgb + p_rf) / 3
    y_hat <- factor(colnames(p_avg)[max.col(p_avg)],
                    levels = levels(y_tr))

    cm <- caret::confusionMatrix(y_hat, y_te)

#    p(sprintf("fold %d/%d  •  finished", fold, length(cvIndex)), amount=1)

    data.frame(Accuracy = cm$overall["Accuracy"],
               MacroF1  = mean(cm$byClass[,"F1"]),
               Fold     = fold)
  } # foreach
  }) # with_progress

#  stopCluster(cl)                               # tidy up
#  registerDoSEQ()                               # back to sequential

  cat(sprintf("CV Accuracy  %.3f ± %.3f\n",
              mean(metrics$Accuracy), sd(metrics$Accuracy)))
  cat(sprintf("CV Macro-F1  %.3f ± %.3f\n\n",
              mean(metrics$MacroF1),  sd(metrics$MacroF1)))

  ## ------------------------------------------------------------------------
  ## 3.  fit final ensemble on all data  (sequential) -----------------------
  ## ------------------------------------------------------------------------
  prep_expr <- function(mat, thr_miss = 0.20) {
    keep <- which(colMeans(is.na(mat)) <= thr_miss)
    mat  <- mat[, keep, drop = FALSE]
#    mat  <- VIM::kNN(mat, k = 5, imp_var = FALSE)
    scale(mat)
  }

  with_progress({
    p2 <- progressor(steps=3)
    p2(sprintf("Starting fit of final ensemble on all data (serial/non-parallel)..."), amount=0)

  levels_all <- sort(unique(APOE_gt))

  X_all <- prep_expr(expr)
  y_all <- factor(APOE_gt, levels=levels_all)

  final_glm <- cv.glmnet(X_all, y_all, family = "multinomial",
                         type.measure = "class", parallel = FALSE)
    p2(sprintf("GLM fit on full data  •  finished"), amount=1)

  xpar_all <- list(eta = 0.1, max_depth = 6, subsample = 0.8,
                   colsample_bytree = 0.8,
                   objective = "multi:softprob",
                   num_class = length(levels(y_all)),
                   nthread = ncores)          # can use all cores now

  d.all    <- xgb.DMatrix(X_all, label = as.numeric(y_all)-1)
  bst_all  <- xgb.train(params=xpar_all, data=d.all,
                        watchlist=list(train=d.all),
                        nrounds = 200, verbose = 0,   #instead of 200, could use best #rounds: bst$best_iteration
                        early_stopping_rounds = 20)
    p2(sprintf("XGboost fit on full data  •  finished"), amount=1)

  rf_all <- ranger(y_all ~ ., data = data.frame(y_all, X_all),
                   probability = TRUE, num.trees = 500,
                   num.threads = ncores)
    p2(sprintf("Random Forest fit on full data  •  finished"), amount=1)

    p2(sprintf("Fit on full data  •  finished"), amount=0)

  ## ------------------------------------------------------------------------
  ## 4.  prediction wrapper -------------------------------------------------
  ## ------------------------------------------------------------------------
  function(new_expr) {
    new_expr <- prep_expr(rbind(expr[1,,drop=FALSE], new_expr))[-1,,drop=FALSE]
    new_expr <- new_expr[, colnames(X_all), drop = FALSE]

    p1 <- predict(final_glm, new_expr, s = "lambda.min",
                  type = "response")[,,1]

    p2 <- matrix(predict(bst_all, new_expr),
                 ncol = length(levels(y_all)), byrow = TRUE)
    colnames(p2) <- levels(y_all)

    p3 <- predict(rf_all, data.frame(new_expr))$predictions
    colnames(p3) <- rf_all$forest$levels

    p <- (p1 + p2 + p3) / 3
    factor(colnames(p)[max.col(p)], levels = levels(y_all))
  }
  }) # with_progress
}
options(future.globals.maxSize= 1024^3)  #1Gb Total size of all global objects that need to be exported - up from 500MB

## Create predict() wrapper  -- basic run (step 1 of 2)
#predict_APOE.b3<-fit_APOE_ensemble(t(cleanDat.4p13b3[,which(!is.na(gt.APOE))]), na.omit(gt.APOE))  # not parallel
#predict_APOE.b3<-fit_APOE_ensemble_par(t(cleanDat.4p13b3[,which(!is.na(gt.APOE))]), na.omit(gt.APOE), ncores=8)

## Impute APOE epsilon 6 genotypes in the 7k samples missing it (step 2 of 2)
#imputed.APOE.b3<-predict_APOE.b3(t(cleanDat.4p13b3[,which(is.na(gt.APOE))]))

## Overcoming technical hurdles - main hurdle is NA values in training data for step 1. Options to overcome:
#option 1 - VIM::kNN on a subset of higher correlation asssays
library(WGCNA)
e4.bicor.to.siteCorr.b3.assays<-bicor(t(cleanDat.4p13b3),numericMeta.reg.b345$APOE4.carrier.imputed, use='p')
# sort by bicor - select top 1250+, 1250- correlate assay names to keep
APOE4.assays.keep<-names(e4.bicor.to.siteCorr.b3.assays[order(unlist(t(e4.bicor.to.siteCorr.b3.assays)),decreasing=TRUE),][c(1:1250,(7335-1250):7335)])
#! training.cleanDat<-VIM::kNN(cleanDat.4p13b3[APOE4.assays.keep, which(!is.na(gt.APOE))], k = 5, imp_var = FALSE)    # still slower than 30 min  - cannot be gracefully killed

##option 2 - predict for imputation
#library(caret)
#pre <- preProcess(t(cleanDat.4p13b3[APOE4.assays.keep, which(!is.na(gt.APOE))]), method = c("medianImpute", "center", "scale"))
#training.cleanDat <- predict(pre, t(cleanDat.4p13b3[APOE4.assays.keep, which(!is.na(gt.APOE))]))          #  no NA, already scaled
#predict_APOE.b3<-fit_APOE_ensemble_par(training.cleanDat, na.omit(gt.APOE), ncores=8)

#option 3 - impute::impute.knn and reduce input expr mat size
training.cleanDat <- impute::impute.knn(cleanDat.4p13b3[, which(!is.na(gt.APOE))])$data   #  no NA
#sample.subset<-sort(sample(length(na.omit(gt.APOE)),750))
#table(na.omit(gt.APOE)[sample.subset])
set.seed(1)
# sample 150 of each genotype for training (n=810 total)
keep_idx<-unlist(lapply(split(seq_along( na.omit(gt.APOE)), na.omit(gt.APOE)), function(i) if (length(i) > 150) sample(i,150) else i), use.names=FALSE)
table(na.omit(gt.APOE)[keep_idx])  # 60 2/2, 150 ea of the other 5 genotypes; took <40 min to complete training on 28/32 cores at 3.2GHz (128 GB RAM rec.)
predict_APOE.b3<-fit_APOE_ensemble_par(t(training.cleanDat[APOE4.assays.keep[c(1:250,2251:length(APOE4.assays.keep))], keep_idx]), na.omit(gt.APOE)[keep_idx], ncores=28)
#CV Accuracy  0.909 ± 0.021
#CV Macro-F1  0.914 ± 0.019

# Impute APOE epsilon 6 genotypes in the 7k samples missing it (step 2 of 2)
imputed.APOE.b3<-predict_APOE.b3(t(impute::impute.knn(cleanDat.4p13b3)$data[APOE4.assays.keep[c(1:250,2251:length(APOE4.assays.keep))], which(is.na(gt.APOE))]))


#option 4 - full data of 13k samples with known APOE genotype having no missing for training (406, site K? are only missing 1 and thrown out, too) - est 10h
rownames(cleanDat.4p13b3)[which(apply(cleanDat.4p13b3[,which(numericMeta.reg.b345$contributor_code=="K")],1,function(x) length(which(is.na(x)))==406))]
#"IRF6|O14896"
training.cleanDat.noNA<-t(na.omit(t(cleanDat.4p13b3[which(!rownames(cleanDat.4p13b3) %in% c("IRF6|O14896")), which(!is.na(gt.APOE))])))
dim(training.cleanDat.noNA)
#  7334 13004   # 7335 13004 with IRF6 in (site K does not have APOE genotypes) -- otherwise it would be 13410 samples for training
training.gt.APOE<-gt.APOE[colnames(training.cleanDat.noNA)]
predict_APOE.b3.all.noNA<-fit_APOE_ensemble_par(t(training.cleanDat.noNA), training.gt.APOE, ncores=14)
# Impute APOE epsilon 6 genotypes in the 7k samples missing it (step 2 of 2)
imputed.APOE.b3.all.noNA<-predict_APOE.b3.all.noNA(t(impute::impute.knn(cleanDat.4p13b3)$data[rownames(training.cleanDat.noNA),which(is.na(gt.APOE))]))

names(imputed.APOE.b3.all.noNA)<-rownames(t(cleanDat.4p13b3[rownames(training.cleanDat.noNA),which(is.na(gt.APOE))]))
saveRDS(predict_APOE.b3.all.noNA,"predict_APOE.b3.all.noNA.RDS")

saveRDS(imputed.APOE.b3.all.noNA,"imputed.APOE.b3.all.noNA.RDS")

predict_APOE.b3.all.noNA<-readRDS("predict_APOE.b3.all.noNA.RDS")
imputed.APOE.b3.all.noNA<-predict_APOE.b3.all.noNA(t(impute::impute.knn(cleanDat.4p13b3)$data[rownames(training.cleanDat.noNA),which(is.na(gt.APOE))]))
table(imputed.APOE.b3.all.noNA)

#UDS.map
UDS.map.APOE<-UDS.map
UDS.map.APOE$APOE.4cohort<-as.numeric(gsub("e","",traits.4cohort$APOE[match(UDS.map[,2],traits.4cohort$LoadedSampleName)]))
UDS.map.APOE$APOE.predict<-imputed.APOE.b3.all.noNA[match(UDS.map.APOE[,1],names(imputed.APOE.b3.all.noNA))]

UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE$APOE.4cohort==22)]<-"e2/e2"
UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE$APOE.4cohort==23)]<-"e2/e3"
UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE$APOE.4cohort==33)]<-"e3/e3"
UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE$APOE.4cohort==24)]<-"e2/e4"
UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE$APOE.4cohort==34)]<-"e3/e4"
UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE$APOE.4cohort==44)]<-"e4/e4"

table(paste0(UDS.map.APOE$APOE.4cohort,".",UDS.map.APOE$APOE.predict))


BH.map.APOE<-BH.map
BH.map.APOE$APOE.4cohort<-as.numeric(gsub("e","",traits.4cohort$APOE[match(BH.map[,2],traits.4cohort$LoadedSampleName)]))
BH.map.APOE$APOE.predict<-imputed.APOE.b3.all.noNA[match(BH.map.APOE[,1],names(imputed.APOE.b3.all.noNA))]

BH.map.APOE$APOE.4cohort[which(BH.map.APOE$APOE.4cohort==22)]<-"e2/e2"
BH.map.APOE$APOE.4cohort[which(BH.map.APOE$APOE.4cohort==23)]<-"e2/e3"
BH.map.APOE$APOE.4cohort[which(BH.map.APOE$APOE.4cohort==33)]<-"e3/e3"
BH.map.APOE$APOE.4cohort[which(BH.map.APOE$APOE.4cohort==24)]<-"e2/e4"
BH.map.APOE$APOE.4cohort[which(BH.map.APOE$APOE.4cohort==34)]<-"e3/e4"
BH.map.APOE$APOE.4cohort[which(BH.map.APOE$APOE.4cohort==44)]<-"e4/e4"

table(paste0(BH.map.APOE$APOE.4cohort,".",BH.map.APOE$APOE.predict))


RM.map.APOE<-RM.map
RM.map.APOE$APOE.4cohort<-as.numeric(gsub("e","",traits.4cohort$APOE[match(RM.map[,2],traits.4cohort$LoadedSampleName)]))
RM.map.APOE$APOE.predict<-imputed.APOE.b3.all.noNA[match(RM.map.APOE[,1],names(imputed.APOE.b3.all.noNA))]

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
# 5715
7469-1754
# 5715


## Rerun ensemble ML with more complete known APOE genotypes for training

training.cleanDat.noNA<-t(na.omit(t(cleanDat.4p13b3[which(!rownames(cleanDat.4p13b3) %in% c("IRF6|O14896")), which(!is.na(gt.APOE))])))
dim(training.cleanDat.noNA)
#  7334 13004   # 7335 13004 with IRF6 in (site K does not have APOE genotypes) -- otherwise it would be 13410 samples for training
training.gt.APOE<-gt.APOE[colnames(training.cleanDat.noNA)]
predict_APOE.b3.allAndMapped.noNA<-fit_APOE_ensemble_par(t(training.cleanDat.noNA), training.gt.APOE, ncores=14)
# Impute APOE epsilon 6 genotypes in the 7k samples missing it (step 2 of 2)
imputed.APOE.b3.allAndMapped.noNA<-predict_APOE.b3.allAndMapped.noNA(t(impute::impute.knn(cleanDat.4p13b3)$data[rownames(training.cleanDat.noNA),which(is.na(gt.APOE))]))
table(imputed.APOE.b3.allAndMapped.noNA)

names(imputed.APOE.b3.allAndMapped.noNA)<-rownames(t(cleanDat.4p13b3[rownames(training.cleanDat.noNA),which(is.na(gt.APOE))]))
saveRDS(predict_APOE.b3.allAndMapped.noNA,"predict_APOE.b3.allAndMapped.noNA.RDS")

saveRDS(imputed.APOE.b3.allAndMapped.noNA,"imputed.APOE.b3.allAndMapped.noNA.RDS")

#predict_APOE.b3.all.noNA<-readRDS("predict_APOE.b3.allAndMapped.noNA.RDS")


## Predict APOE e3/e3 and e4/e4 with very high accuracy in a modified ensemble ML prediction function

## VIM kNN disabled (impute input to the function)
##############################################################################
##  Parallel fit_APOE_ensemble()  -------------------------------------------
##############################################################################
library(caret)
library(glmnet)
library(xgboost)
library(ranger)
library(progressr)
#library(VIM)
library(dplyr)
#library(doParallel)
library(future)
library(doFuture)
library(doRNG)

fit_APOE_ensemble_homo33_44maxAcc <- function(expr, APOE_gt,
                                  nfold = 5, nrep = 5,
                                  ncores = parallel::detectCores() - 1,
                                  seed   = 1) {

  set.seed(seed)

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

  plan(multisession, workers = ncores)  # or multicore on linux/macOS
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
                     .export = "p",            # let workers see 'p'
                     .packages  = c("progressr","glmnet","xgboost","ranger","dplyr")) %dorng% {

    ## helper function - make all rectangular matrices carry all levels/genotypes
    complete_prob_mat <- function(mat, all_levels) {
      miss <- setdiff(all_levels, colnames(mat))
      if (length(miss)) {
        mat <- cbind(mat, matrix(0, nrow(mat), length(miss),
                                 dimnames = list(NULL, miss)))
      }
      mat[ , all_levels, drop = FALSE]           # reorder & keep only wanted
    }

    ## ----------------- announce fold start -----------------------------
    p(sprintf("fold %d/%d  -  started", fold, length(cvIndex)), amount=0)

    set.seed(seed + fold)                        # reproducible inside worker
    tr_idx <- cvIndex[[fold]]
    te_idx <- setdiff(seq_len(nrow(expr)), tr_idx)

    ## ---------- preprocessing  (no kNN, already imputed) --------------
    prep <- function(m) scale(m[, colMeans(is.na(m)) <= .20, drop = FALSE])
    X_tr <- prep(expr[tr_idx, ])
    X_te <- scale(expr[te_idx, colnames(X_tr)],
                  center = attr(X_tr, "scaled:center"),
                  scale  = attr(X_tr, "scaled:scale"))
    X_te[is.na(X_te)] <- 0

    levels_all <- sort(unique(APOE_gt))            # six possible APOE genotypes
    num_class   <- length(levels_all)              # == 6

    y_tr <- factor(APOE_gt[tr_idx], levels = levels_all)
    y_te <- factor(APOE_gt[te_idx], levels = levels_all)

    ## NEW-weights: one vector that all learners can consume
    class_w <- c("e2/e2" = 1, "e2/e3" = 1,
                 "e3/e3" = 2,                # want *very* high precision
                 "e2/e4" = 1, "e3/e4" = 1,
                 "e4/e4" = 8)                # low prevalence -> strong weight

    #For robustness, drop samples and labels whose labels are NA (should not have been input)
    keep_tr <- !is.na(y_tr)
    keep_te <- !is.na(y_te)

    X_tr <- X_tr[keep_tr, , drop = FALSE]
    y_tr <- y_tr[keep_tr]

    X_te <- X_te[keep_te, , drop = FALSE]
    y_te <- y_te[keep_te]

    ## ---------------- glmnet ---------------------------------------------
    w_tr <- class_w[as.character(y_tr)]
    cv_glm <- cv.glmnet(X_tr, y_tr,
                        family = "multinomial",
                        weights = w_tr,
                        type.measure = "class", parallel = FALSE)
#    cv_glm <- cv.glmnet(X_tr, y_tr, family = "multinomial",
#                        type.measure = "class", parallel = FALSE)
    p(sprintf("fold %d/%d  •  glmnet done", fold, length(cvIndex)), amount=1)

    p_glm  <- predict(cv_glm, X_te, s = "lambda.min",
                      type = "response")[,,1]
    p_glm <- complete_prob_mat(p_glm, levels_all)

    ## ---------------- xgboost --------------------------------------------
#    dtr <- xgb.DMatrix(X_tr, label = as.numeric(y_tr) - 1)
    dtr  <- xgb.DMatrix(X_tr, label = as.numeric(y_tr) - 1, weight = w_tr)
    dval <- xgb.DMatrix(X_te, label = as.numeric(y_te) - 1)

    xpar <- list(eta = 0.1, max_depth = 6, subsample = 0.8,
                 colsample_bytree = 0.8,
                 objective = "multi:softprob",
                 num_class = num_class,     # length(levels(y_tr)),
                 nthread = 1)               # <- single thread per worker

    bst  <- xgb.train(params=xpar, data=dtr, watchlist=list(train=dtr, eval=dval), nrounds = 200,
                      verbose = 0, early_stopping_rounds = 20)
    p_xgb <- matrix(predict(bst, X_te),
                    ncol = num_class, byrow = TRUE)
#                    ncol = length(levels(y_tr)), byrow = TRUE)
#    colnames(p_xgb) <- levels(y_tr)
    colnames(p_xgb) <- levels_all                 # <- name first
    p_xgb <- complete_prob_mat(p_xgb, levels_all) # <- then complete
    p(sprintf("fold %d/%d  •  xgboost done", fold, length(cvIndex)), amount=1)

    ## ---------------- random-forest --------------------------------------
#    rf  <- ranger(y_tr ~ ., data = data.frame(y_tr, X_tr),
#                  probability = TRUE, num.trees = 500,
#                  num.threads = 1)          # <- single thread
    rf <- ranger(y_tr ~ ., data = data.frame(y_tr, X_tr),
                 probability   = TRUE,
                 num.trees     = 500, num.threads = 1,
                 class.weights = class_w)

    p(sprintf("fold %d/%d  •  ranger done", fold, length(cvIndex)), amount=1)

    p_rf <- predict(rf, data.frame(X_te))$predictions
    colnames(p_rf) <- rf$forest$levels
    p_rf  <- complete_prob_mat(p_rf, levels_all)

    ## ---------------- ensemble vote --------------------------------------

    ## ----- averaged class-probabilities --------------------------------------
    p_avg <- (p_glm + p_xgb + p_rf) / 3          # n x 6  matrix

    top_class <- apply(p_avg, 1, function(z) names(z)[which.max(z)])
    top_prob  <- apply(p_avg, 1, max)
    second    <- apply(p_avg, 1, function(z) sort(z, TRUE)[2])

    ## confidence rules --------------------------------------------------------
    ##   1. if the winner is e3/e3 or e4/e4  AND  its prob >= 0.80  ? keep it
    ##   2. if the winner is anything else   AND  (prob >= 0.90  *and*
    ##                                             prob-margin >= 0.20) ? keep it
    ##   3. otherwise                                                             ? NA
    keep_homo  <- top_class %in% c("e3/e3", "e4/e4") & top_prob >= 0.80
    keep_other <- !(top_class %in% c("e3/e3", "e4/e4")) &
                  top_prob >= 0.90 & (top_prob - second) >= 0.20

    y_hat <- ifelse(keep_homo | keep_other, top_class, NA_character_)
    y_hat <- factor(y_hat, levels = levels_all)
#    p_avg <- (p_glm + p_xgb + p_rf) / 3
#    y_hat <- factor(colnames(p_avg)[max.col(p_avg)],
#                    levels = levels(y_tr))

    cm <- tryCatch(
             caret::confusionMatrix(y_hat, y_te),
             error = function(e) NULL)

## ---------------- safe extractor --------------------------------------
    safe_metric <- function(cm, cls, what) {
      if (is.null(cm)) return(NA_real_)
      bc <- cm$byClass
      if (is.null(dim(bc))) return(NA_real_)        # 2-class vector - skip
      if (cls %in% rownames(bc) && what %in% colnames(bc))
          return(bc[cls, what])
      NA_real_
    }

    prec_e33 <- safe_metric(cm, "e3/e3", "Precision")
    rec_e33  <- safe_metric(cm, "e3/e3", "Recall")
    prec_e44 <- safe_metric(cm, "e4/e4", "Precision")
    rec_e44  <- safe_metric(cm, "e4/e4", "Recall")

    acc   <- if (!is.null(cm)) cm$overall["Accuracy"] else NA_real_
    macro <- if (!is.null(cm) && !is.null(dim(cm$byClass)))
                 mean(cm$byClass[,"F1"]) else NA_real_

    data.frame(Accuracy = acc,
               MacroF1  = macro,
               Prec_e33 = prec_e33,  Rec_e33 = rec_e33,
               Prec_e44 = prec_e44,  Rec_e44 = rec_e44,
               Fold     = fold)
  } # foreach
  }) # with_progress

#  stopCluster(cl)                               # tidy up
#  registerDoSEQ()                               # back to sequential

  cat(sprintf("CV Accuracy  %.3f ± %.3f\n",
              mean(metrics$Accuracy, na.rm=TRUE), sd(metrics$Accuracy, na.rm=TRUE)))
  cat(sprintf("CV Macro-F1  %.3f ± %.3f\n\n",
              mean(metrics$MacroF1, na.rm=TRUE),  sd(metrics$MacroF1, na.rm=TRUE)))
  cat(sprintf("Precision (e3/e3) %.3f; Recall %.3f\n",
             mean(metrics$Prec_e33, na.rm=TRUE), mean(metrics$Rec_e33, na.rm=TRUE)))
  cat(sprintf("Precision (e4/e4) %.3f; Recall %.3f\n\n",
             mean(metrics$Prec_e44, na.rm=TRUE), mean(metrics$Rec_e44, na.rm=TRUE)))

  ## ------------------------------------------------------------------------
  ## 3.  fit final ensemble on all data  (sequential) -----------------------
  ## ------------------------------------------------------------------------
  prep_expr <- function(mat, thr_miss = 0.20) {
    keep <- which(colMeans(is.na(mat)) <= thr_miss)
    mat  <- mat[, keep, drop = FALSE]
#    mat  <- VIM::kNN(mat, k = 5, imp_var = FALSE)
    scale(mat)
  }

  with_progress({
    p2 <- progressor(steps=3)
    p2(sprintf("Starting fit of final ensemble on all data (serial/non-parallel)..."), amount=0)

  class_w <- c("e2/e2" = 1, "e2/e3" = 1,
               "e3/e3" = 2,                # want *very* high precision
               "e2/e4" = 1, "e3/e4" = 1,
               "e4/e4" = 8)                # low prevalence -> strong weight

  levels_all <- sort(unique(APOE_gt))

  X_all <- prep_expr(expr)
  y_all <- factor(APOE_gt, levels=levels_all)

  final_glm <- cv.glmnet(X_all, y_all, family = "multinomial", weights = class_w[as.character(y_all)],
                         type.measure = "class", parallel = FALSE)
    p2(sprintf("GLM fit on full data  •  finished"), amount=1)

  xpar_all <- list(eta = 0.1, max_depth = 6, subsample = 0.8,
                   colsample_bytree = 0.8,
                   objective = "multi:softprob",
                   num_class = length(levels(y_all)),
                   nthread = ncores)          # can use all cores now

  w_all    <- class_w[as.character(y_all)]
  d.all    <- xgb.DMatrix(X_all, label = as.numeric(y_all)-1, weight = w_all)
  bst_all  <- xgb.train(params=xpar_all, data=d.all,
                        watchlist=list(train=d.all),
                        nrounds = 200, verbose = 0,   #instead of 200, could use best #rounds: bst$best_iteration
                        early_stopping_rounds = 20)
    p2(sprintf("XGboost fit on full data  •  finished"), amount=1)

  rf_all <- ranger(y_all ~ ., data = data.frame(y_all, X_all), class.weights = class_w,
                   probability = TRUE, num.trees = 500,
                   num.threads = ncores)
    p2(sprintf("Random Forest fit on full data  •  finished"), amount=1)

    p2(sprintf("Fit on full data  •  finished"), amount=0)

  ## ------------------------------------------------------------------------
  ## 4.  prediction wrapper -------------------------------------------------
  ## ------------------------------------------------------------------------
  function(new_expr) {
    new_expr <- prep_expr(rbind(expr[1,,drop=FALSE], new_expr))[-1,,drop=FALSE]
    new_expr <- new_expr[, colnames(X_all), drop = FALSE]

    ## helper - makes every probability matrix  n x 6
    complete_prob_mat <- function(mat, all_levels) {
      miss <- setdiff(all_levels, colnames(mat))
      if (length(miss)) {
        mat <- cbind(mat, matrix(0, nrow(mat), length(miss),
                                 dimnames = list(NULL, miss)))
      }
      mat[ , all_levels, drop = FALSE]
    }

    ## ---- individual learners -------------------------------------------------
    p1 <- predict(final_glm, new_expr, s = "lambda.min",
                  type = "response")[,,1]
    p1 <- complete_prob_mat(p1, levels_all)

    p2_mat <- matrix(predict(bst_all, new_expr),
                     ncol = length(levels_all), byrow = TRUE)
    colnames(p2_mat) <- levels_all
    p2_mat <- complete_prob_mat(p2_mat, levels_all)

    p3_mat <- predict(rf_all, data.frame(new_expr))$predictions
    colnames(p3_mat) <- rf_all$forest$levels
    p3_mat <- complete_prob_mat(p3_mat, levels_all)

    ## ---- averaged probabilities & confidence filter --------------------------
    p_avg <- (p1 + p2_mat + p3_mat) / 3      # n x 6

    top_class <- apply(p_avg, 1, function(z) names(z)[which.max(z)])
    top_prob  <- apply(p_avg, 1, max)
    second    <- apply(p_avg, 1, function(z) sort(z, TRUE)[2])

    keep_homo  <- top_class %in% c("e3/e3", "e4/e4") & top_prob >= 0.80
    keep_other <- !(top_class %in% c("e3/e3", "e4/e4")) &
                  top_prob >= 0.90 & (top_prob - second) >= 0.20

    pred <- ifelse(keep_homo | keep_other, top_class, NA_character_)
    factor(pred, levels = levels_all)
    }
  }) # with_progress
}


################[QUICK TEST RUN - 28 THREADS OK ON SMALL TRAINING SET]####################
#option 3 - impute::impute.knn and reduce input expr mat size
training.cleanDat <- impute::impute.knn(cleanDat.4p13b3[, which(!is.na(gt.APOE))])$data   #  no NA
#sample.subset<-sort(sample(length(na.omit(gt.APOE)),750))
#table(na.omit(gt.APOE)[sample.subset])
set.seed(1)
# sample 150 of each genotype for training (n=810 total)
keep_idx<-unlist(lapply(split(seq_along( na.omit(gt.APOE)), na.omit(gt.APOE)), function(i) if (length(i) > 150) sample(i,150) else i), use.names=FALSE)
table(na.omit(gt.APOE)[keep_idx])  # 60 2/2, 150 ea of the other 5 genotypes; took <40 min to complete training on 28/32 cores at 3.2GHz (128 GB RAM rec.)
predict_APOE.b3.allAndMapped.e33e44maxAcc<-fit_APOE_ensemble_homo33_44maxAcc(t(training.cleanDat[APOE4.assays.keep[c(1:250,2251:length(APOE4.assays.keep))], keep_idx]), na.omit(gt.APOE)[keep_idx], ncores=28)

# Impute APOE epsilon 6 genotypes in the 7k samples missing it (step 2 of 2)
imputed.APOE.b3.allAndMapped.e33e44maxAcc<-predict_APOE.b3.allAndMapped.e33e44maxAcc(t(impute::impute.knn(cleanDat.4p13b3)$data[APOE4.assays.keep[c(1:250,2251:length(APOE4.assays.keep))], which(is.na(gt.APOE))]))
################[QUICK TEST RUN - 28 THREADS - 	COMPLETED]####################
table(imputed.APOE.b3.allAndMapped.e33e44maxAcc)
#over prediction of e4/e4, only 419 total predictions non-NA


## Rerun ensemble ML to get very high accuracy calls predicting e4/e4 and e3/e3 in unknown APOE genotype samples

# done above
#training.cleanDat.noNA<-t(na.omit(t(cleanDat.4p13b3[which(!rownames(cleanDat.4p13b3) %in% c("IRF6|O14896")), which(!is.na(gt.APOE))])))
#dim(training.cleanDat.noNA)
##  7334 13004   # 7335 13004 with IRF6 in (site K does not have APOE genotypes) -- otherwise it would be 13410 samples for training
#training.gt.APOE<-gt.APOE[colnames(training.cleanDat.noNA)]
predict_APOE.b3.allAndMapped.e33e44maxAcc<-fit_APOE_ensemble_homo33_44maxAcc(t(training.cleanDat.noNA), training.gt.APOE, ncores=14)
#CV Accuracy  0.972 ± 0.005
#CV Macro-F1  NaN ± NA
#
#Precision (e3/e3) NaN; Recall NaN
#Precision (e4/e4) NaN; Recall NaN
# Impute APOE epsilon 6 genotypes in the 7k samples missing it (step 2 of 2)
imputed.APOE.b3.allAndMapped.e33e44maxAcc<-predict_APOE.b3.allAndMapped.e33e44maxAcc(t(impute::impute.knn(cleanDat.4p13b3)$data[rownames(training.cleanDat.noNA),which(is.na(gt.APOE))]))
table(imputed.APOE.b3.allAndMapped.e33e44maxAcc)
#e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#    0     0     0  2655     0    18
#compared to before
table(imputed.APOE.b3.allAndMapped.noNA)
#e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#   12   578   101  3229  1625   170

names(imputed.APOE.b3.allAndMapped.e33e44maxAcc)<-rownames(t(cleanDat.4p13b3[rownames(training.cleanDat.noNA),which(is.na(gt.APOE))]))
saveRDS(predict_APOE.b3.allAndMapped.e33e44maxAcc,"predict_APOE.b3.allAndMapped.e33e44maxAcc.RDS")

saveRDS(imputed.APOE.b3.allAndMapped.e33e44maxAcc,"imputed.APOE.b3.allAndMapped.e33e44maxAcc.RDS")

#predict_APOE.b3.all.e33e44maxAcc<-readRDS("predict_APOE.b3.allAndMapped.e33e44maxAcc.RDS")


## Below is a minimal-change recipe that:
#gets rid of the last “subscript out of bounds”,
#adds two one-vs-rest (OvR) ensembles - one for e3/e3, one for e4/e4,
##automatically finds the lowest probability cut-off that still gives ~98% precision on the 15 k labelled subjects,
#keeps the current multi-class model for the other four genotypes,
#plugs the new high-confidence logic into the prediction wrapper.

# EDITED:
#remove bias toward e3/e3 (weight=1) so that the 6-way ensemble does not prefer e3/e3 when unsure
#choose the one-vs-rest cutoff by precision rather than a quantile; now a row is promoted only if the training set shows that >=96% of rows above that probability really are the homozygote.
#add a margin for homozygotes (keep_homo); with a 0.10 gap the ensemble must still be clearly "more sure" of a homozygote than of any other class.


## VIM kNN disabled (impute input to the function)
##############################################################################
##  Parallel fit_APOE_ensemble()  -------------------------------------------
##############################################################################
library(caret)
library(glmnet)
library(xgboost)
library(ranger)
library(progressr)
#library(VIM)
library(dplyr)
#library(doParallel)
library(future)
library(doFuture)
library(doRNG)

fit_APOE_ensemble_homo33_44.96Acc <- function(expr, APOE_gt,
                                  nfold = 5, nrep = 5,
                                  ncores = parallel::detectCores() - 1,
                                  target=c("e4/e4"),   # run one genotype vs all others
                                  target_ppv=0.96,
                                  seed   = 1) {

  set.seed(seed)

  memLimit=4*1024^3
  options(future.globals.maxSize= memLimit)  #4GB Total size of all global objects that need to be exported - up from 500MB
  Sys.setenv(R_FUTURE_GLOBALS_MAXSIZE=memLimit) #inherited by workers

  ## Helper function - threshold minimum
  pick_thr <- function(prob, truth, target_ppv = 0.96,
                       min_tp = 30, floor = 0.80) {
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

  plan(multisession, workers = ncores, globals.maxSize=memLimit)  # or multicore on linux/macOS
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
                     .options.future = list(  #expr not exported explicitly; captured only once.
                        globals = list(APOE_gt  = APOE_gt,
                                       pick_thr = pick_thr,
                                       target   = target,
                                       seed     = seed)
                     ),
                     .export = "p",            # let workers see 'p'
                     .packages  = c("progressr","glmnet","xgboost","ranger","dplyr")) %dorng% {

    ## ----------------- announce fold start -----------------------------
    p(sprintf("fold %d/%d  -  started", fold, length(cvIndex)), amount=0)

    set.seed(seed + fold)                        # reproducible inside worker

    tr <- cvIndex[[fold]]
    te <- setdiff(seq_len(nrow(expr)), tr)

    tr_idx <- cvIndex[[fold]]
    te_idx <- setdiff(seq_len(nrow(expr)), tr_idx)

    ## ---------- preprocessing  (no kNN, already imputed) --------------
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
                 nthread = 1)
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
    p_avg <- (p_glm + p_xgb + p_rf) / 3
    thr   <- pick_thr(p_avg, y_bin_tr == "pos", target_ppv)   # use helper

    pred  <- factor(ifelse(p_avg >= thr, target, NA_character_),
                    levels = c(target))

    tp <- sum(pred == target & y_bin_te == "pos", na.rm = TRUE)
    fp <- sum(pred == target & y_bin_te == "neg", na.rm = TRUE)
    fn <- sum(pred != target & y_bin_te == "pos", na.rm = TRUE)

    prec <- if (tp+fp) tp/(tp+fp) else NA
    rec  <- if (tp+fn) tp/(tp+fn) else NA

    data.frame(Precision = prec, Recall = rec, Fold = fold)
  } # foreach
  }) # with_progress

#  stopCluster(cl)                               # tidy up
#  registerDoSEQ()                               # back to sequential

  cat(sprintf("CV %s - Precision %.3f ± %.3f | Recall %.3f ± %.3f\n\n",
              target,
              mean(metrics$Precision, na.rm=TRUE), sd(metrics$Precision, na.rm=TRUE),
              mean(metrics$Recall,    na.rm=TRUE), sd(metrics$Recall,    na.rm=TRUE)))

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
      w_pos  <- if (tg == "e4/e4") 12 else 8           # keep old weights
      w_bin  <- ifelse(y_bin == "pos", w_pos, 1)

      ## ---- (i) glmnet --------------------------------------------------------
      glm_tg <- cv.glmnet(X_all, y_bin, family = "binomial",
                          weights = w_bin, type.measure = "class")
      p2(sprintf("GLM fit on full data  •  finished"), amount=1)

      ## ---- (ii) xgboost ------------------------------------------------------
      d_bin  <- xgb.DMatrix(X_all, label = as.numeric(y_bin) - 1, weight = w_bin)
      xpar_b <- list(objective = "binary:logistic", eta = 0.1,
                     max_depth = 6, subsample = 0.8,
                     colsample_bytree = 0.8, nthread = ncores)
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

  ## ------------------------------------------------------------------------
  ## 4.  prediction wrapper -------------------------------------------------
  ## ------------------------------------------------------------------------
    function(new_expr) {

      new_expr <- prep_expr(rbind(expr[1,,drop=FALSE], new_expr))[-1,,drop=FALSE]
      new_expr <- new_expr[, colnames(X_all), drop = FALSE]

      glm_tg <- ovr[[target]]$glm
      xgb_tg <- ovr[[target]]$xgb
      rf_tg  <- ovr[[target]]$rf
      thr    <- ovr[[target]]$thr

      p_g <- drop(predict(glm_tg, new_expr, s = "lambda.min", type="response"))
      p_x <- drop(predict(xgb_tg, new_expr))
      p_r <- predict(rf_tg, data.frame(new_expr))$predictions[,"pos"]

      p_bin <- (p_g + p_x + p_r) / 3

      factor(ifelse(p_bin >= thr, target, NA_character_), levels = target)
    }
  }) # with_progress
}


################[QUICK TEST RUN - 28 THREADS OK ON SMALL TRAINING SET]####################
#option 3 - impute::impute.knn and reduce input expr mat size
training.cleanDat <- impute::impute.knn(cleanDat.4p13b3[, which(!is.na(gt.APOE))])$data   #  no NA
#sample.subset<-sort(sample(length(na.omit(gt.APOE)),750))
#table(na.omit(gt.APOE)[sample.subset])
set.seed(1)
# sample 400 of each genotype for training (n=810 total)
keep_idx<-unlist(lapply(split(seq_along( na.omit(gt.APOE)), na.omit(gt.APOE)), function(i) if (length(i) > 400) sample(i,400) else i), use.names=FALSE)
table(na.omit(gt.APOE)[keep_idx])  # 66 2/2, now 400 ea of the other 5 genotypes; took <5 min to complete training on 28/32 cores at 3.2GHz (128 GB RAM rec.)
predict_APOE.b3.allAndMapped.e44only.96Acc<-fit_APOE_ensemble_homo33_44.96Acc(t(training.cleanDat[APOE4.assays.keep[c(1:250,2251:length(APOE4.assays.keep))], keep_idx]), na.omit(gt.APOE)[keep_idx], ncores=28, target="e4/e4")
#CV Accuracy  0.976 ± 0.013
#CV Macro-F1  0.977 ± 0.022
#
# e3/e3  -  Precision 0.952 | Recall 0.003
# e4/e4  -  Precision 0.983 | Recall 0.997
#
# Impute APOE epsilon 6 genotypes in the 7k samples missing it (step 2 of 2)
imputed.APOE.b3.allAndMapped.e44only.96Acc<-predict_APOE.b3.allAndMapped.e44only.96Acc(t(impute::impute.knn(cleanDat.4p13b3)$data[APOE4.assays.keep[c(1:250,2251:length(APOE4.assays.keep))], which(is.na(gt.APOE))]))
################[QUICK TEST RUN - 28 THREADS - 	COMPLETED]####################
table(imputed.APOE.b3.allAndMapped.e44only.96Acc)
#over prediction of e4/e4 (445; exptected 170 or less)


## Rerun ensemble ML to get very high accuracy calls predicting e4/e4 (binary) in unknown APOE genotype samples

# done above
#training.cleanDat.noNA<-t(na.omit(t(cleanDat.4p13b3[which(!rownames(cleanDat.4p13b3) %in% c("IRF6|O14896")), which(!is.na(gt.APOE))])))
#dim(training.cleanDat.noNA)
##  7334 13004   # 7335 13004 with IRF6 in (site K does not have APOE genotypes) -- otherwise it would be 13410 samples for training
#training.gt.APOE<-gt.APOE[colnames(training.cleanDat.noNA)]
predict_APOE.b3.allAndMapped.e44only.96Acc<-fit_APOE_ensemble_homo33_44.96Acc(t(training.cleanDat.noNA), training.gt.APOE, ncores=14)
#CV e4/e4 - Precision 0.972 ± 0.025 | Recall 1.000 +/- 0.000

# Impute APOE epsilon 6 genotypes in the 5715 (5k) samples missing it (step 2 of 2)
imputed.APOE.b3.allAndMapped.e44only.96Acc<-predict_APOE.b3.allAndMapped.e44only.96Acc(t(impute::impute.knn(cleanDat.4p13b3)$data[rownames(training.cleanDat.noNA),which(is.na(gt.APOE))]))
table(imputed.APOE.b3.allAndMapped.e44only.96Acc)
#current run
#e4/e4
#   70
#compared to before
table(imputed.APOE.b3.allAndMapped.noNA)
#e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#   12   578   101  3229  1625   170

names(imputed.APOE.b3.allAndMapped.e44only.96Acc)<-rownames(t(cleanDat.4p13b3[rownames(training.cleanDat.noNA),which(is.na(gt.APOE))]))
saveRDS(predict_APOE.b3.allAndMapped.e44only.96Acc,"predict_APOE.b3.allAndMapped.e44only.96Acc.RDS")

saveRDS(imputed.APOE.b3.allAndMapped.e44only.96Acc,"imputed.APOE.b3.allAndMapped.e44only.96Acc.RDS")

#predict_APOE.b3.all.e33e44.96Acc<-readRDS("predict_APOE.b3.allAndMapped.e33e44.96Acc.RDS")


## Rerun ensemble ML to get very high accuracy calls predicting e3/e3 (binary) in unknown APOE genotype samples

# done above
#training.cleanDat.noNA<-t(na.omit(t(cleanDat.4p13b3[which(!rownames(cleanDat.4p13b3) %in% c("IRF6|O14896")), which(!is.na(gt.APOE))])))
#dim(training.cleanDat.noNA)
##  7334 13004   # 7335 13004 with IRF6 in (site K does not have APOE genotypes) -- otherwise it would be 13410 samples for training
#training.gt.APOE<-gt.APOE[colnames(training.cleanDat.noNA)]
predict_APOE.b3.allAndMapped.e33only.96Acc<-fit_APOE_ensemble_homo33_44.96Acc(t(training.cleanDat.noNA), training.gt.APOE, ncores=14, target="e3/e3")
#CV e3/e3 - Precision 0.986 ± 0.006 | Recall 1.000 +/- 0.000

# Impute APOE epsilon 1 genotype in the 5715 (5k) samples missing genotype (step 2 of 2)
imputed.APOE.b3.allAndMapped.e33only.96Acc<-predict_APOE.b3.allAndMapped.e33only.96Acc(t(impute::impute.knn(cleanDat.4p13b3)$data[rownames(training.cleanDat.noNA),which(is.na(gt.APOE))]))
table(imputed.APOE.b3.allAndMapped.e33only.96Acc)
#current run
#e3/e3
# 2840
#compared to before
table(imputed.APOE.b3.allAndMapped.noNA)
#e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#   12   578   101  3229  1625   170

names(imputed.APOE.b3.allAndMapped.e33only.96Acc)<-rownames(t(cleanDat.4p13b3[rownames(training.cleanDat.noNA),which(is.na(gt.APOE))]))
saveRDS(predict_APOE.b3.allAndMapped.e33only.96Acc,"predict_APOE.b3.allAndMapped.e33only.96Acc.RDS")

saveRDS(imputed.APOE.b3.allAndMapped.e33only.96Acc,"imputed.APOE.b3.allAndMapped.e33only.96Acc.RDS")

#predict_APOE.b3.allAndMapped.e33only.96Acc<-readRDS("predict_APOE.b3.allAndMapped.e33only.96Acc.RDS")


# sanity check - cases now predicted as e3/e3 were what in the noNA imputed vector?
table(imputed.APOE.b3.allAndMapped.noNA[which(imputed.APOE.b3.allAndMapped.e33only.96Acc=="e3/e3")])
#e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#    0     0     0  2840     0     0
# and the same question for cases now predicted as e4/e4 in the binary ensemble prediction:
table(imputed.APOE.b3.allAndMapped.noNA[which(imputed.APOE.b3.allAndMapped.e44only.96Acc=="e4/e4")])
#e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#    0     0     0     0     0    70

## Merge the two binary vectors, keeping e4/e4 if both are predicted for the same sample (should not occur)
imputed.APOE.b3.allAndMapped.binary.e33e44.96Acc<-imputed.APOE.b3.allAndMapped.e33only.96Acc
levels(imputed.APOE.b3.allAndMapped.binary.e33e44.96Acc)<-c("e3/e3","e4/e4")
imputed.APOE.b3.allAndMapped.binary.e33e44.96Acc[which(imputed.APOE.b3.allAndMapped.e44only.96Acc=="e4/e4")]<-"e4/e4"
table(imputed.APOE.b3.allAndMapped.binary.e33e44.96Acc) #should be no overlap, so 70 e4/e4 AND the e3/e3 sample count just predicted with very high confidence and Precision.
#e3/e3 e4/e4
# 2840    70

## Decision point: we will not use imputation for e4/e4 and e3/e3 96% enforced accuracy; use ground truth OR all 6 genotypes imputed.
## (But we will do nested fold CV with 20% hold-out to prove we are not overfitting, and also report our top predictive features, and consider removing the bottom least predictive ones
##  from input.)
################################

################################
## First full APOE genotype vector - sanity check

## Assemble predictions with known (ground truth) genotypes for a full nonmissing character/string vector of 22392 genotypes, in order.
APOE.noNA<-gt.APOE
APOE.noNA[names(imputed.APOE.b3.allAndMapped.noNA)]<-as.character(imputed.APOE.b3.allAndMapped.noNA)
table(APOE.noNA)  ## Full population (22392), with predicted genotypes (5715)
#e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#   78  2214   506 11719  6796  1079

## Sanity check using Better histogram leveraging top 6 e4 carrier correlates (3 + and 3 -)
SPC.med=median(cleanDat.4p13b3["SPC25|Q9HBM1",],na.rm=T)
CTF1.med=median(cleanDat.4p13b3["CTF1|Q16619^SL002783@seq.13732.79",],na.rm=T)
NEFL.med=median(cleanDat.4p13b3["NEFL|P07196",],na.rm=T)
S100A13.med=median(cleanDat.4p13b3["S100A13|Q99584",],na.rm=T)
hist.data.SPC=((cleanDat.4p13b3["SPC25|Q9HBM1",] - SPC.med) + (cleanDat.4p13b3["LRRN1|Q6UXK5^SL025922@seq.11293.14",] - LRRN1.med) + (cleanDat.4p13b3["CTF1|Q16619^SL002783@seq.13732.79",] - CTF1.med) - (cleanDat.4p13b3["NEFL|P07196",] - NEFL.med) - (cleanDat.4p13b3["S100A13|Q99584",] - S100A13.med) - (cleanDat.4p13b3["TBCA|O75347",] - TBCA.med))
hist(hist.data.SPC, breaks=100, xlab="log2(abundance ratio):  SPC+LRRN1+CTF1(median-centered) / TBCA+S100A13+NEFL(median-centered)", main="APOE e4 Carrier Status Best Nonmissing Correlate")
#hist(hist.data.SPC[which(numericMeta.reg.b345$APOE==24 | numericMeta.reg.b345$APOE==34)],breaks=60,col="#BBFFBB40",add=T)  #green overlay
#hist(hist.data.SPC[which(numericMeta.reg.b345$APOE==44)],breaks=70,col="#FFBBBB99",add=T)  #red overlay
#hist(hist.data.SPC[which(numericMeta.reg.b345$APOE==33)],breaks=100,col="#BBBBFF99",add=T)  #blue
## Full population (22392), with predicted genotypes (5715)
hist(hist.data.SPC[which(APOE.noNA=="e2/e4" | APOE.noNA=="e3/e4")],breaks=60,col="#BBFFBB40",add=T)  #green overlay
hist(hist.data.SPC[which(APOE.noNA=="e4/e4")],breaks=70,col="#FFBBBB99",add=T)  #red overlay
hist(hist.data.SPC[which(APOE.noNA=="e3/e3")],breaks=100,col="#BBBBFF99",add=T)  #blue

legend("topright",c("e4/e4","e3/e3", "e4 het"),fill=c("#FFBBBB99","#BBBBFF99","#BBFFBB40"))
abline(v=0.625,col="darkgreen",lty=2, lwd=2.2)
abline(v=4.25,col="maroon",lty=2, lwd=2.2)


## Just the imputed (predicted) samples - 5715
SPC.med=median(cleanDat.4p13b3["SPC25|Q9HBM1",],na.rm=T)
CTF1.med=median(cleanDat.4p13b3["CTF1|Q16619^SL002783@seq.13732.79",],na.rm=T)
NEFL.med=median(cleanDat.4p13b3["NEFL|P07196",],na.rm=T)
S100A13.med=median(cleanDat.4p13b3["S100A13|Q99584",],na.rm=T)
hist.data.SPC=((cleanDat.4p13b3["SPC25|Q9HBM1",] - SPC.med) + (cleanDat.4p13b3["LRRN1|Q6UXK5^SL025922@seq.11293.14",] - LRRN1.med) + (cleanDat.4p13b3["CTF1|Q16619^SL002783@seq.13732.79",] - CTF1.med) - (cleanDat.4p13b3["NEFL|P07196",] - NEFL.med) - (cleanDat.4p13b3["S100A13|Q99584",] - S100A13.med) - (cleanDat.4p13b3["TBCA|O75347",] - TBCA.med))
hist(hist.data.SPC, breaks=100, xlab="log2(abundance ratio):  SPC+LRRN1+CTF1(median-centered) / TBCA+S100A13+NEFL(median-centered)", main="[PREDICTED ONLY] APOE e4 Carrier Status Best Nonmissing Correlate")
#hist(hist.data.SPC[which(numericMeta.reg.b345$APOE==24 | numericMeta.reg.b345$APOE==34)],breaks=60,col="#BBFFBB40",add=T)  #green overlay
#hist(hist.data.SPC[which(numericMeta.reg.b345$APOE==44)],breaks=70,col="#FFBBBB99",add=T)  #red overlay
#hist(hist.data.SPC[which(numericMeta.reg.b345$APOE==33)],breaks=100,col="#BBBBFF99",add=T)  #blue
## Only the predicted (imputed) population (n=5715)
hist(hist.data.SPC[c(names(imputed.APOE.b3.allAndMapped.noNA)[which(imputed.APOE.b3.allAndMapped.noNA=="e2/e4")], names(imputed.APOE.b3.allAndMapped.noNA)[which(imputed.APOE.b3.allAndMapped.noNA=="e3/e4")])],breaks=60,col="#BBFFBB40",add=T)  #green overlay
hist(hist.data.SPC[names(imputed.APOE.b3.allAndMapped.noNA)[which(imputed.APOE.b3.allAndMapped.noNA=="e4/e4")] ],breaks=35,col="#FFBBBB99",add=T)  #red overlay
hist(hist.data.SPC[names(imputed.APOE.b3.allAndMapped.noNA)[which(imputed.APOE.b3.allAndMapped.noNA=="e3/e3")] ],breaks=50,col="#BBBBFF99",add=T)  #blue

legend("topright",c("e4/e4","e3/e3", "e4 het"),fill=c("#FFBBBB99","#BBBBFF99","#BBFFBB40"))
abline(v=0.625,col="darkgreen",lty=2, lwd=2.2)
abline(v=4.25,col="maroon",lty=2, lwd=2.2)
#################################


## Return ground truth APOE genotypes to known+mapped only (with 5715 NA)

gt.APOE<-rep(NA,length(numericMeta.reg.b345$APOE))
gt.APOE[which(numericMeta.reg.b345$APOE==22)]<-"e2/e2"
gt.APOE[which(numericMeta.reg.b345$APOE==23)]<-"e2/e3"
gt.APOE[which(numericMeta.reg.b345$APOE==33)]<-"e3/e3"
gt.APOE[which(numericMeta.reg.b345$APOE==24)]<-"e2/e4"
gt.APOE[which(numericMeta.reg.b345$APOE==34)]<-"e3/e4"
gt.APOE[which(numericMeta.reg.b345$APOE==44)]<-"e4/e4"
names(gt.APOE)<-rownames(numericMeta.reg.b345)

#UDS.map
UDS.map.APOE<-UDS.map
UDS.map.APOE$APOE.4cohort<-as.numeric(gsub("e","",traits.4cohort$APOE[match(UDS.map[,2],traits.4cohort$LoadedSampleName)]))
UDS.map.APOE$APOE.predict<-imputed.APOE.b3.all.noNA[match(UDS.map.APOE[,1],names(imputed.APOE.b3.all.noNA))]

UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE$APOE.4cohort==22)]<-"e2/e2"
UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE$APOE.4cohort==23)]<-"e2/e3"
UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE$APOE.4cohort==33)]<-"e3/e3"
UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE$APOE.4cohort==24)]<-"e2/e4"
UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE$APOE.4cohort==34)]<-"e3/e4"
UDS.map.APOE$APOE.4cohort[which(UDS.map.APOE$APOE.4cohort==44)]<-"e4/e4"

table(paste0(UDS.map.APOE$APOE.4cohort,".",UDS.map.APOE$APOE.predict))


BH.map.APOE<-BH.map
BH.map.APOE$APOE.4cohort<-as.numeric(gsub("e","",traits.4cohort$APOE[match(BH.map[,2],traits.4cohort$LoadedSampleName)]))
BH.map.APOE$APOE.predict<-imputed.APOE.b3.all.noNA[match(BH.map.APOE[,1],names(imputed.APOE.b3.all.noNA))]

BH.map.APOE$APOE.4cohort[which(BH.map.APOE$APOE.4cohort==22)]<-"e2/e2"
BH.map.APOE$APOE.4cohort[which(BH.map.APOE$APOE.4cohort==23)]<-"e2/e3"
BH.map.APOE$APOE.4cohort[which(BH.map.APOE$APOE.4cohort==33)]<-"e3/e3"
BH.map.APOE$APOE.4cohort[which(BH.map.APOE$APOE.4cohort==24)]<-"e2/e4"
BH.map.APOE$APOE.4cohort[which(BH.map.APOE$APOE.4cohort==34)]<-"e3/e4"
BH.map.APOE$APOE.4cohort[which(BH.map.APOE$APOE.4cohort==44)]<-"e4/e4"

table(paste0(BH.map.APOE$APOE.4cohort,".",BH.map.APOE$APOE.predict))


RM.map.APOE<-RM.map
RM.map.APOE$APOE.4cohort<-as.numeric(gsub("e","",traits.4cohort$APOE[match(RM.map[,2],traits.4cohort$LoadedSampleName)]))
RM.map.APOE$APOE.predict<-imputed.APOE.b3.all.noNA[match(RM.map.APOE[,1],names(imputed.APOE.b3.all.noNA))]

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
# 5715
7469-1754
# 5715

22392-5715
# 16677


## Output e3/e3 and e4/e4 ground truth samples with regression of site, protecting APOE e4 carrier status, age and sex in output data from within-site 2PAV protein regression

library("doParallel")
parallelThreads=31  #now Windows02  #max is number of processes that can run on your computer at one time
stopCluster(clusterLocal)
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="PSOCK")

registerDoParallel(clusterLocal)


####################### REGRESSION - b4c
  which(!names(gt.APOE)==colnames(cleanDat.unreg.b345))

  cleanDat.unreg.b4c <- cleanDat.unreg.b345[,which(!is.na(gt.APOE))]
  gt.APOE.b4c<-gt.APOE[colnames(cleanDat.unreg.b4c)]

  regvars.b4c<-numericMeta.reg.b345[which(!is.na(gt.APOE)),]
  regvars.b4c$Age<-as.numeric(numericMeta.reg.b345$age_at_visit[match(names(gt.APOE.b4c), rownames(numericMeta.reg.b345))])
  regvars.b4c$Sex.int=as.integer(abs(numericMeta.reg.b345$sex[match(names(gt.APOE.b4c), rownames(numericMeta.reg.b345))] -2))
  regvars.b4c$Sex<-relevel(factor(regvars.b4c$Sex.int), ref="0")
  regvars.b4c$APOE4.carrier<-relevel(factor(ifelse(gt.APOE.b4c %in% c("e2/e4","e3/e4","e4/e4"), 1, 0)), ref="0")
  regvars.b4c$APOE4.int=as.integer(regvars.b4c$APOE4.carrier)
  regvars.b4c$Site <-regvars.b4c$contributor_Fsplit
#  regvars.b4c<-na.omit(regvars.b4c)  # removes samples with missing age, sex
  dim(regvars.b4c)
  # 16677   104

  ##  Run the regression (4p13 b4c) - Site with Age+Sex+APOE.e4 (APOE4.carrier -- only ground truth known APOE genotypes included) protection

  ## covariate names you want to protect
  cov_keep <- c("Age", "Sex", "APOE4.carrier")   # Technically, we cleared out NAs, but one could list protection covariates that have missing values here

  normExpr.reg <- matrix(NA,nrow=nrow(cleanDat.unreg.b4c),ncol=ncol(cleanDat.unreg.b4c), dimnames = dimnames(cleanDat.unreg.b4c))
  good_samp=which(complete.cases(regvars.b4c[, cov_keep]))
  normExpr.reg <-  foreach (i=1:nrow(cleanDat.unreg.b4c), .combine=rbind, .packages="stats") %dopar% {  #** .combine=c, .multicombine=TRUE, .packages="stats", .export="regvars") %dopar% {
    y <- as.numeric(cleanDat.unreg.b4c[i, ])
    ## keep = samples where y and all covariates are present
    keep <- intersect(which(!is.na(y)), good_samp)             # columns that have a value, and regress with variables that have a value (we know Age and Sex are missing some)

    ## initialise result for this row as all NA
    adj <- rep(NA_real_, length(y))

    fit <- tryCatch(
        lm(y[keep]~Age+Sex+APOE4.carrier +Site, data = regvars.b4c[keep, , drop = FALSE]),
        error = function(e) NULL)

    if (!is.null(fit)) {
        coef <- coef(fit)
        ## coefficient[1] + residuals -> adjusted expression
        adj[keep] <- coef["(Intercept)"] + coef["Age"]*regvars.b4c[keep,"Age"] + coef[which(grepl("^Sex", names(coef)))]*regvars.b4c[keep,"Sex.int"] + coef["APOE4.carrier1"]*regvars.b4c[keep,"APOE4.int"] + residuals(fit)
    }

#        if (!is.null(fit)) {
#            ## build a model matrix with Site *zeroed out*
#            X     <- model.matrix(fit)
#            X[ , grep("^Site", colnames(X)) ] <- 0     # drop site contribution
#            y_hat <- as.numeric(X %*% coef(fit))       # fitted without Site
#
#            ## protected expression = y_hat + residuals  ( == y - Site effect )
#            adj[keep] <- y_hat + residuals(fit)
#        }
    adj                                    # returned to foreach
  }
  dimnames(normExpr.reg) <- dimnames(cleanDat.unreg.b4c)

  cat(paste0("Finished Pass 4p13 b4c regression of intersite variance.\n"))

# sanity check: some sites have missing values:
table(apply(normExpr.reg,2,function(x) length(which(is.na(x)))))
#    0  2303  <- # of missing values in column
#14758  2254  <- # of columns with that many NA

saveRDS(normExpr.reg,"4p13b4c.normExpr.reg_sites1-19_Fsplit_knownAPOE_only.RDS")
cleanDat.b4c<-normExpr.reg
#regvars.b4c still valid, ordered for this cleanDat
regvars.b4c$APOE.mapped<-gt.APOE.b4c
saveRDS(regvars.b4c, "4p13b4c.numericMeta_sites1-19_traits_knownAPOE_only.RDS")

## Sanity check
e4.bicor.to.siteCorr.b4c.assays<-bicor(t(cleanDat.b4c),regvars.b4c$APOE4.carrier, use='p')
# sort by bicor
e4.bicor.to.siteCorr.b4c.assays[order(unlist(t(e4.bicor.to.siteCorr.b4c.assays)),decreasing=TRUE),][c(1:5,7331:7335)]
#b4c (current)
#SPC25|Q9HBM1   LRRN1|Q6UXK5^SL025922@seq.11293.14   CTF1|Q16619^SL002783@seq.13732.79            NEFL|P07196        S100A13|Q99584                TBCA|O75347
#   0.7702713                            0.6725630                           0.4673459      ...    -0.5905291            -0.5929931                 -0.6264711

#b4b (previously)
#SPC25|Q9HBM1   LRRN1|Q6UXK5^SL025922@seq.11293.14   CTF1|Q16619^SL002783@seq.13732.79            NEFL|P07196        S100A13|Q99584                TBCA|O75347
#   0.7228271                            0.6734026                           0.4010664      ...    -0.5761990            -0.6053804                 -0.6564871

#b2 (previously)
#SPC25|Q9HBM1   LRRN1|Q6UXK5^SL025922@seq.11293.14   CTF1|Q16619^SL002783@seq.13732.79            NEFL|P07196        S100A13|Q99584                TBCA|O75347
#   0.7568582                            0.6478450                           0.4488890      ...    -0.5650291            -0.5672361                 -0.6022569


## QC (tSNE, VP, histogram)


##################################################
# 4p13b4c. 19 sites (no site U,V,W; site F split to F1, F2, F3) 1xPAVregr+sitewise Regress, Protect Age+Sex+APOE.e4 carrier status (binary, no imputation) QC (tSNE plots)

####################
## Examine tSNE of PLASMA - Human only assays (rows) and samples (columns);  no missing data by row.

numericMeta.plasma<-regvars.b4c
dim(numericMeta.plasma)
#16677  104
exprMat.plasma<-cleanDat.b4c

dim(as.data.frame(na.omit(exprMat.plasma)))
#  5032 16677   # 5031 in H; 7334 in K; 7335 in all others.
#exprMat.plasma<-exprMat.plasma[,which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))]
#numericMeta.plasma<-numericMeta.plasma[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma))))),]

#Group.3mappedCohortsPlusCTimputed alternate
Group=as.factor(regvars.b4c$Group.withCTimputed ) #[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))])


tSNE.list.plasma <- Rtsne::Rtsne(t(as.data.frame(na.omit(exprMat.plasma))),perplexity=20)
dim(tSNE.list.plasma$Y)
#16677   2
tSNE.plasma.xy<-as.data.frame(tSNE.list.plasma$Y)
colnames(tSNE.plasma.xy)<-c('x','y')
tSNE.4p13b4c.plasma.xy<-tSNE.plasma.xy
rownames(tSNE.4p13b4c.plasma.xy)<-rownames(regvars.b4c) #[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))]
library(ggplot2)
#library(ggpubr) - rlang upgrade required, ggplot2 upgrade required. cannot install from older source!
library(ggrepel)

# Get the indices of the first occurrence of each unique value in contributor_Fsplit
first_occurrence_indices.4p13b4c<-first_occurrence_indices <- match(unique(regvars.b4c$contributor_Fsplit), regvars.b4c$contributor_Fsplit)

labels.4p13b4c<-labels<-regvars.b4c$contributor_Fsplit[first_occurrence_indices.4p13b4c]

tSNE.plasma.samples.sites<-ggplot2::ggplot(tSNE.4p13b4c.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=regvars.b4c$contributor_Fsplit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13b4c.plasma.xy[first_occurrence_indices.4p13b4c, ],
                  aes(x=x,y=y, label = labels.4p13b4c), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.age<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$age_at_visit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleMatrix<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_matrix), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleType<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_type), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.apoe4dose<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$APOE4.Dose), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.MMSE<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$MMSE), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.Group<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=Group), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )


pdf(file="4p13b4c.2PAVproteinIntrasiteRegressed+Regress19sites_protectAgeSexAPOE.e4_carrierBinary_Fsplit.tSNE-Plasma(7335x16677)-samples_coloredByTraits.pdf",width=11,height=9)
  print(tSNE.plasma.samples.sites + labs(colour="Site"))
  print(tSNE.plasma.samples.age)
  print(tSNE.plasma.samples.sampleMatrix)
  print(tSNE.plasma.samples.sampleType)
  print(tSNE.plasma.samples.apoe4dose)
  print(tSNE.plasma.samples.MMSE)
  print(tSNE.plasma.samples.Group + labs(colour="Group"))
dev.off()


tSNE.4p13b4c.plasma.samples.sites<-tSNE.plasma.samples.sites


## Plot UNREGRESSED 16677
numericMeta.plasma<-regvars.b4c
dim(numericMeta.plasma)
#16677  104
exprMat.plasma<-cleanDat.unreg.b4c

dim(as.data.frame(na.omit(exprMat.plasma)))
#  5032 16677   # 5031 in H; 7334 in K; 7335 in all others.
#exprMat.plasma<-exprMat.plasma[,which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))]
#numericMeta.plasma<-numericMeta.plasma[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma))))),]

#Group.3mappedCohortsPlusCTimputed alternate
Group=as.factor(regvars.b4c$Group.withCTimputed ) #[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))])


tSNE.list.plasma <- Rtsne::Rtsne(t(as.data.frame(na.omit(exprMat.plasma))),perplexity=20)
dim(tSNE.list.plasma$Y)
#16677   2
tSNE.plasma.xy<-as.data.frame(tSNE.list.plasma$Y)
colnames(tSNE.plasma.xy)<-c('x','y')
tSNE.4p13b4c.plasma.xy<-tSNE.plasma.xy
rownames(tSNE.4p13b4c.plasma.xy)<-rownames(regvars.b4c) #[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))]
library(ggplot2)
#library(ggpubr) - rlang upgrade required, ggplot2 upgrade required. cannot install from older source!
library(ggrepel)

# Get the indices of the first occurrence of each unique value in contributor_Fsplit
first_occurrence_indices.4p13b4c<-first_occurrence_indices <- match(unique(regvars.b4c$contributor_Fsplit), regvars.b4c$contributor_Fsplit)

labels.4p13b4c<-labels<-regvars.b4c$contributor_Fsplit[first_occurrence_indices.4p13b4c]

tSNE.plasma.samples.sites<-ggplot2::ggplot(tSNE.4p13b4c.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=regvars.b4c$contributor_Fsplit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13b4c.plasma.xy[first_occurrence_indices.4p13b4c, ],
                  aes(x=x,y=y, label = labels.4p13b4c), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.age<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$age_at_visit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleMatrix<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_matrix), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleType<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_type), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.apoe4dose<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$APOE4.Dose), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.MMSE<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$MMSE), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.Group<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=Group), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )


pdf(file="4p13b4c.UNREGRESSED.2PAVproteinIntrasiteRegressed+Fsplit.tSNE-Plasma(7335x16677)-samples_coloredByTraits.pdf",width=11,height=9)
  print(tSNE.plasma.samples.sites + labs(colour="Site"))
  print(tSNE.plasma.samples.age)
  print(tSNE.plasma.samples.sampleMatrix)
  print(tSNE.plasma.samples.sampleType)
  print(tSNE.plasma.samples.apoe4dose)
  print(tSNE.plasma.samples.MMSE)
  print(tSNE.plasma.samples.Group + labs(colour="Group"))
dev.off()


tSNE.4p13b4c.UNREGRESSED.plasma.samples.sites<-tSNE.plasma.samples.sites


library("doParallel")
parallelThreads=30  #now Windows02  #max is number of processes that can run on your computer at one time
stopCluster(clusterLocal)
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="PSOCK")

registerDoParallel(clusterLocal)

#############################################################################
## 4p13b4c. Variance Partition regressed (2PAV regression intrasite)+Site regressed; protect Age+Sex+E4 carrier binary status (no NA, NA NOT imputed) (QC)

regvars.vp<-data.frame(regvars.b4c)
regvars.vp$Sex<-factor(abs(regvars.vp$sex -2))
regvars.vp$Age<-as.numeric(regvars.vp$age_at_visit)
regvars.vp$contributor_Fsplit<-factor(regvars.vp$contributor_Fsplit)
regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1)
regvars.vp$RegrBloodPreanalyticFactor.HBZ<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HBZ)
#regvars.vp$APOE.E4carrier.Proxy.LRRN1<-as.numeric(cleanDat.4p13b2["LRRN1|Q6UXK5^SL025922@seq.11293.14",])
regvars.vp$APOE.E4carrier<-factor(regvars.vp$APOE4.carrier)

# too many missing values:
#form <- ~ age_at_visit+(1|Sex)+(1|raceAA)+(1|recruited_control)+(1|ad)+(1|ftd)+(1|pd)+(1|als)+(1|mci_sci)+(1|mi)+(1|C9Orf72)+(1|GRN)+(1|MAPT)+APOE4.Dose+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+MMSE
#form <- ~ (1|contributor_Fsplit) +RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ
form <- ~ Age+(1|Sex)+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+(1|APOE.E4carrier)  #+APOE.E4carrier.Proxy.LRRN1

library(variancePartition)

# (If regressed) REMOVE regressed proteins to avoid "Response variable 4641 has a variance of 0":
#regrProts.idx<-c(which(grepl("^HNRNPA2B1\\|",rownames(cleanDat))),which(grepl("^HBZ\\|",rownames(cleanDat))))
#regrProts.idx
#5569 6087 -- take out
#previously: 6156 4641
#integer(0)  # if already removed

#previously: removed 2 regr proteins in line:  varPart.reg <- fitExtractVarPartModel(impute::impute.knn(cleanDat[-c(4641,6156),])$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = 8, type = "SOCK"))
varPart.b4c <- fitExtractVarPartModel(na.omit(as.matrix(cleanDat.b4c[,])), form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp.b4c <- sortCols(varPart.b4c,FUN=median,last= c("Residuals"))

pdf(file="4p13b4c.contributor_Fsplit_19sites1xPAVregr+SiteRegress_Protect_age+sex+carrierBinaryStatus-VariancePartition-PLASMA-7335x16677.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp.b4c, main="HDS 1.3ms - 4p13b4c - KNOWN APOE 19 sites 2x Regr(2PAV) + Site Regr, Prot age+sex+(e4 binary)" )

	SexSortOrder<-order(vp.b4c[["Sex"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b4c)) { vp.b4c[[i]]<-vp.b4c[[i]][SexSortOrder]; }
	rownames(vp.b4c)<-rownames(vp.b4c)[SexSortOrder]

	plotPercentBars( vp.b4c[1:50,]) + ggtitle( "Top Sex-covariates" )


	AgeSortOrder<-order(vp.b4c[["Age"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b4c)) { vp.b4c[[i]]<-vp.b4c[[i]][AgeSortOrder]; }
	rownames(vp.b4c)<-rownames(vp.b4c)[AgeSortOrder]

	plotPercentBars( vp.b4c[1:50,]) + ggtitle( "Top Age-covariates" )


        BatchSortOrder<-order(vp.b4c[["RegrBloodPreanalyticFactor.HNRNPA2B1"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4c)) { vp.b4c[[i]]<-vp.b4c[[i]][BatchSortOrder]; }
        rownames(vp.b4c)<-rownames(vp.b4c)[BatchSortOrder]

        plotPercentBars( vp.b4c[1:50,]) + ggtitle( "Top Time To Spin (HNRNPA2B1)-covariates" )


        BatchSortOrder<-order(vp.b4c[["RegrBloodPreanalyticFactor.HBZ"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4c)) { vp.b4c[[i]]<-vp.b4c[[i]][BatchSortOrder]; }
        rownames(vp.b4c)<-rownames(vp.b4c)[BatchSortOrder]

        plotPercentBars( vp.b4c[1:50,]) + ggtitle( "Top Preanalytical Factor 2 (HBZ)-covariates" )


        BatchSortOrder<-order(vp.b4c[["contributor_Fsplit"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4c)) { vp.b4c[[i]]<-vp.b4c[[i]][BatchSortOrder]; }
        rownames(vp.b4c)<-rownames(vp.b4c)[BatchSortOrder]

        plotPercentBars( vp.b4c[1:50,]) + ggtitle( "Top Contributor (site)-covariates" )


        BatchSortOrder<-order(vp.b4c[["APOE.E4carrier"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4c)) { vp.b4c[[i]]<-vp.b4c[[i]][BatchSortOrder]; }
        rownames(vp.b4c)<-rownames(vp.b4c)[BatchSortOrder]

        plotPercentBars( vp.b4c[1:50,]) + ggtitle( "Top APOE e4 carrier (Binary, no NA)-covariates" )


#	BatchSortOrder<-order(vp.b4c[["MMSE"]],decreasing=TRUE)
#	#rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
#	for (i in ls(vp.b4c)) { vp.b4c[[i]]<-vp.b4c[[i]][BatchSortOrder]; }
#	rownames(vp.b4c)<-rownames(vp.b4c)[BatchSortOrder]
#
#	plotPercentBars( vp.b4c[1:50,]) + ggtitle( "MMSE Cog. Score-covariates" )

dev.off()

varPart.b4c<-vp.b4c
saveRDS(varPart.b4c,"4p13b4c.contributor_Fsplit_19sites1xPAVregr+SiteRegress_Protect_age+sex+e4_carrierBinaryStatus.na.omit-varPart.b4c.RDS")

#############################################################################
## 4p13b4c. UNREGRESSED for site - Variance Partition (2PAV regression intrasite) (QC)

library(variancePartition)

varPart.b4c.unreg <- fitExtractVarPartModel(na.omit(as.matrix(cleanDat.unreg.b4c[,])), form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp.b4c.unreg <- sortCols(varPart.b4c.unreg,FUN=median,last= c("Residuals"))

pdf(file="4p13b4c.UNREGRESSED.contributor_Fsplit_19sites1xPAVregr-VariancePartition-PLASMA-7335x16677.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp.b4c.unreg, main="HDS 1.3ms - 4p13b4c UNREGRESSED - KNOWN APOE 19 sites 1x Regr(2PAV)" )

	SexSortOrder<-order(vp.b4c.unreg[["Sex"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b4c.unreg)) { vp.b4c.unreg[[i]]<-vp.b4c.unreg[[i]][SexSortOrder]; }
	rownames(vp.b4c.unreg)<-rownames(vp.b4c.unreg)[SexSortOrder]

	plotPercentBars( vp.b4c.unreg[1:50,]) + ggtitle( "Top Sex-covariates" )


	AgeSortOrder<-order(vp.b4c.unreg[["Age"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b4c.unreg)) { vp.b4c.unreg[[i]]<-vp.b4c.unreg[[i]][AgeSortOrder]; }
	rownames(vp.b4c.unreg)<-rownames(vp.b4c.unreg)[AgeSortOrder]

	plotPercentBars( vp.b4c.unreg[1:50,]) + ggtitle( "Top Age-covariates" )


        BatchSortOrder<-order(vp.b4c.unreg[["RegrBloodPreanalyticFactor.HNRNPA2B1"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4c.unreg)) { vp.b4c.unreg[[i]]<-vp.b4c.unreg[[i]][BatchSortOrder]; }
        rownames(vp.b4c.unreg)<-rownames(vp.b4c.unreg)[BatchSortOrder]

        plotPercentBars( vp.b4c.unreg[1:50,]) + ggtitle( "Top Time To Spin (HNRNPA2B1)-covariates" )


        BatchSortOrder<-order(vp.b4c.unreg[["RegrBloodPreanalyticFactor.HBZ"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4c.unreg)) { vp.b4c.unreg[[i]]<-vp.b4c.unreg[[i]][BatchSortOrder]; }
        rownames(vp.b4c.unreg)<-rownames(vp.b4c.unreg)[BatchSortOrder]

        plotPercentBars( vp.b4c.unreg[1:50,]) + ggtitle( "Top Preanalytical Factor 2 (HBZ)-covariates" )


        BatchSortOrder<-order(vp.b4c.unreg[["contributor_Fsplit"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4c.unreg)) { vp.b4c.unreg[[i]]<-vp.b4c.unreg[[i]][BatchSortOrder]; }
        rownames(vp.b4c.unreg)<-rownames(vp.b4c.unreg)[BatchSortOrder]

        plotPercentBars( vp.b4c.unreg[1:50,]) + ggtitle( "Top Contributor (site)-covariates" )


        BatchSortOrder<-order(vp.b4c.unreg[["APOE.E4carrier"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4c.unreg)) { vp.b4c.unreg[[i]]<-vp.b4c.unreg[[i]][BatchSortOrder]; }
        rownames(vp.b4c.unreg)<-rownames(vp.b4c.unreg)[BatchSortOrder]

        plotPercentBars( vp.b4c.unreg[1:50,]) + ggtitle( "Top APOE e4 carrier (Binary, no NA)-covariates" )


#	BatchSortOrder<-order(vp.b4c.unreg[["MMSE"]],decreasing=TRUE)
#	#rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
#	for (i in ls(vp.b4c.unreg)) { vp.b4c.unreg[[i]]<-vp.b4c.unreg[[i]][BatchSortOrder]; }
#	rownames(vp.b4c.unreg)<-rownames(vp.b4c.unreg)[BatchSortOrder]
#
#	plotPercentBars( vp.b4c.unreg[1:50,]) + ggtitle( "MMSE Cog. Score-covariates" )

dev.off()

varPart.b4c.unreg<-vp.b4c.unreg
saveRDS(varPart.b4c.unreg,"4p13b4c.UNREGRESSED.contributor_Fsplit_19sites1xPAVregr.na.omit-varPart.b4c.unreg.RDS")


cleanDat.b4c.NOsiteH<-t(na.omit(t(cleanDat.b4c)))
regvars.b4c.NOsiteH<-regvars.b4c[colnames(cleanDat.b4c.NOsiteH),]

## Better histogram using top 6 correlates (3 + and 3 -)
SPC.med=median(cleanDat.b4c.NOsiteH["SPC25|Q9HBM1",],na.rm=T)
LRRN1.med=median(cleanDat.b4c.NOsiteH["LRRN1|Q6UXK5^SL025922@seq.11293.14",],na.rm=T)
CTF1.med=median(cleanDat.b4c.NOsiteH["CTF1|Q16619^SL002783@seq.13732.79",],na.rm=T)
NEFL.med=median(cleanDat.b4c.NOsiteH["NEFL|P07196",],na.rm=T)
TBCA.med=median(cleanDat.b4c.NOsiteH["TBCA|O75347",],na.rm=T)
S100A13.med=median(cleanDat.b4c.NOsiteH["S100A13|Q99584",],na.rm=T)
hist.data.SPC=((cleanDat.b4c.NOsiteH["SPC25|Q9HBM1",] - SPC.med) + (cleanDat.b4c.NOsiteH["LRRN1|Q6UXK5^SL025922@seq.11293.14",] - LRRN1.med) + (cleanDat.b4c.NOsiteH["CTF1|Q16619^SL002783@seq.13732.79",] - CTF1.med) - (cleanDat.b4c.NOsiteH["NEFL|P07196",] - NEFL.med) - (cleanDat.b4c.NOsiteH["S100A13|Q99584",] - S100A13.med) - (cleanDat.b4c.NOsiteH["TBCA|O75347",] - TBCA.med))
hist(hist.data.SPC, breaks=100, xlab="log2(abundance ratio):  SPC+LRRN1+CTF1(median-centered) / TBCA+S100A13+NEFL(median-centered)", main="APOE e4 Carrier Status Best Nonmissing Correlate")
hist(hist.data.SPC[which(regvars.b4c.NOsiteH$APOE==24 | regvars.b4c.NOsiteH$APOE==34)],breaks=60,col="#BBFFBB40",add=T)  #green overlay
hist(hist.data.SPC[which(regvars.b4c.NOsiteH$APOE==44)],breaks=70,col="#FFBBBB99",add=T)  #red overlay
hist(hist.data.SPC[which(regvars.b4c.NOsiteH$APOE==33)],breaks=100,col="#BBBBFF99",add=T)  #blue
legend("topright",c("e4/e4","e3/e3", "e4 het"),fill=c("#FFBBBB99","#BBBBFF99","#BBFFBB40"))
abline(v=0.625,col="darkgreen",lty=2, lwd=2.2)
abline(v=4.25,col="maroon",lty=2, lwd=2.2)
# saved plot to powerpoint (capture)


######################################
## ANOVA + Volcanoes + DEXstacked Barplots
source("parANOVA.dex.fallback7.25.R")

parallelThreads=20
outFilePrefix="4p13b4c.regr.knownAPOEonly2xRegr"
outFileSuffix="SixDxVolcs.CTimputed.3sitesMappedDx"

cleanDat<-cleanDat.b4c

#Grouping=numericMeta$Group #***
Grouping=regvars.b4c$Group.withCTimputed
Grouping[which(regvars.b4c$Group.withCTimputed=="CI.Other")]<-NA
Grouping[which(regvars.b4c$Group.withCTimputed=="AsymAD")]<-NA
#Grouping[which(regvars.b4c$Group.withCTimputed=="MCI")]<-NA

ANOVAout <- parANOVA.dex()

#*** commented new for 3-genotype only comparisons
##head(ANOVAout)
## plot volcanoes only for comparisons with 2/3 Group components in common.
#comps<-cbind(ANOVAout.colnames=colnames(ANOVAout), as.data.frame(do.call(rbind,strsplit(colnames(ANOVAout),"[.|-]"))))
#lastComp=which(grepl("diff ",comps$ANOVAout.colnames))[1]-1
#lastComp
#155  # 155-2 = 153 comparisons pairwise
## keep comparisons with 2/3 variables the same (n=45 comparisons)
#selectComps=unlist( sapply(3:lastComp,function(x) if(length(which(duplicated(unlist(comps[x,2:7]))))>1) x ) )
#length(selectComps) # selectComps is used directly by plotVolc to choose which comparisons to plot (column #s of p values in ANOVAout)
## Decide which comparison X axis (logFC) to flip
#comps[selectComps,1] # view pairwise comparisons selected for volcano plotting
# None need flipping (more extreme expected phenotype is first (numerator) in each case)
flip=c(3,4,5,6,10,11)
sameScale=FALSE #TRUE
symbolsOnly=TRUE
highlightGeneProducts=c("HNRNPA2B1","HBZ","SPC25","CPLX2","PTN","MAPT","OMG","GDF15","NEFL","NRGN","CTHRC1","HTRA1","TTR","LRP1","NTN1","SFRP1","GPNMB","PAFAH1B3")
labelHighlighted=TRUE      # if true, highlighted spots get text labels with their rownames from ANOVAout
labelTop=5
plotVolc()                 # runs on ANOVAout as input (need not be specified).

#DEXpercentStacked()        # runs on prior function outputs as input; writes stacked bar plot(s) to PDF.


Grouping=as.integer(regvars.b4c$APOE4.carrier)-1
Grouping[Grouping==0]<-"non.carrier"
Grouping[Grouping==1]<-"e4.carrier"
outFileSuffix="e4carrier.vs.NonCarrier"

ANOVAout <- parANOVA.dex()  # current .csv output has labels flipped (forgot to subtract 1 from factor levels above)

#flip=c(3)
flip=c()
colnames(ANOVAout)[c(3:4)]<-c("e4.carrier-non.carrier","diff e4.carrier-non.carrier")
#16 p values are 0! (and FDR too)
# Set to 1e-175 (smaller than 1e-267, first nonzero sorted value)
ANOVAout[which(ANOVAout[,3]==0),3]<- 1e-275
sameScale=FALSE #TRUE
symbolsOnly=TRUE
highlightGeneProducts=c("HNRNPA2B1","HBZ","SPC25","CPLX2","PTN","MAPT","OMG","GDF15","NEFL","NRGN","CTHRC1","HTRA1","TTR","LRP1","NTN1","SFRP1","GPNMB","PAFAH1B3")
labelHighlighted=TRUE      # if true, highlighted spots get text labels with their rownames from ANOVAout
labelTop=20
plotVolc()                 # runs on ANOVAout as input (need not be specified).


####### 4p13b4d
## Output e3/e3 and e4/e4 ground truth samples with regression of site, protecting ONLY age and sex in output data from within-site 2PAV protein regression

library("doParallel")
parallelThreads=31  #now Windows02  #max is number of processes that can run on your computer at one time
stopCluster(clusterLocal)
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="PSOCK")

registerDoParallel(clusterLocal)


####################### SITE REGRESSION - pass b4d - 16,677 samples, protect only Age+Sex
  which(!names(gt.APOE)==colnames(cleanDat.unreg.b345))

#  cleanDat.unreg.b4c <- cleanDat.unreg.b345[,which(!is.na(gt.APOE))]
#  gt.APOE.b4c<-gt.APOE[colnames(cleanDat.unreg.b4c)]
#
#  regvars.b4c<-numericMeta.reg.b345[which(!is.na(gt.APOE)),]
#  regvars.b4c$Age<-as.numeric(numericMeta.reg.b345$age_at_visit[match(names(gt.APOE.b4c), rownames(numericMeta.reg.b345))])
#  regvars.b4c$Sex.int=as.integer(abs(numericMeta.reg.b345$sex[match(names(gt.APOE.b4c), rownames(numericMeta.reg.b345))] -2))
#  regvars.b4c$Sex<-relevel(factor(regvars.b4c$Sex.int), ref="0")
#  regvars.b4c$APOE4.carrier<-relevel(factor(ifelse(gt.APOE.b4c %in% c("e2/e4","e3/e4","e4/e4"), 1, 0)), ref="0")
#  regvars.b4c$APOE4.int=as.integer(regvars.b4c$APOE4.carrier)
#  regvars.b4c$Site <-regvars.b4c$contributor_Fsplit
##  regvars.b4c<-na.omit(regvars.b4c)  # removes samples with missing age, sex
  dim(regvars.b4c)
  # 16677   105

  ##  Run the regression (4p13 b4d) - Site with Age+Sex protected (only ground truth APOE genotyped samples included)

  ## covariate names you want to protect
  cov_keep <- c("Age", "Sex")   # Technically, we cleared out NAs, but one could list protection covariates that have missing values here

  normExpr.reg <- matrix(NA,nrow=nrow(cleanDat.unreg.b4c),ncol=ncol(cleanDat.unreg.b4c), dimnames = dimnames(cleanDat.unreg.b4c))
  good_samp=which(complete.cases(regvars.b4c[, cov_keep]))
  normExpr.reg <-  foreach (i=1:nrow(cleanDat.unreg.b4c), .combine=rbind, .packages="stats") %dopar% {  #** .combine=c, .multicombine=TRUE, .packages="stats", .export="regvars") %dopar% {
    y <- as.numeric(cleanDat.unreg.b4c[i, ])
    ## keep = samples where y and all covariates are present
    keep <- intersect(which(!is.na(y)), good_samp)             # columns that have a value, and regress with variables that have a value (we know Age and Sex are missing some)

    ## initialise result for this row as all NA
    adj <- rep(NA_real_, length(y))

    fit <- tryCatch(
        lm(y[keep]~Age+Sex +Site, data = regvars.b4c[keep, , drop = FALSE]),
        error = function(e) NULL)

    if (!is.null(fit)) {
        coef <- coef(fit)
        ## coefficient[1] + residuals -> adjusted expression
        adj[keep] <- coef["(Intercept)"] + coef["Age"]*regvars.b4c[keep,"Age"] + coef[which(grepl("^Sex", names(coef)))]*regvars.b4c[keep,"Sex.int"] + residuals(fit)
    }

#        if (!is.null(fit)) {
#            ## build a model matrix with Site *zeroed out*
#            X     <- model.matrix(fit)
#            X[ , grep("^Site", colnames(X)) ] <- 0     # drop site contribution
#            y_hat <- as.numeric(X %*% coef(fit))       # fitted without Site
#
#            ## protected expression = y_hat + residuals  ( == y - Site effect )
#            adj[keep] <- y_hat + residuals(fit)
#        }
    adj                                    # returned to foreach
  }
  dimnames(normExpr.reg) <- dimnames(cleanDat.unreg.b4c)

  cat(paste0("Finished Pass 4p13 b4d regression of intersite variance.\n"))

# sanity check: some sites have missing values:
table(apply(normExpr.reg,2,function(x) length(which(is.na(x)))))
#    0  2303  <- # of missing values in column
#14758  2254  <- # of columns with that many NA

saveRDS(normExpr.reg,"4p13b4d.normExpr.reg_sites1-19_Fsplit_ProtAge+Sex-knownAPOE_only.RDS")
cleanDat.b4d<-normExpr.reg
#regvars.b4c still valid, ordered for this cleanDat
#regvars.b4c$APOE.mapped<-gt.APOE.b4c
#same: saveRDS(regvars.b4c, "4p13b4c.numericMeta_sites1-19_traits_knownAPOE_only.RDS")

## Sanity check
e4.bicor.to.siteCorr.b4d.assays<-bicor(t(cleanDat.b4d),regvars.b4c$APOE4.carrier, use='p')
# sort by bicor
e4.bicor.to.siteCorr.b4d.assays[order(unlist(t(e4.bicor.to.siteCorr.b4d.assays)),decreasing=TRUE),][c(1:5,7331:7335)]
#b4d (current Age+Sex protected only)
#SPC25|Q9HBM1   LRRN1|Q6UXK5^SL025922@seq.11293.14   CTF1|Q16619^SL002783@seq.13732.79            NEFL|P07196        S100A13|Q99584                TBCA|O75347
#   0.7636772                            0.6565787                           0.4529732      ...    -0.5720141            -0.5735374                 -0.6058351

#b4c (Age+Sex+Carrier binary status protected, prior)
#SPC25|Q9HBM1   LRRN1|Q6UXK5^SL025922@seq.11293.14   CTF1|Q16619^SL002783@seq.13732.79            NEFL|P07196        S100A13|Q99584                TBCA|O75347
#   0.7702713                            0.6725630                           0.4673459      ...    -0.5905291            -0.5929931                 -0.6264711


##################################################
# 4p13b4d. 19 sites (no site U,V,W; site F split to F1, F2, F3) 1xPAVregr+sitewise Regress, Protect Age+Sex QC (tSNE plots)

####################
## Examine tSNE of PLASMA - Human only assays (rows) and samples (columns);  no missing data by row.

numericMeta.plasma<-regvars.b4c
dim(numericMeta.plasma)
#16677  105
exprMat.plasma<-cleanDat.b4d

dim(as.data.frame(na.omit(exprMat.plasma)))
#  5032 16677   # 5031 in H; 7334 in K; 7335 in all others.
#exprMat.plasma<-exprMat.plasma[,which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))]
#numericMeta.plasma<-numericMeta.plasma[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma))))),]

#Group.3mappedCohortsPlusCTimputed alternate
Group=as.factor(regvars.b4c$Group.withCTimputed ) #[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))])


tSNE.list.plasma <- Rtsne::Rtsne(t(as.data.frame(na.omit(exprMat.plasma))),perplexity=20)
dim(tSNE.list.plasma$Y)
#16677   2
tSNE.plasma.xy<-as.data.frame(tSNE.list.plasma$Y)
colnames(tSNE.plasma.xy)<-c('x','y')
tSNE.4p13b4d.plasma.xy<-tSNE.plasma.xy
rownames(tSNE.4p13b4d.plasma.xy)<-rownames(regvars.b4c) #[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))]
library(ggplot2)
#library(ggpubr) - rlang upgrade required, ggplot2 upgrade required. cannot install from older source!
library(ggrepel)

# Get the indices of the first occurrence of each unique value in contributor_Fsplit
first_occurrence_indices.4p13b4d<-first_occurrence_indices <- match(unique(regvars.b4c$contributor_Fsplit), regvars.b4c$contributor_Fsplit)

labels.4p13b4d<-labels<-regvars.b4c$contributor_Fsplit[first_occurrence_indices.4p13b4d]

tSNE.plasma.samples.sites<-ggplot2::ggplot(tSNE.4p13b4d.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=regvars.b4c$contributor_Fsplit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13b4d.plasma.xy[first_occurrence_indices.4p13b4d, ],
                  aes(x=x,y=y, label = labels.4p13b4d), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.age<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$age_at_visit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleMatrix<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_matrix), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleType<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_type), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.apoe4dose<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$APOE4.Dose), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.MMSE<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$MMSE), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.Group<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=Group), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )


pdf(file="4p13b4d.2PAVproteinIntrasiteRegressed+Regress19sites_protectAgeSex_Fsplit.tSNE-Plasma(7335x16677)-samples_coloredByTraits.pdf",width=11,height=9)
  print(tSNE.plasma.samples.sites + labs(colour="Site"))
  print(tSNE.plasma.samples.age)
  print(tSNE.plasma.samples.sampleMatrix)
  print(tSNE.plasma.samples.sampleType)
  print(tSNE.plasma.samples.apoe4dose)
  print(tSNE.plasma.samples.MMSE)
  print(tSNE.plasma.samples.Group + labs(colour="Group"))
dev.off()


tSNE.4p13b4d.plasma.samples.sites<-tSNE.plasma.samples.sites


library("doParallel")
parallelThreads=31  #now Windows02  #max is number of processes that can run on your computer at one time
stopCluster(clusterLocal)
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="PSOCK")

registerDoParallel(clusterLocal)

#############################################################################
## 4p13b4d. Variance Partition regressed (2PAV regression intrasite)+Site regressed; protect Age+Sex (no NA) (QC)

#regvars.vp<-data.frame(regvars.b4c)
#regvars.vp$Sex<-factor(abs(regvars.vp$sex -2))
#regvars.vp$Age<-as.numeric(regvars.vp$age_at_visit)
#regvars.vp$contributor_Fsplit<-factor(regvars.vp$contributor_Fsplit)
#regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1)
#regvars.vp$RegrBloodPreanalyticFactor.HBZ<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HBZ)
##regvars.vp$APOE.E4carrier.Proxy.LRRN1<-as.numeric(cleanDat.4p13b2["LRRN1|Q6UXK5^SL025922@seq.11293.14",])
#regvars.vp$APOE.E4carrier<-factor(regvars.vp$APOE4.carrier)

# too many missing values:
form <- ~ Age+(1|Sex)+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+(1|APOE.E4carrier)  #+APOE.E4carrier.Proxy.LRRN1

library(variancePartition)

# (If regressed) REMOVE regressed proteins to avoid "Response variable 4641 has a variance of 0":
#regrProts.idx<-c(which(grepl("^HNRNPA2B1\\|",rownames(cleanDat))),which(grepl("^HBZ\\|",rownames(cleanDat))))
#regrProts.idx
#5569 6087 -- take out
#previously: 6156 4641
#integer(0)  # if already removed

#previously: removed 2 regr proteins in line:  varPart.reg <- fitExtractVarPartModel(impute::impute.knn(cleanDat[-c(4641,6156),])$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = 8, type = "SOCK"))
varPart.b4d <- fitExtractVarPartModel(na.omit(as.matrix(cleanDat.b4d[,])), form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp.b4d <- sortCols(varPart.b4d,FUN=median,last= c("Residuals"))

pdf(file="4p13b4d.contributor_Fsplit_19sites1xPAVregr+SiteRegress_Protect_age+sex-VariancePartition-PLASMA-7335x16677.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp.b4d, main="HDS 1.3ms - 4p13b4d - KNOWN APOE 19 sites 2x Regr(2PAV) + Site Regr, Prot age+sex" )

	SexSortOrder<-order(vp.b4d[["Sex"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b4d)) { vp.b4d[[i]]<-vp.b4d[[i]][SexSortOrder]; }
	rownames(vp.b4d)<-rownames(vp.b4d)[SexSortOrder]

	plotPercentBars( vp.b4d[1:50,]) + ggtitle( "Top Sex-covariates" )


	AgeSortOrder<-order(vp.b4d[["Age"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b4d)) { vp.b4d[[i]]<-vp.b4d[[i]][AgeSortOrder]; }
	rownames(vp.b4d)<-rownames(vp.b4d)[AgeSortOrder]

	plotPercentBars( vp.b4d[1:50,]) + ggtitle( "Top Age-covariates" )


        BatchSortOrder<-order(vp.b4d[["RegrBloodPreanalyticFactor.HNRNPA2B1"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4d)) { vp.b4d[[i]]<-vp.b4d[[i]][BatchSortOrder]; }
        rownames(vp.b4d)<-rownames(vp.b4d)[BatchSortOrder]

        plotPercentBars( vp.b4d[1:50,]) + ggtitle( "Top Time To Spin (HNRNPA2B1)-covariates" )


        BatchSortOrder<-order(vp.b4d[["RegrBloodPreanalyticFactor.HBZ"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4d)) { vp.b4d[[i]]<-vp.b4d[[i]][BatchSortOrder]; }
        rownames(vp.b4d)<-rownames(vp.b4d)[BatchSortOrder]

        plotPercentBars( vp.b4d[1:50,]) + ggtitle( "Top Preanalytical Factor 2 (HBZ)-covariates" )


        BatchSortOrder<-order(vp.b4d[["contributor_Fsplit"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4d)) { vp.b4d[[i]]<-vp.b4d[[i]][BatchSortOrder]; }
        rownames(vp.b4d)<-rownames(vp.b4d)[BatchSortOrder]

        plotPercentBars( vp.b4d[1:50,]) + ggtitle( "Top Contributor (site)-covariates" )


        BatchSortOrder<-order(vp.b4d[["APOE.E4carrier"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4d)) { vp.b4d[[i]]<-vp.b4d[[i]][BatchSortOrder]; }
        rownames(vp.b4d)<-rownames(vp.b4d)[BatchSortOrder]

        plotPercentBars( vp.b4d[1:50,]) + ggtitle( "Top APOE e4 carrier (Binary, no NA)-covariates" )


#	BatchSortOrder<-order(vp.b4d[["MMSE"]],decreasing=TRUE)
#	#rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
#	for (i in ls(vp.b4d)) { vp.b4d[[i]]<-vp.b4d[[i]][BatchSortOrder]; }
#	rownames(vp.b4d)<-rownames(vp.b4d)[BatchSortOrder]
#
#	plotPercentBars( vp.b4d[1:50,]) + ggtitle( "MMSE Cog. Score-covariates" )

dev.off()

varPart.b4d<-vp.b4d
saveRDS(varPart.b4d,"4p13b4d.contributor_Fsplit_19sites1xPAVregr+SiteRegress_Protect_age+sex.na.omit-varPart.b4d.RDS")


######################################
## ANOVA + Volcanoes + DEXstacked Barplots (b4d)
source("parANOVA.dex.fallback7.25.R")

parallelThreads=20
outFilePrefix="4p13b4d.regr.knownAPOEonly2xRegr"
outFileSuffix="SixDxVolcs.CTimputed.3sitesMappedDx"

cleanDat<-cleanDat.b4d

#Grouping=numericMeta$Group #***
Grouping=regvars.b4c$Group.withCTimputed
Grouping[which(regvars.b4c$Group.withCTimputed=="CI.Other")]<-NA
Grouping[which(regvars.b4c$Group.withCTimputed=="AsymAD")]<-NA
#Grouping[which(regvars.b4c$Group.withCTimputed=="MCI")]<-NA

ANOVAout <- parANOVA.dex()

#*** commented new for 3-genotype only comparisons
##head(ANOVAout)
## plot volcanoes only for comparisons with 2/3 Group components in common.
#comps<-cbind(ANOVAout.colnames=colnames(ANOVAout), as.data.frame(do.call(rbind,strsplit(colnames(ANOVAout),"[.|-]"))))
#lastComp=which(grepl("diff ",comps$ANOVAout.colnames))[1]-1
#lastComp
#155  # 155-2 = 153 comparisons pairwise
## keep comparisons with 2/3 variables the same (n=45 comparisons)
#selectComps=unlist( sapply(3:lastComp,function(x) if(length(which(duplicated(unlist(comps[x,2:7]))))>1) x ) )
#length(selectComps) # selectComps is used directly by plotVolc to choose which comparisons to plot (column #s of p values in ANOVAout)
## Decide which comparison X axis (logFC) to flip
#comps[selectComps,1] # view pairwise comparisons selected for volcano plotting
# None need flipping (more extreme expected phenotype is first (numerator) in each case)
flip=c(3,4,5,6,10,11)
sameScale=FALSE #TRUE
symbolsOnly=TRUE
highlightGeneProducts=c("HNRNPA2B1","HBZ","SPC25","CPLX2","PTN","MAPT","OMG","GDF15","NEFL","NRGN","CTHRC1","HTRA1","TTR","LRP1","NTN1","SFRP1","GPNMB","PAFAH1B3")
labelHighlighted=TRUE      # if true, highlighted spots get text labels with their rownames from ANOVAout
labelTop=5
plotVolc()                 # runs on ANOVAout as input (need not be specified).

#DEXpercentStacked()        # runs on prior function outputs as input; writes stacked bar plot(s) to PDF.


Grouping=as.integer(regvars.b4c$APOE4.carrier)-1
Grouping[Grouping==0]<-"non.carrier"
Grouping[Grouping==1]<-"e4.carrier"
outFileSuffix="e4carrier.vs.NonCarrier"

ANOVAout <- parANOVA.dex()  # current .csv output has labels flipped (forgot to subtract 1 from factor levels above)

flip=c(3)
#flip=c()
#colnames(ANOVAout)[c(3:4)]<-c("e4.carrier-non.carrier","diff e4.carrier-non.carrier")
#16 p values are 0! (and FDR too)
# Set to 1e-175 (smaller than 1e-267, first nonzero sorted value)
ANOVAout[which(ANOVAout[,3]==0),3]<- 1e-275
sameScale=FALSE #TRUE
symbolsOnly=TRUE
highlightGeneProducts=c("HNRNPA2B1","HBZ","SPC25","CPLX2","PTN","MAPT","OMG","GDF15","NEFL","NRGN","CTHRC1","HTRA1","TTR","LRP1","NTN1","SFRP1","GPNMB","PAFAH1B3")
labelHighlighted=TRUE      # if true, highlighted spots get text labels with their rownames from ANOVAout
labelTop=20
plotVolc()                 # runs on ANOVAout as input (need not be specified).


####### 4p13b4e
## Output e3/e3 and e4/e4 ground truth samples with regression of site, protecting ONLY age and sex in output data from within-site 2PAV protein regression

library("doParallel")
parallelThreads=8  #now Windows02  #max is number of processes that can run on your computer at one time
stopCluster(clusterLocal)
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="PSOCK")

registerDoParallel(clusterLocal)


####################### SITE REGRESSION - pass b4e - 16,677 samples, protect only Age+Sex
  which(!names(gt.APOE)==colnames(cleanDat.unreg.b345))

#  cleanDat.unreg.b4c <- cleanDat.unreg.b345[,which(!is.na(gt.APOE))]
#  gt.APOE.b4c<-gt.APOE[colnames(cleanDat.unreg.b4c)]
#
#  regvars.b4c<-numericMeta.reg.b345[which(!is.na(gt.APOE)),]
#  regvars.b4c$Age<-as.numeric(numericMeta.reg.b345$age_at_visit[match(names(gt.APOE.b4c), rownames(numericMeta.reg.b345))])
#  regvars.b4c$Sex.int=as.integer(abs(numericMeta.reg.b345$sex[match(names(gt.APOE.b4c), rownames(numericMeta.reg.b345))] -2))
#  regvars.b4c$Sex<-relevel(factor(regvars.b4c$Sex.int), ref="0")
#  regvars.b4c$APOE4.carrier<-relevel(factor(ifelse(gt.APOE.b4c %in% c("e2/e4","e3/e4","e4/e4"), 1, 0)), ref="0")
#  regvars.b4c$APOE4.int=as.integer(regvars.b4c$APOE4.carrier)
#  regvars.b4c$Site <-regvars.b4c$contributor_Fsplit
##  regvars.b4c<-na.omit(regvars.b4c)  # removes samples with missing age, sex
  dim(regvars.b4c)
  # 16677   105

  ##  Run the regression (4p13 b4e) - Site with Age+Sex protected (only ground truth APOE genotyped samples included)

  ## covariate names you want to protect
  cov_keep <- c("Age", "Sex")   # Technically, we cleared out NAs, but one could list protection covariates that have missing values here

  normExpr.reg <- matrix(NA,nrow=nrow(cleanDat.unreg.b4c),ncol=ncol(cleanDat.unreg.b4c), dimnames = dimnames(cleanDat.unreg.b4c))
  good_samp=which(complete.cases(regvars.b4c[, cov_keep]))
  normExpr.reg <-  foreach (i=1:nrow(cleanDat.unreg.b4c), .combine=rbind, .packages="stats") %dopar% {  #** .combine=c, .multicombine=TRUE, .packages="stats", .export="regvars") %dopar% {
    y <- as.numeric(cleanDat.unreg.b4c[i, ])
    ## keep = samples where y and all covariates are present
    keep <- intersect(which(!is.na(y)), good_samp)             # columns that have a value, and regress with variables that have a value (we know Age and Sex are missing some)

# ------------------------------------------------------------------------
# ANNOTATION: Run differential-expression/volcano QC contrasts on candidate
# harmonized matrices.
# ------------------------------------------------------------------------

    ## initialise result for this row as all NA
    adj <- rep(NA_real_, length(y))

    fit <- tryCatch(
        lm(y[keep]~Site, data = regvars.b4c[keep, , drop = FALSE]),
        error = function(e) NULL)

    if (!is.null(fit)) {
        coef <- coef(fit)
        ## coefficient[1] + residuals -> adjusted expression
        adj[keep] <- coef["(Intercept)"] + residuals(fit)
    }

#        if (!is.null(fit)) {
#            ## build a model matrix with Site *zeroed out*
#            X     <- model.matrix(fit)
#            X[ , grep("^Site", colnames(X)) ] <- 0     # drop site contribution
#            y_hat <- as.numeric(X %*% coef(fit))       # fitted without Site
#
#            ## protected expression = y_hat + residuals  ( == y - Site effect )
#            adj[keep] <- y_hat + residuals(fit)
#        }
    adj                                    # returned to foreach
  }
  dimnames(normExpr.reg) <- dimnames(cleanDat.unreg.b4c)

  cat(paste0("Finished Pass 4p13 b4e regression of intersite variance.\n"))

# sanity check: some sites have missing values:
table(apply(normExpr.reg,2,function(x) length(which(is.na(x)))))
#    0  2303  <- # of missing values in column
#14758  2254  <- # of columns with that many NA

saveRDS(normExpr.reg,"4p13b4e.normExpr.reg_sites1-19_Fsplit_ProtNothing-knownAPOE_only.RDS")
cleanDat.b4e<-normExpr.reg
#regvars.b4c still valid, ordered for this cleanDat
#regvars.b4c$APOE.mapped<-gt.APOE.b4c
#same: saveRDS(regvars.b4c, "4p13b4c.numericMeta_sites1-19_traits_knownAPOE_only.RDS")

## Sanity check
e4.bicor.to.siteCorr.b4e.assays<-bicor(t(cleanDat.b4e),regvars.b4c$APOE4.carrier, use='p')
# sort by bicor
e4.bicor.to.siteCorr.b4e.assays[order(unlist(t(e4.bicor.to.siteCorr.b4e.assays)),decreasing=TRUE),][c(1:5,7331:7335)]
#b4e (nothing protected, current)
#SPC25|Q9HBM1   LRRN1|Q6UXK5^SL025922@seq.11293.14   CTF1|Q16619^SL002783@seq.13732.79            NEFL|P07196        S100A13|Q99584                TBCA|O75347
#   0.7612027                            0.6574740                           0.4543276      ...    -0.5739026            -0.5788553                 -0.6112647

#b4d (Age+Sex protected only, prior)
#SPC25|Q9HBM1   LRRN1|Q6UXK5^SL025922@seq.11293.14   CTF1|Q16619^SL002783@seq.13732.79            NEFL|P07196        S100A13|Q99584                TBCA|O75347
#   0.7636772                            0.6565787                           0.4529732      ...    -0.5720141            -0.5735374                 -0.6058351

#b4c (Age+Sex+Carrier binary status protected, prior)
#SPC25|Q9HBM1   LRRN1|Q6UXK5^SL025922@seq.11293.14   CTF1|Q16619^SL002783@seq.13732.79            NEFL|P07196        S100A13|Q99584                TBCA|O75347
#   0.7702713                            0.6725630                           0.4673459      ...    -0.5905291            -0.5929931                 -0.6264711


##################################################
# 4p13b4e. 19 sites (no site U,V,W; site F split to F1, F2, F3) 1xPAVregr+sitewise Regress, Protect Age+Sex+APOE.e4 carrier status (binary, imputed NA) QC (tSNE plots)

####################
## Examine tSNE of PLASMA - Human only assays (rows) and samples (columns);  no missing data by row.

numericMeta.plasma<-regvars.b4c
dim(numericMeta.plasma)
#16677  105
exprMat.plasma<-cleanDat.b4e

dim(as.data.frame(na.omit(exprMat.plasma)))
#  5032 16677   # 5031 in H; 7334 in K; 7335 in all others.
#exprMat.plasma<-exprMat.plasma[,which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))]
#numericMeta.plasma<-numericMeta.plasma[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma))))),]

#Group.3mappedCohortsPlusCTimputed alternate
Group=as.factor(regvars.b4c$Group.withCTimputed ) #[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))])


tSNE.list.plasma <- Rtsne::Rtsne(t(as.data.frame(na.omit(exprMat.plasma))),perplexity=20)
dim(tSNE.list.plasma$Y)
#16677   2
tSNE.plasma.xy<-as.data.frame(tSNE.list.plasma$Y)
colnames(tSNE.plasma.xy)<-c('x','y')
tSNE.4p13b4e.plasma.xy<-tSNE.plasma.xy
rownames(tSNE.4p13b4e.plasma.xy)<-rownames(regvars.b4c) #[which(!duplicated(t(as.data.frame(na.omit(exprMat.plasma)))))]
library(ggplot2)
#library(ggpubr) - rlang upgrade required, ggplot2 upgrade required. cannot install from older source!
library(ggrepel)

# Get the indices of the first occurrence of each unique value in contributor_Fsplit
first_occurrence_indices.4p13b4e<-first_occurrence_indices <- match(unique(regvars.b4c$contributor_Fsplit), regvars.b4c$contributor_Fsplit)

labels.4p13b4e<-labels<-regvars.b4c$contributor_Fsplit[first_occurrence_indices.4p13b4e]

tSNE.plasma.samples.sites<-ggplot2::ggplot(tSNE.4p13b4e.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=regvars.b4c$contributor_Fsplit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13b4e.plasma.xy[first_occurrence_indices.4p13b4e, ],
                  aes(x=x,y=y, label = labels.4p13b4e), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.age<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$age_at_visit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleMatrix<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_matrix), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.sampleType<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_type), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.apoe4dose<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$APOE4.Dose), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.MMSE<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$MMSE), size=0.35) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.Group<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=Group), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.plasma.xy[first_occurrence_indices, ],
                  aes(x=x,y=y, label = labels), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )


pdf(file="4p13b4e.2PAVproteinIntrasiteRegressed+Regress19sites_protectNothing_Fsplit.tSNE-Plasma(7335x16677)-samples_coloredByTraits.pdf",width=11,height=9)
  print(tSNE.plasma.samples.sites + labs(colour="Site"))
  print(tSNE.plasma.samples.age)
  print(tSNE.plasma.samples.sampleMatrix)
  print(tSNE.plasma.samples.sampleType)
  print(tSNE.plasma.samples.apoe4dose)
  print(tSNE.plasma.samples.MMSE)
  print(tSNE.plasma.samples.Group + labs(colour="Group"))
dev.off()


tSNE.4p13b4e.plasma.samples.sites<-tSNE.plasma.samples.sites


library("doParallel")
parallelThreads=8  #now Windows02  #max is number of processes that can run on your computer at one time
stopCluster(clusterLocal)
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="PSOCK")

registerDoParallel(clusterLocal)

#############################################################################
## 4p13b4e. Variance Partition regressed (2PAV regression intrasite)+Site regressed; protect Age+Sex (no NA) (QC)

#regvars.vp<-data.frame(regvars.b4c)
#regvars.vp$Sex<-factor(abs(regvars.vp$sex -2))
#regvars.vp$Age<-as.numeric(regvars.vp$age_at_visit)
#regvars.vp$contributor_Fsplit<-factor(regvars.vp$contributor_Fsplit)
#regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1)
#regvars.vp$RegrBloodPreanalyticFactor.HBZ<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HBZ)
##regvars.vp$APOE.E4carrier.Proxy.LRRN1<-as.numeric(cleanDat.4p13b2["LRRN1|Q6UXK5^SL025922@seq.11293.14",])
#regvars.vp$APOE.E4carrier<-factor(regvars.vp$APOE4.carrier)

# too many missing values:
form <- ~ Age+(1|Sex)+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+(1|APOE.E4carrier)  #+APOE.E4carrier.Proxy.LRRN1

library(variancePartition)

# (If regressed) REMOVE regressed proteins to avoid "Response variable 4641 has a variance of 0":
#regrProts.idx<-c(which(grepl("^HNRNPA2B1\\|",rownames(cleanDat))),which(grepl("^HBZ\\|",rownames(cleanDat))))
#regrProts.idx
#5569 6087 -- take out
#previously: 6156 4641
#integer(0)  # if already removed

#previously: removed 2 regr proteins in line:  varPart.reg <- fitExtractVarPartModel(impute::impute.knn(cleanDat[-c(4641,6156),])$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = 8, type = "SOCK"))
varPart.b4e <- fitExtractVarPartModel(na.omit(as.matrix(cleanDat.b4e[,])), form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp.b4e <- sortCols(varPart.b4e,FUN=median,last= c("Residuals"))

pdf(file="4p13b4e.contributor_Fsplit_19sites1xPAVregr+SiteRegress_Protect_nothing-VariancePartition-PLASMA-7335x16677.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp.b4e, main="HDS 1.3ms - 4p13b4e - KNOWN APOE 19 sites 2x Regr(2PAV) + Site Regr, Prot Nothing" )

	SexSortOrder<-order(vp.b4e[["Sex"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b4e)) { vp.b4e[[i]]<-vp.b4e[[i]][SexSortOrder]; }
	rownames(vp.b4e)<-rownames(vp.b4e)[SexSortOrder]

	plotPercentBars( vp.b4e[1:50,]) + ggtitle( "Top Sex-covariates" )


	AgeSortOrder<-order(vp.b4e[["Age"]],decreasing=TRUE)
	#rownames(cleanDat.noNA.unreg)[SexSortOrder][1:50]
	for (i in ls(vp.b4e)) { vp.b4e[[i]]<-vp.b4e[[i]][AgeSortOrder]; }
	rownames(vp.b4e)<-rownames(vp.b4e)[AgeSortOrder]

	plotPercentBars( vp.b4e[1:50,]) + ggtitle( "Top Age-covariates" )


        BatchSortOrder<-order(vp.b4e[["RegrBloodPreanalyticFactor.HNRNPA2B1"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4e)) { vp.b4e[[i]]<-vp.b4e[[i]][BatchSortOrder]; }
        rownames(vp.b4e)<-rownames(vp.b4e)[BatchSortOrder]

        plotPercentBars( vp.b4e[1:50,]) + ggtitle( "Top Time To Spin (HNRNPA2B1)-covariates" )


        BatchSortOrder<-order(vp.b4e[["RegrBloodPreanalyticFactor.HBZ"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4e)) { vp.b4e[[i]]<-vp.b4e[[i]][BatchSortOrder]; }
        rownames(vp.b4e)<-rownames(vp.b4e)[BatchSortOrder]

        plotPercentBars( vp.b4e[1:50,]) + ggtitle( "Top Preanalytical Factor 2 (HBZ)-covariates" )


        BatchSortOrder<-order(vp.b4e[["contributor_Fsplit"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4e)) { vp.b4e[[i]]<-vp.b4e[[i]][BatchSortOrder]; }
        rownames(vp.b4e)<-rownames(vp.b4e)[BatchSortOrder]

        plotPercentBars( vp.b4e[1:50,]) + ggtitle( "Top Contributor (site)-covariates" )


        BatchSortOrder<-order(vp.b4e[["APOE.E4carrier"]],decreasing=TRUE)
        #rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
        for (i in ls(vp.b4e)) { vp.b4e[[i]]<-vp.b4e[[i]][BatchSortOrder]; }
        rownames(vp.b4e)<-rownames(vp.b4e)[BatchSortOrder]

        plotPercentBars( vp.b4e[1:50,]) + ggtitle( "Top APOE e4 carrier (Binary, no NA)-covariates" )


#	BatchSortOrder<-order(vp.b4e[["MMSE"]],decreasing=TRUE)
#	#rownames(cleanDat.noNA.unreg)[BatchSortOrder][1:50]
#	for (i in ls(vp.b4e)) { vp.b4e[[i]]<-vp.b4e[[i]][BatchSortOrder]; }
#	rownames(vp.b4e)<-rownames(vp.b4e)[BatchSortOrder]
#
#	plotPercentBars( vp.b4e[1:50,]) + ggtitle( "MMSE Cog. Score-covariates" )

dev.off()

varPart.b4e<-vp.b4e
saveRDS(varPart.b4e,"4p13b4e.contributor_Fsplit_19sites1xPAVregr+SiteRegress_Protect_nothing.na.omit-varPart.b4e.RDS")


######################################
## ANOVA + Volcanoes + DEXstacked Barplots (b4e)
source("parANOVA.dex.fallback7.25.R")

parallelThreads=8
outFilePrefix="4p13b4e.regr.knownAPOEonly2xRegr"
outFileSuffix="SixDxVolcs.CTimputed.3sitesMappedDx"

cleanDat<-cleanDat.b4e

#Grouping=numericMeta$Group #***
Grouping=regvars.b4c$Group.withCTimputed
Grouping[which(regvars.b4c$Group.withCTimputed=="CI.Other")]<-NA
Grouping[which(regvars.b4c$Group.withCTimputed=="AsymAD")]<-NA
#Grouping[which(regvars.b4c$Group.withCTimputed=="MCI")]<-NA

ANOVAout <- parANOVA.dex()

#*** commented new for 3-genotype only comparisons
##head(ANOVAout)
## plot volcanoes only for comparisons with 2/3 Group components in common.
#comps<-cbind(ANOVAout.colnames=colnames(ANOVAout), as.data.frame(do.call(rbind,strsplit(colnames(ANOVAout),"[.|-]"))))
#lastComp=which(grepl("diff ",comps$ANOVAout.colnames))[1]-1
#lastComp
#155  # 155-2 = 153 comparisons pairwise
## keep comparisons with 2/3 variables the same (n=45 comparisons)
#selectComps=unlist( sapply(3:lastComp,function(x) if(length(which(duplicated(unlist(comps[x,2:7]))))>1) x ) )
#length(selectComps) # selectComps is used directly by plotVolc to choose which comparisons to plot (column #s of p values in ANOVAout)
## Decide which comparison X axis (logFC) to flip
#comps[selectComps,1] # view pairwise comparisons selected for volcano plotting
# None need flipping (more extreme expected phenotype is first (numerator) in each case)
flip=c(3,4,5,6,10,11)
sameScale=FALSE #TRUE
symbolsOnly=TRUE
highlightGeneProducts=c("HNRNPA2B1","HBZ","SPC25","CPLX2","PTN","MAPT","OMG","GDF15","NEFL","NRGN","CTHRC1","HTRA1","TTR","LRP1","NTN1","SFRP1","GPNMB","PAFAH1B3")
labelHighlighted=TRUE      # if true, highlighted spots get text labels with their rownames from ANOVAout
labelTop=5
plotVolc()                 # runs on ANOVAout as input (need not be specified).

#DEXpercentStacked()        # runs on prior function outputs as input; writes stacked bar plot(s) to PDF.


Grouping=as.integer(regvars.b4c$APOE4.carrier)-1
Grouping[Grouping==0]<-"non.carrier"
Grouping[Grouping==1]<-"e4.carrier"
outFileSuffix="e4carrier.vs.NonCarrier"

ANOVAout <- parANOVA.dex()  # current .csv output has labels flipped (forgot to subtract 1 from factor levels above)

flip=c(3)
#flip=c()
#colnames(ANOVAout)[c(3:4)]<-c("e4.carrier-non.carrier","diff e4.carrier-non.carrier")
#16 p values are 0! (and FDR too)
# Set to 1e-175 (smaller than 1e-267, first nonzero sorted value)
ANOVAout[which(ANOVAout[,3]==0),3]<- 1e-275
sameScale=FALSE #TRUE
symbolsOnly=TRUE
highlightGeneProducts=c("HNRNPA2B1","HBZ","SPC25","CPLX2","PTN","MAPT","OMG","GDF15","NEFL","NRGN","CTHRC1","HTRA1","TTR","LRP1","NTN1","SFRP1","GPNMB","PAFAH1B3")
labelHighlighted=TRUE      # if true, highlighted spots get text labels with their rownames from ANOVAout
labelTop=20
plotVolc()                 # runs on ANOVAout as input (need not be specified).


#####################################################
## Output b4c cleanDat (subset, homozygotes of defined diagnosis) for Shijia

b4c.e33.subset.idx=which(regvars.b4c$Group.withCTimputed=="CT" & regvars.b4c$APOE.mapped=="e3/e3")
length(b4c.e33.subset.idx)
#2587

b4c.e44.subset.idx=which(regvars.b4c$Group.withCTimputed %in% c("CT","AsymAD","CI.Other","AD","MCI") & regvars.b4c$APOE.mapped=="e4/e4")
length(b4c.e44.subset.idx)
#510

numericMeta.Shijia<-regvars.b4c[sort(c(b4c.e33.subset.idx,b4c.e44.subset.idx)),  # 2587+510, kept in order of contributor site (A-T)
                                c("person_id","visit","sequential_visit_number","age_at_visit","Sex.int","raceAA","weight_kg","bmi","height_cm","years_of_education",  # demographics (longitudinal); raceAA is binary for African American (1); Sex.int is integer 0/1 F/M
                                  "contributor_code","contributor_Fsplit",       # Contributor Sites (regressed including site F split to F1, F2, F3)
                                  "resting_heart_rate_pulse","systolic_blood_pressure_sitting","diastolic_blood_pressure_sitting",                                     # heart health, numeric parameters
                                  "hypertension","stroke","tia","tbi","diabetes","chf","copd","mi","afib","hyperlipidaemia","depression","anxiety", "alcohol_hx","smoking_hx","total_years_smoked",  # comorbidities (binary status, except years smoked)
                                  "MMSE","cdr",           # Cognitive scores; preferred MMSE imputed from MoCA with education adjustment using Fasnacht et al 2022 rubric
                                  "Group.withCTimputed",  # diagnosis, including those mapped from our Nat Aging cohorts, and Controls with 0 CDR and MMSE >=28 imputed
                                  "C9Orf72","GRN","MAPT", # selected genetic binary mutation status for detrimental mutations
                                  "sample_matrix",        # citrate or EDTA tube-collected plasma
                                  "Lilly.BH.blood.pTau217","UDS.blood.pTau217","AmyloidPositivity.withRM",  # Pathological biomarkers for Nat Aging cohorts
                                  "APOE.mapped")]          # APOE genotype in "e#/e#" string format -- as provided by GNPC HDS v1.3, and mapped for missing values in Nat Aging 3 cohorts [ROSMAP (site R), UDS (site D), and BioHermes (site A)]

cleanDat.Shijia<-cleanDat.b4c[, sort(c(b4c.e33.subset.idx,b4c.e44.subset.idx))]
dim(cleanDat.Shijia)
#[1] 7335 3097
dim(numericMeta.Shijia)
#[1] 3097   41
any(!colnames(cleanDat.Shijia)==rownames(numericMeta.Shijia))
#[1] FALSE

save(cleanDat.Shijia,numericMeta.Shijia,file="APOE.homozygote_matrices-forShijia.RData")  # Z:/EBD/APOE.homozygote_matrices-forShijia.RData


########################################################
## Work toward a Nested CV prediction of 5715 missing APOE genotypes
## 1) get top predictive features for each of 6 genotypes in single binary prediction function, with top feature for the ensemble ranking as output added to the binary predictor function above.
## 2) Include top x features for each genotype in a 6-genotype prediction function, which has 80/20 Nested CV. The inner 5-fold x 5-rep x3-method ensemble -foldCV runs on 80% of data,
##    and the outer 20% uses the ensemble prediction to demonstrate little to no overfitting. 1000 iterations recommended.

##############################################################################
## 1) Binary predictor function with ensemble feature importance ranking output
##############################################################################
library(caret)
library(glmnet)
library(xgboost)
library(ranger)
library(progressr)
library(dplyr)
library(future)
library(doFuture)
library(doRNG)

fit_APOE_binary.95Acc <- function(expr, APOE_gt,
                                  nfold = 5, nrep = 5,
                                  ncores = parallel::detectCores() - 1,
                                  target=c("e4/e4"),   # run one genotype vs all others
                                  target_ppv=0.95,
                                  seed   = 1) {

  set.seed(seed)

  if (!"package:cli" %in% search()) suppressPackageStartupMessages(library(cli))

  memLimit=4*1024^3
  options(future.globals.maxSize= memLimit)  #4GB Total size of all global objects that need to be exported - up from 500MB
  Sys.setenv(R_FUTURE_GLOBALS_MAXSIZE=memLimit) #inherited by workers

  ## Helper function - threshold minimum
  pick_thr <- function(prob, truth, target_ppv = 0.95,
                       min_tp = 30, floor = 0.80) {
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

    ## ---------- preprocessing  (no kNN, already imputed) --------------
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
                 nthread = 1)
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
                     colsample_bytree = 0.8, nthread = ncores)
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
    imp_scaled <- scale(imp_mat)
    imp_mean   <- rowMeans(imp_scaled)
    top_feats  <- sort(imp_mean, decreasing = TRUE) #[1:top_k]

  ## -------- export top_k predictive proteins --------------------------
    if (!exists("rankedProteins", envir=.GlobalEnv)) assign("rankedProteins",list(), envir=.GlobalEnv)
    rankedProteins[[target]] <<- data.frame(feature=names(sort(top_feats, decreasing = TRUE)), importance=sort(top_feats, decreasing = TRUE))
    names(top_feats)<-gsub("\\|","_",names(top_feats))
    print(knitr::kable(data.frame(Protein = names(top_feats)[c(1:10,(top_k-9):top_k)],
                                  Importance = round(top_feats,3)[c(1:10,(top_k-9):top_k)]),
                       caption = sprintf("Top 10 and bottom 10 (of top %d) predictive proteins - (all exported to list rankedProteins)", top_k)))
  ### ------------------------------------------------------------------- ###


  ## ------------------------------------------------------------------------
  ## 4.  prediction wrapper -------------------------------------------------
  ## ------------------------------------------------------------------------
    function(new_expr) {

      new_expr <- prep_expr(rbind(expr[1,,drop=FALSE], new_expr))[-1,,drop=FALSE]
      new_expr <- new_expr[, colnames(X_all), drop = FALSE]

      glm_tg <- ovr[[target]]$glm
      xgb_tg <- ovr[[target]]$xgb
      rf_tg  <- ovr[[target]]$rf
      thr    <- ovr[[target]]$thr

      p_g <- drop(predict(glm_tg, new_expr, s = "lambda.min", type="response"))
      p_x <- drop(predict(xgb_tg, new_expr))
      p_r <- predict(rf_tg, data.frame(new_expr))$predictions[,"pos"]

      p_bin <- (p_g + p_x + p_r) / 3

      factor(ifelse(p_bin >= thr, target, NA_character_), levels = target)
    }
  }) # with_progress
}


################[QUICK TEST RUN - 28 THREADS OK ON SMALL TRAINING SET]####################
#option 3 - impute::impute.knn and reduce input expr mat size
training.cleanDat <- impute::impute.knn(cleanDat.4p13b3[, which(!is.na(gt.APOE))])$data   #  no NA
#sample.subset<-sort(sample(length(na.omit(gt.APOE)),750))
#table(na.omit(gt.APOE)[sample.subset])
set.seed(1)
# sample 400 of each genotype for training (n=810 total)
keep_idx<-unlist(lapply(split(seq_along( na.omit(gt.APOE)), na.omit(gt.APOE)), function(i) if (length(i) > 400) sample(i,400) else i), use.names=FALSE)
table(na.omit(gt.APOE)[keep_idx])  # 66 2/2, now 400 ea of the other 5 genotypes; took <5 min to complete training on 28/32 cores at 3.2GHz (128 GB RAM rec.)

predict_APOE.b3.allAndMapped.e44only.95Acc<-fit_APOE_binary.95Acc(t(training.cleanDat[APOE4.assays.keep[c(1:250,2251:length(APOE4.assays.keep))], keep_idx]), na.omit(gt.APOE)[keep_idx], ncores=28, target="e4/e4")
# CV e4/e4  -  Precision 0.983 ± 0.015 | Recall 1.000 ± 0.000
#
# Impute APOE the the above predict wrapper function output  in the 7k samples missing it (step 2 of 2)
#imputed.APOE.b3.allAndMapped.e44only.95Acc<-predict_APOE.b3.allAndMapped.e44only.95Acc(t(impute::impute.knn(cleanDat.4p13b3)$data[APOE4.assays.keep[c(1:250,2251:length(APOE4.assays.keep))], which(is.na(gt.APOE))]))
################[QUICK TEST RUN - 28 THREADS - 	COMPLETED]####################
table(imputed.APOE.b3.allAndMapped.e44only.95Acc)
#overprediction of e4/e4 (exptected 170 or less)
# 445


## e4/e4 get top 200 features with binary genotype predictor ensemble ML to get very high accuracy calls predicting e4/e4 (binary) in test sets of 5x5 fold CV

# done above:
#training.cleanDat.noNA<-t(na.omit(t(cleanDat.4p13b3[which(!rownames(cleanDat.4p13b3) %in% c("IRF6|O14896")), which(!is.na(gt.APOE))])))
#dim(training.cleanDat.noNA)
##  7334 14758   # 7335 with IRF6 in; without mapping our 3 cohort genotypes: 13004 (site K does not have APOE genotypes)
#training.gt.APOE<-gt.APOE[colnames(training.cleanDat.noNA)]
predict_APOE.b3.allAndMapped.e44only.95Acc<-fit_APOE_binary.95Acc(t(training.cleanDat.noNA), training.gt.APOE, ncores=14, target="e4/e4")
#CV e4/e4 - Precision 0.972 ± 0.025 | Recall 1.000 ± 0.000

# Impute APOE e4/e4 genotype in the 5715 (5k) samples missing it (step 2 of 2)
imputed.APOE.b3.allAndMapped.e44only.95Acc<-predict_APOE.b3.allAndMapped.e44only.95Acc(t(impute::impute.knn(cleanDat.4p13b3)$data[rownames(training.cleanDat.noNA),which(is.na(gt.APOE))]))
table(imputed.APOE.b3.allAndMapped.e44only.96Acc)
#current run
#e4/e4
#   70
#compared to before
table(imputed.APOE.b3.allAndMapped.noNA)
#e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#   12   578   101  3229  1625   170

names(imputed.APOE.b3.allAndMapped.e44only.95Acc)<-rownames(t(cleanDat.4p13b3[rownames(training.cleanDat.noNA),which(is.na(gt.APOE))]))
saveRDS(predict_APOE.b3.allAndMapped.e44only.95Acc,"predict_APOE.b3.allAndMapped.e44only.95Acc.RDS")

saveRDS(imputed.APOE.b3.allAndMapped.e44only.95Acc,"imputed.APOE.b3.allAndMapped.e44only.95Acc.RDS")

#predict_APOE.b3.all.e33e44.95Acc<-readRDS("predict_APOE.b3.allAndMapped.e33e44.95Acc.RDS")


## Run on Windows02 (fails with error:
# e2/e2 and e2/e4 get ranked features with binary genotype predictor ensemble ML in training set of 14758 samples with known APOE genotype
for (this.genotype in c("e2/e2","e2/e4")) {
  this.geno.numeric=gsub("/","",gsub("e","",this.genotype))
  wrapper.out<-fit_APOE_binary.95Acc(t(training.cleanDat.noNA), training.gt.APOE, ncores=8, target=this.genotype)
  assign(paste0("predict_APOE.b3.allAndMapped.e",this.geno.numeric,"only.95Acc"), wrapper.out)

  # Impute APOE single genotype in the 5715 (5k) samples missing it (step 2 of 2)
  predict.geno.out<-wrapper.out(t(impute::impute.knn(cleanDat.4p13b3)$data[rownames(training.cleanDat.noNA),which(is.na(gt.APOE))]))
  assign(paste0("imputed.APOE.b3.allAndMapped.e",this.geno.numeric,"only.95Acc"), predict.geno.out)
  table(predict.geno.out)
} # check rankedProteins and binaryPredictionMetrics lists in global environment.


## Run on Windows03
# 4 genotypes: get ranked features with binary genotype predictor ensemble ML in training set of 14758 samples with known APOE genotype
for (this.genotype in c("e4/e4","e3/e4","e2/e3","e3/e3","e2/e2","e2/e4")) {
  this.geno.numeric=gsub("/","",gsub("e","",this.genotype))
  wrapper.out<-fit_APOE_binary.95Acc(t(training.cleanDat.noNA), training.gt.APOE, ncores=14, target=this.genotype)
  assign(paste0("predict_APOE.b3.allAndMapped.e",this.geno.numeric,"only.95Acc"), wrapper.out)

  # Impute APOE single genotype in the 5715 (5k) samples missing it (step 2 of 2)
  predict.geno.out<-wrapper.out(t(impute::impute.knn(cleanDat.4p13b3)$data[rownames(training.cleanDat.noNA),which(is.na(gt.APOE))]))
  assign(paste0("imputed.APOE.b3.allAndMapped.e",this.geno.numeric,"only.95Acc"), predict.geno.out)
  cat(table(predict.geno.out))
}


#CV e4/e4 - Precision 0.972 ± 0.025 | Recall 1.000 +/- 0.000
#CV e3/e4 - Precision 0.969 ± 0.007 | Recall 1.000 +/- 0.000
#CV e2/e3 - Precision 0.969 ± 0.016 | Recall 1.000 +/- 0.000
#CV e3/e3 - Precision 0.968 ± 0.006 | Recall 1.000 +/- 0.000
#CV e2/e2 - Precision 1.000 ± 0.000 | Recall 1.000 +/- 0.000
#CV e2/e4 - Precision 1.000 ± 0.000 | Recall 1.000 +/- 0.000
# table(predict.geno.out)  or  table(imputed.APOE.b3.allAndMapped.e44only.95Acc)
#e4/e4: 70
#e3/e4: 943
#e3/e3: 2840
#e2/e3: 134
#e2/e2: 0
#e2/e4: 1

# saved data for each binary prediction run:
#binaryPredictionMetrics
#rankedProteins

allImportantFeatures<-unique(unlist(lapply(rankedProteins, function(x) x[which(x[,"importance"]>0),"feature"])))
length(allImportantFeatures)
# 2235

saveRDS(rankedProteins,"rankedProteins.list.RDS")
saveRDS(binaryPredictionMetrics,"binaryPredictionMetrics.list.RDS")

names(imputed.APOE.b3.allAndMapped.e44only.95Acc)<-names(imputed.APOE.b3.allAndMapped.e34only.95Acc)<-names(imputed.APOE.b3.allAndMapped.e23only.95Acc)<-names(imputed.APOE.b3.allAndMapped.e33only.95Acc)<-names(imputed.APOE.b3.allAndMapped.e22only.95Acc)<-names(imputed.APOE.b3.allAndMapped.e24only.95Acc) <- rownames(t(cleanDat.4p13b3[rownames(training.cleanDat.noNA),which(is.na(gt.APOE))]))

saveRDS(predict_APOE.b3.allAndMapped.e44only.95Acc,"predict_APOE.b3.allAndMapped.e44only.95Acc.RDS")
saveRDS(imputed.APOE.b3.allAndMapped.e44only.95Acc,"imputed.APOE.b3.allAndMapped.e44only.95Acc.RDS")

saveRDS(predict_APOE.b3.allAndMapped.e34only.95Acc,"predict_APOE.b3.allAndMapped.e34only.95Acc.RDS")
saveRDS(imputed.APOE.b3.allAndMapped.e34only.95Acc,"imputed.APOE.b3.allAndMapped.e34only.95Acc.RDS")

saveRDS(predict_APOE.b3.allAndMapped.e23only.95Acc,"predict_APOE.b3.allAndMapped.e23only.95Acc.RDS")
saveRDS(imputed.APOE.b3.allAndMapped.e23only.95Acc,"imputed.APOE.b3.allAndMapped.e23only.95Acc.RDS")

saveRDS(predict_APOE.b3.allAndMapped.e33only.95Acc,"predict_APOE.b3.allAndMapped.e33only.95Acc.RDS")
saveRDS(imputed.APOE.b3.allAndMapped.e33only.95Acc,"imputed.APOE.b3.allAndMapped.e33only.95Acc.RDS")

saveRDS(predict_APOE.b3.allAndMapped.e22only.95Acc,"predict_APOE.b3.allAndMapped.e22only.95Acc.RDS")
saveRDS(imputed.APOE.b3.allAndMapped.e22only.95Acc,"imputed.APOE.b3.allAndMapped.e22only.95Acc.RDS")

saveRDS(predict_APOE.b3.allAndMapped.e24only.95Acc,"predict_APOE.b3.allAndMapped.e24only.95Acc.RDS")
saveRDS(imputed.APOE.b3.allAndMapped.e24only.95Acc,"imputed.APOE.b3.allAndMapped.e24only.95Acc.RDS")


#save.image("4p13b4f.ShijiaOutput+6binaryClassifiersAPOE.RData")


##############################################################################
##  2) 6-genotype Parallel fit_APOE_ensemble()  ------------------------------
##############################################################################
library(caret)
library(glmnet)
library(xgboost)
library(ranger)
library(progressr)
#library(VIM)  ## VIM kNN disabled (impute input to the function)
library(dplyr)
#library(doParallel)  ## doRNG for future instead
library(future)
library(doFuture)
library(doRNG)
#library(DMwR2)  # DMwR package had SMOTE function;
library(performanceEstimation)  # For smote function, oversampling rare genotypes with synthetic samples


fit_APOE_ensemble_6way_nestedCV <- function(expr, APOE_gt, confidenceRules=FALSE,
                                  nfold = 5, nrep = 5, outerFolds = 5,
                                  ncores = parallel::detectCores() - 1,
                                  seed   = 1, memLimitGB = 4) {

  start_time <- Sys.time()  # Record start time

  if (!"package:cli" %in% search()) suppressPackageStartupMessages(library(cli))

  memLimit=memLimitGB*1024^3
  options(future.globals.maxSize= memLimit)      # 4GB Total size of all global objects that need to be exported - up from 500MB
#  Sys.setenv(R_FUTURE_GLOBALS_MAXSIZE=memLimit)  # inherited by workers  -- may be sticky and not able to reset without restarting R!

  set.seed(seed)

  # Store original feature names
  original_names <- colnames(expr)
  # Create a data frame with modified names
  modified_names <- colnames(expr) <- make.names(original_names, unique = TRUE)
  # Create a mapping between original and modified names
  name_mapping <- setNames(original_names, modified_names)


#  ## -------- 1 OUTER hold-out ----------------------------------------- ###
#  outer_idx <- createDataPartition(APOE_gt, p = .80, list = FALSE)  ### NEW
#  expr_tr   <- expr[ outer_idx, , drop = FALSE]                     ### NEW
#  expr_hold <- expr[-outer_idx, , drop = FALSE]                     ### NEW
#  y_tr      <- APOE_gt[ outer_idx]
#  y_hold    <- APOE_gt[-outer_idx]

  ## scaling function for outer fold data
  prep_outer <- function(m, ref = NULL) {
    keep <- which(colMeans(is.na(m)) <= .20)
    if (is.null(ref)) {
      x <- scale(m[ , keep, drop = FALSE])
      list(x = x,
           center = attr(x, "scaled:center"),
           scale  = attr(x, "scaled:scale"),
           vars   = colnames(x))
    } else {
      x <- scale(m[ , ref$vars, drop = FALSE],
                 center = ref$center,
                 scale  = ref$scale)
      x[is.na(x)] <- 0
      list(x = x)
    }
  }

  ## -------- Flexible OUTER hold-out (for loop)------------------------- ###
  if (outerFolds <= 1) {  # 1 fold invalid for caret helper; treat same as 0
    outer_folds <- list(integer(0))
  } else {
#    createFolds(APOE_gt, k=outerFolds, returnTrain=FALSE)  # was getting ~7400 samples for holdout set after 50 synthetics added (2x 20%)
    # Create a list to store hold-out indices for each iteration
    outer_folds <- vector("list", outerFolds)
    set.seed(1)  # Ensure reproducibility
#    for (i in seq_len(outerFolds)) {
#      # Generate a unique 20% hold-out sample
#      hold_idx1 <- createDataPartition(APOE_gt, p = 0.2, list = FALSE)
#      # Store the hold-out indices
#      outer_folds[[i]] <- hold_idx1
#    }
    outer_folds <- replicate( outerFolds, createDataPartition(APOE_gt, p=0.20, list=FALSE)[,1], simplify=FALSE)
  }
  ## containers sized to the number of outer iterations
  outer_results   <- vector("list", length(outer_folds))   # collect metrics
  outer_hold_perf <- vector("list", length(outer_folds))   # collect hold scores
#  outer_rankedProteins <- vector("list", length(outer_folds))   # collect features, importances

  for (outer in seq_along(outer_folds)) {      # <-- nested CV OUTER FOLD LOOP
    hold_idx  <- outer_folds[[outer]]
    train_idx <- setdiff(seq_len(nrow(expr)), hold_idx)  # <- full data if hold_idx=0

    expr_tr   <- expr[train_idx , , drop = FALSE]
    expr_hold <- expr[hold_idx  , , drop = FALSE]
    y_tr.inn  <- APOE_gt[train_idx]
    y_hold    <- APOE_gt[hold_idx]

    if (outerFolds > 1) {
      # Check genotype distribution in hold-out data
      genotype_counts <- table(APOE_gt)
      # Identify rare genotypes with fewer than 8 samples
      rare_genotypes <- names(genotype_counts[genotype_counts < 450])  # e2/e2: 66; e2/e4: 405

      if (length(rare_genotypes) > 0) {
        # Apply SMOTE to full dataset
        expr_balanced <- performanceEstimation::smote(y ~ ., data = data.frame(y=APOE_gt,expr), perc.over = 200, k=5, perc.under = 100)

        # Extract synthetic samples for missing genotypes in hold_data
        synthetic_samples<-list()
        for (rare_genotype in rare_genotypes) {
          synthetic_samples[[rare_genotype]] <- expr_balanced[expr_balanced$y %in% rare_genotype, ]
          # keep 15 (2/2) or 80 (2/4) of the synthetic samples to add to holdout data
          keepCount=25 # overwritten in the case of either specific genotype below.
#          if (rare_genotype=="e2/e2") { keepCount=15 } else if (rare_genotype=="e2/e4") { keepCount=80 }
          synthetic_samples[[rare_genotype]] <- synthetic_samples[[rare_genotype]][sample(nrow(synthetic_samples[[rare_genotype]]),keepCount),]
        }
        # Add synthetic samples to hold_data
        synthetic_samples<-do.call(rbind,synthetic_samples)
        expr_hold <- as.matrix(rbind(expr_hold, synthetic_samples[, -1]))  # Remove genotype column; matrix expected by predict functions.
      }
    }
#expr_hold.global<<-expr_hold  # check number of samples (rows) -- should be 20% of input + 50 synthetics.

    # Scale the holdout data identical to training (inner CV) data
    pp <- prep_outer(expr_tr)   # Learn scaling on the outer training part

    X_tr_outer <- pp$x          # Use this for glmnet / xgboost / RF
    X_hd_outer <- prep_outer(expr_hold, pp)$x  # Use this for prediction

    ## ------------------------------------------------------------------------
    ## 0.  create inner CV resampling indices once
    ## ------------------------------------------------------------------------
    cvIndex <- createMultiFolds(y_tr.inn, k = nfold, times = nrep)  #APOE_gt <- y_tr if no outer CV
    nTasks  <- length(cvIndex)
    ## ------------------------------------------------------------------------
    ## 1.  start the cluster
    ## ------------------------------------------------------------------------
    handlers("progress")
    handlers(global=TRUE)  # set handler for progress bar before the cluster

    plan(multisession, workers = ncores)  # or multicore on linux/macOS
    registerDoFuture()

    ## ------------------------------------------------------------------------
    ## 2.  run each fold in a worker  -----------------------------------------
    ## ------------------------------------------------------------------------

    ## Helper function - catch GLMnet errors
       safe_glmnet <- function(x, y, w, ...) {
        # need >=2 obs per *present* class
        if (any(table(y) < 2)) return(NULL)
        tryCatch(
          cv.glmnet(x, y, family = "multinomial",
                    weights = w, type.measure = "class", parallel=FALSE, ...),
          error = function(e) NULL)
      }

    ## helper function - make all rectangular matrices carry all levels/genotypes
    complete_prob_mat <- function(mat, all_levels) {
      miss <- setdiff(all_levels, colnames(mat))
      if (length(miss)) {
        mat <- cbind(mat, matrix(0, nrow(mat), length(miss),
                                 dimnames = list(NULL, miss)))
      }
      mat[ , all_levels, drop = FALSE]           # reorder & keep only wanted
    }

    calc_prec <- function(pred, ref, cls) {
      tp <- sum(pred == cls & ref == cls, na.rm = TRUE)
      fp <- sum(pred == cls & ref != cls, na.rm = TRUE)
      if (tp + fp == 0) return(NA_real_)
      tp / (tp + fp)
    }

    calc_rec  <- function(pred, ref, cls) {
      tp <- sum(pred == cls & ref == cls, na.rm = TRUE)
      fn <- sum(pred != cls & ref == cls, na.rm = TRUE)
      if (tp + fn == 0) return(NA_real_)
      tp / (tp + fn)
    }

    # Genotype counts from training data
#    genotype_counts <- c("e2/e2" = 66, "e2/e3" = 1636,
#                         "e2/e4" = 405, "e3/e3" = 8490,
#                         "e3/e4" = 5171, "e4/e4" = 909)
    genotype_counts <- as.vector(table(APOE_gt))
    # Square root weighting - balanced, preventing excessive bias toward rare classes
    genotype_counts <- 1 / sqrt(genotype_counts)

    # Min-max scaling function (scaling between 1 and 15)
    min_val <- min(genotype_counts)
    max_val <- max(genotype_counts)

    scaled_weights <- 1 + (genotype_counts - min_val) / (max_val - min_val) * (15 - 1)

    # Assign scaled weights
#    class_w <- scaled_weights
    # Assign equal weights
    class_w <- rep(1,length(names(table(APOE_gt))))

    #print(class_w)  # if scaled:
    #15.00000  2.730127  5.844367  1.000000  1.380867  3.783445
    names(class_w) <- names(table(APOE_gt))


    with_progress({                              # << all progress lives here
      n_sub <- 3   # Progress bar increments 3x per fold
      p <- progressor(steps=length(cvIndex) * n_sub )     #along = cvIndex)          # one step per fold

      p(sprintf("Initializing %d workers...", ncores), amount=0)


    metrics <- foreach(fold = seq_along(cvIndex),
                       .combine   = rbind,
                       .export = "p",            # let workers see 'p'
                       .packages  = c("progressr","glmnet","xgboost","ranger","dplyr")) %dorng% {

      ## ----------------- announce fold start -----------------------------
      p(sprintf("fold %d/%d  -  started", fold, length(cvIndex)), amount=0)

      set.seed(seed + fold)                        # reproducible inside worker
      tr_idx <- cvIndex[[fold]]
      te_idx <- setdiff(seq_len(nrow(expr_tr)), tr_idx)

      ## ---------- preprocessing  (no kNN, already imputed) --------------
      prep <- function(m) { scale(m[, colMeans(is.na(m)) <= 0.20, drop = FALSE]) }
      X_tr <- prep(expr_tr[tr_idx, ])
      X_te <- scale(expr_tr[te_idx, colnames(X_tr)],
                    center = attr(X_tr, "scaled:center"),
                    scale  = attr(X_tr, "scaled:scale"))
      X_te[is.na(X_te)] <- 0

      levels_all <- sort(unique(APOE_gt))            # six possible APOE genotypes in FULL data
      num_class   <- length(levels_all)              # == 6

      y_tr <- factor(y_tr.inn[tr_idx], levels = levels_all)
      y_te <- factor(y_tr.inn[te_idx], levels = levels_all)

#      ## NEW-weights: one vector that all learners can consume
#      class_w <- c("e2/e2" = 15, "e2/e3" = 3,
#                   "e3/e3" = 1,                # 2, if one may want *very* high precision
#                   "e2/e4" = 12, "e3/e4" = 2,
#                   "e4/e4" = 8)                # low prevalence -> strong weight

      #For robustness, drop samples and labels whose labels are NA (should not have been input)
      keep_tr <- !is.na(y_tr)
      keep_te <- !is.na(y_te)

      X_tr <- X_tr[keep_tr, , drop = FALSE]
      y_tr <- y_tr[keep_tr]

      X_te <- X_te[keep_te, , drop = FALSE]
      y_te <- y_te[keep_te]

      ## ---------------- glmnet ---------------------------------------------
      w_tr <- class_w[as.character(y_tr)]
      if (any(table(y_tr) < 2)) message(sprintf("Fold %d - skipped glmnet (rare class)", fold))
      cv_glm <- safe_glmnet(X_tr, y_tr, w_tr)
      #          cv.glmnet(X_tr, y_tr,
      #                    family = "multinomial",
      #                    weights = w_tr,
      #                    type.measure = "class", parallel = FALSE)
  #    cv_glm <- cv.glmnet(X_tr, y_tr, family = "multinomial",
  #                        type.measure = "class", parallel = FALSE)
      p(sprintf("fold %d/%d  •  glmnet done", fold, length(cvIndex)), amount=1)

      if (is.null(cv_glm)) {
        ## give a zero-probability matrix when glmnet is absent
        p_glm <- matrix(0, nrow(X_te), num_class,
                        dimnames = list(NULL, levels_all))
      } else {
        p_glm <- predict(cv_glm, X_te, s = "lambda.min",
                         type = "response")[,,1]
        p_glm <- complete_prob_mat(p_glm, levels_all)
      }

      ## ---------------- xgboost --------------------------------------------
  #    dtr <- xgb.DMatrix(X_tr, label = as.numeric(y_tr) - 1)
      dtr  <- xgb.DMatrix(X_tr, label = as.numeric(y_tr) - 1, weight = w_tr)
      dval <- xgb.DMatrix(X_te, label = as.numeric(y_te) - 1)

      xpar <- list(eta = 0.1, max_depth = 6, subsample = 0.8,
                   colsample_bytree = 0.8,
                   objective = "multi:softprob",
                   num_class = num_class,     # length(levels(y_tr)),
                   nthread = 1,               # <- single thread per worker
                   eval_metric = "mlogloss")  # explicitly set eval metric to avoid warning

      bst  <- xgb.train(params=xpar, data=dtr, watchlist=list(train=dtr, eval=dval), nrounds = 200,
                        verbose = 0, early_stopping_rounds = 20)
      p_xgb <- matrix(predict(bst, X_te),
                      ncol = num_class, byrow = TRUE)
  #                    ncol = length(levels(y_tr)), byrow = TRUE)
  #    colnames(p_xgb) <- levels(y_tr)
      colnames(p_xgb) <- levels_all                 # <- name first
      p_xgb <- complete_prob_mat(p_xgb, levels_all) # <- then complete
      p(sprintf("fold %d/%d  •  xgboost done", fold, length(cvIndex)), amount=1)

      ## ---------------- random-forest --------------------------------------
  #    rf  <- ranger(y_tr ~ ., data = data.frame(y_tr, X_tr),
  #                  probability = TRUE, num.trees = 500,
  #                  num.threads = 1)          # <- single thread
      rf <- ranger(y_tr ~ ., data = data.frame(y_tr, X_tr, check.names=TRUE),
                   probability   = TRUE,
                   num.trees     = 500, num.threads = 1,
                   class.weights = class_w, verbose=FALSE)

      p(sprintf("fold %d/%d  •  ranger done", fold, length(cvIndex)), amount=1)

      p_rf <- predict(rf, data.frame(X_te, check.names=TRUE))$predictions
      colnames(p_rf) <- rf$forest$levels
      p_rf  <- complete_prob_mat(p_rf, levels_all)

      ## ---------------- ensemble vote --------------------------------------

      ## ----- averaged class-probabilities --------------------------------------
      probs_list <- list(p_xgb, p_rf)           # always available
      if (!is.null(cv_glm)) probs_list <- c(list(p_glm), probs_list)
      p_avg <- Reduce(`+`, probs_list) / length(probs_list)
      #p_avg <- (p_glm + p_xgb + p_rf) / 3          # n x 6  matrix

      top_class <- apply(p_avg, 1, function(z) names(z)[which.max(z)])
      top_prob  <- apply(p_avg, 1, max)
      second    <- apply(p_avg, 1, function(z) sort(z, TRUE)[2])

      ## confidence rules --------------------------------------------------------
      if(confidenceRules) {
        ##   1. if the winner is e3/e3 or e4/e4  AND  its prob >= 0.80  ? keep it
        ##   2. if the winner is anything else   AND  (prob >= 0.90  *and*
        ##                                             prob-margin >= 0.20) ? keep it
        ##   3. otherwise                                                             ? NA
        keep_homo  <- top_class %in% c("e3/e3", "e4/e4") & top_prob >= 0.80
        keep_other <- !(top_class %in% c("e3/e3", "e4/e4")) &
                      top_prob >= 0.90 & (top_prob - second) >= 0.20

        y_hat <- ifelse(keep_homo | keep_other, top_class, NA_character_)
      } else {
        y_hat = top_class
      }
        y_hat <- factor(y_hat, levels = levels_all)
  #    p_avg <- (p_glm + p_xgb + p_rf) / 3
  #    y_hat <- factor(colnames(p_avg)[max.col(p_avg)],
  #                    levels = levels(y_tr))

      cm <- tryCatch(
               caret::confusionMatrix(y_hat, y_te),
               error = function(e) NULL)

  ## ---------------- safe extractor --------------------------------------

      prec_e22 <- calc_prec (y_hat, y_te, "e2/e2")
      rec_e22  <- calc_rec  (y_hat, y_te, "e2/e2")
      prec_e23 <- calc_prec (y_hat, y_te, "e2/e3")
      rec_e23  <- calc_rec  (y_hat, y_te, "e2/e3")
      prec_e24 <- calc_prec (y_hat, y_te, "e2/e4")
      rec_e24  <- calc_rec  (y_hat, y_te, "e2/e4")
      prec_e33 <- calc_prec (y_hat, y_te, "e3/e3")
      rec_e33  <- calc_rec  (y_hat, y_te, "e3/e3")
      prec_e34 <- calc_prec (y_hat, y_te, "e3/e4")
      rec_e34  <- calc_rec  (y_hat, y_te, "e3/e4")
      prec_e44 <- calc_prec (y_hat, y_te, "e4/e4")
      rec_e44  <- calc_rec  (y_hat, y_te, "e4/e4")

      acc   <- if (!is.null(cm)) cm$overall["Accuracy"] else NA_real_
      macro <- if (!is.null(cm) && !is.null(dim(cm$byClass)))
                   mean(cm$byClass[,"F1"]) else NA_real_

      data.frame(Accuracy = acc,
                 MacroF1  = macro,
                 Prec_e33 = prec_e33,  Rec_e33 = rec_e33,
                 Prec_e44 = prec_e44,  Rec_e44 = rec_e44,
                 Prec_e23 = prec_e23,  Rec_e23 = rec_e23,
                 Prec_e34 = prec_e34,  Rec_e34 = rec_e34,
                 Prec_e24 = prec_e24,  Rec_e24 = rec_e24,
                 Prec_e22 = prec_e22,  Rec_e22 = rec_e22,
                 Fold     = fold)
    } # foreach
    p(sprintf("Finished evaluating inner fold set %d  •  Inner %d x %d folds DONE.", outer, nfold, nrep), amount=0)
    cat(sprintf("Finished evaluating inner fold set %d  •  Inner %d x %d folds DONE. Summary stats:\n\n", outer, nfold, nrep))
    }) # with_progress

  #  stopCluster(cl)                               # tidy up
  #  registerDoSEQ()                               # back to sequential

    cat(sprintf("CV Accuracy  %.3f ± %.3f\n",
                mean(metrics$Accuracy, na.rm=TRUE), sd(metrics$Accuracy, na.rm=TRUE)))
    cat(sprintf("CV Macro-F1  %.3f ± %.3f\n\n",
                mean(metrics$MacroF1, na.rm=TRUE),  sd(metrics$MacroF1, na.rm=TRUE)))
    cat(sprintf("Precision (e2/e2) %.3f; Recall %.3f\n",
               mean(metrics$Prec_e22, na.rm=TRUE), mean(metrics$Rec_e22, na.rm=TRUE)))
    cat(sprintf("Precision (e3/e3) %.3f; Recall %.3f\n",
               mean(metrics$Prec_e33, na.rm=TRUE), mean(metrics$Rec_e33, na.rm=TRUE)))
    cat(sprintf("Precision (e4/e4) %.3f; Recall %.3f\n",
               mean(metrics$Prec_e44, na.rm=TRUE), mean(metrics$Rec_e44, na.rm=TRUE)))
    cat(sprintf("Precision (e2/e3) %.3f; Recall %.3f\n",
               mean(metrics$Prec_e23, na.rm=TRUE), mean(metrics$Rec_e23, na.rm=TRUE)))
    cat(sprintf("Precision (e2/e4) %.3f; Recall %.3f\n",
               mean(metrics$Prec_e24, na.rm=TRUE), mean(metrics$Rec_e24, na.rm=TRUE)))
    cat(sprintf("Precision (e3/e4) %.3f; Recall %.3f\n",
               mean(metrics$Prec_e34, na.rm=TRUE), mean(metrics$Rec_e34, na.rm=TRUE)))
    cat(sprintf("--------------------------------------\n\nEvaluating outer fold %d unseen hold-out data...", outer))

    ## ---------------- evaluate once on the untouched hold-out ---------- ###
    levels_all <- sort(unique(APOE_gt))            # six possible APOE genotypes in FULL data
    num_class   <- length(levels_all)              # == 6

    if (length(hold_idx)) {                                    # (normal case, outerFolds>1)

      ## ---- fit an ensemble on expr_tr / y_tr.inn  (training part of outer fold)

#      ## NEW-weights: one vector that all learners can consume
#      class_w <- c("e2/e2" = 15, "e2/e3" = 3,
#                   "e3/e3" = 1,                # 2, if one may want *very* high precision
#                   "e2/e4" = 12, "e3/e4" = 2,
#                   "e4/e4" = 8)                # low prevalence -> strong weight
      w_tr=class_w[as.character(y_tr.inn)]

      y_tr.inn.factor<-factor(y_tr.inn, levels=levels_all)

      if (any(table(y_tr.inn) < 2)) message(sprintf("Outer Fold %d - skipped glmnet (rare class)", outer))
      outer_glm <- safe_glmnet(X_tr_outer, y_tr.inn.factor, w_tr)

      dtr.outer  <- xgb.DMatrix(X_tr_outer, label = as.numeric(y_tr.inn.factor) - 1, weight = w_tr)
#      dval.outer <- xgb.DMatrix(X_te, label = as.numeric(y_te) - 1)

      xpar.outer <- list(eta = 0.1, max_depth = 6, subsample = 0.8,
                   colsample_bytree = 0.8,
                   objective = "multi:softprob",
                   num_class = num_class,          # length(levels(y_tr)),
                   nthread = ncores,               # <- all threads for the outer evaluation
                   eval_metric = "mlogloss")       # explicitly set eval metric to avoid warning

#      bst  <- xgboost::xgb.train(params=xpar, data=dtr, watchlist=list(train=dtr, eval=dval), nrounds = 200,
#                        verbose = 0, early_stopping_rounds = 20)

      bst_all   <- xgboost::xgb.train(params=xpar.outer, data = dtr.outer, # expr_tr,
#                                    label = as.numeric(factor(y_tr.inn)) - 1,
                                    nrounds = 200,
                                    verbose = 0)

      rf_all    <- ranger(y_tr.inn.factor ~ ., data = data.frame(X_tr_outer, y_tr.inn.factor, check.names=TRUE),
                          probability = TRUE, num.trees = 500, class.weights = class_w, verbose=FALSE)

#      prep <- function(m) { scale(m[, colMeans(is.na(m)) <= 0.20, drop = FALSE]) }
#      new_expr <- prep(expr_hold)
#      new_expr <- X_hd_outer

      if (is.null(outer_glm)) {
        ## give a zero-probability matrix when glmnet is absent
        p_g <- matrix(0, nrow(X_hd_outer), num_class,
                        dimnames = list(NULL, levels_all))
      } else {
        p_g <- predict(outer_glm, X_hd_outer, s = "lambda.min",
                         type = "response")[,,1]
        p_g <- complete_prob_mat(p_g, levels_all)
      }
      p_x <- drop(matrix(predict(bst_all, X_hd_outer), ncol=num_class, byrow=TRUE))
      p_x <- complete_prob_mat(p_x, levels_all)
      p_r <- predict(rf_all, data.frame(X_hd_outer, check.names=TRUE))$predictions
      p_r <- complete_prob_mat(p_r, levels_all)
      #p_avg <- (p_g + p_x + p_r)/3
      probs_list <- list(p_x, p_r)           # always available
      if (!is.null(outer_glm)) probs_list <- c(list(p_g), probs_list)
      p_avg <- Reduce(`+`, probs_list) / length(probs_list)

      ## Apply ensemble voting
    #  factor(ifelse(p_bin >= ovr[[target]]$thr, target, NA_character_),
    #         levels = target)       # using a minimum accuracy threshold percentile target.
#      hold_pred <- factor(colnames(p_avg)[max.col(p_avg)], levels=levels_all)  # not using a minimum accuracy threshold percentile target.
      hold_pred <- apply(p_avg, 1, function(z) names(z)[which.max(z)])
      hold_pred <- factor(hold_pred, levels=levels_all)
### ready to remove later
      tp <- fp <- fn <- prec_vec <- rec_vec <- setNames(numeric(length(levels_all)), levels_all)
      ## --------------------------------  TP / FP / FN  ---------------------------
      for (cls in levels_all) {
        tp[cls]   <- sum(hold_pred == cls & y_hold == cls)
        fp[cls]   <- sum(hold_pred == cls & y_hold != cls)
        fn[cls]   <- sum(hold_pred != cls & y_hold == cls)
        prec_vec[cls] <- calc_prec(hold_pred, y_hold, cls)  #if (tp[cls] + fp[cls]) tp[cls]/(tp[cls]+fp[cls]) else 0 #NA
        rec_vec [cls] <- calc_rec (hold_pred, y_hold, cls)  #if (tp[cls] + fn[cls]) tp[cls]/(tp[cls]+fn[cls]) else 0 #NA
      }
###
      # Compute precision and recall using safe extraction functions
      prec_e22 <- calc_prec(hold_pred, y_hold, "e2/e2")
      rec_e22  <- calc_rec(hold_pred, y_hold, "e2/e2")
      prec_e23 <- calc_prec(hold_pred, y_hold, "e2/e3")
      rec_e23  <- calc_rec(hold_pred, y_hold, "e2/e3")
      prec_e24 <- calc_prec(hold_pred, y_hold, "e2/e4")
      rec_e24  <- calc_rec(hold_pred, y_hold, "e2/e4")
      prec_e33 <- calc_prec(hold_pred, y_hold, "e3/e3")
      rec_e33  <- calc_rec(hold_pred, y_hold, "e3/e3")
      prec_e34 <- calc_prec(hold_pred, y_hold, "e3/e4")
      rec_e34  <- calc_rec(hold_pred, y_hold, "e3/e4")
      prec_e44 <- calc_prec(hold_pred, y_hold, "e4/e4")
      rec_e44  <- calc_rec(hold_pred, y_hold, "e4/e4")


      accuracy <- sum(hold_pred == y_hold) / length(y_hold)

      ## Macro-F1 calc
      #safe_prec <- function(tp, fp) if (tp + fp == 0) 0 else tp / (tp + fp)
      #safe_rec  <- function(tp, fn) if (tp + fn == 0) 0 else tp / (tp + fn)
      #safe_f1   <- function(p , r ) if (p +  r  == 0) 0 else 2 * p * r / (p + r)
      safe_f1 <- function(p, r) if (is.na(p) || is.na(r) || (p + r) == 0) { NA_real_ } else { 2 * p * r / (p + r) }
#      prec_vec <- mapply(safe_prec, tp, fp)
#      rec_vec  <- mapply(safe_rec , tp, fn)
      f1_vec   <- mapply(safe_f1  , prec_vec, rec_vec)

      macro_f1 <- mean(f1_vec, na.rm=TRUE)          # always six numbers -> no accidental NA-drop


      ## store them together with the already-reported metrics
      outer_hold_perf[[outer]] <- data.frame(
#        Genotype  = levels_all,
        Prec_e33  = prec_e33,  Rec_e33 = rec_e33,
        Prec_e44  = prec_e44,  Rec_e44 = rec_e44,
        Prec_e23  = prec_e23,  Rec_e23 = rec_e23,
        Prec_e34  = prec_e34,  Rec_e34 = rec_e34,
        Prec_e24  = prec_e24,  Rec_e24 = rec_e24,
        Prec_e22  = prec_e22,  Rec_e22 = rec_e22,
        Accuracy  = accuracy,
        MacroF1   = macro_f1,
        OuterFold = outer, stringsAsFactors=FALSE
      )

      # Print results for each genotype
      cat(paste0("\n---- 20 % hold-out (outer fold ", outer, ") ----\n"))
      for (cls in levels_all) {
        cat(sprintf("%-5s  Precision %.3f | Recall %.3f\n", cls, calc_prec(hold_pred, y_hold, cls), calc_rec(hold_pred, y_hold, cls)))
      }
    } else {                                                     # (inner-only case, outerFolds=0, or 1)
#      outer_hold_perf[[outer]] <- data.frame(Genotype=levels_all, TP=NA_real_, FP=NA_real_, FN=NA_real_, Precision = NA_real_, Recall = NA_real_, OuterFold = outer)
     outer_hold_perf[[outer]] <- data.frame(
#        Genotype  = levels_all,
        Prec_e33  = NA_real_,  Rec_e33 = NA_real_,
        Prec_e44  = NA_real_,  Rec_e44 = NA_real_,
        Prec_e23  = NA_real_,  Rec_e23 = NA_real_,
        Prec_e34  = NA_real_,  Rec_e34 = NA_real_,
        Prec_e24  = NA_real_,  Rec_e24 = NA_real_,
        Prec_e22  = NA_real_,  Rec_e22 = NA_real_,
        Accuracy  = NA_real_,
        MacroF1   = NA_real_,
        OuterFold = outer
      )
    }

    outer_results[[outer]] <- metrics
#    outer_rankedProteins[[outer]] <- rankedProteins
    # Save performance values of all folds completed to global environment, in case of early termination.
    outer_results.global<<-outer_results
    outer_hold_perf.global<<-outer_hold_perf
  } # close outer CV loop

  ## Pool the results from the outer fold(s) for the final report
  metrics <- dplyr::bind_rows(outer_results)
  hold_df <- dplyr::bind_rows(outer_hold_perf, .id=NULL)

#  rankedProteins[["sixGeno"]] <<- outer_rankedProteins

#  cat(sprintf("\n========== multi-fold outer CV summary ==========\n"))
#  cat(sprintf("Outer-fold Precision %.3f ± %.3f | Recall %.3f ± %.3f\n\n",
#              mean(hold_df$Precision, na.rm = TRUE), sd(hold_df$Precision, na.rm = TRUE),
#              mean(hold_df$Recall   , na.rm = TRUE), sd(hold_df$Recall   , na.rm = TRUE)))

  outer_macro <- hold_df |>
    dplyr::group_by(.data$OuterFold) |>
    dplyr::summarise(MacroPrec_e22 = mean(Prec_e22, na.rm = TRUE),
                     MacroRec_e22  = mean(Rec_e22,    na.rm = TRUE),
                     MacroPrec_e23 = mean(Prec_e23, na.rm = TRUE),
                     MacroRec_e23  = mean(Rec_e23,    na.rm = TRUE),
                     MacroPrec_e24 = mean(Prec_e24, na.rm = TRUE),
                     MacroRec_e24  = mean(Rec_e24,    na.rm = TRUE),
                     MacroPrec_e33 = mean(Prec_e33, na.rm = TRUE),
                     MacroRec_e33  = mean(Rec_e33,    na.rm = TRUE),
                     MacroPrec_e34 = mean(Prec_e34, na.rm = TRUE),
                     MacroRec_e34  = mean(Rec_e34,    na.rm = TRUE),
                     MacroPrec_e44 = mean(Prec_e44, na.rm = TRUE),
                     MacroRec_e44  = mean(Rec_e44,    na.rm = TRUE),
                     Accuracy      = mean(Accuracy, na.rm = TRUE),
                     MacroF1       = mean(MacroF1, na.rm = TRUE), .groups="drop")

  cat(sprintf(
    "\n\n========== multi-fold outer CV summary ==========\n"   ),
    sprintf("e2/e2 Precision  %.3f ± %.3f | Recall %.3f ± %.3f\n e2/e3 Precision  %.3f ± %.3f | Recall %.3f ± %.3f\n e2/e4 Precision  %.3f ± %.3f | Recall %.3f ± %.3f\n e3/e3 Precision  %.3f ± %.3f | Recall %.3f ± %.3f\n e3/e4 Precision  %.3f ± %.3f | Recall %.3f ± %.3f\n e4/e4 Precision  %.3f ± %.3f | Recall %.3f ± %.3f\n",
            mean(outer_macro$MacroPrec_e22, na.rm = TRUE),
            sd  (outer_macro$MacroPrec_e22, na.rm = TRUE),
            mean(outer_macro$MacroRec_e22 , na.rm = TRUE),
            sd  (outer_macro$MacroRec_e22 , na.rm = TRUE),
            mean(outer_macro$MacroPrec_e23, na.rm = TRUE),
            sd  (outer_macro$MacroPrec_e23, na.rm = TRUE),
            mean(outer_macro$MacroRec_e23 , na.rm = TRUE),
            sd  (outer_macro$MacroRec_e23 , na.rm = TRUE),
            mean(outer_macro$MacroPrec_e24, na.rm = TRUE),
            sd  (outer_macro$MacroPrec_e24, na.rm = TRUE),
            mean(outer_macro$MacroRec_e24 , na.rm = TRUE),
            sd  (outer_macro$MacroRec_e24 , na.rm = TRUE),
            mean(outer_macro$MacroPrec_e33, na.rm = TRUE),
            sd  (outer_macro$MacroPrec_e33, na.rm = TRUE),
            mean(outer_macro$MacroRec_e33 , na.rm = TRUE),
            sd  (outer_macro$MacroRec_e33 , na.rm = TRUE),
            mean(outer_macro$MacroPrec_e34, na.rm = TRUE),
            sd  (outer_macro$MacroPrec_e34, na.rm = TRUE),
            mean(outer_macro$MacroRec_e34 , na.rm = TRUE),
            sd  (outer_macro$MacroRec_e34 , na.rm = TRUE),
            mean(outer_macro$MacroPrec_e44, na.rm = TRUE),
            sd  (outer_macro$MacroPrec_e44, na.rm = TRUE),
            mean(outer_macro$MacroRec_e44 , na.rm = TRUE),
            sd  (outer_macro$MacroRec_e44 , na.rm = TRUE) ))
  cat(sprintf(
        "=================================================\n"   ),
      sprintf("Accuracy %.3f ± %.3f | Macro-F1 %.3f ± %.3f\n\n",
            mean(outer_macro$Accuracy, na.rm = TRUE),
            sd  (outer_macro$Accuracy, na.rm = TRUE),
            mean(outer_macro$MacroF1 , na.rm = TRUE),
            sd  (outer_macro$MacroF1 , na.rm = TRUE) ))

#  hold_df.byGeno <- hold_df |>
#    dplyr::group_by(Genotype) |>
#    dplyr::summarise(across(c(TP,FP,FN), sum, na.rm = TRUE))
#
#  cat("===== aggregated TP / FP / FN over all outer folds =====\n")
#  print(colSums(hold_df.byGeno[, c("TP","FP","FN")], na.rm = TRUE))
  cat("========================================================\n")

  ## ------------------------------------------------------------------------
  ## 3.  fit final ensemble on all data  (sequential) -----------------------
  ## ------------------------------------------------------------------------
  prep_expr <- function(mat, thr_miss = 0.20) {
    keep <- which(colMeans(is.na(mat)) <= thr_miss)
    mat  <- mat[, keep, drop = FALSE]
#    mat  <- VIM::kNN(mat, k = 5, imp_var = FALSE)
    scale(mat)
  }

  with_progress({
    p2 <- progressor(steps=3)
    p2(sprintf("Starting fit of final ensemble on all data (serial/non-parallel)..."), amount=0)

  levels_all <- sort(unique(APOE_gt))

  X_all <- prep_expr(expr)
  y_all <- factor(APOE_gt, levels=levels_all)

#  if (any(table(y_all) < 2)) message(sprintf("Final ensemble - skipped glmnet (rare class)"))
#  final_glm <- safe_glmnet(X_all, y_all, class_w[as.character(y_all)])
   final_glm <- cv.glmnet(X_all, y_all, family = "multinomial", weights = class_w[as.character(y_all)],
                          type.measure = "class", parallel = FALSE)
    p2(sprintf("GLM fit on full data  •  finished"), amount=1)

  xpar_all <- list(eta = 0.1, max_depth = 6, subsample = 0.8,
                   colsample_bytree = 0.8,
                   objective = "multi:softprob",
                   num_class = length(levels(y_all)),
                   nthread = ncores,          # can use all cores now
                   eval_metric = "mlogloss")       # explicitly set eval metric to avoid warning

  w_all    <- class_w[as.character(y_all)]
  d.all    <- xgb.DMatrix(X_all, label = as.numeric(y_all)-1, weight = w_all)
  bst_all  <- xgb.train(params=xpar_all, data=d.all,
                        watchlist=list(train=d.all),
                        nrounds = 200, verbose = 0,   #instead of 200, could use best #rounds: bst$best_iteration
                        early_stopping_rounds = 20)
    p2(sprintf("XGboost fit on full data  •  finished"), amount=1)

  rf_all <- ranger(y_all ~ ., x=data.frame(X_all, check.names=TRUE), y=y_all, class.weights = class_w,
                   probability = TRUE, num.trees = 500,
                   importance = "impurity",
                   num.threads = ncores, verbose=FALSE)
    p2(sprintf("Random Forest fit on full data  •  finished"), amount=1)

    p2(sprintf("Fit on full data  •  finished"), amount=0)

  ## Final models have been fit.


  ## ---------------- feature importance -------------------------------- ###
  importances.6geno.list<-coef(final_glm, s="lambda.min")
  importances.6geno.df<-as.data.frame(do.call(cbind, lapply(importances.6geno.list, function(x) abs(as.matrix(x)[-1,1]))))
  colnames(importances.6geno.df)<-names(coef(final_glm))
  importances.6geno.df$Sum<-apply(importances.6geno.df,1,sum)
  imp_glm <- importances.6geno.df  # abs(as.matrix(coef(final_glm, s="lambda.min"))[-1,1])
#final_glm.global<<-final_glm

  imp_xgb <- {
       g <- xgb.importance(model = bst_all)
       setNames(g$Gain, g$Feature)
  }
  imp_rf  <- rf_all$variable.importance
##  names(imp_rf)<-sub("\\.","|",names(imp_rf))
#imp_xgb.global<<-imp_xgb
#imp_rf.global<<-imp_rf
#imp_glm.global<<-imp_glm

  # Rename importance scores back to original feature names
  rownames(imp_glm)<-name_mapping[rownames(imp_glm)]
  names(imp_xgb) <-name_mapping[names(imp_xgb)]
  names(imp_rf) <- name_mapping[names(imp_rf)]

  # Ensure all features exist in each importance list (original column names from input, so check.names=FALSE)
  imp_glm[setdiff(name_mapping, rownames(imp_glm)), "Sum"] <- 0
  imp_xgb[setdiff(name_mapping, names(imp_xgb))] <- 0
  imp_rf[setdiff(name_mapping, names(imp_rf))] <- 0

  ## put everything on the same scale and average
  all_feats <- union(rownames(imp_glm), union(names(imp_xgb), names(imp_rf)))
  top_k=length(all_feats)

  imp_mat   <- cbind(
     glm = imp_glm[all_feats,"Sum"], xgb = imp_xgb[all_feats], rf = imp_rf[all_feats])
#imp_mat.global<<-imp_mat
  imp_mat[is.na(imp_mat)] <- 0
  imp_scaled <- scale(imp_mat)
  imp_mean   <- rowMeans(imp_scaled)
  top_feats  <- sort(imp_mean, decreasing = TRUE) #[1:top_k]

  ## -------- export top_k predictive proteins --------------------------
  if (!exists("rankedProteins", envir=.GlobalEnv)) { assign("rankedProteins",list(), envir=.GlobalEnv) } else { rankedProteins[["sixGeno"]] <- NULL }
  rankedProteins[["sixGeno"]] <<- data.frame(feature=names(sort(top_feats, decreasing = TRUE)), importance=sort(top_feats, decreasing = TRUE))

  names(top_feats)<-gsub("\\|","_",names(top_feats))
  if (Sys.getenv("RSTUDIO") == 1) {  # inside RStudio, avoid using print()
    cat(sprintf("Top 10 and bottom 10 (of all %d) predictive proteins - (all exported to list rankedProteins)\n", top_k))
    print(data.frame(Protein = names(top_feats)[c(1:10,(top_k-9):top_k)],
                     Importance = round(top_feats,3)[c(1:10,(top_k-9):top_k)], check.names=FALSE))
  } else {
  print(knitr::kable(data.frame(Protein = names(top_feats)[c(1:10,(top_k-9):top_k)],
                                Importance = round(top_feats,3)[c(1:10,(top_k-9):top_k)]),
                     caption = sprintf("Top 10 and bottom 10 (of all %d) predictive proteins - (all exported to list rankedProteins)", top_k)))
  }
  ### ------------------------------------------------------------------- ###

  end_time <- Sys.time()    # Record end time
  elapsed_time <- difftime(end_time, start_time, units="secs")  # Compute elapsed time
  if (elapsed_time<60) { message("Execution time: ", round(elapsed_time,1), " seconds") } else if (elapsed_time < 3600) { message("Execution time: ", round(elapsed_time / 60, 1), " minutes") } else { message("Execution time: ", round(elapsed_time / 3600, 1), " hours") }

  ## ------------------------------------------------------------------------
  ## 4.  prediction wrapper -------------------------------------------------
  ## ------------------------------------------------------------------------
  fn.out<-function(new_expr, confidenceRules=FALSE) {
#    colnames(new_expr)<-make.names(colnames(new_expr),unique=TRUE)
    new_expr <- prep_expr(rbind(expr[1,,drop=FALSE], new_expr))[-1,,drop=FALSE]
    new_expr <- new_expr[, colnames(X_all), drop = FALSE]

    ## helper - makes every probability matrix  n x 6
    complete_prob_mat <- function(mat, all_levels) {
      miss <- setdiff(all_levels, colnames(mat))
      if (length(miss)) {
        mat <- cbind(mat, matrix(0, nrow(mat), length(miss),
                                 dimnames = list(NULL, miss)))
      }
      mat[ , all_levels, drop = FALSE]
    }

    ## ---- individual learners -------------------------------------------------
    #p1 <- predict(final_glm, new_expr, s = "lambda.min",
    #              type = "response")[,,1]
    #p1 <- complete_prob_mat(p1, levels_all)
    if (!is.null(final_glm)) {
      p1 <- predict(final_glm, new_expr, s="lambda.min",
                    type="response")[,,1]
      p1 <- complete_prob_mat(p1, levels_all)
    }

    p2_mat <- matrix(predict(bst_all, new_expr),
                     ncol = length(levels_all), byrow = TRUE)
    colnames(p2_mat) <- levels_all
    p2_mat <- complete_prob_mat(p2_mat, levels_all)

    p3_mat <- predict(rf_all, data.frame(new_expr, check.names=TRUE))$predictions
    colnames(p3_mat) <- rf_all$forest$levels
    p3_mat <- complete_prob_mat(p3_mat, levels_all)

    if (!is.null(final_glm)) {
      p_list <- list(p1, p2_mat, p3_mat)
    } else {
      p_list <- list(p2_mat, p3_mat)          # skip glmnet
    }

    ## ---- averaged probabilities & confidence filter --------------------------
    #p_avg <- (p1 + p2_mat + p3_mat) / 3      # n x 6
    p_avg <- Reduce(`+`, p_list) / length(p_list)

    top_class <- apply(p_avg, 1, function(z) names(z)[which.max(z)])
    top_prob  <- apply(p_avg, 1, max)
    second    <- apply(p_avg, 1, function(z) sort(z, TRUE)[2])

    if(confidenceRules) {
      keep_homo  <- top_class %in% c("e3/e3", "e4/e4") & top_prob >= 0.80
      keep_other <- !(top_class %in% c("e3/e3", "e4/e4")) &
                    top_prob >= 0.90 & (top_prob - second) >= 0.20

      pred <- ifelse(keep_homo | keep_other, top_class, NA_character_)
    } else {
      pred <- top_class
    }
    factor(pred, levels = levels_all)
    }
  return(fn.out)
  }) # with_progress
}

#option 3 - not run
#predict_APOE.b3.allAndMapped.6geno.nestedCV.smallTest<-fit_APOE_ensemble_6way_nestedCV(t(training.cleanDat[APOE4.assays.keep[c(1:250,2251:length(APOE4.assays.keep))], keep_idx]), na.omit(gt.APOE)[keep_idx], ncores=8, outerFolds=2)

# data with no missing values to train on (na.omit by sample), and IRF6 (single row missing in a cohort) removed first
dim(training.cleanDat.noNA)
# 7334 14758
length(training.gt.APOE)
# 14758
22392-5715
# out of 16677 possible


# On VM Windows02:  2235 features with Positive Importance in any of the binary genotype models (allImportantFeatures)
start_time <- Sys.time()  # Record start time
#  predict_APOE.b3.allAndMapped.6geno.nestedCV.smallTest<-fit_APOE_ensemble_6way_nestedCV(t(training.cleanDat[APOE4.assays.keep[c(1:250,2251:length(APOE4.assays.keep))], keep_idx]), na.omit(gt.APOE)[keep_idx], ncores=8, nfold=3, nrep=3, outerFolds=2)
#  501 features run on full set of 400x5+33 e2/e2 - for quick test run
  predict_APOE.b3.allAndMapped.6geno.nestedCV<-fit_APOE_ensemble_6way_nestedCV(t(training.cleanDat.noNA[allImportantFeatures,]), training.gt.APOE, outerFolds=3, ncores=8)  #20 cores ok on VM 03
#  #2235 features with importance scaled mean > 0 in any of 6 binary prediction runs
end_time <- Sys.time()    # Record end time
elapsed_time <- difftime(end_time, start_time, units="secs")  # Compute elapsed time
if (elapsed_time<60) { message("Execution time: ", round(elapsed_time,1), " seconds") } else if (elapsed_time < 3600) { message("Execution time: ", round(elapsed_time / 60, 1), " minutes") } else { message("Execution time: ", round(elapsed_time / 3600, 1), " hours") }


# On VM Windows03: 7334 (all but 1) features
  predict_APOE.b3.allAndMapped.6geno.nestedCV<-fit_APOE_ensemble_6way_nestedCV(t(training.cleanDat.noNA), training.gt.APOE, outerFolds=3, ncores=14, memLimitGB=Inf)  #20 cores ok on VM 03 for 2235 features.
  # all 7334 features used for training.  9.6 hours run time
# Impute APOE 6 genotypes using the weighted 3ML wrapper in the 5715 (5k) samples missing it (step 2 of 2)
# fails due to naming being nonstandard... cannot use output wrapper.
#imputed.APOE.b3.allAndMapped.6geno.nestedCV.uneqWeight<-predict_APOE.b3.allAndMapped.6geno.nestedCV(t(impute::impute.knn(cleanDat.4p13b3)$data[rownames(training.cleanDat.noNA),which(is.na(gt.APOE))]))
#table(imputed.APOE.b3.allAndMapped.6geno.nestedCV.uneqWeight)

# Used to troubleshoot
predict_APOE.b3.allAndMapped.6geno.nestedCV.smallTest<-fit_APOE_ensemble_6way_nestedCV(t(training.cleanDat[APOE4.assays.keep[c(1:250,2251:length(APOE4.assays.keep))], keep_idx]), na.omit(gt.APOE)[keep_idx], ncores=28, nfold=3, nrep=3, outerFolds=2)
test<-predict_APOE.b3.allAndMapped.6geno.nestedCV.smallTest(t(impute::impute.knn(cleanDat.4p13b3)$data[APOE4.assays.keep[c(1:250,2251:length(APOE4.assays.keep))], which(is.na(gt.APOE))]))
table(test) #e4/e4 overpredicted by >10x

# 2235 features with importance scaled mean > 0 in any of 6 binary prediction runs
predict_APOE.b3.allAndMapped.6geno.nestedCV.2235feat<-fit_APOE_ensemble_6way_nestedCV(t(training.cleanDat.noNA[allImportantFeatures,]), training.gt.APOE, outerFolds=0, ncores=30, memLimitGB=Inf)  #30 cores ok on VM 03 for 2235 features.
#2.6 hours
imputed.APOE.b3.allAndMapped.6geno.nestedCV.eqWt.2235feat<-predict_APOE.b3.allAndMapped.6geno.nestedCV.2235feat( t(impute::impute.knn(cleanDat.4p13b3)$data[allImportantFeatures, which(is.na(gt.APOE))]) )
table(imputed.APOE.b3.allAndMapped.6geno.nestedCV.eqWt.2235feat)
#e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#   13   579   103  3247  1603   170
top1000features<-rankedProteins$sixGeno[1:1000,]

# Top 1000 ranked features from above used as input in next round;
predict_APOE.b3.allAndMapped.6geno.nestedCV.1000feat<-fit_APOE_ensemble_6way_nestedCV(t(training.cleanDat.noNA[top1000features$feature,]), training.gt.APOE, outerFolds=0, ncores=30, memLimitGB=Inf)  #30 cores ok on VM 03 for 1000 features.
#2.5 hours
imputed.APOE.b3.allAndMapped.6geno.nestedCV.eqWt.1000feat<-predict_APOE.b3.allAndMapped.6geno.nestedCV.1000feat( t(impute::impute.knn(cleanDat.4p13b3)$data[top1000features$feature, which(is.na(gt.APOE))]) )
table(imputed.APOE.b3.allAndMapped.6geno.nestedCV.eqWt.1000feat)
#e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#   16   575   107  3192  1658   167
top500features<-rankedProteins$sixGeno[1:500,]

# Top 500 ranked features from above used as input in next round;
predict_APOE.b3.allAndMapped.6geno.nestedCV.500feat<-fit_APOE_ensemble_6way_nestedCV(t(training.cleanDat.noNA[top500features$feature,]), training.gt.APOE, outerFolds=0, ncores=30, memLimitGB=Inf)  #30 cores ok on VM 03 for 500 features.
#1.0 hours
imputed.APOE.b3.allAndMapped.6geno.nestedCV.eqWt.500feat<-predict_APOE.b3.allAndMapped.6geno.nestedCV.500feat( t(impute::impute.knn(cleanDat.4p13b3)$data[top500features$feature, which(is.na(gt.APOE))]) )
table(imputed.APOE.b3.allAndMapped.6geno.nestedCV.eqWt.500feat)
#e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#   15   582   104  3181  1659   174
top100features<-rankedProteins$sixGeno[1:100,]

# Top 100 ranked features from above used as input in next round;
predict_APOE.b3.allAndMapped.6geno.nestedCV.100feat<-fit_APOE_ensemble_6way_nestedCV(t(training.cleanDat.noNA[top100features$feature,]), training.gt.APOE, outerFolds=0, ncores=30, memLimitGB=Inf)  #30 cores ok on VM 03 for 500 features.
#8.4 minutes
imputed.APOE.b3.allAndMapped.6geno.nestedCV.eqWt.100feat<-predict_APOE.b3.allAndMapped.6geno.nestedCV.100feat( t(impute::impute.knn(cleanDat.4p13b3)$data[top100features$feature, which(is.na(gt.APOE))]) )
table(imputed.APOE.b3.allAndMapped.6geno.nestedCV.eqWt.100feat)
#e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#   16   579   113  3190  1653   164

# Top 25 ranked features from above used as input in next round; Outer Folds=1000 for statistics collection on hold out data
predict_APOE.b3.allAndMapped.6geno.nestedCV.25feat.1000outerFold<-fit_APOE_ensemble_6way_nestedCV(t(training.cleanDat.noNA[top100features$feature[1:25],]), training.gt.APOE, outerFolds=1000, ncores=30, memLimitGB=Inf)  #30 cores ok on VM 03 for 25 features.
#6.9 minutes for 2 outerFolds

#outer_hold_perf.global  # outer stats, if we want to summarize from a run that does not complete 1000 outer folds, the completed fold stats are here.
#outer_results.global    # inner stats

imputed.APOE.b3.allAndMapped.6geno.nestedCV.eqWt.25feat.1000outerFold<-predict_APOE.b3.allAndMapped.6geno.nestedCV.25feat.1000outerFold( t(impute::impute.knn(cleanDat.4p13b3)$data[top100features$feature[1:25], which(is.na(gt.APOE))]) )
table(imputed.APOE.b3.allAndMapped.6geno.nestedCV.eqWt.25feat.1000outerFold)


#outer_hold_perf.global.25eachRare
#outer_results.global.25eachRare

# Top 25 features missingness inthe 5715
data.frame(Feature=top100features$feature[1:25], Count.missed=apply(cleanDat.4p13b3[top100features$feature[1:25], which(is.na(gt.APOE))],1,function(x) length(which(is.na(x)))))

#outer_hold_perf.global
 hold_df <- dplyr::bind_rows(outer_hold_perf.global, .id=NULL)

  outer_macro <- hold_df |>
    dplyr::group_by(.data$OuterFold) |>
    dplyr::summarise(MacroPrec_e22 = mean(Prec_e22, na.rm = TRUE),
                     MacroRec_e22  = mean(Rec_e22,    na.rm = TRUE),
                     MacroPrec_e23 = mean(Prec_e23, na.rm = TRUE),
                     MacroRec_e23  = mean(Rec_e23,    na.rm = TRUE),
                     MacroPrec_e24 = mean(Prec_e24, na.rm = TRUE),
                     MacroRec_e24  = mean(Rec_e24,    na.rm = TRUE),
                     MacroPrec_e33 = mean(Prec_e33, na.rm = TRUE),
                     MacroRec_e33  = mean(Rec_e33,    na.rm = TRUE),
                     MacroPrec_e34 = mean(Prec_e34, na.rm = TRUE),
                     MacroRec_e34  = mean(Rec_e34,    na.rm = TRUE),
                     MacroPrec_e44 = mean(Prec_e44, na.rm = TRUE),
                     MacroRec_e44  = mean(Rec_e44,    na.rm = TRUE),
                     Accuracy      = mean(Accuracy, na.rm = TRUE),
                     MacroF1       = mean(MacroF1, na.rm = TRUE), .groups="drop")

  cat(sprintf(
    "\n\n========== multi-fold outer CV summary ==========\n"   ),
    sprintf("e2/e2 Precision  %.3f ± %.3f | Recall %.3f ± %.3f\n e2/e3 Precision  %.3f ± %.3f | Recall %.3f ± %.3f\n e2/e4 Precision  %.3f ± %.3f | Recall %.3f ± %.3f\n e3/e3 Precision  %.3f ± %.3f | Recall %.3f ± %.3f\n e3/e4 Precision  %.3f ± %.3f | Recall %.3f ± %.3f\n e4/e4 Precision  %.3f ± %.3f | Recall %.3f ± %.3f\n",
            mean(outer_macro$MacroPrec_e22, na.rm = TRUE),
            sd  (outer_macro$MacroPrec_e22, na.rm = TRUE),
            mean(outer_macro$MacroRec_e22 , na.rm = TRUE),
            sd  (outer_macro$MacroRec_e22 , na.rm = TRUE),
            mean(outer_macro$MacroPrec_e23, na.rm = TRUE),
            sd  (outer_macro$MacroPrec_e23, na.rm = TRUE),
            mean(outer_macro$MacroRec_e23 , na.rm = TRUE),
            sd  (outer_macro$MacroRec_e23 , na.rm = TRUE),
            mean(outer_macro$MacroPrec_e24, na.rm = TRUE),
            sd  (outer_macro$MacroPrec_e24, na.rm = TRUE),
            mean(outer_macro$MacroRec_e24 , na.rm = TRUE),
            sd  (outer_macro$MacroRec_e24 , na.rm = TRUE),
            mean(outer_macro$MacroPrec_e33, na.rm = TRUE),
            sd  (outer_macro$MacroPrec_e33, na.rm = TRUE),
            mean(outer_macro$MacroRec_e33 , na.rm = TRUE),
            sd  (outer_macro$MacroRec_e33 , na.rm = TRUE),
            mean(outer_macro$MacroPrec_e34, na.rm = TRUE),
            sd  (outer_macro$MacroPrec_e34, na.rm = TRUE),
            mean(outer_macro$MacroRec_e34 , na.rm = TRUE),
            sd  (outer_macro$MacroRec_e34 , na.rm = TRUE),
            mean(outer_macro$MacroPrec_e44, na.rm = TRUE),
            sd  (outer_macro$MacroPrec_e44, na.rm = TRUE),
            mean(outer_macro$MacroRec_e44 , na.rm = TRUE),
            sd  (outer_macro$MacroRec_e44 , na.rm = TRUE) ))
  cat(sprintf(
        "=================================================\n"   ),
      sprintf("Accuracy %.3f ± %.3f | Macro-F1 %.3f ± %.3f\n\n",
            mean(outer_macro$Accuracy, na.rm = TRUE),
            sd  (outer_macro$Accuracy, na.rm = TRUE),
            mean(outer_macro$MacroF1 , na.rm = TRUE),
            sd  (outer_macro$MacroF1 , na.rm = TRUE) ))

rm(hold_df)
rm(outer_macro)


# 6 genotype ensemble classifier prediction frequencies in 5715 unknowns:

table(imputed.APOE.b3.allAndMapped.noNA)  # 7334 features input
#e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#   12   578   101  3229  1625   170
#...6 binary classifier runs were used to get all features with any absolute importance to any genotype prediction (n=2235)
table(imputed.APOE.b3.allAndMapped.6geno.nestedCV.eqWt.2235feat)
#e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#   13   579   103  3247  1603   170
top1000features<-rankedProteins$sixGeno[1:1000,]

table(imputed.APOE.b3.allAndMapped.6geno.nestedCV.eqWt.1000feat)
#e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#   16   575   107  3192  1658   167
top500features<-rankedProteins$sixGeno[1:500,]

table(imputed.APOE.b3.allAndMapped.6geno.nestedCV.eqWt.500feat)
#e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#   15   582   104  3181  1659   174
top100features<-rankedProteins$sixGeno[1:100,]

table(imputed.APOE.b3.allAndMapped.6geno.nestedCV.eqWt.100feat)
#e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#   16   579   113  3190  1653   164
#top 25 of 100 features directly used without reranking recursively.

table(imputed.APOE.b3.allAndMapped.6geno.nestedCV.eqWt.25feat.0outerFold)
#e2/e2 e2/e3 e2/e4 e3/e3 e3/e4 e4/e4
#   18   575   112  3302  1548   160


  ## Run the regression (4p13 b3) - Site with Age+Sex protection
  normExpr.reg <- matrix(NA,nrow=nrow(cleanDat.unreg),ncol=ncol(cleanDat.unreg))
  #coefmat <- matrix(NA,nrow=nrow(cleanDat.unreg),ncol=9)  #ncol(regvars)+1) ## change this to ncol(regvars)+2 when condition has 2 levels if BOOT=TRUE, +1 if BOOT=FALSE

  #RNG seed set for reproducibility
  set.seed(1234567);
  normExpr.reg <-  foreach (i=1:nrow(cleanDat.unreg), .combine=rbind, .packages="stats") %dopar% {  #** .combine=c, .multicombine=TRUE, .packages="stats", .export="regvars") %dopar% {
        set.seed(1234567)
        options(stringsAsFactors=FALSE)
        tryCatch({
          lmmod1 <- lm(as.numeric(cleanDat.unreg[i,])~Age+Sex +Site,data=regvars)
          #** list(list(coef=coef(lmmod1),residuals=lmmod1$residuals))  #return a list of (1) vector of coefficients; and (2) vector of residuals length nrow(cleanDat.unreg)
          ##datpred <- predict(object=lmmod1,newdata=regvars)
          coef <- coef(lmmod1)
#          coef[1] + coef[2]*regvars[,"HNRNPA2B1.TimeToSpin"] + coef[3]*regvars[,"HBZ.TimeToDecant"] + lmmod1$residuals ## The full data + with undesired covariates
          coef[1] + coef[2]*regvars[,"Age"] + coef[3]*abs(Sex-1) + lmmod1$residuals ## The full data - the undesired covariates left out  (sex=0 Female mean effect is represented by the coefficient, since males were reference level)
        }, error=function(e) { rep(NA_real_,ncol(cleanDat.unreg)) })
  }
  rownames(normExpr.reg) <- rownames(cleanDat.unreg)
  colnames(normExpr.reg) <- colnames(cleanDat.unreg)

  cat(paste0("Finished Pass 4p13 b2 regression of intersite variance."\n"))


saveRDS(normExpr.reg,"4p13b2.normExpr.reg_sites1-19_Fsplit.RDS")

cleanDat.4p13b2<-normExpr.reg
#numericMeta.reg still valid, ordered for this cleanDat


#.#.#


cleanDat.backup<-cleanDat

######################################
## ANOVA + Volcanoes + DEXstacked Barplots
source("parANOVA.dex.fallback7.25.R")

parallelThreads=7
outFilePrefix="4p12b.regr.6SV+medianStrip"
outFileSuffix="SixDxVolcs.CTimputed"

cleanDat<-cleanDat.zeroSiteMedian

#Grouping=numericMeta$Group #***
Grouping=numericMeta$Group.withCTimputed
Grouping[which(numericMeta$Group.withCTimputed=="CI.Other")]<-NA
Grouping[which(numericMeta$Group.withCTimputed=="AsymAD")]<-NA
#Grouping[which(numericMeta$Group.withCTimputed=="MCI")]<-NA

ANOVAout <- parANOVA.dex()

#*** commented new for 3-genotype only comparisons
##head(ANOVAout)
## plot volcanoes only for comparisons with 2/3 Group components in common.
#comps<-cbind(ANOVAout.colnames=colnames(ANOVAout), as.data.frame(do.call(rbind,strsplit(colnames(ANOVAout),"[.|-]"))))
#lastComp=which(grepl("diff ",comps$ANOVAout.colnames))[1]-1
#lastComp
#155  # 155-2 = 153 comparisons pairwise
## keep comparisons with 2/3 variables the same (n=45 comparisons)
#selectComps=unlist( sapply(3:lastComp,function(x) if(length(which(duplicated(unlist(comps[x,2:7]))))>1) x ) )
#length(selectComps) # selectComps is used directly by plotVolc to choose which comparisons to plot (column #s of p values in ANOVAout)
## Decide which comparison X axis (logFC) to flip
#comps[selectComps,1] # view pairwise comparisons selected for volcano plotting
# None need flipping (more extreme expected phenotype is first (numerator) in each case)
flip=c(3,4,5,6,7,8)
sameScale=FALSE #TRUE
symbolsOnly=TRUE
highlightGeneProducts=c("HNRNPA2B1","HBZ","SPC25","CPLX2","PTN","MAPT","OMG","GDF15","NEFL","NRGN","CTHRC1","HTRA1","TTR","LRP1","NTN1","SFRP1","GPNMB","PAFAH1B3")
labelHighlighted=TRUE      # if true, highlighted spots get text labels with their rownames from ANOVAout
labelTop=5
plotVolc()                 # runs on ANOVAout as input (need not be specified).

#DEXpercentStacked()        # runs on prior function outputs as input; writes stacked bar plot(s) to PDF.


numericMeta.plasma<-numericMeta #[which(numericMeta$sample_matrix=="CSF"),]
dim(numericMeta.plasma)
library(ggplot2)
library(ggrepel)

## tSNE plot for Each site alone - colored by Age (19 pages)
tSNE.plasma.samples.age.1site<-Age<-list()
for(site in names(table(numericMeta.plasma$contributor_code))) {

  Age[[site]]=numericMeta.plasma$age_at_visit[which(numericMeta.plasma$contributor_code==site)]

  tSNE.plasma.samples.age.1site[[site]]<-ggplot2::ggplot(tSNE.plasma.xy[which(numericMeta.plasma$contributor_code==site),],label=1:length(which(numericMeta.plasma$contributor_code==site))) + geom_point(aes(x=x,y=y, color=Age[[site]]), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
    geom_text_repel(data = tSNE.plasma.xy[which(numericMeta.plasma$contributor_code==site)[1], ],
                    aes(x=x,y=y, label = site), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
    theme_minimal() +  # Minimal theme
    theme(
      panel.background = element_blank(),  # Remove plot area color
  #    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
      legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
      legend.key = element_rect(fill = "white"),  # Keep legend keys clean
      axis.title.x = element_text(size = 28),  # Double x-axis label text size
      axis.title.y = element_text(size = 28)
    )
}


pdf(file="4p12b.regr.6SVregr+medianStrip-tSNE-Plasma(7335x22392)_protectAgeSex-separate19sites-samples_coloredByTraits.pdf",width=11,height=9)
  for(site in names(table(numericMeta.plasma$contributor_code))) print(tSNE.plasma.samples.age.1site[[site]])
dev.off()


## Redo-color by Dx group (with CT imputed and 3 Emory study sites updated Dx group.pathCog
tSNE.plasma.samples.age.1site<-GroupDx<-list()
for(site in names(table(numericMeta.plasma$contributor_code))) {

  #Age[[site]]=numericMeta.plasma$age_at_visit[which(numericMeta.plasma$contributor_code==site)]
  GroupDx[[site]]=as.factor(numericMeta$Group.withCTimputed[which(numericMeta.plasma$contributor_code==site)])

  tSNE.plasma.samples.age.1site[[site]]<-ggplot2::ggplot(tSNE.plasma.xy[which(numericMeta.plasma$contributor_code==site),],label=1:length(which(numericMeta.plasma$contributor_code==site))) + geom_point(aes(x=x,y=y, color=GroupDx[[site]]), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
    geom_text_repel(data = tSNE.plasma.xy[which(numericMeta.plasma$contributor_code==site)[1], ],
                    aes(x=x,y=y, label = site), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
    theme_minimal() +  # Minimal theme
    theme(
      panel.background = element_blank(),  # Remove plot area color
  #    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
      legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
      legend.key = element_rect(fill = "white"),  # Keep legend keys clean
      axis.title.x = element_text(size = 28),  # Double x-axis label text size
      axis.title.y = element_text(size = 28)
    )
}

pdf(file="4p12b.regr.6SVregr+medianStrip-tSNE-Plasma(7335x22392)_protectAgeSex-separate19sites-samples_coloredByImputedDxGroup.pdf",width=11,height=9)
  for(site in names(table(numericMeta.plasma$contributor_code))) print(tSNE.plasma.samples.age.1site[[site]])
dev.off()


######################################
## ANOVA + Volcanoes - without the Site F (1600 CT after imputation)
source("parANOVA.dex.fallback7.25.R")

parallelThreads=7
outFilePrefix="4p12b.regr.6SV+medianStrip[noSiteF]"
outFileSuffix="SixDxVolcs.CTimputed[1600siteF_CT_out]"

cleanDat<-cleanDat.zeroSiteMedian[,which(numericMeta$contributor_code!="F")]  #(already have cleanDat.backup)

numericMeta.noF<-numericMeta[which(numericMeta$contributor_code!="F"),]

#Grouping=numericMeta$Group #***
Grouping=numericMeta.noF$Group.withCTimputed
Grouping[which(numericMeta.noF$Group.withCTimputed=="CI.Other")]<-NA
Grouping[which(numericMeta.noF$Group.withCTimputed=="AsymAD")]<-NA
#Grouping[which(numericMeta$Group.withCTimputed=="MCI")]<-NA

ANOVAout <- parANOVA.dex()

#*** commented new for 3-genotype only comparisons
##head(ANOVAout)
## plot volcanoes only for comparisons with 2/3 Group components in common.
#comps<-cbind(ANOVAout.colnames=colnames(ANOVAout), as.data.frame(do.call(rbind,strsplit(colnames(ANOVAout),"[.|-]"))))
#lastComp=which(grepl("diff ",comps$ANOVAout.colnames))[1]-1
#lastComp
#155  # 155-2 = 153 comparisons pairwise
## keep comparisons with 2/3 variables the same (n=45 comparisons)
#selectComps=unlist( sapply(3:lastComp,function(x) if(length(which(duplicated(unlist(comps[x,2:7]))))>1) x ) )
#length(selectComps) # selectComps is used directly by plotVolc to choose which comparisons to plot (column #s of p values in ANOVAout)
## Decide which comparison X axis (logFC) to flip
#comps[selectComps,1] # view pairwise comparisons selected for volcano plotting
# None need flipping (more extreme expected phenotype is first (numerator) in each case)
flip=c(3,4,5,6,7,8)
sameScale=FALSE #TRUE
symbolsOnly=TRUE
highlightGeneProducts=c("HNRNPA2B1","HBZ","SPC25","CPLX2","PTN","MAPT","OMG","GDF15","NEFL","NRGN","CTHRC1","HTRA1","TTR","LRP1","NTN1","SFRP1","GPNMB","PAFAH1B3")
labelHighlighted=TRUE      # if true, highlighted spots get text labels with their rownames from ANOVAout
labelTop=5
plotVolc()                 # runs on ANOVAout as input (need not be specified).

#DEXpercentStacked()        # runs on prior function outputs as input; writes stacked bar plot(s) to PDF.


######################################
## ANOVA + Volcanoes - without the Site F (1600 CT after imputation) - no Median stripping
source("parANOVA.dex.fallback7.25.R")

parallelThreads=7
outFilePrefix="4p12b.regr.6SV_NOmedianStrip[noSiteF]"
outFileSuffix="SixDxVolcs.CTimputed[1600siteF_CT_out]"

cleanDat<-cleanDat.backup[,which(numericMeta$contributor_code!="F")]  #(already have cleanDat.backup)

numericMeta.noF<-numericMeta[which(numericMeta$contributor_code!="F"),]

#Grouping=numericMeta$Group #***
Grouping=numericMeta.noF$Group.withCTimputed
Grouping[which(numericMeta.noF$Group.withCTimputed=="CI.Other")]<-NA
Grouping[which(numericMeta.noF$Group.withCTimputed=="AsymAD")]<-NA
#Grouping[which(numericMeta$Group.withCTimputed=="MCI")]<-NA

ANOVAout <- parANOVA.dex()

#*** commented new for 3-genotype only comparisons
##head(ANOVAout)
## plot volcanoes only for comparisons with 2/3 Group components in common.
#comps<-cbind(ANOVAout.colnames=colnames(ANOVAout), as.data.frame(do.call(rbind,strsplit(colnames(ANOVAout),"[.|-]"))))
#lastComp=which(grepl("diff ",comps$ANOVAout.colnames))[1]-1
#lastComp
#155  # 155-2 = 153 comparisons pairwise
## keep comparisons with 2/3 variables the same (n=45 comparisons)
#selectComps=unlist( sapply(3:lastComp,function(x) if(length(which(duplicated(unlist(comps[x,2:7]))))>1) x ) )
#length(selectComps) # selectComps is used directly by plotVolc to choose which comparisons to plot (column #s of p values in ANOVAout)
## Decide which comparison X axis (logFC) to flip
#comps[selectComps,1] # view pairwise comparisons selected for volcano plotting
# None need flipping (more extreme expected phenotype is first (numerator) in each case)
flip=c(3,4,5,6,7,8)
sameScale=FALSE #TRUE
symbolsOnly=TRUE
highlightGeneProducts=c("HNRNPA2B1","HBZ","SPC25","CPLX2","PTN","MAPT","OMG","GDF15","NEFL","NRGN","CTHRC1","HTRA1","TTR","LRP1","NTN1","SFRP1","GPNMB","PAFAH1B3")
labelHighlighted=TRUE      # if true, highlighted spots get text labels with their rownames from ANOVAout
labelTop=5
plotVolc()                 # runs on ANOVAout as input (need not be specified).

#DEXpercentStacked()        # runs on prior function outputs as input; writes stacked bar plot(s) to PDF.


## Visualize 4p12 intrasite median-zeroed tSNE site B (calculate tSNE for all points in all sites

## tSNE plot for (site B) alone - colored by Age - print to console
tSNE.plasma.samples.age.1site<-Age<-list()
for(site in c("B")) { #names(table(numericMeta.plasma$contributor_code))) {

  Age[[site]]=numericMeta.plasma$age_at_visit[which(numericMeta.plasma$contributor_code==site)]

  tSNE.plasma.samples.age.1site[[site]]<-ggplot2::ggplot(tSNE.plasma.xy[which(numericMeta.plasma$contributor_code==site),],label=1:length(which(numericMeta.plasma$contributor_code==site))) + geom_point(aes(x=x,y=y, color=Age[[site]]), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
    geom_text_repel(data = tSNE.plasma.xy[which(numericMeta.plasma$contributor_code==site)[1], ],
                    aes(x=x,y=y, label = site), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
    theme_minimal() +  # Minimal theme
    theme(
      panel.background = element_blank(),  # Remove plot area color
  #    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
      legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
      legend.key = element_rect(fill = "white"),  # Keep legend keys clean
      axis.title.x = element_text(size = 28),  # Double x-axis label text size
      axis.title.y = element_text(size = 28)
    )
}
print(tSNE.plasma.samples.age.1site[[site]])


## Remove tight outlier cluster of samples represented by points in the tSNE plot for isolated samples from site B

## Equation of the line through (20,20) and (30,10):
##   slope  m  = (10 - 20)/(30 - 20) = -1
##   so     y = -x + 40
## Indexes of points lying **above** that line
#B.outCluster.idx <- intersect( which(tSNE.plasma.xy$y > ((-1) * tSNE.plasma.xy$x + 40)),
#                               which(numericMeta$contributor_code=="B") )


## line through (-20,-20) and (-30,-10):
#m <- (-10 - (-20)) / (-30 - (-20))   # slope  = -1
#b <- -20 - m * (-20)                 # intercept = -40   (so y = -x - 40)

## indexes of points **below** that line
B.outCluster.idx <- intersect( which(tSNE.plasma.xy$y < (-1) * tSNE.plasma.xy$x - 40),
                               which(numericMeta$contributor_code=="B") )
length(B.outCluster.idx)
#417 samples are removed from B out of
length(which(numericMeta$contributor_code=="B"))
#1228


## Repeat tSNE on data unregressed for the 6SVs used to get (4p12) data
numericMeta.plasma<-numericMeta #[which(numericMeta$sample_matrix=="CSF"),]
dim(numericMeta.plasma)
#22392   91   #previously (16 sites with Calibrators):  18739   93
exprMat.plasma<-cleanDat.unreg #cleanDat.zeroSiteMedian  #exprMat0[,match(rownames(numericMeta.plasma),colnames(exprMat0))]


unreg.tSNE.list.plasma <- Rtsne::Rtsne(t(as.data.frame(na.omit(exprMat.plasma))),perplexity=20)
dim(unreg.tSNE.list.plasma$Y)
unreg.tSNE.plasma.xy<-as.data.frame(unreg.tSNE.list.plasma$Y)
colnames(unreg.tSNE.plasma.xy)<-c('x','y')


## unreg tSNE plot for (site B) alone - colored by outlier sample cluster - print to console
unreg.tSNE.plasma.samples.age.1site<-OLcluster<-list()
for(site in c("B")) { #names(table(numericMeta.plasma$contributor_code))) {

  OLcluster[[site]]=rep("spreadOut",nrow(unreg.tSNE.plasma.xy))
  OLcluster[[site]][B.outCluster.idx]="4p12.tightCluster"
  OLcluster[[site]]<-factor(OLcluster[[site]][which(numericMeta.plasma$contributor_code==site)])

  unreg.tSNE.plasma.samples.age.1site[[site]]<-ggplot2::ggplot(unreg.tSNE.plasma.xy[which(numericMeta.plasma$contributor_code==site),],label=1:length(which(numericMeta.plasma$contributor_code==site))) + geom_point(aes(x=x,y=y, color=OLcluster[[site]]), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
    geom_text_repel(data = unreg.tSNE.plasma.xy[which(numericMeta.plasma$contributor_code==site)[1], ],
                    aes(x=x,y=y, label = site), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
    theme_minimal() +  # Minimal theme
    theme(
      panel.background = element_blank(),  # Remove plot area color
  #    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
      legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
      legend.key = element_rect(fill = "white"),  # Keep legend keys clean
      axis.title.x = element_text(size = 28),  # Double x-axis label text size
      axis.title.y = element_text(size = 28)
    )
}
print(unreg.tSNE.plasma.samples.age.1site[[site]])


## unreg tSNE plot for (site F) alone - colored by Age - print to console
unreg.tSNE.plasma.samples.age.1site<-Age<-list()
for(site in c("F")) { #names(table(numericMeta.plasma$contributor_code))) {

  Age[[site]]=numericMeta.plasma$age_at_visit[which(numericMeta.plasma$contributor_code==site)]

  unreg.tSNE.plasma.samples.age.1site[[site]]<-ggplot2::ggplot(unreg.tSNE.plasma.xy[which(numericMeta.plasma$contributor_code==site),],label=1:length(which(numericMeta.plasma$contributor_code==site))) + geom_point(aes(x=x,y=y, color=Age[[site]]), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE DImension 2") +  # Axis labels
    geom_text_repel(data = unreg.tSNE.plasma.xy[which(numericMeta.plasma$contributor_code==site)[1], ],
                    aes(x=x,y=y, label = site), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
    theme_minimal() +  # Minimal theme
    theme(
      panel.background = element_blank(),  # Remove plot area color
  #    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
      legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
      legend.key = element_rect(fill = "white"),  # Keep legend keys clean
      axis.title.x = element_text(size = 28),  # Double x-axis label text size
      axis.title.y = element_text(size = 28)
    )
}
print(unreg.tSNE.plasma.samples.age.1site[[site]])


#<STOP 4p12 here>
# NOT RUN BELOW HERE 04/25/2025
#
# ANNOTATION: Removed 1801 archival lines that were explicitly marked as not run below here in the source script.
# The original upload is unchanged; this cleaned manuscript-facing copy stops at the recorded execution boundary.
