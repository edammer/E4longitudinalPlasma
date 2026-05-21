##############################################################################
# Pipeline annotation header: 3+4.(2ndRegr.2)4p13c1.PLASMA_known+predictedAPOE_e4carrierRegressedSite(secondRegr)-WGCNA.R
# Manuscript code section(s): 3 / 4
#
# Purpose:
# Use the finalized APOE-aware regression matrix for WGCNA network
# construction, module membership refinement, module-trait visualization,
# differential-expression volcanoes, and enrichment/network summaries.
#
# Principal inputs:
#   - saved.image-genotype_prediction_finalized.RData
#   - geneListFET_customLabels-fixedScale+15thPlotNominalP.R
#   - parANOVA.dex.fallback7.25.R
#   - GOparallel-FET.R
#   - buildIgraphs.R
#
# Principal outputs:
#   - 4p13c1.regr.2PAV+Site(protAgeSexE4).PLASMA_sft.22392.RDS
#   - 4p13c1.GlobalNetworkPlots*.pdf
#   - 4p13c1.ModuleAssignments*.txt
#   - 22392_sample_net_12MEs.RDS
#   - saved.image-genotype_prediction_finalized.RData
#
# Step overview:
#   1. Load the finalized genotype prediction/regression image and confirm
#      residual APOE-e4 proxy associations.
#   2. Optionally perform outlier checks while preserving the final
#      sample/protein matrix.
#   3. Select a signed bicor WGCNA power and construct modules with
#      blockwiseModules.
#   4. Remove the two PAV regression assays, compute module eigengenes, signed
#      kME, and refined module membership assignments.
#   5. Generate global network plots, enrichment summaries, ANOVA/volcano
#      outputs, and module graph visualizations.
#
# Notes added during manuscript-code alignment:
#   - Original executable code was preserved except for whitespace/encoding cleanup
#     and the explicitly marked non-run archive described below when applicable.
#   - Comments inserted by this pass are prefixed with "ANNOTATION:".
##############################################################################

# ------------------------------------------------------------------------
# ANNOTATION: Load the finalized APOE-aware harmonized data and set up the
# WGCNA workspace.
# ------------------------------------------------------------------------
##################################################
# 4. Plasma network exploration

rootdir="z:/EBD/grid/4p13b3forAPOEpredict+2ndRegrAgain/"
#rootdir="c:/Users/edammer_4ceb6ff/Downloads/"
setwd(rootdir)
#load("2.saved.image_trait+human_cleanup_nm+em0_V1_3ms_03-27-25.RData")
   #("3.saved.image_trait+human_cleanup2_2fluidSplit+person_ids_CSFandPlasma_(both).RData")


#cleanDat<-readRDS("4p13b4c.normExpr.reg_sites1-19_Fsplit_knownAPOE_only.RDS")
##regvars.b4c still valid, ordered for this cleanDat
##regvars.b4c$APOE.mapped<-gt.APOE.b4c
#numericMeta<-readRDS("4p13b4c.numericMeta_sites1-19_traits_knownAPOE_only.RDS")

load("saved.image-genotype_prediction_finalized.RData") # now includes all QC of unreg and c1 regressed 2x data, and cv_folds
numericMeta<-regvars.c1  # overwrites 22547 (155 samples of unknown age/sex, old)

library(WGCNA)

## Sanity check
e4.bicor.to.siteCorr.b4c.assays<-bicor(t(cleanDat),numericMeta$APOE4.carrier, use='p')
# sort by bicor
e4.bicor.to.siteCorr.b4c.assays[order(unlist(t(e4.bicor.to.siteCorr.b4c.assays)),decreasing=TRUE),][c(1:5,7331:7335)]
#c1 (current 22392 2x regr)
#c1 (current, 22392, age, sex, apoe e4 status protected)
#SPC25|Q9HBM1   LRRN1|Q6UXK5^SL025922@seq.11293.14   CTF1|Q16619^SL002783@seq.13732.79            NEFL|P07196        S100A13|Q99584                TBCA|O75347
#   0.7834346                            0.6711973                           0.4755896      ...    -0.5858044            -0.5916719                 -0.6225905

#b4c (previous 16677 2x regr)
#SPC25|Q9HBM1   LRRN1|Q6UXK5^SL025922@seq.11293.14   CTF1|Q16619^SL002783@seq.13732.79            NEFL|P07196        S100A13|Q99584                TBCA|O75347
#   0.7702713                            0.6725630                           0.4673459      ...    -0.5905291            -0.5929931                 -0.6264711


library(purrr)

table(numericMeta$sample_type)
#Sample
# 22392

#previously:
#Calibrator     Sample
#      1115      27850


# ------------------------------------------------------------------------
# ANNOTATION: Run sanity checks for residual APOE-e4 correlations and
# sample/site counts.
# ------------------------------------------------------------------------
table(numericMeta$is_somalogic)
#    1
#22392

table(numericMeta$contributor_code) #--without MS
# here (PLASMA)
#   A    B    C    D    E    F    G    H    I    J    K    L    M    N    P    Q    R    S    T    U    V    W
#1058 1303 2120  788  678 4076 1330 2254 1411  814  406 1191  991  839 1596 1455  710  215  322 4212 1061  135  with calibrator samples
# 983 1228 2000  743  638 3966 1330 2254 1331  814  406 1101  931  789 1491 1370  670  200  302 4212  996   95  without calibrator samples (missed is_somalogic==1 setting for Calibrator samples)
# 983 1228 2000  743  648 3966 1330 2210 1333  827  406 1101  931  789 1491 1370  670  200  302 4215 1040   95  previously without calibrators but some extra samples compared to now
# 977 1228 1847  715  593 3928 1212 1919 1191           1024                1169  670       204			Prior 16677
# 977 1228 2000  743  638 3966 1330 2254 1192  808  406 1100  931  788 1491 1369  670  200  301


library("doParallel")
#stopCluster(clusterLocal)  #return parallel processing to local workstation if previously set
parallelThreads=32  #set to # of threads on your computer
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="PSOCK")
registerDoParallel(clusterLocal)
enableWGCNAThreads(nThreads=parallelThreads) #speeds the pickSoftThreshold function, outlier removal
  # or allowWGCNAThreads() depending on platform


#<SKIP>
#=============================#

# ------------------------------------------------------------------------
# ANNOTATION: Optional robust outlier-screening code block retained for
# reproducibility but not used for final removal.
# ------------------------------------------------------------------------
#  Check and Remove Outliers  #
#=============================#

numericMeta.withOutliers<-numericMeta


if(!exists("numericMeta")) numericMeta<-traits

sdout=2.5 #Z.k SD fold for outlier threshold
outliers.noOLremoval<-outliers.All<-vector()
cleanDat.noOLremoval<-cleanDat
targets.All=numericMeta

# Define a function for robust standardization (z-scores) using median and MAD.
robust_z <- function(x) {
  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, constant = 1, na.rm = TRUE)
  if(mad_val == 0) return(rep(0, length(x)))
  (x - med) / mad_val
}


for (repeated in 1:20) {
# slow:  normadj <- (0.5+0.5*bicor(cleanDat,use="pairwise.complete.obs")^2)
# Use an approximation instead:

# Apply robust standardization to each column (sample).
cleanDat_std <- apply(cleanDat, 2, robust_z)
# Compute the Pearson correlation on the standardized data.
#cor_mat <- cor(cleanDat_std, use = "pairwise.complete.obs")
cor_mat <- foreach(i = 1:ncol(cleanDat), .combine = cbind) %dopar% {  # ...=cbind, .packages = "WGCNA")
  cor(cleanDat_std[, i, drop = FALSE], cleanDat_std, use = "pairwise.complete.obs")
}

# Create the normalized adjacency matrix.
normadj <- 0.5 + 0.5 * matrix(cor_mat,nrow=ncol(cleanDat_std),ncol=ncol(cleanDat_std),byrow=FALSE)^2


## Calculate connectivity
#netsummary <- fundamentalNetworkConcepts(normadj)
ku <- rowSums(normadj) - diag(normadj)  #much more efficient and equivalent to:  # netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))  #corrected, moved parenthesis open to leading position
## Declare as outliers those samples which are more than sdout sd above the mean connectivity based on the chosen measure
outliers <- (z.ku < mean(z.ku)-sdout*sd(z.ku))  #previously had | z.ku > mean(z.ku)+sdout*sd(z.ku)) to remove outliers on high connectivity end...
print(paste0("There are ",sum(outliers)," outlier samples based on a bicor distance sample network connectivity standard deviation above ",sdout,".  [Round ",repeated,"]"))

cleanDat <- cleanDat[,!outliers]
numericMeta <- targets <- targets.All[!outliers,]
outliers.All<-c(outliers.All,outliers)

if (sum(outliers)==0) break;
} #repeat up to 20 times

#All outliers removed
print(paste0("There are ",sum(outliers.All)," total outlier samples removed in ",repeated," iterations:"))
names(which(outliers.All))
outliersRemoved<-names(which(outliers.All))
#Note outliers as comment below, copied from R session.

#at sdout>2.5
#[1] "There are 0 total outlier samples removed in 1 iterations:"


## Enforce <50% missingness (1 less than half of cleanDat columns (or round down half if odd number of columns))
LThalfSamples<-length(colnames(cleanDat))/2
LThalfSamples<-LThalfSamples - if ((length(colnames(cleanDat)) %% 2)==1) { 0.5 } else { 1.0 }

## If operating on log2(FPKM) data, remove rows with >=50% originally 0 FPKM values (only if there are some rows to be removed)
#IndexHighMissing<-rowsRemoved<-zeroVarRows<-vector()
#temp2<-data.frame(ThrowOut=apply(cleanDat,1,function(x) length(x[x==log2(0+0.05)])>LThalfSamples))
#cleanDat<-cleanDat[!temp2$ThrowOut,]
#dim(cleanDat) #still have x genes, now for y total samples

## If working on log2(protein abundance or ratio) with NA missing values; Enforce <50% missingness (1 less than half of cleanDat columns (or round down half if odd number of columns))
#remove rows with >=50% missing values (only if there are some rows to be removed)
IndexHighMissing<-rowsRemoved<-zeroVarRows<-vector()
temp2<-as.data.frame(cleanDat[which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples),])
#handle condition if temp2 is for one row of cleandat (a vector instead of a data frame)
if (ncol(temp2)==1) {
  temp2<-t(temp2)
  rownames(temp2)=rownames(cleanDat)[which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples)]
}

if (nrow(temp2)>0) { IndexHighMissing=which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples); rowsRemoved<-rownames(cleanDat)[IndexHighMissing]; cleanDat<-cleanDat[-IndexHighMissing,]; }

dim(cleanDat)
#[1] 7334 23568  # no rows removed:

rownames(temp2)
# none


## Write filtered, outlier checked/removed log2(rel abun) matrix to csv file
#write.csv(cleanDat,file="4p3.unreg.Plasma19sites-cleanDat_postTAMPOR_noOL.csv")
#write.csv(numericMeta,file="4p3.unreg.Plasma19sites-numericMeta_postTAMPOR_noOL.csv")


#<END SKIP>


powers <- seq(4,14,by=1)  #initial power check -- try to get SFT.R.sq to go > 0.80
sft <- pickSoftThreshold(t(cleanDat),blockSize=nrow(cleanDat)+1000,   #always calculate power within a single block (blockSize > # of rows in cleanDat)
                         powerVector=powers,
                         corFnc="bicor",networkType="signed")


# ------------------------------------------------------------------------
# ANNOTATION: Select WGCNA soft-threshold power using signed bicor scale-
# free topology diagnostics.
# ------------------------------------------------------------------------
#saveRDS(sft,file="4p4.PLASMA_sft.noGIS.RDS")
saveRDS(sft,file="4p13c1.regr.2PAV+Site(protAgeSexE4).PLASMA_sft.22392.RDS")


##Paste/replace output below and annotate based on graphical plot that follows
## current 22382 2x Regressed data
# ANNOTATION: pasted console output header below was commented to keep the script parseable.
# Power SFT.Rsq  slope truncated.R2 mean.k. median.k.  max.k
#   4   0.811 -2.73         0.917   635.0     538.00  1240
#   5   0.829 -2.07         0.940   388.0     297.00   990
#   6   0.837 -1.70         0.941   250.0     167.00   824  << power=6
#   7   0.847 -1.49         0.929   170.0      95.90   710

#Previous 16677 2x Regressed data
#   4   0.829 -2.90         0.931   621.0     535.00  1180
#   5   0.832 -2.21         0.931   376.0     293.00   940
#   6   0.856 -1.81         0.935   239.0     164.00   785  << power=6
#   7   0.866 -1.56         0.927   161.0      93.90   675

## 4p4b. now TAMPOR mode 1 on calibrator samples; 3415 sampleMedianRows, ***Calibrator samples removed after TAMPOR; power series for SFT:
# ANNOTATION: pasted console output header below was commented to keep the script parseable.
# Power SFT.Rsq  slope truncated.R2 mean.k. median.k.  max.k
#<...>
#   8   0.717 -1.070        0.945   209.0      93.40   949
#   9   0.738 -1.000        0.944   168.0      60.80   875
#  10   0.772 -0.965        0.937   139.0      40.30   811
#  11   0.795 -0.929        0.923   117.0      27.00   757  << power=11
#  12   0.824 -0.916        0.908   101.9      18.40   708
#  13   0.844 -0.896        0.908    88.2      12.60   665


## 4p4. prev TAMPOR mode 1 on calibrator samples; 3415 sampleMedianRows, ***Calibrator samples left in after TAMPOR; power series for SFT:
# ANNOTATION: pasted console output header below was commented to keep the script parseable.
# Power SFT.Rsq  slope truncated.R2 mean.k. median.k.  max.k
#<...>
#   8   0.699 -1.070        0.949   207.0      93.7    944
#   9   0.710 -1.000        0.945   166.0      61.1    870
#  10   0.760 -0.964        0.950   137.0      40.4    806
#  11   0.791 -0.928        0.935   116.0      27.0    752  << power=11
#  12   0.811 -0.905        0.926    99.9      18.3    704
#  13   0.829 -0.891        0.911    87.3      12.6    661

## previous mode4 TAMPOR 3922 sampleMedianRows, power series for SFT:
#<...>
#   8 0.8534091 -1.930777     0.9694948 120.06727  96.988946  425.5853
#<...>


#plot initial SFT.R.sq vs. power curve
tableSFT<-sft[[2]]
plot(tableSFT[,1],tableSFT[,2],xlab="Power (Beta)",ylab="SFT R^2")


#Remove 2 regressed rows with ~0 variance:
badRow.idx<-which(grepl("^HBZ",rownames(cleanDat)) | grepl("^HNRNPA2B1",rownames(cleanDat)))
badRow.idx
#5569 6087
cleanDat<-cleanDat[-badRow.idx,]
dim(cleanDat)
# 7333 22392

# ------------------------------------------------------------------------
# ANNOTATION: Remove the two PAV proxy assay rows before final WGCNA
# construction because their variance was regressed out.
# ------------------------------------------------------------------------


#choose power at elbow of SFT R^2 curve approaching asymptote near or ideally above 0.80
power=6  # (sft$powerEstimate=14)
enforceMMS=FALSE

## Run an automated network analysis (ds=4 and mergeCutHeight=0.07, more liberal)
# choose parameters deepSplit and mergeCutHeight to get respectively more modules and more stringency sending more low connectivity genes to grey (not in modules).
net <- blockwiseModules(t(cleanDat),power=power,deepSplit=4,minModuleSize=10,

# ------------------------------------------------------------------------
# ANNOTATION: Construct a signed bicor WGCNA network with blockwiseModules
# using the selected parameters.
# ------------------------------------------------------------------------
                        mergeCutHeight=0.07,TOMDenom="mean", #detectCutHeight=0.9999,                        #TOMDenom="mean" may get more small modules here.  NOTE: CAPITAL "D"
                        corType="bicor",networkType="signed",pamStage=TRUE,pamRespectsDendro=TRUE,
                        verbose=3,saveTOMs=FALSE,maxBlockSize=nrow(cleanDat)+1000,reassignThresh=0.05)       #maxBlockSize always more than the number of rows in cleanDat
#blockwiseModules can take 30 min+ for large numbers of gene products/proteins (10000s of rows); much quicker for smaller proteomic data sets

nModules<-length(table(net$colors))-1
modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
as.data.frame(cbind(orderedModules,Size=modules))

##copy R session output;
#22392 deepSplit=4; minModSize=10 TOMDenom="mean" (current)
#<...>
#  M12            tan   24

#16677 deepSplit=4; minModSize=10 TOMDenom="mean" (current)
#<...>
#  M12            tan   24 (current 16677 2x regressed, site, protect Age, sex, APOE e4 carrier +/-)

#deepSplit=4; minModSize=10 TOMDenom="mean" (previous)
#<...>
#  M27          white   10  (blood plasma post TAMPOR mode 1) - no outliers found at 2.5 SD (Calibrator samples out)

#<...>
#  M28        skyblue   10  (blood plasma post TAMPOR mode 1) - no outliers found at 2.5 SD (Calibrator samples n=1075 left in)

#<...>
#  M17         grey60   17  (blood plasma post TAMPOR mode 4)

#  M68     orangered3   19  (previously, CSF)

#deepsplit=4; minModSize=10 TOMdenom=.. ("min")
#<...>
#  M41 lightsteelblue1   18  (previously, last module)
#  M48   darkslateblue   19  (now, with (7334x3233)_4403sampleMedianRows)

#net.ds4.mms10<-net  #with TOMDenom="mean"
net.ds4.mms10.TOMDenomMean.2xRegr.122392samp<-net


#we will explore the blockwiseModules() function-built network with parameter deepSplit=4, minimum (initial) module size of 10 and TOMDenom="min"
net<-net.ds4.mms10.TOMDenomMean.2xRegr.16677samp


minModSize=10
# If necessary, return module members of small modules below size minSize=X to grey
if (enforceMMS) {
  removedModules<-orderedModules[which(modules<minModSize),"Color"]
  for(i in removedModules) { net$colors[net$colors==i] <- "grey" }
  for(i in removedModules) { net$MEs[,paste0("ME",i)] <- NULL }

  nModules<-length(table(net$colors))-1
  modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
  orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
  modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
  as.data.frame(cbind(orderedModules,Size=modules))
}
minModSize=10


#calculate kME table up front, in case we need to correct color assignments
MEs<-tmpMEs<-data.frame()
MEList = moduleEigengenes(t(cleanDat), colors = net$colors)
MEs = orderMEs(MEList$eigengenes)
net$MEs <- MEs
colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
rownames(MEs)<-rownames(numericMeta)

tmpMEs <- MEs #net$MEs
colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
MEs[,"grey"] <- NULL
tmpMEs[,"MEgrey"] <- NULL

# ------------------------------------------------------------------------
# ANNOTATION: Calculate module eigengenes and signed kME values used for
# module membership refinement.
# ------------------------------------------------------------------------

kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc="bicor")


table(net$colors)["grey"]
# 2436  #  2497  #  1708  #  1463  #  1618  #  1243  # 1084   # previously 1883
paste0(round(table(net$colors)["grey"]/nrow(cleanDat)*100,2),"% grey")
# 33.22% grey  #  34.05% grey  #  23.29% grey  (mode1 TAMPOR calibrators out) #  19.95% grey (plasma mode1 TAMPOR calibrators left in)  #  22.06% grey (plasma mode4 TAMPOR)  #  16.95% grey  # 14.78% grey  # previously (4a) 17.85% grey


##ITERATIVE until condition met that all module membes are at least 0.28 kMEintramodule.
#Go back and do final algorithm fix of module colors (remove kMEintramodule<0.28 members, reassign grey with kMEintramodule>0.35; max difference from kMEmax<0.10)

retry=TRUE;
kMEmaxDiff=0.1
reassignIfGT=0.30
greyIfLT=0.30
iter=1;
while (retry) {
  cat(paste0("\nkME table Cleanup, processing iteration ",iter,"..."))
  colorVecFixed<-colorVecBackup<-net$colors
  orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
  kMEintramoduleVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]) )  #all sig digits (no rounding), so max will be unique.
  colorVecFixed[kMEintramoduleVector<greyIfLT]<-"grey"
  kMEmaxVec<-apply( as.data.frame(kMEdat),1,function(x) max(x) )
  kMEmaxColorsVec<-apply( as.data.frame(cbind(kMEmaxVec,kMEdat)),1, function(x) gsub("kME","",colnames(kMEdat)[which(x==x[1])[2]-1]) )
  kMEintramoduleVector<-unlist(lapply(kMEintramoduleVector,function(x) if(length(x)==0) { 1 } else { x }))   #grey will be ignored in checking for kMEmaxVec-kMEintramoduleVector difference max
  kMEmaxDiffTooBig<-(kMEmaxVec-kMEintramoduleVector) >= kMEmaxDiff
  colorVecFixed[which( (colorVecFixed=="grey" & kMEmaxVec>reassignIfGT) | kMEmaxDiffTooBig )] <- kMEmaxColorsVec[which( (colorVecFixed=="grey" & kMEmaxVec>reassignIfGT) | kMEmaxDiffTooBig )]
  net$colors<-colorVecFixed

#  table(net$colors)["grey"]  #decreased to x


# Are colors still in rank order? -- put them in order by recoloring modules that changed rank
  sort(table(net$colors),decreasing=TRUE)[!names(sort(table(net$colors),decreasing=TRUE))=="grey"]

  oldcolors <- names(sort(table(net$colors),decreasing=TRUE)[!names(sort(table(net$colors),decreasing=TRUE))=="grey"])
  for (i in 1:length(oldcolors)) {
    net$colors[net$colors==oldcolors[i]]<-paste0("proxy",labels2colors(i))
  }
  for (i in 1:length(oldcolors)) {
    net$colors[net$colors==paste0("proxy",labels2colors(i))]<-labels2colors(i)
  }

# one can check that colors are in order by size now
  #sort(table(net$colors),decreasing=TRUE)[!names(sort(table(net$colors),decreasing=TRUE))=="grey"]

# recalculate kME table, since we have corrected color assignments
  MEs<-tmpMEs<-data.frame()
  MEList = moduleEigengenes(t(cleanDat), colors = net$colors, verbose=0)
  MEs = orderMEs(MEList$eigengenes)
  net$MEs <- MEs
  colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
  rownames(MEs)<-rownames(numericMeta)

  tmpMEs <- MEs #net$MEs
  colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
  MEs[,"grey"] <- NULL
  tmpMEs[,"MEgrey"] <- NULL

  kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc="bicor")

# recheck min kMEintramodule and max diff from kMEmax
  nModules<-length(table(net$colors))-1
  modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
  orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
  orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
  kMEsIntramoduleVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) if(!x[1]=="grey") { paste0(round(as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]),4)) } else { 1 } ) #grey proteins set to dummy value of 1 (ignore)

  kMEmaxVec<-apply( as.data.frame(kMEdat),1,function(x) max(x) )
  kMEintramoduleVector<-unlist(lapply(kMEintramoduleVector,function(x) if(length(x)==0) { 1 } else { x }))   #grey will be ignored in checking for kMEmaxVec-kMEintramoduleVector difference max
  kMEmaxDiffCalc<- kMEmaxVec-kMEintramoduleVector
  if (min(kMEsIntramoduleVector)>=greyIfLT & max(kMEmaxDiffCalc)<=kMEmaxDiff) { cat(paste0("\nkME table 'clean' in ",iter," iterations.")); retry=FALSE; }
  iter=iter+1
  if (iter>30) break; #**
}
#** breaks after iteration 30 if did not reach criteria.


nModules<-length(table(net$colors))-1
modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
as.data.frame(cbind(orderedModules,Size=modules))


# Final modules
#cleanDat, pwr=6, ds=4 mms=10 TOMDenom="mean"; 2x regr 22392; clean in 12 iterations:
#  turquoise     M1   turquoise 1838
#       blue     M2        blue  804
#      brown     M3       brown  777
#     yellow     M4      yellow  669
#      green     M5       green  288
#        red     M6         red  269
#      black     M7       black  153
#       pink     M8        pink  133
#    magenta     M9     magenta   81
#     purple    M10      purple   80
#greenyellow    M11 greenyellow   78
#        tan    M12         tan   56


#cleanDat, pwr=6, ds=4 mms=10 TOMDenom="mean"; 2x regr 16677; clean in 8 iterations:
#  turquoise     M1   turquoise 1743
#       blue     M2        blue  821
#      brown     M3       brown  813
#     yellow     M4      yellow  682
#      green     M5       green  280
#        red     M6         red  257
#      black     M7       black  158
#       pink     M8        pink  126
#    magenta     M9     magenta   81
#     purple    M10      purple   80
#greenyellow    M11 greenyellow   79
#        tan    M12         tan   57


#cleanDat, pwr=11 ds=4 mms=10 TOMDenom="mean; TAMPOR mode 1;
#<...>
#  M27          white   15  mode1 TAMPOR calibrators out (unreg)

#<...>
#  M28       skyblue3   15  mode1 TAMPOR calibrators in (unreg)

#cleanDat, pwr=6 ds=4 mms=10 TOMDenom="mean"; 9 iter cleanup complete
#<...>
#  M17         grey60   62  (now, blood plasma post TAMPOR)


table(net$colors)["grey"]
# 2107  #  2156  #  527  #  458  #  1126  now-plasma  #  463  # 381  # previously (4a) 647
paste0(round(table(net$colors)["grey"]/nrow(cleanDat)*100,2),"% grey")
# 28.73% grey  #  29.4% grey  #  7.18% grey  #  6.24% grey  #  15.35% grey  #  6.31% grey  # 5.19% grey  # previously (4a) "6.13% grey"

projectFilesOutputTag="mms10.ds4.pwr6"
## saved image of R session after running and finalizing blockwiseModules() function WGCNA output (now includes net data structure)
#not saved: save.image(paste0("4.saved.image.Plasma16677.WGCNA.",projectFilesOutputTag,".RData"))  #overwrites
#load("4b.saved.image.Plasma16677.WGCNA.mms10.ds4.pwr6.RData")


rootdir="z:/EBD/grid/4p13b3forAPOEpredict+2ndRegrAgain/"
outputfigs<-outputtabs<- rootdir
# setwd(rootdir)


## Traits finalization/breakout


# ------------------------------------------------------------------------
# ANNOTATION: Generate global network trait correlation figures and order
# modules/eigenproteins.
# ------------------------------------------------------------------------

#numericMeta.Shijia<-regvars.b4c[sort(c(b4c.e33.subset.idx,b4c.e44.subset.idx)),  # 2587+510, kept in order of contributor site (A-T)
#                                c("person_id","visit","sequential_visit_number","age_at_visit","Sex.int","raceAA","weight_kg","bmi","height_cm","years_of_education",  # demographics (longitudinal); raceAA is binary for African American (1); Sex.int is integer 0/1 F/M
#                                  "contributor_code","contributor_Fsplit",       # Contributor Sites (regressed including site F split to F1, F2, F3)
#                                  "resting_heart_rate_pulse","systolic_blood_pressure_sitting","diastolic_blood_pressure_sitting",                                     # heart health, numeric parameters
#                                  "hypertension","stroke","tia","tbi","diabetes","chf","copd","mi","afib","hyperlipidaemia","depression","anxiety", "alcohol_hx","smoking_hx","total_years_smoked",  # comorbidities (binary status, except years smoked)
#                                  "MMSE","cdr",           # Cognitive scores; preferred MMSE imputed from MoCA with education adjustment using Fasnacht et al 2022 rubric
#                                  "Group.withCTimputed",  # diagnosis, including those mapped from our Nat Aging cohorts, and Controls with 0 CDR and MMSE >=28 imputed
#                                  "C9Orf72","GRN","MAPT", # selected genetic binary mutation status for detrimental mutations
#                                  "sample_matrix",        # citrate or EDTA tube-collected plasma
#                                  "Lilly.BH.blood.pTau217","UDS.blood.pTau217","AmyloidPositivity.withRM",  # Pathological biomarkers for Nat Aging cohorts
#                                  "APOE.mapped")]          # APOE genotype in "e#/e#" string format -- as provided by GNPC HDS v1.3, and mapped for missing values in Nat Aging 3 cohorts [ROSMAP (site R), UDS (site D), and BioHermes (site A)]


# Amyloid Positivity (already set)
#numericMeta$AmyloidPositivity.01<-NA
#numericMeta$AmyloidPositivity.01[which(numericMeta$AmyloidPositivity.withRM=="NEGATIVE")]<-0
#numericMeta$AmyloidPositivity.01[which(numericMeta$AmyloidPositivity.withRM=="POSITIVE")]<-1


## MMSE values above 30 are only remaining trait not ready for Global Network Plots
length(which(numericMeta$MMSE>30))
# 0; previously 24
#numericMeta$MMSE[which(numericMeta$MMSE>30)]<-NA


## Output GlobalNetworkPlots and kMEtable
####################################################################################################################
FileBaseName=paste0("Plasma_(7333x22392)_WGCNA.",projectFilesOutputTag)


library(Cairo)
CairoPDF(file=paste0(outputfigs,"4p13c1.GlobalNetworkPlots-part1-",FileBaseName,".pdf"),width=16,height=14)

## Plot dendrogram with module colors and trait correlations
MEs<-tmpMEs<-data.frame()
MEList = moduleEigengenes(t(cleanDat), colors = net$colors)
MEs = orderMEs(MEList$eigengenes)
colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
rownames(MEs)<-rownames(numericMeta)

#numericIndices<-sort(unique(c( which(!is.na(apply(numericMeta,2,function(x) sum(as.numeric(x))))), which(!(apply(numericMeta,2,function(x) sum(as.numeric(x),na.rm=T)))==0) )))
##Warnings OK; This determines which traits are numeric and if forced to numeric values, non-NA values do not sum to 0

# Specify manually numericIndices -- many are too sparse to be useful for heatmap:
#numericNames<-c("age_at_visit","sex","raceAA","years_of_education","height_cm","weight_kg","bmi","resting_heart_rate_pulse","systolic_blood_pressure_sitting","diastolic_blood_pressure_sitting","alcohol_hx","smoking_hx","total_years_smoked","bmi","stroke","tia","tbi","recruited_control","ad","ftd","pd","als","mci_sci","copd","mi","afib","depression","anxiety","cdr","MMSE","MoCA","visit","sequential_visit_number","C9Orf72","GRN","MAPT","APOE4.Dose")
numericNames<-c("age_at_visit","sex","raceAA","years_of_education","height_cm","weight_kg","bmi","resting_heart_rate_pulse","systolic_blood_pressure_sitting","diastolic_blood_pressure_sitting","hypertension","alcohol_hx","smoking_hx","total_years_smoked","tia","tbi","stroke","chf","mi","afib","angina","hyperlipidaemia","recruited_control","ad","ftd","pd","als","mci_sci","mi","depression","anxiety","cdr","MMSE","MoCA","sequential_visit_number","C9Orf72","GRN","MAPT","APOE4.Dose","Lilly.BH.blood.pTau217","UDS.blood.pTau217","AmyloidPositivity.01","RegrBloodPreanalyticFactor.HNRNPA2B1","RegrBloodPreanalyticFactor.HBZ","TimeToSpin","TimeToDecant","TimeToFreeze","FedFastedTime","FreezeThawCycles")
numericIndices<-match(numericNames,colnames(numericMeta))   #c(2,4:6,11,12,14,16,18,20,22,24,26,28,30:32,34:53,55,56,62,63,66)]


geneSignificance <- cor(sapply(numericMeta[,numericIndices],as.numeric),t(cleanDat),use="pairwise.complete.obs")
rownames(geneSignificance) <- colnames(numericMeta)[numericIndices]
geneSigColors <- t(numbers2colors(t(geneSignificance),,signed=TRUE,lim=c(-1,1),naColor="black"))
rownames(geneSigColors) <- colnames(numericMeta)[numericIndices]

par(mar=c(2,12,2,3))
plotDendroAndColors(dendro=net$dendrograms[[1]],
                    colors=t(rbind(net$colors,geneSigColors)),
                    cex.dendroLabels=1.2,addGuide=TRUE,
                    dendroLabels=FALSE,
                    marAll=c(2,12,2,2),
                    groupLabels=c("Module Colors",colnames(numericMeta)[numericIndices]))

## Plot eigengene dendrogram/heatmap - using bicor
tmpMEs <- MEs #net$MEs
colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
MEs[,"grey"] <- NULL
tmpMEs[,"MEgrey"] <- NULL

plotEigengeneNetworks(tmpMEs, "Eigengene Network", marHeatmap = c(3,4,2,2), marDendro = c(0,4,2,0),plotDendrograms = TRUE, xLabelsAngle = 90,heatmapColors=blueWhiteRed(50))


######################
## Find differences between Groups (as defined in Traits input); Finalize Grouping of Samples for ANOVA

table(numericMeta$Group)
#  AD  ALS   CT  FTD  MCI   PD
#1173  190 1156   78 1422  648

#Set a vector of strings that represent each sample in order, calling out each sample as a member of named groups (used by GlobalNetworkPlot boxplots, and later, ANOVA DiffEx)
Grouping<-numericMeta$Group  #here, we will calculate the P value for one-way ANOVA, controlling age and sex; typically there is a column "Group" loaded as a column in the traits.csv file
#Grouping[numericMeta$Diagnosis==0]<-"Normal"  #only necessary if Group was numerically encoded; does nothing if the Grouping vector has no numeric values
#Grouping[numericMeta$Diagnosis==1]<-"AD"


# This gets one-way ANOVA (if ranked, Kruskal-Wallis) nonparametric p-values for groupwise comparison of interest.
# look at numericMeta (traits data) and choose traits to use for linear model-determination of p value
head(numericMeta)
# Change below line to point to a factored trait, which will define groups for ANOVA
regvars <- data.frame(as.factor( numericMeta$Group ) ) #, as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
colnames(regvars) <- c("Group") #,"Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group,data=regvars) # any second or later variable effects are removed by the linear model

pvec.group <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.group[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.group) <- colnames(MEs)


# Group, age+sex adjusted
regvars <- data.frame(as.factor( numericMeta$Group ) , as.numeric(numericMeta$age_at_visit), as.factor(abs(numericMeta$sex-2)))
colnames(regvars) <- c("Group","Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group +Age+Sex, data=regvars) # any second or later variable effects are removed by the linear model

pvec.group.ageSexAdj <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.group.ageSexAdj[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.group.ageSexAdj) <- colnames(MEs)


# Group as Batch
regvars <- data.frame(as.factor( numericMeta$contributor_code ) ) #, as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
colnames(regvars) <- c("Group") #,"Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group,data=regvars) # any second or later variable effects are removed by the linear model

pvec.batch <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.batch[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.batch) <- colnames(MEs)


# Group as C9Orf72 mutation status
regvars <- data.frame(as.factor( numericMeta$C9Orf72 ) ) #, as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
colnames(regvars) <- c("Group") #,"Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group,data=regvars) # any second or later variable effects are removed by the linear model

pvec.C9Orf72 <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.C9Orf72[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.C9Orf72) <- colnames(MEs)


# Group as GRN mutation status
regvars <- data.frame(as.factor( numericMeta$GRN ) ) #, as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
colnames(regvars) <- c("Group") #,"Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group,data=regvars) # any second or later variable effects are removed by the linear model

pvec.GRN <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.GRN[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.GRN) <- colnames(MEs)


# Group as MAPT mutation status
regvars <- data.frame(as.factor( numericMeta$MAPT ) ) #, as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
colnames(regvars) <- c("Group") #,"Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group,data=regvars) # any second or later variable effects are removed by the linear model

pvec.MAPT <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.MAPT[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.MAPT) <- colnames(MEs)


# Group as TIA status
regvars <- data.frame(as.factor( numericMeta$tia ) ) #, as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
colnames(regvars) <- c("Group") #,"Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group,data=regvars) # any second or later variable effects are removed by the linear model

pvec.tia <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.tia[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.tia) <- colnames(MEs)


# Group as TBI status
regvars <- data.frame(as.factor( numericMeta$tbi ) ) #, as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
colnames(regvars) <- c("Group") #,"Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group,data=regvars) # any second or later variable effects are removed by the linear model

pvec.tbi <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.tbi[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.tbi) <- colnames(MEs)


# Group as MI myocardial infarct status
regvars <- data.frame(as.factor( numericMeta$mi ) ) #, as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
colnames(regvars) <- c("Group") #,"Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group,data=regvars) # any second or later variable effects are removed by the linear model

pvec.mi <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.mi[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.mi) <- colnames(MEs)


# Group as hyperlipidaemia status
regvars <- data.frame(as.factor( numericMeta$hyperlipidaemia ) ) #, as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
colnames(regvars) <- c("Group") #,"Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group,data=regvars) # any second or later variable effects are removed by the linear model

pvec.hyperlipidaemia <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.hyperlipidaemia[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.hyperlipidaemia) <- colnames(MEs)


# Group as hypertension status
regvars <- data.frame(as.factor( numericMeta$hypertension ) ) #, as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
colnames(regvars) <- c("Group") #,"Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group,data=regvars) # any second or later variable effects are removed by the linear model

pvec.hypertension <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.hypertension[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.hypertension) <- colnames(MEs)


# Group as APOE e4 Dose mutation status
regvars <- data.frame(as.factor( numericMeta$APOE4.Dose ) ) #, as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
colnames(regvars) <- c("Group") #,"Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group,data=regvars) # any second or later variable effects are removed by the linear model

pvec.APOE4.Dose <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.APOE4.Dose[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.APOE4.Dose) <- colnames(MEs)


##[ pvec calculation, for 2 new boxplots ]############################
#Set a vector of strings that represent each sample in order, calling out each sample as a member of named groups (used by GlobalNetworkPlot boxplots, and later, ANOVA DiffEx)

# Groups including 3 mapped cohorts' Group (1, RM); or Group.pathCog (2, BH, UDS)
regvars <- data.frame(as.factor( numericMeta$Group.3mappedCohorts ) ) #, as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
colnames(regvars) <- c("Group") #,"Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group,data=regvars) # any second or later variable effects are removed by the linear model

pvec.group.3mappedCohorts <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.group.3mappedCohorts[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.group.3mappedCohorts) <- colnames(MEs)


# Groups including 3 mapped cohorts' Group (1, RM); or Group.pathCog (2, BH, UDS), age+sex adjusted
regvars <- data.frame(as.factor( numericMeta$Group.3mappedCohorts ) , as.numeric(numericMeta$age_at_visit), as.factor(abs(numericMeta$sex-2)))
colnames(regvars) <- c("Group","Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group +Age+Sex, data=regvars) # any second or later variable effects are removed by the linear model

pvec.group.3mappedCohorts.ageSexAdj <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.group.3mappedCohorts.ageSexAdj[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.group.3mappedCohorts.ageSexAdj) <- colnames(MEs)


# Groups including 3 mapped cohorts' Group (1, RM); or Group.pathCog (2, BH, UDS)
regvars <- data.frame(as.factor( numericMeta$Group.pathCog.mapped ) ) #, as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
colnames(regvars) <- c("Group") #,"Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group,data=regvars) # any second or later variable effects are removed by the linear model

pvec.group.pathCog.mapped <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.group.pathCog.mapped[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.group.pathCog.mapped) <- colnames(MEs)


# Groups including 3 mapped cohorts' Group (1, RM); or Group.pathCog (2, BH, UDS), age+sex adjusted
regvars <- data.frame(as.factor( numericMeta$Group.pathCog.mapped ) , as.numeric(numericMeta$age_at_visit), as.factor(abs(numericMeta$sex-2)))
colnames(regvars) <- c("Group","Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group +Age+Sex, data=regvars) # any second or later variable effects are removed by the linear model

pvec.group.pathCog.mapped.ageSexAdj <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.group.pathCog.mapped.ageSexAdj[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.group.pathCog.mapped.ageSexAdj) <- colnames(MEs)


# Groups including 3 mapped cohorts' Group (1, RM); or Group.pathCog (2, BH, UDS) and CT imputed
regvars <- data.frame(as.factor( numericMeta$Group.withCTimputed ) ) #, as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
colnames(regvars) <- c("Group") #,"Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group,data=regvars) # any second or later variable effects are removed by the linear model

pvec.group.withCTimputed <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.group.withCTimputed[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.group.withCTimputed) <- colnames(MEs)


# Groups including 3 mapped cohorts' Group (1, RM); or Group.pathCog (2, BH, UDS) and CT imputed, age+sex adjusted
regvars <- data.frame(as.factor( numericMeta$Group.withCTimputed ) , as.numeric(numericMeta$age_at_visit), as.factor(abs(numericMeta$sex-2)))
colnames(regvars) <- c("Group","Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group +Age+Sex, data=regvars) # any second or later variable effects are removed by the linear model

pvec.group.withCTimputed.ageSexAdj <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.group.withCTimputed.ageSexAdj[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.group.withCTimputed.ageSexAdj) <- colnames(MEs)


# Amyloid Positivity 3 cohorts curated trait
regvars <- data.frame(as.factor( numericMeta$AmyloidPositivity.withRM ) ) #, as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
colnames(regvars) <- c("Group") #,"Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group,data=regvars) # any second or later variable effects are removed by the linear model

pvec.AmyloidPositivity.withRM <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.AmyloidPositivity.withRM[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.AmyloidPositivity.withRM) <- colnames(MEs)


# Apoe genotypes (n=6) - all cases
regvars <- data.frame(as.factor( numericMeta$APOE.mapped.predicted ) ) #, as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
colnames(regvars) <- c("Group") #,"Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~Group,data=regvars) # any second or later variable effects are removed by the linear model

pvec.APOEgeno.allSamp <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.APOEgeno.allSamp[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.APOEgeno.allSamp) <- colnames(MEs)


# Apoe genotypes (n=6) - AD+MCI only
regvars <- data.frame(as.factor( numericMeta$APOE.mapped.predicted[which(Grouping=="AD" | Grouping=="MCI")] ) ) #, as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
colnames(regvars) <- c("Group") #,"Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs[which(Grouping=="AD" | Grouping=="MCI"),])~Group,data=regvars) # any second or later variable effects are removed by the linear model

pvec.APOEgeno.ADorMCIonly <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.APOEgeno.ADorMCIonly[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.APOEgeno.ADorMCIonly) <- colnames(MEs)


# Apoe genotypes (n=6) - Control only
regvars <- data.frame(as.factor( numericMeta$APOE.mapped.predicted[which(Grouping=="CT")] ) ) #, as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
colnames(regvars) <- c("Group") #,"Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
##aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs[which(Grouping=="CT"),])~Group,data=regvars) # any second or later variable effects are removed by the linear model

pvec.APOEgeno.CTonly <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.APOEgeno.CTonly[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.APOEgeno.CTonly) <- colnames(MEs)


######################
## Get sigend kME values
kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc="bicor")


######################
## Plot eigengene-trait correlations - using p value of bicor for heatmap scale
library(RColorBrewer)
MEcors <- bicorAndPvalue(MEs,numericMeta[,numericIndices])
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p


textMatrix = apply(moduleTraitCor,2,function(x) signif(x, 2))
#textMatrix = paste(signif(moduleTraitCor, 2), " (",
#  signif(moduleTraitPvalue, 1), ")", sep = "");
#dim(textMatrix) = dim(moduleTraitCor)
# par(mfrow=c(1,1))
# par(mar = c(6, 8.5, 3, 3));

par(mar=c(6, 12, 3, 3) )
par(mfrow=c(2,1))

## Display the p value heatmap with text correlation values #modules on x, transposed - same as next bicor heatmap
cexx <- if(nModules>75) { 0.8 } else { 1 }
#rowMin(moduleTraitPvalue) # if we want to resort rows by min P value in the row
xlabAngle <- if(nModules>75) { 90 } else { 45 }
cex.text <- if(ncol(MEs)>50) { 0.3 } else { 0.5 }

labelMat<-matrix(nrow=(length(names(MEs))), ncol=2,data=c(rep(1:(length(names(MEs)))),labels2colors(1:(length(names(MEs))))))
labelMat<-labelMat[match(names(MEs),labelMat[,2]),]
for (i in 1:(length(names(MEs)))) { labelMat[i,1]<-paste("M",labelMat[i,1],sep="") }
for (i in 1:length(names(MEs))) { labelMat[i,2]<-paste("ME",labelMat[i,2],sep="") }

colvec <- rep("white",1500)
colvec[1:500] <- colorRampPalette(rev(brewer.pal(8,"BuPu")[2:8]))(500)
colvec[501:1000]<-colorRampPalette(c("white",brewer.pal(8,"BuPu")[2]))(3)[2] #interpolated color for 0.05-0.1 p
labeledHeatmap(Matrix = t(apply(moduleTraitPvalue,2,as.numeric)),
               xLabels = labelMat[,2], # paste0("ME",names(MEs)),
               yLabels = colnames(numericMeta)[numericIndices],
               xSymbols = labelMat[,1], # names(MEs),
               xColorLabels = TRUE,
               xLabelsAngle = xlabAngle,
               colors = colvec,
               textMatrix = t(textMatrix),
               setStdMargins = FALSE,
               cex.text = cex.text,
               cex.lab.x= cexx,
               zlim = c(0,0.15),
               main = paste("Module-trait relationships\n bicor r-value shown as text\nHeatmap scale: Student correlation p value"),
               cex.main=0.8)


######################
## Plot eigengene-trait heatmap custom - using bicor color scale

numericMetaCustom<-numericMeta[,numericIndices]
MEcors <- bicorAndPvalue(MEs,numericMetaCustom)
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p

moduleTraitPvalue.txt<-signif(moduleTraitPvalue, 1)
moduleTraitPvalue.txt[moduleTraitPvalue.txt > as.numeric(0.05)]<-as.character("")

textMatrix = moduleTraitPvalue.txt; #paste(signif(moduleTraitCor, 2), " / (", moduleTraitPvalue, ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
                                #textMatrix = gsub("()", "", textMatrix,fixed=TRUE)

# par(mar=c(16, 12, 3, 3) )
# par(mfrow=c(1,1))

bw<-colorRampPalette(c("#0058CC", "white"))
wr<-colorRampPalette(c("white", "#CC3300"))

colvec<-c(bw(50),wr(50))

labeledHeatmap(Matrix = t(moduleTraitCor)[,],
               yLabels = colnames(numericMetaCustom),
               xLabels = labelMat[,2],
               xSymbols = labelMat[,1],
               xColorLabels=TRUE,
               colors = colvec,
               textMatrix = t(textMatrix)[,],
               setStdMargins = FALSE,
               cex.text = cex.text,
               cex.lab.x = cexx,
               xLabelsAngle = xlabAngle,
               verticalSeparator.x=c(rep(c(1:length(colnames(MEs))),as.numeric(ncol(MEs)))),
               verticalSeparator.col = 1,
               verticalSeparator.lty = 1,
               verticalSeparator.lwd = 1,
               verticalSeparator.ext = 0,
               horizontalSeparator.y=c(rep(c(1:ncol(numericMetaCustom)),ncol(numericMetaCustom))),
               horizontalSeparator.col = 1,
               horizontalSeparator.lty = 1,
               horizontalSeparator.lwd = 1,
               horizontalSeparator.ext = 0,
               zlim = c(-1,1),
               main = "Module-trait Relationships\n Heatmap scale: signed bicor r-value", # \n (Signif. p-values shown as text)"),
               cex.main=0.8)
dev.off()


## Plot annotated heatmap - annotate all the metadata, plot the eigengenes!
CairoPDF(file=paste0(outputfigs,"4p13c1.GlobalNetworkPlots-part2-",FileBaseName,".pdf"),width=16,height=12)


# This is where we will first use the Grouping vector of string group descriptions we set above.
toplot <- MEs

colnames(toplot) <- colnames(MEs)
rownames(toplot) <- rownames(MEs)
toplot <- t(toplot)

# Windsorize ME data for heatmap
minExtreme=min(abs(range(toplot,na.rm=T)))
toplot[toplot> minExtreme]<- minExtreme
toplot[toplot< -minExtreme]<- -minExtreme

colnames(toplot)<-NULL

pvec <- pvec.group[match(names(pvec.group),rownames(toplot))]
#rownames(toplot) <- paste(rownames(toplot),"\np = ",signif(pvec,2),sep="")
rownames(toplot) <- paste(orderedModules[match(colnames(MEs),orderedModules[,2]),1]," ",rownames(toplot),"  |  AOV p=",signif(pvec,2),sep="")

# add any traits of interest you want to be in the legend
Gender=as.numeric(numericMeta$sex)
Gender[Gender==2]<-"Female"
Gender[Gender==1]<-"Male"

African=as.numeric(numericMeta$raceAA)
African[African==0]<-"Caucasian"
African[African==1]<-"African"

metdat=data.frame(Group.all=numericMeta$Group.withCTimputed,Group=numericMeta$Group, Site=numericMeta$contributor_code, Age=as.numeric(numericMeta$age_at_visit), Sex=Gender, MMSE=as.numeric(numericMeta$MMSE), TwoRaces=African)

# set colors for the traits in the legend
heatmapLegendColors=list('Group.all'=c('red','greenyellow','yellow','thistle2','blue','darkorange','hotpink','darkviolet'),  # AD ALS AsymAD CI.Other CT FTD MCI PD
                         'Group'=c('red','greenyellow','blue','darkorange','hotpink','darkviolet'),  # AD ALS CT FTD MCI PD
                         'Site'=labels2colors(141:156),  #:161),  #c('red','pink','hotpink','white','lightyellow'), #
                         'Age'=c('white','darkviolet'), # low to high
                         'Sex'=c("pink","dodgerblue"),  # Female, Male
                         'MMSE'=c('lightgreen','darkslateblue'),
                         'TwoRaces'=c('chocolate3','paleturquoise'),
                         'Modules'=sort(colnames(MEs)))

library(NMF)
# par(mfrow=c(2,1))
#par(mar=c(3, 12, 1, 3) )  # changed bottom and top margins from above page output
layout(matrix(c(1,1, 2,3), nrow = 2, ncol = 2, byrow=TRUE),
       heights = c(0.95,1.3), # Heights of the rows
       widths = c(0.88,0.12)) # Widths of the columns  -- the distance to squash the second plot to the left, (because we do not duplicate legends)
# sapply(c(1:3),layout.show)

# a way to increase right margin of lower plot, since we are not plotting annLegend again.

aheatmap(x=toplot, ## Numeric Matrix
         main="Plot of Eigengene-Trait Relationships - SAMPLES IN ORIG. BATCH ORDER",
         annCol=metdat,
         annRow=data.frame(Modules=colnames(MEs)),
         annColors=heatmapLegendColors,
         border=list(matrix = TRUE),
         scale="row",
         distfun="correlation",hclustfun="average", ## Clustering options
         cexRow=0.8, ## Character sizes
         cexCol=0.8,
         labCol=NA,  # for this data, labels are not informative
         col=blueWhiteRed(100), ## Color map scheme
         treeheight=80,
         Rowv=TRUE, Colv=NA) ## Do not cluster columns - keep given order


### THIS PLOT WILL NOT COMPLETE PLOTTING IN 12+ HOURS -- skipped
##source("../samePage.aheatmap.below.R")
##aheatmap.noAnnLegend.sameLayout(... low-level functions not accessible/exported
#aheatmap(x=toplot, ## Numeric Matrix
#         main="Plot of Eigengene-Trait Relationships - SAMPLES CLUSTERED",
#         annCol=metdat,
#         annRow=data.frame(Modules=colnames(MEs)),
#         annColors=heatmapLegendColors,
#         annLegend=FALSE,  #second plot on same age, same trait metadata tracks.
#         legend=FALSE,
#         border=list(matrix = TRUE),
#         scale="row",
#         distfun="correlation",hclustfun="average", ## Clustering options
#         cexRow=0.8, ## Character sizes
#         cexCol=0.8,
#         labCol=NA,  # for this data, labels are not informative
#         col=blueWhiteRed(100), ## Color map scheme
#         treeheight=80,
#         Rowv=TRUE,Colv=TRUE) ## Cluster columns

dev.off()


#deconvProportions<-read.csv(file="../Deconvolution/2.Ensemble-13cellTypes-bulkRNA_Raj59-proportionEstimates.csv",header=TRUE,row.names=1,check.names=F)
#deconvProportions<-deconvProportions[match(numericMeta$RNAsample,rownames(deconvProportions)),]
#rownames(deconvProportions)<-rownames(numericMeta)
#MEpropCors <- bicorAndPvalue(MEs,deconvProportions)
#modulePropCor <- MEpropCors$bicor
#modulePropPvalue <- MEpropCors$p


######################################
## INDIVIDUAL ME BOXPLOTS AND SCATTERPLOTS - 1 page / ME

## Change the below code in the for loop using the following session output

toplot <- MEs

colnames(toplot) <- colnames(MEs)
rownames(toplot) <- rownames(MEs)
toplot <- t(toplot)


#numericIndices.more<-c(which(colnames(numericMeta) %in% c("RegrBloodPreanalyticFactor.HNRNPA2B1","RegrBloodPreanalyticFactor.HBZ")))
#MEcors.more <- bicorAndPvalue(MEs,numericMeta[,numericIndices.more])
#moduleTraitCor.more <- MEcors.more$bicor
#moduleTraitPvalue.more <- MEcors.more$p


CairoPDF(file=paste0(outputfigs,"4p13c1.GlobalNetworkPlots-part3-",FileBaseName,"-Cairo.pdf"),width=18,height=18)


#These are your numerically coded traits:
colnames(numericMeta)[numericIndices] #choose traits for correlation scatterplots (verboseScatterplot functions below)

#These are your ANOVA sample groups and the number of samples in each
table(Grouping) #alphabetically ordered, you choose the order of groups in the boxplot function by typing them in

## Make changes after checking output on console for the above 2 lines
#CairoPDF(file=paste0(outputfigs,"/GlobalNetPlots(BoxPlots)_",FileBaseName,"-CAIRO.pdf"),width=18,height=11.25)
##pdf(file=paste0(outputfigs,"/GlobalNetPlots(BoxPlots)_",FileBaseName,".pdf"),width=18,height=11.25)

#par(mfrow=c(4,6))
par(mar=c(6.5,6,4.5,1.5))

layout(matrix(c(1,1,2,3,4,4, 5,6,7,8,9,10, 11,12,13,14,15,16, 17,18,19,20,21,22, 23,23,24,24,25,25, 26,27,28,29,30,31, 32,33,34,35,36,37), nrow = 7, ncol = 6, byrow=TRUE),
       heights = c(0.9,0.9,0.9,0.9), # Heights of the four rows
       widths = c(1,1,1,1,1,1)) # Widths of the 5 columns
##sapply(c(1:24),layout.show)

library(beeswarm)
library(gplots)

for (i in 1:(nrow(toplot))) {  # grey already excluded, no -1
  titlecolor<-if(signif(pvec.group,2)[i] <0.05) { "red" } else { "black" }
  boxplot(toplot[i,]~factor(Grouping,names(table(Grouping))),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nGroup 1 var AOV p = ",signif(pvec.group,2)[i]," | age+sex adj AOV p = ",signif(pvec.group.ageSexAdj,2)[i]),xlab=NULL,col.main=titlecolor)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  beeswarm(toplot[i,]~factor(Grouping,names(table(Grouping))),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism

  verboseScatterplot(x=abs(numericMeta[,"sex"]-2),y=toplot[i,],xlab="Sex (1=male)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"sex"],2),", p=",signif(moduleTraitPvalue[i,"sex"],2),"\n"),col.main=if(moduleTraitPvalue[i,"sex"]<0.05) { "red" } else { "black" }, las=1.5)
  verboseScatterplot(x=numericMeta[,"age_at_visit"],y=toplot[i,],xlab="Age",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"age_at_visit"],2),", p=",signif(moduleTraitPvalue[i,"age_at_visit"],2),"\n"),col.main=if(moduleTraitPvalue[i,"age_at_visit"]<0.05) { "red" } else { "black" }, las=1.5)

  titlecolor<-if(signif(pvec.batch,2)[i] <0.05) { "red" } else { "black" }
  boxplot(toplot[i,]~factor(numericMeta$contributor_code,names(table(numericMeta$contributor_code))),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nBatch 1 var AOV p = ",signif(pvec.batch,2)[i]),xlab=NULL,col.main=titlecolor)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  beeswarm(toplot[i,]~factor(numericMeta$contributor_code,names(table(numericMeta$contributor_code))),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism


#  verboseScatterplot(x=numericMeta[,"Age.CT"],y=toplot[i,],xlab="Age of Control Indiv.",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"Age.CT"],2),", p=",signif(moduleTraitPvalue[i,"Age.CT"],2),"\n"),col.main=if(moduleTraitPvalue[i,"Age.CT"]<0.05) { "red" } else { "black" }, las=1.5)

  verboseScatterplot(x=numericMeta[,"height_cm"],y=toplot[i,],xlab="Height (cm)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"height_cm"],2),", p=",signif(moduleTraitPvalue[i,"height_cm"],2),"\n"),col.main=if(moduleTraitPvalue[i,"height_cm"]<0.05) { "red" } else { "black" }, las=1.5)
  verboseScatterplot(x=numericMeta[,"bmi"],y=toplot[i,],xlab="BMI",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"bmi"],2),", p=",signif(moduleTraitPvalue[i,"bmi"],2),"\n"),col.main=if(moduleTraitPvalue[i,"bmi"]<0.05) { "red" } else { "black" }, las=1.5)
  verboseScatterplot(x=numericMeta[,"resting_heart_rate_pulse"],y=toplot[i,],xlab="Resting HR (pulse)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"resting_heart_rate_pulse"],2),", p=",signif(moduleTraitPvalue[i,"resting_heart_rate_pulse"],2),"\n"),col.main=if(moduleTraitPvalue[i,"resting_heart_rate_pulse"]<0.05) { "red" } else { "black" }, las=1.5)
  verboseScatterplot(x=numericMeta[,"systolic_blood_pressure_sitting"],y=toplot[i,],xlab="Systolic BP",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"systolic_blood_pressure_sitting"],2),", p=",signif(moduleTraitPvalue[i,"systolic_blood_pressure_sitting"],2),"\n"),col.main=if(moduleTraitPvalue[i,"systolic_blood_pressure_sitting"]<0.05) { "red" } else { "black" }, las=1.5)
  verboseScatterplot(x=numericMeta[,"diastolic_blood_pressure_sitting"],y=toplot[i,],xlab="Diastolic BP",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"diastolic_blood_pressure_sitting"],2),", p=",signif(moduleTraitPvalue[i,"diastolic_blood_pressure_sitting"],2),"\n"),col.main=if(moduleTraitPvalue[i,"diastolic_blood_pressure_sitting"]<0.05) { "red" } else { "black" }, las=1.5)
#  verboseScatterplot(x=numericMeta[,"mi"],y=toplot[i,],xlab="Myocardial Infarct (MI)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"mi"],2),", p=",signif(moduleTraitPvalue[i,"mi"],2),"\n"),col.main=if(moduleTraitPvalue[i,"mi"]<0.05) { "red" } else { "black" }, las=1.5)
  titlecolor<-if(signif(pvec.mi,2)[i] <0.05) { "red" } else { "black" }
  boxplot(toplot[i,]~factor(numericMeta$mi,names(table(numericMeta$mi))),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nMyocardial Infarct (=1)\n1 var AOV p = ",signif(pvec.mi,2)[i]),xlab=NULL,col.main=titlecolor)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  beeswarm(toplot[i,]~factor(numericMeta$mi,names(table(numericMeta$mi))),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism


#  verboseScatterplot(x=numericMeta[,"APOE4.Dose"],y=toplot[i,],xlab="APOE e4 Dose (0/1/2)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"APOE4.Dose"],2),", p=",signif(moduleTraitPvalue[i,"APOE4.Dose"],2),"\n"),col.main=if(moduleTraitPvalue[i,"APOE4.Dose"]<0.05) { "red" } else { "black" }, las=1.5)
  titlecolor<-if(signif(pvec.APOE4.Dose,2)[i] <0.05) { "red" } else { "black" }
  boxplot(toplot[i,]~factor(numericMeta$APOE4.Dose,names(table(numericMeta$APOE4.Dose))),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nAPOE e4 Dose\n1 var AOV p = ",signif(pvec.APOE4.Dose,2)[i]),xlab=NULL,col.main=titlecolor)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  beeswarm(toplot[i,]~factor(numericMeta$APOE4.Dose,names(table(numericMeta$APOE4.Dose))),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism

#  verboseScatterplot(x=numericMeta[,"tia"],y=toplot[i,],xlab="Transient Ischemic Attack (TIA)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"tia"],2),", p=",signif(moduleTraitPvalue[i,"tia"],2),"\n"),col.main=if(moduleTraitPvalue[i,"tia"]<0.05) { "red" } else { "black" }, las=1.5)
  titlecolor<-if(signif(pvec.tia,2)[i] <0.05) { "red" } else { "black" }
  boxplot(toplot[i,]~factor(numericMeta$tia,names(table(numericMeta$tia))),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nTransient Ischemic Attack\n1 var AOV p = ",signif(pvec.tia,2)[i]),xlab=NULL,col.main=titlecolor)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  beeswarm(toplot[i,]~factor(numericMeta$tia,names(table(numericMeta$tia))),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism

#  verboseScatterplot(x=numericMeta[,"tbi"],y=toplot[i,],xlab="Traumatic Brain Injury (TBI)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"tbi"],2),", p=",signif(moduleTraitPvalue[i,"tbi"],2),"\n"),col.main=if(moduleTraitPvalue[i,"tbi"]<0.05) { "red" } else { "black" }, las=1.5)
  titlecolor<-if(signif(pvec.tbi,2)[i] <0.05) { "red" } else { "black" }
  boxplot(toplot[i,]~factor(numericMeta$tbi,names(table(numericMeta$tbi))),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nTraumatic Br Injury\n1 var AOV p = ",signif(pvec.tbi,2)[i]),xlab=NULL,col.main=titlecolor)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  beeswarm(toplot[i,]~factor(numericMeta$tbi,names(table(numericMeta$tbi))),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism

  verboseScatterplot(x=numericMeta[,"MMSE"],y=toplot[i,],xlab="MMSE Cog. Score",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"MMSE"],2),", p=",signif(moduleTraitPvalue[i,"MMSE"],2),"\n"),col.main=if(moduleTraitPvalue[i,"MMSE"]<0.05) { "red" } else { "black" }, las=1.5)
  verboseScatterplot(x=numericMeta[,"MoCA"],y=toplot[i,],xlab="MoCA Cog. Score",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"MoCA"],2),", p=",signif(moduleTraitPvalue[i,"MoCA"],2),"\n"),col.main=if(moduleTraitPvalue[i,"MoCA"]<0.05) { "red" } else { "black" }, las=1.5)
  verboseScatterplot(x=numericMeta[,"sequential_visit_number"],y=toplot[i,],xlab="Sequential Visit #",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"sequential_visit_number"],2),", p=",signif(moduleTraitPvalue[i,"sequential_visit_number"],2),"\n"),col.main=if(moduleTraitPvalue[i,"sequential_visit_number"]<0.05) { "red" } else { "black" }, las=1.5)

# UNUSED numericMeta columns for these correlation plots:
#"raceAA","years_of_education",,"weight_kg","alcohol_hx", "recruited_control","ad","ftd","pd","als","mci_sci"
#"depression","anxiety","visit"

## No non-NA values in the 16677 for rare mutation statuses.
#  titlecolor<-if(signif(pvec.C9Orf72,2)[i] <0.05) { "red" } else { "black" }
#  boxplot(toplot[i,]~factor(numericMeta$C9Orf72,c(0,1)),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nC9Orf72 Mutated (=1)\n1 var AOV p = ",signif(pvec.C9Orf72,2)[i]),xlab=NULL,col.main=titlecolor)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
#  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
#  beeswarm(toplot[i,]~factor(numericMeta$C9Orf72,c(0,1)),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism
#
#  titlecolor<-if(signif(pvec.GRN,2)[i] <0.05) { "red" } else { "black" }
#  boxplot(toplot[i,]~factor(numericMeta$GRN,c(0,1)),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nGRN Mutated (=1)\n1 var AOV p = ",signif(pvec.GRN,2)[i]),xlab=NULL,col.main=titlecolor)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
#  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
#  beeswarm(toplot[i,]~factor(numericMeta$GRN,c(0,1)),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism
#
#  titlecolor<-if(signif(pvec.MAPT,2)[i] <0.05) { "red" } else { "black" }
#  boxplot(toplot[i,]~factor(numericMeta$MAPT,c(0,1)),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nMAPT Mutated (=1)\n1 var AOV p = ",signif(pvec.MAPT,2)[i]),xlab=NULL,col.main=titlecolor)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
#  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
#  beeswarm(toplot[i,]~factor(numericMeta$MAPT,c(0,1)),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism

# replacement plots 17, 18, 19: 6 APOE genotypes, in all or AD+MCI only, CT only
  titlecolor<-if(signif(pvec.APOEgeno.allSamp,2)[i] <0.05) { "red" } else { "black" }
  boxplot(toplot[i,]~factor(numericMeta$APOE.mapped.predicted,names(table(numericMeta$APOE.mapped.predicted))),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nAPOE genotype (All samp)\n1 var AOV p = ",signif(pvec.APOEgeno.allSamp,2)[i]),xlab=NULL,col.main=titlecolor, las=2)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  beeswarm(toplot[i,]~factor(numericMeta$APOE.mapped.predicted,names(table(numericMeta$APOE.mapped.predicted))),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism

  titlecolor<-if(signif(pvec.APOEgeno.CTonly,2)[i] <0.05) { "red" } else { "black" }
  boxplot(toplot[i,which(Grouping=="CT")]~factor(numericMeta$APOE.mapped.predicted[which(Grouping=="CT")],names(table(numericMeta$APOE.mapped.predicted))),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nAPOE genotype (CT only)\n1 var AOV p = ",signif(pvec.APOEgeno.CTonly,2)[i]),xlab=NULL,col.main=titlecolor)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  beeswarm(toplot[i,which(Grouping=="CT")]~factor(numericMeta$APOE.mapped.predicted[which(Grouping=="CT")],names(table(numericMeta$APOE.mapped.predicted))),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism

  titlecolor<-if(signif(pvec.APOEgeno.ADorMCIonly,2)[i] <0.05) { "red" } else { "black" }
  boxplot(toplot[i,which(Grouping=="AD" | Grouping=="MCI")]~factor(numericMeta$APOE.mapped.predicted[which(Grouping=="AD" | Grouping=="MCI")],names(table(numericMeta$APOE.mapped.predicted))),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nAPOE genotype (AD+MCI only)\n1 var AOV p = ",signif(pvec.APOEgeno.ADorMCIonly,2)[i]),xlab=NULL,col.main=titlecolor)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  beeswarm(toplot[i,which(Grouping=="AD" | Grouping=="MCI")]~factor(numericMeta$APOE.mapped.predicted[which(Grouping=="AD" | Grouping=="MCI")],names(table(numericMeta$APOE.mapped.predicted))),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism


  verboseScatterplot(x=numericMeta[,"cdr"],y=toplot[i,],xlab="Clinical Dem Rating (CDR)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"cdr"],2),", p=",signif(moduleTraitPvalue[i,"cdr"],2),"\n"),col.main=if(moduleTraitPvalue[i,"cdr"]<0.05) { "red" } else { "black" }, las=1.5)
  verboseScatterplot(x=numericMeta[,"total_years_smoked"],y=toplot[i,],xlab="Total Years Smoked",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"total_years_smoked"],2),", p=",signif(moduleTraitPvalue[i,"total_years_smoked"],2),"\n"),col.main=if(moduleTraitPvalue[i,"total_years_smoked"]<0.05) { "red" } else { "black" }, las=1.5)

  verboseScatterplot(x=numericMeta[,"RegrBloodPreanalyticFactor.HNRNPA2B1"],y=toplot[i,],xlab="Preanalytical Factor\nProxy 1 (HNRNPA2B1)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"RegrBloodPreanalyticFactor.HNRNPA2B1"],2),", p=",signif(moduleTraitPvalue[i,"RegrBloodPreanalyticFactor.HNRNPA2B1"],2),"\n"),col.main=if(moduleTraitPvalue[i,"RegrBloodPreanalyticFactor.HNRNPA2B1"]<0.05) { "red" } else { "black" }, las=1.5)


  # New plot1 - all groups including 3 cohorts best mapping and CT imputed
  titlecolor<-if(signif(pvec.group.withCTimputed,2)[i] <0.05) { "red" } else { "black" }
  boxplot(toplot[i,]~factor(numericMeta$Group.withCTimputed,names(table(numericMeta$Group.withCTimputed))),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nGroup+RM.BH.UDSmapped+CTimputed 1 var AOV p = ",signif(pvec.group.withCTimputed,2)[i],"\nage+sex adj AOV p = ",signif(pvec.group.withCTimputed.ageSexAdj,2)[1]),xlab=NULL,col.main=titlecolor,las=2)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  beeswarm(toplot[i,]~factor(numericMeta$Group.withCTimputed,names(table(numericMeta$Group.withCTimputed))),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism

  # New plot - all groups including 3 cohorts best mapping
  titlecolor<-if(signif(pvec.group.3mappedCohorts,2)[i] <0.05) { "red" } else { "black" }
  boxplot(toplot[i,]~factor(numericMeta$Group.3mappedCohorts,names(table(numericMeta$Group.3mappedCohorts))),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nGroup with 3 mapped cohorts 1 var AOV p = ",signif(pvec.group.3mappedCohorts,2)[i],"\nage+sex adj AOV p = ",signif(pvec.group.3mappedCohorts.ageSexAdj,2)[i]),xlab=NULL,col.main=titlecolor,las=2)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  beeswarm(toplot[i,]~factor(numericMeta$Group.3mappedCohorts,names(table(numericMeta$Group.3mappedCohorts))),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism

  # New plot2 - only Group.pathCog from BH and UDS
  titlecolor<-if(signif(pvec.group.pathCog.mapped,2)[i] <0.05) { "red" } else { "black" }
  boxplot(toplot[i,]~factor(numericMeta$Group.pathCog.mapped,names(table(numericMeta$Group.pathCog.mapped))),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nGroup.pathCog (BH+UDS only) 1 var AOV p = ",signif(pvec.group.pathCog.mapped,2)[i],"\nage+sex adj AOV p = ",signif(pvec.group.pathCog.mapped.ageSexAdj,2)[1]),xlab=NULL,col.main=titlecolor,las=2)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  beeswarm(toplot[i,]~factor(numericMeta$Group.pathCog.mapped,names(table(numericMeta$Group.pathCog.mapped))),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism


# new plots
  verboseScatterplot(x=numericMeta[,"Lilly.BH.blood.pTau217"],y=toplot[i,],xlab="Lilly pTau217 (BH)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"Lilly.BH.blood.pTau217"],2),", p=",signif(moduleTraitPvalue[i,"Lilly.BH.blood.pTau217"],2),"\n"),col.main=if(moduleTraitPvalue[i,"Lilly.BH.blood.pTau217"]<0.05) { "red" } else { "black" }, las=1.5)
  verboseScatterplot(x=numericMeta[,"UDS.blood.pTau217"],y=toplot[i,],xlab="pTau217 (UDS)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"UDS.blood.pTau217"],2),", p=",signif(moduleTraitPvalue[i,"UDS.blood.pTau217"],2),"\n"),col.main=if(moduleTraitPvalue[i,"UDS.blood.pTau217"]<0.05) { "red" } else { "black" }, las=1.5)
#  verboseScatterplot(x=numericMeta[,"AmyloidPositivity.withRM"],y=toplot[i,],xlab="Age",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"AmyloidPositivity.withRM"],2),", p=",signif(moduleTraitPvalue[i,"AmyloidPositivity.withRM"],2),"\n"),col.main=if(moduleTraitPvalue[i,"AmyloidPositivity.withRM"]<0.05) { "red" } else { "black" }, las=1.5)
  titlecolor<-if(signif(pvec.AmyloidPositivity.withRM,2)[i] <0.05) { "red" } else { "black" }
  boxplot(toplot[i,]~factor(numericMeta$AmyloidPositivity.withRM,names(table(numericMeta$AmyloidPositivity.withRM))),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nAmy. Positivity (3 cohorts)\n1 var AOV p = ",signif(pvec.AmyloidPositivity.withRM,2)[i]),xlab=NULL,col.main=titlecolor)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  beeswarm(toplot[i,]~factor(numericMeta$AmyloidPositivity.withRM,names(table(numericMeta$AmyloidPositivity.withRM))),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism

  titlecolor<-if(signif(pvec.hypertension,2)[i] <0.05) { "red" } else { "black" }
  boxplot(toplot[i,]~factor(numericMeta$hypertension,names(table(numericMeta$hypertension))),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nHypertension\n1 var AOV p = ",signif(pvec.hypertension,2)[i]),xlab=NULL,col.main=titlecolor)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  beeswarm(toplot[i,]~factor(numericMeta$hypertension,names(table(numericMeta$hypertension))),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism

  titlecolor<-if(signif(pvec.hyperlipidaemia,2)[i] <0.05) { "red" } else { "black" }
  boxplot(toplot[i,]~factor(numericMeta$hyperlipidaemia,names(table(numericMeta$hyperlipidaemia))),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nHyperlipidemia\n1 var AOV p = ",signif(pvec.hyperlipidaemia,2)[i]),xlab=NULL,col.main=titlecolor)  #rotate x labs: ,las=2  #no outliers: ,outline=FALSE)
  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  beeswarm(toplot[i,]~factor(numericMeta$hyperlipidaemia,names(table(numericMeta$hyperlipidaemia))),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism

  verboseScatterplot(x=numericMeta[,"RegrBloodPreanalyticFactor.HBZ"],y=toplot[i,],xlab="Preanalytical Factor\nProxy 2 (HBZ)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"RegrBloodPreanalyticFactor.HBZ"],2),", p=",signif(moduleTraitPvalue[i,"RegrBloodPreanalyticFactor.HBZ"],2),"\n"),col.main=if(moduleTraitPvalue[i,"RegrBloodPreanalyticFactor.HBZ"]<0.05) { "red" } else { "black" }, las=1.5)

  verboseScatterplot(x=numericMeta[,"TimeToSpin"],y=toplot[i,],xlab="Time to Spin\n(SomaSignal)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"TimeToSpin"],2),", p=",signif(moduleTraitPvalue[i,"TimeToSpin"],2),"\n"),col.main=if(moduleTraitPvalue[i,"TimeToSpin"]<0.05) { "red" } else { "black" }, las=1.5)
  verboseScatterplot(x=numericMeta[,"TimeToDecant"],y=toplot[i,],xlab="Time to Decant\n(SomaSignal)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"TimeToDecant"],2),", p=",signif(moduleTraitPvalue[i,"TimeToDecant"],2),"\n"),col.main=if(moduleTraitPvalue[i,"TimeToDecant"]<0.05) { "red" } else { "black" }, las=1.5)
  verboseScatterplot(x=numericMeta[,"TimeToFreeze"],y=toplot[i,],xlab="Time to Freeze\n(SomaSignal)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"TimeToFreeze"],2),", p=",signif(moduleTraitPvalue[i,"TimeToFreeze"],2),"\n"),col.main=if(moduleTraitPvalue[i,"TimeToFreeze"]<0.05) { "red" } else { "black" }, las=1.5)
  verboseScatterplot(x=numericMeta[,"FedFastedTime"],y=toplot[i,],xlab="Fed-fasted Time\n(SomaSignal)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"FedFastedTime"],2),", p=",signif(moduleTraitPvalue[i,"FedFastedTime"],2),"\n"),col.main=if(moduleTraitPvalue[i,"FedFastedTime"]<0.05) { "red" } else { "black" }, las=1.5)
  verboseScatterplot(x=numericMeta[,"FreezeThawCycles"],y=toplot[i,],xlab="Freeze-thaw Cycles\n(SomaSignal)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"FreezeThawCycles"],2),", p=",signif(moduleTraitPvalue[i,"FreezeThawCycles"],2),"\n"),col.main=if(moduleTraitPvalue[i,"FreezeThawCycles"]<0.05) { "red" } else { "black" }, las=1.5)

  verboseScatterplot(x=numericMeta[,"RegrBloodPreanalyticFactor.HNRNPA2B1"],y=toplot[i,],xlab="Preanalytical Factor\nProxy 1 (HNRNPA2B1)",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"RegrBloodPreanalyticFactor.HNRNPA2B1"],2),", p=",signif(moduleTraitPvalue[i,"RegrBloodPreanalyticFactor.HNRNPA2B1"],2),"\n"),col.main=if(moduleTraitPvalue[i,"RegrBloodPreanalyticFactor.HNRNPA2B1"]<0.05) { "red" } else { "black" }, las=1.5)

}


##outputs sample-by-sample eigenprotein barplots (not useful for large number of samples)
#while(!par('page')) plot.new()
#for (i in 1:nrow(toplot)) {
# barplot(height=rev(toplot[i,]),width=5,col=colnames(MEs)[i],xlab=paste(colnames(MEs)[i]," Eigenprotein Relative Expression"),main=rownames(toplot)[i],ylab=NULL,las=2,space=0.4,horiz=TRUE) #las=2 for rotated 90° X-axis labels  main=rownames(toplot)[i]
## text(bargr,par("usr")[3] - 0.025, srt=45, adj =1, labels= c(colnames(toplot)),xpd=TRUE,font=2) # bargr <- barplot(... above; gives rotated 45° x-axis labels but overwrites on top of existing ones
#}

dev.off() #finishes and closes writing of globalNetworkPlots PDF


########################################
#Write Module Membership/kME table
orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
kMEtableSortVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) if(!x[1]=="grey") { paste0(paste(orderedModulesWithGrey[match(x[1],orderedModulesWithGrey[,2]),],collapse=" "),"|",round(as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]),4)) } else { paste0("grey|AllKmeAvg:",round(mean(as.numeric(x[-1],na.rm=TRUE)),4)) } )
kMEtable=cbind(c(1:nrow(cleanDat)),rownames(cleanDat),net$colors,kMEdat,kMEtableSortVector)[order(kMEtableSortVector,decreasing=TRUE),]
write.table(kMEtable,file=paste0(outputtabs,"/4p13c1.ModuleAssignments-",FileBaseName,".txt"),sep="\t",row.names=FALSE)
#(load above file in excel and apply green-yellow-red conditional formatting heatmap to the columns with kME values); then save as excel.

## saved image of R session
#save.image(paste0("./4p.saved.image.",FileBaseName,".Rdata"))  #overwrites
#load("4p.saved.image.PLASMA_TAMPORmode4.209iter.(7334x27783)_3414sampleMedianRows.WGCNA.mms10.ds4.pwr6.Rdata")


## SOMA 19 Gene Lists - FET
#####################################################################################

#rootdir="f:/OneDrive - Emory/SOMAplasmaMultibatch/FETs.Figure5/Up.Down.All.SeparatePages/"
#setwd(rootdir)


# ------------------------------------------------------------------------
# ANNOTATION: Export the final module assignment/kME table.
# ------------------------------------------------------------------------
#load("f:/OneDrive - Emory/SOMAplasmaMultibatch/Nets/BH+RM+UDS+EMT.mode1TAMPOR/4b.4cohort.TAMPORmode1.27mods_regrHNRNPA2B1+HBZ_forORA.Rdata")
##contains:
#net<-net4cohort.27mod
#cleanDat<-cleanDat4cohort.27mod
#numericMeta<-numericMeta4cohort.27mod

#setwd("./FET/")


######################################
## Cell Type Enrichment
#source("./geneListFET.R")

#geneListFET(FileBaseName="5.4CohortBloodNet_FET_to_SelectedGeneLists",
#            heatmapTitle="4 Cohort Blood Plasma Network 27 Module Overlaps with Selected Gene Lists",
#            modulesInMemory=TRUE,categorySpeciesCode="hsapiens",  # use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
#            refDataFiles=c("FETinputs-MainFig5.csv"),
#            speciesCode=c("hsapiens"),refDataDescription="mainFET_GeneLists",  # file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?
#            paletteColors=c("Spectral"),
#            verticalCompression=1)  # colors from RColorBrewer::display.brewer.all()


### Repeat with Descriptive labels for 17 modules in 15 split tracks
#MElabels<-read.csv("./plasma_21cohort_17moduleOntologies_unregr+color.csv",header=TRUE)
#
#we do not have labels yet, just make placeholder table of M#s
nModules<-length(table(net$colors))-1
MElabels<-data.frame( Mnum=paste0("M",c(1:nModules)),Module.color.Description=paste0("M",c(1:nModules)," unannot."),Module.color=paste0("M",c(1:nModules)," ",labels2colors(c(1:nModules))),Description=rep("unannot.",nModules),Color=labels2colors(c(1:nModules)) )

#
#source("./geneListFET_customLabels.R")
#geneListFET(FileBaseName="Figure5.21CohortBloodNet_rectangularFETs_15pages_separateFDR",
#            heatmapTitle="4 Cohort Blood Plasma Network 27 Module Overlaps with Selected Gene Lists",
#            modulesInMemory=TRUE,categorySpeciesCode="hsapiens",  # use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
#            refDataFiles=c("1.Amyloidosis.csv","2.CogFn.FDRsig.csv",
#            "3.CogFn.AmyAdj.pSig.csv","4.ROSMAP.Cogfn.NonPathAssoc.csv",
#            "5.Amyloidosis.APOEdep.csv","6.Tau.NFT.Tangles.Union.csv","7a.Tau.NFT.csv","7b.Tau.Tangles.csv",
#            "8.DLBdx.csv","9.CI.mct.csv","10.CI.gct.csv","11.Arteriosclerosis.csv","12.CAA.4gp.csv","13.CVDA.4gp2.csv","14.TDP.st4.csv","15.CoxHR-Intersection.3way.csv"),
#            speciesCode=c("hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens"),refDataDescription="mainFET_GeneLists",  # file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?
#            paletteColors=c("Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral"),
#            verticalCompression=5)  # colors from RColorBrewer::display.brewer.all()
#
#
## repeat for unadj p value scale - all plots
#geneListFET(FileBaseName="Figure5.21CohortBloodNet_rectangularFETs_15pages_separate(NO_FDR)",
#            heatmapTitle="4 Cohort Blood Plasma Network 27 Module Overlaps with Selected Gene Lists",
#            modulesInMemory=TRUE,categorySpeciesCode="hsapiens",  # use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
#            heatmapScale="p.unadj",
#            refDataFiles=c("1.Amyloidosis.csv","2.CogFn.FDRsig.csv",
#            "3.CogFn.AmyAdj.pSig.csv","4.ROSMAP.Cogfn.NonPathAssoc.csv",
#            "5.Amyloidosis.APOEdep.csv","6.Tau.NFT.Tangles.Union.csv","7a.Tau.NFT.csv","7b.Tau.Tangles.csv",
#            "8.DLBdx.csv","9.CI.mct.csv","10.CI.gct.csv","11.Arteriosclerosis.csv","12.CAA.4gp.csv","13.CVDA.4gp2.csv","14.TDP.st4.csv","15.CoxHR-Intersection.3way.csv"),
#            speciesCode=c("hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens","hsapiens"),refDataDescription="mainFET_GeneLists",  # file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?
#            paletteColors=c("Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral","Spectral"),
#            verticalCompression=5)  # colors from RColorBrewer::display.brewer.all()


# repeat with fixed scale -logFDR 0-20, and 15th plot uses nominal -log10(p)
source("geneListFET_customLabels-fixedScale+15thPlotNominalP.R")
geneListFET(FileBaseName="4p13c1.Figure5.MultiCohortBloodNet_rectangularFETs_15pages_separate(oneScale0-20)_lastPlotNOMINALp_v02-20pp",
            heatmapTitle=paste0("HDS Blood Plasma Network ",nModules," Module Overlaps with Selected Gene Lists"),
            modulesInMemory=TRUE,categorySpeciesCode="hsapiens",  # use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
#            heatmapScale="p.unadj",
            refDataFiles=c("1.Amyloidosis.csv","2.CogFn.FDRsig.csv",
            "3.CogFn.AmyAdj.pSig.csv","4.ROSMAP.Cogfn.NonPathAssoc.csv",
            "5.Amyloidosis.APOEdep.csv","6.Tau.NFT.Tangles.Union.csv","7a.Tau.NFT.csv","7b.Tau.Tangles.csv",
            "8.DLBdx.csv","9.CI.mct.csv","10.CI.gct.csv","11.Arteriosclerosis.csv","12.CAA.4gp.csv","13.CVDA.4gp2.csv","14.TDP.st4.csv","15.CoxHR-Intersection.3way.csv",
            "16.PlaqDplusN.csv","17.PlaqD.csv","18.PlaqN.csv","19.Brain-linked.csv"),
            speciesCode=rep("hsapiens",20),refDataDescription="mainFET_GeneLists",  # file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?
            paletteColors=rep("Spectral",20),
            verticalCompression=1.1)  # colors from RColorBrewer::display.brewer.all()
            #changed page height from 15 to 3.5, verticalCompression from 5 to 1.1


#setwd("../")

cleanDat.full<-cleanDat
numericMeta.full<-numericMeta

length(unique(numericMeta.full$person_id))
#17484 (in 22392)  # was 13783 (in 16677)

numericMeta<-numericMeta.full

# ------------------------------------------------------------------------
# ANNOTATION: Run module-level enrichment and custom-labeled plots.
# ------------------------------------------------------------------------
numericMeta$sample_id<-rownames(numericMeta)
library(dplyr)
numericMeta <- numericMeta %>%
  group_by(person_id) %>%
  filter(sequential_visit_number == max(sequential_visit_number)) %>%
  slice_head(n = 1)
numericMeta<-as.data.frame(numericMeta)
rownames(numericMeta)<-numericMeta$sample_id

dim(numericMeta)
# 17484   107
cleanDat<-cleanDat.full[,match(rownames(numericMeta),colnames(cleanDat.full))]


######################################
## ANOVA + Volcanoes + DEXstacked Barplots
source("parANOVA.dex.fallback7.25.R")

parallelThreads=8
outFilePrefix="4p13c1"
outFileSuffix="SixDxVolcs.CTimputed"

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

# ------------------------------------------------------------------------
# ANNOTATION: Run parANOVA-based phenotype contrasts used for module/assay
# volcano outputs.
# ------------------------------------------------------------------------
# None need flipping (more extreme expected phenotype is first (numerator) in each case)
flip=c(3:8,15:17)
sameScale=FALSE #TRUE
symbolsOnly=TRUE
highlightGeneProducts=c("HNRNPA2B1","HBZ","SPC25","CPLX2","PTN","MAPT","OMG","GDF15","NEFL","NRGN","CTHRC1","HTRA1","TTR","LRP1","NTN1","SFRP1","GPNMB","PAFAH1B3","CRIP1")
labelHighlighted=TRUE      # if true, highlighted spots get text labels with their rownames from ANOVAout
labelTop=5
plotVolc()                 # runs on ANOVAout as input (need not be specified).

DEXpercentStacked()        # runs on prior function outputs as input; writes stacked bar plot(s) to PDF.

#outputs moved to folder 4b.DiffEx_Ctx


######################################
## Alternative Correlation (to linear trait), stats table with volcanoes

# These parameters are specific to trait correlation statistics generation; traits are provided as columns of the data frame stored in the provided example RData as the variable numericMeta.
cor.traits=c("age_at_visit","MMSE", "Lilly.BH.blood.pTau217","UDS.blood.pTau217", "RegrBloodPreanalyticFactor.HNRNPA2B1","RegrBloodPreanalyticFactor.HBZ",
             "TimeToSpin","TimeToDecant","TimeToFreeze","FedFastedTime","FreezeThawCycles")                 # Molecular and quantitative Traits to correlate to in numericMeta columns (colnames)
filter.trait="Group.withCTimputed"             # Trait on which to subset case samples
filter.trait.subsets=c("ALL","AD","CT") # Subsets of case samples will be used for correlation to the cell type proportion estimates
                                    # (4 separate cor.traits x 4 sample subsets = 16 total p and R value columns to generate)
corFn="bicor"                       #'bicor'; other options are 'kendall', 'spearman', ...anything else will cause Pearson (cor) to be used


#source("./parANOVA.dex.R")
CORout <- trait.corStat()                      # runs on cleanDat and Grouping variables as required input.
# Correlation p + R table calculations complete. If you want to use the table with plotVolc(), set the variable corVolc=TRUE and use variable CORout to store the table generated.

dim(CORout)
#[1] 7333   69


outFileSuffix="BICORvolcs"
corVolc=TRUE        # changes the behavior of plotVolc, DEXpercentStacked, and GOparallel functions later in the pipeline, to use CORout
useNETcolors=TRUE
sameScale=TRUE
#highlightGeneProducts=rownames(CORout)[which(CORout$Filter.ALLsamples.FavoriteCorrStat.Sig)]  # Specifies which spots should be large
#labelHighlighted=TRUE

#Change values less than 1e-50 to 1e-50, so same scale volcanoes will not be stretched to -log10(p) of 200!
plateauMinP=1e-50
CORout[,3:(which(grepl("bicor ",colnames(CORout)))[1]-1)]<-apply(CORout[,3:(which(grepl("bicor ",colnames(CORout)))[1]-1)],2,function(x) { x1<-x; x1[x1<plateauMinP]<-plateauMinP; x1; })
plotVolc()          # Plots PDFs and HTMLs

#DEXpercentStacked(CORout)


## not run
#source("./GOparallel-FET.R")  #Available from the repository file https://github.com/edammer/GOparallel/blob/main/GOparallel-FET.R
#ANOVAgroups=TRUE
#parallelThreads=8
#outFilename="4p.CORsig_Volcs.GO"
#GOparallel(CORout)


# Going forward, do not use correlation statistics, use ANOVA groupwise stats.
corVolc=FALSE

## This section's outputs can be moved to subfolder: "corStats.volcanoes"


cleanDat<-cleanDat.full
numericMeta<-numericMeta.full

rm(numericMeta.full)
rm(cleanDat.full)

############################################################################################
# iGRAPHs (Multiple Toggle Options, e.g. BioGRID interactome overlap) // CONNECTIVITY PLOT #
############################################################################################
## Configuration section  ##
## - check these settings ##
############################################################################################
PPIedges=FALSE #*** matrix too large for adjacency calc in 12 hr       # TRUE will be somewhat slower, and if few parallelThreads specified below...
myHumanBioGrid.tsvFile = "BIOGRID-ORGANISM-Homo_sapiens-4.4.235.mitab.HUMANsimple.txt"  # "nonexistent.file"  # ONLY NEEDED IF PPIedges=TRUE
parallelThreads=8      # needed if PPIedges=TRUE; set to # of threads on your computer
#Symbol to add to iGraphs that don't include it, if checking for protein-protein interaction edges (add 4, one for each corner; can edit output later to remove unwanted ones.)
symbols2Add=c()     # interactions to these extra nodes (not necessarily in the modules) will be drawn for every module.
                    # There must be 4 symbols exactly. These nodes appear in the 4 corners of each plot.
#c()          # if you do not want to add these supplemental nodes to corners! (or do not set variable)
showAllPPIs=FALSE       # if there are 4 symbols2Add, unless this flag is true, only PPIs to the 4 corner nodes will be drawn.
#####################
GOIlist<- c("HNRNPA2B1","HBZ","SPC25","CPLX2","PTN","MAPT","OMG","GDF15","NEFL","NRGN","CTHRC1","HTRA1","TTR","LRP1","NTN1","SFRP1","GPNMB","PAFAH1B3","CRIP1")  # unique(unlist(as.list(read.csv(file="RNAbindingProtein.FET.lists.csv",header=TRUE))))
# Genes in modules to highlight -- could be all RNA binding proteins, or a list of MAGMA-significant genes for a disease...
# ...they get highlighted yellow, or cyan if the module is yellow.
# The above file of ontology-style organized RNA-binding and co-aggregating proteins is from Guo Q et al, 2021: https://www.frontiersin.org/articles/10.3389/fnmol.2021.623659/full
vertexsize=16            # 8 for regular, 16 for large balls
netSpecies="human"       # current option "mouse" will convert bioGRID to mouse symbols before drawing PPI edges.
showTOMhubs=FALSE        # calculates TOM from cleanDat, so that hubs can be found and shown on a second graph for each module.
                         # the second graph uses Fruchterman-Reingold layout and only nodes with an edge TOM value calculated by WGCNA
                         # that is equal to or greater than the Nth top such edge. Default N=150. See full parameter list to change.
outFilePrefix="4p13c1"
outFileSuffix="PlasmaModules.22392samp"

recalcMEs=FALSE   # save 10 minutes here!
############################################################################################
power  # will be used for adjacency/TOM recalculation
# 6

source("buildIgraphs.R")
buildIgraphs()


save.image("saved.image-genotype_prediction_finalized.RData")   #"4p13c1.saved.image.PLASMA_(7333x22392)_WGCNA.mms10.ds4.pwr6.RData")  #overwrites
#load("saved.image-genotype_prediction_finalized.RData")

#note: loads in 5 min, 27GB footprint in RAM


##############################################################################
## GOparallel - v1.2 - with hitLists included in .csv outputs.

modulesInMemory=TRUE            # uses cleanDat, net, and kMEdat from pipeline already in memory
ANOVAgroups=FALSE               # if true, modulesInMemory ignored. Volcano pipeline code should already have been run!
outFilename="4p13c1.PlasmaModules.GO"

GMTdatabaseFile="z:/EBD/Human_GO_AllPathways_noPFOCR_with_GO_iea_March_01_2025_symbol.gmt"
GO.OBOfile="z:/EBD/go.obo"
source("GOparallel-FET.R")
GOparallel()


#modulesInMemory=FALSE            # uses cleanDat, net, and kMEdat from pipeline already in memory
#ANOVAgroups=TRUE                 # if true, modulesInMemory ignored. Volcano pipeline code should already have been run!
#testIndexMasterList=c(3,4,5,6,7)
#outFilename="4b.Plasma_DEX_lists_byGroup_TAMPORnoRegr.GO"   #added _byGroup to better specify test in folder name
#

# ------------------------------------------------------------------------
# ANNOTATION: Build module graph visualizations for selected proteins/genes
# of interest.
# ------------------------------------------------------------------------
#GOparallel()


setwd(rootdir)

saveRDS(MEs,"22392_sample_net_12MEs.RDS")  #Z:/EBD/grid/FullSet_APOEfinalization2ndRegr/22392_sample_et_12MEs.RDS
#load("z:/EBD/APOE.homozygote_matrices-forShijia.RData") # contains cleanDat.Shijia, numericMeta.Shijia

MEs.3177<-MEs[match(rownames(numericMeta.3177),rownames(MEs)),]
save(MEs.3177,file="Z:/EBD/APOE.homozygote_MEs.3177.RData")

MEs.4199<-MEs[match(rownames(numericMeta.4199),rownames(MEs)),]
save(MEs.4199,file="Z:/EBD/APOE.homozygote_MEs.4199.RData")
