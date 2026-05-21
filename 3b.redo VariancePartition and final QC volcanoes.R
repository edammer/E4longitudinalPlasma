##############################################################################
# Pipeline annotation header: 3b.redo VariancePartition and final QC volcanoes.R
# Manuscript code section(s): 3 QC
#
# Purpose:
# Re-run final variance partitioning and manuscript-facing QC volcano/t-SNE
# visualizations on the APOE-aware regressed plasma matrix.
#
# Principal inputs:
#   - saved.image-genotype_prediction_finalized.RData
#   - parANOVA.dex.fallback7.25.R
#   - ANOVA CSV outputs from prior parANOVA runs
#
# Principal outputs:
#   - 4p13c1...VariancePartition...pdf
#   - EffectSize_Plots...pdf
#   - final t-SNE QC PDF
#
# Step overview:
#   1. Define variance-partition model terms for age, sex, contributor_Fsplit,
#      and APOE e4 carrier status.
#   2. Fit and plot variancePartition results for the final regression matrix.
#   3. Generate volcano plots for APOE e4 carrier, sex, AD-vs-control, and age
#      associations.
#   4. Re-render t-SNE plots colored by site, age, and APOE e4 dose for final
#      QC figures.
#
# Notes added during manuscript-code alignment:
#   - Original executable code was preserved except for whitespace/encoding cleanup
#     and the explicitly marked non-run archive described below when applicable.
#   - Comments inserted by this pass are prefixed with "ANNOTATION:".
##############################################################################

# ------------------------------------------------------------------------
# ANNOTATION: Load the finalized harmonized matrix and regression metadata
# for final QC visualization.
# ------------------------------------------------------------------------
##################################################
# 4. Plasma network exploration  - return to VP final, QC Volcanoes (bicor for age, sex, AD/CT single)


rootdir="z:/EBD/grid/4p13b3forAPOEpredict+2ndRegrAgain/"
#rootdir="c:/Users/edammer_4ceb6ff/Downloads/"
setwd(rootdir)

load("saved.image-genotype_prediction_finalized.RData") # now includes all QC of unreg and c1 regressed 2x data, and cv_folds


## 4p13c1. Variance Partition regressed (2PAV regression intrasite)+Site regressed; protect Age+Sex+E4 carrier binary status (no NA, NA NOT imputed) (QC)

regvars.vp<-data.frame(regvars.c1)
regvars.vp$Sex<-factor(abs(regvars.vp$sex -2))

# ------------------------------------------------------------------------
# ANNOTATION: Construct the variancePartition model using biological and
# technical terms retained for final QC.
# ------------------------------------------------------------------------
regvars.vp$Age<-as.numeric(regvars.vp$age_at_visit)
regvars.vp$contributor_Fsplit<-factor(regvars.vp$contributor_Fsplit)
regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HNRNPA2B1)
regvars.vp$RegrBloodPreanalyticFactor.HBZ<-as.numeric(regvars.vp$RegrBloodPreanalyticFactor.HBZ)
#regvars.vp$APOE.E4carrier.Proxy.LRRN1<-as.numeric(cleanDat.4p13b2["LRRN1|Q6UXK5^SL025922@seq.11293.14",])
regvars.vp$APOE.E4carrier<-factor(regvars.vp$APOE4.carrier)

# too many missing values:
##form <- ~ age_at_visit+(1|Sex)+(1|raceAA)+(1|recruited_control)+(1|ad)+(1|ftd)+(1|pd)+(1|als)+(1|mci_sci)+(1|mi)+(1|C9Orf72)+(1|GRN)+(1|MAPT)+APOE4.Dose+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+MMSE
##form <- ~ (1|contributor_Fsplit) +RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ
#form <- ~ Age+(1|Sex)+(1|contributor_Fsplit)+RegrBloodPreanalyticFactor.HNRNPA2B1+RegrBloodPreanalyticFactor.HBZ+(1|APOE.E4carrier)  #+APOE.E4carrier.Proxy.LRRN1
form <- ~ Age+(1|Sex)+(1|contributor_Fsplit)+(1|APOE.E4carrier)

library(variancePartition)

# (If regressed) REMOVE regressed proteins to avoid "Response variable 4641 has a variance of 0":
#regrProts.idx<-c(which(grepl("^HNRNPA2B1\\|",rownames(cleanDat))),which(grepl("^HBZ\\|",rownames(cleanDat))))
#regrProts.idx
#5569 6087 -- take out
#previously: 6156 4641
#integer(0)  # if already removed

#previously: removed 2 regr proteins in line:  varPart.reg <- fitExtractVarPartModel(impute::impute.knn(cleanDat[-c(4641,6156),])$data, form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = 8, type = "SOCK"))
varPart.b4c <- fitExtractVarPartModel(na.omit(as.matrix(cleanDat.c1[,])), form, regvars.vp, BPPARAM=BiocParallel::SnowParam(workers = parallelThreads, type = "SOCK"))


vp.b4c <- sortCols(varPart.b4c,FUN=median,last= c("Residuals"))

pdf(file="4p13c1.contributor_Fsplit_19sites1xPAVregr+SiteRegress_Protect_age+sex+carrierBinaryStatus-VariancePartition-PLASMA-7335x22392.pdf", width=15,height=11)
par(mfrow=c(1,1))

plotVarPart( vp.b4c, main="HDS 1.3ms - 4p13c1 - KNOWN APOE 19 sites 2x Regr(2PAV) + Site Regr, Prot age+sex+(e4 binary)" )

# ------------------------------------------------------------------------
# ANNOTATION: Fit and plot final variance-partition estimates for the APOE-
# aware matrix.
# ------------------------------------------------------------------------

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


cleanDat.full<-cleanDat
numericMeta.full<-numericMeta

length(unique(numericMeta.full$person_id))
#17484 (in 22392)  # was 13783 (in 16677)

numericMeta<-numericMeta.full
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

parallelThreads=31
outFilePrefix="4p13c1.QC."
outFileSuffix="AD.CT_volc-ANOVA_tTest.CTimputed"

#Grouping=numericMeta$Group #***
Grouping=numericMeta$Group.withCTimputed
Grouping[which(numericMeta$Group.withCTimputed=="CI.Other")]<-NA

# ------------------------------------------------------------------------
# ANNOTATION: Load the high-precision parANOVA/volcano plotting code used
# for final contrast summaries.
# ------------------------------------------------------------------------
Grouping[which(numericMeta$Group.withCTimputed=="AsymAD")]<-NA
## Different - here we will keep only AD and CT - for a t test volcano
Grouping[which(numericMeta$Group.withCTimputed=="MCI")]<-NA
Grouping[which(numericMeta$Group.withCTimputed=="PD")]<-NA
Grouping[which(numericMeta$Group.withCTimputed=="ALS")]<-NA
Grouping[which(numericMeta$Group.withCTimputed=="FTD")]<-NA


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
flip=c(3)  #c(3:8,15:17)
sameScale=FALSE #TRUE
symbolsOnly=TRUE
highlightGeneProducts=c()  #"HNRNPA2B1","HBZ","SPC25","CPLX2","PTN","MAPT","OMG","GDF15","NEFL","NRGN","CTHRC1","HTRA1","TTR","LRP1","NTN1","SFRP1","GPNMB","PAFAH1B3","CRIP1")
labelHighlighted=TRUE      # if true, highlighted spots get text labels with their rownames from ANOVAout
labelTop=30
plotVolc()                 # runs on ANOVAout as input (need not be specified).

#DEXpercentStacked()        # runs on prior function outputs as input; writes stacked bar plot(s) to PDF.


Grouping<-rep("Female",nrow(numericMeta))
Grouping[numericMeta$Sex==1]<-"Male"
table(Grouping)
#Female   Male
#  9000   8484

outFileSuffix="Sex_volc-ANOVA_tTest"
ANOVAout <- parANOVA.dex()

flip=c(3)  # Female Up on right
plotVolc()                 # runs on ANOVAout as input (need not be specified).


Grouping<-rep("Noncarrier",nrow(numericMeta))
Grouping[numericMeta$APOE4.carrier.imputed==1]<-"e4+ Carrier"
table(Grouping)
#e4+ Carrier  Noncarrier
#       6535       10949

outFileSuffix="e4carrier_volc-ANOVA_tTest"
ANOVAout <- parANOVA.dex()

flip=c(3)  # Female Up on right
plotVolc()                 # runs on ANOVAout as input (need not be specified).


#outputs moved to folder 4p13c1.QC


plotVolcanoBeeswarm <- function(ANOVAout, maxPoints=100, censorRange=0, contrast="", flip=FALSE, xAxisLabel=as.expression(bquote('log'[2] ~ 'Fold Change'))) {
  # ANOVAout: data.frame with p-values in col 3, log2FC in col 4, colors in col 5
  # maxPoints: number of top positive and negative rows to keep
  # censorRange: proportion of effect size range around zero to exclude
  require(beeswarm)

  # Extract relevant columns
  pvals   <- ANOVAout[[3]]
  log2FC  <- ANOVAout[[4]]
  colors  <- ANOVAout[[5]]

  # Compute -log10(p)
  pvals[pvals==0]<-min(1e-320,min(pvals[pvals!=0],na.rm=T))
  logp <- -log10(pvals)

  # Order by significance separately for positive and negative effect sizes
  posIdx <- which(log2FC > 0)
  negIdx <- which(log2FC < 0)

  # tie-aware ordering
  ordPos <- order(pvals[posIdx], -abs(log2FC[posIdx]))
  topPos <- posIdx[ordPos][1:min(maxPoints, length(posIdx))]
  ordNeg <- order(pvals[negIdx], -abs(log2FC[negIdx]))
  topNeg <- negIdx[ordNeg][1:min(maxPoints, length(negIdx))]

  selIdx <- c(topPos, topNeg)

  df <- data.frame(log2FC=log2FC[selIdx],
                   logp=logp[selIdx],
                   color=colors[selIdx])
  rownames(df) <- rownames(ANOVAout)[selIdx]
  if (flip==TRUE) { df$log2FC = (-1)*df$log2FC }

  # Scale point sizes between min and max logp
  this.logpRange=range(df$logp, na.rm=TRUE)
  sizeRange <- if (this.logpRange[1]==this.logpRange[2] & this.logpRange[2]>1.3) { c(1.3,this.logpRange[2]) } else { this.logpRange }
  sizeRange.cex <- c(1.7,6)
  df$size <- scales::rescale(df$logp, to=sizeRange.cex, from=sizeRange)  # adjust min/max sizes

  # Optionally censor a range around zero
  if(censorRange > 0) {
    maxAbs <- max(abs(df$log2FC))
    censor <- censorRange * maxAbs
    df$log2FC[df$log2FC > -censor & df$log2FC < censor] <- NA
  }

  # Plot setup
  op <- par(mar=c(2,1,1,8) + 0.1, xpd=TRUE)  # add extra space on the right

  plot(NA, xlim=range(df$log2FC, na.rm=TRUE)*1.01, ylim=c(0,1),
       xlab="", ylab="", yaxt="n", cex.axis=1.5, type="n") #,
  #     main=paste0("Top Significant Effects n=",nrow(df)))
  #mtext(as.expression(bquote('log'[2] ~ 'Fold Change \n(' * .(contrast) * ')')), side=1, line=1.5, font=2, cex=1.2)
  mtext(xAxisLabel, side=1, line=1.5, font=2, cex=1.2)
  #y axis title:   mtext(contrast, side=2, line=1.5, font=2, cex=1.2)

  # Draw horizontal line at y=0.5
#  abline(h=0.5, col="black", lwd=2)
  usr <- par("usr")   # horizontal line only across xlim
  segments(x0=usr[1], x1=usr[2], y0=0.5, y1=0.5, col="black", lwd=2)
  # Vertical hatch at 0
  segments(x0=0, x1=0, y0=0.65, y1=0.35, col="black", lwd=2.5)


  # Beeswarm jittering
  beeswarm(log2FC~factor(rep("thisGroup",nrow(df))), horizontal=TRUE,
           at=0.5,
           corralWidth=0.9,
           corral="gutter",
           data=df,
           pch=16,
#           col="black",
           pwcol=df$color,
           pwcex=df$size,
           add=TRUE)

  # --- Labeling with simple repel logic ---
  # jitter y positions to reduce overlap
  yjitter <- jitter(rep(0.55, nrow(df)), amount=0.35)
  yjitterFlip <- runif(length(yjitter))
  for (i in 1:length(yjitter)) yjitter[i] = ifelse(yjitterFlip[i] < 0.5, 0.5-(yjitter[i]-0.5), yjitter[i])

  for(i in seq_len(nrow(df))) {
    text(df$log2FC[i], yjitter[i], labels=strsplit(rownames(df)[i],"[|]")[[1]][1], cex=0.95, pos=if(yjitter[i]>0.5) { 3 } else { 1 })
    segments(df$log2FC[i], 0.5, df$log2FC[i], yjitter[i], col="grey30")
  }

  # Legend for point sizes
  legendSizes <- seq(sizeRange[1], sizeRange[2], length.out=5)
  legendCex   <- scales::rescale(legendSizes, to=sizeRange.cex, from=sizeRange)
  legend(x=par("usr")[2] + 0.005 * diff(par("usr")[1:2]),  # a bit to the right of plot
       y=0.95, title.cex=1.4, title.font=2,                # y, near top
       legend=round(legendSizes,1), y.intersp = 2.2, x.intersp = 1.75,
       pt.cex=legendCex, pch=21, col="black", bg="grey",
       title="-log10(p)", bty="n")
  legend(x=usr[1], y=usr[4],legend="", xjust=0,yjust=1,title.font=2,title.cex=2,title=paste0(contrast,"  Top sig. effects: ",length(topNeg)," (-) | ",length(topPos)," (+)"), bty="n")
}


pdf("EffectSize_Plots-17484lastVisits_of22392_volcanoTop20+20-.pdf",width=22,height=12)
par(mfrow=c(3,1))

ANOVAout.e4<-read.csv("4p13c1.QC..ANOVA_diffEx-ALL-e4carrier_volc-ANOVA_tTest.csv",header=TRUE,row.names=1,check.names=FALSE)
plotVolcanoBeeswarm(ANOVAout.e4,20,0,contrast="APOE e4+ vs e4-",flip=TRUE)

ANOVAout.sex<-read.csv("4p13c1.QC..ANOVA_diffEx-ALL-Sex_volc-ANOVA_tTest.csv",header=TRUE,row.names=1,check.names=FALSE)
#ANOVAout.sex[which(ANOVAout.sex[,3]==0),3]<-1e-318  #now handled in function
plotVolcanoBeeswarm(ANOVAout.sex,20,0,contrast="Female vs Male",flip=TRUE)


# ------------------------------------------------------------------------
# ANNOTATION: Plot top effect-size points for APOE e4, sex, and AD-vs-
# control contrasts.
# ------------------------------------------------------------------------
ANOVAout.ADvsCT<-read.csv("4p13c1.QC..ANOVA_diffEx-ALL-AD.CT_volc-ANOVA_tTest.CTimputed.csv",header=TRUE,row.names=1,check.names=FALSE)
plotVolcanoBeeswarm(ANOVAout.ADvsCT,20,0,contrast="AD vs CT",flip=TRUE)

dev.off()

#ESplot<-recordPlot()
#pdf("EffectSize_Plots-17484lastVisits_of22392_volcanoTop20+20-.pdf",width=11,height=8.5)
#  print(ESplot)
#dev.off()


######################################
## Alternative Correlation (to linear trait), stats table with volcanoes

# These parameters are specific to trait correlation statistics generation; traits are provided as columns of the data frame stored in the provided example RData as the variable numericMeta.
cor.traits=c("Age") #
           #c("age_at_visit","MMSE", "Lilly.BH.blood.pTau217","UDS.blood.pTau217", "RegrBloodPreanalyticFactor.HNRNPA2B1","RegrBloodPreanalyticFactor.HBZ",
           #  "TimeToSpin","TimeToDecant","TimeToFreeze","FedFastedTime","FreezeThawCycles")                 # Molecular and quantitative Traits to correlate to in numericMeta columns (colnames)
filter.trait="Group.withCTimputed"             # Trait on which to subset case samples
filter.trait.subsets=c("ALL","AD","CT") # Subsets of case samples will be used for correlation to the cell type proportion estimates
                                    # (4 separate cor.traits x 4 sample subsets = 16 total p and R value columns to generate)
corFn="bicor"                       #'bicor'; other options are 'kendall', 'spearman', ...anything else will cause Pearson (cor) to be used


#source("./parANOVA.dex.R")
CORout <- trait.corStat()                      # runs on cleanDat and Grouping variables as required input.
# Correlation p + R table calculations complete. If you want to use the table with plotVolc(), set the variable corVolc=TRUE and use variable CORout to store the table generated.

dim(CORout)
#[1] 7333   9


outFileSuffix="BICORvolc"
corVolc=TRUE        # changes the behavior of plotVolc, DEXpercentStacked, and GOparallel functions later in the pipeline, to use CORout
useNETcolors=TRUE
sameScale=TRUE
#highlightGeneProducts=rownames(CORout)[which(CORout$Filter.ALLsamples.FavoriteCorrStat.Sig)]  # Specifies which spots should be large
#labelHighlighted=TRUE

##Change values less than 1e-50 to 1e-50, so same scale volcanoes will not be stretched to -log10(p) of 200!
#plateauMinP=1e-50
#CORout[,3:(which(grepl("bicor ",colnames(CORout)))[1]-1)]<-apply(CORout[,3:(which(grepl("bicor ",colnames(CORout)))[1]-1)],2,function(x) { x1<-x; x1[x1<plateauMinP]<-plateauMinP; x1; })
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


## Output bicor plot of volc top to PDF

pdf("EffectSize_Plots-17484lastVisits_of22392_bicor(Age)volcanoTop20+20-.pdf",width=22,height=12)
par(mfrow=c(3,1))

plotVolcanoBeeswarm(CORout[,c(1:3,6,9)],20,0,contrast="Age (bicor)",flip=FALSE,xAxisLabel="Bicor (rho)")
dev.off()

#ESplot<-recordPlot()
#pdf("EffectSize_Plots-17484lastVisits_of22392.pdf",width=26,height=12)
#  print(ESplot)
#dev.off()


# ------------------------------------------------------------------------
# ANNOTATION: Plot the age bicor volcano/effect-size view.
# ------------------------------------------------------------------------

## undo last visit cull
cleanDat<-cleanDat.full
numericMeta<-numericMeta.full

rm(numericMeta.full)
rm(cleanDat.full)


## replot tSNE for 4p13c1 regressed final data

tSNE.4p13c1.plasma.xy.UNREGRESSED<-tSNE.4p13c1.plasma.xy
tSNE.4p13c1.plasma.xy<-tSNE.4p13c1.plasma.samples.sites$data


library(ggplot2)
#library(ggpubr) - rlang upgrade required, ggplot2 upgrade required. cannot install from older source!
library(ggrepel)
library(viridisLite)

## Get the indices of the first occurrence of each unique value in contributor_Fsplit
#first_occurrence_indices.4p13c1<-first_occurrence_indices <- match(unique(regvars.c1$contributor_Fsplit), regvars.c1$contributor_Fsplit)
#
#labels.4p13c1<-labels<-regvars.c1$contributor_Fsplit[first_occurrence_indices.4p13c1]

tSNE.plasma.samples.sites<-ggplot2::ggplot(tSNE.4p13c1.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=regvars.c1$contributor_Fsplit), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE Dimension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13c1.plasma.xy[first_occurrence_indices.4p13c1, ],
                  aes(x=x,y=y, label = labels.4p13c1), size = 6, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +  # Labels for top 10 positive & negative y points
  theme_minimal() +  # Minimal theme
  theme(
    panel.background = element_blank(),  # Remove plot area color
#    legend.position = c(0.095, 0.85),  # Position legend at top right within plot area
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Add background to legend for clarity
    legend.key = element_rect(fill = "white"),  # Keep legend keys clean
    axis.title.x = element_text(size = 28),  # Double x-axis label text size
    axis.title.y = element_text(size = 28)
  )

tSNE.plasma.samples.age<-ggplot2::ggplot(tSNE.4p13c1.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=as.numeric(numericMeta.plasma$age_at_visit)), size=0.35) + scale_color_viridis_c(option = "plasma", name = "Age") +  guides(color = guide_colorbar(barwidth = unit(0.5, "cm"), barheight = unit(4, "cm"))) + labs(x = "tSNE Dimension 1", y = "tSNE Dimension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13c1.plasma.xy[first_occurrence_indices, ],
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

tSNE.plasma.samples.sampleMatrix<-ggplot2::ggplot(tSNE.4p13c1.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_matrix), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE Dimension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13c1.plasma.xy[first_occurrence_indices, ],
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

tSNE.plasma.samples.sampleType<-ggplot2::ggplot(tSNE.4p13c1.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta.plasma$sample_type), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE Dimension 2") +  # Axis labels
  geom_text_repel(data = tSNE.4p13c1.plasma.xy[first_occurrence_indices, ],
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

tSNE.plasma.samples.apoe4dose<-ggplot2::ggplot(tSNE.4p13c1.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=as.character(numericMeta$APOE4.Dose)), size=0.35) + scale_color_manual(name="APOE4 Dose",
    values = c("0" = "lightgreen",
               "1" = "darkorange",
               "2" = "maroon"), labels = c("0 copies", "1 copy", "2 copies")) + labs(x = "tSNE Dimension 1", y = "tSNE Dimension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.4p13c1.plasma.xy[first_occurrence_indices, ],
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

tSNE.plasma.samples.MMSE<-ggplot2::ggplot(tSNE.4p13c1.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta$MMSE), size=0.35) + scale_color_viridis_c(option = "mako", name = "MMSE") +  guides(color = guide_colorbar(barwidth = unit(0.5, "cm"), barheight = unit(4, "cm"))) + labs(x = "tSNE Dimension 1", y = "tSNE Dimension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
  geom_text_repel(data = tSNE.4p13c1.plasma.xy[first_occurrence_indices, ],
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

tSNE.plasma.samples.Group<-ggplot2::ggplot(tSNE.plasma.xy,label=1:ncol(exprMat.plasma)) + geom_point(aes(x=x,y=y, color=numericMeta$Group.withCTimputed), size=0.35) + guides(color = guide_legend(override.aes = list(size = 3.5))) + labs(x = "tSNE Dimension 1", y = "tSNE Dimension 2") +  # Axis labels   #+ guides(color = guide_legend(override.aes = list(size = 3.5))) +
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


pdf(file="4p13c1.2PAVproteinIntrasiteRegressed+Regress19sites_protectAgeSexAPOE.e4_carrierBinary_Fsplit.tSNE-Plasma(7335x22392)-samples_coloredByTraits-RedoneByECBJreq.pdf",width=11,height=9)
  print(tSNE.plasma.samples.sites + labs(colour="Site"))
  print(tSNE.plasma.samples.age)
  print(tSNE.plasma.samples.sampleMatrix)
  print(tSNE.plasma.samples.sampleType)
  print(tSNE.plasma.samples.apoe4dose)
  print(tSNE.plasma.samples.MMSE)
  print(tSNE.plasma.samples.Group + labs(colour="Group"))
dev.off()
