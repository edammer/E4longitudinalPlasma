## =============================================================================
## 6. Endophenotypes_FETSplitLists - FET enrichment of sliding-window /
##    before-after-50 significant-assay lists against pathology &
##    cognition-associated protein ("endophenotype") marker lists
## =============================================================================
##
## PURPOSE
##   Tests whether the proteins identified as significant in each 5-year EYO
##   sliding window (and in simple before-/after-age-50 splits) are enriched
##   for proteins previously reported as associated with neuropathology or
##   cognitive function endophenotypes (Nat Aging cohort marker lists), using
##   two-sided Fisher's Exact Test (FET) overlap enrichment via the lab's
##   general-purpose geneListFET() wrapper (a GOparallel-family FET tool).
##
## STEP-BY-STEP PIPELINE
##   1. Read the 67-window x assay sliding-window significant-hit list
##      (67_5year_sliding_windowSigHits.csv) and remove any non-protein
##      "seq.####" control-assay entries from every column, left-shifting
##      remaining entries upward and padding with blanks
##      (remove_seq_and_shift()), writing the cleaned list back out as
##      67_5year_sliding_windowSigHits_noNonProteinAssays.csv.
##   2. Define the FET run configuration: heatmap color scale/titles, the
##      3 candidate Nat Aging reference marker-list files
##      (refDataFiles[1:3]) and their species codes, and the 2 category list
##      files to be tested (the cleaned 67-window hit list, and a separate
##      before/after-age-50 split list).
##   3. Source geneListFET(), the Seyfried/Emory-lab general FET enrichment
##      function (a GOparallel-family wrapper around hypergeometric
##      enrichment, heatmap plotting, and PDF/CSV export).
##   4. Run geneListFET() once per (scale x legend-scale) combination for the
##      67-window hit list against the ROSMAP pathology-hits reference list
##      (refDataFiles[1]): unadjusted-p and FDR scales, each in both natural
##      and -log legend display.
##   5. Repeat the same 4 scale/legend combinations for the 67-window hit
##      list against the cognitive-function-association FDR<0.05 reference
##      list (refDataFiles[3]). (refDataFiles[2], the uncorrected p<0.05
##      cognitive-function list, is defined for reference/reuse but not run
##      in this script.)
##
## REQUIRED INPUTS
##   - ./geneListFET.R                          (source()'d FET/heatmap function)
##   - 67_5year_sliding_windowSigHits.csv        (per-window significant-assay
##                                                lists; output of the prefix-7
##                                                ExtractSigAssays scripts)
##   - SOMAbkgr_7289.csv                         (full SomaScan assay
##                                                background gene list for FET)
##   - SplitList_BEFORE_AFTER_50yo.csv           (before/after age-50 assay
##                                                split category list)
##   - ROSMAP_Path_Hits(NatAgingS8).csv          (Nat Aging pathology marker
##                                                reference list)
##   - META_CogFnS7_AmyloidosisNoE4.S2-pLT0.05.csv      (Nat Aging cognition
##                                                marker reference list, p<0.05)
##   - META_CogFnS7_AmyloidosisNoE4.S2-FDR_LT0.05.csv   (Nat Aging cognition
##                                                marker reference list, FDR<0.05)
##
## MAJOR OUTPUTS
##   - 67_5year_sliding_windowSigHits_noNonProteinAssays.csv (cleaned hit list)
##   - One heatmap PDF + companion FET-statistics CSV/XLSX per geneListFET()
##     call, named with the "#1a/#1b/#1c/#1d" and "#2a/#2b/#2c/#2d" prefixes
##     used in FileBaseName1 above, summarizing overlap enrichment between
##     the sliding-window hit list and each Nat Aging marker reference list.
## =============================================================================

######  no "seq.(...)" entries as gene symbols - 06-17-2026

############ PREPARATION OF IN-MEMORY VARIABLES (optional) ##########################
rootdir="F:/OneDrive - Emory/Legacy/e4_homozygoteStudy/DL/3.Five_yr_slidingWindow/FET/GitHub_test/"
setwd(rootdir)

# ---- STEP 1. Strip non-protein "seq.####" control-assay entries out of the
# sliding-window significant-hit list and write a cleaned copy. ----
remove_seq_and_shift <- function(infile) {

    # 1. Read CSV (header = TRUE, no row.names)
    df <- read.csv(infile, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

    # 2. Process each column independently
    df.cleaned <- as.data.frame(
        lapply(df, function(col) {

            # Keep only values NOT starting with "seq."
            keep <- col[!grepl("^seq\\.", col)]

            # Number of padding rows needed
            pad.n <- length(col) - length(keep)

            # Return shifted + padded column
            c(keep, rep("", pad.n))
        }),
        stringsAsFactors = FALSE
    )

    # 3. Build output filename
    outfile <- sub("\\.csv$", "_noNonProteinAssays.csv", infile)

    # 4. Write cleaned data frame
    write.csv(df.cleaned, outfile, row.names = FALSE)

    return(outfile)
}

remove_seq_and_shift("67_5year_sliding_windowSigHits.csv")


# ---- STEP 2. FET run configuration: heatmap appearance, reference marker
# list files (Nat Aging pathology/cognition associations), and the two
# category (assay-list) files to be tested for enrichment against them. ----
#################### CONFIGURATION PARAMETERS FULL LIST #############################
##             WITH SAMPLE VALUES GIVEN FOR SAMPLE DATA PROVIDED                   ##

heatmapScale="minusLogFDR"                        					# Accepted options are "p.unadj" or "minusLogFDR"
heatmapTitle1="EYO 5 year Window Hits Overlap with Nat Aging Pathology/Cognition Associated Protein Lists"	# What are your categories (or WGCNA) list of lists based on?
heatmapTitle2="Before/After 50 Marker Overlap with Nat Aging Pathology/Cognition Associated Protein Lists"	# What are your categories (or WGCNA) list of lists based on?
                                                                                        # And What gene lists are your reference lists?
paletteColors="RdPu"                                          # See valid palettes using RColorBrewer::display.brewer.all() or viridisLite:: functions
                                                                # Can be a vector if there are more than 1 refDataFiles (heatmaps to generate)

FileBaseName1="EYO67x5yrWindows_FET_to_NatAgingMarkerLists_7289"
FileBaseName2="Before+After50_FET_to_NatAgingMarkerLists_7289"
refDataDescription="NatAgingMarkerLists"				# One Description of reference Data list(s) specified in PDF file name below
# File Names of Reference List(s): You will get one output PDF page per file
refDataFiles <- c(      "ROSMAP_Path_Hits(NatAgingS8).csv",
			"META_CogFnS7_AmyloidosisNoE4.S2-pLT0.05.csv",
			"META_CogFnS7_AmyloidosisNoE4.S2-FDR_LT0.05.csv")
speciesCode=c("hsapiens","hsapiens","hsapiens") 				# species code(s) for biomaRt (one for each refDataFile)#one for each .csv in refDataFiles


# Use Modules in memory OR a .csv file with your input gene lists.
modulesInMemory=FALSE                              		# Load modules as categories? (If TRUE, categoriesFile not used, but you need cleanDat, net[["colors"]] and numericMeta variables)
categoriesFile1="67_5year_sliding_windowSigHits_noNonProteinAssays.csv"	# File Name of Categories (Lists of Fly genes), only loaded if modulesInMemory=FALSE
								# NOTE this file format has a column for official gene symbols of each module or cluster, with the cluster name/ID as column names in row 1
categorySpeciesCode="hsapiens"				# What species are the gene sybmols in categoriesFile?
categoriesFile2="SplitList_BEFORE_AFTER_50yo.csv"

# Other Options
allowDuplicates=TRUE				# Allow duplicate symbols across different lists for overlap?
						# (should be true if you have general cell type lists and e.g. disease-associated phenotype cell type lists)
resortListsDecreasingSize=FALSE			# resort categories/modules and reference data lists? (decreasing size order)
barOption=FALSE					# draw bar charts for each list overlap instead of a heatmap.
adjustFETforLookupEfficiency=FALSE		# adjust p FET input for cross-species lookup inefficiency/loss of list member counts?
verticalCompression=3				# Plot(s) are squeezed into 1 row out of this many in each PDF page, compressing the heatmap tracks vertically (or the bar chart heights) for each reference list)
reproduceHistoricCalc=FALSE			# should be FALSE unless trying to reproduce exact calculations of prior publications listed.
#####################################################################################


## Generate Sample Outputs

# ---- STEP 3. Load the Seyfried/Emory FET enrichment wrapper geneListFET(). ----
# Load Seyfried/Emory pipeline FET as function geneListFET() having all the parameters described above, many with defaults used.
source("./geneListFET.R")

# ---- STEP 4. Run FET enrichment of the 67-window sliding-window hit list
# against the ROSMAP pathology marker list (refDataFiles[1]); 4 calls cover
# unadjusted-p vs FDR color scale, each with natural vs -log legend scale. ----
geneListFET(FileBaseName=paste0("#1a.p.unadj_-log-UNT+fullBG-",FileBaseName1), verticalCompression=1, asterisksOnly=TRUE, heatmapScale="p.unadj",
            heatmapTitle=paste0("",heatmapTitle1), paletteColors=paletteColors, maxPcolor=0.05, bkgrFileForCategories="SOMAbkgr_7289.csv", strictSymmetry=TRUE,legendScale="minusLog",
            modulesInMemory=FALSE, categoriesFile=categoriesFile1, categorySpeciesCode=categorySpeciesCode,  #use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
            refDataFiles=refDataFiles[1],speciesCode=speciesCode,refDataDescription=refDataDescription)  #file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?

geneListFET(FileBaseName=paste0("#1b.p.unadj_unlog-UNT+fullBG-",FileBaseName1), verticalCompression=1, asterisksOnly=TRUE, heatmapScale="p.unadj",
            heatmapTitle=paste0("",heatmapTitle1), paletteColors=paletteColors, maxPcolor=0.05, bkgrFileForCategories="SOMAbkgr_7289.csv", strictSymmetry=TRUE,legendScale="unlog",
            modulesInMemory=FALSE, categoriesFile=categoriesFile1, categorySpeciesCode=categorySpeciesCode,  #use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
            refDataFiles=refDataFiles[1],speciesCode=speciesCode,refDataDescription=refDataDescription)  #file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?

geneListFET(FileBaseName=paste0("#1c.FDR_-log-UNT+fullBG-",FileBaseName1), verticalCompression=1, asterisksOnly=TRUE, heatmapScale="minusLogFDR",legendScale="minusLog",
            heatmapTitle=paste0("",heatmapTitle1), paletteColors=paletteColors, maxPcolor=0.05, bkgrFileForCategories="SOMAbkgr_7289.csv", strictSymmetry=TRUE,
            modulesInMemory=FALSE, categoriesFile=categoriesFile1, categorySpeciesCode=categorySpeciesCode,  #use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
            refDataFiles=refDataFiles[1],speciesCode=speciesCode,refDataDescription=refDataDescription)  #file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?

geneListFET(FileBaseName=paste0("#1d.FDR_unlog-UNT+fullBG-",FileBaseName1), verticalCompression=1, asterisksOnly=TRUE, heatmapScale="minusLogFDR",legendScale="unlog",
            heatmapTitle=paste0("",heatmapTitle1), paletteColors=paletteColors, maxPcolor=0.05, bkgrFileForCategories="SOMAbkgr_7289.csv", strictSymmetry=TRUE,
            modulesInMemory=FALSE, categoriesFile=categoriesFile1, categorySpeciesCode=categorySpeciesCode,  #use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
            refDataFiles=refDataFiles[1],speciesCode=speciesCode,refDataDescription=refDataDescription)  #file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?


# ---- STEP 5. Repeat the same 4 scale/legend combinations, this time against
# the cognition-association FDR<0.05 marker list (refDataFiles[3]). ----
geneListFET(FileBaseName=paste0("#2a.p.unadj_-log-UNT+fullBG-",FileBaseName1), verticalCompression=1, asterisksOnly=TRUE, heatmapScale="p.unadj",
            heatmapTitle=paste0("",heatmapTitle1), paletteColors=paletteColors, maxPcolor=0.05, bkgrFileForCategories="SOMAbkgr_7289.csv", strictSymmetry=TRUE,legendScale="minusLog",
            modulesInMemory=FALSE, categoriesFile=categoriesFile1, categorySpeciesCode=categorySpeciesCode,  #use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
            refDataFiles=refDataFiles[3],speciesCode=speciesCode,refDataDescription=refDataDescription)  #file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?

geneListFET(FileBaseName=paste0("#2b.p.unadj_unlog-UNT+fullBG-",FileBaseName1), verticalCompression=1, asterisksOnly=TRUE, heatmapScale="p.unadj",
            heatmapTitle=paste0("",heatmapTitle1), paletteColors=paletteColors, maxPcolor=0.05, bkgrFileForCategories="SOMAbkgr_7289.csv", strictSymmetry=TRUE,legendScale="unlog",
            modulesInMemory=FALSE, categoriesFile=categoriesFile1, categorySpeciesCode=categorySpeciesCode,  #use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
            refDataFiles=refDataFiles[3],speciesCode=speciesCode,refDataDescription=refDataDescription)  #file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?

geneListFET(FileBaseName=paste0("#2c.FDR_-log-UNT+fullBG-",FileBaseName1), verticalCompression=1, asterisksOnly=TRUE, heatmapScale="minusLogFDR",legendScale="minusLog",
            heatmapTitle=paste0("",heatmapTitle1), paletteColors=paletteColors, maxPcolor=0.05, bkgrFileForCategories="SOMAbkgr_7289.csv", strictSymmetry=TRUE,
            modulesInMemory=FALSE, categoriesFile=categoriesFile1, categorySpeciesCode=categorySpeciesCode,  #use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
            refDataFiles=refDataFiles[3],speciesCode=speciesCode,refDataDescription=refDataDescription)  #file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?

geneListFET(FileBaseName=paste0("#2d.FDR_unlog-UNT+fullBG-",FileBaseName1), verticalCompression=1, asterisksOnly=TRUE, heatmapScale="minusLogFDR",legendScale="unlog",
            heatmapTitle=paste0("",heatmapTitle1), paletteColors=paletteColors, maxPcolor=0.05, bkgrFileForCategories="SOMAbkgr_7289.csv", strictSymmetry=TRUE,
            modulesInMemory=FALSE, categoriesFile=categoriesFile1, categorySpeciesCode=categorySpeciesCode,  #use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
            refDataFiles=refDataFiles[3],speciesCode=speciesCode,refDataDescription=refDataDescription)  #file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?
