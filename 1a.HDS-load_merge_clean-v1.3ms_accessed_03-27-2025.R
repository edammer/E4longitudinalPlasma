##############################################################################
# Pipeline annotation header: 1a.HDS-load_merge_clean-v1.3ms_accessed_03-27-2025.R
# Manuscript code section(s): 0 / 1
#
# Purpose:
# Download GNPC HDS v1.3ms clinical, SomaScan, genetics, lookup, person-map,
# and proteomics-metadata tables; construct a raw log2-ready proteomic
# matrix and harmonized sample metadata table.
#
# Principal inputs:
#   - PostgreSQL DSN "PostgreSQL" on GNPC ADDI Workbench
#   - ClinicalV1_3ms, SomalogicAnalyteInfoV1_3ms, Somalogic01-12V1_3ms,
#     GeneticsV1_3ms, GNPCLookupV1_3ms, PersonMappingV1_3ms,
#     SomalogicMetaV1_3ms
#
# Principal outputs:
#   - ClinicalV1_3MS_032725.csv
#   - AnalytesV1_3MS_032725.csv
#   - loadedV1_3ms_03-27-25.RData
#
# Step overview:
#   1. Connect to GNPC PostgreSQL and enumerate/source data tables.
#   2. Create unique SomaScan assay row identifiers using gene symbol, UniProt
#      ID, SomaID, and aptamer name as needed for replicated targets.
#   3. Merge the twelve Somalogic matrices by sample_id, remove duplicated
#      metadata columns, transpose to assays x samples, and replace sentinel
#      or invalid RFU values with NA.
#   4. Map APOE and other genetic annotations into sample metadata, then
#      append proteomics metadata, clinical traits, and person identifiers.
#   5. Split the merged data into plasma/serum and CSF subsets for downstream
#      analyses.
#
# Notes added during manuscript-code alignment:
#   - Original executable code was preserved except for whitespace/encoding cleanup
#     and the explicitly marked non-run archive described below when applicable.
#   - Comments inserted by this pass are prefixed with "ANNOTATION:".
##############################################################################
rootdir="/home/workspace/files/EBD/"
setwd(rootdir)


# ------------------------------------------------------------------------
# ANNOTATION: Connect to the GNPC PostgreSQL database on the ADDI Workbench
# and load packages used for table extraction and high-volume matrix
# assembly.
# ------------------------------------------------------------------------
#New Connection (menu), postgreSQL (establish connection to DB)
library(DBI)
con <- dbConnect(odbc::odbc(), dsn="PostgreSQL", timeout = 10)

library(xaputils)
library(data.table)


# ------------------------------------------------------------------------
# ANNOTATION: Inventory and download the HDS v1.3ms source tables used
# throughout the downstream pipeline.
# ------------------------------------------------------------------------
#to receive a list of all available tables
tbls <- xap.list_tables()
#to read the data from a single table
clinical <- xap.read_table("ClinicalV1_3ms")
#previously 2/12/2025: date error, fixed 2/13/2025
dim(clinical)
#[1] 31083    56  in v1.3ms  previously:  31111    52
write.csv(clinical,"ClinicalV1_3MS_032725.csv")



# ------------------------------------------------------------------------
# ANNOTATION: Prepare SomaScan analyte metadata and construct unique assay
# identifiers suitable for row names in the expression matrix.
# ------------------------------------------------------------------------
analytes<-xap.read_table("SomalogicAnalyteInfoV1_3ms")
analytes$origNames<-gsub("\\.","_",analytes$apt_name)

## Get simple Unique ID feature names to set columns of dat.df
# fix missing UniProt and Symbols for non-proteins to SomaId and AptamerName
analytes$uni_prot[!analytes$type=="Protein"]<-analytes$soma_id[!analytes$type=="Protein"]
analytes$entrez_gene_symbol[!analytes$type=="Protein"]<-analytes$apt_name[!analytes$type=="Protein"]
# fix missing Uniprot and Symbols for odd proteins
#analytes[which(analytes$EntrezGeneSymbol==""),] #all are for Fc_MOUSE
analytes$entrez_gene_symbol[analytes$uni_prot=="Q99LC4"]<-analytes$apt_name[analytes$uni_prot=="Q99LC4"]
#analytes$entrez_gene_symbol[na.omit(which(analytes$uni_prot=="Q99LC4"))]<-analytes$apt_name[na.omit(which(analytes$uni_prot=="Q99LC4"))]

uniqueIDs<-paste0(as.data.frame(do.call("rbind",strsplit(analytes$entrez_gene_symbol,"[|]")))[,1], "|", as.data.frame(do.call("rbind",strsplit(analytes$uni_prot,"[|]")))[,1] )
#analytes[5414:5434,]  # rows of featData without Uniprot or Enrtrez Symbol are Spuriomer, Non-Biotin, or Non-Cleavable type control Somamers
#uniqueIDs[5414:5434]
#analytes[c(9,12,18,77),]
#uniqueIDs[c(9,12,18,77)]

length(uniqueIDs)
#[1] 7644 (in v1.3ms)  previously: 7596
length(unique(uniqueIDs))
#[1] 6717 (in v1.3ms)  previously: 6669

techRepUniqueIDs<-uniqueIDs[duplicated(uniqueIDs)]
#techRepUniqueIDs

length(techRepUniqueIDs)
#[1] 927 (same in v1.3ms)

nonUniqueIDsAtSomamer<-paste0(uniqueIDs[which(uniqueIDs %in% techRepUniqueIDs)], "^", analytes$soma_id[which(uniqueIDs %in% techRepUniqueIDs)], "@", analytes$apt_name[which(uniqueIDs %in% techRepUniqueIDs)])   #^SomaId @AptName
length(nonUniqueIDsAtSomamer)
#1726 now unique IDs for all replicated measures (different somamers)
uniqueIDs[which(uniqueIDs %in% techRepUniqueIDs)]<-nonUniqueIDsAtSomamer

which(duplicated(uniqueIDs))  #requires AptName to be unique...
#integer(0)

analytes$uniqueID<-uniqueIDs

write.csv(analytes,"AnalytesV1_3MS_032725.csv")


# ------------------------------------------------------------------------
# ANNOTATION: Read all twelve Somalogic abundance table partitions before
# reassembling them into one sample-by-assay table.
# ------------------------------------------------------------------------
p01<-xap.read_table("Somalogic01V1_3ms")
p02<-xap.read_table("Somalogic02V1_3ms")
p03<-xap.read_table("Somalogic03V1_3ms")
p04<-xap.read_table("Somalogic04V1_3ms")
p05<-xap.read_table("Somalogic05V1_3ms")
p06<-xap.read_table("Somalogic06V1_3ms")
p07<-xap.read_table("Somalogic07V1_3ms")
p08<-xap.read_table("Somalogic08V1_3ms")
p09<-xap.read_table("Somalogic09V1_3ms")
p10<-xap.read_table("Somalogic10V1_3ms")
p11<-xap.read_table("Somalogic11V1_3ms")
p12<-xap.read_table("Somalogic12V1_3ms")

# Grab other tables from SQL DB:
genetics <- xap.read_table("GeneticsV1_3ms")
dim(genetics)
#[1] 21332    4  previously: [1] 22287     4
GNPClookup <- xap.read_table("GNPCLookupV1_3ms")
dim(GNPClookup)
#[1] 205  4  previously: [1] 177   4
personMapping <- xap.read_table("PersonMappingV1_3ms")
dim(personMapping)
#[1] 31083    4  previously: 31111     2
proteomicsMeta <- xap.read_table("SomalogicMetaV1_3ms")  #previously ProteomicsMetaV1_3
dim(proteomicsMeta)
#[1] 33848   21  previously: 31111    20

#interim
save.image("loadedV1_3ms_03-27-25.RData")  #previously: loadedV1_3.RData


dat<-cbind(as.data.frame(p01),as.data.frame(p02[match(p01$sample_id,p02$sample_id),]))
dat<-cbind(dat,as.data.frame(p03[match(p01$sample_id,p03$sample_id),]))

# ------------------------------------------------------------------------
# ANNOTATION: Column-bind Somalogic partitions by sample_id so all assay
# values are aligned per sample.
# ------------------------------------------------------------------------
dat<-cbind(dat,as.data.frame(p04[match(p01$sample_id,p04$sample_id),]))
dat<-cbind(dat,as.data.frame(p05[match(p01$sample_id,p05$sample_id),]))
dat<-cbind(dat,as.data.frame(p06[match(p01$sample_id,p06$sample_id),]))
dat<-cbind(dat,as.data.frame(p07[match(p01$sample_id,p07$sample_id),]))
dat<-cbind(dat,as.data.frame(p08[match(p01$sample_id,p08$sample_id),]))
dat<-cbind(dat,as.data.frame(p09[match(p01$sample_id,p09$sample_id),]))
dat<-cbind(dat,as.data.frame(p10[match(p01$sample_id,p10$sample_id),]))
dat<-cbind(dat,as.data.frame(p11[match(p01$sample_id,p11$sample_id),]))
dat<-cbind(dat,as.data.frame(p12[match(p01$sample_id,p12$sample_id),]))

head(dat[,duplicated(colnames(dat))])
#-- same values in all duplicate columns

# assign final uniqueIDs to analytes.
dat2<-dat[,which(!duplicated(colnames(dat)))]

# ------------------------------------------------------------------------
# ANNOTATION: Remove duplicated metadata columns, rename assay columns to
# unique identifiers, and transpose to assays x samples.
# ------------------------------------------------------------------------
oldIDs.dat2<-colnames(dat2)[5:ncol(dat2)]
unknownIDs.idx<-which(is.na(match(oldIDs.dat2,analytes$origNames))) +4
dat2.unknownAnalytes<-dat2[,c(1:4, unknownIDs.idx)]

which(duplicated(dat2$sample_id))
#integer(0)   previously:  [1] 18911 18918 18919 20287 20289 20430 21933 21953
sampleIDs.rownames<-make.unique(dat2$sample_id)

dat3<-dat2[,c(5:ncol(dat2))[which(!5:ncol(dat2) %in% unknownIDs.idx)] ]
oldIDs.dat3<-colnames(dat3)
colnames(dat3)<-analytes$uniqueID[match(oldIDs.dat3,analytes$origNames)]
rownames(dat3)<-sampleIDs.rownames

dat3.unknownIDs<-dat2[,c(5:ncol(dat2))[which(5:ncol(dat2) %in% unknownIDs.idx)] ]
rownames(dat3.unknownIDs)<-sampleIDs.rownames

numericMeta0<-dat2[,2:4]
rownames(numericMeta0)<-sampleIDs.rownames

#dat4<-rbind(as.vector(dat3[,1]),as.vector(dat3[,2]))
#dat4<-for(i in 3:ncol(dat3)) dat4<-rbind(dat4,dat3[,i])
exprMat0<-t(cbind(as.matrix(dat3),as.matrix(dat3.unknownIDs)))
colnames(exprMat0)<-sampleIDs.rownames
dim(exprMat0)
#[1]  7745 33848   previously: 7744 31119
dim(dat2)
#[1] 33848  7749   previously: 31119  7748


## CLEAR MEMORY SPACE  (cut 8GB from 14.x GB footprint)
rm(list=c("p01","p02","p03","p04","p05","p06","p07","p08","p09","p10","p11","p12","con"))

# ------------------------------------------------------------------------
# ANNOTATION: Free large source objects before sample/assay value cleaning
# to reduce memory use in the Workbench session.
# ------------------------------------------------------------------------


# all samples with -1 values
sapply(1:7745, function(x) if(length(which(exprMat0[,x]==-1))>0) x )

# ------------------------------------------------------------------------
# ANNOTATION: Audit Somalogic sentinel values; -1 and invalid non-positive
# RFU values are converted to NA before log2 handling.
# ------------------------------------------------------------------------
#7745 indices! (all)  previously 7744

#total samples with -1 values
length(sapply(1:7745, function(x) if(length(which(exprMat0[,x]==-1))>0) x))
#7745

# total -1 values
minus1.count<-sum(sapply(1:7745, function(x) length(which(exprMat0[,x]==-1)) ))
minus1.count
#3213531  previously: 1153856
total.count<-nrow(exprMat0)*ncol(exprMat0)
total.count
#262152760  previously: 240985536
minus1.count/total.count
#0.01225824  previously: 0.004788072
sum(sapply(1:7745, function(x) length(which(is.na(exprMat0[,x]))) ))
#0

exprMat0[exprMat0==-1]<-NA

length(which(exprMat0==0))
#[1] 0
length(which(exprMat0<=0))
#[1] 2
exprMat0[which(exprMat0<=0)]
#[1] -0.06452822 -0.05013292  current v1.3ms
#[1] -0.05013292 -0.06452822  previously

exprMat0[which(exprMat0<0)]<-NA

# population of standard 7k aptamers across all samples
hist(log2(exprMat0[1:7644,]),breaks=150)  #previously 7595
#population centers at 10, left tail fattens, falls off before reaching 1

# population of measurements of 5k custom aptamers without names
hist(log2(exprMat0[7645:nrow(exprMat0),]),breaks=150,col="#FF999999", add=FALSE)
#left peak much higher than pop center at 0 for unknown aptamers.


## Grab other tables from SQL DB:  (now grabbed above followed by interim save on RStudio in Workbench)
#genetics <- xap.read_table("GeneticsV1_3ms")
#dim(genetics)
##[1] 21332    4  previously: [1] 22287     4
#GNPClookup <- xap.read_table("GNPCLookupV1_3ms")
#dim(GNPClookup)
##[1] 205  4  previously: [1] 177   4
#personMapping <- xap.read_table("PersonMappingV1_3ms")
#dim(personMapping)
##[1] 31083    4  previously: 31111     2
#proteomicsMeta <- xap.read_table("SomalogicMetaV1_3ms")  #previously ProteomicsMetaV1_3
#dim(proteomicsMeta)
##[1] 33848   21  previously: 31111    20

#interim
save.image("loadedV1_3ms_03-27-25.RData")  #previously: loadedV1_3.RData


# Examine key variables in proteomicsMeta, genetics tables
table(proteomicsMeta$sample_matrix)

# ------------------------------------------------------------------------
# ANNOTATION: Summarize proteomics metadata, genetics fields, and lookup
# codes for reproducibility/QC.
# ------------------------------------------------------------------------
## v1.3ms now:
#        -1 Citrate Plasma            CSF    EDTA Plasma        Serum
#      2765            638           3233          23000         4212
## previously:
# Citrate Plasma            CSF    EDTA Plasma
#           648           3233          27230

table(proteomicsMeta$sample_type)
## v1.3ms now:
#    Buffer Calibrator         QC     Sample
#       781       1240        744      31083
## previously (v1.3)
#    -1 Sample
#  4696  26415

table(proteomicsMeta$row_check)
## v1.3ms now:
#  -1  FLAG  PASS
#3410   815 29623
## previously (v1.3)
#   -1  FLAG  PASS
# 3366   810 26935
head(genetics)
## v1.3ms now:
#  contributor_code                            sample_id gene variant
#1                G 2a1ac6bc-2a6f-4f7c-91db-2b1bead63fbd APOE      34
#2                W 5d12a257-14b7-473e-a829-2b46ddb41083 APOE      44
#3                F e55867f7-d1f0-47eb-b595-2b74c696fe4b APOE      24
#4                L d42c2793-0ead-4c0e-8fc5-2b858e881396 APOE      33
#5                C 8060f61d-386d-4d6f-9fc9-2c32d2284ff3 APOE      33
#6                N aadd8c24-e5ad-4ec3-ad2c-2c5746ea2dc0 MAPT       1
## previously (v1.3)
#  contributor_code                            sample_id gene variant
#1                F 8ec1ccc9-c387-4fca-bd22-000048ff66b3 APOE      24
#2                G 70dde1a5-f19c-4102-9947-00063ee7fcac APOE      34
#3                Q fb56cab2-fd4e-4f07-a1e1-000b7224f985 APOE      33
#4                Q 058f0c61-56d5-4af7-a838-cc36d70399ee APOE      33
#5                D 6a9ad566-ea63-4791-ae03-863b02ac70a8 APOE      -1
#6                F cefb9d4a-fec3-4d3f-a3d2-00335e939178 APOE      34
table(genetics$gene)
#    ANG    APOE      C9 C9ORF72     FUS     GRN    MAPT PROGRAN    SETX    SOD1     TAU   TDP43    VAPB     VCP
#      8   18929     576     341      10     689     726       3       8      19       2       5       7       9
## identical to prior
#    ANG    APOE      C9 C9ORF72     FUS     GRN    MAPT PROGRAN    SETX    SOD1     TAU   TDP43    VAPB     VCP
#      8   19884     576     341      10     689     726       3       8      19       2       5       7       9
table(genetics$variant)
## v1.3ms now:
#  -1    0    1   22   23   24   32   33   34   44
#1445 1864  539   65 1550  404  147 8938 5430  950
## previously (v1.3)
#  -1    0    1   22   23   24   32   33   34   44
#2381 1864  539   66 1551  403  147 8938 5441  957

# Integrate genetic metadata into numericMeta data frame
numericMeta0$APOE <- NA

# ------------------------------------------------------------------------
# ANNOTATION: Map APOE and other genetic variants from the genetics table
# onto the sample-level metadata.
# ------------------------------------------------------------------------
length(match(rownames(numericMeta0),genetics$sample_id))
#[1] 33848  previously: 31119
genetics.apoe.idx<-intersect(match(rownames(numericMeta0),genetics$sample_id),which(genetics$gene=="APOE"))
length(genetics.apoe.idx)
#[1] 18929  previously: 19884
genetics.apoe<-genetics$variant[genetics.apoe.idx]
length(unique(genetics$sample_id))
#[1] 20110  previously: 21065
names(genetics.apoe)<-genetics$sample_id[genetics.apoe.idx]
length(unique(names(genetics.apoe)))
#[1] 18929  previously: 19884
numericMeta0$APOE[which(rownames(numericMeta0) %in% names(genetics.apoe))]<-genetics.apoe[na.omit(match(rownames(numericMeta0),names(genetics.apoe)))]

genetics.c9.idx<-intersect(match(rownames(numericMeta0),genetics$sample_id),which(genetics$gene=="C9" | genetics$gene=="C9ORF72"))
genetics.c9<-genetics$variant[genetics.c9.idx]
names(genetics.c9)<-genetics$sample_id[genetics.c9.idx]
length(unique(names(genetics.c9)))
#565 (previously 609)  -- less than 576+341 -- same cases under both?
numericMeta0$C9Orf72<-NA
numericMeta0$C9Orf72[which(rownames(numericMeta0) %in% names(genetics.c9))]<-genetics.c9[na.omit(match(rownames(numericMeta0),names(genetics.c9)))]
table(genetics.c9)
## v1.3ms now:
#  0   1
#300 265
## previously
#  0   1
#345 264

genetics.grn.idx<-intersect(match(rownames(numericMeta0),genetics$sample_id),which(genetics$gene=="GRN" | genetics$gene=="PROGRAN"))
genetics.grn<-genetics$variant[genetics.grn.idx]
names(genetics.grn)<-genetics$sample_id[genetics.grn.idx]
length(unique(names(genetics.grn)))
#295  previously: 241
numericMeta0$GRN<-NA
numericMeta0$GRN[which(rownames(numericMeta0) %in% names(genetics.grn))]<-genetics.grn[na.omit(match(rownames(numericMeta0),names(genetics.grn)))]
table(genetics.grn)
## v1.3ms now:
#  0   1
#184 111
## previously
#  0   1
#130 111

genetics.mapt.idx<-intersect(match(rownames(numericMeta0),genetics$sample_id),which(genetics$gene=="MAPT" | genetics$gene=="TAU"))
genetics.mapt<-genetics$variant[genetics.mapt.idx]
names(genetics.mapt)<-genetics$sample_id[genetics.mapt.idx]
length(unique(names(genetics.mapt)))
#304  previously: 316
numericMeta0$MAPT<-NA
numericMeta0$MAPT[which(rownames(numericMeta0) %in% names(genetics.mapt))]<-genetics.mapt[na.omit(match(rownames(numericMeta0),names(genetics.mapt)))]
table(genetics.mapt)
## v1.3ms now:
#  0   1
#154 150
## previously
#  0   1
#166 150

# Integrate proteomic metadata select traits into numericMeta data frame
length(na.omit(match(rownames(numericMeta0),proteomicsMeta$sample_id)))

# ------------------------------------------------------------------------
# ANNOTATION: Append proteomics metadata, clinical variables, and person IDs
# to build the first numericMeta object.
# ------------------------------------------------------------------------
#[1] 33848  -- all rows match
length(unique(na.omit(match(rownames(numericMeta0),proteomicsMeta$sample_id))))
#33848  previously: 31111
numericMeta.proteomicsMeta<-proteomicsMeta[match(rownames(numericMeta0),proteomicsMeta$sample_id),c("contributor_code","visit","units","row_check","sample_matrix","sample_type")]
#units has log2(RFU) for 1156 samples, others are unlogged.
table(proteomicsMeta$units)
#      -1 log2 RFU
#   32692     1156  previously: 29955     1156

# Integrate clinical metadata select traits into numericMeta data frame
numericMeta.clinical<-clinical[match(rownames(numericMeta0),clinical$sample_id),]
numericMeta.clinical[numericMeta.clinical<0]<-NA

# Integrate person mapping into numericMeta data frame
numericMeta.personMap<-personMapping[match(rownames(numericMeta0),personMapping$sample_id),]

## Finalize numericMeta
numericMeta<-as.data.frame(cbind(numericMeta0,numericMeta.proteomicsMeta,numericMeta.clinical,numericMeta.personMap))
dim(numericMeta)
#[1] 33848    73  previously: 31119    67

#split plasma from CSF sample_matrix
numericMeta.plasma<-numericMeta[which(numericMeta$sample_matrix=="Citrate Plasma" | numericMeta$sample_matrix=="EDTA Plasma" | numericMeta$sample_matrix=="Serum"),]

# ------------------------------------------------------------------------
# ANNOTATION: Create plasma/serum and CSF-specific subsets for separate
# downstream analysis branches.
# ------------------------------------------------------------------------
exprMat0.plasma<-exprMat0[,match(rownames(numericMeta.plasma),colnames(exprMat0))]

numericMeta.CSF<-numericMeta[which(numericMeta$sample_matrix=="CSF"),]
exprMat0.CSF<-exprMat0[,match(rownames(numericMeta.CSF),colnames(exprMat0))]

rm(list=c("dat","dat2","dat3","xap.conn"))

save.image("loadedV1_3ms_03-27-25.RData")  #previously: loadedV1_3.RData


# ------------------------------------------------------------------------
# ANNOTATION: Save the full raw extraction and merged metadata workspace for
# subsequent trait cleaning.
# ------------------------------------------------------------------------

GNPClookup
# CURRENT FULL CONTENTS (v1.3ms):
       table_name                             column_name     key_lookup                               key_description
1        Clinical                      computed_age_range             -1                                           N/A
2        Clinical                      computed_age_range              0                                   Z=-1 to Z=1
3        Clinical                      computed_age_range              1                   Z=-2 to Z=-1 and Z=1 to Z=2
4        Clinical                      computed_age_range              2                   Z=-3 to Z=-2 and Z=2 to Z=3
5        Clinical                      computed_age_range              3                                   Z<3 and Z>3
6        Clinical                                     sex             -1                                    Not Stated
7        Clinical                                     sex              1                                          Male
8        Clinical                                     sex              2                                        Female
9        Clinical                                    race             -1                                    Not Stated
10       Clinical                                    race              1                 American_Indian/Alaska_Native
11       Clinical                                    race              2                     Black_or_African_American
12       Clinical                                    race              3     Native_Hawaiian_or_Other_Pacific_Islander
13       Clinical                                    race              4                                         Asian
14       Clinical                                    race              5                                         White
15       Clinical                                    race             66                               Other(Group_C)
16       Clinical                                    race             77               Non_Caucasian_or_Black(Group_B)
17       Clinical                                    race             88                        Non_Caucasian(Group_A)
18       Clinical          computed_years_education_range             -1                                           N/A
19       Clinical          computed_years_education_range              0                                   Z=-1 to Z=1
20       Clinical          computed_years_education_range              1                   Z=-2 to Z=-1 and Z=1 to Z=2
21       Clinical          computed_years_education_range              2                   Z=-3 to Z=-2 and Z=2 to Z=3
22       Clinical          computed_years_education_range              3                                   Z<3 and Z>3
23       Clinical                   computed_height_range             -2                             Incompatible data
24       Clinical                   computed_height_range             -1                                           N/A
25       Clinical                   computed_height_range              0                                   Z=-1 to Z=1
26       Clinical                   computed_height_range              1                   Z=-2 to Z=-1 and Z=1 to Z=2
27       Clinical                   computed_height_range              2                   Z=-3 to Z=-2 and Z=2 to Z=3
28       Clinical                   computed_height_range              3                                   Z<3 and Z>3
29       Clinical                   computed_weight_range             -1                                           N/A
30       Clinical                   computed_weight_range              0                                   Z=-1 to Z=1
31       Clinical                   computed_weight_range              1                   Z=-2 to Z=-1 and Z=1 to Z=2
32       Clinical                   computed_weight_range              2                   Z=-3 to Z=-2 and Z=2 to Z=3
33       Clinical                   computed_weight_range              3                                   Z<3 and Z>3
34       Clinical                      computed_bmi_range             -1                                           N/A
35       Clinical                      computed_bmi_range              0                                   Z=-1 to Z=1
36       Clinical                      computed_bmi_range              1                   Z=-2 to Z=-1 and Z=1 to Z=2
37       Clinical                      computed_bmi_range              2                   Z=-3 to Z=-2 and Z=2 to Z=3
38       Clinical                      computed_bmi_range              3                                   Z<3 and Z>3
39       Clinical                    computed_pulse_range             -1                                           N/A
40       Clinical                    computed_pulse_range              0                                   Z=-1 to Z=1
41       Clinical                    computed_pulse_range              1                   Z=-2 to Z=-1 and Z=1 to Z=2
42       Clinical                    computed_pulse_range              2                   Z=-3 to Z=-2 and Z=2 to Z=3
43       Clinical                    computed_pulse_range              3                                   Z<3 and Z>3
44       Clinical  computed_systolic_blood_pressure_range             -1                                           N/A
45       Clinical  computed_systolic_blood_pressure_range              0                                   Z=-1 to Z=1
46       Clinical  computed_systolic_blood_pressure_range              1                   Z=-2 to Z=-1 and Z=1 to Z=2
47       Clinical  computed_systolic_blood_pressure_range              2                   Z=-3 to Z=-2 and Z=2 to Z=3
48       Clinical  computed_systolic_blood_pressure_range              3                                   Z<3 and Z>3
49       Clinical computed_diastolic_blood_pressure_range             -2                             Incompatible data
50       Clinical computed_diastolic_blood_pressure_range             -1                                           N/A
51       Clinical computed_diastolic_blood_pressure_range              0                                   Z=-1 to Z=1
52       Clinical computed_diastolic_blood_pressure_range              1                   Z=-2 to Z=-1 and Z=1 to Z=2
53       Clinical computed_diastolic_blood_pressure_range              2                   Z=-3 to Z=-2 and Z=2 to Z=3
54       Clinical computed_diastolic_blood_pressure_range              3                                   Z<3 and Z>3
55       Clinical                              alcohol_hx             -1                                    Not Stated
56       Clinical                              alcohol_hx              0                                         Never
57       Clinical                              alcohol_hx              1                                       Current
58       Clinical                              alcohol_hx              2                                      Previous
59       Clinical                              alcohol_hx              3                             Never_or_Previous
60       Clinical                              smoking_hx             -1                                    Not Stated
61       Clinical                              smoking_hx              0                                         Never
62       Clinical                              smoking_hx              1                                       Current
63       Clinical                              smoking_hx              2                                      Previous
64       Clinical                              smoking_hx              3                             Never_or_Previous
65       Clinical                   computed_years_smoked             -2                             Incompatible data
66       Clinical                   computed_years_smoked             -1                                           N/A
67       Clinical                   computed_years_smoked              0                                   Z=-1 to Z=1
68       Clinical                   computed_years_smoked              1                   Z=-2 to Z=-1 and Z=1 to Z=2
69       Clinical                   computed_years_smoked              2                   Z=-3 to Z=-2 and Z=2 to Z=3
70       Clinical                   computed_years_smoked              3                                   Z<3 and Z>3
71       Clinical                                  stroke             -1                                    Not Stated
72       Clinical                                  stroke              0                                            No
73       Clinical                                  stroke              1                                           Yes
74       Clinical                                  stroke              2                                        Unsure
75       Clinical                                     tia             -1                                    Not Stated
76       Clinical                                     tia              0                                            No
77       Clinical                                     tia              1                                           Yes
78       Clinical                                     tia              2                                        Unsure
79       Clinical                                     tbi             -1                                    Not Stated
80       Clinical                                     tbi              0                                            No
81       Clinical                                     tbi              1                                           Yes
82       Clinical                                     tbi              2                                        Unsure
83       Clinical                       recruited_control             -1                                    Not Stated
84       Clinical                       recruited_control              0                                            No
85       Clinical                       recruited_control              1                                           Yes
86       Clinical                       recruited_control              2                                       Pending
87       Clinical                                      ad             -1                                    Not Stated
88       Clinical                                      ad              0                                            No
89       Clinical                                      ad              1                                           Yes
90       Clinical                                      ad              2                                        Unsure
91       Clinical                                     ftd             -1                                    Not Stated
92       Clinical                                     ftd              0                                            No
93       Clinical                                     ftd              1                                           Yes
94       Clinical                                     ftd              2                                        Unsure
95       Clinical                                      pd             -1                                    Not Stated
96       Clinical                                      pd              0                                            No
97       Clinical                                      pd              1                                           Yes
98       Clinical                                      pd              2                                        Unsure
99       Clinical                                     als             -1                                    Not Stated
100      Clinical                                     als              0                                            No
101      Clinical                                     als              1                                           Yes
102      Clinical                                     als              2                                        Unsure
103      Clinical                                 mci_sci             -1                                    Not Stated
104      Clinical                                 mci_sci              0                                            No
105      Clinical                                 mci_sci              1                                           Yes
106      Clinical                                 mci_sci              2                                        Unsure
107      Clinical                                  cancer             -1                                    Not Stated
108      Clinical                                  cancer              0                                            No
109      Clinical                                  cancer              1                                           Yes
110      Clinical                                  cancer              2                                        Unsure
111      Clinical                                diabetes             -1                                    Not Stated
112      Clinical                                diabetes              0                                            No
113      Clinical                                diabetes              1                                           Yes
114      Clinical                                diabetes              2                                        Unsure
115      Clinical                                     chf             -1                                    Not Stated
116      Clinical                                     chf              0                                            No
117      Clinical                                     chf              1                                           Yes
118      Clinical                                     chf              2                                        Unsure
119      Clinical                                    copd             -1                                    Not Stated
120      Clinical                                    copd              0                                            No
121      Clinical                                    copd              1                                           Yes
122      Clinical                                    copd              2                                        Unsure
123      Clinical                                      mi             -1                                    Not Stated
124      Clinical                                      mi              0                                            No
125      Clinical                                      mi              1                                           Yes
126      Clinical                                      mi              2                                        Unsure
127      Clinical                                    afib             -1                                    Not Stated
128      Clinical                                    afib              0                                            No
129      Clinical                                    afib              1                                           Yes
130      Clinical                                    afib              2                                        Unsure
131      Clinical                                  angina             -1                                    Not Stated
132      Clinical                                  angina              0                                            No
133      Clinical                                  angina              1                                           Yes
134      Clinical                                  angina              2                                        Unsure
135      Clinical                         hyperlipidaemia             -1                                    Not Stated
136      Clinical                         hyperlipidaemia              0                                            No
137      Clinical                         hyperlipidaemia              1                                           Yes
138      Clinical                         hyperlipidaemia              2                                        Unsure
139      Clinical                            hypertension             -1                                    Not Stated
140      Clinical                            hypertension              0                                            No
141      Clinical                            hypertension              1                                           Yes
142      Clinical                            hypertension              2                                        Unsure
143      Clinical                              depression             -1                                    Not Stated
144      Clinical                              depression              0                                            No
145      Clinical                              depression              1                                           Yes
146      Clinical                              depression              2                                        Unsure
147      Clinical                 depression_test_battery           GADS
148      Clinical                 depression_test_battery            GDS
149      Clinical                 depression_test_battery           HADS
150      Clinical                 depression_test_battery            N/A
151      Clinical                                 anxiety             -1                                    Not Stated
152      Clinical                                 anxiety              0
153      Clinical                                 anxiety              1
154      Clinical                                 anxiety              2
155      Clinical                                     cdr             -1                                    Not Stated
156      Clinical                                     cdr              0
157      Clinical                                     cdr            0.5
158      Clinical                                     cdr              1
159      Clinical                                     cdr              2
160      Clinical                                     cdr              3
161      Clinical           computed_cognitive_test_score             -6             Verbal Refusal (Did not complete)
162      Clinical           computed_cognitive_test_score             -5              Other problem (Did not complete)
163      Clinical           computed_cognitive_test_score             -4 Cognitive/behavior problem (Did not complete)
164      Clinical           computed_cognitive_test_score             -3           Physical problem (Did not complete)
165      Clinical           computed_cognitive_test_score             -2                             Incompatible Data
166      Clinical           computed_cognitive_test_score             -1                                    Not Stated
167      Clinical           computed_cognitive_test_score              1                      Test Score within limits
168      Clinical                  cognitive_test_battery         ALSFRS
169      Clinical                  cognitive_test_battery            AH4
170      Clinical                  cognitive_test_battery            CAR     ALS Cognitive Behavioral Screen (ALS-CBS)
171      Clinical                  cognitive_test_battery           MMSE
172      Clinical                  cognitive_test_battery           MOCA
173      Clinical                  cognitive_test_battery            N/A
174      Clinical             computed_clinical_diagnosis             -1                                    Not Stated
175      Clinical             computed_clinical_diagnosis              0                                        Normal
176      Clinical             computed_clinical_diagnosis              1                                      Dementia
177      Clinical             computed_clinical_diagnosis              2                                           MCI
178      Clinical           computed_cognitive_impairment             -1                                    Not Stated
179      Clinical           computed_cognitive_impairment              0                                  Not Impaired
180      Clinical           computed_cognitive_impairment              1                                      Impaired
181      Clinical                            is_neuropath             -1                                    Not Stated
182      Clinical                            is_neuropath              0                                            No
183      Clinical                            is_neuropath              1                                           Yes
184      Clinical                            is_biomarker             -1                                    Not Stated
185      Clinical                            is_biomarker              0                                            No
186      Clinical                            is_biomarker              1                                           Yes
187 SomalogicMeta                                   units             -1                                           N/A
188 SomalogicMeta                                   units       log2 RFU
189 SomalogicMeta                               row_check             -1                                           N/A
190 SomalogicMeta                               row_check           FLAG
191 SomalogicMeta                               row_check           PASS
192 SomalogicMeta                           sample_matrix             -1                                           N/A
193 SomalogicMeta                           sample_matrix    EDTA Plasma
194 SomalogicMeta                           sample_matrix Citrate Plasma
195 SomalogicMeta                           sample_matrix            CSF
196 SomalogicMeta                           sample_matrix          Serum
197 SomalogicMeta                             sample_type         Sample
198 SomalogicMeta                             sample_type         Buffer
199 SomalogicMeta                             sample_type     Calibrator
200 SomalogicMeta                             sample_type             QC
201 PersonMapping                            is_somalogic              0                                            N0
202 PersonMapping                            is_somalogic              1                                           Yes
203 PersonMapping                            is_mass_spec              0                                            N0
204 PersonMapping                            is_mass_spec              1                                           Yes
205      Clinical             computed_clinical_diagnosis              3                                      Other ND

## PREVIOUSLY:
#         table_name                             column_name     key_lookup                               key_description
# 1         Clinical                      computed_age_range             -1                                           N/A
# 2         Clinical                      computed_age_range              0                                   Z=-1 to Z=1
# 3         Clinical                      computed_age_range              1                   Z=-2 to Z=-1 and Z=1 to Z=2
# 4         Clinical                      computed_age_range              2                   Z=-3 to Z=-2 and Z=2 to Z=3
# 5         Clinical                      computed_age_range              3                                   Z<3 and Z>3
# 6         Clinical                                     sex             -1                                    Not Stated
# 7         Clinical                                     sex              1                                          Male
# 8         Clinical                                     sex              2                                        Female
# 9         Clinical                                    race             -1                                    Not Stated
# 10        Clinical                                    race              1                 American_Indian/Alaska_Native
# 11        Clinical                                    race              2                     Black_or_African_American
# 12        Clinical                                    race              3     Native_Hawaiian_or_Other_Pacific_Islander
# 13        Clinical                                    race              4                                         Asian
# 14        Clinical                                    race              5                                         White
# 15        Clinical                                    race             66                               Other(Group_C)
# 16        Clinical                                    race             77               Non_Caucasian_or_Black(Group_B)
# 17        Clinical                                    race             88                        Non_Caucasian(Group_A)
# 18        Clinical                                    race             99                                       Unknown
# 19        Clinical          computed_years_education_range             -1                                           N/A
# 20        Clinical          computed_years_education_range              0                                   Z=-1 to Z=1
# 21        Clinical          computed_years_education_range              1                   Z=-2 to Z=-1 and Z=1 to Z=2
# 22        Clinical          computed_years_education_range              2                   Z=-3 to Z=-2 and Z=2 to Z=3
# 23        Clinical          computed_years_education_range              3                                   Z<3 and Z>3
# 24        Clinical                   computed_height_range             -2                             Incompatible data
# 25        Clinical                   computed_height_range             -1                                           N/A
# 26        Clinical                   computed_height_range              0                                   Z=-1 to Z=1
# 27        Clinical                   computed_height_range              1                   Z=-2 to Z=-1 and Z=1 to Z=2
# 28        Clinical                   computed_height_range              2                   Z=-3 to Z=-2 and Z=2 to Z=3
# 29        Clinical                   computed_height_range              3                                   Z<3 and Z>3
# 30        Clinical                   computed_weight_range             -1                                           N/A
# 31        Clinical                   computed_weight_range              0                                   Z=-1 to Z=1
# 32        Clinical                   computed_weight_range              1                   Z=-2 to Z=-1 and Z=1 to Z=2
# 33        Clinical                   computed_weight_range              2                   Z=-3 to Z=-2 and Z=2 to Z=3
# 34        Clinical                   computed_weight_range              3                                   Z<3 and Z>3
# 35        Clinical                      computed_bmi_range             -1                                           N/A
# 36        Clinical                      computed_bmi_range              0                                   Z=-1 to Z=1
# 37        Clinical                      computed_bmi_range              1                   Z=-2 to Z=-1 and Z=1 to Z=2
# 38        Clinical                      computed_bmi_range              2                   Z=-3 to Z=-2 and Z=2 to Z=3
# 39        Clinical                      computed_bmi_range              3                                   Z<3 and Z>3
# 40        Clinical                    computed_pulse_range             -1                                           N/A
# 41        Clinical                    computed_pulse_range              0                                   Z=-1 to Z=1
# 42        Clinical                    computed_pulse_range              1                   Z=-2 to Z=-1 and Z=1 to Z=2
# 43        Clinical                    computed_pulse_range              2                   Z=-3 to Z=-2 and Z=2 to Z=3
# 44        Clinical                    computed_pulse_range              3                                   Z<3 and Z>3
# 45        Clinical  computed_systolic_blood_pressure_range             -1                                           N/A
# 46        Clinical  computed_systolic_blood_pressure_range              0                                   Z=-1 to Z=1
# 47        Clinical  computed_systolic_blood_pressure_range              1                   Z=-2 to Z=-1 and Z=1 to Z=2
# 48        Clinical  computed_systolic_blood_pressure_range              2                   Z=-3 to Z=-2 and Z=2 to Z=3
# 49        Clinical  computed_systolic_blood_pressure_range              3                                   Z<3 and Z>3
# 50        Clinical computed_diastolic_blood_pressure_range             -2                             Incompatible data
# 51        Clinical computed_diastolic_blood_pressure_range             -1                                           N/A
# 52        Clinical computed_diastolic_blood_pressure_range              0                                   Z=-1 to Z=1
# 53        Clinical computed_diastolic_blood_pressure_range              1                   Z=-2 to Z=-1 and Z=1 to Z=2
# 54        Clinical computed_diastolic_blood_pressure_range              2                   Z=-3 to Z=-2 and Z=2 to Z=3
# 55        Clinical computed_diastolic_blood_pressure_range              3                                   Z<3 and Z>3
# 56        Clinical                              alcohol_hx             -1                                    Not Stated
# 57        Clinical                              alcohol_hx              0                                         Never
# 58        Clinical                              alcohol_hx              1                                       Current
# 59        Clinical                              alcohol_hx              2                                      Previous
# 60        Clinical                              alcohol_hx              3                             Never_or_Previous
# 61        Clinical                              smoking_hx             -1                                    Not Stated
# 62        Clinical                              smoking_hx              0                                         Never
# 63        Clinical                              smoking_hx              1                                       Current
# 64        Clinical                              smoking_hx              2                                      Previous
# 65        Clinical                              smoking_hx              3                             Never_or_Previous
# 66        Clinical                   computed_years_smoked             -2                             Incompatible data
# 67        Clinical                   computed_years_smoked             -1                                           N/A
# 68        Clinical                   computed_years_smoked              0                                   Z=-1 to Z=1
# 69        Clinical                   computed_years_smoked              1                   Z=-2 to Z=-1 and Z=1 to Z=2
# 70        Clinical                   computed_years_smoked              2                   Z=-3 to Z=-2 and Z=2 to Z=3
# 71        Clinical                   computed_years_smoked              3                                   Z<3 and Z>3
# 72        Clinical                                  stroke             -1                                    Not Stated
# 73        Clinical                                  stroke              0                                            No
# 74        Clinical                                  stroke              1                                           Yes
# 75        Clinical                                  stroke              2                                        Unsure
# 76        Clinical                                     tia             -1                                    Not Stated
# 77        Clinical                                     tia              0                                            No
# 78        Clinical                                     tia              1                                           Yes
# 79        Clinical                                     tia              2                                        Unsure
# 80        Clinical                                     tbi             -1                                    Not Stated
# 81        Clinical                                     tbi              0                                            No
# 82        Clinical                                     tbi              1                                           Yes
# 83        Clinical                                     tbi              2                                        Unsure
# 84        Clinical                       recruited_control             -1                                    Not Stated
# 85        Clinical                       recruited_control              0                                            No
# 86        Clinical                       recruited_control              1                                           Yes
# 87        Clinical                       recruited_control              2                                       Pending
# 88        Clinical                                      ad             -1                                    Not Stated
# 89        Clinical                                      ad              0                                            No
# 90        Clinical                                      ad              1                                           Yes
# 91        Clinical                                      ad              2                                        Unsure
# 92        Clinical                                     ftd             -1                                    Not Stated
# 93        Clinical                                     ftd              0                                            No
# 94        Clinical                                     ftd              1                                           Yes
# 95        Clinical                                     ftd              2                                        Unsure
# 96        Clinical                                      pd             -1                                    Not Stated
# 97        Clinical                                      pd              0                                            No
# 98        Clinical                                      pd              1                                           Yes
# 99        Clinical                                      pd              2                                        Unsure
# 100       Clinical                                     als             -1                                    Not Stated
# 101       Clinical                                     als              0                                            No
# 102       Clinical                                     als              1                                           Yes
# 103       Clinical                                     als              2                                        Unsure
# 104       Clinical                                 mci_sci             -1                                    Not Stated
# 105       Clinical                                 mci_sci              0                                            No
# 106       Clinical                                 mci_sci              1                                           Yes
# 107       Clinical                                 mci_sci              2                                        Unsure
# 108       Clinical                                  cancer             -1                                    Not Stated
# 109       Clinical                                  cancer              0                                            No
# 110       Clinical                                  cancer              1                                           Yes
# 111       Clinical                                  cancer              2                                        Unsure
# 112       Clinical                                diabetes             -1                                    Not Stated
# 113       Clinical                                diabetes              0                                            No
# 114       Clinical                                diabetes              1                                           Yes
# 115       Clinical                                diabetes              2                                        Unsure
# 116       Clinical                                     chf             -1                                    Not Stated
# 117       Clinical                                     chf              0                                            No
# 118       Clinical                                     chf              1                                           Yes
# 119       Clinical                                     chf              2                                        Unsure
# 120       Clinical                                    copd             -1                                    Not Stated
# 121       Clinical                                    copd              0                                            No
# 122       Clinical                                    copd              1                                           Yes
# 123       Clinical                                    copd              2                                        Unsure
# 124       Clinical                                      mi             -1                                    Not Stated
# 125       Clinical                                      mi              0                                            No
# 126       Clinical                                      mi              1                                           Yes
# 127       Clinical                                      mi              2                                        Unsure
# 128       Clinical                                    afib             -1                                    Not Stated
# 129       Clinical                                    afib              0                                            No
# 130       Clinical                                    afib              1                                           Yes
# 131       Clinical                                    afib              2                                        Unsure
# 132       Clinical                                  angina             -1                                    Not Stated
# 133       Clinical                                  angina              0                                            No
# 134       Clinical                                  angina              1                                           Yes
# 135       Clinical                                  angina              2                                        Unsure
# 136       Clinical                         hyperlipidaemia             -1                                    Not Stated
# 137       Clinical                         hyperlipidaemia              0                                            No
# 138       Clinical                         hyperlipidaemia              1                                           Yes
# 139       Clinical                         hyperlipidaemia              2                                        Unsure
# 140       Clinical                            hypertension             -1                                    Not Stated
# 141       Clinical                            hypertension              0                                            No
# 142       Clinical                            hypertension              1                                           Yes
# 143       Clinical                            hypertension              2                                        Unsure
# 144       Clinical                              depression             -1                                    Not Stated
# 145       Clinical                              depression              0                                            No
# 146       Clinical                              depression              1                                           Yes
# 147       Clinical                              depression              2                                        Unsure
# 148       Clinical                 depression_test_battery           GADS
# 149       Clinical                 depression_test_battery            GDS
# 150       Clinical                 depression_test_battery           HADS
# 151       Clinical                 depression_test_battery            N/A
# 152       Clinical                                 anxiety             -1                                    Not Stated
# 153       Clinical                                 anxiety              0                                            No
# 154       Clinical                                 anxiety              1                                           Yes
# 155       Clinical                                 anxiety              2                                        Unsure
# 156       Clinical           computed_cognitive_test_score             -6             Verbal Refusal (Did not complete)
# 157       Clinical           computed_cognitive_test_score             -5              Other problem (Did not complete)
# 158       Clinical           computed_cognitive_test_score             -4 Cognitive/behavior problem (Did not complete)
# 159       Clinical           computed_cognitive_test_score             -3           Physical problem (Did not complete)
# 160       Clinical           computed_cognitive_test_score             -2                             Incompatible Data
# 161       Clinical           computed_cognitive_test_score             -1                                    Not Stated
# 162       Clinical           computed_cognitive_test_score              1                      Test Score within limits
# 163       Clinical                  cognitive_test_battery         ALSFRS
# 164       Clinical                  cognitive_test_battery            CAR
# 165       Clinical                  cognitive_test_battery           MMSE
# 166       Clinical                  cognitive_test_battery           MOCA
# 167       Clinical                  cognitive_test_battery            N/A
# 168 ProteomicsMeta                                   units             -1                                           N/A
# 169 ProteomicsMeta                                   units       log2 RFU
# 170 ProteomicsMeta                               row_check             -1                                           N/A
# 171 ProteomicsMeta                               row_check           FLAG
# 172 ProteomicsMeta                               row_check           PASS
# 173 ProteomicsMeta                           sample_matrix    EDTA Plasma
# 174 ProteomicsMeta                           sample_matrix Citrate Plasma
# 175 ProteomicsMeta                           sample_matrix            CSF
# 176 ProteomicsMeta                             sample_type             -1                                           N/A
# 177 ProteomicsMeta                             sample_type         Sample
