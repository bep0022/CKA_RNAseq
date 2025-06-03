######## Script 4: Normalizing the HiSAT2 Mapping Counts for analysis #####

## PURPOSE: to create .csv files of averaged normalized mapping counts to assess the effects of age and sex
##  the .csv files will be assessing the number and identity of genes between age and sex groups.

## INPUT FILES: 
##  the gene count matrix .csv file generated after running the mapping and stringtie scripts... 
##    in which the different sample IDs are in columns and genes in rows (genexID)
##  a .txt metadata file containing the corresponding age/sex information for the given sample
##  IDs... in which the columns are age and sex and rows are the samples IDs(IDxAge/Sex)

## OUTPUT FILES: 
##  a .csv file containing the averaged normalized mapping counts for all of the individuals of the same sex within an age 

## SCRIPT VARIABLE DESCRIPTIONS:
##  Counts - table containing the DGE count data for your samples (genes x samples)
##  d0 - object containing the information from the DGE Counts variable
##  d - object containing the information from the DGE Counts variable with normalization factors calculated
##    d is the object that will be used for creating plots and downstream analyses 
##  phenodata - metadata file containing the treatment information for the corresponding samples
##    IN THE SAME ORDER! i.e. the order in the metadata correctly matches the order of the Counts data file


#### LOAD PACKAGES AND INPUT FILES ####

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")

library(edgeR)
library(limma)
library(gplots)
library(ggplot2)

### Read in the input .csv counts table from HiSAT2xStringtie output
###   row.names = 1 acknowledges that the rows have names that are stored in the 1st column 
###   Rows = gene names ; ex. gene-LOC124189989
###   Columns = individual samples ex. CKA_004
counts <- read.csv("Output_Data/HiSAT2/gene_count_matrix.csv", row.names = 1)
head(counts)
dim(counts) # rows x columns    ex. 30359 genes x 18 samples 


####### SET VARIABLES FOR DATA AND METADATA FILES #####

## Create DGEList object
#   DGEList - Creates a DGEList object from a table of counts (rows=features, columns=samples), 
#     group indicator for each row, library size (optional) and a table of feature annotation (optional).
d0 <- DGEList(counts)
d0

## Calculate normalization factors for use downstream (does not actually normalize the data at this step, just stores the 
##    normalization factors calculated based on the TMM method in the norm.factors column within the d0 DGElist variable)
d0 <- calcNormFactors(d0)
d0

#   Ex. BEFORE calculation of norm.factors: 
#           group   lib.size    norm.factors
#CKA_004     1      52111082            1

#   Ex. AFTER calculation of norm.factors 
#           group    lib.size   norm.factors
# CKA_004     1     52111082      0.7760032 


## Filter low-expressed genes
cutoff <- 1 
drop <- which(apply(cpm(d0), 1, max) < cutoff) # 1 indicates rows ; takes the max value of the cpm for a given gene (row) 
# across the different samples and sets those that are less than 1 to be input into variable "drop"
d <- d0[-drop,] 
dim(d0) # number of genes prior to filtering out low-expressed genes
dim(d) # number of genes left after filtering low-expressed genes (< 1 cpm)
d$samples
head(d$counts)
## Sample names    
snames <- colnames(counts)
snames
# Ex. "CKA_004" "CKA_007" "CKA_015" "CKA_018" "CKA_024" "CKA_025"


## After normalizing the counts we want to average all of the counts for the 
#   individuals within the same sex for each age.
## Input the phenotype data
# Make sure the individual names match between the count data and the metadata!
pheno <-(read.csv("Phenodata.csv", header=TRUE, row.names=1))
dim(pheno)
head(pheno)


#### SUBSET COUNTS DATA ####
# 7 month females
counts_7F_subset <- counts[, c("CKA_113", "CKA_109", "CKA_106", "CKA_111", "CKA_110")]
counts_7F_subset$`7F_Avg` <- rowMeans(counts_7F_subset, na.rm = TRUE)
nrow(counts_7F_subset)

# 7 month males
counts_7M_subset <- counts[, c("CKA_107", "CKA_122", "CKA_125", "CKA_108", "CKA_112")]
counts_7M_subset$`7M_Avg` <- rowMeans(counts_7M_subset, na.rm = TRUE)
nrow(counts_7M_subset)

# 18 month females
counts_18F_subset <- counts[, c("CKA_079", "CKA_081", "CKA_077", "CKA_078", "CKA_065")]
counts_18F_subset$`18F_Avg` <- rowMeans(counts_18F_subset, na.rm = TRUE)
nrow(counts_18F_subset)

# 18 month males
counts_18M_subset <- counts[, c("CKA_070", "CKA_073", "CKA_075", "CKA_063", "CKA_076")]
counts_18M_subset$`18M_Avg` <- rowMeans(counts_18M_subset, na.rm = TRUE)
nrow(counts_18M_subset)

# 36 month females
counts_36F_subset <- counts[, c("CKA_018", "CKA_033", "CKA_015", "CKA_004", "CKA_007")]
counts_36F_subset$`36F_Avg` <- rowMeans(counts_36F_subset, na.rm = TRUE)
nrow(counts_36F_subset)

# 36 month males
counts_36M_subset <- counts[, c("CKA_027", "CKA_024", "CKA_026", "CKA_028", "CKA_025")]
counts_36M_subset$`36M_Avg` <- rowMeans(counts_36M_subset, na.rm = TRUE)
nrow(counts_36M_subset)

# 48 month females
counts_48F_subset <- counts[, c("CKA_051", "CKA_047", "CKA_048", "CKA_046", "CKA_055")]
counts_48F_subset$`48F_Avg` <- rowMeans(counts_48F_subset, na.rm = TRUE)
nrow(counts_48F_subset)

# 48 month males
counts_48M_subset <- counts[, c("CKA_053", "CKA_043", "CKA_058", "CKA_049", "CKA_054")]
counts_48M_subset$`48M_Avg` <- rowMeans(counts_48M_subset, na.rm = TRUE)
nrow(counts_48M_subset)

# 60 month females
counts_60F_subset <- counts[, c("CKA_135", "CKA_130", "CKA_136", "CKA_138", "CKA_132")]
counts_60F_subset$`60F_Avg` <- rowMeans(counts_60F_subset, na.rm = TRUE)
nrow(counts_60F_subset)

# 60 month males
counts_60M_subset <- counts[, c("CKA_131", "CKA_127", "CKA_128")]
counts_60M_subset$`60M_Avg` <- rowMeans(counts_60M_subset, na.rm = TRUE)
nrow(counts_60M_subset)


#### MERGE DATAFRAMES ####
# Merge all of the dataframes that were subset by age x sex 
merged_counts <- cbind(counts_60F_subset, 
                       counts_60M_subset, 
                       counts_48F_subset, 
                       counts_48M_subset,
                       counts_36F_subset,
                       counts_36M_subset,
                       counts_18F_subset,
                       counts_18M_subset,
                       counts_7F_subset,
                       counts_7M_subset)
View(merged_counts)


#### FINALIZE DATAFRAME ####
# Remove all of the columns starting with "CKA" - the non-averaged counts
merged_counts_no_CKA <- merged_counts[, !grepl("^CKA", colnames(merged_counts))]
View(merged_counts_no_CKA)
nrow(merged_counts_no_CKA)

write.csv(merged_counts_no_CKA, "avg_normalized_gene_count_matrix.csv", row.names = TRUE)




