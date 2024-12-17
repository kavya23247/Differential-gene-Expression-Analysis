# Install BiocManager if not installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("airway")

BiocManager::install('DESeq2')

install.packages('tidyverse')

library(DESeq2)
library(tidyverse)
library(airway)

#-----------------Step 0: Downloading Data-----------------------------------
# Load the airway dataset
data(airway)
airway

# Extract sample information and modify
sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')

# Save sample information to a CSV file
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = TRUE, row.names = TRUE, quote = FALSE)

# Extract counts data and save to a CSV file
countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = TRUE, row.names = TRUE, quote = FALSE)

#------------------step 1: reading data--------------------------------------
## reading Downloaded counts data
counts_data <- read.csv('/content/counts_data.csv')
head(counts_data)
dim(counts_data)

## reading Downloaded sample info data
sample_info <- read.csv('/content/sample_info.csv')
head(sample_info)
dim(sample_info)

## Checking if column names in count data are matching to the row names in sample info data
all(colnames(counts_data) %in% rownames(sample_info))

## Checking if column names in count data are in same order to the row names in sample info data
all(colnames(counts_data) == rownames(sample_info))

#-------------------------step 2: Construct DESeq2 Dataset-------------------
## constructing the dataset
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = sample_info,
                              design= ~ dexamethasone)
dds

## prefiltering data : removing rows with low gene counts
## keeping genes with atleast 10 reads
## here we will only keep those genes which are having sum of counts in row greater than 10
keep <- rowSums(counts(dds))>=10
keep

dds <- dds[keep,]
dds

## in dexamethasone setting reference level as untreated
dds$dexamethasone <- relevel(dds$dexamethasone,ref = "untreated")
dds$dexamethasone

## building expression analysis using DESeq function
dds <- DESeq(dds)
dds

## saving results to res variable
res <- results(dds)
res

'''here we take p value as 0.05 which indicates that 5% of ur genes are not really differentially expressed,
but they are just occured due to random chance that means drug is having no real effect on 5% of genes.
that means from our dataset of ~22000 genes ~1100 genes are false positive(they are not differentially expressed 
but due to chance they are shown as differentially expressed'''
"in DESeq2 there are several several methods to adjust p value"

#--------------------step 3: Exploring Results------------------------------------

## getting summary of our results
summary(res) ## here adjusted p-value < 0.1 is considered

## we can adjust it by the following way
res.0.01 <- results(dds,alpha = 0.01)
summary(res.0.01)

## knowing what is compared against what
resultsNames(dds)

## here we cann compare results of dexamethasone untreated with different parameters
results(dds, contrast = c('dexamethasone','treated_4hrs','untreated'))

#--------------------step 4: Vizualizing Data------------------------------------

## to vizualize data we can make MA plot
plotMA(res) ## genes in blue are differentially expressed
