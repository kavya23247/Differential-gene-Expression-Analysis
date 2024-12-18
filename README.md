# RNA-Seq Differential Gene Expression Analysis Workflow

## Step 0: Downloading and Preparing Data

### 1. Installing and Loading Libraries
**BiocManager**: A package manager for Bioconductor packages.
**DESeq2**: Primary package for differential gene expression (DGE) analysis.
**tidyverse**: Utilities for data manipulation and visualization.
**airway**: A real RNA-Seq dataset used for practice.

### 2. Loading the Dataset
The airway dataset contains gene expression count data from airway smooth muscle cells.
The experiment involves treated and untreated samples with dexamethasone (a corticosteroid).

### 3. Extracting and Modifying Metadata
Extract sample metadata including treatment condition and cell line.
Rename "trt" and "untrt" to "treated" and "untreated" for clarity.
Save the metadata as sample_info.csv.

### 4. Extracting Count Data
Retrieve raw count data of genes across samples.
Save the counts matrix as counts_data.csv.

---

## Step 1: Reading Data

### 1. Reading Counts and Metadata
Load the previously saved counts data and sample information.

### 2. Verifying Matching Columns
Ensure that the column names in the counts matrix match the row names in the sample metadata.
Proper alignment is necessary for DESeq2 to perform DGE analysis.

---

## Step 2: Constructing the DESeq2 Dataset

### 1. Creating the DESeq2 Dataset
Combine gene expression counts and metadata into a DESeq2 object.
Specify the experimental design to compare "treated" vs "untreated" conditions.

### 2. Prefiltering Low Counts
Remove genes with low counts (less than 10 reads) across all samples.
Improves statistical power by reducing noise from non-informative genes.

### 3. Setting Reference Level
Set "untreated" as the reference level for comparison.
Results are reported as "treated relative to untreated."

### 4. Running DESeq
Run the DESeq2 pipeline:
  - Normalize for sequencing depth.
  - Estimate gene expression variability (dispersions).
  - Fit a Negative Binomial Generalized Linear Model (GLM).

---

## Step 3: Exploring Results

### 1. Extracting Results
Extract the results table with metrics such as:
  - **log2FoldChange**: Indicates up/down-regulation of genes.
  - **p-value**: Tests statistical significance.
  - **Adjusted p-value**: Corrects for multiple testing.

### 2. Adjusting P-value Threshold
Apply a stricter threshold (e.g., 0.01) for identifying significantly expressed genes.

### 3. Comparing Results
Specify custom contrasts for comparison between "treated" and "untreated" conditions.

---

## Step 4: Visualizing Results

### 1. MA Plot
Visualize gene expression changes:
  - X-axis: Average expression (mean).
  - Y-axis: log2FoldChange.
  - Highlight significantly differentially expressed genes.

---

## Key Insights from the Workflow

1. **Normalization**: Accounts for sequencing depth differences to ensure fair comparisons.
2. **Statistical Testing**: Uses robust Negative Binomial GLM for significance testing.
3. **Prefiltering**: Removes low-count genes to reduce noise and improve analysis quality.
4. **Visualization**: MA plots provide an intuitive representation of gene expression changes.

This workflow enables accurate identification of genes that respond to a treatment (e.g., dexamethasone) and ensures reproducibility in RNA-Seq experiments.
