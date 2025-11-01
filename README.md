## BIO36101
### [Lecture] A single-cell sequencing guide for immunologists - BIO36101 Cell biology & Genetics Laboratory course @ UNIST
## scRNA_seq_analysis (Tutorial)
### 1. scRNA_seq_with_Seurat.R
* Retrieved & Modified R script from Seurat v5 PBMC 3k tutorial.
* To perform single-sample analysis.

### 2. scRNA_seq_with_Scanpy.ipynb
* Retrieved & Modified .ipynb file from Scanpy PBMC 3k tutorial.
* To perform single-sample analysis.

## scRNA_seq_analysis_multi (Report)
## R
### 1. scRNA_seq_multi_with_Harmony.R
* Main downstream analysis using Seurat.
* Doublet removal using scDblFinder.
* Normalization using SCTransform v2 with the Gamma-Poisson model fitting (glmGamPoi).
* Using 3000 highly variable genes (SCTransform default).
* Data scaling with regressing out percent.mt.
* Data integration using Harmony.
* Using 15 harmony-corrected dimensions.
* Using the Euclidean distance metric.
* To perform multi-sample analysis.

### 2. scDblFinder.R
* Running scDblFinder to perform doublet removal.

## Python
### 1. scRNA_seq_multi_with_Harmony.ipynb
* Main downstream analysis using Scanpy.
* Doublet removal using scDblFinder.
* Log1p-Normalization.
* Using 2000 highly variable genes.
* Data scaling with regressing out total_counts & pct_counts_mt.
* Data integration using Harmony.
* Using 10 harmony-corrected dimensions.
* Using the Euclidean distance metric.
* To perform multi-sample analysis.

### 2. command_scDblFinder_Python.txt
* To install scDblFinder package.

## References
### 1. (Tutorial) Data retrieved from 10X PBMC 3k.
### 2. (Report) Data retrieved from GSE179633.
