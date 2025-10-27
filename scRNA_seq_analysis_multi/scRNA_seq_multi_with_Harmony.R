library(dplyr)
library(Seurat)
library(patchwork)
library(DoubletFinder)
library(ggplot2)
library(cowplot)
library(sets)
library(ggpubr)
library(SeuratDisk)

options(Seurat.object.assay.version = "v5")

# data.dir = directory for 10X output files (Please refer to Seurt vignette)
seu1.data <- Read10X(data.dir = "")
seu2.data <- Read10X(data.dir = "")
seu3.data <- Read10X(data.dir = "")

seu.1 <- CreateSeuratObject(counts = seu1.data, project = "Normal_whole", min.cells = 3, min.features = 200)
seu.2 <- CreateSeuratObject(counts = seu2.data, project = "Normal_tomato_pos", min.cells = 3, min.features = 200)
seu.3 <- CreateSeuratObject(counts = seu3.data, project = "NASH_whole", min.cells = 3, min.features = 200)
seu.4 <- CreateSeuratObject(counts = seu4.data, project = "NASH_tomato_pos", min.cells = 3, min.features = 200)

### next step
seu <- merge(seu.1, y = c(seu.2, seu.3, seu.4), add.cell.ids = c("nowh", "notp", "nawh", "natp"), project = "SeuratProject")

seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
Layers(seu[["RNA"]])

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seu <- subset(seu, subset = nFeature_RNA > 300 & percent.mt < 10)

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seu.list <- SplitObject(seu, split.by = "orig.ident")

### Doublet removal using scDblFinder
# scDblFinder
seu.list <- lapply(X = seu.list, FUN = function(x) {
  x <- DoubletDetect(x)
})

### next step
seu <- merge(x = seu.list[[1]], y = seu.list[2:length(seu.list)], project = "SeuratProject")

seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
Layers(seu[["RNA"]])

seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)

##### Integration part
### Harmony with SCTransform
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
seu <- SCTransform(seu, vars.to.regress = "percent.mt", vst.flavor = "v2", method = "glmGamPoi")

seu <- RunPCA(seu)

seu.combined <- IntegrateLayers(
  object = seu, method = HarmonyIntegration, 
  orig.reduction = "pca", new.reduction = "harmony", 
  normalization.method = "SCT"
)

### UMAP with Harmony
sig.dims <- SelectOptimalDimension(seu.combined, jackstraw = FALSE)

seu.combined <- RunUMAP(seu.combined, reduction = "harmony", dims = 1:sig.dims)
seu.combined <- FindNeighbors(seu.combined, reduction = "harmony", dims = 1:sig.dims)
seu.combined <- FindClusters(seu.combined, resolution = 0.3)

# If you want to perform clustering with leiden algorithm, 
#library(reticulate)
#reticulate::install_python(version = "3.10")
#reticulate::py_install("leidenalg")
#reticulate::py_install("pandas")
seu.combined <- FindClusters(seu.combined, resolution = 0.5, algorithm = 4)
#

DimPlot(seu.combined, reduction = "umap", label = TRUE)
DimPlot(seu.combined, reduction = "umap", group.by = "orig.ident", label = FALSE)

### Find marker genes 
# with SCTransform
DefaultAssay(seu.combined) <- "SCT"
seu.combined <- PrepSCTFindMarkers(seu.combined)

Idents(seu.combined) <- "seurat_clusters"

markers <- FindAllMarkers(seu.combined, assay = "SCT", test.use = "wilcox", 
                          min.pct = 0.25, only.pos = TRUE)

### Celltype annotation


### Differentially expressed genes

