library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(sets)
library(ggpubr)
library(SeuratDisk)

options(Seurat.object.assay.version = "v5")

sample_list <- c("HC1120E", "DLE1113E", "SLE1116E")
seu.list <- as.list(sample_list)
names(seu.list) <- sample_list

for(id_sample in sample_list) {
  data_dir <- paste0("data_GSE179633_SLE/", id_sample, "/")
  seu.temp.data <- Read10X(data.dir = data_dir)
  seu.temp <- CreateSeuratObject(counts = seu.temp.data, project = id_sample, min.cells = 3, min.features = 200)
  seu.temp <- RenameCells(object = seu.temp, add.cell.id = id_sample)
  seu.list[[id_sample]] <- seu.temp
  print(id_sample)
}

seu <- merge(x = seu.list[[1]], y = seu.list[2:length(seu.list)], project = "SeuratProject")

### next step
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
Layers(seu[["RNA"]])

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mt < 10)

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
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- SCTransform(seu, vars.to.regress = "percent.mt", vst.flavor = "v2", method = "glmGamPoi")

seu <- RunPCA(seu)

seu.combined <- IntegrateLayers(
  object = seu, method = HarmonyIntegration, 
  orig.reduction = "pca", new.reduction = "harmony", 
  normalization.method = "SCT"
)

ElbowPlot(seu.combined)

### UMAP with Harmony
seu.combined <- RunUMAP(seu.combined, reduction = "harmony", dims = 1:15)
seu.combined <- FindNeighbors(seu.combined, reduction = "harmony", dims = 1:15)
seu.combined <- FindClusters(seu.combined, resolution = 0.5)

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
# (Write your code)

### Differentially expressed genes
# (Write your code)
