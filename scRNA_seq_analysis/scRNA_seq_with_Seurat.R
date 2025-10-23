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

### data.dir = directory for 10X output files (Please refer to Seurt vignette)
seu.data <- Read10X(data.dir = "")
seu <- CreateSeuratObject(counts = seu.data, project = "pbmc3k", min.cells = 3, min.features = 200)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# with SCTransform
seu <- SCTransform(seu, vars.to.regress = "percent.mt", vst.flavor = "v2", method = "glmGamPoi")

# with standard log1p-normalization
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, vars.to.regress = "percent.mt")

seu <- RunPCA(seu)

seu <- RunUMAP(seu, dims = 1:10)
seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 0.5)

##### If you want to perform clustering with leiden algorithm, 
#library(reticulate)
#reticulate::install_python(version = "3.10")
#reticulate::py_install("leidenalg")
#reticulate::py_install("pandas")
seu <- FindClusters(seu, resolution = 0.5, algorithm = 4)
#####

DimPlot(seu, reduction = "umap", label = TRUE)
DimPlot(seu, reduction = "umap", group.by = "orig.ident", label = FALSE)

# with SCTransform
DefaultAssay(seu) <- "SCT"
Idents(seu) <- "seurat_clusters"

markers <- FindAllMarkers(seu, assay = "SCT", test.use = "wilcox", 
                          min.pct = 0.25, only.pos = TRUE)

# with standard log1p-normalization
DefaultAssay(seu) <- "RNA"
Idents(seu) <- "seurat_clusters"

markers <- FindAllMarkers(seu, test.use = "wilcox", 
                          min.pct = 0.25, only.pos = TRUE)

# Cell type annotation
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

seu$annotation <- Idents(seu)

# differentially expressed genes
# with SCTransform
deg_list <- FindMarkers(seu, assay = "SCT", test.use = "wilcox", 
                        ident.1 = "Memory CD4 T", 
                        ident.2 = "Naive CD4 T")

# with standard log1p-normalization
deg_list <- FindMarkers(seu, assay = "RNA", test.use = "wilcox", 
                        ident.1 = "Memory CD4 T", 
                        ident.2 = "Naive CD4 T")

deg_list <- deg_list[which(deg_list$p_val_adj < 0.05), ]
deg_list <- deg_list[which(abs(deg_list$avg_log2FC) > 0.5), ]