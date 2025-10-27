DoubletDetect <- function(Object) {
  library(scDblFinder)
  sce <- as.SingleCellExperiment(Object)
  sce <- scater::logNormCounts(sce)
  
  sce$cluster <- fastcluster(sce)
  sce <- scDblFinder(sce, clusters = "cluster", nfeatures = 2000)
  
  seu <- as.Seurat(sce, counts = "counts", data = NULL)
  seu <- subset(seu, subset = scDblFinder.class == "singlet")
  
  seu <- CreateSeuratObject(counts = seu@assays$RNA@counts, project = as.character(Object$orig.ident)[1], min.cells = 3, min.features = 200)
  seu$orig.ident <- rep(as.character(Object$orig.ident)[1], length(seu$orig.ident))
  
  return(seu)
}

