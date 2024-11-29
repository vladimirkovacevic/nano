library(Seurat)
library(scCustomize)
library("ggplot2")
options(future.globals.maxSize = 8000 * 1024^2)

data.ref <- readRDS("seurat_object.rds")

data.ref <- NormalizeData(data.ref, normalization.method = "LogNormalize", scale.factor = 1e3)
data.ref <- FindVariableFeatures(data.ref)
data.ref <- ScaleData(data.ref)
data.ref <- RunPCA(data.ref)
data.ref <- RunUMAP(data.ref, dims = 1:30, return.model = TRUE)

samples_path = "RNA"
samples = list.files(path=samples_path, pattern="^Uzorak_.*")
for (sample in samples){
  path = file.path(samples_path, sample, "uzorak")
  
  sampleData <- Read10X(data.dir = path)
  data.query <- CreateSeuratObject(counts = sampleData, project = sample, min.cells = 3, min.features = 200)
  
  data.query <- NormalizeData(data.query, normalization.method = "LogNormalize", scale.factor = 1e3)
  #data.query <- FindVariableFeatures(data.query)
  data.query <- ScaleData(data.query)
  data.query <- FindVariableFeatures(data.query)
  data.query <- RunPCA(data.query)
  data.query <- RunUMAP(data.query, dims = 1:30, return.model = TRUE)
  data.query <- FindNeighbors(data.query, dims = 1:30)
  data.query <- FindClusters(data.query, resolution = 0.3, verbose = FALSE)
  data.query <- FindVariableFeatures(data.query)
  data.query <- ScaleData(data.query)
  print(path)
  anchors <- FindTransferAnchors(reference = data.ref, query = data.query, reference.reduction = "pca")

  predictions <- TransferData(anchorset = anchors, refdata = data.ref$CellType, prediction.assay = FALSE, weight.reduction = data.query[["pca"]], dims = 1:50)

  data.query <- AddMetaData(data.query, metadata = predictions)

  table(data.query$predicted.id)
  DimPlot(data.query, reduction = 'umap', group.by="predicted.id")

  data.query <- MapQuery(anchorset = anchors, reference = data.ref, query = data.query, refdata = list(celltype = "CellType"), reference.reduction = "pca", reduction.model = "umap")

  p1 <- DimPlot(data.ref, reduction = "umap", group.by = "CellType", label = TRUE, label.size = 3, repel = TRUE) + ggtitle("Reference annotations")
  p2 <- DimPlot(data.query, reduction = "ref.umap", group.by = "predicted.id", label = TRUE, label.size = 3, repel = TRUE) + ggtitle("Query transferred labels")
  p1 + p2

  options(max.print = 1e7)

  table(data.query$predicted.id)
  write.csv(as.matrix(data.query$'predicted.celltype'), file=paste(sample, ".csv", sep=""))
  as.anndata(x = data.query, file_path = getwd(), file_name = paste(sample, ".h5ad", sep=""))
}
