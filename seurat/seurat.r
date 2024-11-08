library(sceasy)                                         
library(reticulate)
library(Seurat)
library(sctransform)

sceasy::convertFormat('4K/Mouse_brain_SC.h5ad', from="anndata", to="seurat", outFile='filename.rds')

ref <- readRDS("/home/ubuntu/bgi_projects/filename.rds")
ref

# file_path <- "annotation.txt"
file_path <- "/home/ubuntu/bgi_projects/annotation.txt"
lines <- readLines(file_path)
lines <- gsub("\n", "", lines)

Idents(ref) <- lines
Idents(ref)

# sceasy::convertFormat("PROBA2.h5ad", from="anndata", to="seurat", outFile='adata_st.rds')
sceasy::convertFormat("mouse_brain_visium_cell2location.h5ad", from="anndata", to="seurat", outFile='adata_st.rds')

st <- readRDS("adata_st.rds")
st

st <- SCTransform(st, assay = "RNA", ncells = 3000, verbose = FALSE)
st <- RunPCA(st)
st <- RunUMAP(st, dims = 1:30)
st <- FindNeighbors(st, dims = 1:30)
st <- FindClusters(st, resolution = 0.3, verbose = FALSE)
st <- FindVariableFeatures(st)
st <- ScaleData(st)

ref <- SCTransform(ref)
ref <- FindVariableFeatures(ref)
ref <- ScaleData(ref)

anchors <- FindTransferAnchors(reference = ref, query = st)

ref$celltypes <- Idents(ref)

predictions.assay <- TransferData(anchorset = anchors, refdata = ref$celltypes, prediction.assay = TRUE, weight.reduction = st[["pca"]], dims = 1:50)
st[["predictions"]] <- predictions.assay



DefaultAssay(st) <- "predictions"

options(max.print = 1e7)

sink("output.txt")
st[['predictions']]$'data'
sink()
