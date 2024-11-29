library(Seurat)


file_path <- "C:\\Users\\zivic\\Desktop\\scrna\\SCP424\\metadata\\meta.txt"

data <- read.delim(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

colnames(data) <- as.character(data[1, ])
column_types <- as.character(data[2, ])

data <- data[-c(1, 2), ]

#colnames(data) <- as.character(data[1, ])

data <- data[-1, ]

for (i in 1:ncol(data)) {
  col_type <- column_types[i]
  if (col_type == "integer") {
    data[[i]] <- as.integer(data[[i]])
  } else if (col_type == "numeric") {
    data[[i]] <- as.numeric(data[[i]])
  } else if (col_type == "character") {
    data[[i]] <- as.character(data[[i]])
  }else if (col_type == "group") {
    data[[i]] <- as.factor(data[[i]])
  }
}

expression_matrix <- Read10X(data.dir = "C:\\Users\\zivic\\Desktop\\scrna\\SCP424\\uzorak")

seurat_object <- CreateSeuratObject(
  counts = expression_matrix,
  project = "MyProject",
  meta.data = data
)


seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)

saveRDS(seurat_object, file = "seurat_object.rds")