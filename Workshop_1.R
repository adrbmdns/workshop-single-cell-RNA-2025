library(Seurat)
library(ggplot2)
library(cowplot)

getwd()

# Read in matrix file
Rep1_ICBdT_data <- Read10X_h5("data/Rep1_ICBdT-sample_filtered_feature_bc_matrix.h5")

# Create a Seurat object
Rep1_ICBdT_data_seurat_obj <- CreateSeuratObject(counts = Rep1_ICBdT_data,
                                                 project = "Rep1_ICBdT",
                                                 min.cells = 10,
                                                 min.features = 100)

# Do the same for sample Rep1_ICB
Rep1_ICB_data <- Read10X_h5("data/Rep1_ICB-sample_filtered_feature_bc_matrix.h5")
Rep1_ICB_data_seurat_obj <- CreateSeuratObject(counts = Rep1_ICB_data,
                                               project = "Rep1_ICB",
                                               min.cells = 10,
                                               min.features = 100)

# print the Seurat object
print(Rep1_ICB_data_seurat_obj)
Assays(Rep1_ICB_data_seurat_obj)
Layers(Rep1_ICB_data_seurat_obj)

Rep1_ICB_data_seurat_obj[["RNA"]]
Rep1_ICB_data_seurat_obj[["RNA"]]$counts

head(Rep1_ICB_data_seurat_obj@meta.data)
head(Rep1_ICB_data_seurat_obj[[]])

head(Rep1_ICB_data_seurat_obj@meta.data$orig.ident)
head(Rep1_ICB_data_seurat_obj[['orig.ident']])

Rep1_ICB_data_seurat_obj[["percent.mt"]] <-
  PercentageFeatureSet(Rep1_ICB_data_seurat_obj,
                       pattern = "^mt-",
                       assay = "RNA")

# violin plot to look at quality of data
p1 <- VlnPlot(Rep1_ICB_data_seurat_obj,
              layer = "counts",
              features = c("nCount_RNA"),
              pt.size = 0) +
  scale_y_continuous(breaks = c(0,5000,10000,15000,20000,25000,30000))
p1

p2 <- VlnPlot(Rep1_ICB_data_seurat_obj,
              layer = "counts",
              features = c("nFeature_RNA"),
              pt.size = 0) +
  scale_y_continuous(breaks = c(0,1000,2000,3000,4000,5000,6000))
p2

p3 <- VlnPlot(Rep1_ICB_data_seurat_obj,
              layer = "counts",
              features = c("percent.mt"),
              pt.size = 0) +
  scale_y_continuous(breaks = c(0,5,10,15,20))
p3

p <- plot_grid(p1, p2, p3, ncol = 3)
p

sample_names <- c("Rep1_ICBdT", "Rep3_ICBdT", "Rep5_ICBdT",
                  "Rep1_ICB", "Rep3_ICB", "Rep5_ICB")

sample.data <- list()

for (sample in sample_names) {
  cat("\nWorking on sample: ", sample, "\n")
  path <- paste("data/", sample, "-sample_filtered_feature_bc_matrix.h5", sep="")
  cat("Input file path: ", path, "\n")
  data <- Read10X_h5(path)
  seurat_obj <- CreateSeuratObject(counts = data, project = sample,
                                   min.cells = 10, min.features = 100)
  sample.data[[sample]] <- seurat_obj
}
