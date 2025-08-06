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
  cat("Working on sample:", sample, "\n")
  path <- paste("data/", sample, "-sample_filtered_feature_bc_matrix.h5", sep="")
  cat("Input file path:", path, "\n")
  data <- Read10X_h5(path)
  seurat_obj <- CreateSeuratObject(counts = data, project = sample,
                                   min.cells = 10, min.features = 100)
  sample.data[[sample]] <- seurat_obj
}

############ Workshop 2 ##################

# add mitochondrial RNA percentage
for (sample in sample_names) {
  sample.data[[sample]][['percent.mt']] <-
    PercentageFeatureSet(sample.data[[sample]], pattern = "^mt-", assay = "RNA")
}

# create violin plots for all samples
for (sample in sample_names) {
  jpeg(sprintf("%s/%s_unfilteredQC.jpg", "figures", sample),
       width = 16, height = 5, units = "in", res = 150)
  p1 <- VlnPlot(sample.data[[sample]], features = c("nCount_RNA"),
                pt.size = 0)
  p2 <- VlnPlot(sample.data[[sample]], features = c("nFeature_RNA"),
                pt.size = 0) +
    scale_y_continuous(breaks = c(0, 300, 500, 1000, 2000, 4000))
  p3 <- VlnPlot(sample.data[[sample]], features = c("percent.mt"),
                pt.size = 0) +
    scale_y_continuous(breaks = c(0, 12.5, 25, 50))
  p <- plot_grid(p1, p2, p3, ncol = 3)
  print(p)
  dev.off()
}

# write for loop to add labels
for (sample in sample_names) {
  sample.data[[sample]][["keep_cell_percent.mt"]] <-
    ifelse(sample.data[[sample]][["percent.mt"]] <= 12, TRUE, FALSE)

  sample.data[[sample]][["keep_cell_nFeature_RNA"]] <-
    ifelse(sample.data[[sample]][["nFeature_RNA"]] > 1000, TRUE, FALSE)
}

# print out number of cells after filtering
for (sample in sample_names) {
  print(sample)
  print(sum(sample.data[[sample]][["keep_cell_percent.mt"]] &
              sample.data[[sample]][["keep_cell_nFeature_RNA"]]))
}

# merge samples
unfiltered_merged <- merge(x = sample.data[["Rep1_ICBdT"]],
                           y = c(sample.data[["Rep3_ICBdT"]],
                                 sample.data[["Rep5_ICBdT"]],
                                 sample.data[["Rep1_ICB"]],
                                 sample.data[["Rep3_ICB"]],
                                 sample.data[["Rep5_ICB"]]),
                           add.cell.ids = sample_names)

merged <- subset(unfiltered_merged, nFeature_RNA > 1000 & percent.mt <=12) |>
  JoinLayers()

# run normalisation
merged <- NormalizeData(merged, assay = "RNA",
                        normalization.method = "LogNormalize",
                        scale.factor = 10000)

# scale data
merged <- ScaleData(merged, verbose = TRUE)
