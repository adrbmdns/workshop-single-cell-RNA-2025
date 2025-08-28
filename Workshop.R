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

################## workshop 3 ########################
merged <- FindVariableFeatures(merged,
                               assay = "RNA",
                               nfeatures = 2000,
                               mean.cutoff = c(0.1, 8),
                               dispersion.cutoff = c(1, Inf))

top10 <- head(VariableFeatures(merged), 10)

VariableFeaturePlot(object = merged)

HoverLocator(VariableFeaturePlot(object = merged))

# PCA
merged <- RunPCA(merged, npcs = 50, assay = "RNA")

# Create elbow plot
elbow <- ElbowPlot(merged, ndims = 30)
jpeg("figures/Elbow.jpg", width = 8, height = 6, units = 'in', res = 150)
print(elbow)
dev.off()

# create heatmap
DimHeatmap(merged, dims = 1, cells = 500, nfeatures = 30,
           balanced = TRUE, fast = FALSE)

# PC 1 to 12
jpeg("figures/DimHm1_12.jpg", width = 10, height = 20,
     units = 'in', res = 150)
DimHeatmap(merged, dims = 1:12, balanced = TRUE, cells = 500)
dev.off()

# PC 13 to 24
jpeg("figures/DimHm13_24.jpg", width = 10, height = 20,
     units = 'in', res = 150)
DimHeatmap(merged, dims = 13:24, balanced = TRUE, cells = 500)
dev.off()

# PC 25 to 36
jpeg("figures/DimHm25_36.jpg", width = 10, height = 20,
     units = 'in', res = 150)
DimHeatmap(merged, dims = 25:36, balanced = TRUE, cells = 500)
dev.off()

# access SD values
length(which(merged@reductions$pca@stdev > 2))

# cell clustering
PC <- 22
merged <- FindNeighbors(merged, dims = 1:PC)
merged <- FindClusters(merged, resolution = 1.2,
                       cluster.name = "seurat_clusters_res1.2")

# Run UMAP
merged <- RunUMAP(merged, dims = 1:PC)
jpeg("figures/UMAP.jpg", width = 5, height = 4,
     units = "in", res = 150)
DimPlot(merged, label = TRUE, group.by = 'seurat_clusters_res1.2')
dev.off()

# color cells by sample

jpeg("figures/UMAP_by_sample.jpg", width = 5, height = 4, units = 'in',
     res = 150)
DimPlot(merged, label = FALSE, group.by = "orig.ident")
dev.off()

# color by one sample
highlighted_cells <- WhichCells(merged, expression = orig.ident == "Rep1_ICBdT")
DimPlot(merged, group.by = 'orig.ident', cells.highlight = highlighted_cells)

# compare resolution
merged <- FindClusters(merged, resolution = 0.8,
                       cluster.name = 'seurat_clusters_res0.8')
merged <- FindClusters(merged, resolution = 0.5,
                       cluster.name = 'seurat_clusters_res0.5')

# visualisation of different res
jpeg("figures/UMAP_compare_res.jpg", width = 20, height = 5, units = "in",
     res = 150)
DimPlot(merged, label = TRUE, group.by = 'seurat_clusters_res0.5') +
  DimPlot(merged, label = TRUE, group.by = 'seurat_clusters_res0.8') +
  DimPlot(merged, label = TRUE, group.by = 'seurat_clusters_res1.2')
dev.off()

# cell cycle scoring
library(gprofiler2)

s.genes <- gorth(cc.genes.updated.2019$s.genes, source_organism = 'hsapiens',
                 target_organism = 'mmusculus')$ortholog_name
g2m.genes <- gorth(cc.genes.updated.2019$g2m.genes, source_organism = 'hsapiens',
                   target_organism = 'mmusculus')$ortholog_name

merged <- CellCycleScoring(object = merged, s.features = s.genes,
                           g2m.features = g2m.genes)

FeaturePlot(merged, features = c("S.Score", "G2M.Score")) +
  DimPlot(merged, group.by = "Phase")

# save Seurat object
saveRDS(merged, file = "rep135_clustered.rds")

################# Workshop 4 #####################
library(SingleR)
library(celldex)
library(Seurat)
library(cowplot)

merged <- readRDS("rep135_clustered.rds")
ref_immgen <- ImmGenData()

# predict cell type with main labels
predictions_main <- SingleR(test = GetAssayData(merged),
                            ref = ref_immgen,
                            labels = ref_immgen$label.main)

predictions_fine <- SingleR(test = GetAssayData(merged),
                            ref = ref_immgen,
                            labels = ref_immgen$label.fine)

# plot delta distribution
plotDeltaDistribution(predictions_main)
plotDeltaDistribution(predictions_fine)
plotScoreHeatmap(predictions_main)

# add labels to Seurat metadata table
merged[["immgen_singler_main"]] <- rep('NA', ncol(merged))
merged$immgen_singler_main[rownames(predictions_main)] <-
  predictions_main$labels

merged[["immgen_singler_fine"]] <- rep('NA', ncol(merged))
merged$immgen_singler_fine[rownames(predictions_fine)] <-
  predictions_fine$labels

# visualise cell type composition within sample
library(viridis)
library(ggplot2)

ggplot(merged[[]], aes(x = orig.ident, fill = immgen_singler_main)) +
  geom_bar(position = "fill") +
  scale_fill_viridis(discrete = TRUE)

ggplot(merged[[]], aes(x = immgen_singler_main, fill = orig.ident)) +
  geom_bar(position = "fill") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_viridis(discrete = TRUE)

# map cell annotation to cluster
DimPlot(merged, group.by = c("immgen_singler_fine"), label = TRUE)

# use a different reference
ref_mouserna <- celldex::MouseRNAseqData()

# save file
saveRDS(merged, file = "preprocessed_object.rds")

##################### Workshop 5 ########################
library(Seurat)
library(dplyr)
library(EnhancedVolcano)
library(presto)

merged <- readRDS("preprocessed_object.rds")

# view cells with Epcam expression
FeaturePlot(merged, features = "Epcam")

# view clusters together
DimPlot(merged, group.by = "seurat_clusters_res0.8", label = TRUE) +
  FeaturePlot(merged, features = "Epcam") +
  DimPlot(merged, group.by = "immgen_singler_main")

# create a violin plot to visualise which cluster express Epcam
VlnPlot(merged, group.by = "seurat_clusters_res0.8", features = "Epcam")

# subset the Seurat object to only include clusters 8 and 13
merged <- SetIdent(merged, value = "seurat_clusters_res0.8")
merged_epithelial <- subset(merged, idents = c("8", "13"))

# plot the original and subset side-by-side
DimPlot(merged, group.by = "seurat_clusters_res0.8", label = TRUE) +
  DimPlot(merged_epithelial, group.by = "seurat_clusters_res0.8", label = TRUE)

# count the number of cells
table(merged$seurat_clusters_res0.8)
table(merged_epithelial$seurat_clusters_res0.8)

# find differential expressed genes in these clusters
merged_epithelial <- SetIdent(merged_epithelial, value = "seurat_clusters_res0.8")
epithelial_de <- FindMarkers(merged_epithelial, ident.1 = "8", ident.2 = "13",
                             min.pct = 0.25, logfc.threshold = 0.1)

# subset DE table with genes have a p-value less than 0.001
epithelial_de_sig <- epithelial_de[epithelial_de$p_val_adj < 0.001,]
epithelial_de_sig_top20 <- epithelial_de_sig |> top_n(n=20, wt=abs(avg_log2FC))

# get the top 20 gene names
epithelial_de_sig_top20_genes <- rownames(epithelial_de_sig_top20)
VlnPlot(merged_epithelial, features = epithelial_de_sig_top20_genes,
        group.by = "seurat_clusters_res0.8", ncol=5, pt.size=0)

# plot top 20 genes on UMAP
FeaturePlot(merged_epithelial, features = epithelial_de_sig_top20_genes, ncol=5)

# plot top 20 genes on a Dot Plot
DotPlot(merged_epithelial, features = epithelial_de_sig_top20_genes,
        group.by = "seurat_clusters_res0.8") +
  RotatedAxis()

# enhanced volcano
EnhancedVolcano(epithelial_de, lab = rownames(epithelial_de), x = "avg_log2FC",
                y = "p_val_adj", title = "Cluster 8 wrt 13", pCutoff = 0.05,
                FCcutoff = 0.5, pointSize = 3, labSize = 5, colAlpha = 0.3)

# re-run FindMarkers
epithelial_de_gsea <- FindMarkers(merged_epithelial, ident.1 = "8", ident.2 = "13",
                                  min.pct = 0.25, logfc.threshold = 0)
epithelial_de_gsea <- tibble::rownames_to_column(epithelial_de_gsea, var = "gene")
write.table(x = epithelial_de_gsea, file = "epithelial_de_all.tsv", sep = "\t",
            row.names = FALSE)

############### Pseudobulk DE Analysis ##################
pb_epithelial <- AggregateExpression(merged_epithelial, assays = 'RNA',
                                     return.seurat = TRUE,
                                     group.by = c('orig.ident',
                                                  'seurat_clusters_res0.8'))

# Have a look of the log-normalised table
pb_epithelial[['RNA']]$data

# change Ident to seurat clusters
Idents(pb_epithelial) <- "seurat_clusters_res0.8"

# find marker genes for comparing cluster 8 and 13
pb_epithelial_de <- FindMarkers(object = pb_epithelial, test.use = "DESeq2",
                                ident.1 = "8", ident.2 = "13")

# remove NA values
pb_epithelial_de <- na.omit(pb_epithelial_de)

# significant genes
pb_epithelial_de_sig <- pb_epithelial_de[pb_epithelial_de$p_val_adj < 0.001, ]

# compare result between single-cell and pseudobulk
# significant genes both present in single-cell and pseudobulk DE analysis
pb_and_sc_genes <- intersect(rownames(epithelial_de_sig),
                             rownames(pb_epithelial_de_sig))

# genes only significant in single-cell DE
only_sc_genes <- setdiff(rownames(epithelial_de_sig),
                         rownames(pb_epithelial_de_sig))

# genes only significant in pseudobulk DE
only_pb_genes <- setdiff(rownames(pb_epithelial_de_sig),
                         rownames(epithelial_de_sig))

################## DE Analysis on T cells ####################
# check all annotation types
unique(merged$immgen_singler_main)
t_celltypes_names <- c("T cells", "NKT", "Tgd")

# subset T cells
merged <- SetIdent(merged, value = "immgen_singler_main")
merged_tcells <- subset(merged, idents = t_celltypes_names)

# confirm the subset is correct
DimPlot(merged, group.by = "immgen_singler_main") +
  DimPlot(merged_tcells, group.by = "immgen_singler_main")

# create a new column in metadata table to store conditions
merged_tcells@meta.data$experimental_condition <- NA
merged_tcells@meta.data$experimental_condition[
  merged_tcells@meta.data$orig.ident %in% c("Rep1_ICB", "Rep3_ICB", "Rep5_ICB")
] <- "ICB"
merged_tcells@meta.data$experimental_condition[
  merged_tcells@meta.data$orig.ident %in% c("Rep1_ICBdT", "Rep3_ICBdT", "Rep5_ICBdT")
] <- "ICBdT"

# DE analysis between two groups
merged_tcells <- SetIdent(merged_tcells, value = "experimental_condition")
tcells_de <- FindMarkers(merged_tcells, ident.1 = "ICBdT", ident.2 = "ICB")
tcells_de_sig <- tcells_de[tcells_de$p_val_adj < 0.001, ]
tcells_de_sig |> top_n(n=5, wt=-avg_log2FC)
tcells_de_sig |> top_n(n=5, wt=avg_log2FC)

# inspect expression patterns for CD8/4 genes
VlnPlot(merged_tcells, features = c("Cd8a", "Cd8b1", "Cd4"))
merged_cd8tcells <- subset(merged_tcells, subset = Cd8b1>1 & Cd8a>1 & Cd4<0.1)

merged_cd8tcells <- SetIdent(merged_cd8tcells, value = "experimental_condition")
cd8tcells_de <- FindMarkers(merged_cd8tcells, ident.1 = "ICBdT", ident.2 = "ICB",
                            min.pct = 0.25)

cd8tcells_de_sig <- cd8tcells_de[cd8tcells_de$p_val_adj < 0.001, ]

cd8tcells_de_sig_top20 <- cd8tcells_de_sig |> top_n(n=20, wt=abs(avg_log2FC))
cd8tcells_de_sig_top20_genes <- rownames(cd8tcells_de_sig_top20)

# compare top 20 genes
VlnPlot(merged_cd8tcells, features = cd8tcells_de_sig_top20_genes[1:10],
        group.by = "experimental_condition", ncol=5, pt.size=0)
VlnPlot(merged_cd8tcells, features = cd8tcells_de_sig_top20_genes[11:20],
        group.by = "experimental_condition", ncol=5, pt.size=0)

FeaturePlot(merged_cd8tcells, features = cd8tcells_de_sig_top20_genes[1:10],
            ncol=5)
FeaturePlot(merged_cd8tcells, features = cd8tcells_de_sig_top20_genes[11:20],
            ncol=5)

DotPlot(merged_cd8tcells, features = cd8tcells_de_sig_top20_genes,
        group.by = "experimental_condition") + RotatedAxis()

EnhancedVolcano(cd8tcells_de, lab=rownames(cd8tcells_de), x='avg_log2FC',
                y='p_val_adj', title="ICBdT wrt ICB", pCutoff=0.05,
                FCcutoff = 0.5, pointSize = 3, labSize = 5, colAlpha = 0.3)

# save DE result to TSV files
cd8tcells_de_gsea <- FindMarkers(merged_cd8tcells, ident.1 = "ICBdT",
                                 ident.2 = "ICB", min.pct=0.25,
                                 logfc.threshold=0)
cd8tcells_de_gsea <- tibble::rownames_to_column(cd8tcells_de_gsea, var="gene")
write.table(cd8tcells_de_gsea, file="cd8tcells_de_gsea.tsv", sep='\t',
            row.names = FALSE)






