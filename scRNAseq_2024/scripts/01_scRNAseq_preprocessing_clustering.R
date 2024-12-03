###########################################################
# This script is intended to recapitulate the working     #
# working session in course: "scRNAseq data analysis"     #
# which was held on Nov/ 04-05, 2024 by Juan Jovel        #
#                                                         #
# Use it at your own risk                                 #
###########################################################

# load required libraries
library(Seurat)    # main library for the analysis of scRNAseq data
library(tximport)  # used to import text files produced by salmon alevin
library(tidyverse) # includes several R packages needed
library(patchwork)

# set the working directory
setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/juanJovel/courses/2024/scRNAseq/scRNAseq_course2024_materials/alevin_quants')

# define the path of the alevin directory
alevin_dir = "pbmc3k_alevin_output"

# File path to the 'quants_mat.gz file, where the CB and UMIs lists are too
files <- file.path(alevin_dir, "alevin", "quants_mat.gz")
files

# import nedeed files with tximport
txi <- tximport(files, type = "alevin")
  
# extract counts data
data <- txi$counts
  
# Create a Seurat object per dataset
obj <- CreateSeuratObject(counts = data, project = "pbmc3k")

## QC ###
# Calculate nFeature_RNA and nCount_RNA if they don't exist
obj <- AddMetaData(obj, metadata = Matrix::colSums(obj@assays$RNA$counts > 0), col.name = "nFeature_RNA")
obj <- AddMetaData(obj, metadata = Matrix::colSums(obj@assays$RNA$counts), col.name = "nCount_RNA")

# Calculate the percentage of mitochondrial reads
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# Function to generate QC violin plots
generateQCPlot <- function(obj, barcode, suffix = "") {
  filename <- paste0(barcode, '_seurat_QCplot', suffix, '.png')
  png(filename, width = 12, height = 9, units = 'in', pointsize = 24, res = 300)
  p <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, cols = 'dodgerblue')
  print(p)
  dev.off()
}

# Function to generate feature relationship scatter plots
generateFeatureRelationshipPlot <- function(obj, barcode, suffix = "") {
  filename <- paste0(barcode, '_seurat_feature-relationship', suffix, '.png')
  png(filename, width = 12, height = 9, units = 'in', pointsize = 24, res = 300)
  plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot1 + plot2)
  dev.off()
}

# Generate initial QC plot and feature relationship plot
generateQCPlot(obj, "pbmc3k")
generateFeatureRelationshipPlot(obj, "pbmc3k")
  
# Filter cells based on QC criteria
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 & 
                  nCount_RNA > 500 & nCount_RNA < 10000)
  
# Generate post-filtering QC plot and feature relationship plot
generateQCPlot(obj, "pbmc3k", suffix = "_AF")
generateFeatureRelationshipPlot(obj, "pbmc3k", suffix = "_AF")
  
# Normalize counts in each of the datasets
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)
  
# Function to save plots
saveSinglePlot <- function(plot, filename){
  png(filename, width = 8, height = 6, units = 'in',  res = 300)
  print(plot)
  dev.off()
}

plot <- VizDimLoadings(obj, dims = 1:2, reduction = "pca")
saveSinglePlot(plot, 'VizDimLoadings_pbmc3k.png')

plot <- DimHeatmap(obj, dims = 1:8, cells = 500, balanced = TRUE)
ggsave(filename = "DimHeatmap_pbmc3k.png", plot = plot, width = 10, height = 8, dpi = 300)

plot <- ElbowPlot(obj)
saveSinglePlot(plot, 'ElbowPlot_pbmc3k.png')

########## UMAP ##########
obj <- RunUMAP(obj, reduction = "pca", dims = 1:15)

########## tSNE ##########
# Remove duplicates
obj2 <- obj[["pca"]]@cell.embeddings[,1:13]
duplicate_cells <- duplicated(obj2)
obj_unique <- obj[,!duplicate_cells]

# Now run t-SNE on the deduplicated data
obj_unique <- RunTSNE(obj,
                           dims = 1:13,
                           perplexity = 50,
                           max_iter = 1000,
                           seed.use = 42)

# Clustering
obj <- FindNeighbors(obj, dims = 1:15)
obj <- FindClusters(obj, resolution = 1)

obj_unique <- FindNeighbors(obj_unique, dims = 1:15)
obj_unique <- FindClusters(obj_unique, resolution = 1)


# Visual comparison plots
p1 <- DimPlot(obj, reduction = "pca", group.by = "seurat_clusters", 
              dims = c(1,2)) + ggtitle("PCA")
p2 <- DimPlot(obj_unique, reduction = "tsne",  group.by = "seurat_clusters", 
              label = TRUE, label.size = 3) + ggtitle("tSNE")
p3 <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", 
              label = TRUE, label.size = 3) + ggtitle("UMAP")

# Visualize UMAP
filename <- "PCA_tSNE_UMAP_integrated-object.png"
png(filename, width = 18, height = 7, units = 'in', pointsize = 24, res = 300)
combined_plot <- p1 + p2 + p3 + plot_layout(ncol = 3)
print(combined_plot)
dev.off()

# Save in RDS format
saveRDS(obj, 'pbmc3k.rds')

