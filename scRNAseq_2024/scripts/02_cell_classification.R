##########################################################
# This script is intended to recapitulate the working     #
# working session in course: "scRNAseq data analysis"     #
# which was held on Nov/ 04-05, 2024 by Juan Jovel        #
#                                                         #
# Use it at your own risk                                 #
###########################################################

# load required libraries
library(Seurat)     # main library for the analysis of scRNAseq data
library(Azimuth)    # to automatically classify cells
library(tidyverse)  # multi-package for data manipulation and deployment
library(SeuratData) # to visualize profiles availble in Azimuth
library(ggrepel)    # to add labels that do not overlap
library(cowplot)    #  Accesory package to arrange multiple ggplot2 plots

setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/me/courses/2024/scRNAseq/scRNAseq_course2024_materials/alevin_quants')

# Load your h5Seurat file
pbmc3k <- readRDS("pbmc3k.rds")

# Let's recapitulate the Seurat clusters derived in our previous session
DimPlot(pbmc3k, reduction = "umap", group.by = "RNA_snn_res.1",
        label = TRUE, label.size = 3) + ggtitle("UMAP")

# Which references are available in Azimuth
available_data <- AvailableData()

# Run Azimuth
pbmc3k <- RunAzimuth(pbmc3k, reference = 'pbmcref')

DimPlot(pbmc3k, reduction = "umap", group.by = "predicted.celltype.l1",
        label = TRUE, label.size = 3) + ggtitle("UMAP")

DimPlot(pbmc3k, reduction = "umap", group.by = "predicted.celltype.l2",
        label = TRUE, label.size = 3) + ggtitle("UMAP")

library(pals)

cell_types <- unique(pbmc3k@meta.data$predicted.celltype.l2)
cell_types
num_cell_types <- length(cell_types)
num_cell_types
# Generate a polychrome palette with as many colors as there are cell types
palette <- polychrome(n = num_cell_types)

# Assign names to the palette for clarity (optional but recommended)
names(palette) <- cell_types

# Create the initial DimPlot without labels
dim_plot <- DimPlot(pbmc3k, 
                    reduction = "umap", 
                    group.by = "predicted.celltype.l2",
                    label = FALSE,  # Disable default labels
                    cols = palette) + 
  ggtitle("UMAP - Predicted Cell Types") +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title

# Extract UMAP embeddings and cluster assignments
umap_embeddings <- Embeddings(pbmc3k, reduction = "umap")
clusters <- pbmc3k$predicted.celltype.l2

# Combine embeddings with cluster information
umap_df <- as.data.frame(umap_embeddings)
umap_df$cluster <- clusters

# Calculate the mean UMAP coordinates for each cluster
cluster_centers <- umap_df %>%
  group_by(cluster) %>%
  summarize(umap_1 = mean(umap_1), umap_2 = mean(umap_2))

# Add labels beside the clusters using ggrepel
dim_plot <- dim_plot +
  geom_text_repel(data = cluster_centers, 
                  aes(x = umap_1, y = umap_2, label = cluster),
                  size = 3,                # Adjust text size as needed
                  nudge_x = 0.5,           # Adjust horizontal position
                  nudge_y = 0.5,           # Adjust vertical position
                  box.padding = 0.35,      # Padding around labels
                  point.padding = 0.5,     # Padding around data points
                  segment.color = 'grey50',# Color of the line connecting label to cluster
                  show.legend = FALSE)     # Hide legend for labels

# Display the enhanced plot
print(dim_plot)

# Save the plot as a PNG file
ggsave(filename = "pbmc3k_UMAP_custom_labels.png", 
       plot = dim_plot, 
       width = 10, 
       height = 8, 
       dpi = 300)

########## Manual classification of clusters ##########

# This step requires manual inspection of marker genes
markers <- FindAllMarkers(pbmc3k, only.pos = TRUE, min.pct = 0.25, 
                          logfc.threshold = 0.25)
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Visualize top markers using DoHeatmap
file.name <- "pbmc3k_markers_heatmap.png"
png(file.name, width = 600, height = 900)
DoHeatmap(pbmc3k, features = top_markers$gene)
dev.off()

markers <- FindAllMarkers(pbmc3k, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top3_markers <- markers %>%
   group_by(cluster) %>%
   arrange(desc(avg_log2FC)) %>%
   slice_head(n = 3) %>%
   ungroup()
marker_genes_top3 <- top3_markers$gene

print(marker_genes_top3)

# Create a unique list of marker genes
unique_marker_genes <- unique(marker_genes_top3)

# Verify the unique genes
print(unique_marker_genes)

# Create the DotPlot
dotplot <- DotPlot(pbmc3k, 
                   features = unique_marker_genes, 
                   group.by = "predicted.celltype.l2",
                   cols = c("lightgrey", "blue")) +  # Customize colors as needed
  RotatedAxis() +                       # Rotate x-axis labels for better readability
  ggtitle("DotPlot of Top 3 Marker Genes per Cluster") +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title

# Display the plot
file.name <- "pbmc3k_markers_dotplot.png"
png(file.name, width = 900, height = 600)
print(dotplot)
dev.off()

dotplot <- DotPlot(pbmc3k, 
                   features = unique_marker_genes, 
                   group.by = "seurat_clusters",
                   cols = c("lightgrey", "blue")) +  # Customize colors as needed
  RotatedAxis() +                       # Rotate x-axis labels for better readability
  ggtitle("DotPlot of Top 3 Marker Genes per Cluster") +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title

# Display the plot
file.name <- "pbmc3k_markers_dotplot.png"
png(file.name, width = 900, height = 600)
print(dotplot)
dev.off()

# Create the FeaturePlot
gene_of_interest <- c("IGHM")
feature_plot <- FeaturePlot(pbmc3k, 
                            features = gene_of_interest, 
                            reduction = "umap", 
                            cols = c("lightgrey", "blue"),   # Customize color gradient
                            pt.size = 1.5,                   # Adjust point size
                            blend = FALSE) +                 # Disable blending (for single gene)
  ggtitle(paste("FeaturePlot of", gene_of_interest)) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))  # Center and style the title

# Display the FeaturePlot
file.name <- "pbmc3k_markers_feature-plot.png"
png(file.name, width = 900, height = 600)
print(feature_plot)
dev.off()

# Determine the number of clusters
num_clusters <- length(unique(pbmc3k$seurat_clusters))
print(paste("Number of unique clusters:", num_clusters))

# Generate a color palette with sufficient colors
# Example using pals' polychrome
palette_colors <- polychrome(n = num_clusters)
names(palette_colors) <- levels(pbmc3k$seurat_clusters)  # Assign cluster names to colors

# Define the gene of interest
gene_of_interest <- "IGHM"

# Create the ViolinPlot without auto-printing
violin_plot <- VlnPlot(pbmc3k, 
                                 features = gene_of_interest, 
                                 group.by = "seurat_clusters",  # Adjust based on your metadata
                                 pt.size = 0.1,                 # Adjust point size for individual cells
                                 cols = palette_colors) +      # Use the generated color palette
  ggtitle(paste("Violin Plot of", gene_of_interest)) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))  # Center and style the title

# Explicitly print the ViolinPlot
file.name <- "pbmc3k_markers_violin-plot.png"
png(file.name, width = 900, height = 600)
print(violin_plot)
dev.off()

# Create the RidgePlot
ridge_plot <- RidgePlot(pbmc3k, 
                        features = gene_of_interest, 
                        group.by = "seurat_clusters",
                        cols = palette_colors) +      
  ggtitle(paste("Ridge Plot of", gene_of_interest)) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

# Display the RidgePlot
file.name <- "pbmc3k_markers_ridge-plot.png"
png(file.name, width = 900, height = 600)
print(ridge_plot)
dev.off()
