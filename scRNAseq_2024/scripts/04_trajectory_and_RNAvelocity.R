###########################################################
# This script is intended to recapitulate the working     #
# working session in course: "scRNAseq data analysis"     #
# which was held on Nov/ 04-05, 2024 by Juan Jovel        #
#                                                         #
# Use it at your own risk                                 #
###########################################################

library(Seurat)  # For Monocle3 and scVelo compatibility
library(monocle3)
library(Matrix)
library(SeuratWrappers)
library(magrittr)
library(patchwork)

setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/me/courses/2024/scRNAseq/initial_quants/alevin_output_3p_ACDA')

pbmc4k <- readRDS('pbmc4k.rds')

# Monocle requires that NFSR is run before implamentation

# Convert the Seurat object to a Monocle3 CellDataSet.
cds <- as.cell_data_set(pbmc4k)

# Preprocess with 20 principal components
cds <- preprocess_cds(cds, num_dim = 20)

# Perform UMAP dimensionality reduction to visualize the data.
cds <- reduce_dimension(cds, reduction_method = "UMAP")

# Cluster the cells to identify distinct cell groups.
cds <- cluster_cells(cds, resolution = 1e-3)

# Construct the trajectory graph to enable pseudotime analysis.
cds <- learn_graph(cds)

# Identify a root node, or let Monocle3 pick the best one
cds <- order_cells(cds, reduction_method = "UMAP")

# Plot Trajectories and Pseudotime
plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = TRUE)

# Plot partitions
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = TRUE)

# Plot pseudo-time
plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)

##############


