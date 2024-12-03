###########################################################
# This script is intended to recapitulate the working     #
# working session in course: "scRNAseq data analysis"     #
# which was held on Nov/ 04-05, 2024 by Juan Jovel        #
#                                                         #
# Use it at your own risk                                 #
###########################################################

# load required libraries
library(Seurat)          # main library for the analysis of scRNAseq data
library(EnhancedVolcano) # Produces high quality, customizable, volcano plots

# Set working directory
setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/me/courses/2024/scRNAseq/scRNAseq_course2024_materials')

#  Import covid RDS file
covid <- readRDS('real_covid_10Ksample.rds')

# Calculate nFeature_RNA and nCount_RNA if they don't exist
covid <- AddMetaData(covid, metadata = Matrix::colSums(covid@assays$RNA$counts > 0), col.name = "nFeature_RNA")
covid <- AddMetaData(covid, metadata = Matrix::colSums(covid@assays$RNA$counts), col.name = "nCount_RNA")

# Calculate the percentage of mitochondrial reads
covid[["percent.mt"]] <- PercentageFeatureSet(covid, pattern = "^MT-")



# Filter cells based on QC criteria
covid <- subset(covid, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10 &
                nCount_RNA > 500 & nCount_RNA < 10000)

# Ensure 'disease' is set as the active identity class
Idents(covid) <- "disease"

table(covid@meta.data$disease)

# Perform differential expression analysis between COVID-19 and normal
de_results <- FindMarkers(
  object = covid,
  ident.1 = "COVID-19",
  ident.2 = "normal",
  assay = "RNA", # Specify the RNA assay layer
  slot = "data", # Use normalized data (change to "counts" if raw counts are desired)
  test.use = "wilcox" # Choose the test (e.g., Wilcoxon by default)
)

# View the top differential expression results
head(de_results)

sign_res <- subset(de_results, abs(avg_log2FC) > 0.25 &  p_val_adj < 0.05)
sign_res$gene <- rownames(sign_res)
sign_res <- sign_res[, c('gene', colnames(sign_res)[1:(ncol(sign_res) - 1)])]

# Export results
write.table(sign_res, "DE_results_wilcox.tsv", sep = '\t', quote = F)

library(EnhancedVolcano)

# Generate the volcano plot
EnhancedVolcano(
  de_results,
  lab       = rownames(de_results),
  x         = 'avg_log2FC', # Column for log fold change
  y         = 'p_val_adj',  # Column for adjusted p-values
  xlab      = 'Log2 Fold Change',
  ylab      = '-Log10 Adjusted P-value',
  title     = 'Differential Expression COVID-19 vs. Normal',
  pCutoff   = 0.05,         # Adjust this threshold as needed
  FCcutoff  = 0.25,         # Adjust this fold-change threshold as needed
  pointSize = 3.0,
  labSize   = 3.0
)

########## What about cells? ##########
# Get list of cell types with counts >= 20
cell_types <- names(which(table(covid@meta.data$predicted.celltype.l2) >= 20))

# Initialize list to store DE results
de_results_list <- list()

for (cell in cell_types) {
  # Subset cells of the current cell type
  subset_cells <- subset(covid, subset = predicted.celltype.l2 == cell)
  
  # Ensure there are at least two disease groups (COVID-19 and normal) to proceed with DE analysis
  if (length(unique(subset_cells@meta.data$disease)) < 2) {
    message(paste("Skipping", cell, "- not enough groups"))
    next
  }
  
  # Perform differential expression analysis
  de_markers <- FindMarkers(
    object          = subset_cells,
    ident.1         = "COVID-19",
    ident.2         = "normal",
    logfc.threshold = 0.25,
    min.pct         = 0.1,
    test.use        = "wilcox",
    latent.vars     = "donor_id"
  )
  
  # Store results
  de_markers$cell_type    <- cell
  de_results_list[[cell]] <- de_markers
}

# Combine results
de_results_cell_types <- bind_rows(de_results_list, .id = "cell_type")
sign_res <- subset(de_results_cell_types, abs(avg_log2FC) > 0.25 &  p_val_adj < 0.05)
sign_res$gene <- rownames(sign_res)

# Move the new "gene" column to the beginning of the dataframe
sign_res <- sign_res[, c("gene", colnames(sign_res)[1:(ncol(sign_res) - 1)])]

write.table(sign_res, "DE_genes_normal_vs_COVID19_cell-type.tsv", 
            sep='\t', row.names = F, quote = F)
