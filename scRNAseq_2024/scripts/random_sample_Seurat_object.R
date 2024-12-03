load your libraries and Seurat object
...

# First approach
# Define target sample size
sample_size <- 100000

# Subset and sample each class
set.seed(123) # for reproducibility
seurat_object.g_equal <- subset(seurat_object, idents = c("COVID-19", "normal"))
seurat_object.g_equal <- subset(seurat_object_equal, cells = unlist(lapply(levels(Idents(seurat_object_equal)), 
                                                         function(x) sample(Cells(subset(seurat_object.g_equal, idents = x)), sample_size))))
seurat_object.g_equal <- RunAzimuth(seurat_object_equal, reference = "pbmcref")

# Verify the cell numbers
table(Idents(seurat_object.g_equal))

bulk <- AggregateExpression(seurat_object.g_equal,
                            return.seurat = TRUE,
                            assays = "RNA",
                            group.by = c("predicted.celltype.l2", 
                                         "donor_id", "disease"))

bulk <- subset(bulk, subset = disease %in% c("normal", "COVID-19"))
bulk <- subset(bulk, subset = predicted.celltype.l2 != "Doublet")
bulk$disease <- factor(bulk$disease, levels = c("normal", "COVID-19"))

# Verify the cell numbers
table(Idents(bulk))

saveRDS(bulk, "real_seurat_object.g_10Ksample.rds")



