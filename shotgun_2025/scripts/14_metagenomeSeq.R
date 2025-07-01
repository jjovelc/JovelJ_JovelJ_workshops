# This script reads kraken2 results and conducts
# normalization, statistical comparissons and
# plotting.

rm(list = ls())

# Make sure you have all packages.
# Otherwise:
# BiocManager::install("metagenomeSeq")
# install.packages("tidyverse")
# install.packages("reshape2")
# install.packages("RColorBrewer")
# install.packages("pheatmap")

library(metagenomeSeq)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(pheatmap)

# ===== DATA LOADING ===== #
mouse <- loadMeta("inflammation_all_samples_kraken2-counts.tsv")
meta <- loadPhenoData("metadata.txt", tran = TRUE)

# Get the counts and taxa
counts_matrix <- mouse$counts
taxa_data     <- mouse$taxa

# Ensure consistent naming
sample_names <- colnames(counts_matrix)
feature_names <- rownames(counts_matrix)

# Make sure metadata has the right row names
if (!all(sample_names %in% rownames(meta))) {
  meta <- meta[sample_names, , drop = FALSE]
}

# Make sure taxa data has the right row names
if (nrow(taxa_data) == nrow(counts_matrix)) {
  rownames(taxa_data) <- feature_names
}

# Create AnnotatedDataFrame objects
inflammData <- AnnotatedDataFrame(meta)
featureData <- AnnotatedDataFrame(taxa_data)

# Create MRexperiment object
mouseMR <- newMRexperiment(
  counts = counts_matrix,
  phenoData = inflammData,
  featureData = featureData
)

print("MRexperiment object created successfully!")
print(mouseMR)

# ===== DATA EXPLORATION ===== #
# Basic statistics
cat("\n === BASIC DATA SUMMARY ===\n")
cat("Number of features:", nrow(mouseMR), "\n")
cat("Number of samples:", ncol(mouseMR), "\n")

# Check sparsity
sparsity <- sum(MRcounts(mouseMR) == 0) / (nrow(MRcounts(mouseMR)) * ncol(MRcounts(mouseMR)))
cat("Data sparsity:", round(sparsity * 100, 2), "%\n")

# ===== FILTERING ===== #
# Filter low abundance features
# Keep features present in at least 2 samples with at least 2 counts
# filterData returns the filtered MRexperiment object directly
mouseMR_filtered <- filterData(mouseMR, present = 2, depth = 2)

cat("Original features:", nrow(mouseMR), "\n")
cat("Features after filtering:", nrow(mouseMR_filtered), "\n")

# If you want to see which features were kept, you can compare:
original_features <- rownames(MRcounts(mouseMR))
filtered_features <- rownames(MRcounts(mouseMR_filtered))

cat("Features removed:", length(original_features) - length(filtered_features), "\n")

# Check if filtering worked
if (nrow(mouseMR_filtered) > 0) {
  cat("Filtering successful!\n")
} else {
  cat("No features passed the filter, trying less stringent criteria...\n")
  # Try more lenient filtering
  mouseMR_filtered <- filterData(mouseMR, present = 2, depth = 1)
  cat("Features with present=2, depth=1:", nrow(mouseMR_filtered), "\n")
}

cat("\n=== AFTER FILTERING ===\n")
cat("Features remaining:", nrow(mouseMR_filtered), "\n")
cat("Samples remaining:", ncol(mouseMR_filtered), "\n")

# ===== NORMALIZATION =====
# Calculate normalization factors using Cumulative Sum Scaling (CSS)
mouseMR_filtered <- cumNorm(mouseMR_filtered, p = cumNormStatFast(mouseMR_filtered))

# Get normalized counts
normalized_counts <- MRcounts(mouseMR_filtered, norm = TRUE, log = TRUE)

cat("\n=== NORMALIZATION COMPLETE ===\n")
cat("Normalization factors:\n")
print(normFactors(mouseMR_filtered))

# ===== STATISTICAL ANALYSIS =====
# Set up the model for differential abundance testing
# Using group as the main factor
mod <- model.matrix(~group, data = pData(mouseMR_filtered))
settings <- zigControl()

# Perform zero-inflated Gaussian model fitting
mouseFit <- fitZig(obj = mouseMR_filtered, mod = mod, control = settings)
cat("\n=== MODEL FITTING COMPLETE ===\n")

# Extract results
zigResults <- MRcoefs(mouseFit, by = 2, coef = 2)

# Check the structure of zigResults
cat("Column names in zigResults:\n")
print(colnames(zigResults))

cat("\nFirst few rows of zigResults:\n")
print(head(zigResults))

# Add taxa information to results
if (nrow(zigResults) > 0) {
  zigResults$taxa <- taxa_data[rownames(zigResults), "taxa"]

  # Sort by adjusted p-value
  zigResults <- zigResults[order(zigResults$adjPvalues), ]

  # Filter significant results (adjust threshold as needed)
  significant_results <- zigResults[zigResults$adjPvalues < 0.05, ]

  cat("\n=== DIFFERENTIAL ABUNDANCE RESULTS ===\n")
  cat("Total features tested:", nrow(zigResults), "\n")
  cat("Significant features (padj < 0.05):", nrow(significant_results), "\n")

  if (nrow(significant_results) > 0) {
    cat("\nTop 10 significant features:\n")
    print(head(significant_results[, c("taxa", "group02_DSS", "pvalues", "adjPvalues")], 10))
  } else {
    cat("No significantly different features found at padj < 0.05\n")
    cat("Top 10 features by p-value:\n")
    print(head(zigResults[, c("taxa", "group02_DSS", "pvalues", "adjPvalues")], 10))
  }
}

# ===== VISUALIZATION ===== #


# 1. PCA plot
pca_data <- prcomp(t(normalized_counts), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca_data$x[, 1],
  PC2 = pca_data$x[, 2],
  Group = pData(mouseMR_filtered)$group,
  Sample = rownames(pca_data$x)
)

variance_explained <- round(100 * pca_data$sdev^2 / sum(pca_data$sdev^2), 1)

p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  geom_text(aes(label = Sample), vjust = -0.5, size = 3) +
  theme_minimal() +
  labs(
    title = "PCA of Normalized Counts",
    x = paste0("PC1 (", variance_explained[1], "%)"),
    y = paste0("PC2 (", variance_explained[2], "%)")
  ) +
  scale_color_brewer(type = "qual", palette = "Set1")

print(p2)

# 2. Heatmap of top variable features
# Get top 50 most variable features
var_features <- apply(normalized_counts, 1, var)
top_var_indices <- order(var_features, decreasing = TRUE)[1:min(50, length(var_features))]
top_var_data <- normalized_counts[top_var_indices, ]

# Create annotation for samples
annotation_col <- data.frame(
  Group = pData(mouseMR_filtered)$group,
  row.names = colnames(top_var_data)
)

# Create heatmap
pheatmap(top_var_data,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  annotation_col = annotation_col,
  show_rownames = FALSE,
  main = "Heatmap of Top 50 Most Variable Features",
  color = colorRampPalette(c("blue", "white", "red"))(100)
)

# 3. Volcano plot (if we have differential abundance results)
if (exists("zigResults") && nrow(zigResults) > 0) {
  volcano_data <- data.frame(
    logFC = zigResults$group02_DSS,
    neglog10p = -log10(zigResults$pvalues),
    adjPvalues = zigResults$adjPvalues,
    taxa = zigResults$taxa
  )
  
  volcano_data$significance <- "Not significant"
  volcano_data$significance[volcano_data$adjPvalues < 0.05] <- "Significant"
  volcano_data$significance[volcano_data$adjPvalues < 0.01] <- "Highly significant"
  
  # SOLUTION 1: Add jitter to separate overlapping points
  p4_jitter <- ggplot(volcano_data, aes(x = logFC, y = neglog10p, color = significance)) +
    geom_jitter(alpha = 0.8, size = 3, width = 0.05, height = 1) +  # Small jitter
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(
      title = "Volcano Plot with Jitter (Solution 1)",
      subtitle = "Slight random offset to separate overlapping points",
      x = "Log2 Fold Change (DSS vs Control)",
      y = "-log10(p-value)"
    ) +
    scale_color_manual(values = c(
      "Not significant" = "grey",
      "Significant" = "orange", 
      "Highly significant" = "red"
    ))
  
  # SOLUTION 2: Use geom_point with position_jitter
  p4_position_jitter <- ggplot(volcano_data, aes(x = logFC, y = neglog10p, color = significance)) +
    geom_point(alpha = 0.8, size = 3, position = position_jitter(width = 0.05, height = 1)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(
      title = "Volcano Plot with Position Jitter (Solution 2)",
      x = "Log2 Fold Change (DSS vs Control)",
      y = "-log10(p-value)"
    ) +
    scale_color_manual(values = c(
      "Not significant" = "grey",
      "Significant" = "orange", 
      "Highly significant" = "red"
    ))
  
  # SOLUTION 3: Use ggrepel to add labels and spread points
  library(ggrepel)
  
  # Create short labels for taxa (first few characters)
  volcano_data$short_label <- paste0("T", 1:nrow(volcano_data))
  
  p4_repel <- ggplot(volcano_data, aes(x = logFC, y = neglog10p, color = significance)) +
    geom_point(alpha = 0.8, size = 3) +
    geom_text_repel(aes(label = short_label), 
                    size = 3, 
                    box.padding = 0.5,
                    point.padding = 0.5,
                    force = 2,
                    max.overlaps = Inf) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(
      title = "Volcano Plot with ggrepel Labels (Solution 3)",
      subtitle = "Labels help identify individual points",
      x = "Log2 Fold Change (DSS vs Control)",
      y = "-log10(p-value)"
    ) +
    scale_color_manual(values = c(
      "Not significant" = "grey",
      "Significant" = "orange", 
      "Highly significant" = "red"
    ))
  
  # SOLUTION 4: Original plot with larger, semi-transparent points
  p4_large <- ggplot(volcano_data, aes(x = logFC, y = neglog10p, color = significance)) +
    geom_point(alpha = 0.6, size = 5, stroke = 1.5) +  # Larger, more transparent
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(
      title = "Volcano Plot with Large Transparent Points (Solution 4)",
      subtitle = "Overlapping points create darker areas",
      x = "Log2 Fold Change (DSS vs Control)",
      y = "-log10(p-value)"
    ) +
    scale_color_manual(values = c(
      "Not significant" = "grey",
      "Significant" = "orange", 
      "Highly significant" = "red"
    ))
  
  # Show all solutions
  print(p4_jitter)
  print(p4_position_jitter)
  
  # Only show ggrepel plot if package is available
  if(require(ggrepel, quietly = TRUE)) {
    print(p4_repel)
  } else {
    cat("Install ggrepel package for labeled solution: install.packages('ggrepel')\n")
  }
  
  print(p4_large)
  
  # Print summary
  cat("\n=== OVERLAP ANALYSIS ===\n")
  cat("Points with identical coordinates:\n")
  duplicated_coords <- duplicated(paste(volcano_data$logFC, volcano_data$neglog10p))
  if(any(duplicated_coords)) {
    cat("Found", sum(duplicated_coords), "points with identical coordinates\n")
  } else {
    cat("No exactly identical coordinates, but points may be very close\n")
  }
}
  
# 4. Box plots for top significant features (if available)
if (exists("significant_results") && nrow(significant_results) > 0) {
  # Get top 6 significant features for plotting
  top_features <- rownames(significant_results)[1:min(6, nrow(significant_results))]

  plot_data <- list()
  for (i in seq_along(top_features)) {
    feature_counts <- normalized_counts[top_features[i], ]
    plot_data[[i]] <- data.frame(
      Sample = names(feature_counts),
      Count = as.numeric(feature_counts),
      Group = pData(mouseMR_filtered)$group,
      Feature = paste0(
        "Feature_", i, ": ",
        substr(taxa_data[top_features[i], "taxa"], 1, 30)
      )
    )
  }

  combined_plot_data <- do.call(rbind, plot_data)

  p5 <- ggplot(combined_plot_data, aes(x = Group, y = Count, fill = Group)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.7) +
    facet_wrap(~Feature, scales = "free_y", ncol = 2) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 8)
    ) +
    labs(
      title = "Top Significant Features by Group",
      x = "Group", y = "Normalized Count (log)"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set1")

  print(p5)
}

# ===== SAVE RESULTS ===== #
# Save significant results to CSV
if (exists("zigResults")) {
  write.csv(zigResults, "differential_abundance_results.csv", row.names = TRUE)
  cat("\nResults saved to 'differential_abundance_results.csv'\n")
}

# Save normalized counts
write.csv(normalized_counts, "normalized_counts.csv", row.names = TRUE)
cat("Normalized counts saved to 'normalized_counts.csv'\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
