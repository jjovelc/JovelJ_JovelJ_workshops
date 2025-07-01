# This script takes as input a series of kraken2 (mpa style) generated
# by running: python parse_taxa.py kraken2_mpa-style_report.tsv
# e.g. python parse_taxa.py inflammation_all_samples_kraken2-counts.tsv
# The data will be normalized, and alpha and beta diversity analyses
# conducted, and plots will be generated

# Make sure you have the required libraries, otherwise:
# BiocManager::install("phyloseq")
# install.packages("tidyverse")
# install.packages("remotes")
# remotes::install_github("vegandevs/vegan")
rm(list = ls())
library(phyloseq)
library(tidyverse)
library(vegan)

# Change directory to the directory where you have the input data
setwd("/Users/juanjovel/Library/CloudStorage/OneDrive-UniversityofCalgary/jj/UofC/data_analysis/me/courses/2025/shotgun")

# Define the levels you want to iterate over
levels <- 1:7

# Define the taxonomic ranks
taxonomic_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Function to extract terminal leaf from taxonomic string
extract_terminal_leaf <- function(taxa_name) {
  # Split by | and get the last non-empty element
  parts <- strsplit(taxa_name, "\\|")[[1]]
  parts <- parts[parts != "" & !is.na(parts)]
  if (length(parts) > 0) {
    terminal <- parts[length(parts)]
    # Remove prefixes like k__, p__, c__, etc.
    terminal <- gsub("^[kpcofgs]__", "", terminal)
    return(terminal)
  } else {
    return("Unknown")
  }
}

# Iterate over each level
for (level in levels) {
  # Generate file name based on level
  infile <- sprintf("inflammation_all_samples_kraken2-counts_level_%d.tsv", level)
  
  # Extract level suffix and number
  level_suffix <- sprintf('level_%d', level)
  taxonomic_rank_name <- taxonomic_ranks[level]
  rank_name_for_file <- gsub(" ", "_", tolower(taxonomic_rank_name))
  
  # Read the OTU table from the file
  shotgun_table <- read.csv(infile, sep = '\t', header = TRUE, check.names = FALSE, row.names = 1)
  
  # Ensure taxa names in OTU table are consistent and formatted correctly
  taxa_names_in_otu_table <- rownames(shotgun_table)
  
  # For diversity indices calculations requiring raw counts
  otu_table_obj_raw <- otu_table(as.matrix(shotgun_table), taxa_are_rows = TRUE)
  
  # For parts of analysis where relative abundances are appropriate
  shotgun_table_relative_abundance <- apply(shotgun_table, 2, function(x) x / sum(x))
  otu_table_obj_relative <- otu_table(as.matrix(shotgun_table_relative_abundance), taxa_are_rows = TRUE)
  
  # Select the taxonomic ranks to use based on the level
  ranks_to_use <- taxonomic_ranks[1:level]
  
  # Create taxa information
  taxa_info <- data.frame(Taxa = taxa_names_in_otu_table) %>%
    separate(Taxa, into = ranks_to_use, sep = "\\|", fill = "right") %>%
    as.data.frame()
  
  row.names(taxa_info) <- taxa_names_in_otu_table
  
  tax_table_obj <- tax_table(as.matrix(taxa_info))
  
  samples <- colnames(shotgun_table)
  
  # Metadata
  sample_metadata_df <- data.frame(
    SampleID = samples,
    Condition = c("Control", "Control", "Control", "Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment", "Treatment")
  )
  
  row.names(sample_metadata_df) <- sample_metadata_df$SampleID
  
  sample_metadata <- sample_data(sample_metadata_df)
  
  # Create phyloseq object for raw counts
  physeq_raw <- phyloseq(otu_table_obj_raw, tax_table_obj, sample_metadata)
  
  # Create phyloseq object for relative abundances
  # This object can be used for beta diversity analysis or other analyses where relative abundances are appropriate
  physeq_relative <- phyloseq(otu_table_obj_relative, tax_table_obj, sample_metadata)
  
  # Generate stacked bar plots for levels 2-7 (Phylum through Species)
  if (level >= 2) {
    cat(sprintf("\nGenerating stacked bar plot for %s (Level %d)...\n", taxonomic_rank_name, level))
    
    # Calculate mean relative abundance across all samples for each taxon
    mean_abundance <- rowMeans(shotgun_table_relative_abundance)
    
    # Get top 20 most abundant taxa (or max available if less than 20)
    n_taxa <- min(20, length(mean_abundance))
    top_taxa_indices <- order(mean_abundance, decreasing = TRUE)[1:n_taxa]
    top_taxa_names <- names(mean_abundance)[top_taxa_indices]
    
    # Subset data to top taxa
    top_taxa_data <- shotgun_table_relative_abundance[top_taxa_names, , drop = FALSE]
    
    # Extract terminal leaves for legend
    terminal_leaves <- sapply(rownames(top_taxa_data), extract_terminal_leaf)
    
    # Create data frame for plotting
    plot_data <- top_taxa_data %>%
      as.data.frame() %>%
      rownames_to_column("Taxa") %>%
      mutate(Terminal_Leaf = terminal_leaves) %>%
      pivot_longer(cols = -c(Taxa, Terminal_Leaf), names_to = "Sample", values_to = "Relative_Abundance") %>%
      left_join(sample_metadata_df, by = c("Sample" = "SampleID"))
    
    # Create discrete color palette with support for >20 colors
    n_colors <- length(unique(plot_data$Terminal_Leaf))
    
    # Define a comprehensive palette combining multiple discrete palettes
    discrete_palette <- c(
      # RColorBrewer Set3 (12 colors)
      RColorBrewer::brewer.pal(12, "Set3"),
      # RColorBrewer Paired (12 colors)
      RColorBrewer::brewer.pal(12, "Paired"),
      # RColorBrewer Set1 (9 colors)
      RColorBrewer::brewer.pal(9, "Set1"),
      # RColorBrewer Set2 (8 colors)  
      RColorBrewer::brewer.pal(8, "Set2"),
      # RColorBrewer Dark2 (8 colors)
      RColorBrewer::brewer.pal(8, "Dark2"),
      # RColorBrewer Accent (8 colors)
      RColorBrewer::brewer.pal(8, "Accent"),
      # Additional pastel colors
      "#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", "#BAE1FF", "#E6BAFF",
      "#FFCCCB", "#FFDAB9", "#E0E0E0", "#B0E0E6", "#F0E68C", "#DDA0DD",
      "#98FB98", "#F5DEB3", "#D3D3D3", "#AFEEEE", "#FAFAD2", "#THISTLE"
    )
    
    # Remove duplicates and select needed colors
    discrete_palette <- unique(discrete_palette)
    colors <- discrete_palette[1:min(n_colors, length(discrete_palette))]
    
    # If we still need more colors, supplement with generated pastels
    if (n_colors > length(discrete_palette)) {
      additional_needed <- n_colors - length(discrete_palette)
      hues <- seq(0, 1, length.out = additional_needed + 1)[1:additional_needed]
      additional_colors <- hsv(h = hues, s = 0.35, v = 0.85)
      colors <- c(colors, additional_colors)
    }
    
    # Create stacked bar plot
    p_stacked <- ggplot(plot_data, aes(x = Sample, y = Relative_Abundance, fill = Terminal_Leaf)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = colors) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      ) +
      labs(
        x = "Sample",
        y = "Relative Abundance",
        fill = taxonomic_rank_name,
        title = paste("Top", n_taxa, "Most Abundant", taxonomic_rank_name, "- Stacked Bar Plot")
      ) +
      facet_wrap(~ Condition, scales = "free_x", nrow = 1)
    
    # Save stacked bar plot
    stacked_plot_filename <- paste("stacked_barplot_", rank_name_for_file, "_top", n_taxa, ".png", sep = "")
    ggsave(stacked_plot_filename, plot = p_stacked, width = 12, height = 8, dpi = 300)
    
    cat(sprintf("Stacked bar plot saved as: %s\n", stacked_plot_filename))
  }
  
  # Alpha diversity calculations should use raw counts
  indices <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
  for (index in indices) {
    p <- plot_richness(physeq_raw, x = "Condition", measures = index) + 
      geom_boxplot(aes(color = Condition)) +
      ggtitle(paste(level_suffix, index, sep = "_"))
    
    # Save plot
    plot_filename <- paste("alpha_diversity_", rank_name_for_file, "_", index, ".png", sep = "")
    ggsave(plot_filename, plot = p, width = 10, height = 8)
  }
  
  # Beta diversity analysis can use the relative abundance data
  distances <- c("bray", "jaccard", "jsd")
  for (dist in distances) {
    distance_matrix <- phyloseq::distance(physeq_relative, method = dist)
    pcoa_results <- cmdscale(distance_matrix, eig = TRUE)
    pcoa_df <- as.data.frame(pcoa_results$points)
    pcoa_df$condition <- sample_metadata_df$Condition
    
    # Conduct PERMANOVA using adonis2
    permanova_results <- adonis2(distance_matrix ~ Condition, data = sample_metadata_df)
    # Conduct ANOSIM
    anosim_results <- anosim(distance_matrix, sample_metadata_df$Condition)
    
    # Extract p-value from PERMANOVA results
    permanova_pvalue <- permanova_results$`Pr(>F)`[1]
    
    # Add error checking for p-value extraction
    if (is.null(permanova_pvalue) || is.na(permanova_pvalue)) {
      permanova_pvalue <- "NA"
      cat("Warning: Could not extract PERMANOVA p-value\n")
    }
    
    # Print results
    cat(sprintf("\n%s %s PERMANOVA Results:\n", level_suffix, dist))
    print(permanova_results)
    cat(sprintf("\n%s %s ANOSIM Results:\n", level_suffix, dist))
    print(anosim_results)
    
    # Plot PCoA with corrected p-value access
    p <- ggplot(pcoa_df, aes(x = V1, y = V2, color = condition)) +
      geom_point(size = 5) +
      geom_text(aes(label = row.names(pcoa_df)), nudge_x = 0.02, nudge_y = 0.02, check_overlap = TRUE) +
      theme_bw() +
      theme(legend.position = 'top', legend.title = element_blank()) +
      ggtitle(paste(level_suffix, dist, "PERMANOVA p =", format(permanova_pvalue, digits = 4), "ANOSIM R =", format(anosim_results$statistic, digits = 4)))
    
    # Save plot with results
    plot_filename <- paste("pcoa_", rank_name_for_file, "_", dist, "_with_stats.png", sep = "")
    ggsave(plot_filename, plot = p, width = 10, height = 8)
  }
}
