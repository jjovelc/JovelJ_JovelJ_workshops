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
library(phyloseq)
library(tidyverse)
library(vegan)

# Change directory to the directory where you have the input data
setwd("/Users/juanjovel/Library/CloudStorage/OneDrive-UniversityofCalgary/jj/UofC/data_analysis/me/courses/2025/shotgun")

# Define the levels you want to iterate over
levels <- 1:7

# Define the taxonomic ranks
taxonomic_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

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