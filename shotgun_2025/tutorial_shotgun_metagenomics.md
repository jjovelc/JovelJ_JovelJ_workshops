# Shotgun Metagenomics Analysis Pipeline Tutorial

## Overview

This tutorial provides a complete, step-by-step guide for analyzing shotgun metagenomics data. The pipeline takes you from raw sequencing reads through quality control, taxonomic classification, functional profiling, and statistical analysis. All commands can be copied and pasted directly into your terminal. When appropriate, clear instructions are given to execute Python or R scritps.

### What You'll Learn
- Quality assessment and control of metagenomic reads
- Host DNA removal from microbiome samples
- Taxonomic classification and abundance estimation
- Functional profiling of microbial communities
- Diversity analysis and statistical comparisons
- Data visualization and interpretation

### Prerequisites
- Basic command line knowledge
- Access to a high-performance computing cluster with SLURM
- Conda/Mamba package manager installed

## Pipeline Overview

The complete pipeline consists of 14 steps:

1. **Quality Control** - Initial assessment with FastQC/MultiQC
2. **Adapter Removal** - Remove sequencing adapters with Cutadapt
3. **Low Quality Filtering** - Filter poor quality sequences with Prinseq
4. **Host Removal** - Remove host DNA contamination with BMTagger
5. **Post-filtering QC** - Quality assessment after filtering
6. **Taxonomic Classification** - Classify reads with Kraken2/Bracken
7. **Functional Profiling** - Metabolic pathway analysis with HUMAnN
8. **HUMAnN Post-processing** - Normalize and regroup functional data
9. **Data Merging** - Combine sample results
10. **Data Filtering** - Remove low-abundance features
11. **Taxa Parsing** - Extract taxonomic levels
12. **Diversity Analysis** - Alpha and beta diversity with R
13. **Differential Abundance Analysis (DESeq2)** - Statistical testing with R
14. **Statistical Analysis (metagenomeSeq)** - Alternative statistical approach

---

## Environment Setup

### 1. Create Conda Environment

First, create the main analysis environment:

```bash
# Create environment from YAML file
mamba env create -f shotgun2025.yaml

# Activate environment
conda activate shotgun2025
```
NOTE: It is possible that, depending on the configuration of your user in ARC, some software might fail to install during the previous step. If that happens, you can create an additional mamba/conda environment for that app. Let's assume that cutadapt failed to install. Do the following:

```bash
mamba create -n cutadapt cutadapt
```

Make sure you activate the environment, ideally right before you call the software.

### 2. Directory Structure

Set up your working directory:

```bash
# Create main project directory
mkdir -p shotgun2025
cd shotgun2025

# Create subdirectories
mkdir -p {data,scripts}

cp /work/vetmed_shared_dbs/shotgun_workshop_2025/data/*fq ./data

cp /work/vetmed_shared_dbs/shotgun_workshop_2025/scripts/* ./scripts/

```

---

## Step 1: Initial Quality Control

### Purpose
Assess the quality of raw sequencing reads to identify potential issues before processing.

###### Script: `01_raw_fastqc_multiqc.sh`

### How to Run

```bash

# Submit to SLURM
sbatch scripts/01_raw_fastqc_multiqc.sh

# Check job status
squeue -u $USER

# View results
ls output/01_raw_fastqc/
ls output/01_raw_multiqc/
```

### Expected Output

- Individual FastQC reports for each sample
- Aggregated MultiQC report showing quality metrics across all samples
- HTML report viewable in web browser

Figure 1. Typical fastqc and multiqc results files.

![fastqc and multiqc results](images/01_fastqc_multiqc_results.png)

If you don't get the expected results, you can inspect the corresponding log files.

Figure 2. Log files corresponding to the run of fastqc and multiqc.

![fastqc and multiqc log files](images/01_logFiles_fastqc_multiqc.png)

---

## Step 2: Adapter Removal

### Purpose
Remove sequencing adapters and perform quality trimming to improve data quality.

###### Script: `02_cutadapt.sh`

### How to Run

```bash
# Submit job
sbatch scripts/02_cutadapt.sh

# Monitor progress
tail -f logs/cutadapt_sbatch_job.*.out
```

---

## Step 3: Low Complexity Filtering

### Purpose
Remove low-complexity sequences and apply additional quality filters using Prinseq.

###### Script: `03_prinseq.sh`

### How to Run

```bash
# Make sure prinseq-lite.pl is in the scripts directory

# Submit job
sbatch scripts/03_prinseq.sh
```
---

## Step 4: Host DNA Removal

### Purpose
Remove host contamination (mouse reads) from microbiome samples using BMTagger.

###### Script: `04_bmtagger.sh`

### How to Run

```bash
sbatch scripts/04_bmtagger.sh
```

### Database Requirements
BMTagger requires pre-built host reference databases. For mouse samples, you need:
- `mice_reference.bitmask` - BMFilter index
- `mice_reference.srprism` - SRPrism index

---

## Step 5: Post-Host-Removal Quality Control

### Purpose
Assess data quality after host removal to ensure successful filtering.

### Script: `05_bmtagger_fastqc_multiqc.sh`

```bash
#!/usr/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=7-00:00:00
#SBATCH --mem=4G
#SBATCH --output=logs/QC_bmtagged_sbatch_job.%A.out
#SBATCH --error=logs/QC_bmtagged_sbatch_job.%A.err

# Initialize conda
conda init bash &> /dev/null
source ~/.bashrc &> /dev/null

echo "Starting post-BMTagger quality control..."
echo "Started at: `date`"

conda activate shotgun2025

# Input directory
bmtagger_out_dir="$(pwd)/output/04_bmtagger/"
forward_reads="_bmtagged_1.fastq"
reverse_reads="_bmtagged_2.fastq"
bmtagger_fastqc_dir="$(pwd)/output/05_bmtagger_fastqc/"

mkdir -p ${bmtagger_fastqc_dir}

# Run FastQC on BMTagger output
for sample in "${bmtagger_out_dir}"/*"$forward_reads"; do
    SAMPLE=$(basename "${sample}")
    SAMPLE=$(echo "$SAMPLE" | cut -d'_' -f1)
    echo "Processing sample: ${SAMPLE}"

    echo "${bmtagger_out_dir}${SAMPLE}${forward_reads}"
    echo "${bmtagger_out_dir}${SAMPLE}${reverse_reads}"

    # Run FastQC
    fastqc -o ${bmtagger_fastqc_dir} -t 20 \
        ${bmtagger_out_dir}${SAMPLE}${forward_reads} \
        ${bmtagger_out_dir}${SAMPLE}${reverse_reads}
done

echo "Running MultiQC..."

bmtagger_multiqc_dir="$(pwd)/output/05_bmtagger_multiqc"
mkdir -p ${bmtagger_multiqc_dir}

conda activate multiqc

# Run MultiQC
multiqc -f ${bmtagger_fastqc_dir} -o ${bmtagger_multiqc_dir} -n multiqc_report.html

echo "All samples processed."
echo "Finished at: `date`"
```

### How to Run

```bash
sbatch scripts/05_bmtagger_fastqc_multiqc.sh

# View results
firefox output/05_bmtagger_multiqc/multiqc_report.html
```

---

## Step 6: Taxonomic Classification

### Purpose
Classify reads taxonomically and estimate abundance using Kraken2 and Bracken.

### Script: `06_kraken_bracken_5.sh`

```bash
#!/usr/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=7-00:00:00
#SBATCH --mem=128G
#SBATCH --output=logs/annotation_sbatch_job.%A.out
#SBATCH --error=logs/annotation_sbatch_job.%A.err

echo "Starting taxonomic classification..."
echo "Started at: `date`"

# Initialize conda
eval "$(conda shell.bash hook)"
conda activate shotgun2025

# Directory paths
bmtagger_out_dir="$(pwd)/output/04_bmtagger/"
Kraken_DIR="$(pwd)/output/06_kraken2"
Bracken_DIR="$(pwd)/output/06_bracken2"

# Create output directories
mkdir -p "${Kraken_DIR}" "${Bracken_DIR}"

# Path to Kraken2 database
DB="/work/vetmed_shared_dbs/kraken2_dbs/kraken2_NCBI_Oct22/"

echo "First Kraken step to produce MPA-style reports..."

# First Kraken2 run - MPA style reports
for sample in "${bmtagger_out_dir}"/*_bmtagged_1.fastq; do
    SAMPLE=$(basename "${sample}" "_bmtagged_1.fastq")
    echo "${SAMPLE}: First alignment"

    R1=${SAMPLE}_bmtagged_1.fastq
    R2=${SAMPLE}_bmtagged_2.fastq

    # Run Kraken2 with MPA style output
    kraken2 \
        --db "${DB}" \
        --threads 48 \
        --memory-mapping \
        --confidence 0.1 \
        --report "${Kraken_DIR}/${SAMPLE}_mpa.tax" \
        --paired "${bmtagger_out_dir}/${R1}" "${bmtagger_out_dir}/${R2}" \
        --output "${Kraken_DIR}/${SAMPLE}.krk2"
done

echo "Second Kraken step for Bracken compatibility..."

# Second Kraken2 run - standard reports for Bracken
for sample in "${bmtagger_out_dir}"/*_bmtagged_1.fastq; do
    SAMPLE=$(basename "${sample}" "_bmtagged_1.fastq")

    R1=${SAMPLE}_bmtagged_1.fastq
    R2=${SAMPLE}_bmtagged_2.fastq

    echo "[INFO] Processing sample: ${SAMPLE}"

    kraken2 \
        --db "${DB}" \
        --threads 48 \
        --memory-mapping \
        --confidence 0.1 \
        --paired "${bmtagger_out_dir}/${R1}" "${bmtagger_out_dir}/${R2}" \
        --report "${Kraken_DIR}/${SAMPLE}.krk2rpt" \
        --output "${Kraken_DIR}/${SAMPLE}.krk2"

    echo "[INFO] Finished sample: ${SAMPLE}"
done

echo "Running Bracken for abundance estimation..."

# Run Bracken on Kraken reports
for sample in "${Kraken_DIR}"/*krk2rpt; do
    SAMPLE=$(basename "${sample}")
    echo "Processing ${SAMPLE}"

    bracken \
        -d $DB \
        -i "${Kraken_DIR}/$SAMPLE" \
        -r 300 \
        -t 100 \
        -o "${Bracken_DIR}"/${SAMPLE/.krk2rpt/.bck} \
        -w "${Bracken_DIR}"/${SAMPLE/.krk2rpt/.bckrpt}
done

echo "Generating visualization files..."

Kraken_bracken_figures="$(pwd)/output/06_Kr-Br-html_reports"
mkdir -p "${Kraken_bracken_figures}"

# Convert to Krona format
for sample in "${Bracken_DIR}"/*bckrpt; do 
    kreport2krona.py -r $sample -o ${sample/.bckrpt/_krona.txt} --no-intermediate-ranks
done

# Generate individual Krona plots
for sample in "${Bracken_DIR}"/*krona.txt; do 
    ktImportText $sample -o "${Kraken_bracken_figures}"/${sample/txt/html}
done

# Generate combined Krona plot
ktImportText "${Bracken_DIR}"/*krona.txt -o "${Kraken_bracken_figures}"/all_krona_plots.html

echo "All samples processed."
echo "Finished at: `date`"
```

### How to Run

```bash
sbatch scripts/06_kraken_bracken_5.sh

# View Krona plots
firefox output/06_Kr-Br-html_reports/all_krona_plots.html
```

### Database Requirements
Requires pre-built Kraken2 database. The script assumes access to a shared database, but you can build your own:

```bash
# Download and build standard database (requires ~100GB disk space)
kraken2-build --standard --threads 32 --db kraken2_standard_db
```

---

## Step 7: Functional Profiling

### Purpose
Analyze metabolic pathways and gene families using HUMAnN.

### Script: `07_humann_multiple_samples.sh`

```bash
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=0-10:00:00
#SBATCH --mem=32G
#SBATCH --error=logs/humann_run.%J.err
#SBATCH --output=logs/humann_run.%J.out
#SBATCH --array=0-9 # job array index

# Initialize conda
eval "$(conda shell.bash hook)"
conda activate shotgun2025

bmtagger_out_dir="$(pwd)/output/04_bmtagger/"

# Path to HUMAnN databases
DB="/work/vetmed_shared_dbs/humann3_dbs"
nuc_db="${DB}/chocophlan"
prot_db="${DB}/uniref"

forward_reads="_bmtagged_1.fastq"
reverse_reads="_bmtagged_2.fastq"

# Get sample names array
names=($(for FILE in *_R1_100K.fq; do echo "${FILE/_R1_100K.fq/}"; done))

SAMPLE=${names[${SLURM_ARRAY_TASK_ID}]}
echo "Processing sample: ${SAMPLE}"

merged_dir="$(pwd)/output/07_merged_reads/output_${SAMPLE}"
mkdir -p "${merged_dir}"

echo "Starting to merge reads"

# Construct file paths
input1="${bmtagger_out_dir}/${SAMPLE}${forward_reads}"
input2="${bmtagger_out_dir}/${SAMPLE}${reverse_reads}"
output="${merged_dir}/${SAMPLE}_merged.fastq"

# Merge paired-end reads
cat ${input1} ${input2} > ${output}

echo "Reads merged successfully"
echo "Starting HUMAnN analysis"

# HUMAnN output directory
humann_out_dir="$(pwd)/output/07_humann/"
mkdir -p "${humann_out_dir}"

echo "Started at: `date`"

# Run HUMAnN
humann --input ${merged_dir}/${SAMPLE}_merged.fastq \
       --output ${humann_out_dir}/${SAMPLE}_humann3_res \
       --threads 6 \
       --protein-database ${prot_db} \
       --nucleotide-database ${nuc_db} \
       --metaphlan-options=" --index mpa_vJun23_CHOCOPhlAnSGB_202307 --bowtie2db /bulk/IMCbinf_bulk/sbagheri/Projects_IMC/databases/metaphlan/vJun23_version"

echo "Finished at: `date`"
```

### How to Run

```bash
# Submit array job (processes multiple samples in parallel)
sbatch scripts/07_humann_multiple_samples.sh

# Check array job status
squeue -u $USER
```

### Database Requirements
HUMAnN requires several databases:
- ChocoPhlAn nucleotide database
- UniRef protein database
- MetaPhlAn marker gene database

Download instructions:
```bash
# Download HUMAnN databases (requires significant disk space)
humann_databases --download chocophlan full /path/to/databases/
humann_databases --download uniref uniref90_diamond /path/to/databases/
```

---

## Step 8: HUMAnN Post-processing

### Purpose
Normalize and regroup HUMAnN output for downstream analysis.

### Script: `08_postprocessing_humann3.sh`

```bash
#!/usr/bin/bash

echo "Starting HUMAnN post-processing..."

# Navigate to HUMAnN output directory
cd output/07_humann/

echo "Normalizing gene families to CPM..."
# Normalize gene families to copies per million (CPM)
for FILE in *genefamilies.tsv; do
    echo "Processing $FILE"
    humann_renorm_table --input $FILE \
        --output ${FILE/.tsv/_cpm.tsv} --units cpm --update-snames
done

echo "Regrouping gene families to reactions..."

# Regroup normalized data to reactions
echo "Processing normalized data"
for FILE in *_cpm.tsv; do
    echo "Regrouping $FILE"
    humann_regroup_table --input $FILE \
        --output ${FILE/.tsv/_rxn.tsv} --groups uniref90_rxn
done

# Regroup raw data to reactions
echo "Processing raw data"
for FILE in *_genefamilies.tsv; do
    echo "Regrouping $FILE"
    humann_regroup_table --input $FILE \
        --output ${FILE/.tsv/_rxn.tsv} --groups uniref90_rxn
done

echo "Renaming reactions with MetaCyc names..."

# Rename normalized reactions
echo "Processing normalized reaction data"
for FILE in *_cpm_rxn.tsv; do
    echo "Renaming $FILE"
    humann_rename_table --input $FILE \
        --output ${FILE/.tsv/_ren.tsv} --names metacyc-rxn
done

# Rename raw reactions
echo "Processing raw reaction data"
for FILE in *_genefamilies_rxn.tsv; do
    echo "Renaming $FILE"
    humann_rename_table --input $FILE \
        --output ${FILE/.tsv/_ren.tsv} --names metacyc-rxn
done

echo "HUMAnN post-processing complete!"
```

### How to Run

```bash
# Make script executable and run
chmod +x scripts/08_postprocessing_humann3.sh
./scripts/08_postprocessing_humann3.sh
```

### Output Files Explanation
- `*_genefamilies.tsv`: Raw gene family abundances
- `*_genefamilies_cpm.tsv`: Normalized gene families (CPM)
- `*_genefamilies_rxn.tsv`: Gene families regrouped to reactions
- `*_genefamilies_cpm_rxn_ren.tsv`: Final normalized, regrouped, and renamed files

---

## Step 9: Data Merging

### Purpose
Combine individual sample results into unified tables for analysis.

### Script: `09_merge.py`

The merge script is provided and should be used as follows:

```bash
# Navigate to Kraken results directory
cd output/06_kraken2/

# Merge MPA-style taxonomic profiles
python ../../scripts/09_merge.py *_mpa.tax > ../../merged_taxonomy_table.tsv

# Navigate to HUMAnN results
cd ../07_humann/

# Merge gene family tables (example)
python ../../scripts/09_merge.py *_genefamilies_cpm.tsv > ../../merged_genefamilies_cpm.tsv

# Merge pathway abundance tables
python ../../scripts/09_merge.py *_pathabundance_cpm.tsv > ../../merged_pathways_cpm.tsv

# Return to main directory
cd ../../
```

### How to Run

```bash
# Make sure you're in the correct directory with results
# Run merging for different data types

# For taxonomic data
python scripts/09_merge.py output/06_kraken2/*_mpa.tax > merged_taxonomy_table.tsv

# For gene families
python scripts/09_merge.py output/07_humann/*_genefamilies_cpm.tsv > merged_genefamilies_cpm.tsv

# For pathways
python scripts/09_merge.py output/07_humann/*_pathabundance_cpm.tsv > merged_pathways_cpm.tsv
```

---

## Step 10: Data Filtering

### Purpose
Remove low-abundance features to focus on meaningful biological signals.

### Script: `10_filter_by_average.sh`

```bash
#!/bin/bash

# Usage: ./filter_by_average.sh input_file threshold

infile=$1
threshold=$2

# Validate input
if [[ -z "$infile" || -z "$threshold" ]]; then
    echo "Usage: $0 <input_file> <threshold>"
    echo "Example: $0 merged_taxonomy_table.tsv 0.01"
    exit 1
fi

# Filter features by average abundance across samples
awk -v threshold="$threshold" '{
    for (i = 2; i <= NF; i++) {
        sum += $i
    } 
    if ((sum / (NF - 1) >= threshold) || NF == 1) {
        print
    } 
    sum = 0
}' < "$infile"
```

### How to Run

```bash
# Filter taxonomy table (keep features with >0.01% average abundance)
./scripts/10_filter_by_average.sh merged_taxonomy_table.tsv 0.01 > filtered_taxonomy_table.tsv

# Filter gene families (keep features with >10 average CPM)
./scripts/10_filter_by_average.sh merged_genefamilies_cpm.tsv 10 > filtered_genefamilies_cpm.tsv

# Filter pathways (keep features with >1 average CPM)
./scripts/10_filter_by_average.sh merged_pathways_cpm.tsv 1 > filtered_pathways_cpm.tsv
```

---

## Step 11: Taxa Parsing

### Purpose
Extract specific taxonomic levels from the merged taxonomy table for level-specific analysis.

### Script: `11_parse_taxa.py`

```bash
# Usage with the inflammation dataset
python scripts/11_parse_taxa.py inflammation_all_samples_kraken2-counts.tsv
```

This script will generate separate files for each taxonomic level:
- `inflammation_all_samples_kraken2-counts_level_1.tsv` (Kingdom)
- `inflammation_all_samples_kraken2-counts_level_2.tsv` (Phylum)
- `inflammation_all_samples_kraken2-counts_level_3.tsv` (Class)
- `inflammation_all_samples_kraken2-counts_level_4.tsv` (Order)
- `inflammation_all_samples_kraken2-counts_level_5.tsv` (Family)
- `inflammation_all_samples_kraken2-counts_level_6.tsv` (Genus)
- `inflammation_all_samples_kraken2-counts_level_7.tsv` (Species)

### How to Run

```bash
# Parse your filtered taxonomy table
python scripts/11_parse_taxa.py filtered_taxonomy_table.tsv

# Or use the provided inflammation dataset
python scripts/11_parse_taxa.py inflammation_all_samples_kraken2-counts.tsv
```

---

## Step 12: Alpha and Beta Diversity Analysis

### Purpose
Analyze microbial diversity within samples (alpha) and between samples (beta), with statistical testing and visualization.

### Script: `12_alpha-beta_diversity_norm.R`

This comprehensive R script performs:
- Alpha diversity calculations (Shannon, Simpson, Chao1, etc.)
- Beta diversity analysis (PCoA with multiple distance metrics)
- Statistical testing (PERMANOVA, ANOSIM)
- Stacked bar plots of most abundant taxa
- Analysis across all taxonomic levels

### Required R packages

```r
# Install required packages if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")
install.packages(c("tidyverse", "vegan", "RColorBrewer"))

# If vegan installation fails, try:
# remotes::install_github("vegandevs/vegan")
```

### How to Run

```bash
# First, ensure you have the required input file
# The script expects: inflammation_all_samples_kraken2-counts.tsv
# This should be your merged and filtered taxonomy count table

# Create metadata file
cat > metadata.txt << 'EOF'
SampleID	group	label
CA1E	Control	CA1E
CA2E	Control	CA2E
CA3E	Control	CA3E
CA4E	Control	CA4E
CA5E	Control	CA5E
TA1E	Treatment	TA1E
TA2E	Treatment	TA2E
TA3E	Treatment	TA3E
TA4E	Treatment	TA4E
TA5E	Treatment	TA5E
EOF

# Run the R script
Rscript scripts/12_alpha-beta_diversity_norm.R
```

### Key Features of the Script

1. **Taxonomic Level Processing**: Analyzes all 7 taxonomic levels automatically
2. **Alpha Diversity Metrics**:
   - Observed species richness
   - Chao1 estimator
   - ACE estimator
   - Shannon diversity
   - Simpson diversity
   - Inverse Simpson

3. **Beta Diversity Analysis**:
   - Bray-Curtis dissimilarity
   - Jaccard distance
   - Jensen-Shannon divergence
   - PCoA visualization with statistical testing

4. **Statistical Testing**:
   - PERMANOVA for testing group differences
   - ANOSIM for alternative group testing
   - Automated p-value reporting

5. **Visualization**:
   - Box plots for alpha diversity metrics
   - PCoA plots with group coloring and statistics
   - Stacked bar plots showing top 20 most abundant taxa

### Expected Output Files

The script generates multiple files for each taxonomic level and analysis:

**Alpha Diversity Plots**:
- `alpha_diversity_kingdom_Observed.png`
- `alpha_diversity_kingdom_Shannon.png`
- `alpha_diversity_phylum_Observed.png`
- etc.

**Beta Diversity Plots**:
- `pcoa_kingdom_bray_with_stats.png`
- `pcoa_phylum_jaccard_with_stats.png`
- etc.

**Stacked Bar Plots**:
- `stacked_barplot_phylum_top20.png`
- `stacked_barplot_genus_top20.png`
- etc.

### Understanding the Output

1. **Alpha Diversity**: Higher values generally indicate more diverse communities
2. **Beta Diversity**: PCoA plots show sample relationships - closer points are more similar
3. **Statistical Results**: 
   - PERMANOVA p < 0.05 indicates significant group differences
   - ANOSIM R values: 0 = no difference, 1 = complete separation
4. **Stacked Bar Plots**: Show relative abundance of most abundant taxa across samples

---

## Step 13: Differential Abundance Analysis with DESeq2

### Purpose
Identify significantly different microbial features between experimental groups using DESeq2's robust statistical framework.

### Script: `13_DESeq2_microbiome.R`

This script performs:
- Differential expression analysis adapted for microbiome data
- Hierarchical clustering and PCoA visualization
- Volcano plot generation
- Excel export of results

### Required R Packages

```r
# Install DESeq2 and dependencies
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

# Install other required packages
install.packages(c("ecodist", "RColorBrewer", "pheatmap", "vegan", "ggplot2"))
```

### Input Files Required

The script expects two input files:

1. **Count Data**: `inflammation_all_samples_kraken2-counts.tsv`
   - Rows = microbial features (taxa)
   - Columns = samples
   - Values = raw counts

2. **Metadata**: `metadata.txt`
   - Must contain columns: sample names, group information
   - Tab-separated format

### Create the Metadata File

```bash
# Create metadata file for the inflammation study
cat > metadata.txt << 'EOF'
	group	label
CA1E	01_Control	CA1E
CA2E	01_Control	CA2E
CA3E	01_Control	CA3E
CA4E	01_Control	CA4E
CA5E	01_Control	CA5E
TA1E	02_DSS	TA1E
TA2E	02_DSS	TA2E
TA3E	02_DSS	TA3E
TA4E	02_DSS	TA4E
TA5E	02_DSS	TA5E
EOF
```

### How to Run

```bash
# Make sure input files are in the working directory
ls inflammation_all_samples_kraken2-counts.tsv metadata.txt

# Run DESeq2 analysis
Rscript scripts/13_DESeq2_microbiome.R
```

### Key Features of the Script

1. **Data Preprocessing**:
   - Converts to integer counts (required by DESeq2)
   - Adds pseudocount of 1 to handle zeros
   - Creates DESeq2 data object

2. **Exploratory Analysis**:
   - Regularized log transformation (rlog) for visualization
   - Hierarchical clustering with Bray-Curtis distance
   - Principal Coordinates Analysis (PCoA)

3. **Differential Analysis**:
   - Uses geometric mean for size factor estimation
   - Fits negative binomial model
   - Multiple testing correction (Benjamini-Hochberg)

4. **Visualization**:
   - Sample distance heatmap
   - PCoA plot with group coloring
   - Volcano plot highlighting significant features

5. **Results Export**:
   - Multiple significance thresholds (p < 0.1, 0.05, 0.01)
   - Normalized count data included
   - Tab-separated output files

### Expected Output Files

**Statistical Results**:
- `inflammation_all_samples_kraken2-counts_p0.1_DE_features.tsv`
- `inflammation_all_samples_kraken2-counts_p0.05_DE_features.tsv`
- `inflammation_all_samples_kraken2-counts_p0.01_DE_features.tsv`
- `inflammation_all_samples_kraken2-counts_allRes_normData.tsv`

**Visualization Files**:
- `inflammation_all_samples_kraken2-counts_HCheatmap.pdf`
- `inflammation_all_samples_kraken2-counts_PCoA.pdf`
- `inflammation_all_samples_kraken2-counts_volcanoPlot.pdf`

### Interpreting Results

1. **Log2 Fold Change**: 
   - Positive values = higher in Treatment group
   - Negative values = higher in Control group
   - |log2FC| > 1 indicates 2-fold change

2. **Adjusted p-values (padj)**:
   - < 0.05 = statistically significant
   - Controls for multiple testing

3. **Volcano Plot Colors**:
   - Red: Significantly higher in treatment (padj < 0.05, log2FC > 1)
   - Green: Significantly lower in treatment (padj < 0.05, log2FC < -1)
   - Blue: Significant but small fold change
   - Gray: Large fold change but not significant
   - Black: Not significant

### Example Results Interpretation

```bash
# View top significant results
head -20 inflammation_all_samples_kraken2-counts_p0.05_DE_features.tsv

# Count significant features
wc -l inflammation_all_samples_kraken2-counts_p0.05_DE_features.tsv
```

---

## Step 14: Statistical Analysis with metagenomeSeq

### Purpose
Alternative statistical analysis using metagenomeSeq's zero-inflated Gaussian model, specifically designed for sparse microbiome data.

### Script: `14_metagenomeSeq.R`

This script provides:
- CSS (Cumulative Sum Scaling) normalization
- Zero-inflated Gaussian model fitting
- Advanced visualization options
- Multiple statistical approaches for microbiome data

### Required R Packages

```r
# Install metagenomeSeq and dependencies
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("metagenomeSeq")

# Install other packages
install.packages(c("tidyverse", "reshape2", "RColorBrewer", "pheatmap", "ggrepel"))
```

### How to Run

```bash
# Ensure input files are available
ls inflammation_all_samples_kraken2-counts.tsv metadata.txt

# Run metagenomeSeq analysis
Rscript scripts/14_metagenomeSeq.R
```

### Key Features of the Script

1. **Data Import and Setup**:
   - Loads count and metadata files
   - Creates MRexperiment object (metagenomeSeq format)
   - Validates data consistency

2. **Data Exploration**:
   - Basic statistics (features, samples, sparsity)
   - Data quality assessment

3. **Filtering**:
   - Removes low-abundance features
   - Configurable presence and depth thresholds
   - Reports filtering statistics

4. **Normalization**:
   - CSS (Cumulative Sum Scaling) normalization
   - Calculates and reports normalization factors
   - Generates normalized count matrix

5. **Statistical Analysis**:
   - Zero-inflated Gaussian (ZIG) model
   - Differential abundance testing
   - Multiple testing correction

6. **Advanced Visualization**:
   - PCA with variance explained
   - Heatmap of most variable features
   - Multiple volcano plot solutions for overlapping points
   - Box plots of significant features

### Unique Features of metagenomeSeq

1. **CSS Normalization**: Designed specifically for microbiome data sparsity
2. **Zero-Inflation Handling**: Explicitly models excess zeros common in microbiome data
3. **Robust Statistics**: Less sensitive to outliers than standard methods

### Expected Output Files

**Statistical Results**:
- `differential_abundance_results.csv` - Complete results table
- `normalized_counts.csv` - CSS-normalized count matrix

**Visualization Outputs**:
- PCA plot with group separation
- Heatmap of top 50 most variable features
- Multiple volcano plot variations
- Box plots of top significant features

### Volcano Plot Solutions

The script provides four different approaches to handle overlapping points in volcano plots:

1. **Jitter Method**: Adds small random offsets
2. **Position Jitter**: Uses ggplot2's position_jitter
3. **ggrepel Labels**: Adds non-overlapping labels
4. **Large Transparent Points**: Uses transparency to show density

### Interpreting metagenomeSeq Results

1. **Log2 Fold Change**: Similar to DESeq2 interpretation
2. **Adjusted p-values**: Multiple testing corrected p-values
3. **Taxa Information**: Includes full taxonomic classification
4. **Normalization Factors**: Shows sample-specific scaling factors

### Comparison with DESeq2

| Feature | DESeq2 | metagenomeSeq |
|---------|--------|---------------|
| Model | Negative Binomial | Zero-Inflated Gaussian |
| Normalization | Geometric Mean | CSS |
| Zero Handling | Implicit | Explicit |
| Best For | RNA-seq adapted | Microbiome-specific |

### Running Both Analyses

```bash
# Run both statistical approaches for comprehensive analysis
Rscript scripts/13_DESeq2_microbiome.R
Rscript scripts/14_metagenomeSeq.R

# Compare results
echo "DESeq2 significant features:"
wc -l inflammation_all_samples_kraken2-counts_p0.05_DE_features.tsv

echo "metagenomeSeq results available in:"
ls differential_abundance_results.csv
```

---

## Complete Pipeline Summary

### Quick Start Commands

```bash
# 1. Setup environment
conda activate shotgun2025
mkdir -p shotgun_analysis/{data,scripts,output,logs}
cd shotgun_analysis

# 2. Quality control
sbatch scripts/01_raw_fastqc_multiqc.sh

# 3. Read processing
sbatch scripts/02_cutadapt.sh
sbatch scripts/03_prinseq.sh
sbatch scripts/04_bmtagger.sh
sbatch scripts/05_bmtagger_fastqc_multiqc.sh

# 4. Taxonomic classification
sbatch scripts/06_kraken_bracken_5.sh

# 5. Functional profiling
sbatch scripts/07_humann_multiple_samples.sh
./scripts/08_postprocessing_humann3.sh

# 6. Data merging and filtering
python scripts/09_merge.py output/06_kraken2/*_mpa.tax > merged_taxonomy.tsv
./scripts/10_filter_by_average.sh merged_taxonomy.tsv 0.01 > filtered_taxonomy.tsv
python scripts/11_parse_taxa.py filtered_taxonomy.tsv

# 7. Statistical analysis
Rscript scripts/12_alpha-beta_diversity_norm.R
Rscript scripts/13_DESeq2_microbiome.R
Rscript scripts/14_metagenomeSeq.R
```

### Troubleshooting Common Issues

1. **Database Path Errors**: Update database paths in scripts to match your system
2. **Memory Issues**: Increase memory allocation in SLURM headers
3. **Conda Environment**: Ensure all required packages are installed
4. **File Permissions**: Make scripts executable with `chmod +x`
5. **Array Jobs**: Adjust array indices based on your sample count

### Key Output Files for Publication

1. **Quality Metrics**: MultiQC reports from steps 1 and 5
2. **Taxonomic Composition**: Krona plots from step 6
3. **Diversity Analysis**: All plots from step 12
4. **Statistical Results**: Significant features from steps 13-14
5. **Functional Profiles**: HUMAnN pathway abundances from steps 7-8

This complete pipeline provides a robust framework for shotgun metagenomics analysis, from raw reads to publication-ready results. Each step builds upon the previous ones, creating a comprehensive analysis workflow suitable for various microbiome research questions.
