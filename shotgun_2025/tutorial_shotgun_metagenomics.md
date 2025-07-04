# Shotgun Metagenomics Analysis Pipeline Tutorial

## Overview

This tutorial provides a complete, step-by-step guide for analyzing shotgun metagenomics data. The pipeline takes you from raw sequencing reads through quality control, taxonomic classification, functional profiling, and statistical analysis. All commands can be copied and pasted directly into your terminal.

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
13. **Differential Abundance (DESeq2)** - Statistical testing with R
14. **Statistical Analysis (metagenomeSeq)** - Alternative statistical approach

---

## Environment Setup

### 1. Create Conda Environment

First, create the main analysis environment:

```bash
# Create environment from YAML file
conda env create -f shotgun2025.yaml

# Activate environment
conda activate shotgun2025
```

### 2. Directory Structure

Set up your working directory:

```bash
# Create main project directory
mkdir -p shotgun_analysis
cd shotgun_analysis

# Create subdirectories
mkdir -p {data,scripts,output,logs,dbs}
mkdir -p output/{01_raw_fastqc,01_raw_multiqc,02_cutadapt,02_cutadaptQC,03_prinseq,04_bmtagger,05_bmtagger_fastqc,05_bmtagger_multiqc,06_kraken2,06_bracken2,06_Kr-Br-html_reports,07_merged_reads,07_humann}
```

### 3. Required Software Installation

Install required tools in your conda environment:

```bash
# Activate environment
conda activate shotgun2025

# Install core tools
conda install -c conda-forge -c bioconda fastqc=0.12.1
conda install -c conda-forge -c bioconda cutadapt
conda install -c conda-forge -c bioconda kraken2
conda install -c conda-forge -c bioconda bracken
conda install -c conda-forge -c bioconda humann
conda install -c conda-forge -c bioconda bmtagger

# Install Python packages
pip install multiqc

# For R analysis
conda install -c conda-forge -c bioconda r-base r-essentials
```

---

## Step 1: Initial Quality Control

### Purpose
Assess the quality of raw sequencing reads to identify potential issues before processing.

### Script: `01_raw_fastqc_multiqc.sh`

```bash
#!/usr/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=7-00:00:00
#SBATCH --mem=4G
#SBATCH --output=logs/QC_sbatch_job.%A.out
#SBATCH --error=logs/QC_sbatch_job.%A.err

# Initialize conda
conda init bash &> /dev/null
source ~/.bashrc &> /dev/null

echo "Starting FastQC analysis..."
echo "Started at: `date`"

conda activate shotgun2025

# Define paths
data_dir="./data/"
forward_reads="_R1_100K.fq"
reverse_reads="_R2_100K.fq"
fastqc_out_dir="$(pwd)/output/01_raw_fastqc"

mkdir -p ${fastqc_out_dir}

# Run FastQC on all samples
for sample in "${data_dir}"/*"$forward_reads"; do
    SAMPLE=$(basename "${sample}")
    SAMPLE=$(echo "$SAMPLE" | cut -d'_' -f1)
    echo "Processing sample: ${SAMPLE}"

    echo "${data_dir}${SAMPLE}${forward_reads}"
    echo "${data_dir}${SAMPLE}${reverse_reads}"

    # Run FastQC
    fastqc -o ${fastqc_out_dir} -t 20 ${data_dir}${SAMPLE}${forward_reads} ${data_dir}${SAMPLE}${reverse_reads}
done

echo "Running MultiQC to aggregate results..."

multiqc_out_dir="$(pwd)/output/01_raw_multiqc"
mkdir -p ${multiqc_out_dir}

conda activate multiqc

# Run MultiQC
multiqc -f ${fastqc_out_dir} -o ${multiqc_out_dir} -n multiqc_report.html

echo "All samples processed."
echo "Finished at: `date`"
```

### How to Run

```bash
# Copy script to scripts directory
cp 01_raw_fastqc_multiqc.sh scripts/

# Submit to SLURM
sbatch scripts/01_raw_fastqc_multiqc.sh

# Check job status
squeue -u $USER

# View results
firefox output/01_raw_multiqc/multiqc_report.html
```

### Expected Output
- Individual FastQC reports for each sample
- Aggregated MultiQC report showing quality metrics across all samples
- HTML report viewable in web browser

---

## Step 2: Adapter Removal

### Purpose
Remove sequencing adapters and perform quality trimming to improve data quality.

### Script: `02_cutadapt.sh`

```bash
#!/usr/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=7-00:00:00
#SBATCH --mem=4G
#SBATCH --output=logs/cutadapt_sbatch_job.%A.out
#SBATCH --error=logs/cutadapt_sbatch_job.%A.err

# Initialize conda
conda init bash &> /dev/null
source ~/.bashrc &> /dev/null

echo "Starting adapter removal with Cutadapt..."
echo "Started at: `date`"

conda activate shotgun2025

# Parameters
minlength=60
maxn=15
cpus=7

# Define directories
data_dir="./data/"
forward_reads="_R1_100K.fq"
reverse_reads="_R2_100K.fq"
cutadapt_out_dir="$(pwd)/output/02_cutadapt"

mkdir -p ${cutadapt_out_dir}

# Loop over forward read files
for file in "${data_dir}"/*"${forward_reads}"; do
    filename=$(basename "${file}")
    sample_name="${filename%${forward_reads}}"

    echo "Processing sample: ${sample_name}"

    # Construct file paths
    input1="${data_dir}/${sample_name}${forward_reads}"
    input2="${data_dir}/${sample_name}${reverse_reads}"
    output1="${cutadapt_out_dir}/${sample_name}_trimmed_R1.fq"
    output2="${cutadapt_out_dir}/${sample_name}_trimmed_R2.fq"

    echo "Input1: ${input1}"
    echo "Input2: ${input2}"
    echo "Output1: ${output1}"
    echo "Output2: ${output2}"

    # Run Cutadapt for adapter removal
    cutadapt -O 19 -e 0.15 -m ${minlength} --max-n ${maxn} \
        -g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
        -G GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
        -j ${cpus} -o ${output1} -p ${output2} ${input1} ${input2}
done

echo "Adapter removal finished at: `date`"

# Quality trimming step
echo "Starting quality trimming..."

cutadapt_Qtrim_dir="$(pwd)/output/02_cutadaptQC"
forward_reads="_trimmed_R1.fq"
reverse_reads="_trimmed_R2.fq"

mkdir -p ${cutadapt_Qtrim_dir}

# Loop over trimmed files for quality trimming
for file in "${cutadapt_out_dir}"/*"${forward_reads}"; do
    filename=$(basename "${file}")
    sample_name="${filename%${forward_reads}}"

    echo "Quality trimming sample: ${sample_name}"

    # Construct file paths
    input1="${cutadapt_out_dir}/${sample_name}${forward_reads}"
    input2="${cutadapt_out_dir}/${sample_name}${reverse_reads}"
    output1="${cutadapt_Qtrim_dir}/${sample_name}${forward_reads}"
    output2="${cutadapt_Qtrim_dir}/${sample_name}${reverse_reads}"

    # Run quality trimming
    cutadapt --nextseq-trim=20 --poly-a -m ${minlength} \
        -o ${output1} -p ${output2} ${input1} ${input2}
done

echo "All samples processed."
echo "Quality trimming finished at: `date`"
```

### How to Run

```bash
# Submit job
sbatch scripts/02_cutadapt.sh

# Monitor progress
tail -f logs/cutadapt_sbatch_job.*.out
```

### Parameters Explained
- `-O 19`: Minimum overlap between adapter and read
- `-e 0.15`: Maximum error rate (15%)
- `-m 60`: Minimum read length after trimming
- `--max-n 15`: Maximum number of N bases allowed
- `--nextseq-trim=20`: Quality trimming from 3' end

---

## Step 3: Low Complexity Filtering

### Purpose
Remove low-complexity sequences and apply additional quality filters using Prinseq.

### Script: `03_prinseq.sh`

```bash
#!/usr/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=7-00:00:00
#SBATCH --mem=4G
#SBATCH --output=logs/prinseq_sbatch_job.%A.out
#SBATCH --error=logs/prinseq_sbatch_job.%A.err

# Initialize conda
conda init bash &> /dev/null
source ~/.bashrc &> /dev/null

echo "Starting Prinseq filtering..."
echo "Started at: `date`"

conda activate shotgun2025

# Input directory
cutadapt_Qtrim_dir="$(pwd)/output/02_cutadaptQC"
forward_reads="_trimmed_R1.fq"
reverse_reads="_trimmed_R2.fq"
prinseq_out_dir="$(pwd)/output/03_prinseq"

mkdir -p ${prinseq_out_dir}

# Prinseq parameters
lc_method="dust"
lc_threshold=7
minlength=60
maxn=15

# Loop over quality-trimmed files
for file in ${cutadapt_Qtrim_dir}/*${forward_reads}; do
    filename=$(basename "${file}")
    sample_name="${filename%${forward_reads}}"

    echo "Processing sample: ${sample_name}"

    # Construct file paths
    input1="${cutadapt_Qtrim_dir}/${sample_name}${forward_reads}"
    input2="${cutadapt_Qtrim_dir}/${sample_name}${reverse_reads}"
    output="${prinseq_out_dir}/${sample_name}_filtered"

    echo "-------------------------------------"
    echo "Sample:      ${sample_name}"
    echo "Forward:     ${input1}"
    echo "Reverse:     ${input2}"
    echo "-------------------------------------"

    # Run Prinseq
    perl $(pwd)/scripts/prinseq-lite.pl -fastq ${input1} -fastq2 ${input2} \
        -out_good ${output} -out_bad null \
        -lc_method ${lc_method} -lc_threshold ${lc_threshold} \
        -derep 1 -min_len ${minlength} -ns_max_n ${maxn}

    echo "Prinseq finished for sample: ${sample_name}"
    echo
done

echo "All samples processed."
echo "Script finished at: `date`"
```

### How to Run

```bash
# Make sure prinseq-lite.pl is in scripts directory
cp prinseq-lite.pl scripts/

# Submit job
sbatch scripts/03_prinseq.sh
```

### Parameters Explained
- `lc_method dust`: Low complexity detection method
- `lc_threshold 7`: Threshold for complexity filtering
- `derep 1`: Remove exact duplicates
- `min_len 60`: Minimum sequence length
- `ns_max_n 15`: Maximum N bases allowed

---

## Step 4: Host DNA Removal

### Purpose
Remove host contamination from microbiome samples using BMTagger.

### Script: `04_bmtagger.sh`

```bash
#!/usr/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=7-00:00:00
#SBATCH --mem=24G
#SBATCH --output=logs/bmtagger_sbatch_job.%A.out
#SBATCH --error=logs/bmtagger_sbatch_job.%A.err

# Initialize conda
conda init bash &> /dev/null
source ~/.bashrc &> /dev/null

echo "Starting host removal with BMTagger..."
echo "Started at: `date`"

conda activate shotgun2025

shotgun_dir="/work/vetmed_shared_dbs/shotgun_workshop_2025"

# BMTagger database paths
bmfilter_ref="${shotgun_dir}/bmtagger_DB/mice_reference.bitmask"
srprism_ref="${shotgun_dir}/bmtagger_DB/mice_reference.srprism"

# Input/output directories
prinseq_out_dir="$(pwd)/output/03_prinseq"
forward_reads="_filtered_1.fastq"
reverse_reads="_filtered_2.fastq"
bmtagger_out_dir="$(pwd)/output/04_bmtagger"

mkdir -p ${bmtagger_out_dir}

# Loop over filtered files
for file in ${prinseq_out_dir}/*${forward_reads}; do
    filename=$(basename "${file}")
    sample_name="${filename%${forward_reads}}"

    echo "Processing sample: ${sample_name}"

    # Construct file paths
    input1="${prinseq_out_dir}/${sample_name}${forward_reads}"
    input2="${prinseq_out_dir}/${sample_name}${reverse_reads}"
    output="${bmtagger_out_dir}/${sample_name}_bmtagged"

    echo "-------------------------------------"
    echo "Sample:      ${sample_name}"
    echo "Forward:     ${input1}"
    echo "Reverse:     ${input2}"
    echo "Output:      ${output}"
    echo "-------------------------------------"

    # Run BMTagger
    bmtagger.sh \
        -b "${bmfilter_ref}" \
        -x "${srprism_ref}" \
        -q 1 \
        -1 "${input1}" \
        -2 "${input2}" \
        -o "${output}" \
        -X

    echo "BMTagger finished for sample: ${sample_name}"
    echo
done

echo "All samples processed."
echo "Script finished at: `date`"
```

### How to Run

```bash
# Note: Requires BMTagger database - update paths as needed
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
