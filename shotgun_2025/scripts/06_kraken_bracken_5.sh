#!/usr/bin/bash


#SBATCH --partition=synergy,cpu2019,cpu2021,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=7-00:00:00
#SBATCH --mem=200G
#SBATCH --output=logs/annotation_sbatch_job.%A.out
#SBATCH --error=logs/annotation_sbatch_job.%A.err



echo "Started at: `date`" && echo -e "\n"

#installing required tools
#conda install bioconda::kraken2
#conda install bioconda::bracken
#bioconda::krakentools

# please update as needed
eval "$(conda shell.bash hook)"
conda activate shotgun2025


# Run Kraken2 on all samples, but generating a metaphlan-like report (--use-mpa-style), which is desirable
# for some downstream applications.

# Directory paths (modify as needed)
bmtagger_out_dir="$(pwd)/output/04_bmtagger/"

Kraken_DIR="$(pwd)/output/06_kraken2"
Bracken_DIR="$(pwd)/output/06_bracken2"

# Ensure output directory exists
mkdir -p "${Kraken_DIR}"

# Path to the kraken2 database (modify accordingly)
DB="/bulk/IMCshared_bulk/sbagheri/workshops/spring_workshop_2025/final_files/kraken2_NCBI_Oct22"


echo -e "\n\n****** First Kraken step to produce files with mpa.tax suffix ******\n\n"


# Loop through R1 files in the data directory
for sample in "${bmtagger_out_dir}"/*_bmtagged_1.fastq; do

    # Extract the sample name by removing the path and the _R1.fq suffix
    SAMPLE=$(basename "${sample}" "_bmtagged_1.fastq")

    echo "${SAMPLE}: First alignment"

    R1=${SAMPLE}_bmtagged_1.fastq
    R2=${SAMPLE}_bmtagged_2.fastq
    echo "DB is set to: $DB"
    echo "Sample is set to: $SAMPLE"

    echo "R1 is set to: $R1"
    echo "R2 is set to: $R2"

    # Run kraken2
    kraken2 \
      --db "${DB}" \
      --threads 48 \
      --memory-mapping \
      --confidence 0.1 \
      --report "${Kraken_DIR}/${SAMPLE}_mpa.tax" \
      --paired "${bmtagger_out_dir}/${R1}" "${bmtagger_out_dir}/${R2}" \
      --output "${Kraken_DIR}/${SAMPLE}.krk2"

done



echo -e "\n\n****** Second Kraken step to produce files with krk2rpt suffix ******\n\n"


# Run Kraken2 on all samples, but this time generating a typical Kraken2 report that is compatible with
# Bracken (for quantification). This command is identical to the previous one, except for lacking the
# flag '--use-mpa-style'

for sample in "${bmtagger_out_dir}"/*_bmtagged_1.fastq; do

    # Extract sample name by removing the path and `_R1.fq` from the filename
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



echo -e "\n\n****** First Bracken step to produce files with bckrpt suffix ******\n\n"

# Ensure output directory exists
mkdir -p "${Bracken_DIR}"


# Run bracken on kraken reports
for sample in "${Kraken_DIR}"/*krk2rpt; do

    # Extract the sample name by removing the path and the _R1.fq suffix
    SAMPLE=$(basename "${sample}")
    echo "${SAMPLE}"

    bracken \
    -d $DB \
    -i "${Kraken_DIR}/$SAMPLE" \
    -r 300 \
    -t 100 \
    -o "${Bracken_DIR}"/${SAMPLE/.krk2rpt/.bck} \
    -w "${Bracken_DIR}"/${SAMPLE/.krk2rpt/.bckrpt}

done

echo ${Bracken_DIR}



echo -e "\n\n****** Visualization of Bracken results ******\n\n"

Kraken_bracken_figures="$(pwd)/output/06_Kr-Br-html_reports"


# Ensure output directory exists
mkdir -p "${Kraken_bracken_figures}"


## Run kreport2krona.py
for sample in "${Bracken_DIR}"/*bckrpt; do kreport2krona.py -r $sample -o ${sample/.bckrpt/_krona.txt} --no-intermediate-ranks; done

## Produce individual Krona plots for individual samples
for sample in "${Bracken_DIR}"/*krona.txt; do ktImportText $sample -o "${Kraken_bracken_figures}"/${sample/txt/html}; done

## Produce s single HTML document containing Krona plots for all samples
ktImportText "${Bracken_DIR}"/*krona.txt -o "${Kraken_bracken_figures}"/all_krona_plots.html


echo "All samples processed."

echo "Finished at: `date`" && echo -e "\n"

