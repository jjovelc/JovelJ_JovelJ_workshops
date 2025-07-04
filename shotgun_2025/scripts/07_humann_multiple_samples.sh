#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=0-10:00:00
#SBATCH --mem=32G
#SBATCH --error=logs/humann_run.%J.err
#SBATCH --output=logs/humann_run.%J.out
#SBATCH --array=0-9 # job array index



# please update as needed
eval "$(conda shell.bash hook)"
conda activate shotgun2025

bmtagger_out_dir="$(pwd)/output/04_bmtagger/"

# Path to the humann database (modify accordingly)
DB="/work/vetmed_shared_dbs/humann3_dbs"

threads=20
nuc_db="${DB}/chocophlan"
prot_db="${DB}/uniref"


forward_reads="_bmtagged_1.fastq"
reverse_reads="_bmtagged_2.fastq"

####### Merging reads #######
names=$(for FILE in *_R1_100K.fq; do echo "${FILE/_R1_100K.fq/}"; done)


SAMPLE=${names[${SLURM_ARRAY_TASK_ID}]}
echo ${SAMPLE}

merged_dir="$(pwd)/output/07_merged_reads/output_${SAMPLE}"

# Ensure output directory exists
mkdir -p "${merged_dir}"


echo "Starting to merge reads"

# Construct full paths for the input files using the sample name
input1="${bmtagger_out_dir}/${SAMPLE}${forward_reads}"
input2="${bmtagger_out_dir}/${SAMPLE}${reverse_reads}"

output="${merged_dir}/${SAMPLE}_merged.fastq"
cat ${input1} ${input2} > ${output}

echo "All samples were merged"

echo "Starting to run humann"


# Directory paths (modify as needed)
humann_out_dir="$(pwd)/output/07_humann/"

mkdir -p "${humann_out_dir}"

echo "Started at: `date`" && echo -e "\n"


humann --input ${merged_dir}/${SAMPLE}_merged.fastq --output ${SAMPLE} --threads 6 \
       --protein-database ${prot_db} \
       --nucleotide-database ${nuc_db} \
       --metaphlan-options=" --index mpa_vJun23_CHOCOPhlAnSGB_202307 --bowtie2db /bulk/IMCbinf_bulk/sbagheri/Projects_IMC/databases/metaphlan/vJun23_version" \
       --output ${humann_out_dir}/${SAMPLE}_humann3_res
:


echo "Finished at: `date`" && echo -e "\n"

