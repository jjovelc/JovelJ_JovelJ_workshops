#!/usr/bin/bash


#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=7-00:00:00
#SBATCH --mem=4G
#SBATCH --output=logs/prinseq_sbatch_job.%A.out
#SBATCH --error=logs/prinseq_sbatch_job.%A.err


#Innitalize conda
conda init bash &> /dev/null

#refresh shell environment after conda init
source ~/.bashrc &> /dev/null


echo "Started at: `date`" && echo -e "\n"


conda activate shotgun2025


echo "##########################"
echo "#        PRINSEQ         #"
echo "##########################"


cutadapt_Qtrim_dir="$(pwd)/output/02_cutadaptQC"

forward_reads="_trimmed_R1.fq"
reverse_reads="_trimmed_R2.fq"

prinseq_out_dir="$(pwd)/output/03_prinseq"

mkdir -p ${prinseq_out_dir}



### Prinseq parameters
lc_method="dust"
lc_threshold=7

# Minimum length for reads
minlength=60

# Maximum number of N bases allowed
maxn=15


# Loop over forward read files in the input directory
for file in ${cutadapt_Qtrim_dir}/*${forward_reads}; do
    # Extract the basename (filename only)
    filename=$(basename "${file}")

    # Remove the suffix to get the sample name (e.g., "dfdrt" from "dfdrt_R1.fq")
    sample_name="${filename%${forward_reads}}"

    echo "Sample name: ${sample_name}"  # Example output: dfdrt

    echo "Processing sample: ${sample_name}"

    # Construct full paths for the input files using the sample name
    input1="${cutadapt_Qtrim_dir}/${sample_name}${forward_reads}"
    input2="${cutadapt_Qtrim_dir}/${sample_name}${reverse_reads}"

    # Construct full paths for the output files (note the underscore added before bmtagged)
    output="${prinseq_out_dir}/${sample_name}_filtered"

    echo "-------------------------------------"
    echo "Sample:      ${sample_name}"
    echo "Forward:     ${input1}"
    echo "Reverse:     ${input2}"
    echo "-------------------------------------"


    perl $(pwd)/scripts/prinseq-lite.pl -fastq ${input1} -fastq2 ${input2} \
    -out_good ${output} -out_bad null -lc_method ${lc_method} -lc_threshold ${lc_threshold} \
    -derep 1 -min_len ${minlength} -ns_max_n ${maxn}

    echo "Prinseq finished for sample: ${sample_name}"
    echo

done

echo "All samples processed."

echo "Script finished at: `date`" && echo -e "\n"


