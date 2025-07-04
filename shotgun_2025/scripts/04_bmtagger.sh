#!/usr/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=7-00:00:00
#SBATCH --mem=24G
#SBATCH --output=logs/bmtagger_sbatch_job.%A.out
#SBATCH --error=logs/bmtagger_sbatch_job.%A.err

#Innitalize conda
conda init bash &> /dev/null

#refresh shell environment after conda init
source ~/.bashrc &> /dev/null

echo "Started at: `date`" && echo -e "\n"

conda activate shotgun2025

shotgun_dir="/work/vetmed_shared_dbs/shotgun_workshop_2025"

# Index for bmfilter (part of bmtagger), bitmask file
bmfilter_ref="${shotgun_dir}/bmtagger_DB/mice_reference.bitmask"

# Index for srprism (part of bmtagger)
srprism_ref="${shotgun_dir}/bmtagger_DB/mice_reference.srprism"

prinseq_out_dir="$(pwd)/output/03_prinseq"
forward_reads="_filtered_1.fastq"
reverse_reads="_filtered_2.fastq"
bmtagger_out_dir="$(pwd)/output/04_bmtagger"

mkdir -p ${bmtagger_out_dir}

# Loop over forward read files in the data directory
for file in ${prinseq_out_dir}/*${forward_reads}; do
    # Extract the basename (filename only)
    filename=$(basename "${file}")

    # Remove the suffix to get the sample name (e.g., "dfdrt" from "dfdrt_R1.fq")
    sample_name="${filename%${forward_reads}}"

    echo "Sample name: ${sample_name}"  # Example output: dfdrt

    echo "Processing sample: ${sample_name}"

    # Construct full paths for the input files using the sample name
    input1="${prinseq_out_dir}/${sample_name}${forward_reads}"
    input2="${prinseq_out_dir}/${sample_name}${reverse_reads}"

    # Construct full paths for the output files (note the underscore added before bmtagged)
    output="${bmtagger_out_dir}/${sample_name}_bmtagged"

    echo "-------------------------------------"
    echo "Sample:      ${sample_name}"
    echo "Forward:     ${input1}"
    echo "Reverse:     ${input2}"
    echo "Output: ${output}"
    echo "-------------------------------------"

    # Run BMTagger
    bmtagger.sh \
        -b "${bmfilter_ref}" \
        -x "${srprism_ref}" \
        -q 1 \
        -1 "${input1}" \
        -2 "${input2}" \
        -o "${output}" \
        -X;

    echo "BMTagger finished for sample: ${sample_name}"
    echo
done

echo "All samples processed."

echo "Script finished at: `date`" && echo -e "\n"


