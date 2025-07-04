#!/usr/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=7-00:00:00
#SBATCH --mem=4G
#SBATCH --output=logs/cutadapt_sbatch_job.%A.out
#SBATCH --error=logs/cutadapt_sbatch_job.%A.err


#Innitalize conda
conda init bash &> /dev/null

#refresh shell environment after conda init
source ~/.bashrc &> /dev/null


echo "Started at: `date`" && echo -e "\n"

conda activate shotgun2025

minlength=60
maxn=15
cpus=7

# Define directories and filename parts
data_dir="/work/vetmed_shared_dbs/shotgun_workshop_2025/data/"
forward_reads="_R1_100K.fq"
reverse_reads="_R2_100K.fq"
cutadapt_out_dir="$(pwd)/output/02_cutadapt"


mkdir -p ${cutadapt_out_dir}

# Loop over forward read files in the data directory
for file in "${data_dir}"/*"${forward_reads}"; do
    # Extract the basename (filename only)
    filename=$(basename "${file}")

    # Remove the suffix to get the sample name (e.g., "dfdrt" from "dfdrt_R1.fq")
    sample_name="${filename%${forward_reads}}"

    echo "Sample name: ${sample_name}"  # Example output: dfdrt

    echo "Processing sample: ${sample_name}"

    # Construct full paths for the input files using the sample name
    input1="${data_dir}/${sample_name}${forward_reads}"
    input2="${data_dir}/${sample_name}${reverse_reads}"

    # Construct full paths for the output files (note the underscore added before trimmed/untrimmed)
    output1="${cutadapt_out_dir}/${sample_name}_trimmed_R1.fq"
    output2="${cutadapt_out_dir}/${sample_name}_trimmed_R2.fq"

    echo "Input1: ${input1}"
    echo "Input2: ${input2}"
    echo "Output1: ${output1}"
    echo "Output2: ${output2}"

    cutadapt -O 19 -e 0.15 -m ${minlength} --max-n ${maxn} -g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -G GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -j ${cpus} -o ${output1} -p ${output2} ${input1} ${input2}
done

echo "Adapter removal finished at: `date`" && echo -e "\n"


cutadapt_Qtrim_dir="$(pwd)/output/02_cutadaptQC"

forward_reads="_trimmed_R1.fq"
reverse_reads="_trimmed_R2.fq"

mkdir -p ${cutadapt_Qtrim_dir}


# Loop over forward read files in the data directory
for file in "${cutadapt_out_dir}"/*"${forward_reads}"; do
    # Extract the basename (filename only)
    filename=$(basename "${file}")

    # Remove the suffix to get the sample name (e.g., "dfdrt" from "dfdrt_R1.fq")
    sample_name="${filename%${forward_reads}}"

    echo "Sample name: ${sample_name}"  # Example output: dfdrt

    echo "Processing sample: ${sample_name}"

    # Construct full paths for the input files using the sample name
    input1="${cutadapt_out_dir}/${sample_name}${forward_reads}"
    input2="${cutadapt_out_dir}/${sample_name}${reverse_reads}"

    # Construct full paths for the output files (note the underscore added before trimmed/untrimmed)
    output1="${cutadapt_Qtrim_dir}/${sample_name}${forward_reads}"
    output2="${cutadapt_Qtrim_dir}/${sample_name}${reverse_reads}"

    echo "Input1: ${input1}"
    echo "Input2: ${input2}"
    echo "Output1: ${output1}"
    echo "Output2: ${output2}"

    cutadapt --nextseq-trim=20 --poly-a -m ${minlength} -o ${output1} -p ${output2} ${input1} ${input2}

done

echo "All samples processed."

echo "Quality trimming finished at: `date`" && echo -e "\n"
