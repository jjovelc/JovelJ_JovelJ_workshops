#!/usr/bin/bash


#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=7-00:00:00
#SBATCH --mem=4G
#SBATCH --output=logs/QC_bmtagged_sbatch_job.%A.out
#SBATCH --error=logs/QC_bmtagged_sbatch_job.%A.err


#Innitalize conda
conda init bash &> /dev/null

#refresh shell environment after conda init
source ~/.bashrc &> /dev/null


##Required tools to install
#conda install -c conda-forge -c bioconda fastqc==0.12.1
#conda install pip
#pip install multiqc


echo "#########################"
echo "#        FASTQC         #"
echo "#########################"



echo "Started at: `date`" && echo -e "\n"

conda activate shotgun2025

bmtagger_out_dir="$(pwd)/output/04_bmtagger/"

forward_reads="_bmtagged_1.fastq"
reverse_reads="_bmtagged_2.fastq"

bmtagger_fastqc_dir="$(pwd)/output/05_bmtagger_fastqc/"

mkdir -p ${bmtagger_fastqc_dir}



# Run bracken on kraken reports
for sample in "${bmtagger_out_dir}"/*"$forward_reads"; do

    # Extract the sample name by removing the path and the _R1.fq suffix
    SAMPLE=$(basename "${sample}")
    SAMPLE=$(echo "$SAMPLE" | cut -d'_' -f1)
    echo "${SAMPLE}"

    echo ${bmtagger_out_dir}${SAMPLE}${forward_reads}
    echo ${bmtagger_out_dir}${SAMPLE}${reverse_reads}

    ### Running fastqc
    fastqc -o ${bmtagger_fastqc_dir} -t 20 ${bmtagger_out_dir}${SAMPLE}${forward_reads} ${bmtagger_out_dir}${SAMPLE}${reverse_reads}

done


echo "##########################"
echo "#        MULTIQC         #"
echo "##########################"


bmtagger_multiqc_dir="$(pwd)/output/05_bmtagger_multiqc"


mkdir -p ${bmtagger_multiqc_dir}

conda activate multiqc

### Running multiQC
multiqc -f ${bmtagger_fastqc_dir} -o ${bmtagger_multiqc_dir} -n multiqc_report.html

echo "All samples processed."

echo "Finished at: `date`" && echo -e "\n"

