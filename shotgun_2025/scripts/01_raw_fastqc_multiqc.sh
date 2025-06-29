#!/usr/bin/bash


#SBATCH --partition=synergy,cpu2019,cpu2021,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=7-00:00:00
#SBATCH --mem=200G
#SBATCH --output=logs/QC_sbatch_job.%A.out
#SBATCH --error=logs/QC_sbatch_job.%A.err


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

data_dir="/work/vetmed_shared_dbs/shotgun_workshop_2025/shotgun_inflammation_libraries/100K/"
forward_reads="_R1_100K.fq"
reverse_reads="_R2_100K.fq"
fastqc_out_dir="$(pwd)/output/01_raw_fastqc"

mkdir -p ${fastqc_out_dir}



# Run bracken on kraken reports
for sample in "${data_dir}"/*"$forward_reads"; do

    # Extract the sample name by removing the path and the _R1.fq suffix
    SAMPLE=$(basename "${sample}")
    SAMPLE=$(echo "$SAMPLE" | cut -d'_' -f1)
    echo "${SAMPLE}"

    echo ${data_dir}${SAMPLE}${forward_reads}
    echo ${data_dir}${SAMPLE}${reverse_reads}

    ### Running fastqc
    fastqc -o ${fastqc_out_dir} -t 20 ${data_dir}${SAMPLE}${forward_reads} ${data_dir}${SAMPLE}${reverse_reads}

done


echo "##########################"
echo "#        MULTIQC         #"
echo "##########################"


multiqc_out_dir="$(pwd)/output/01_raw_multiqc"


mkdir -p ${multiqc_out_dir}

### Running multiQC
multiqc -f ${fastqc_out_dir} -o ${multiqc_out_dir} -n multiqc_report.html

echo "All samples processed."

echo "Finished at: `date`" && echo -e "\n"

