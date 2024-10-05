#!/usr/bin/bash

kallisto_path=$(which kallisto)

# If you need to download indices available from the kallisto repo
# available indices up to 240117 are:
# human, mouse, dog, monkey, zebrafish
# For example, for mouse:
kb ref -d mouse -i index.idx -g t2g.txt --kallisto $kallisto_path

# If an index other than the five listed above needs to be generated
# the following command can be used (with the relevant modifications)
kb ref -i output_index_kallisto.idx -g output_t2g.txt -f1 output_transcriptome.fa genome.fa  genes.gtf

# To generate text output format just remove --h5ad 
# --h5ad can also be replaced by 
# (data.h5ad file will be in the subdirectory counts_unfiltered)
kb count -i index.idx -g t2g.txt -x 10xv2 --h5ad -t 4 -o kb_out_h5ad SRR8599150_R1.fastq.gz SRR8599150_R2.fastq.gz --kallisto $kallisto_path



