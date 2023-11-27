#!/usr/bin/bash

module load kallisto/0.46.1
INDEX=Homo_sapiens.GRCh38.cdna.all.idx

for FILE in *_1.fq.gz
do
	kallisto quant -i "$INDEX" -o ${FILE/_*/}_kallisto --bias $FILE ${FILE/_1.fq/_2.fq}

done	
