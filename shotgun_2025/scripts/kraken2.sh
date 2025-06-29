#!/usr/bin/bash

start=`date +%s`

# please update as needed
eval "$(conda shell.bash hook)"
conda activate xxxx

# I can pass you the database, or you try to compile
# your own one.
DB=/home/dayhoff/sdd/juan/db/kraken2_NCBI_Oct22

# Erase previous symbolic liks if they exist
rm /dev/shm/hash.k2d
rm /dev/shm/opts.k2d
rm /dev/shm/taxo.k2d

# Create symbolic links to the new Kraken2 database files
ln -s ${DB}/hash.k2d /dev/shm/hash.k2d
ln -s ${DB}/opts.k2d /dev/shm/opts.k2d
ln -s ${DB}/taxo.k2d /dev/shm/taxo.k2d

# Run Kraken2 on all samples, but generating a metaphlan-like report (--use-mpa-style), which is desirable
# for some downstream applications.
for FILE in *_R1_001_q30_5M.fastq ; do echo "$FILE: First alignment"; kraken2 --db $DB --threads 48 --use-mpa-style --memory-mapping --confidence 0.1 --report ${FILE/_R1*/}_mpa.tax $FILE ${FILE/_R1/_R2} --output ${FILE/_R1*/}.krk2; done

# Remove Kraken2 standard output fron the previous run since it will run again.
rm *_kraken2.txt

# Run Kraken2 on all samples, but this time generating a typical Kraken2 report that is compatible with
# Bracken (for quantification). This command is identical to the previous one, except for lacking the
# flag '--use-mpa-style'
for FILE in *_R1_001_q30_5M.fastq; do echo "$FILE: Second alignment"; kraken2 --db $DB --threads 48 --memory-mapping --confidence 0.1 --report ${FILE/_R1*/}.krk2rpt --paired ${FILE} ${FILE/_R1/_R2} --output ${FILE/_R1*/}.krk2; done

rm -rf k2
mkdir k2_reports
mkdir k2_outputs

mv *krk2rpt k2_reports
mv *krk2 k2_outputs
mv *_mpa.tax k2_outputs 

cd k2_reports

# Run bracken on kraken reports
for FILE in *krk2rpt; do bracken -d $DB -i $FILE -r 300 -t 100 -o ${FILE/.krk2rpt/.bck} -w ${FILE/.krk2rpt/.bckrpt}; done

## Run kreport2krona.py
for FILE in *bckrpt; do kreport2krona.py -r $FILE -o ${FILE/.bckrpt/_krona.txt} --no-intermediate-ranks; done

## Produce individual Krona plots for individual samples
for FILE in *krona.txt; do ktImportText $FILE -o ${FILE/txt/html}; done

## Produce s single HTML document containing Krona plots for all samples
ktImportText *krona.txt -o all_krona_plots.html

end=`date +%s`
runtime=$((end-start))
echo $runtime
