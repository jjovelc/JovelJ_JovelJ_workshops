#!/usr/bin/bash

# Start a for cycle that iterates along all kallisto result folders
for DIR in *_kallisto
do
# Grab the first column of each abundance file and save it in a file called names
# This will be done 20 times but only the last iteration would be kept.
	cut -f 1 "${DIR}"/abundance.tsv > names
# Grab est.counts from abundance.tsv files, delete the title and the add a new title
# containing the name of the sample = ${DIR/_kallisto/}
	cut -f 4 "${DIR}"/abundance.tsv | sed 1d | sed '1s/^/'${DIR/_kallisto/}'\n/' > ${DIR/_kallisto/_counts}
	cut -f 5 "${DIR}"/abundance.tsv | sed 1d | sed '1s/^/'${DIR/_kallisto/}'\n/' > ${DIR/_kallisto/_tpms}
done

# Pasta names plus all counts
paste names *counts > all_samples_counts.tsv
# Pasta names plus all tpms
paste names *tpms > all_samples_tpms.tsv

# Remoce temporary files
rm names *counts *tpms

# Inform users where results have been saved
echo "Your count results are in file all_samples_counts.tsv"
echo "Your TPM results are in file all_samples_tpms.tsv"

