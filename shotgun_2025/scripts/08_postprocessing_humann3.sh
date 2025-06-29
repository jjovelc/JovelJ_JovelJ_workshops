#!/usr/bin/bash

for FILE in *genefamilies.tsv
do
	humann_renorm_table --input $FILE \
	    --output ${FILE/.tsv/_cpm.tsv} --units cpm --update-snames
done



echo "______________________________________"
echo " Regrouping"
echo " Processing normalized data"
for FILE in *_cpm.tsv
do
	humann_regroup_table --input $FILE \
	    --output ${FILE/.tsv/_rxn.tsv} --groups uniref90_rxn
done

echo "__________"

echo " Processing raw data"
for FILE in *_genefamilies.tsv
do
	humann_regroup_table --input $FILE \
	    --output ${FILE/.tsv/_rxn.tsv} --groups uniref90_rxn
done



echo "______________________________________"
echo "Renaming"
echo "Processing normalized data"

for FILE in *_cpm_rxn.tsv
do
	 humann_rename_table --input $FILE \
		     --output ${FILE/.tsv/_ren.tsv} --names metacyc-rxn
done

echo "__________"

echo "Processing raw data"

for FILE in *_genefamilies_rxn.tsv
do
	 humann_rename_table --input $FILE \
		     --output ${FILE/.tsv/_ren.tsv} --names metacyc-rxn
done

