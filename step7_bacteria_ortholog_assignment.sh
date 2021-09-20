# Script to use dfast-core to uniformly annotate bacterial genomes
# Written by TES in June, 2021


## NOTES
# User must download the following protein annotations from the Ensembl website:
# Bacteriodes fragilis YCH46 (ASM992v1), Escherichia coli MG1655 (ASM584v2), 
# Buchnera from A. pisum (ASM960v1), B. pistaciae (ASM772v1),
# Schizaphis graminum (ASM736v1), and S. symbiotica from A. pisum str. Tucson (ASM18648v1)
# Place these files in a folder called bacteria_ensembl_proteomes


# User must also install orthofinder by following the instructions at:
# https://github.com/davidemms/OrthoFinder

# Next, utilize this script within the same folder using the following command:
# bash step7_bacteria_ortholog_assignment.sh


## Extract protein sequences from the dfast folders of each species
# Copy bacteria dfast protein annotations to a new folder
mkdir bacteria_dfast_protein 
for path in `find -name protein.faa -print`
do
species=$( echo $path | cut -f 2 -d "/" | cut -f 1 -d "_") 
cp $path bacteria_dfast_protein/
mv bacteria_dfast_protein/protein.faa bacteria_dfast_protein/${species}.faa
done


## Ensure protein files have unique sequence headers
cd bacteria_dfast_protein/
# Append the species name to the front of each protein sequence header
for file in *.faa
do
species=$( echo $file | cut -f 1 -d ".")
# In all current headers, replace | with -
# then insert species name at the start of each sequence header.
# Also replace any [ or ] in headers so as not to confuse grep with ranges
sed -i -e "s/|/-/g" -e "s|>|>${species}-|g" -e "s/\[/_/g" -e "s/\]/_/g" $file
done


## Remove sequences with stop codons and/or frameshifts
# First, extract sequences of genes with stop codons and/or frameshifts
for file in *.faa
do
species=$( echo $file | cut -f 1 -d ".")
grep -A 1 --no-group-separator \
-f ../${species}_dfast/genes_with_internal_stop_codons.txt $file > \
${species}_temp.faa 
grep -A 1 --no-group-separator \
-f ../${species}_dfast/genes_with_frameshifts.txt $file >> \
${species}_temp.faa
# Remove duplicate lines from pseudogenes file
awk '!visited[$0]++' ${species}_temp.faa > ${species}_pseudogenes.faa
# Second, remove these sequences from each protein.faa file
grep -v -f ${species}_pseudogenes.faa $file > ${species}_filt.faa 
done 
rm *_temp.faa


## Pre-processing of ensmbl proteomes
cd ../bacteria_ensembl_proteomes/
# Trim headers to include gene names and unique identifiers only
for file in *.pep.all.fa
do
species=$( echo $file | cut -f 1 -d ".")
if echo "$file" | grep -q "ych46"; then
ID=`echo "Bafr_Ensmbl"`
elif echo "$file" | grep -q "Buchnera_aphidicola_bcc"; then
ID=`echo "BuapCice_Ensmbl"`
elif echo "$file" | grep -q "acyrthosiphon_pisum"; then
ID=`echo "BuapAcpi_Ensmbl"`
elif echo "$file" | grep -q "baizongia_pistaciae"; then
ID=`echo "BuapBapi_Ensmbl"`
elif echo "$file" | grep -q "schizaphis_graminum"; then
ID=`echo "BuapScgr_Ensmbl"`
elif echo "$file" | grep -q "mg1655"; then
ID=`echo "Esco_Ensmbl"`
elif echo "$file" | grep -q "str_tucson"; then
ID=`echo "SesyAcpi_Ensmbl"`
fi
if grep -q "gene_symbol" $file; then
sed -r "s/^>(\w+)\s.*gene_symbol:(\w+)\s.*/>\2_\1/g" $file > ${ID}_filt.faa
else
sed -r "s/^>(\w+)\s.*/>\1/g" $file > ${ID}_filt.faa
fi
done
mv *_filt.faa ../bacteria_dfast_protein/


## Run OrthoFinder on filtered dfast and Ensembl annotations
cd ../bacteria_dfast_protein/
# Move all input files into their own directory
mkdir no_pseudogenes
mv *_filt.faa  no_pseudogenes/
# Run OrthoFinder 
ulimit -n 31076
orthofinder -t 36 -a 36 -I 10 -o ../bacteria_orthofinder/ -f no_pseudogenes/
