# Script to use dfast-core to uniformly annotate bacterial genomes
# Written by TES in May, 2021


## NOTES
# User must first download the bacterial genomes listed in Table S9,
# then add unique headers to each file name 
# (i.e. change "GCF_000009605.1_ASM960v1_genomic.fna" to
# "BuapAcpi_GCF_000009605.1_ASM960v1_genomic.fna").

# User must also install dfast by following the instructions at:
# https://github.com/nigyta/dfast_core#installation

# Next, utilize this script within the same folder using the following command:
# bash step6_bacteria_genome_annotations.sh


## Modify headers for each fasta file
# Append file/assembly name to every header so that
# each sequence has a unique identifier out of all genome sequences
for file in *.f*
do
# Assign a custom species identifier based on file name
ID=$( echo $file | cut -f 1 -d "_")
# Insert identifier at the start of each sequence header,
# then remove everything following the first space
sed "s|>|>${ID}_|g" $file | cut -f 1 -d " " > ${ID}.header.fna
done


## Annotate genomes with dfast
# Perform annotation on the intermediate .header.fna files
for file in *.header.fna
# Extract species name from file name
do
species=$( echo $file | cut -f 1 -d ".")
dfast -g $file -o ${species}_dfast --force
done


## Identify and remove pseudogenes
# For each species, identify genes determined by dfast to contain
# internal stop codons
for file in `find -name genome.gff -print`
do
species=$( echo $file | cut -f 2 -d "/" | cut -f 1 -d "_" )
grep "internal stop codon" $file |
sed -r "s/.*(ID=\w+).*(locus_tag=\w+).*/\1-\2/g" | 
sed "s/ID=/${species}-/g" | sed "s/locus_tag=//g" > \
${species}_dfast/genes_with_internal_stop_codons.txt
done
# For each species, identify determined by dfast to contain
# frameshift mutations
for file in `find -name genome.gff -print`
do
species=$( echo $file | cut -f 2 -d "/" | cut -f 1 -d "_" ) 
grep "frameshift" $file |
sed -r "s/.*(ID=\w+).*(locus_tag=\w+).*/\1-\2/g" | 
sed "s/ID=/${species}-/g" | sed "s/locus_tag=//g" > \
${species}_dfast/genes_with_frameshifts.txt
done
