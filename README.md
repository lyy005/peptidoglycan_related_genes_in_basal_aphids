# Analytical pipeline of the 

This document is a walkthrough of the methods and code used to analyze the peptidoglycan related genes in aphid genomes and Bunchnera genomes.

## 1 - Aphid genome assembly

Raw data can be downloaded from NCBI BioProject (PRJNA758084). 

### 1.1 - Quality control: 

We cleaned both genomic data and RNA-seq data using Trimmomatic version 0.38. For example, we used the following commands for the Geopemphigus genome: 

        java -jar trimmomatic-0.38.jar PE -phred33 Geo_1.fastq.gz Geo_2.fastq.gz Geo_pe.1.fq.gz Geo_se.1.fq.gz Geo_pe.2.fq.gz Geo_se.2.fq.gz ILLUMINACLIP:/adapters/combined.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  

### 1.2 - Genome assembly using SPADES version 3.14.1

We first assembled draft genomes based on the genomic data using SPADES: 

  spades.py -1 Geo_pe.1.fq.gz -2 Geo_pe.2.fq.gz --s1 Geo_se.1.fq.gz --s2 Geo_se.2.fq.gz -o run1_spades_out --threads 48 --memory 300 -k 21,33,55,77,99,127

Then we use RNA-seq data to scaffold the draft genome assemblies with P_RNA_scaffolder (https://github.com/CAFS-bioinformatics/P_RNA_scaffolder): 
  
  hisat2 -x scaffolds.fasta -1 Geo_RNA_pe.1.fq.gz -2 Geo_RNA_pe.2.fq.gz -k 3 -p 10 --pen-noncansplice 1000000 -S input.sam

  bash /home/P_RNA_scaffolder.sh -d /stor/work/Ochman/yli19/projects/tom_aphid_genomes/bin/P_RNA_scaffolder/P_RNA_scaffolder/ -i input.sam -j discovar.fasta -F Geo_RNA_pe.1.fq.gz -R Geo_RNA_pe.2.fq.gz -o ./run1/
  
The resulting assemblies are named **P_RNA_scaffold.fasta** under the ./run1/ folder. 

### 1.3 - Contamination removal using BlobTools version 1.1.1

BLAST final assemblies against nt database
  # download NT database locally: 
  /home/BLAST/ncbi-blast-2.9.0+/bin/update_blastdb.pl --decompress --blastdb_version 5 nt
  /home/ncbi-blast-2.11.0+/bin/blastn -task megablast -query P_RNA_scaffold.mask.fasta -db /home/nt_database/nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -culling_limit 5 -num_threads 40 -evalue 1e-25 -out assembly_vs_nt.cul5.1e25.megablast.out



  
## 2 - Aphid genome annotation

## 3 - Aphid orthologous assignment

## 4 - Aphid phylogeny construction

## 5 - Buchnera genome assembly

## 6 - Buchnera genome annotation



# 


# Citation
Smith T.E., Li Y., Moran N.A. Comparative genomics of eight aphid subfamilies reveals variable relationships between host horizontally-transferred genes and symbiont peptidoglycan metabolism

