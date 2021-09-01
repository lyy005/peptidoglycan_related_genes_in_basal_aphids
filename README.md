# The distribution of peptidoglycan related genes in aphids and Buchnera 

This document is a walkthrough of the methods and code used to analyze the peptidoglycan related genes in aphid genomes and Bunchnera genomes.

## 1 - Aphid genome assembly

Raw data can be downloaded from NCBI BioProject (PRJNA758084). 

### 1.1 - Quality control: 

We cleaned both genomic data and RNA-seq data using Trimmomatic version 0.38. For example, we used the following commands for the Geopemphigus genome: 

        java -jar trimmomatic-0.38.jar PE -phred33 Geo_1.fastq.gz Geo_2.fastq.gz Geo_pe.1.fq.gz Geo_se.1.fq.gz Geo_pe.2.fq.gz Geo_se.2.fq.gz ILLUMINACLIP:/adapters/combined.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  

### 1.2 - Genome assembly

We first assembled draft genomes based on the genomic data using SPADES version 3.14.1: 

        spades.py -1 Geo_pe.1.fq.gz -2 Geo_pe.2.fq.gz --s1 Geo_se.1.fq.gz --s2 Geo_se.2.fq.gz -o run1_spades_out --threads 48 --memory 300 -k 21,33,55,77,99,127

Then we use RNA-seq data to scaffold the draft genome assemblies with P_RNA_scaffolder (https://github.com/CAFS-bioinformatics/P_RNA_scaffolder): 
  
        hisat2 -x scaffolds.fasta -1 Geo_RNA_pe.1.fq.gz -2 Geo_RNA_pe.2.fq.gz -k 3 -p 10 --pen-noncansplice 1000000 -S input.sam
        bash /home/P_RNA_scaffolder.sh -d /stor/work/Ochman/yli19/projects/tom_aphid_genomes/bin/P_RNA_scaffolder/P_RNA_scaffolder/ -i input.sam -j discovar.fasta -F Geo_RNA_pe.1.fq.gz -R Geo_RNA_pe.2.fq.gz -o ./run1/
  
The resulting assemblies are named **P_RNA_scaffold.fasta** under the ./run1/ folder. You can rename the file to Geo.genome.fasta using:
        mv P_RNA_scaffold.fasta Geo.genome.fasta

### 1.3 - Contamination removal

We explored potential contamination using BlobTools version 1.1.1 (https://blobtools.readme.io/docs/blobplot)

BLAST final assemblies against nt database
        
        # download NT database locally 
        /home/BLAST/ncbi-blast-2.9.0+/bin/update_blastdb.pl --decompress --blastdb_version 5 nt
        
        /home/ncbi-blast-2.11.0+/bin/blastn -task megablast -query Geo.genome.fasta -db nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -culling_limit 5 -num_threads 40 -evalue 1e-25 -out assembly_vs_nt.cul5.1e25.megablast.out

Estimating sequencing depth: 
        bwa mem -t 50 -o Geo_pe.sam Geo.genome.fasta Geo_pe.1.fq.gz Geo_pe.2.fq.gz

        samtools view -h -b -S Geo_pe.sam -o Geo_pe.bam --threads 20
        samtools sort Geo_pe.bam -o Geo_pe.sorted.bam --threads 20
        samtools index Geo_pe.sorted.bam

Create Blobplots:

        blobtools create -i Geo.genome.fasta -b Geo_pe.sorted.bam -t discovar_vs_nt.cul5.1e25.megablast.out -o blob_out
        blobtools view -i blob_out.blobDB.json -r all -o blob_taxonomy
        blobtools plot -i blob_out.blobDB.json -r genus
        blobtools plot -i blob_out.blobDB.json -r order
        
For *Chaitophorus*, an extra filtering was used to map the assembly to the Pemphigus assembly
        blastn -query Geo.genome.fasta -db Pem.genome.fasta -out Cha.scaf2genome.out -evalue 1e-10 -outfmt "6 qlen slen qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads 32

Assign scaffolds as Pemphigus contaminations (1) if the scaffold > 1000bp, has 1000 bp alignment and >= 95% identity (2) or if shorter than 1000bp, with >= 95% identity

        perl filter_scaffolds.alignment_only.pl Cha.scaf2genome.out ../blob_taxonomy.blob_out.blobDB.table.txt Cha.scaf2genome.out.list 
        cat Cha.scaf2genome.out.list.Blast Cha.scaf2genome.out.list.Blast.short | awk '{print $3}' | uniq > Cha.PemContam.list

Lastly, any 


  
## 2 - Aphid genome annotation

## 3 - Aphid orthologous assignment

## 4 - Aphid phylogeny construction

## 5 - Buchnera genome assembly

## 6 - Buchnera genome annotation



# 


# Citation
Smith T.E., Li Y., Moran N.A. Comparative genomics of eight aphid subfamilies reveals variable relationships between host horizontally-transferred genes and symbiont peptidoglycan metabolism

