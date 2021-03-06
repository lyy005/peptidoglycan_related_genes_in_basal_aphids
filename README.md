# The distribution of peptidoglycan related genes in aphids and Buchnera 

This document is a walkthrough of the methods and code used to analyze the peptidoglycan related genes in aphid genomes and Bunchnera genomes.

## 1 - Aphid genome assembly

Raw data can be downloaded from NCBI BioProject (PRJNA758084). 

### 1.1 - Quality control: 

We cleaned both genomic data and RNA-seq data using Trimmomatic version 0.38. For example, we used the following commands for the Geopemphigus genome: 

        # For genomic data
        java -jar trimmomatic-0.38.jar PE -phred33 Geo_1.fastq.gz Geo_2.fastq.gz Geo_pe.1.fq.gz Geo_se.1.fq.gz Geo_pe.2.fq.gz Geo_se.2.fq.gz ILLUMINACLIP:/adapters/combined.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        # For RNA-seq data
        java -jar trimmomatic-0.38.jar PE -phred33 Geo_RNA_1.fastq.gz Geo_RNA_2.fastq.gz Geo_RNA_pe.1.fq.gz Geo_RNA_se.1.fq.gz Geo_RNA_pe.2.fq.gz Geo_RNA_se.2.fq.gz ILLUMINACLIP:/home/Trimmomatic-0.38/adapters/combined.fa.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  

### 1.2 - Genome assembly

We first assembled draft genomes based on the genomic data using SPADES version 3.14.1: 

        spades.py -1 Geo_pe.1.fq.gz -2 Geo_pe.2.fq.gz --s1 Geo_se.1.fq.gz --s2 Geo_se.2.fq.gz -o run1_spades_out --threads 48 --memory 300 -k 21,33,55,77,99,127

The resulting assemblies are named as **scaffolds.fasta**.

Then we use RNA-seq data (mapped with HISAT2 version 2.2.1) to scaffold the draft genome assemblies with P_RNA_scaffolder (https://github.com/CAFS-bioinformatics/P_RNA_scaffolder): 
  
        hisat2 -x scaffolds.fasta -1 Geo_RNA_pe.1.fq.gz -2 Geo_RNA_pe.2.fq.gz -k 3 -p 10 --pen-noncansplice 1000000 -S input.sam
        bash /home/P_RNA_scaffolder/P_RNA_scaffolder.sh -d /home/P_RNA_scaffolder/ -i input.sam -j scaffolds.fasta -F Geo_RNA_pe.1.fq.gz -R Geo_RNA_pe.2.fq.gz -o ./run1/
  
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

Assign scaffolds as *Pemphigus* contaminations (1) if the scaffold > 1000bp, has 1000 bp alignment and >= 95% identity (2) or if shorter than 1000bp, with >= 95% identity

        perl filter_scaffolds.alignment_only.pl Cha.scaf2genome.out ../blob_taxonomy.blob_out.blobDB.table.txt Cha.scaf2genome.out.list 
        cat Cha.scaf2genome.out.list.Blast Cha.scaf2genome.out.list.Blast.short | awk '{print $3}' | uniq > Cha.PemContam.list

Lastly, any scaffolds that are identified as contamination will be removed: 
        
        # Extract scaffold IDs with primate or bacteria origin
        grep 'Primate\|Homo' blob_taxonomy.blob_out.blobDB.table.txt > Primate.blob.list
        grep "Bacteria" blob_taxonomy.blob_out.blobDB.table.txt > Bacteria.blob.list
        perl pick_sequences_NOT_on_list.pl Geo.genome.fasta Bacteria.list Chaitophorus_decon1.decontam.mask.fasta
        
        # For Chaitophorus, contaminated scaffolds from Pemphigus are also removed
        perl pick_sequences_NOT_on_list.pl Cha.genome.fasta Cha.PemContam.list Chaitophorus_decon1.decontam.mask.fasta

## 2 - Aphid genome annotation

### 2.1 - General gene annotation

Repeat masking using Tantan implemented in Funannotate https://github.com/nextgenusfs/funannotate

        funannotate-docker mask -i P_RNA_scaffold.sort.fasta -o P_RNA_scaffold.mask.fasta
        
Gene annotation using BRAKER version 2.1.5 https://github.com/Gaius-Augustus/BRAKER

        time perl /home/BRAKER/scripts/braker.pl --genome P_RNA_scaffold.mask.fasta --bam combined.sorted.bam --softmasking --cores 48 --workingdir run1 --species cha_spadesrna_prot --GENEMARK_PATH=/home/GeneMark/gmes_linux_64/ --prot_seq=protein_evidence_for_annotation.fasta --etpmode --PROTHINT_PATH=/home/ProtHint/ProtHint/bin/

### 2.2 - Peptidoglycan (PGN) related gene annotation

Following genome assembly and annotation, insect PGN genes were identified by BLAST search using known PGN-associated genes as queries. User must first download the following files into separate directories: 

1) Insect genome annotations (Table S1) into ~/insect_annotations/

2) Query protein sequences (Table S7; separate fasta file per query) into ~/queries/

3) Insect genome assemblies (Table S1) into ~/insect_genomes/ 

Next, utilize this script within the parent directory using the following command:

        bash step2.2_insect_blast_searches.sh

### 2.3 - Multiple sequence alignment

Protein alignments performed with MUSCLE in SeaView version 5.0.4, then, for RlpA sequences, manually adjusted. SeaView is available for download at http://doua.prabi.fr/software/seaview

### 2.4 - Alignment trimming

Alignments were trimmed using trimAl such that only 20% of sequences may contain a gap at each position. TrimAl is available at http://trimal.cgenomics.org

        ~/trimal-trimAl/source/trimal -in ~/alignment.fasta -out alignment_trimal.fasta -htmlout alignment_trimal.html -gt 0.8

### 2.5 - Gene tree construction

Protein trees were generated from trimmed alignments using IQ-TREE 2, available at http://www.iqtree.org, by first identifying the best-fit substitution model using ModelFinder, then using Transfer Bootstrap Estimates with the best-fit model (VTR5 in our case) to determine the final tree.

        iqtree -s alignment_trimal.fasta -m MFP -T auto -pre alignment_trimal_MFP
        iqtree -s alignment_trimal.fasta -m VT+R5 -pre alignment_trimal_VTR5-TBE -T auto -b 100 --tbe

### 2.6 - Test for positive selection

Following step 2.2, tests for positive selection were performed for aphid horizontally-transferred genes.

This step utilizes the output of the step 2.2, which includes a fasta file for each gene family containing the protein sequences of the resulting BLAST hits (i.e. ldcA.fasta). Nucleotide sequences for each protein sequence were retrieved manually by searching the respective insect coding sequence annotation file for each protein accession number (specified in the sequence header). Coding sequence files are available from the genome sources described in Table S1. 

Next, nucleotide sequences were codon-aligned to match amino acid sequence alignments using PAL2NAL, available at http://www.bork.embl.de/pal2nal. Finally, tests for positive selection were performed using BUSTED operated through Datamonkey, available at 
https://www.datamonkey.org/busted.

## 3 - Aphid ortholog assignment

Assign genes to orthologous groups using OrthoFinder version 2.5.2 (https://github.com/davidemms/OrthoFinder), including 23 genomes (available on Zenodo: https://zenodo.org/record/5517159). 
        
        orthofinder -t 36 -a 36 -f amino_acids_from_23_genomes/


## 4 - Aphid phylogeny construction

Reconstruct the phylogenetic tree with 1000 ultrafast bootstraps using IQ-Tree version 2.0.6. Bemisia tabaci (Beta), Diaphorina citri (Dici), Pachypsylla venusta (Pave) were used as outgroups. 

        /home/IQTree/iqtree-2.0.6-Linux/bin/iqtree2 -s Orthogroups.GeneCount.tsv.0.80.concate.fasta -o Beta,Dici,Pave -T 50 -B 1000
        

## 5 - Sitobion miscanthi Buchnera genome assembly

To assemble the Buchnera genome from Sitobion miscanthi, we downloaded the Illumina reads from NCBI using fastq-dump in NCBI SRA-tools (https://github.com/ncbi/sra-tools). Reads were from a previous Sitobion miscanthi genome sequencing project: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA532495

        fastq-dump --split-files SRR8988435
        
Genome assembly using SPADES

        java -jar /home/trimmomatic-0.38.jar PE -phred33 Sit.1.fq.gz Sit.2.fq.gz Sit_pe.1.fq.gz Sit_se.1.fq.gz Sit_pe.2.fq.gz Sit_se.2.fq.gz ILLUMINACLIP:/home/Trimmomatic-0.38/adapters/combined.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -threads 72

        /home/SPAdes-3.15.2-Linux/bin/spades.py -1 Sit_pe.1.fq.gz -2 Sit_pe.2.fq.gz --s1 Sit_se.1.fq.gz --s2 Sit_se.2.fq.gz -o run1_spades --threads 100 --memory 900 -k 21,33,55,77,99,127
        
## 6 - Buchnera genome annotation

User must first download the bacterial genomes listed in Table S9, then add unique headers to each file name (i.e. change "GCF_000009605.1_ASM960v1_genomic.fna" to "BuapAcpi_GCF_000009605.1_ASM960v1_genomic.fna").

User must also install dfast by following the instructions at: https://github.com/nigyta/dfast_core#installation

Next, utilize this script within the same directory using the following command:

        bash step6_bacteria_genome_annotations.sh

The output are the same annotated bacterial genomes available at https://zenodo.org/record/5484415

## 7 - Buchnera ortholog assignment

User must install orthofinder by following the instructions at: https://github.com/davidemms/OrthoFinder

Next, utilize this script within the same directory using the following command: 
        bash step7_bacteria_ortholog_assignment.sh

The output are the input for Step 8.

## 8 - Functional assignment of Buchnera orthologs

User must install R and/or RStudio along with the following R packages: tidyverse (v1.3.1), broom (v0.7.8), pheatmap (v1.0.12)

Additionally, the User must include the following files in the parent directory: "All_instances_of_Genes_in_Escherichia_coli_K-12_substr._MG1655.txt" and "aphid_PGN_count.csv". The former was originally downloaded from EcoCyc.org and contains current gene names for E. coli. The latter is the manually created gene counts file from blastp and tblastn searches of insect PGN genes. Both are available at https://zenodo.org/record/5484415

Next, utilize this script within the parent directory using the following command: 

        Rscript step8_host_symbiont_gene_count.R

The output are the heatmaps used to construct Figure 3.

# Citation
Smith T.E., Li Y., Moran N.A. Comparative genomics of eight aphid subfamilies reveals variable relationships between host horizontally-transferred genes and symbiont peptidoglycan metabolism

