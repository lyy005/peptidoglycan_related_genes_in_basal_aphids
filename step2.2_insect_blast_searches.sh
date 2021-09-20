# Script to identify aphid PGN genes by BLAST
# Written by TES in May, 2021
# Last updated by TES in Sep, 2021

## NOTES
# Following genome assembly and annotation, insect PGN genes were identified by BLAST 
# search using known PGN-associated genes as queries. User must first download the 
# following files into separate directories: 
# 1) Insect genome annotations (Table S1) into ~/insect_annotations/
# 2) Query protein sequences (Table S7; separate fasta file per query) into ~/queries/
# 3) Insect genome assemblies (Table S1) into ~/insect_genomes/ 

# Next, utilize this script within the parent directory using the following command:
# $ bash step6_bacteria_genome_annotations.sh


## Format annotation files
# Append header name to every sequence so that
# each line has a unique identifier after concatenation
cd ~/insect_annotations/
for file in *
do
# Assign a custom species identifier based on file name
if [[ $file == Aphis_glycines_4.v2.1.scaffolds.fa.gff.pep.fa ]]
then ID=`echo "Apgl"`
elif [[ $file == GCF_004010815.1_ASM401081v1_translated_cds.faa ]]
then ID=`echo "Apgo"`
elif [[ $file == Pentalonia_nigronervosa.v1.scaffolds.gff.pep.fa ]]
then ID=`echo "Peni"`
elif [[ $file == Rpadi_v2_pep.fa ]]
then ID=`echo "Rhpa"`
elif [[ $file == GCF_003676215.2_ASM367621v3_translated_cds.faa ]]
then ID=`echo "Rhma"`
elif [[ $file == GCF_005508785.1_pea_aphid_22Mar2018_4r6ur_translated_cds.faa ]]
then ID=`echo "Acpi"`
elif [[ $file == GCF_001186385.1_Dnoxia_1.0_translated_cds.faa ]]
then ID=`echo "Dino"`
elif [[ $file == GCF_001856785.1_MPER_G0061.0_translated_cds.faa ]]
then ID=`echo "Mype"`
elif [[ $file == M.cerasi.V1.protein.fasta ]]
then ID=`echo "Myce"`
elif [[ $file == Chr_genome_final_gene.gff3.pep.fa ]]
then ID=`echo "Simi"`
elif [[ $file == Hormaphis_cornu.pep.fasta ]]
then ID=`echo "Hoco"`
elif [[ $file == cinced3A.pep.fa ]]
then ID=`echo "Cice"`
elif [[ $file == Eriosoma_lanigerum_v1.0.scaffolds.gff.pep.fa ]]
then ID=`echo "Erla"`
elif [[ $file == Geopemphigus.aa ]]
then ID=`echo "Gesp"`
elif [[ $file == Stegophylla.aa ]]
then ID=`echo "Stsp"`
elif [[ $file == Chaitophorus.aa ]]
then ID=`echo "Chvi"`
elif [[ $file == Pemphigus.aa ]]
then ID=`echo "Peob"`
elif [[ $file == OGS3.1_20170928_protein.fasta ]]
then ID=`echo "Davi"`
elif [[ $file == coding_gene_translated_sequences.pep ]]
then ID=`echo "Erpe"`
elif [[ $file == HiC_genome_protein.fasta ]]
then ID=`echo "Phso"`
elif [[ $file == GCF_001854935.1_ASM185493v1_translated_cds.faa ]]
then ID=`echo "Beta"`
elif [[ $file == GCF_000475195.1_Diaci_psyllid_genome_assembly_version_1.1_translated_cds.faa ]]
then ID=`echo "Dici"`
elif [[ $file == augustus.hints.fa ]]
then ID=`echo "Pave"`
fi
# Make sure there are no periods or asterisks at the end of each sequence
sed -e "/^>/! s/\.//g" -e "/^>/! s/\*//g" $file > ${ID}.clean.aa.fasta
# Insert identifier at the start of each sequence header,
# then remove everything following the first space
sed "s/>/>${ID}_/g" ${ID}.clean.aa.fasta | cut -f 1 -d " " > ${ID}.header.aa.fasta
done

# Concatentate all insect annotations together
cat *.header.aa.fasta > insect_annotations.fasta 

# Finally, make the blast database
makeblastdb -dbtype prot -in insect_annotations.fasta -parse_seqids


## Next, perform blastp searches using the insect_annotation blastdb
# For each query sequence, perform a blastp search and save the output
cd ~/queries/
mkdir insect_blastp
for file in *.fasta
do
query=$(echo $file | cut -f 1 -d ".")
blastp -query $file -db ~/insect_annotations/insect_annotations.fasta \
-out ~/queries/insect_blastp/blastp_${query}.txt -evalue 1e-10 \
-outfmt "6 sseqid sstart send sseq"
done

# Produce a linearized insect annotation file such that each sequence is on a single line,
# making it easy to use grep -A 1 with a query that matches the sequence headers
cd ~/insect_annotations/
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' \
insect_annotations.fasta > insect_annotations_linearized.fasta

# For each resulting sequence file, retrieve full protein sequences from insect annotations
cd ~/queries/insect_blastp
for file in *.txt
do
query=$(echo $file | cut -f 1 -d "." | cut -f 3- -d "_")
# First, extract the resulting protein accession numbers
proteins=$(cut -f 1 -d $'\t' $file)
# Then, search insect annotations for these proteins and return sequences plus headers
# and save them as a fasta file
grep -A 1 --no-group-separator "$proteins" \
~/insect_annotations/insect_annotations_linearized.fasta > ${query}.fasta
done

# Sequence alignments were performed in SeaView using muscle
# Manual adjustment was necessary for some alignments, especially RlpA proteins.
# The manually aligned RlpA alignment is available at https://zenodo.org/record/5484415#.YTob355KgWp


## Format genome files
# Append header name to every sequence so that
# each line has a unique identifier after concatenation
cd ~/insect_genomes/
for file in *
do
# Assign a custom species identifier based on file name
if [[ $file == Aphis_glycines_4.v2.1.scaffolds.fa ]]
then ID=`echo "Apgl"`
elif [[ $file == GCF_004010815.1.fasta ]]
then ID=`echo "Apgo"`
elif [[ $file == GCA_014851325.1.fasta ]]
then ID=`echo "Peni"`
elif [[ $file == R_padi_v2.fasta ]]
then ID=`echo "Rhpa"`
elif [[ $file == GCF_003676215.2.fasta ]]
then ID=`echo "Rhma"`
elif [[ $file == GCF_005508785.1.fasta ]]
then ID=`echo "Acpi"`
elif [[ $file == GCF_001186385.1_Dnoxia_1.0_genomic.fna ]]
then ID=`echo "Dino"`
elif [[ $file == GCF_001856785.1.fasta ]]
then ID=`echo "Mype"`
elif [[ $file == Myzus.cerasi_genome.v1.1.fasta ]]
then ID=`echo "Myce"`
elif [[ $file == GCA_008086715.1.fasta ]]
then ID=`echo "Simi"`
elif [[ $file == hormaphis_cornu_26Sep2017_PoQx8.fasta ]]
then ID=`echo "Hoco"`
elif [[ $file == cinced3.scaffolds.fa ]]
then ID=`echo "Cice"`
elif [[ $file == GCA_013282895.1.fasta ]]
then ID=`echo "Erla"`
elif [[ $file == Geopemphigus.decontam.mask.fasta ]]
then ID=`echo "Gesp"`
elif [[ $file == Stegophylla.decontam.mask.fasta ]]
then ID=`echo "Stsp"`
elif [[ $file == Chaitophorus_decon2.decontam.mask.fasta ]]
then ID=`echo "Chvi"`
elif [[ $file == Pemphigus.decontam.mask.fasta ]]
then ID=`echo "Peob"`
elif [[ $file == Dv_genome_V3.1.fa ]]
then ID=`echo "Davi"`
elif [[ $file == assembly.fa ]]
then ID=`echo "Erpe"`
elif [[ $file == HiC_genome.fasta ]]
then ID=`echo "Phso"`
elif [[ $file == GCF_001854935.1_ASM185493v1_genomic.fna ]]
then ID=`echo "Beta"`
elif [[ $file == GCF_000475195.1_Diaci_psyllid_genome_assembly_version_1.1_genomic.fna ]]
then ID=`echo "Dici"`
elif [[ $file == GCA_012654025.1_Pven_dovetail_genomic.fna ]]
then ID=`echo "Pave"`
fi
# Insert identifier at the start of each sequence header,
# then remove everything following the first space
sed "s|>|>${ID}_|g" $file | cut -f 1 -d " " > ${ID}.header.fasta
done

# Concatentate all insect genomes together
cat *.header.fasta > insect_genomes.fasta  

# Finally, make the blast database
makeblastdb -dbtype nucl -in insect_genomes.fasta -parse_seqids


## Next, perform tblastn searches using the insect_genomes blastdb
# For the bLys query, perform a tblastn search and save the output
cd ~/queries/
mkdir insect_tblastn
for file in Acpi_bLys_lysozyme.fasta
do
query=$(echo $file | rev | cut -f 1 -d "/" | rev | cut -f 1 -d ".")
tblastn -query $file -db ~/insect_genomes/insect_genomes.fasta -out \
~/queries/insect_tblastn/tblastn_${query}.txt \
-evalue 1e-10 -outfmt "6 sseqid sstart send sseq"
done

# For each resulting sequence file, convert format to fasta
for file in *.txt
do
query=$(echo $file | cut -f 1 -d "." | cut -f 3- -d "_")
# Add > to start of each line and replace tab delimiters
sed -e 's/^/>/g' \
-e 's/\t/_/' \
-e 's/\t/-/' \
-e 's/\t/\n/' $file > ${query}.fasta
done

