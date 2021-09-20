# R script for identifying and counting genes per Buchnera genome
# from OrthoFinder output
# Written by Thomas E. Smith in June 2021
# Last updated by TES September 2021

#### NOTE
# This script assumes the User is in the same directory "that
# step7_bacteria_ortholog_assignment.sh ended in; bacteria_dfast_protein/

# Additionally, the User must place the following files 
# in the parent directory:
# "All_instances_of_Genes_in_Escherichia_coli_K-12_substr._MG1655.txt"
# and "aphid_PGN_count.csv"
# The former was originally downloaded from EcoCyc.org and contains
# current gene names for E. coli. 
# The latter is the manually created gene counts file from blastp and tblastn 
# searches of insect PGN genes and is available at 
# https://zenodo.org/record/5484415#.YTob355KgWp


## Setting the stage
# Load necessary packages 
library("tidyverse")
library("broom")
# Change into parent directory
#setwd("..")
# Import smart table from EcoCyc.org containing current MG1655 E. coli 
# gene names and synonyms
gene.names <- readLines("All_instances_of_Genes_in_Escherichia_coli_K-12_substr._MG1655.txt", 
                        warn = F)[-1] %>% 
  tibble(names.to.split = .) %>%
  separate(names.to.split, into = c("gene", "synonyms"), sep = "\t")
# Import results of aphid gene presence/absence data 
## based on blastp and tblastn searches.
df.aphid <- read.table("aphid_PGN_count.csv", header = T, sep = ",") %>%
  as_tibble()
# Establish working directories
# Change into the bacteria_orthofinder directory
setwd("bacteria_orthofinder")
# Locate the most recent OrthoFinder results folder
res.fold <- row.names(file.info(list.files()))[which.max(file.info(list.files())$mtime)]
wd <- paste0(getwd(), "/", res.fold, "/")
setwd(paste0(wd))


## Import OrthoFinder data
# Import both assigned and unassigned orthogroups from Orthogroups 
# directory and combine together
df.assigned <- read.table("Orthogroups/Orthogroups.tsv", header = T, sep = "\t")
df.unassigned <- read.table("Orthogroups/Orthogroups_UnassignedGenes.tsv", 
                            header = T, sep = "\t")
df <- rbind(df.assigned, df.unassigned)


## Convert data frame to tibble, reformat, and calculate gene copies per
## orthogroup per species
df.format <- df %>% 
  as_tibble() %>%
  # remove _filt tag of each species column header
  rename_with( ~ str_remove(.x, "_filt") ) %>%
  # create two new columns of gene names derived from either 
  # Esco_Ensmbl or BuapAcpi_Ensmbl columns
  rowwise() %>%
  mutate(Esco.gene.names = if_else(any(str_detect(Esco_Ensmbl, ",")),
                                   str_replace_all(Esco_Ensmbl, ", ", "\\."),
                                   Esco_Ensmbl) ) %>%
  mutate(Esco.gene.names = str_replace_all(Esco.gene.names, 
                                           "_[A-Z]+\\d+", "") ) %>%
  mutate(BuapAcpi.gene.names = if_else(any(str_detect(BuapAcpi_Ensmbl, ",")),
                                       str_replace_all(BuapAcpi_Ensmbl, ", ", "\\."),
                                       BuapAcpi_Ensmbl) ) %>%
  mutate(BuapAcpi.gene.names = str_replace_all(BuapAcpi.gene.names, 
                                               "_[A-Z]+\\d+", "") ) %>%
  ungroup() %>%
  # convert each cell to a count of the number of genes described per cell
  mutate(across(2:(length(.)-2), ~ str_count(.x, pattern = "\\w+(-\\w+){0,2}") )) %>%
  # replace NAs from any unassigned orthogroups with 0's
  mutate(across(everything(), ~replace(., is.na(.), 0) )) %>%
  # convert columns to numeric
  mutate(across(2:(length(.)-2), ~as.numeric(.)) )


## Ensure genes are uniformly named by identifying and changing any gene name
## synonyms to the current preferred EcoCyc.org gene name
df.names <- df.format
gene.columns <- c("Esco.gene.names", "BuapAcpi.gene.names")
for (h in 1:length(gene.columns)) {
  gene.column <- gene.columns[h]
  for (i in 1:nrow(df.names)) {
    genes.i <- df.names[i, gene.column]
    if (str_count(genes.i) > 0) {
      if (str_count(genes.i, "\\.") > 0) {
        genes.i <- unlist(str_split(genes.i, "\\."))
        new.genes <- c()
        for (j in 1:length(genes.i)) {
          genes.j <- genes.i[j]
          syn.check <- grepl(paste0(genes.j, "(?![a-zA-Z]{1}\\.*)"), gene.names$synonyms, perl = T)
          con.check <- grepl(paste0(genes.j, "(?![a-zA-Z]{1}\\.*)"), gene.names$gene, perl = T)
          # is this gene not a synonym?
          if (!(any(syn.check))) { 
            new.gene <- genes.j
            new.genes <- c(new.genes, new.gene) 
            next }
          # is this gene name a synonym?
          if (any(syn.check)) {
            # is this gene also a consensus gene name?
            if (any(con.check)) { 
              new.gene <- genes.j
              new.genes <- c(new.genes, new.gene) 
              next } 
            # this gene is only a synonym
            else { new.gene <- gene.names$gene[ syn.check ]
            new.genes <- c(new.genes, new.gene) }
          }
          # is this gene envC? envC is the synonymous gene name for acrE, 
          # so must specifically avoid renaming envC
          #if (grepl("envC", genes.j)) { new.gene <- "envC" }
          new.genes <- unique(new.genes)
          print(paste(paste(genes.i, collapse = "."), " in column ", gene.column, " row ", i, " to ", paste(new.genes, collapse = ".")))
          df.names[i, gene.column] <- paste(new.genes, collapse = ".")
        }
      } 
      else if (str_count(genes.i, "\\.") == 0) {
        syn.check <- grepl(paste0(genes.i, "(?![a-zA-Z]{1}\\.*)"), gene.names$synonyms, perl = T)
        con.check <- grepl(paste0(genes.i, "(?![a-zA-Z]{1}\\.*)"), gene.names$gene, perl = T)
        # is this gene not a synonym?
        if (!(any(syn.check))) { next }
        # is this gene name a synonym?
        if (any(syn.check)) {
          # is this gene also a consensus gene name?
          if (any(con.check)) { next } 
          # this gene is only a synonym
          else { new.gene <- gene.names$gene[ syn.check ] }
          # if only a three-letter code is given, grep might return multiple 
          # gene names. Defer to E. coli gene name in this instance
          if (length(new.gene) > 1) { new.gene <- df.names[ i, "Esco.gene.names" ] }
        }
        # is this gene envC?
        #if (grepl("envC", genes.i)) { new.gene <- "envC" }
        print(paste(genes.i, " in column ", gene.column, " row ", i, " to ", new.gene))
        df.names[i, gene.column] <- new.gene
      }
    }
  }
}


## Determine whether any genes between Esco and BuapAcpi Ensmbl 
## annotations fall into different orthogroups and merge these orthogroups
# Extract both Esco and BuapAcpi genes
Esco.genes <- c()
BuapAcpi.genes <- c()
for (i in 1:nrow(df.names)) {
  if ( str_count(df.names$Esco.gene.names[i], ".") > 0 ) {
    x <- df.names$Esco.gene.names[i]
    if ( grepl("\\.", x) ) {
      Esco.genes <- c(Esco.genes, unlist(strsplit(x, "\\.")))
    } else {
      Esco.genes <- c(Esco.genes, x)
    }
  }
  if ( str_count(df.names$BuapAcpi.gene.names[i], ".") > 0 ) {
    y <- df.names$BuapAcpi.gene.names[i]
    if ( grepl("\\.", y) ) {
      BuapAcpi.genes <- c(BuapAcpi.genes, unlist(strsplit(y, "\\.")))
    } else {
      BuapAcpi.genes <- c(BuapAcpi.genes, y)
    }
  }
}
# Combine Esco and BuapAcpi genes, remove duplicates and unnamed genes, and
# gather gene Orthogroup information
all.genes <- tibble(genes = c(Esco.genes, BuapAcpi.genes)) %>%
  distinct(genes) %>%
  filter( str_count(genes, ".") > 0 ) %>%
  # For each species, retrieve the orthogroup to which each gene 
  # belongs
  rowwise() %>%
  mutate( Esco.OG = if_else( length(df.names$Orthogroup[grep(paste0(genes, "(?![a-zA-Z]{1}\\.*)"), df.names$Esco.gene.names, perl = T)]) > 1,
                             df.names$Orthogroup[grep(paste0(genes, "(?![a-zA-Z]{1}\\.*)"), df.names$Esco.gene.names, perl = T)][match(genes, unlist(str_split(df.names$Esco.gene.names[grep(paste0(genes, "(?![a-zA-Z]{1}\\.*)"), df.names$Esco.gene.names, perl = T)], "\\.")) )],
                             df.names$Orthogroup[grep(paste0(genes, "(?![a-zA-Z]{1}\\.*)"), df.names$Esco.gene.names, perl = T)][1] ),
          BuapAcpi.OG = if_else( length(df.names$Orthogroup[grep(paste0(genes, "(?![a-zA-Z]{1}\\.*)"), df.names$BuapAcpi.gene.names, perl = T)]) > 1,
                                 df.names$Orthogroup[grep(paste0(genes, "(?![a-zA-Z]{1}\\.*)"), df.names$BuapAcpi.gene.names, perl = T)][match(genes, unlist(str_split(df.names$BuapAcpi.gene.names[grep(paste0(genes, "(?![a-zA-Z]{1}\\.*)"), df.names$BuapAcpi.gene.names, perl = T)], "\\.")) )],
                                 df.names$Orthogroup[grep(paste0(genes, "(?![a-zA-Z]{1}\\.*)"), df.names$BuapAcpi.gene.names, perl = T)][1] ) ) %>%
  ungroup()
# For which genes don't the Esco and BuapAcpi orthogroups match?
nonmatching.OGs <- all.genes %>%
  mutate( OG.match = if_else( Esco.OG == BuapAcpi.OG, T, F) ) %>%
  filter( OG.match == F ) %>%
  select( -OG.match ) %>%
  rowid_to_column( )
# Extract the orthofinder data for these genes...
merged.OGs <- df.names %>%
  rowwise() %>%
  mutate( merge.num = ifelse(any(grepl(Orthogroup, nonmatching.OGs)),
                             nonmatching.OGs$rowid[grep(Orthogroup, unlist(nonmatching.OGs[grepl(Orthogroup, nonmatching.OGs)]) )],
                             0) ) %>%
  ungroup() %>%
  filter( merge.num > 0 ) %>%
  # ...then merge all data columns
  group_by( merge.num ) %>%
  summarise( family = paste0(Orthogroup, collapse = "."),
             genes = paste0(Esco.gene.names, BuapAcpi.gene.names, collapse = "."),
             across( 2:(length(.)-3), ~ sum(.) )) %>%
  select( -merge.num )
# Remove these orthogroup rows from the full data then add the merged data
df.bacteria <- df.names %>%
  filter( !(Orthogroup %in% unlist(nonmatching.OGs[,3:4])) ) %>%
  rename( family = Orthogroup,
          genes = Esco.gene.names ) %>%
  select( -BuapAcpi.gene.names ) %>%
  rbind( merged.OGs ) %>%
  # remove ensembl columns, as they have fulfilled their purpose
  select(!ends_with("_Ensmbl")) %>%
  # reorder columns for readability
  relocate(c(family, genes, 2:(length(.)-1)))


## Extract each gene group of interest from df.bacteria
# Define genes of interest for each pathway/group
pgn.genes <- c("glmS", "glmM", "glmU", "ispU", "murA", "murB", "murC", "murI", 
               "murD","murE", "alr", "dadX", "metC", "glyA", "ddlA", "ddlB", 
               "murF", "mraY", "murG", "murJ", "bacA", "pgpB", "ybjG", # lipid II biosynthesis
               "mrcA", "mrcB", "mtgA", "mrdA", "ftsI", "mrdB", "ftsW", "dacA", 
               "dacC", "dacD", "ldtE", "ldtD","ldtA", "ldtB", "ldtC", "lpp",
               "lpoA", "lpoB", "pbpC",  # cell wall synthesis
               "amiA", "amiB", "amiC", "amiD", "slt", "mltA", "mltB", "mltC", 
               "mltD", "emtA", "mltF", "mltG", "rlpA", "mepA", "mepS", "mepM", 
               "mepH", "ampH", "mepK", "envC", "nlpD", "dacB", "yfeW", "pbpG", # cell wall remodeling
               "ampG", "mppA", "oppB", "oppC", "oppD", "oppF", "ampD", "nagZ", 
               "ldcA", "mpl", "anmK", "murP", "murQ", "nagK", "nagA", "nagB") # PGN recycling
div.genes <- c("minC", "minD", "minE", #minCDE
               "ftsZ", "ftsA", "zipA", "ftsE", "ftsX", "ftsQ", "ftsL", "ftsB", 
               "ftsN", "ftsK", "ftsP", "zapA", "zapC", "zapD", "zapB", "sulA", # divisome
               "mreB", "mreC", "mreD", "rodZ") # elongasome
fbb.genes <- c("fliF", "flgI", "flgH", # basal body
               "fliE", "flgB", "flgC", "flgF", "flgG", # rod
               "fliK", "flgD", "flgE", "flgK", "flgL", # hook
               "fliC", "fliD", # filament
               "flhA", "flhB", "fliO", "fliP", "fliQ", "fliR", "fliH", "fliI",
               "fliJ", "flgJ", # export apparatus
               "fliG", "fliM", "fliN", # rotor
               "motA", "motB", # stator
               "flgN", "flgA", # chaperones
               "cheY", # response regulator
               "fliA") # sigma factor
# 10 genes of df.pgn shared between pgn and div pathways: mrcA, lpoA, mrcB, lpoB, amiB, amiC, rlpA, ftsW, ftsI, mrdB
# 1 gene of df.fbb shared between pgn and fbb pathways: flgJ
# 0 genes shared between div and fbb pathways

# For peptidoglycan genes...
df.bacteria.pgn <- df.bacteria %>%
  left_join(x = tibble(genes = pgn.genes), by = "genes") %>%
  relocate(c(2, 1, 3:length(.)))
# Some genes of interest are part of orthogroups containing multiple 
# genes, such that the match function of left_join doesn't find them.
# Extract genes for which data was not successfully imported through join
multigeneOG.genes <- df.bacteria.pgn$genes[is.na(df.bacteria.pgn$family)]
# Extract rows of df.bacteria containing these genes
df.multigeneOGs.genes <- sapply(multigeneOG.genes, function(x) 
{ df.bacteria[grepl(x, df.bacteria$genes), ] }) %>%
  unlist() %>%
  matrix(nrow = length(.)/ncol(df.bacteria.pgn),
         ncol = length(.)/length(multigeneOG.genes),
         byrow = T) %>%
  as_tibble(.name_repair = "minimal")
names(df.multigeneOGs.genes) <- names(df.bacteria.pgn)
# Delete the empty rows from df.bacteria.pgn, then insert the rows from
# df.multigeneOGs.genes and place them in the proper order
df.bacteria.pgn <- df.bacteria.pgn %>%
  filter(!is.na(family)) %>%
  rbind(df.multigeneOGs.genes) %>%
  # Remove duplicate rows to eliminate two equivalent rows for dacA/dacC
  distinct(family, .keep_all = T)

# For cell division genes, follow the same procedure as for pgn genes
df.bacteria.div <- df.bacteria %>%
  left_join(x = tibble(genes = div.genes), by = "genes") %>%
  relocate(c(2, 1, 3:length(.)))
multigeneOG.genes <- df.bacteria.div$genes[is.na(df.bacteria.div$family)]
df.multigeneOGs.genes <- sapply(multigeneOG.genes, function(x) 
{ df.bacteria[grepl(x, df.bacteria$genes), ] }) %>%
  unlist() %>%
  matrix(nrow = length(.)/ncol(df.bacteria.div), 
         ncol = length(.)/length(multigeneOG.genes), 
         byrow = T) %>%
  as_tibble(.name_repair = "minimal")
names(df.multigeneOGs.genes) <- names(df.bacteria.div)
# Swap empty rows for filled ones
df.bacteria.div <- df.bacteria.div %>%
  filter(!is.na(family)) %>%
  rbind(df.multigeneOGs.genes) %>%
  # Remove duplicate rows to eliminate any equivalent rows
  distinct(family, .keep_all = T)

# For FBB genes, follow the same procedure as for pgn genes
df.bacteria.fbb <- df.bacteria %>%
  left_join(x = tibble(genes = fbb.genes), by = "genes") %>%
  relocate(c(2, 1, 3:length(.)))
multigeneOG.genes <- df.bacteria.fbb$genes[is.na(df.bacteria.fbb$family)]
df.multigeneOGs.genes <- sapply(multigeneOG.genes, function(x) 
{ df.bacteria[grepl(x, df.bacteria$genes), ] }) %>%
  unlist() %>%
  matrix(nrow = length(.)/ncol(df.bacteria.fbb), 
         ncol = length(.)/length(multigeneOG.genes), 
         byrow = T) %>%
  as_tibble(.name_repair = "minimal")
names(df.multigeneOGs.genes) <- names(df.bacteria.fbb)
# Swap empty rows for filled ones
df.bacteria.fbb <- df.bacteria.fbb %>%
  filter(!is.na(family)) %>%
  rbind(df.multigeneOGs.genes) %>%
  # Remove duplicate rows to eliminate any equivalent rows
  distinct(family, .keep_all = T)

# Combine the three data subsets together
df.bacteria.gois <- rbind(df.bacteria.pgn, 
                          df.bacteria.div, 
                          df.bacteria.fbb)

# Extract species for which >2 copies of any gene are present, 
# then manually investigate whether this is true and edit the 
# data frame at these specific coordinates
inspect <- which(df.bacteria.gois > 1, arr.ind = T)
# omit genes and family columns
inspect <- inspect[inspect[,2] > 2, ]
df.inspect <- tibble(gene = df.bacteria.gois$genes[inspect[,1]],
                     species = names(df.bacteria.gois)[inspect[,2]]) %>%
  # Omit all columns that aren't Buchnera or Skilesia
  filter(grepl("Buap", species)|grepl("Skal", species)) %>%
  rowwise() %>%
  mutate(orthogroup = df.bacteria.gois$family[str_detect(gene, df.bacteria.gois$genes)] )
# Output manual corrections to be made
print("BuapApna MGA_328 and MGA_329 are two pieces of a pseudogenized flgN gene. Set 0 for flgN (OG0000565). E. coli gene is OG0002125.")
print("BuapErla MGA_228 through MGA_231 are pieces of a pseudogenized fliF gene. Set to 0 copies for fliF (OG0000341).")
print("Both BuapErla MGA_205 and MGA_547 BLAST to Buchnera fliQ and are located on small contigs (< 10 kb). They are identical in sequence, so set to 1 copy for fliQ (OG0000366).")
print("BuapErla MGA_235 and MGA_236 are two parts of a pseudogenized fliH gene. Set to 0 for fliH (OG0000405).")
print("BuapErla MGA_208 and MGA_549 BLAST to fliN from Pseudidomarina homiensis and Pusillimonas harenae, but neighboring genes on the same contig BLAST to Buchnera homologs. They are identical in  sequence, so set to 1 copy for fliN (OG0000590). OG0000571 is also fliN.")
print("BuapHoco MGA-361 and MGA-362 are the N- and C-terminal regions of a pseudogenized flgB gene. Set to 0 for flgB (OG0000447).")
print("BuapHoco MGA_473 and MGA_475 correspond to the C- and N-terminal regions of flhA, an are separated by a sequence gap. Likely a complete gene judging from longer sequence length of each part. Set to 1 for flhA (OG0000344).")
print("SkalGesp MGA_380 is ftsW, while MGA_1063-LOCUS_10590 is mrdB (rodA), which all Buchnera lack. Part of elongasome. Set to 1 copy for both mrdB (OG0000706) and ftsW (OG0000373).")
# Make manual corrections
df.bacteria.gois <- df.bacteria.gois %>%
  mutate(BuapApna = if_else(str_detect(genes, "flgN"), "1", BuapApna),
         BuapErla = if_else(str_detect(genes, "fliF"), "0", BuapErla),
         BuapErla = if_else(str_detect(genes, "fliQ"), "1", BuapErla),
         BuapErla = if_else(str_detect(genes, "fliH"), "0", BuapErla),
         BuapErla = if_else(str_detect(genes, "fliN"), "1", BuapErla),
         BuapHoco = if_else(str_detect(genes, "flgB"), "0", BuapHoco),
         BuapHoco = if_else(str_detect(genes, "flhA"), "1", BuapHoco),
         SkalGesp = if_else(str_detect(genes, "ftsW"), "1", SkalGesp),
         SkalGesp = if_else(str_detect(genes, "mrdB"), "1", SkalGesp) )


## Format aphid gene presence/absence data similarly to bacterial data
# Dici amiD gene was acquired from a different source than aphid 
# amiD gene, so manually change it to 0
df.aphid <- df.aphid %>%
  mutate(Dici = ifelse(family == "AmiD", 0, Dici))
# Rename columns and remove outgroup columns
df.aphid <- df.aphid %>%
  rename(genes = family) %>%
  mutate(family = paste0("INSECT", 1:11)) %>%
  relocate(25, 1:24)
df.aphid.filt <- df.aphid %>%
  select(-c(20:25))  %>%
  # Reorder remaining columns by order in phylogenetic tree
  relocate(c(1:14, Cice, Erla, Gesp, Pesp, Hoco))


## Select only species for which host and aphid genomes are available
host.species <- colnames(df.aphid.filt)[3:19]
df.bacteria.gois.filt <- df.bacteria.gois %>%
  select(c(family, genes, contains(host.species)) ) %>%
  select(-contains("Sesy"))
names(df.bacteria.gois.filt)[3:19] <- host.species
# Combine bacterial and aphid data together
df.gois <- rbind(df.aphid.filt, df.bacteria.gois.filt) 
# Remove genes that no species contain
df.gois.filt <- df.gois %>%
  # Convert data to numeric
  mutate(across(3:19, ~as.numeric(.)) ) %>%
  rowwise() %>%
  mutate(row.sum = sum(c_across(3:19)) ) %>%
  filter(row.sum > 0) %>%
  select(-row.sum)
# Add pathway information to facilitate data parsing
df.gois.filt <- df.gois.filt %>%
  mutate(role = case_when(match(family, df.gois.filt$family) %in% 1:9 ~ "host",
                          match(family, df.gois.filt$family) %in% 10:27 ~ "lipid II biosynthesis",
                          match(family, df.gois.filt$family) %in% 28:32 ~ "cell wall synthesis",
                          match(family, df.gois.filt$family) %in% 33:38 ~ "cell wall remodeling",
                          match(family, df.gois.filt$family) %in% 39 ~ "cell wall synthesis",
                          match(family, df.gois.filt$family) %in% 40 ~ "cell wall remodeling",
                          match(family, df.gois.filt$family) %in% 41:51 ~ "cell division",
                          match(family, df.gois.filt$family) %in% 52:nrow(.) ~ "flagellar basal body") ) %>%
  relocate(c(1:2, 20, 3:19))
# Define each pathway of interest in a vector
pathways <- c("host", "lipid II biosynthesis", "cell wall synthesis",
              "cell wall remodeling", "peptidoglycan recycling",
              "cell division", "flagellar basal body")
# Assign colors to gene absence and presence for each pathway
pathway.colors <- data.frame(role = pathways,
                             absence = c("honeydew",
                                         "cornsilk",
                                         "lemonchiffon",
                                         "mistyrose",
                                         "ghostwhite",
                                         "lavenderblush",
                                         "aliceblue"),
                             presence = c("palegreen2",
                                          "goldenrod2",
                                          "gold2",
                                          "palevioletred2",
                                          "gray60",
                                          "orchid2",
                                          "steelblue2"))


## Construct heatmaps depicting the peptidoglycan gene 
## presence/absence in each species, including both host and 
## symbiont genes
library("pheatmap")
for (i in 1:length(pathways)) {
  pathway <- pathways[i]
  if (any(pathway %in% df.gois.filt$role)) {
    
    # Format the data matrix
    matrix.gois <- df.gois.filt %>%
      filter(role == pathway)
    # store the gene names in a vector before removing them
    gois <- matrix.gois$genes
    # convert data frame to matrix
    matrix.gois <- matrix.gois %>%
      select(-family, -genes, -role) %>%
      as.matrix()
    
    # Set binary colors for the potentially non-binary range 
    # of data values
    matrix.colors <- c(pathway.colors$absence[grep(pathway, pathway.colors$role)], 
                       rep(pathway.colors$presence[grep(pathway, pathway.colors$role)],
                           range(matrix.gois)[2]) )
    
    # Establish whether gene copy numbers should be shown within
    # each cell or not
    show.numbers <- case_when(pathway == "host" ~ T,
                              pathway != "host" ~ F)
    # Cells with numbers should be slightly larger then those without
    cell.height <- case_when(pathway == "host" ~ 8,
                             pathway != "host" ~ 6)
    # Cluster genes by presence/absence pattern except for host genes
    cluster.rows <- case_when(pathway == "host" ~ F,
                              pathway != "host" ~ T)
    
    # Create and export the heatmap
    file.name <- paste0("heatmap_", pathway, ".pdf")
    pdf(file.name)
    pheatmap(matrix.gois, 
             color = matrix.colors,
             legend = F,
             labels_row = gois, labels_col = host.species,
             cluster_rows = cluster.rows, cluster_cols = F,
             treeheight_row = 0,
             fontsize = 7,
             cellheight = cell.height, cellwidth = 12,
             show_rownames = T, 
             display_numbers = show.numbers, fontsize_number = 5, 
             number_color = "black", number_format = "%.0f",
             border_color = "goldenrod4")
    dev.off()
  }
}