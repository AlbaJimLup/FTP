## *************************** ##
## Load Libraries
## *************************** ##

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(readr)
library(tidyr)
library(dplyr)
library(patchwork)
library(cowplot)
library(reshape2)


## *************************** ##
## Input Files                 ##
## *************************** ##

# Read command line arguments
# Rscript AF_sample_fasta_plot.R ${out_dir}/merged_mutation_catalogue.simplified.out    ${out_dir}/Copy_Number.Matrix.out   ${plot_out_dir}/${name}.pdf
# variants_AF_Counts dataframe will have AF of each allele in all samples
#args <- commandArgs(trailingOnly = TRUE)
#mutation_data <- read.table(args[1],header=TRUE)
#CN.data.matrix <- read.table(args[2],header=TRUE) 
#output_file <- args[3]

inpath <- "~/Desktop/ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/T2T/01_simulated_reads/Miguel_sim_out/var/jobs/SyntheticRDNA/00_data_for_analysis/"
mutation_data <-  read.table(paste0(inpath, "merged_mutation_catalogue.simplified.out"), header=TRUE)
CN.data.matrix <- read.table(paste0(inpath, "Copy_Number.Matrix.out"),header=TRUE)
blacklisted <- read.table("~/Desktop/ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/T2T/scripts/data/47S_pre-rRNA.blacklisted_repetitive_seq.bed")
regions<-read.table("~/Desktop//ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/T2T/scripts/data/pre-rRNA_47S.regions.bed")

catalogue_size <- 50
number_of_samples<-100

## ******************************* ##
## Functions                       ##
## ******************************* ##

get_box_stats <- function(y, upper_limit ) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "n =", length(y), "\n"
    )
  ))
}


get_box_sum <- function(y, upper_limit ) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "n =", sum(y), "\n"
    )
  ))
}

# ---------------------------------------------------------------------------#
# modify_sample_names:
# This function modifies the sample names by removing file extension.
#> if data_type="CN_matrix" then 
#> this function modifies the row names of matrix 
#> else it modifies sample_name column
# Output: Modified CN_data with updated sample names i.e sample_0.987ncvey9ye --> sample_0
# ---------------------------------------------------------------------------#

modify_sample_names<-function(CN_data,data_type="CN_matrix"){
  if (data_type=="CN_matrix"){
    rownames(CN_data) <- sapply(rownames(CN_data), function(x) strsplit(x, "_")[[1]][1])
  }
  else {
    CN_data$sample_name <- sapply(CN_data$sample_name, function(x) strsplit(x, "_")[[1]][1])
  }
  return(CN_data)
}

# ---------------------------------------------------------------------------#
# mask_positions:
# This function masks positions in one dataframe based on ranges given another dataframe.
# Input:
#   df1: Df with  blacklisted regions with column 2 containg start position and column 3 containing stop position for region to be masked
#   df2: Df with positions to be masked
# Output:
#  Data frame df2 with an additional column "class" indicating if positions are blacklisted or callable.
# ---------------------------------------------------------------------------#
mask_positions <- function(df1, df2) {
  # Create an empty column named "class" in df2
  df2$class <- NA
  
  # use POS column in df2 or pos_in_ref
  if ("POS" %in% colnames(df2)){ 
    df2$POS<-as.numeric(df2$POS)
  } else {
    df2$POS<-as.numeric(df2$pos_in_ref)
  }
  #df2$pos_in_ref<-as.numeric(df2$pos_in_ref) 
  df1$V2<-as.numeric(df1$V2)
  df1$V3<-as.numeric(df1$V3)
  # Loop through each variant position in df2
  for (i in 1:nrow(df2)) {
    pos <- df2$POS[i]
    blacklisted <- "callable"
    # Loop through each row in df1 to check if the pos is in the blacklisted range 
    for (j in 1:nrow(df1)) {
      if (pos >= df1$V2[j] && pos <= df1$V3[j]) {
        blacklisted <- "blacklisted"
        break
      }
    }
    # Assign blacklisted or callable to the class column of df2 for the current row
    df2$class[i] <- blacklisted
  }
  
  return(df2)
}


# ---------------------------------------------------------------------------#
# add_variant_type:
#  This function adds a "type" column to a data frame based on the size of 
#  reference and alternate alleles, categorizing them as SNP, Indel, or MIXED.
# ---------------------------------------------------------------------------#

add_variant_type <- function(data) {
  # Calculate the size of the reference and alternate alleles
  data$ref_size <- nchar(data$ref_allele)
  data$alt_size <- nchar(data$alt_allele)
  # Add a new column "type" based on allele size
  data$type <- ifelse(data$ref_size == 1 & data$alt_size == 1, "SNP", "indel")
  # Check if "mixed" annotation is required
  if ("SNP" %in% data$type & "indel" %in% data$type) {
    # Get the position of "SNP" and "indel" rows
    snp_pos <- data[data$type=="SNP",]$pos_in_ref
    indel_pos <- data[data$type=="indel",]$pos_in_ref
    # Check if any positions have both "SNP" and "indel" rows
    mixed_pos <- intersect(snp_pos, indel_pos)
    # Annotate the mixed positions as "mixed"
    if(length(mixed_pos)>0){
      data[data$pos_in_ref %in% mixed_pos,]$type <- "mixed"
    }
  }
  # Remove the intermediate columns "ref_size" and "alt_size"
  data <- subset(data, select = -c(ref_size, alt_size))
  # Return the modified dataframe
  return(data)
}


# ---------------------------------------------------------------------------#
# add_regions_chrR:
# Input: Data frame with chromosome coordinates (POS or pos_in_ref column)
# ---------------------------------------------------------------------------#

add_regions_chrR <- function(data){
  if ("POS" %in% colnames(data)){ 
    POS<-as.numeric(data$POS)
  } else {
    POS<-as.numeric(data$pos_in_ref)
  }
  # <= only in IGS at start but all other should be > only but then variants are 0 bed based
  #so variant at the end of region be annotated to be in wrong region so we should use >=
  data$region[POS<=9338] <-"IGS"
  data$region[POS>9338] <-"5'ETS"
  data$region[POS>12995] <-"18S"
  data$region[POS>14864] <-"ITS1"
  data$region[POS>15934] <-"5.8S"
  data$region[POS>16091] <-"ITS2"
  data$region[POS>17258] <-"28S"
  data$region[POS>22309] <-"3'ETS"
  data$region[POS>22670] <-"IGS"
  return(data)}

# ---------------------------------------------------------------------------#
# predefined theme for ggplot2 plots, to have all plots with same theme
# ---------------------------------------------------------------------------#

common_theme2<-theme_bw()+
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12,face = "bold", color = "black"),
        axis.title = element_text(size = 14,face = "bold", color = "black"),
        strip.text.x = element_text(size = 10),
        plot.title = element_text(size = 15), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.text=element_text(size=13),
        legend.title = element_blank())




## *************************** ##
##                             ##
##    Preparing Dataframes     ##
##                             ##
## *************************** ##

##------------------------------##
## CN data matrix
##------------------------------##

# Modify sample names in CN.data.matrix eg sample_0.98bfgt876 --> sample_0
CN.data.matrix<-modify_sample_names(CN.data.matrix)
# Replace "Sample" with "sample_" in row names of CN.data.matrix 
rownames(CN.data.matrix) <- gsub("Sample", "sample_", rownames(CN.data.matrix))
# keep only the synthetic morph name and remove info about from which morph out of 24 T2T morphs it originated
colnames(CN.data.matrix) <- sapply(colnames(CN.data.matrix), function(x) strsplit(x, "\\.")[[1]][1])
colnames(CN.data.matrix) <- gsub("synth_", "rDNA_synth_seq_", colnames(CN.data.matrix))

# reorder matrix by sample numbers and synthetic morph numbers
CN.data.matrix <- CN.data.matrix[paste0("sample_", 0:(number_of_samples-1)), paste0("rDNA_synth_seq_", 0:(catalogue_size-1))] ## UPDATED

# start numbering samples and morphs from 1 instead of 0
rownames(CN.data.matrix) <- paste0("sample_", 1:number_of_samples)  
colnames(CN.data.matrix) <- paste0("rDNA_synth_seq_", 1:catalogue_size) ## UPDATED
head(CN.data.matrix)
dim(CN.data.matrix)
saveRDS(CN.data.matrix, "~/Desktop/simulated_WGS_data.sample_rDNA_sequences.rds")
simulated_WGS_data.sample_rDNA_sequences <- CN.data.matrix
##------------------------------##
## merged mutations dataframe
##------------------------------##

# Add 8838 to start coordinates 

mutation_data$mutation_added<-paste0(
  mutation_data$ref_morph_name,":",mutation_data$pos_in_ref+8838, 
  ":", mutation_data$ref_allele, ">",mutation_data$alt_allele
)

# Update 'pos_in_ref' by adding 8838 to each value
mutation_data$pos_in_ref <- as.numeric(mutation_data$pos_in_ref) + 8838

# Modify 'synth_morph_name' by replacing "synth_" with "rDNA_synth_seq_" 
# synth_morph_name should be same as columnnames of CN matrix
mutation_data$synth_morph_name <- gsub("synth_", 
                                       "rDNA_synth_seq_", 
                                       mutation_data$synth_morph_name )

# Split 'synth_morph_name' using "." and take the first part to remove ".chrR"
# mutation_data$synth_morph_name <- gsub(".chrR", "", mutation_data$synth_morph_name)

# Split 'synth_morph_name' using "." and remove everything after . i.e info about origion from 24 T2T morphs
mutation_data$synth_morph_name <-  sapply(mutation_data$synth_morph_name, 
                                          function(x) strsplit(x, "\\.")[[1]][1])
# start synthetic morph numbering from 1 instead of 0 because we updated them in CN matrix too
mutation_data$synth_morph_name <- sapply(mutation_data$synth_morph_name, function(i)
  paste0( paste(unlist(strsplit(i, split = "_"))[1:3], collapse = "_"),
          "_",
          as.numeric( unlist(strsplit(i, split = "_"))[[4]] ) + 1
  ))
sort(unique(mutation_data$synth_morph_name)) ## just to check

head(mutation_data)
dim(mutation_data)

mutation_data<-add_regions_chrR(mutation_data)
mutation_data <- mutation_data[mutation_data$region!="IGS",]
saveRDS(mutation_data, "~/Desktop/variant_catalogue.rds")
variant_catalogue <- mutation_data

length(unique(mutation_data$mutation_added))
##-----------------------------------------------------------##
## metadata dataframe i.e merging CN and mutations information
##------------------------------------------------------------##

# Calculate the variant counts based on copy numbers and morph names
# converting matrix to long-format dataframe
# variants count <-- long format of the matrix <-- CN count of each morph in each sample
metadata <- CN.data.matrix %>%
  tibble::rownames_to_column(var = "sample_name") %>%
  tidyr::gather(synth_morph_name, count, -sample_name)
# total contigs per sample : <-- required to calculate allele frequencies
metadata<-metadata %>% group_by(sample_name) %>% mutate(Total_contigs=sum(count))
# merging CN and mutation data
metadata <- merge(metadata, mutation_data, 
                  by = "synth_morph_name",
                  all.y = TRUE)
# remove rows showing 0 count for an allele in a sample
metadata<-metadata %>% filter(!count==0)

# some checks
head(metadata)
sort(unique(metadata$sample_name))

## add columns for variant type, rDNA region and masked positions
metadata<-add_variant_type(metadata)
metadata<-add_regions_chrR(metadata)
table(metadata$region)
table(metadata$region, metadata$type)
unique( metadata[metadata$region=="IGS", "pos_in_ref"])
metadata_IGS <- metadata
metadata <- metadata[metadata$region != "IGS",]
metadata<-mask_positions(df1 = blacklisted, df2=metadata)
head(metadata)
nrow(metadata)
nrow(mutation_data)
head(mutation_data)
head(metadata)
saveRDS(metadata, "~/Desktop/sample_catalogue.rds")
sample_catalogue <- metadata

##-----------------------------------------------------------##
## AF counts dataframe
##------------------------------------------------------------##


## Count of each variant in each sample
variants_AF_Counts<-metadata %>% 
  group_by(sample_name,
           mutation_added,
           pos_in_ref, 
           ref_allele,
           alt_allele,
           Total_contigs,
           type,
           region,
           class) %>% 
  summarise(count=sum(count))

## calculate allele freq from counts and total contigs
variants_AF_Counts<-variants_AF_Counts %>% 
  mutate(freq=count/Total_contigs)

# remove unecessary information from variant_id i.e "chr_name:10".
# :10 is the padding size
variants_AF_Counts$mutation_added<-gsub("chrR:8839_23170:|:10", "", 
                                        variants_AF_Counts$mutation_added)


##-----------------------------------------------------------##
## Number of mutations shared
##------------------------------------------------------------##

### num of samples carrying each mutation : 
df_Shared_Mutations  <- variants_AF_Counts %>% 
  group_by(mutation_added,pos_in_ref,type,class) %>% 
  summarize(num_samples=n())

###############################
##### RAQUELS Dataframes ######
###############################
# --> this is based on mutation df
variant_data<-add_variant_type(mutation_data)
variant_data<-add_regions_chrR(variant_data)
variant_data<-mask_positions(df1 = blacklisted, df2=variant_data)
colnames(variant_data) <- c("rDNA_copy_id", "variant_id", "seq", "position", "ref", "alt","region", "type","class","POS") 
variant_data$variant_id<-gsub("chrR:8839_23170:", "", variant_data$variant_id)
head( variant_data )
length(unique(variant_data$variant_id))
length(unique(variant_data$position))


rDNA_copy_summary <- variant_data %>% group_by(rDNA_copy_id, type,class) %>% summarize(count=n()) 
head(rDNA_copy_summary)


# --> this is based on CN matrix
cn_data <- CN.data.matrix


# --> sample rDNA configurations based on metadata df
sample_metadata <- metadata[,-ncol(metadata)]
colnames(sample_metadata) <- c("rDNA_copy_id",
                               "sample_id",
                               "copy_cn",
                               "total_cn",
                               "variant_id",
                               "seq",
                               "position",
                               "ref",
                               "alt",
                               "type",
                               "region",
                               "class")

sample_metadata$type <- gsub("Indel", "indel", sample_metadata$type)
sample_metadata$type <- gsub("MIXED", "mixed", sample_metadata$type)


# based on variant AF counts
variants_AF <- variants_AF_Counts
colnames(variants_AF) <- c("sample_id","variant_id", "position", "ref", "alt", "total_copies",
                           "type", "region","class", "count", "AF")

# --> based on shared_mutations df
variant_sharing <- df_Shared_Mutations
colnames(variant_sharing) <- c("variant_id", "position", "type", "class", "number_of_samples")

