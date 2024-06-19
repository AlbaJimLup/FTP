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
regions<-read.table("~/Desktop/ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/T2T/scripts/data/pre-rRNA_47S.regions.bed")

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

mutation_data<-add_regions_chrR(mutation_data)

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


###########################################
## PLOTS
###########################################

# total number of variants (complete catalogue) stratified by type and class ----
# create summary df
variant_summary <- variant_sharing %>%
  ungroup() %>%
  select(class, type, variant_id) %>%
  #distinct(position, .keep_all = TRUE) %>% # with this line reports variants (not variant alleles)
  ungroup() %>%
  group_by(type, class) %>%
  summarise(number_of_variant_positions = n()) 
variant_summary_df <- as.data.frame(variant_summary)
#variant_summary_df$type <- c("indel","indel","mixed","SNP","SNP")  ### I AM MODIFYING THIS LINE BECAUSE WE HAVE TWO ROWS WITH "MIXED"
## CORRECTION
variant_summary_df <- variant_summary_df %>%
  mutate(type = case_when(
    type == "Indel" ~ "indel",
    type == "MIXED" ~ "mixed",
    type == "snp" ~ "SNP",
    TRUE ~ type  # Keep other values unchanged
  ))
## END CORRECTION
variant_summary_df$type <- factor(variant_summary_df$type,
                                  levels = c("SNP", "indel", "mixed"),
                                  order = T)
variant_summary_df$class <- factor(variant_summary_df$class,
                                   levels = c("callable", "blacklisted"),
                                   order = T)
# plot
p1 <- ggplot(variant_summary_df,
             aes(y = number_of_variant_positions, x = type, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +  # Use position = "dodge" to create a dodge plot
  labs(y = "number of variant alleles", x = "") +
  geom_text(aes(label = number_of_variant_positions),
            position = position_dodge(width = 0.9),  # Adjust the width for dodge positioning
            vjust = -0.5,  # Adjust the vjust value for label positioning within dodges
            size = 3) +  # Adjust the size of the labels
  theme_bw() +
  scale_fill_brewer(palette = "Pastel2",
                    guide = guide_legend(reverse=TRUE)) +
  theme(legend.position = "bottom", 
        legend.title = element_blank()) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

# number of callable variants per region  ----
# define the desired order of the bars
regions <- c("5'ETS","18S", "ITS1" ,"5.8S", "ITS2","28S","3'ETS","IGS") ## I ADDED IGS else it add <NA> in callable_variant_summary_df
callable_variant_summary_df <- as.data.frame(variants_AF %>% ungroup() %>% group_by(variant_id, position, type, region, class) %>% filter(class=="callable") %>%  summarize(num_samples=n()) %>%
                                               ungroup() %>%
                                               select(class, type, region, position) %>%
                                               #distinct(position, .keep_all = TRUE) %>%  # with this line reports variants (not variant alleles) (255 variant positions)
                                               ungroup() %>%
                                               group_by(type,region=factor(region,levels=regions)) %>%
                                               summarise(number_of_variant_positions = n()) ) 
sum(callable_variant_summary_df$number_of_variant_positions) # 269
callable_variant_summary_df$type <- gsub("Indel", "indel", callable_variant_summary_df$type)
callable_variant_summary_df$type <- gsub("MIXED", "mixed", callable_variant_summary_df$type)
callable_variant_summary_df$type <- factor(callable_variant_summary_df$type,
                                           levels = rev(c("SNP", "indel", "mixed")),
                                           order = T)
callable_variant_summary_df$region <- factor(callable_variant_summary_df$region,
                                             levels = regions,
                                             order = T)
p2 <- ggplot(callable_variant_summary_df,
             aes(y = number_of_variant_positions, fill = type, x = region)) +
  geom_bar(stat = "identity") +  # Use position = "dodge" to create a dodge plot
  labs(y = "number of callable variant alleles", x = "", fill = "region") +
  geom_text(aes(label = number_of_variant_positions),
            position = position_stack(vjust = 0.5),  # Adjust the width for dodge positioning
            
            size = 3) +  # Adjust the size of the labels
  theme_grey() +
  theme_bw() +
  scale_fill_brewer(palette = "Pastel1",
                    guide = guide_legend(reverse=TRUE)) +
  theme(legend.position = "bottom",
        legend.title = element_blank())


# callable variant alleles per rDNA synthetic sequence ----
rDNA_copy_summary <- variant_data %>% group_by(rDNA_copy_id, type, class) %>% 
  summarize(count=n()) %>% 
  filter(class == "callable")
rDNA_copy_summary$rDNA_copy_id <- factor(rDNA_copy_summary$rDNA_copy_id,
                                         levels = paste0("rDNA_synth_seq_",1:catalogue_size),
                                         order = T)

rDNA_copy_summary$type <- factor(rDNA_copy_summary$type,
                                 levels = rev(c("SNP", "indel", "mixed")), order = T)

sapply(paste0("rDNA_synth_seq_", 1:catalogue_size), function(i) nrow( variant_data %>% filter(rDNA_copy_id == i, ) ) )
rDNA_copy_summary[rDNA_copy_summary$rDNA_copy_id == paste0("rDNA_synth_seq_", 1),]

p3 <- ggplot(rDNA_copy_summary,
             aes(y = count, fill = type, x = rDNA_copy_id)) +
  geom_bar(stat = "identity") +  # Use position = "dodge" to create a dodge plot
  labs(y = "number of variant alleles", x = "", fill = "type") +
  geom_text(aes(label = count),
            position = position_stack(vjust = 0.5), 
            size = 3) +  # Adjust the size of the labels
  theme_bw() +
  scale_fill_brewer(palette = "Pastel1",
                    guide = guide_legend(reverse=TRUE)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

# sharing of variants between rDNA synthetic sequences ----
variant_sharing <- merge(variant_sharing, variant_data %>% group_by(variant_id) %>% summarise(number_of_rDNA_synth_seq = n()), by = "variant_id")
p3_prime <- ggplot(variant_sharing %>% filter(class == "callable"),
                   aes(x = "", y = number_of_rDNA_synth_seq)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.5) +
  ylab("number of rDNA synthetic sequences")  +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
  ) +
  stat_summary(
    fun.data = get_box_stats,
    fun.args = list(upper_limit = max(variant_sharing$number_of_rDNA_synth_seq) * 1.15),
    geom = "text",
    hjust = 0.5, vjust = 0.9, size = 3) 

# p1, p2, p3 & p3_prime ----
plot_grid(p1, p2, p3, p3_prime,
          rel_widths = c(0.15, 0.25, 0.5, 0.1),
          nrow = 1)

# distribution of the number of rDNA copies per sample ----
head( cn_data )
CN <- data.frame(number_of_copies = rowSums(cn_data))
p4 <- ggplot(CN, aes(x = "", y = number_of_copies)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.5) +
  ylab("number of total number of rDNA copies")  +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
  ) +
  stat_summary(
    fun.data = get_box_stats,
    fun.args = list(upper_limit = max(CN$number_of_copies) * 1.15),
    geom = "text",
    hjust = 0.5, vjust = 0.9, size = 3) 

# distribution of the number of unique rDNA copies per sample ----
head(cn_data)
CN$sample <- rownames(CN)

# barplot with total number of copies
CN$sample <- factor(CN$sample,
                    levels = paste0("sample_", 1:number_of_samples),
                    order = T)
p5.1 <- ggplot(CN,
               aes(x = sample, y = number_of_copies)) +
  geom_bar(stat= "identity") +
  geom_text(aes(label = number_of_copies),
            position = position_stack(vjust = 0.5), 
            size = 3, angle = 90) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("total number of rDNA copies") + xlab("") 

# boxplot wiht number of copies of each unique rDNA copy
cn_data$sample <- rownames(cn_data)
cn_data_df <- melt(cn_data)
cn_data_df <- cn_data_df[cn_data_df$value > 0,]
cn_data_df$sample <- factor(cn_data_df$sample, 
                            levels = paste0("sample_", 1:number_of_samples),
                            order = T)
cn_data_df$variable <- factor(cn_data_df$variable,
                              levels = paste0('rDNA_synth_seq_', 1:catalogue_size),
                              order = T)
p5.2 <- ggplot(cn_data_df[cn_data_df$sample %in% paste0("sample_", 1:number_of_samples),],
               aes(x = sample, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.5)   +
  ylab("number of rDNA copies")  + xlab("") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  stat_summary(
    fun.data = get_box_stats,
    fun.args = list(upper_limit = max(cn_data_df$value) * 1.15),
    geom = "text",
    hjust = 1, vjust = 0.75, size = 2, angle = 90) 

p5 <- plot_grid(p5.1, p5.2,
                nrow = 2,
                rel_heights = c(0.4, 0.6))  

# p4 & p5 ----
plot_grid(p4, p5,
          rel_widths = c(0.1, 0.9))

# distribution of the number of callable variant alleles per sample ----
variants_AF$variant_position <- sapply(variants_AF$variant_id, function(i) unlist(strsplit(i, split = ":"))[[2]])
num_variant_alleles_per_sample <- variants_AF %>% filter(class=="callable") %>% group_by(sample_id) %>% summarize(num_variants=n()) 
# one variant can be present more than once , this is counting number of unique variants

p6 <- ggplot(num_variant_alleles_per_sample, 
             aes(x = "", y = num_variants)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.5) +
  theme_bw() +
  ylab("number of variant alleles") +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
  ) +
  stat_summary(
    fun.data = get_box_stats,
    fun.args = list(upper_limit = max(num_variant_alleles_per_sample$num_variants) * 1.15),
    geom = "text",
    hjust = 0.5, vjust = 0.9, size = 3) 


# number of variants (variant alleles) per sample ----
variant_sample_summary <- variants_AF %>% filter(class == "callable") %>% group_by(sample_id, type)  %>%  summarise(count = n())  
variant_sample_summary$type <- gsub("Indel", "indel", variant_sample_summary$type)
variant_sample_summary$type <- gsub("MIXED", "mixed", variant_sample_summary$type)
variant_sample_summary$type <- factor(variant_sample_summary$type,
                                      levels = rev(c("SNP", "indel", "mixed")),
                                      order = T)
variant_sample_summary$sample_id <- factor(variant_sample_summary$sample_id,
                                           levels = paste0("sample_", 1:number_of_samples),
                                           order = T)
p7.1 <- ggplot(variant_sample_summary,
               aes(y = count, fill = type, x = sample_id)) +
  geom_bar(stat = "identity") +  # Use position = "dodge" to create a dodge plot
  labs(y = "number of variant alleles", x = "", fill = "type") +
  geom_text(aes(label = count),
            position = position_stack(vjust = 0.5), 
            size = 3, angle = 90) +  # Adjust the size of the labels
  theme_bw() +
  scale_fill_brewer(palette = "Pastel1",
                    guide = guide_legend(reverse=TRUE)) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
#axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

# distribution of AF per sample ----
variants_AF$sample_id <- factor(variants_AF$sample_id,
                                   levels = paste0("sample_", 1:number_of_samples),
                                   order = T)
p7.2 <- ggplot(variants_AF[variants_AF$class=="callable",],
               aes(x = sample_id, y = AF)) +
  geom_boxplot(outlier.size = 0.5) +
  #geom_jitter(size = 0.1)   +
  ylab("variant allele frequency")  + xlab("") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  stat_summary(
    fun.data = get_box_stats,
    fun.args = list(upper_limit = 1.15),
    geom = "text",
    hjust = 1, vjust = 0.75, size = 2, angle = 90) 


# p6 & p7 ----
p7 <- plot_grid(p7.1, p7.2,
                rel_widths = c(0.4, 0.6),
                nrow = 2)

# sharing of variants across samples ----
head( variant_sharing )
p8 <- ggplot(variant_sharing %>% filter(class == "callable"), aes(x = "", y = number_of_samples)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.5) +
  ylab("number of samples")  +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
  ) +
  stat_summary(
    fun.data = get_box_stats,
    fun.args = list(upper_limit = max(variant_sharing$number_of_samples) * 1.15),
    geom = "text",
    hjust = 0.5, vjust = 0.9, size = 3) 

plot_grid(p6, p7, p8,
          rel_widths = c(0.1, 0.8, 0.1), nrow = 1)


# distribution of AF per variant across samples for bottom and top variant AF ----
# calculate median AG per variany across samples
variant_median_AF_df <- as.data.frame(variants_AF %>% filter(class == "callable") %>% group_by(variant_id)  %>%  summarise(median_AF = median(AF)) )
variant_median_AF_df <-  variant_median_AF_df[order(variant_median_AF_df$median_AF, decreasing = F),] 
# select the bottom 10 and top 10 variants with lowest and highest median AF (leats and most frequent variants across samples)
subset_variants <- variant_median_AF_df$variant_id[c(1:10, (nrow(variant_median_AF_df)-9):nrow(variant_median_AF_df))]
variants_AF$variant_id <- factor(variants_AF$variant_id,
                                 levels = variant_median_AF_df$variant_id) 

pp1 <- ggplot(variants_AF %>% filter(variant_id %in% subset_variants, class == "callable"),
              aes(x = variant_id, y = count)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.6) +
  labs(y = "number of rDNA copies", x = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") + 
  geom_text(data = variant_sharing %>% filter(variant_id %in% subset_variants, class == "callable"),
            aes(x = variant_id , label = paste0("n = ", number_of_samples), y = 400),size = 2.5, hjust = 0.5) +
  geom_vline(xintercept = 10.5, lty = 2)

pp2 <- ggplot(variants_AF %>% filter(variant_id %in% subset_variants, class == "callable"),
              aes(x = variant_id, y = AF)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.6) +
  labs(y = "allele frequency", x = "") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") + 
  geom_text(data = variant_sharing %>% filter(variant_id %in% subset_variants, class == "callable"),
            aes(x = variant_id , label = paste0("n = ", number_of_samples), y = 1),size = 2.5, hjust = 0.5) +
  geom_vline(xintercept = 10.5, lty = 2)

plot_grid(pp1, pp2, nrow = 1)
ggplot(variant_sharing,
       aes( x = 1, y = number_of_samples)) +
  geom_boxplot()
head(variant_sharing)




