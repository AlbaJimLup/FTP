# libraries ----
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(cowplot)
library(reshape2)

# 47S pre-rRNA sequence similarity ----
dist <- read.csv("~/bsc83_Projects/ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/T2T/data/rDNA_45S/alin/chrR_and_CHM13_it10/chrR_and_CHM13.dist.tab",
                 sep = "\t", header = F, row.names = 1)
rownames(dist)[1] <- "chrR_-0_copies_reference"
colnames(dist) <- rownames(dist)
dist2 <- round(dist, 2)
dist3 <- as.matrix(apply(dist2, 2, function(i) gsub("99.", ".", as.character(i))))
dist3 <- as.matrix(apply(dist3, 2, function(i) gsub("100", "", as.character(i))))
rownames(dist3) <- rownames(dist)
colnames(dist3) <- colnames(dist)
head(dist3)

chr_cols <- c("#36C5F0", "#E01E5A", "#ECB223", "#2EB67D", "#6977EC", "white")
names(chr_cols) <- c("chr13", "chr14", "chr15", "chr21", "chr22","reference")
hrowanno <- HeatmapAnnotation("number of copies" = anno_barplot( as.numeric(sapply(colnames(dist), function(i) 
  unlist(strsplit(unlist(strsplit(i, split = "-"))[[2]] , split = "_"))[[1]] )),
  border = FALSE, 
  add_numbers = TRUE, 
  gp = (gpar( fill = "gray", col = "gray" )),
  axis_param = list(labels_rot = 90)),
                  "chromosome" =  as.factor(sapply(colnames(dist), function(i) 
                    unlist(strsplit(i, split = "_"))[[4]])),
  col = list(chromosome = chr_cols),
  annotation_name_rot = 90,
  annotation_name_gp = gpar(fontsize = 8),
  which = "row")

Heatmap(as.matrix( dist ),
        name = "% identity",
        right_annotation = hrowanno,
        col = brewer.pal(9, "Blues"),
        row_labels = gsub("rDNA", "sequence", sapply(rownames(dist), function(i) 
          unlist(strsplit(i, split = "-"))[[1]]) ),
        column_labels = gsub("rDNA", "sequence", sapply(colnames(dist), function(i) 
          unlist(strsplit(i, split = "-"))[[1]]) ),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(dist3[i, j], x, y,  gp = gpar(fontsize = 8))
        }
        )

# CHM13 MSA variants ----
CHM13_variants <- read.delim("~/bsc83_Projects/ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/T2T/data/rDNA_45S/alin/chrR_and_CHM13_it10/CHM13_variants.msa.tsv")
CHM13_variants$chr <- "chrR"
colnames(CHM13_variants)[1:4] <- c("position","reference", "alternative", "chrR")
CHM13_variants <- CHM13_variants[,c(ncol(CHM13_variants),1:(ncol(CHM13_variants)-1))]
head(CHM13_variants)

# number of copies per rDNA sequence ----
CHM13_copy_number <- cbind.data.frame("rDNA_sequence" = sapply(colnames(CHM13_variants)[6:ncol(CHM13_variants)], function(i) unlist(strsplit(i, split = "\\."))[[1]]),
                 "n_copies" = sapply(colnames(CHM13_variants)[6:ncol(CHM13_variants)], function(i) 
                   as.numeric(unlist(strsplit(unlist(strsplit(i, split = "\\."))[[2]], split = "_"))[[1]])),
                 "chromosome" = sapply(colnames(CHM13_variants)[6:ncol(CHM13_variants)], function(i) unlist(strsplit(i, split = "_"))[[3]])
)
rownames(CHM13_copy_number) <- NULL
sum(CHM13_copy_number$n_copies)

# assign region
assign_region <- function(position){
  # <= only in IGS at start but all other should be > only but then variants are 0 bed based
  #so variant at the end of region be annotated to be in wrong region so we should use >=
  if(position <= 9338){
    "IGS"
  }else if(position <= 12995){
    "5'ETS"
  }else if(position <= 14864){
    "18S"
  }else if(position <= 15934){
    "ITS1"
  }else if(position <= 16091){
    "5.8S"
  }else if(position <= 17258){
    "ITS2"
  }else if(position <= 22309){
    "28S"
  }else if(position <= 22670){
    "3'ETS"
  }else{
    "IGS"
  }
}
CHM13_variants$region <- sapply(CHM13_variants$position, function(i) assign_region(i))

# callable or blacklisted
repeat_info <- read.delim("~/bsc83_Projects/ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/T2T/data/rDNA_45S/alin/chrR_and_CHM13_it10/CHM13_variants.msa.repeat_info.tsv", header = F)
CHM13_variants$class <- ifelse(repeat_info$V6 == "chrR", "blacklisted", "callable")
ncol(CHM13_variants)

CHM13_variants <- CHM13_variants[,c(1:4,30,31,5:29)]
head(CHM13_variants)

# fix 28S 'callable' variants
CHM13_variants[CHM13_variants$position==22309, "class"] <- "blacklisted"
head(CHM13_variants)
# write.table(CHM13_variants,
#             "~/Desktop/CHM13_variants.msa.tab",
#             col.names = T, row.names = F, quote = F, sep = "\t")

# plots ----

# fix reference alleles so they are chrR (GATK reported) reference alleles ----
CHM13_variants[CHM13_variants$position==15060,]
CHM13_variants[CHM13_variants$position==18040,]
CHM13_variants[CHM13_variants$position==18120,]
CHM13_variants[CHM13_variants$position==20575,]
CHM13_variants[CHM13_variants$position==22309,]

# number of callable and blacklisted variants  ----
regions <- c("5'ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3'ETS")

df <- CHM13_variants %>% 
  group_by(region, class)   %>% 
  #group_by(TYPE, region, class)   %>% 
  summarize(variants=n())
df
#df$TYPE <- factor(df$TYPE, levels = c("SNP", "INDEL", "MIXED", "MNP"), order = T)
df$region <- factor(df$region, levels = regions, order = T)
df$class <- factor(df$class, levels = c("callable", "blacklisted"), order = T)
df[df$region=="5'ETS",]
ggplot(df, 
       aes( x = region, y = variants, fill = class)) +
       #aes( x = region, y = variants, fill = TYPE)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  #facet_wrap(~class, nrow = 1) +
  geom_text(aes(label = variants),
            position = position_stack(),
            #position = position_dodge(width = 0.9),  # Adjust the width for dodge positioning
            vjust = 1,  # Adjust the vjust value for label positioning within dodges
            size = 3) +
  scale_fill_brewer(palette = "Pastel1",
                    guide = guide_legend(reverse=FALSE), 
                      )

# to do -> classify variants as SNPs, indels or mixed ----
# head( CHM13_variants )
# CHM13_variants$number_of_alleles <- sapply(CHM13_variants$alt, function(i) length(unlist(strsplit(i, split = '/')))+1)
# table(CHM13_variants$number_of_alleles)
# df <- CHM13_variants %>% 
#   group_by(number_of_alleles, region, class)   %>% 
#   summarize(variants=n())
# df
# df$number_of_alleles <- factor(df$number_of_alleles, levels = unique(df$number_of_alleles), order = T)
# df$region <- factor(df$region, levels = regions, order = T)
# df$class <- factor(df$class, levels = c("callable", "blacklisted"), order = T)
# df[df$region=="5'ETS",]
# ggplot(df, 
#        aes( x = region, y = variants, fill = number_of_alleles)) +
#   geom_bar(stat = "identity") +
#   theme_bw() +
#   facet_wrap(~class, nrow = 1) +
#   geom_text(aes(label = variants),
#             position = position_stack(),
#             #position = position_dodge(width = 0.9),  # Adjust the width for dodge positioning
#             vjust = 1,  # Adjust the vjust value for label positioning within dodges
#             size = 3) +
#   scale_fill_brewer(palette = "Pastel2",
#                     guide = guide_legend(reverse=FALSE))

nrow(CHM13_variants)
table(CHM13_variants$class)

# number of variants in each sequence ----
head(CHM13_variants)
colnames(CHM13_variants)[7] <- "chrR"
variants_per_sequence <- cbind.data.frame("sequence" = paste0("sequence_", 1:24),
                  "variants" = sapply(8:31, function(i) sum(CHM13_variants$chrR != CHM13_variants[,i]))
)
variants_per_sequence$sequence <- factor(variants_per_sequence$sequence,
                                         levels  = paste0("sequence_", 1:24),
                                         order = T)
ggplot(variants_per_sequence,
       aes( x = sequence, y = variants)) +
  geom_bar(stat = "identity", fill = "gray") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab('') + ylab("number of variants") +
  geom_text(aes(label = variants),
            position = position_identity(),
            #position = position_dodge(width = 0.9),  # Adjust the width for dodge positioning
            vjust = 1,  # Adjust the vjust value for label positioning within dodges
            size = 3)

# sharing of variants ----
barplot(table(sapply(1:nrow(CHM13_variants), function(row)
  sum(CHM13_variants[row,8:31]==CHM13_variants[row,'chrR'])
  )),
  main = "sharing of referece allele")


# distribution of alt allele freq per region and class ----
CHM13_AF <- readRDS("~/Downloads/R_freq_CHM13_variants.rds")
CHM13_AF <- CHM13_AF[order(CHM13_AF$position),]
CHM13_AF$region <- factor(CHM13_AF$region, levels = regions, order = T)
CHM13_AF$class <- factor(CHM13_AF$class, levels = c("callable", "blacklisted"), order = T)
head(CHM13_AF)

nrow(CHM13_AF)
length(unique(CHM13_AF$position))
365-162
nrow(CHM13_AF[CHM13_AF$allele != CHM13_AF$reference,])

head(CHM13_AF)
ggplot(CHM13_AF[CHM13_AF$allele != CHM13_AF$reference,],
       aes( x = region, y = allelefreq, fill = class)) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw()

ggplot(CHM13_AF[CHM13_AF$allele != CHM13_AF$reference,],
       aes( x = region, y = allelefreq)) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw() +
  facet_wrap(~class)
table(table(CHM13_AF$position))

library(dplyr)
nrow(CHM13_AF %>%
  group_by(position) %>%
  summarise(sum_value = sum(allelefreq)))


head(CHM13_AF)
nrow(CHM13_AF)

# Alba computed allele frequencies
allele_freq <- readRDS("~/Downloads/R_freq_CHM13_variants.rds")
allele_freq <- allele_freq[order(allele_freq$position),]

table(allele_freq$position)
head( allele_freq )
nrow(allele_freq)
sum(CHM13_variants$number_of_alleles)

############################### 
# Getting first table
###############################
#d <- readRDS("CHM13_variants.rds")
d <- CHM13_msa_variants

t <-  data.frame("position" = d$position,
                 #"type" = d$TYPE, 
                 "reference" = d$chrR,
                 "allele" = d$ref, 
                 "region" = d$region,
                 "class" = d$class, 
                 "allelefreq" = rep(0, nrow(d)),
                 "rDNAs" = I(rep(list(NULL), nrow(d))),
                 "rDNAsUNITs" = rep(0, nrow(d)),
                 "copycount" = rep(0, nrow(d)))

get_rDNAs_freq <- function(d, row, nucl ){
  
  rDNAs <- list()
  counts <- 0 
  # Get frequencies
  for ( column in  10:33){
    allele <-  d[row, column] # Get allele in copy
    
    if (allele == nucl ){ 
      rdna <-  colnames(d[column])
      rDNAs <- c(rDNAs,  sub("\\..*", "", rdna) ) #name rDNA
      counts <- counts +  as.numeric(gsub(".*\\.(\\d+)_copies.*", "\\1", rdna))# num copies
    }} 
  
  return(list(rDNAs = rDNAs, 
              freqs = counts/219, 
              rDNAunits = length(rDNAs), 
              copies = counts))
}


for ( row in 1:nrow(d)){
  # For REF
  stats <- get_rDNAs_freq(d, row,  d$ref[row])
  t$rDNAs[row] <- paste0(stats$rDNAs, collapse = ", ")
  t$allelefreq[row] <-  stats$freqs 
  t$rDNAsUNITs[row] <- stats$rDNAunits
  t$copycount[row] <- stats$copies
  # For ALTs
  alts <- unlist(strsplit(unique(d$alt[row]), ","))
  
  # more than one alternative
  if (length(alts) > 1 ){
    for (alt in alts){
      # Get frequencies
      stats <- get_rDNAs_freq(d, row, alt )
      
      newrow<-  data.frame("position" = d$position[row],
                           "type" = d$TYPE[row], 
                           "reference" = d$chrR[row],
                           "allele" = alt, 
                           "region" = d$region[row],
                           "class" = d$class[row], 
                           "allelefreq" = stats$freqs,
                           "rDNAs"= paste0(stats$rDNAs, collapse = ", "),
                           "rDNAsUNITs" = stats$rDNAunits,
                           "copycount" = stats$copies)
      t <- rbind(t,newrow) 
    }
  }
  else {
    alt <- d$ALT[row]
    # Get frequencies
    stats <- get_rDNAs_freq(d, row, alt )
    
    newrow<-   data.frame("position" = d$position[row],
                          #"type" = d$TYPE[row], 
                          "reference" = d$chrR[row],
                          "allele" = alt, 
                          "region" = d$region[row],
                          "class" = d$class[row], 
                          "allelefreq" = stats$freqs,
                          "rDNAs"= paste0(stats$rDNAs, collapse = ", "),
                          "rDNAsUNITs" = stats$rDNAunits,
                          "copycount" = stats$copies)
    
    t <- rbind(t,newrow) # add reference line to current data
  }}

saveRDS(t, "freq_CHM13_variants.rds")
