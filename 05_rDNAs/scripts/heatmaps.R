################### libraries #######################
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(cowplot)
library(reshape2)
#unique_sequences
library("Biostrings")
########## Get donor and donor ID ###############
args <- commandArgs(trailingOnly = TRUE)
#GATK version either "4.5.0.0" or "4.3.0.0"
donor=args[1] 
# STEP of the pipeline this is for gVCF files so either "VC" or"filtering"
donor_id <- c("GWHDOOG000000", "GWHDQZJ000000")
replacement <- c("YAO_pat_chr", "YAO_mat_chr")

get_names<- function(names){
  names <- gsub(donor_id[1],replacement[1], names)
  names <- gsub(donor_id[2],replacement[2], names)
  for (n in 1:length(names)){
    if (grepl("rDNA45S_", names[n])) { 
      x <- sub("rDNA45S_", "", names[n]) 
      names[n]<- paste0("CHM13_", x)
    }else{
      names[n] <- gsub(":|-.*", "_", names[n])
      names[n] <- sub("_[^_]*$", "", names[n])
    }
  }
  return(names)
}
###### Directories #######
dir ="~/Desktop/ribosomal_RNAs/Alba/05_rDNA/"

dirI = paste0(dir, "/", "Heatmaps/")
dirF = paste0(dir, "/")
dir = paste0(dir, "/MSA/")
############################ Create Heatmaps ############################
# initialize row order to later filter unique sequences
row_cluster <-  c()

for (i in c(1,25)){
  
  dist <-  data.matrix(read.table(paste0(dir,"All_it", i, "/All_45S.dist"), fill=T, row.names=1, header =F, skip=1))
  
  names<- get_names(rownames(dist))
  
  rownames(dist) <-  names
  colnames(dist) <-  names
  dist2 <- round(dist, 2)
  dist3 <- as.matrix(apply(dist2, 2, function(i) gsub("99.", ".", as.character(i))))
  dist3 <- as.matrix(apply(dist3, 2, function(i) gsub("100", "", as.character(i))))
  rownames(dist3) <-names
  colnames(dist3) <-names
  # head(dist3)
  chromosome_names <- ifelse(grepl("^CHM13", names), "CHM13", 
                             ifelse(grepl("^YAO_mat", names), "YAO_mat", 
                                    ifelse(grepl("^YAO_pat", names), "YAO_pat", "chrR")))

  chr_col_names <- c( "CHM13", "YAO_mat", "YAO_pat", "chrR")
  chr_cols <- c( "#ECB22A", "#36C5F0","#6977EC", "#E01E5A")
  names(chr_cols) <- chr_col_names
  
  hrowanno <- HeatmapAnnotation( "Sample" = as.factor(chromosome_names),
                                  col = list(Sample = chr_cols[chromosome_names]),
                                  annotation_name_rot = 90,
                                  annotation_name_gp = gpar(fontsize = 0),
                                  which = "row" )
                              
  png(file=paste0(dirI, "Heatmap_45S", i, ".png"), res=180, width=3000, height=2000)
  draw(    Heatmap(  as.matrix(dist),
                         column_title = paste("Distance matrix for the 45S", i, "iterations"),
                         name = "% identity",
                         col = brewer.pal(9, "Blues"),
                         row_labels = gsub("rDNA", "sequence", rownames(dist)),
                         row_names_gp = gpar(fontsize = 10),
                         column_labels = gsub("rDNA", "sequence", colnames(dist)),
                         column_names_gp = gpar(fontsize = 10),
                         column_names_rot = 55,
                         cell_fun = function(j, i, x, y, w, h, col) {
                           grid.text(dist3[i, j], x, y,  gp = gpar(fontsize = 8)) }

    ) + hrowanno
  )
  dev.off()
  #################################
  if (i==25) row_cluster <-  row_order(heatmap)
}