# libraries ----
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(cowplot)
library(reshape2)


dist <-  data.matrix(read.table("/45S/hp/out/YAO_it1/YAO_45S_hp.dist", 
                                fill=T, row.names=1, header =F, skip=1))
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
