### Load data ######
dir <- "~/Desktop/ribosomal_RNAs/Alba/04_T2T-YAO/stats/"
dirF <- "~/Desktop/ribosomal_RNAs/Alba/04_T2T-YAO/BED/"
dirI <- "~/Desktop/ribosomal_RNAs/Alba/04_T2T-YAO/IMAGES/BED_rRNA_STATS/"
# mat <- read.gff(paste0(dir, "YAO_mat_rRNA.gff3"), GFF3 = TRUE) 
# pat <- read.gff(paste0(dir, "YAO_pat_rRNA.gff3"), GFF3 = TRUE) 
# hp <-  read.gff(paste0(dir, "YAO_hp_rRNA.gff3"), GFF3 = TRUE) 
### libraries #####
library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(ggrepel)
#### Palettes #####
# palette1 <- c("#c57c3d",   "#a265c2", "#72a553",    "#ca5572",    "#6097ce")
palette2<- c( "18S"= "#a7c9c4","5.8S"= "#6bac97", "28S"="#00493e",
             "5S"="#ffcb39",
            "Mt_rRNA"= "#bd0142", "Mt"= "#bd0142",
            "rRNA"="#72a553",
            #"rRNA_pseudogene"="#6097ce", "ψ"="#6097ce",   "#ca5572",    "#6097ce")
            "rRNA_pseudogene"="#6097ce", "rDNA_pseudogene"="#6097ce",   "#ca5572",    "#6097ce")

palette3<- c("S5" ="#a797ca", "Pseudogenes"="#ffd07b", "S45"="#90ccf9")

### Functions #####
# Per chromosomes plot of counts used for pseudogenes, 45S, 5S
get_plot <- function(data, title) {
  
      return (ggplot(data, aes(x = factor(Var1, levels = c(1:22,"X","Y")), y = Freq)) +
              geom_bar(position=position_stack(reverse=TRUE),stat = "identity",  aes( fill = Var2)) +
              geom_text(aes(label=ifelse(Freq !=0, Freq, ""), y=Freq), size=4.5, fontface="bold", color="#18332a",
                        position=position_stack(vjust=1.1 ,reverse=FALSE)) +
              scale_fill_manual(values = palette2)+
                labs(title = title, y = "Counts", x = "",  fill="Gene Type") +
              theme_minimal() +
              theme(axis.text.x = element_text(vjust=0, hjust=1, size=12),
                    legend.position = "bottom"))
}
# Plot with all the distributions TO CHANGE
get_plot_genes <- function(data, title) {
  return (  ggplot(data, aes(x=factor(chr, levels = c(1:22, "X", "Y")),y=Freq, fill=reorder(Type, -Freq)))+
              geom_bar(stat = "identity", position = "stack") +
              geom_text(aes(label=ifelse(Freq !=0, Freq, "")), position = position_stack(.5), size=4.5, fontface="bold", color="#18332a") +
              scale_fill_manual(values = palette2)+
              labs(title =title, y="Counts", x="", fill="Gene")  +
              theme_minimal()+
              theme(axis.text.x = element_text( vjust = 1, hjust=1, size = 11), legend.position = "right"))
}
# Plot table with frequency of terms
plot_table <-  function( table, title, angle, size, subtitle){
  return(ggplot(table, aes(x=factor(Var1, levels = c("rDNA_pseudogene", "Mt",  "18S",  "5.8S", "28S", "5S", "rRNA" )),
                           y=Freq, fill=Var1, label=Freq))+
           geom_bar(stat ="identity")+
           geom_text(position = position_stack(), vjust=-0.25,size=size)+
           labs(title = title,subtitle = subtitle,  y="number of sequences", x="")+
           theme_minimal()+
           scale_x_discrete( labels =    c(paste0('\n' ,"rDNA", '\n' ,"pseudogenes"), "Mt",  "18S",  "5.8S", "28S", "5S"))+
           coord_cartesian(ylim = c(22, 510))+
           theme(axis.text.x = element_text(  size=9, angle=angle, vjust = 0.6), legend.position = "None", 
                 axis.text.y = element_text(size=8), axis.title = element_text(size=10), 
                 # legend.margin=margin(-100,0,-100,0),
                 axis.line = element_line(color = "black", size=0.35))+
           scale_fill_manual(values=palette2))
}
# Get label for each annotation file
get_label <- function(name) {
  if (name == "mat") {
    return("for the maternal annotations")
  } else if (name == "pat") {
    return("for the paternal annotations")
  } else if (name == "hp") {
    return("for the haplotype annotations")
  }
}
################--------------- Get al plots & filtered gff ------------------------------##############
for (name in c("mat", "pat", "hp")){

  d <- read.gff(paste0(dir, "YAO_", name, "_rRNA.gff3"), GFF3 = TRUE) 
  from <-  get_label(name)
  #### -------------------------------- Parsing & Data prep ------------------------------------ ##################
  
  # Only keep 45S genes & 5S genes, discard transcripts and exons: 
  d <- d[d$type=="gene" | d$type == "rRNA",]
  ## Parsing through the attributes column 
  d$gene_type <-  ""
  # d$transcript_type <- ""  # Transcript type only if type has exon or transcripts in it
  d$name <- ""
  d$key <-  ""
  d$product <- ""
  
  for (i in 1:nrow(d)){
      # Raw no parsed attributes separated with a ";"
      att<-     d$attributes[i]
      # Get gene type
      g <-  regmatches(att, regexpr("(?<=gene_biotype=).*?(?=;)", att, perl=TRUE)) # check other ways it is showing
      if (length(g) == 0)   g <- regmatches(att, regexpr("(?<=gene_type=).*?(?=;)", att, perl=TRUE))
      
      # # Transcript type only if type has exon or transcripts in it
      # t <-  regmatches(att, regexpr("(?<=transcript_biotype=).*?(?=;)", att, perl=TRUE))
      # if (length(t) == 0) t <- regmatches(att, regexpr("(?<=transcript_type=).*?(?=;)", att, perl=TRUE))
    
      # Gene Names 
      n <-  regmatches(att, regexpr("(?<=gene_name=).*?(?=;)", att, perl=TRUE))
      if (length(n) == 0){
        n <- regmatches(att, regexpr("(?<=Name=).*?(?=;)", att, perl=TRUE))
        if(length(n) == 0){
          n <- regmatches(att, regexpr("(?<=source_gene_common_name=).*?(?=;)", att, perl=TRUE))
        }
      }
      k <-  regmatches(att, regexpr("(?<=gbkey=).*?(?=;)", att, perl=TRUE))
      p <-  regmatches(att, regexpr("(?<=product=).*?(?=;)", att, perl=TRUE))
      # print(i)
      if(length(g) != 0) d$gene_type[i] <- g
      # if(length(t) != 0) d$transcript_type[i] <- t
      if(length(n) != 0) d$name[i] <- sub("^RNA", "", n) # remove the RNA pattern of the name
      if(length(k) != 0) d$key[i] <- k
      if(length(p) != 0) d$product[i] <- gsub(" ribosomal RNA", "", p)
  }
  # Remove the attributes column 
  d <- d[, - 9]
  
  # 5S and 45S genes
  genes <- d[! d$gene_type %in% c("rRNA_pseudogene", "Mt_rRNA"),] 
  # INFO and Context 
  # Mat genes dim:  350  12  (remove pseudos)  332  12  
  # Pat genes dim:  536  12  (remove pseudos) 518  12 
  #
  # unique(genes$name) there are some names that containg a "P" this stands for pseudogene:
  # [1] "RNA5SP533" "RNA5SP529" "RNA5SP536" "RNA5SP534" ...

  ## Remove pseudogenes: those with a P on the name: 18 rows less 
  more_pseudo <- genes[ grepl("P", genes$name), ] # keeping pseudogenes
  genes <- genes[ !grepl("P", genes$name), ]  
  # Change 5S and 5.8S genes to be all the same:
  genes$name<- gsub("^5S.*", "5S", genes$name) # 5S  after only keeping no pseudogenes
  genes$name<- gsub("5_8S_rRNA", "5.8S", genes$name) #5.8
  genes$name<- ifelse(genes$name == "", genes$product, genes$name) # morph them into one
  genes$gene_type<- ifelse(genes$gene_type == "", genes$key, genes$gene_type) # morph them into one
  
  # Get 5S genes 
  S5 <- genes[genes$name =="5S",]
  # Get 45S genes
  S45 <- genes[genes$name !="5S",]
  # Pseudogenes:
  pseudos <-  d[d$gene_type=="rRNA_pseudogene",]
  pseudos <- rbind(pseudos, more_pseudo)
  # We only find two types of pseudogenes: 
  pseudos$name<- gsub("^5S.*", "5S", pseudos$name)
  pseudos$name<- gsub("5-8S.*", "5.8S", pseudos$name) #5.8
  ############--------------------------  Chr prep Plotting------------------------------------- ####################
  
  # Table pseudogenes distribution
  table_chr_pse<- as.data.frame(table(pseudos$seqid, pseudos$name))
  table_chr_pse$Var1 <- gsub("chr(\\w+).*", "\\1", table_chr_pse$Var1)
  table_chr_pse$Var1 <- gsub("_hap[12]", "", table_chr_pse$Var1)
  # Table 5S distribution
  table_chr_S5<- as.data.frame(table(S5$seqid, S5$name))
  table_chr_S5$Var1 <- gsub("chr(\\w+).*", "\\1", table_chr_S5$Var1)
  table_chr_S5$Var1 <- gsub("_hap[12]", "", table_chr_S5$Var1)
  # Table 45S distribution
  table_chr_S45<- as.data.frame(table(S45$seqid, S45$name))
  table_chr_S45$Var1 <- gsub("chr(\\w+).*", "\\1", table_chr_S45$Var1)
  table_chr_S45$Var1 <- gsub("_hap[12]", "", table_chr_S45$Var1)
  # Remove mitchondrial in case any didn't get filtered out:
  table_chr_pse <- table_chr_pse[table_chr_pse$Var1 !="M",]
  table_chr_S5 <- table_chr_S5[table_chr_S5$Var1 !="M",]
  table_chr_S45 <- table_chr_S45[table_chr_S45$Var1 !="M",]
  # Table rRNA genes (no pseudogenes) terms found per chromosomes found 
  chr <-  data.frame(matrix(nrow =   length(table_chr_S5$Var1)+length(table_chr_S45$Var1),
                            ncol=0))
  chr$chr <-  c( table_chr_S5$Var1, table_chr_S45$Var1)
  chr$Type <- c( table_chr_S5$Var2, table_chr_S45$Var2)
  chr$Freq <- c( table_chr_S5$Freq, table_chr_S45$Freq) 
 
  chr$chr <- factor(chr$chr, levels = c(1:22, "X", "Y"))
  chr <-  arrange(chr, desc(Freq))
  
  # missing_data <- expand.grid(chr = c(1:22, "X", "Y"), Type = c("S5", "S45"), Value = 0)
  # missing_data$chr <- factor(missing_data$chr, levels = c(1:22, "X", "Y"))
  # long_genes_chr <- bind_rows(long_genes_chr, missing_data)
  # long_genes_chr <- distinct(long_genes_chr, chr, Type, .keep_all = TRUE)

  ############----------------------------  Plotting Chr  -------------------------------------- ####################
  # get plot function used here
  png(file=paste0(dirI, "chr_pseudogenes_stats_",  name,".png"), res=220, width=2300, height=1300);
  print(get_plot(data = table_chr_pse, title=paste("Found pseudogenes terms per chromosome",from)  ))
  dev.off();

  # Summary plot of all genes
  png(file=paste0(dirI, "chr_genes_distr_", name,".png"),  res=200, width=2300, height=1300);
  print(get_plot_genes(data=chr, title=paste("Found gene terms per chromosome", from)))
  dev.off();
  
  # png(file=paste0(dirI, name, "_chr_5s_stats.png"), res=200, width=2300, height=1300);
  # print(get_plot(data = table_chr_S5, title="Found pseudogenes terms per chromosome"))
  # dev.off();
  # 
  # png(file=paste0(dirI, name, "_chr_45s_stats.png"),  res=200, width=2300, height=1300);
  # print(get_plot(data = table_chr_S45, title="Found  45s terms per chromosome"))
  # dev.off();
  
  ############---------------------------- Other Plotting -------------------------------------- ####################
  # Table gene terms
  gt <- as.data.frame(table(d$gene_type[d$gene_type!=""]))
  kt <- as.data.frame(table(d$key[d$key!=""]))
  types <-  merge(gt, kt, by="Var1", all=T)
  types$Freq <- rowSums(types[, -1], na.rm = TRUE)
  # All gene types: 
  nt <- as.data.frame(table(genes$name))

  png(file=paste0(dirI, "product_type_stats_", name,".png"), res=320, width=1300, height=2300);
  print( plot_table(table = nt, title="Found rRNA terms per product type", subtitle=from, angle =0 , size = 3))
  dev.off();

  png(file=paste0(dirI, "raw_rRNA_types_", name,".png"), res=300, width=1300, height=1300);
  print(plot_table(table = types, title="Raw rRNA terms found",subtitle=from, angle =15, size = 6 ))
  dev.off();
  
  plt = data.frame(Var1 = c(nt$Var1, types$Var1), 
                   Freq =  c(nt$Freq, types$Freq))
  plt <-  plt[!plt$Var1=="rRNA" ,]
  plt$Var1 <- as.character(plt$Var1)
  #plt$Var1[plt$Var1=="rRNA_pseudogene"] <- "ψ"
  plt$Var1[plt$Var1=="rRNA_pseudogene"] <- "rDNA_pseudogene"
  plt$Var1[plt$Var1=="Mt_rRNA"] <- "Mt"
  
  # factor(plt$Var1, levels = c("ψ", "Mt",  "18S",  "5.8S", "28S", "5S", "rRNA" ))
  factor(plt$Var1, levels = c("rDNA_pseudogene", "Mt",  "18S",  "5.8S", "28S", "5S", "rRNA" ))
  ##############    IMPORTANT PLOT  ##############
  # 
  png(file=paste0(dirI, "raw_annotations_RNA_", name,".png"), res=300, width=1200, height=1200);
  print(plot_table(table = plt, title="",subtitle="", angle =15, size = 4 ))
  dev.off();
  
  png(file=paste0("~/Desktop/ribosomal_RNAs/Alba/IMAGES TFG/raw_annotations_RNA_", name,".png"), 
      res=330, width=1000, height=1400);
  print(plot_table(table = plt, title="",subtitle="", angle =0, size = 4 ))
  dev.off();
  ########## ---------------------------- Exporting BED file -----------------------------######
   
  bed <- genes[, 1:8]
  bed$score <- rep(".", nrow(bed))
  bed$phase <- rep(".", nrow(bed))
  # bed$attributes <- paste0("gene_type=", genes$gene_type, ";name=", genes$name, sep = "")
  bed$name <-  genes$name
  # 0 based
  bed$start <- bed$start - 1 
  
  write.table(bed, paste0(dirF, "YAO_rRNA_", name, ".bed"), row.names = F,col.names = F, sep="\t", quote=FALSE)
}


# > length(genes[genes$product =="18S",])
# [1] 12
# > length(genes[genes$product =="28S",])
# [1] 12
# > length(genes[genes$product =="5.8S",])
# [1] 12



