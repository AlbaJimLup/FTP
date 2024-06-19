### Load paths ######
dir <- "~/Desktop/ribosomal_RNAs/Alba/04_T2T-YAO/stats/"
dirF <- "~/Desktop/ribosomal_RNAs/Alba/04_T2T-YAO/BED/"
dirI <- "~/Desktop/ribosomal_RNAs/Alba/04_T2T-YAO/IMAGES/STATS/"
### libraries #####
library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(ggrepel)
### Palettes #####

palette1 <- c("13"="#00a278","14"="#7fbf83","15"="#c3df7a", "21"="#5d8fff","22"="#a3b7dc" )
palette1<- c( "13"= "#74999c","14"= "#a7c9c4", "15"="#6bac97",  "21"="#ffcb39","22"="#8ab1d8")
palette1<- c( "13"= "#6097ce","14"= "#ffcb39", "15"="#a7c9c4",  "21"="#ce586c","22"="#6bac97")
palette2<- c("28S"="#7d928b", "18S"= "#54b06c", "5S"="#a797ca","5.8S"= "#91e977")

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
  return (  ggplot(data, aes(x=factor(Var1, levels = c(1:22, "X", "Y")),y=Freq, fill=reorder(Var2, -Freq)))+
              geom_bar(stat = "identity", position = "stack") +
              geom_text(aes(label=ifelse(Freq !=0, Freq, "")), position = position_stack(.5), size=4.5, fontface="bold", color="#18332a") +
              scale_fill_manual(values = palette2)+
              labs(title =title, y="Counts", x="", fill="Gene")  +
              theme_minimal()+
              theme(axis.text.x = element_text( vjust = 1, hjust=1, size = 11), legend.position = "right"))
}
# Plot table with frequency of terms
plot_nt <-  function( table, title, angle, size, subtitle){
  return(ggplot(table, aes(x=type,y=counts, fill=type, label=counts))+
           geom_bar(stat ="identity")+
           geom_text(position = position_stack(1.03), size= size)+
           labs(title = title,subtitle = subtitle,  y="Counts", x="")+
           theme_minimal()+
           theme(axis.text.x = element_text(  size=12, angle=angle), legend.position = "None", 
                 axis.text.y = element_text(size=13), axis.title = element_text(size=14))+
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
# Get expected CN:
get_CN <- function(name) {
  if (name == "mat") {
    return(79)
  } else if (name == "pat") {
    return(149)
  } else if (name == "hp") {
    return(98)
  }
}

get_ID <-  function(name){
  if (name == "mat") {
    return("GWHDQZJ000000")
  } else if (name == "pat") {
    return("GWHDOOG000000")
  } else if (name == "hp") {
    return("GWHDQZI000000")
  }
}

get_chroms <-  function(name){
  
  id <-  get_ID(name)
  
  chroms <-  c(paste0(id, "0", 1:9), paste0(id, 10:23) )
}

###### Initialize variables ######

colnames= c("chr", "source", "type","start","end","score","strand","phase","gene")
# 13 14 15 21 and 22:
CN <- data.frame(matrix(ncol=0, nrow=(5*3)))

######## For each bed file from stats_raw ######## 
for (name in c("mat", "pat", "hp")){

######################## 1 Remove extra copies that are not in continuous operon #####################
    d <- read.table(paste0(dirF, "YAO_rRNA_", name,".bed")) # Import output sata_raw.R
    from <-  get_label(name)
    colnames(d) <- colnames # Change colnames
    #Rename chr col
    IDs <- d$chr
    d$chr <- gsub("chr(\\w+).*", "\\1", d$chr)
    d$chr <-   gsub("_hap[12]", "", d$chr) # only chrN as chr 
    
    ############-------------------------- Context -------------------------########
    #table(m$chr, m$gene)
    #       18S 28S 5.8S 5S
    # chr1    0   0    0 59
    # chr13  45  44   44  0     ! TO REMOVE
    # chr14  31  30   31  0     ! 2 TO REMOVE
    # chr15  17  16   16  0     ! TO REMOVE
    # chr2    0   0    0  2
    # chr20   0   0    3  0
    # chr21   8   7    7  0     ! TO REMOVE
    # chr22  53  52   52  0     ! TO REMOVE
    # chr6    0   0    0  1
        
    # table(p$chr, p$gene)
    #       18S 28S 5.8S 5S
    # chr1    0   0    0 59
    # chr13  45  44   44  0    ! TO REMOVE
    # chr14  31  30   31  0    ! 2 TO REMOVE
    # chr15  17  16   16  0    ! TO REMOVE
    # chr2    0   0    0  2
    # chr20   0   0    3  0
    # chr21   8   7    7  0    ! TO REMOVE
    # chr22  53  52   52  0    ! TO REMOVE
    # chr6    0   0    0  1
    
    # table(h$chr, h$gene)
    #         18S 28S 5.8S 5S
    # chr1    0   0    0 80
    # chr13  14  13   13  0   ! TO REMOVE
    # chr14  11  10   11  0   ! 2 TO REMOVE
    # chr15  17  16   16  0   ! TO REMOVE
    # chr2    0   0    0  2
    # chr20   0   0    3  0
    # chr21   8   7    7  0   ! TO REMOVE
    # chr22  53  52   52  0   ! TO REMOVE
    # chr6    0   0    0  1
    ############----------------------  REMOVE 18S & EXTRA 5.8S ---------------------#############
    ## We remove extra 18 and une 5.8 S:
    for (chr in c("13","14", "15", "21", "22")){
      s18 <- d$start[ d$chr==chr & d$gene =="18S"]
      s58 <- d$start[ d$chr==chr & d$gene =="5.8S"]
      s28 <- d$start[ d$chr==chr & d$gene =="28S"]
      
      if (length(s18)!=length(s28)){
        # deleted <-  rbind(deleted,d[( d$start == max(s18) ),] )
        d <-  d[!( d$start == max(s18) ),]
      }
      if(length(s58)!=length(s28)){ # chr 14 also has an extra 5.8S
        # deleted <-  rbind(deleted,d[( d$start == max(s58) ),] )
        d <-  d[!( d$start == max(s58) ),]
      }
    }
    
    # We want to remove also the 3 extra 5.8S:
    d <- d[ !d$chr == "20", ]
    # Do we want to remove the 5S copies from no chr 1 regions?
    d <- d[ !(d$chr != "1" & d$gene =="5S"), ]
    
    ############----------------------------  Plotting Chr  ---------------------#############
    df <- as.data.frame(table(d$chr, d$gene))
    
    missing <- expand.grid(Var1 = c(1:22, "X", "Y"), Var2 = c("S5", "S18", "S28","S5.8" ), Freq = 0)
    missing$Var1 <- factor(missing$Var1, levels = c(1:22, "X", "Y"))
    
    df <-  bind_rows(df, missing)
    df <-  distinct(df, Var1, Var2, .keep_all = TRUE)
    
    png(file=paste0(dirI,"chr_continous_operons_stats_", name,".png"), res=220, width=2300, height=1300);
    print(get_plot_genes(df, title = paste("Selecting continous operons", from)))
    dev.off();
    
    ############---------------------------  Other Plotting ---------------------#############
    # df$file <- name
    # dfs <-  rbind(dfs, df)
    # 
    df <- df[df$Freq !=0, ]
    
    nt <-  data.frame(type = unique(df$Var2), counts = NA)
    nt$counts <-aggregate(Freq ~ Var2, data = df, FUN = sum)$Freq
    
    png(file=paste0(dirI,"genes_stats_",  name,".png"), res=320, width=1300, height=2300);
    print( plot_nt(table = nt, title="Found rRNA terms per product type", subtitle=from, angle =0 , size = 5))
    dev.off();
    
########################## 2 Transform into continuous operon  #########################
    ############---------------- Merge rows continous 18S 5.8S 28S -------------------#############
    operons <-  d[d$gene !="5S",]
    d <-  d[d$gene =="5S",]
    
      if ((nrow(operons)  %% 3)==0 ){
        for (i in seq(1, nrow(operons), by=3)){
          
          row <- operons[i,]
          row$end <-  operons$end[i+2]
          row$gene <-  "45S"
          
          d <-  rbind(d, row)
        }
      }else{
        print("There are more sparsed 45S copies not within an operon")
      }
    # Get bed file with the genes to extract fasta files
    w<-  d
    w$chr <- paste0("chr", w$chr,  "_", gsub("chr\\w+_(\\w+)", "\\1", IDs[1])) 
    write.table(w, paste0(dirF,"YAO_rRNA_45S_5S_", name,".bed"), row.names = F,col.names = F, sep="\t", quote=FALSE)    
  
    # # Get  bed files with matching chr names
    # w<-  d
    # chr <-  paste0(get_ID(name), "0", w$chr[as.numeric(w$chr)<10])
    # chr <- c(chr,  paste0(get_ID(name), w$chr[as.numeric(w$chr)>10]))
    # w$chr<-  chr
    # write.table(w, paste0(dirF,"YAO_rRNA_45S_5S_", name,".bed"), row.names = F,col.names = F, sep="\t", quote=FALSE)    
    
    
    ############---------------- Plot expected vs found values -------------------#############
    operons <-  d[d$gene !="5S",]
    t <-   as.data.frame( table(operons$chr))
    
    df <-  data.frame(name = name, 
                      expect = get_CN(name), 
                      found = dim(operons)[1],
                      chr  =  t$Var1,
                      counts = t$Freq)
    
    CN <- rbind(CN, df)
    
    # # Save only 45S: 13 14 15 21 and 22:
    # w <- operons
    # w$chr <-  paste0("chr", w$chr,  "_", gsub("chr\\w+_(\\w+)", "\\1", IDs[1])) 
    # write.table(w, paste0(dirF,"YAO_rRNA_45S_", name,".bed"), row.names = F,col.names = F, sep="\t", quote=FALSE)    

    # Get  bed files with matching chr names
    w <- operons
    chr <- paste0(get_ID(name), w$chr[as.numeric(w$chr)>10])
    w$chr<-  chr
    # Save only 45S: 13 14 15 21 and 22:
    write.table(w, paste0(dirF,"YAO_rRNA_45S_", name,".bed"), row.names = F,col.names = F, sep="\t", quote=FALSE)    
    
    
 }

########------------------------- Plotting expected vs got ------------------------------#########

expected <- unique(CN$expect)

png(file=paste0(dirI,"continous_operons_results.png"), res=500, width=2000, height=2000);
ggplot(CN, aes(x = name, y = counts, fill = factor(chr))) +
  geom_segment(aes(x="hp",  y = expected[1],  xend = "pat", yend = expected[1]+0.5), linetype = "dashed", color = "red") +
  geom_segment(aes(x=2.5, y = expected[2],    xend = 3.5,     yend = expected[2]+0.5), linetype = "dashed", color = "red") +
  geom_segment(aes(x=0.5,     y = expected[3],xend =1.5, yend = expected[3]+0.5), linetype = "dashed", color = "red") +
  annotate("text", x = c("mat", "pat", "hp"), y = expected, label = paste(expected, "Reported"),
           color = "red", size = 3.8, vjust = -1) +
  geom_bar(position = "stack", stat = "identity") +
  geom_text(aes(label = counts), color = "black", size = 4, position = position_stack(0.5)) +
  scale_fill_manual(values = palette1) +
  labs(title = "",x = "", y="Number of 45S sequences", fill = "Chromosome") +
  scale_x_discrete(labels=c("hp"="Offspring", "mat"="Maternal", "pat"="Paternal"))+
  coord_cartesian(xlim =c(1, 3), ylim = c(0, 155))+
  theme_minimal()+
  theme(axis.text.x = element_text(size=11), 
        axis.text.y = element_text(size=11),
        axis.title= element_text(size=10),
        legend.position = "bottom", 
        legend.title = element_text(size=10),
        legend.margin=margin(-10,-10,0,-15))
dev.off();

#########  IMPORTANT PLOT ########

png(file=paste0("~/Desktop/ribosomal_RNAs/Alba/IMAGES TFG/","continous_operons_results.png"), res=500, width=1750, height=2000);
ggplot(CN, aes(x = name, y = counts, fill = factor(chr))) +
  geom_segment(aes(x="hp",  y = expected[1],  xend = "pat", yend = expected[1]+0.5), linetype = "dashed", color = "#ff4b92") +
  geom_segment(aes(x=2.5, y = expected[2],    xend = 3.5,     yend = expected[2]+0.5), linetype = "dashed", color = "#ff4b92") +
  geom_segment(aes(x=0.5,     y = expected[3],xend =1.5, yend = expected[3]+0.5), linetype = "dashed", color = "#ff4b92") +
  annotate("text", x = c("mat", "pat", "hp"), y = expected, label = paste(expected, "Reported"),
           color = "#ff0065", size = 3.8, vjust = -1) +
  geom_bar(position = "stack", stat = "identity") +
  geom_text(aes(label = counts), color = "black", size = 4, position = position_stack(0.5)) +
  scale_fill_manual(values = palette1) +
  labs(title = "",x = "", y="number of 45S sequences", fill = "Chromosome") +
  scale_x_discrete(labels=c("hp"="F1", "mat"="Maternal", "pat"="Paternal"))+
  coord_cartesian(xlim =c(1, 3), ylim = c(7, 155))+
  scale_y_continuous(breaks = seq(0, 150, by = 50))+
  theme_minimal()+
  theme(axis.text.x = element_text(size=10), 
        axis.text.y  = element_text(size=7), 
        axis.title= element_text(size=10),
        legend.position = "bottom", 
        legend.title = element_text(size=9),
        legend.margin=margin(-10,-10,0,-25),
        legend.key.size = unit(3.5, "mm" ),
        axis.line = element_line(color = "black", size=0.3))
dev.off();

