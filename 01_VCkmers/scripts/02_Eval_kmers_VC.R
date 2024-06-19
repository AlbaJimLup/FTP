###### GET DATA #######
dir="~/Desktop/ribosomal_RNAs/Alba/01_VCkmers/"
#
s30 <- read.table(paste0(dir, "x30stats_evaluationGATKvsMSA.tbl"), 
                    header = TRUE, row.names = NULL, sep = ",")
s100 <- read.table(paste0(dir, "x100stats_evaluationGATKvsMSA.tbl"), 
                    header = TRUE, row.names = NULL, sep = ",")
s1000 <- read.table(paste0(dir, "x1000stats_evaluationGATKvsMSA.tbl"), 
                    header = TRUE, row.names = NULL, sep = ",")
library(forcats)
library(ggplot2)
library(dplyr)

#### PREP #####
s30$coverage <- rep("30x", nrow(s30))
s100$coverage <- rep("100x", nrow(s100))
s1000$coverage <- rep("1000x", nrow(s1000))


counts_ev <-  rbind(s30, s100, s1000)

counts_ev$evaluation[counts_ev$evaluation =="same variant"] <-  "TP"
counts_ev$evaluation[counts_ev$evaluation =="missing variant"] <- "FN"
counts_ev$evaluation[counts_ev$evaluation =="missing allele"] <-"FN allele" 
counts_ev$evaluation[counts_ev$evaluation =="different alleles"] <- "FP allele"

counts_ev <- counts_ev %>% arrange(counts) 

###### COLOR ########

counts_ev$evaluation <- factor(counts_ev$evaluation, 
                               levels = c( "FP allele", "FP", "FN allele","FN", "TP"))


counts_ev$class[counts_ev$class=="callable"]<- "Non-repetitive regions"
counts_ev$class[counts_ev$class=="blacklisted"]<- "Repetitive regions"


palette <- c("TP" ="#6bac65", "TP_missing_allele"="#7fc3a1",  
              "FN"="#ffc60e","FN allele"="#ecd189",
              "FP"= "#bd0752","FP_sample"="#ff4f8d", "FP allele"="#ed8c96","err"= "#a6a6a6") 

###### PLOT ########
p <- ggplot( counts_ev, aes(x = factor(coverage, levels = c("30x", "100x", "1000x")),  y=freqs,    fill = factor(evaluation))) +
  geom_bar( stat = "identity", position = "stack")+
  geom_text(aes(label = counts),
            position = position_stack(vjust = 0.3), size = 3.3, fontface = "bold", color="black") +
  theme_classic()+
  labs(y="Frequency", x="Set of kmers")+
  scale_fill_manual(values=palette, name = "")+
  labs(fill ="")+
  facet_grid(.~factor(class))+
  coord_cartesian(ylim=c(0, 1))+
  theme(legend.text= element_text(size=9), 
        axis.title = element_text(size=9),
        axis.text.y = element_text(size=9),
        legend.position = "bottom",
        legend.margin=margin(-5,-10,0,0),
        legend.key.size = unit(4, "mm" )); p

png(file=paste0("~/Desktop/ribosomal_RNAs/Alba/IMAGES TFG/kmers_GATKvsMSAalleles.png"), 
    res=400, width=1750, height=1300);
p
dev.off()

png(file=paste0(dir, "IMAGES/kmers_GATKvsMSAalleles.png"), res=400, width=1600, height=1300);
p
dev.off()