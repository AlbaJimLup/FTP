# setwd("~/Documents/Practicals docs/Procedures/simulated_reads")
dir <-  "~/Desktop/ribosomal_RNAs/Alba/02_simulated_reads/v1/"
t <-  read.table(paste0(dir, "pipeline_reads_counts.tbl"), header = T)
# t <-  data.frame( sample=s$sample, 
#                   coverage = s$coverage,
#                   rDNA_like_reads = s$num_seqs_fq, #  initial number beggining pipeline
#                   rRNA_reads = s$RNA_reads, # kept after mummer step
#                   rDNA_not_rRNA = s$RNA_like, # discarted in mummer step
#                   uniquely_mapped_reads = s$unikely_mapped, # kept after BWA step
#                   rRNA_not_uniq = s$RNA ) # discarted after BWA step 

for (i in 1:nrow(t)){
  t$rDNA_not_pseudogenes[i] <- t$rRNA_reads[i]  / t$rDNA_like_reads[i] # fraction kept in mummer step "pseudogenes"
  t$uniquely_mapped_rDNA[i] <- t$uniquely_mapped_reads[i] /t$rRNA_reads[i] # fraction kept in BWA step "multimapped"
}


library(tidyr)
t_long <- pivot_longer(t, cols = c(rDNA_not_pseudogenes, uniquely_mapped_rDNA),
                       names_to = "RNA_Reads", values_to = "freq")
t_long$percent <- paste0(as.character(round(t_long$freq*100, 1)), "%")
# write.table(t_long, "simulated_reads_stats_long.tbl", sep = "\t", quote = FALSE, row.names = FALSE)

t_long$RNA_Reads[t_long$RNA_Reads=="rDNA_not_pseudogenes"]<- "Not pseudogenes"
t_long$RNA_Reads[t_long$RNA_Reads=="uniquely_mapped_rDNA"]<- "Uniquely mapped"
#### PLOTTING#########
library(ggplot2)
library(scales)

palette = c("Not pseudogenes"="#bd0142","Uniquely mapped"="#66c865" )

png(file=paste0(dir,"IMAGES/MUMmer_scatter_percent.png"), res=150, width=3000, height=2000);
ggplot(t_long, aes(x= factor(sample, levels =  paste0("Sample_", 1:100)), y = freq, color = RNA_Reads, label = percent) )+
  geom_point(size=2)+
  scale_color_manual(values = palette) +
  geom_text( size=2.5,color="black",  nudge_y = 0.009)+
  labs(title="", x="", y="Number of sequences", color="rDNA Reads")+
  scale_y_continuous(limit =c(0.68, 1.005), expand = c(0, 0)) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90))+
  coord_flip()
dev.off();



palette = c("#bd0142","#66c865" )
t_long$RNA_Reads[t_long$RNA_Reads=="Not pseudogenes"]<- "rDNA-like \n not pseudogenes "
t_long$RNA_Reads[t_long$RNA_Reads=="Uniquely mapped"]<- "rDNA \n   uniquely mapped"


t_long$RNA_Reads<- factor(t_long$RNA_Reads , levels=c("rDNA-like \n not pseudogenes " , "rDNA \n   uniquely mapped"))

png(file=paste0(dir,"IMAGES/MUMmer_boxplot_percent.png"), res=550, width=1150, height=1500);
ggplot(t_long, aes( x=RNA_Reads, y= freq,color = RNA_Reads) )+
  geom_jitter(size=1.2,  width=0.5, alpha =0.6 )+
  geom_boxplot(alpha =0.3, color="black", width=0.5, outlier.shape = NA)+
  scale_color_manual(values = palette) +
  scale_y_continuous(limit =c(0.68, 1.005), expand = c(0, 0)) +
  labs(title="", x="", y="Fraction of reads kept per sample", color="")+
  guides(color="none")+
  theme_minimal()+
  theme(legend.position="bottom", 
        axis.ticks.x = element_blank(),
        axis.title= element_text(size=8), 
        axis.text.y = element_text(size=7), 
        axis.text.x = element_text(size=6.5, color="black"), 
        axis.line = element_line(color = "black", size=0.2))
dev.off();


png(file=paste0("~/Desktop/ribosomal_RNAs/Alba/IMAGES TFG/MUMmer_boxplot_percent.png"), res=550, width=1300, height=1500);
ggplot(t_long, aes( x=RNA_Reads, y= freq,color = RNA_Reads) )+
  geom_jitter(size=1.2,  width=0.5, alpha =0.6 )+
  geom_boxplot(alpha =0.3, color="black", width=0.5, outlier.shape = NA)+
  scale_color_manual(values = palette) +
  scale_y_continuous(limit =c(0.68, 1.005), expand = c(0, 0)) +
  labs(title="", x="", y="Fraction of reads kept per sample", color="")+
  guides(color="none")+
  theme_minimal()+
  theme(legend.position="bottom", 
        axis.ticks.x = element_blank(),
        axis.title= element_text(size=8), 
        axis.text.y = element_text(size=7), 
        axis.text.x = element_text(size=6.5, color="black"), 
        axis.line = element_line(color = "black", size=0.2))
dev.off();

# ggplot(t_long, aes( x=factor(RNA_Reads, levels=c( "rDNA_not_pseudogenes", "uniquely_mapped_rDNA")) , y= freq,color = RNA_Reads) )+
#   geom_jitter(size=1.5, width=0.4,alpha =0.65)+
#   geom_boxplot(alpha =0.2, color="black", width=0.5, outlier.shape = NA)+
#   scale_color_manual(values = palette, labels = c("Not pseudogenes","Uniquely mapped")) +

#   labs(x="rDNA reads", y="Fraction", color="")+#   theme_minimal()+
#   coord_cartesian(ylim = c(0.7,1))+
#   guides(color = guide_legend(override.aes = list(shape = NA)))+
#   theme(legend.position="bottom", 
#         axis.text.x = element_blank(), 
#         axis.ticks.x = element_blank(),
#         axis.title= element_text(size=8),
#         axis.title.x = element_text(size=8, margin = margin(t = 15)),
#         axis.text.y = element_text(size=7.5), 
#         legend.text = element_text(size=7,  margin = margin(r=0, unit = "pt")),
#         legend.box.spacing = unit(0, "pt"),
#         legend.margin=margin(0,-11,0,0),
#         legend.box.margin=margin(-23,5,10,-10),
#         legend.key.size = unit(2, "mm" ),
#         axis.line = element_line(color = "black", size=0.2))

