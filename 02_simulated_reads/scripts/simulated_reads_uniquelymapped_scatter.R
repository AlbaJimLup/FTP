# setwd("~/Documents/Practicals docs/Procedures/simulated_reads")
dir="~/Desktop/ribosomal_RNAs/Alba/02_simulated_reads/v1/"
t <-  read.table(paste0(dir, "pipeline_reads_counts.tbl"), header = T)
# write.table(CN, "~/Documents/Practicals\ docs/Procedures/simulated_reads/CN.tbl", sep = "\t", quote = FALSE, row.names = FALSE)
CN <- read.table(paste0(dir, "CN.tbl"), header=T)
#### PLOTTING#########
library(ggplot2)
library(scales)


t$CN <-CN$number_of_copies

palette = setNames( c("#677ad1","#66c865" ), c( "freqrRNA_reads", "freq_uniq_mapp"))


png(file=paste0(dir,"IMAGES/STATS_numseqs_MUMmer_coverage.png"), res=400, width=3200, height=2000);
ggplot(t, aes(x= coverage, y = rDNA_like_reads, size=CN))+
  geom_point(alpha=0.6, color="darkblue" )+
  theme_minimal()+
  labs(title ="NEAT",y= "Number of rDNA like reads", x ="Coverage")+
  scale_y_continuous(labels = comma) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 17),
        title = element_blank(),
        legend.title = element_text(size = 18),
        legend.text = element_text(size=15))
dev.off();

  

png(file=paste0(dir,"IMAGES/STATS_numseqs_MUMmer_coverage_flip.png"), res=400, width=3200, height=2000);
ggplot(t, aes(x= coverage, y = rDNA_like_reads, size=CN))+
  geom_point(alpha=0.6, color="darkblue" )+
  theme_minimal()+
  labs(title ="NEAT",y= "Number of rDNA like reads", x ="Coverage")+
  scale_y_continuous(labels = comma) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 17),
        title = element_blank(),
        legend.title = element_text(size = 18),
        legend.text = element_text(size=15))
  coord_flip()
dev.off();
