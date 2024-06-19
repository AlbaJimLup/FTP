# libraries ----
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(cowplot)
library(reshape2)
library(ggrepel)

# granges-> collapse 


################# GET AND PREP FILES #################
dir="~/Desktop/ribosomal_RNAs/Alba/"
dirI=paste0(dir, "05_rDNA/IMAGES/")

chrR<- read.table(paste0(dir, "data/chrR.bed"))
chrR<-  chrR[c(2,3,4)]
names(chrR)<- c("start", "end",  "region")
chrR$length<- chrR$end-(chrR$start)
# bed<- data.frame()
# for (region in c("5_ETS","18S","ITS1","5.8S","ITS2","28S","3_ETS")){
#   t <- read.table(paste0(dir, "data/pre-rRNA_47S.included_intervals.", region, ".bed"))
# 
#   bed<- rbind(bed, t)
# }
bed<-read.table(paste0(dir, "data/chrR_overlap_with_callable_regions.tsv"), header = T)
# bed<-read.table(paste0(dir, "data/pre-rRNA_47S.included_intervals_per_region.bed"))
# names(bed)<- c("chr", "start", "end",  "region")
# bed$region[bed$region=="5'_ETS"]<- "5_ETS"
# # bed$region[bed$region=="3'_ETS"]<- "3_ETS"
# bed$start<- bed$start +1
# bed$length=bed$end - bed$start # window of callable length

################# GET LENGTHS CALLABLE AND BLACKLISTED #################
t<- data.frame(matrix(ncol=4, nrow=7))
names(t)<- c("region", "length","callable", "blacklisted")
t$region<- c("5_ETS","18S","ITS1","5.8S","ITS2","28S","3_ETS")
# t$length<- c(12995-9338, # "5_ETS"
#              14864-12995, # "18S"
#              15934-14864, # ITS1
#              16091-15934, #5.8S
#              17258-16091, # ITS2
#              22309-17258, # 28S
#              22670-22309) # 3_ETS

for (region in t$region){
  t$length[t$region==region]<- chrR$length[chrR$region==region] # total lengths region
  t$callable[t$region==region]<- sum(bed$Overlap[bed$Name==region]) #callable lengths
}
# length blacklisted
t$blacklisted<- t$length-t$callable

################# FREQUENCY #################
t$freq_c<-rep(0, nrow(t))
t$freq_b<- t$freq_c
for (i in 1:nrow(t)){
  t$freq_c[i]<- t$callable[i]/t$length[i]
  t$freq_b[i]<- 1-t$freq_c[i]
}
################# PREP PLOTTING #################
t$region[t$region=="5_ETS"]<- "5'ETS"
t$region[t$region=="3_ETS"]<- "3'ETS"
levels(t$region)<- c("5'ETS","18S","ITS1","5.8S","ITS2","28S","3'ETS")
names(t)<- c("region", "lengths", "Non repetitive regions", "Repetitive regions", "freq_c", "freq_b")

t_long <- pivot_longer(t, cols = c( "Non repetitive regions", "Repetitive regions"), names_to = "class", values_to = "counts")
freq_long <- pivot_longer(t, cols = c( "freq_c", "freq_b"), names_to = "class", values_to = "counts")


t_long$freq<- freq_long$counts
levels(t_long$region)<- c("5'ETS","18S","ITS1","5.8S","ITS2","28S","3'ETS")
levels(t_long$class)<- c("Non repetitive regions", "Repetitive regions")
palette<- c("Non repetitive regions"="#81b5e7","Repetitive regions"="#aabbc6")
#
#################  PLOT LONG FORMAT  ##################
png(file=paste0(dirI, "callable_vs_blacklisted_region.png"), res=600, width=2500, height=2000);
ggplot(t_long, aes(x=factor(region, levels=c("5'ETS","18S","ITS1","5.8S","ITS2","28S","3'ETS")), y= counts, fill=class, label=counts))+
  geom_bar(stat="identity")+
  geom_text_repel(position = position_stack(vjust = 0.5), force = 0.001, fontface = "bold", size=3) +
  theme_minimal()+
  labs(y="Number of positions", x="Region", fill="Masked regions")+
  scale_fill_manual(values=palette)+
  theme(axis.text.x = element_text(color="black"),legend.margin=margin(-15,-5,-5,-15),
        legend.text = element_text(size = 8.5),legend.key.size = unit(3, "mm") )
dev.off();


png(file=paste0(dirI, "callable_vs_blacklisted_norm_region.png"), res=600, width=2500, height=2000);
ggplot(t_long, aes(x=factor(region, levels=c("5'ETS","18S","ITS1","5.8S","ITS2","28S","3'ETS")), y= freq, fill=class, label=counts))+
  geom_bar(stat="identity")+
  geom_text_repel(position = position_stack(vjust = 0.5), force = 0.001, fontface = "bold", size=3) +
  theme_minimal()+
  labs(y="Number of positions", x="Region", fill="Masked regions")+
  scale_fill_manual(values=palette)+
  theme(axis.text.x = element_text(color="black"),legend.margin=margin(-15,-5,-5,-10),
        legend.text = element_text(size = 8.5),legend.key.size = unit(3, "mm") )
dev.off();

############ Do total percentage ##########

t$lengths<- as.numeric(t$lengths)
t$`Repetitive regions`<- as.numeric(t$`Repetitive regions`)
t$`Non repetitive regions` <- as.numeric(t$`Non repetitive regions`)


d<- data.frame(matrix(ncol = ncol(t), nrow=1))
names(d)<- names(t)
d[1,]<- c("Average", sum(t$lengths), 
               sum(t$`Non repetitive regions`),
               sum(t$`Repetitive regions`), 
               sum(t$`Non repetitive regions`)/ sum(t$lengths), 
               sum(t$`Repetitive regions`)/ sum(t$lengths))

d_long <- pivot_longer(d, cols = c( "Non repetitive regions", "Repetitive regions"), names_to = "class", values_to = "counts")
freq_long <- pivot_longer(d, cols = c( "freq_c", "freq_b"), names_to = "class", values_to = "counts")

d_long$freq<- freq_long$counts
levels(d_long$region)<- c("Average")

png(file=paste0(dirI, "Average_callable_vs_blacklisted_region.png"), res=800, width=2000, height=3000);
ggplot(d_long, aes(x=factor(region, levels=c("Average")), y= as.numeric(freq), fill=class, label=paste(round(as.numeric(freq)*100,2), "%")))+
  geom_bar(stat="identity")+
  geom_text_repel(position = position_stack(vjust = 0.5), force = 0.001, fontface = "bold", size=3) +
  theme_minimal()+
  labs(y="Number of positions", x="Region", fill="Masked regions")+
  scale_fill_manual(values=palette)+
  theme(axis.text.x = element_text(color="black"),legend.margin=margin(-15,-5,-5,-15),
        legend.text = element_text(size = 8.5),legend.key.size = unit(3, "mm") )
dev.off();



##############################

## !! plot main
d_long$counts<- paste0(round(as.numeric(freq_long$counts)*100,2), "%")
levels(t_long$region)<- c("5'ETS","18S","ITS1","5.8S","ITS2","28S","3'ETS", "Average")
all<- rbind(t_long, d_long)

all$level<- c(rep("ChrR region", 14), rep("Total", 2))

png(file=paste0(dirI, "all_callable_vs_blacklisted.png"), res=700, width=3400, height=2000);
ggplot(all, aes(x=factor(region, levels=c("5'ETS","18S","ITS1","5.8S","ITS2","28S","3'ETS", "Average")),
                y=  as.numeric(freq), fill=reorder(class, as.numeric(freq)), label=counts))+
  geom_bar(stat="identity")+
  geom_text(position = position_stack(vjust = 0.7),  fontface = "bold", size=2.8) +
  theme_classic()+
  facet_grid(.~level, space="free",  scales="free_x" )+
  # guides(fill="none")+
  labs(y="Fraction", x="", fill="Masked regions")+
  scale_fill_manual(values=palette)+
  theme(axis.text.x = element_text(color="black"),legend.margin=margin(-15,0,-5,-10),
        legend.title = element_text(size=10),
        legend.text = element_text(size = 7, margin = margin(r = -5)),
        legend.key.size = unit(3, "mm"), 
        panel.grid.major = element_line(colour = "#e6e6e6"),
        panel.grid.minor = element_line(colour = "#e6e6e6") )
dev.off();




t_long$counts<- paste0(round(as.numeric(t_long$freq)*100,2), "%")
d_long$counts<- paste0(round(as.numeric(d_long$freq)*100,2), "%")

levels(t_long$region)<- c("5'ETS","18S","ITS1","5.8S","ITS2","28S","3'ETS", "Average")
all<- rbind(t_long, d_long)

all$level<- c(rep("ChrR region", 14), rep("Total", 2))

png(file=paste0(dirI, "all_callable_vs_blacklisted_percent.png"), res=700, width=3000, height=2000);
ggplot(all, aes(x=factor(region, levels=c("5'ETS","18S","ITS1","5.8S","ITS2","28S","3'ETS", "Average")),
                y=  as.numeric(freq), fill=class, label=counts))+
  geom_bar(stat="identity")+
  geom_text(position = position_stack(vjust = 0.7),  fontface = "bold", size=2.5) +
  theme_classic()+
  facet_grid(.~level, space="free",  scales="free_x" )+
  guides(fill="none")+
  labs(y="Number of positions", x="Region", fill="Masked regions")+
  scale_fill_manual(values=palette)+
  theme(axis.text.x = element_text(color="black"),legend.margin=margin(-15,-5,-5,-10),
        legend.text = element_text(size = 8.5),legend.key.size = unit(3, "mm"),   
        panel.grid.major = element_line(colour = "#e6e6e6"),
        panel.grid.minor = element_line(colour = "#e6e6e6") )
dev.off();
