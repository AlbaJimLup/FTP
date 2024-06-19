library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)
####### GET DATA  ######
dir = "~/Desktop/ribosomal_RNAs/Alba/07_Performance/"
dirI<- paste0(dir, "IMAGES/")

# Get all tables
s<-  read.table(paste0(dir, "Performance_sim_errors_vc_filtered_sample.tbl"),  sep=",",header = T)
s<- s[s$sample %in% c("2samples", "15samples"),]
#### --------------------------------------  LOAD DATA ---------------------------------############
raw<- data.frame()
for (v in 2:4){
  vc<-  read.table( paste0(dir,"Performance_sim_errors_vc_filtered_v", v, ".tbl"), sep=",")
  raw<- rbind(raw, vc[vc$V1 %in% c(50, 100), ])
}
raw<- raw[, 1:4]
raw$ploidy=raw$V1
raw$sample="No filter"
raw$qual=c("QUAL30","QUAL30", "QUAL400", "QUAL400", "QUAL1000","QUAL1000")
raw$type<- "vc"
names(raw)<- names(s)
raw$set<- paste(raw$ploidy, raw$qual)

s<- rbind(s, raw)
s$precision<- as.numeric(s$precision)
s$recall<- as.numeric(s$recall)
s$F1_score<- as.numeric(s$F1_score)

best<- data.frame()
for (p in unique(s$ploidy)){
  
  best<- rbind(best, s[s$ploidy==p & s$F1_score==max(s$F1_score[s$sample=="2samples" & s$ploidy==p]),])
}
write.table(best, paste0(dir,  "Performance_sim_errors_vc_sample_best.tbl"), sep=",", quote=F, row.names=F, col.names=T)

best<- data.frame()
for (q in unique(s$qual)){
  
  best<- rbind(best, s[s$qual==q & s$F1_score==max(s$F1_score[s$qual==q]),])
}


#####################   Palettes    #########################
palette2<- c("QUAL1000"="#588571", "QUAL400"="#61A9ED", "QUAL30"="#ffaa2c", "50"="black", "100"="lightgray")
# palette2<- c("No filter"="#588571", "2samples"="#61A9ED", "10samples"="#ff006a", "15samples"="#ff006a","50"="black", "100"="lightgray")
##############   PLOTTING precision recall  ########################
png(file=paste0(dirI, "Final_Sample_Precision_vs_Recall_p50&100.png"), res=470, width=1300, height=800);
ggplot(s, aes(x = recall, y = precision, shape = factor(sample, levels =c("No filter","2samples", "15samples")))) +
  geom_point(aes(color = factor(qual, levels = c("QUAL30", "QUAL400", "QUAL1000"))), size=4, alpha=0.7) +
  geom_point(aes(color = factor(ploidy)), shape=3, size=3, alpha=0.6) +
  labs(x = "Recall", y = "Precision", color = "", shape = "") +
  scale_shape_manual( values= c( 0,18, 20))+
  scale_color_manual(values=palette2)+
  theme_minimal()+
  guides(color = guide_legend(override.aes = list(size=2))) +   
  guides(shape = guide_legend(override.aes = list(size=1.5)))+
  # coord_cartesian(xlim =c(0.9, 0.93), ylim = c(0.85, 0.95))+
  theme(legend.key.size = unit(0, 'lines'), 
        legend.key.height = unit(3, 'mm'), #change legend key height
        legend.key.width = unit(0, 'cm'),
        legend.spacing.x = unit(1, 'mm'),
        legend.text = element_text(size=6, margin = margin(l = 0.7, unit = "mm")), 
        legend.title = element_blank(), 
        axis.text= element_text(size=6),
        axis.title = element_text(size=7),
        axis.line = element_line(color = "black", size=0.2),
        legend.margin=margin(-5,-2,0,0),
        legend.box.margin=margin(5,-2,0,0))
dev.off();



##############   PLOTTING F1 score  ########################
palette <- c(  "50 QUAL30"="#e5729b", "50 QUAL400"="#ce3f7a", "50 QUAL1000"="#9e0142",
               "50 2samples QUAL30"="#FF9999",  "50 2samples QUAL400"="#FF6666",  "50 2samples QUAL1000"="#FF3333",
               "50 15samples QUAL30"="#ffc185", "50 15samples QUAL400"="#ffaa2c", "50 15samples QUAL1000"="#f09000",
               
               "100 QUAL30"="#a9a3c5", "100 QUAL400"="#6658a4", "100 QUAL1000"="#261477",
               "100 2samples QUAL30"="#BaCFFF",  "100 2samples QUAL400"="#6666FF",  "100 2samples QUAL1000"="#3333FF",
               "100 15samples QUAL30"="#B3E6E6", "100 15samples QUAL400"="#66CCCC", "100 15samples QUAL1000"="#339999")

levels(s$set)<- names(palette)
png(file=paste0(dirI, "Final_Sample_F1_score_callable_vc_filters_label_p50&100_v2.png"), res=180, width=1100, height=200);
ggplot(s, aes(x=F1_score, y = factor(set, levels = names(palette)) , fill = factor(set))) +
  coord_cartesian(xlim=c(0.75, 0.99))+
  geom_bar(stat="identity")+
  geom_text(aes(label=round(F1_score, 4)), hjust=1.4,  size=3, color="white", fontface="bold")+
  scale_fill_manual(values=palette)+
  guides(fill="none")+
  theme_minimal()+
  theme( axis.text.x = element_text(size=7,vjust=1, hjust=1),  axis.title.x = element_text(size=8),
         axis.title.y = element_blank(), axis.text.y = element_text(size=7))
dev.off();

# Get all tables
s<-  read.table(paste0(dir, "Performance_sim_errors_vc_filtered_sample.tbl"),  sep=",",header = T)
s<- s[s$sample %in% c("2samples", "5samples", "10samples", "15samples"),]
#### --------------------------------------  LOAD DATA ---------------------------------############
raw<- data.frame()
for (v in 2:4){
  vc<-  read.table( paste0(dir,"Performance_sim_errors_vc_filtered_v", v, ".tbl"), sep=",")
  raw<- rbind(raw, vc[vc$V1 %in% c(50, 100), ])
}
raw<- raw[, 1:4]
raw$ploidy=raw$V1
raw$sample="No filter"
raw$qual=c("QUAL30","QUAL30", "QUAL400", "QUAL400", "QUAL1000","QUAL1000")
raw$type<- "vc"
names(raw)<- names(s)
raw$set<- paste(raw$ploidy, raw$qual)

s<- rbind(s, raw)
s$precision<- as.numeric(s$precision)
s$recall<- as.numeric(s$recall)
s$F1_score<- as.numeric(s$F1_score)

s<- s[s$qual=="QUAL1000", ]

palette <- c(  "50 15samples QUAL1000" = "#710000",  "50 10samples QUAL1000" = "#A90000",
  "50 5samples QUAL1000" = "#d73027",  "50 2samples QUAL1000" = "#EA7000",
  "50 QUAL1000" = "#f09000",  "100 15samples QUAL1000" = "#00045A",
  "100 10samples QUAL1000" = "#313695",  "100 5samples QUAL1000" = "#4C6EA2",
  "100 2samples QUAL1000" = "#598AA8",  "100 QUAL1000" = "#66A6AE")


levels(s$set)<- names(palette)
png(file=paste0(dirI, "Final_Sample_F1_score_callable_vc_filters_label_QUAL1000.png"), res=180, width=1100, height=400);
ggplot(s, aes(x=F1_score, y = factor(set, levels = names(palette)) , fill = factor(set))) +
  coord_cartesian(xlim=c(0.75, 0.95))+
  geom_bar(stat="identity")+
  geom_text(aes(label=round(F1_score, 4)), hjust=1.4,  size=3, color="white", fontface="bold")+
  scale_fill_manual(values=palette)+
  guides(fill="none")+
  theme_minimal()+
  theme( axis.text.x = element_text(size=7,vjust=1, hjust=1),  axis.title.x = element_text(size=8),
         axis.title.y = element_blank(), axis.text.y = element_text(size=7))
dev.off();