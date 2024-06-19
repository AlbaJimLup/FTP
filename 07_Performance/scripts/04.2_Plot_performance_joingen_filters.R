library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)
####### GET DATA  ######
dir = "~/Desktop/ribosomal_RNAs/Alba/07_Performance/"
dirI<- paste0(dir, "IMAGES/")
# Get all tables
s1<-  read.table(paste0(dir, "Performance_sim_errors_postjoin.tbl"), sep=",", header = T)
s2<-  read.table(paste0(dir, "Performance_sim_errors_jgeontype_filtered_v1.tbl"),sep=",", header = T)
s3<-  read.table(paste0(dir, "Performance_sim_errors_jgeontype_filtered_v2.tbl"), sep=",",header = T)
s4<-  read.table(paste0(dir, "Performance_sim_errors_jgeontype_filtered_v3.tbl"), sep=",",header = T)


s<-  rbind(s1,s2,s3, s4)
# s$eval_level <- c( rep("NoFilter", 6), rep("Filter1", 6),  rep("Filter2", 6), rep("Filter3", 6))
s$eval_level <-factor(s$eval_level, levels= c("NoFilter","Filter1", "Filter2", "Filter3" ))

##############    Get best ########################
s$file<- factor(paste0(s$eval_level, "_p", s$set), levels = paste0(s$eval_level, "_p", s$set))

best<- data.frame()

for (ploidy in unique(s$set)){
  # keep filter which best performance
  best<- rbind(best, s[s$F1_score == max(s$F1_score[s$set == ploidy]),])
}

best<- best[!best$eval_level=="Filter1",]

write.table(best, paste0(dir,  "Performance_sim_errors_joingen_filtered_best.tbl"), sep=",", quote=F, row.names=F, col.names=T)

####################### Palette ###################################
# palette <- c("NoFilter" ="#9e2142","Filter1" ="#f8b95d", "Filter2" ="#50917c", "Filter3"  ="#61a9ed" )
palette1<- c( "Filter1"="#9a1142", "Filter2"="#ffaa2c","Filter3"="royalblue", "Filter4"="#377459")

palette<- c("#ff994a", "#9a1142", "#e767c5", "#7066bc", "#93c4ff", "#56ae6c")

palette2<-c("NoFilter_p2"="#ECA7AE","NoFilter_p5"="#E8436A","NoFilter_p7"="#B13254","NoFilter_p10"="#741434","NoFilter_p50"="#510016", "NoFilter_p100"="#3B0010",
            "Filter1_p2"="#ffc672","Filter1_p5"="#ffaa2c","Filter1_p7"="#f09000","Filter1_p10"="#bf7200","Filter1_p50"="#6D4100","Filter1_p100"="#4E2C00",
            "Filter2_p2"="#bae1d4","Filter2_p5"="#99e2ca","Filter2_p7"="#69bfa3","Filter2_p10"="#488571","Filter2_p50"="#184e3c","Filter2_p100"="#003e2a",
            "Filter3_p2"="#CBD8FF","Filter3_p5"="#8BBCF4","Filter3_p7"="#61A9ED","Filter3_p10"="#3C659C","Filter3_p50"="#2F3777","Filter3_p100"="#17224B")
##############   PLOTTING precision recall  ########################
png(file=paste0(dirI, "Joingen_Precision_vs_Recall_filters.png"), res=350, width=1100, height=800);
ggplot(s, aes(x = recall, y = precision, color = factor(set), shape = factor(eval_level))) +
  geom_point(size=2.5, alpha=0.7) +
  labs(x = "Recall", y = "Precision", color = "Ploidy", shape = "Filter") +
  scale_shape_manual( values= c(18,  2, 20, 3))+
  scale_color_manual(values=palette)+
  theme_minimal()+
  guides(color = guide_legend(override.aes = list(size=2))) +   
  guides(shape = guide_legend(override.aes = list(size=1.5)))+
  # geom_point(shape=3,size=2, alpha=0.6) +
  coord_cartesian(xlim =c(0.55, 1), ylim = c(0.55, 1))+
  theme(legend.key.size = unit(0, 'lines'), 
        legend.key.height = unit(3, 'mm'), #change legend key height
        legend.key.width = unit(0, 'cm'),
        legend.spacing.x = unit(1, 'mm'),
        legend.text = element_text(size=6, margin = margin(l = 0.7, unit = "mm")), 
        legend.title = element_text(size=7.5), 
        axis.text.x = element_text(size=7),  
        axis.text.y = element_text(size=7),
        axis.title = element_text(size=8),
        axis.line = element_line(color = "black", size=0.2),
        legend.margin=margin(-15,0,5,-7),
        legend.box.margin=margin(25,0,0,0))
dev.off();

##############   PLOTTING F1 score  ########################

png(file=paste0(dirI, "Joingen_F1_score_callable_filters_label.png"), res=180, width=1100, height=650);
ggplot(s, aes(x=F1_score, y = file , fill = factor(file))) +
  coord_cartesian(xlim=c(0.3, 0.98))+
  geom_bar(stat="identity")+
  geom_text(aes(label=round(F1_score, 3)), hjust=1.4,  size=3, color="white", fontface="bold")+
  scale_fill_manual(values=palette2)+
  guides(fill="none")+
  theme_minimal()+
  theme( axis.text.x = element_text(size=7,vjust=1, hjust=1),  axis.title.x = element_text(size=8),
         axis.title.y = element_blank(), axis.text.y = element_text(size=7))
dev.off();
##############   PLOTTING BEST ########################
png(file=paste0(dirI, "Joingen_F1_score_callable_joingen_label_best.png"), res=180, width=900, height=300);
ggplot(best, aes(x=F1_score, y = factor(set) , fill = factor(eval_level))) +
  coord_cartesian(xlim=c(0.15, 1))+
  geom_bar(stat="identity")+
  geom_text(aes(label=round(F1_score, 3)), hjust=1.4,  size=3.5, color="white", fontface="bold")+
  scale_fill_manual(values=palette1)+
  labs(fill="", y="Ploidy")+
  theme_minimal()+
  theme( axis.text.x = element_text(size=9,vjust=1, hjust=1),  axis.title = element_text(size=10),
         axis.text.y = element_text(size=10), legend.margin=margin(-5,-4,-2,-22),  legend.key.width = unit(2, 'mm'))
dev.off();

png(file=paste0(dirI, "Joingen_Precision_vs_Recall_filters_zoom_best.png"), res=450, width=1200, height=800);
ggplot(best, aes(x = recall, y = precision, color = factor(set), shape = factor(eval_level))) +
  geom_point(size=2.5, alpha=0.7) +
  labs(x = "Recall", y = "Precision", color = "Ploidy", shape = "Filter") +
  scale_shape_manual( values= c( 16))+
  scale_color_manual(values=palette)+
  theme_minimal()+
  guides(color = guide_legend(override.aes = list(size=2))) +   
  guides(shape = guide_legend(override.aes = list(size=1.5)))+
  # geom_point(shape=3,size=2, alpha=0.6) +
  coord_cartesian(xlim =c(0.55, 1), ylim = c(0.55, 1))+
  theme(legend.key.size = unit(0, 'lines'), 
        legend.key.height = unit(3, 'mm'), #change legend key height
        legend.key.width = unit(0, 'cm'),
        legend.spacing.x = unit(1, 'mm'),
        legend.text = element_text(size=5, margin = margin(l = 0.7, unit = "mm")), 
        legend.title = element_text(size=6), 
        axis.text.x = element_text(size=7),  
        axis.text.y = element_text(size=7),
        axis.title = element_text(size=8),
        axis.line = element_line(color = "black", size=0.2),
        legend.margin=margin(-5,-4,-2,-2),
        legend.box.margin=margin(25,0,0,0))
dev.off();
