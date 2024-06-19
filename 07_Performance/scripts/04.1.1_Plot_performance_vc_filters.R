library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)
####### GET DATA  ######
dir = "~/Desktop/ribosomal_RNAs/Alba/07_Performance/"
dirI<- paste0(dir, "IMAGES/")
# Get all tables
s1<-  read.table(paste0(dir, "Performance_sim_errors_vc_all.tbl"), sep=",", header = T)
s2<-  read.table(paste0(dir, "Performance_sim_errors_vc_filtered_v1.tbl"),sep=",", header = T)
s3<-  read.table(paste0(dir, "Performance_sim_errors_vc_filtered_v2.tbl"), sep=",",header = T)
s4<-  read.table(paste0(dir, "Performance_sim_errors_vc_filtered_v3.tbl"), sep=",",header = T)
s5<-  read.table(paste0(dir, "Performance_sim_errors_vc_filtered_v4.tbl"), sep=",",header = T)


s<-  rbind(s1,s2,s3, s4, s5)

s$eval_level <- c( rep("NoFilter", 9), rep("QUAL1", 9),  rep("QUAL30", 9), rep("QUAL400", 9), rep("QUAL1000", 9))
s$eval_level <-factor(s$eval_level, levels= c("NoFilter","QUAL1", "QUAL30", "QUAL400" , "QUAL1000" ))

##############    Get best ########################
s$file<- factor(paste0(s$eval_level, "_ploidy", s$set), levels = paste0(s$eval_level, "_ploidy", s$set))


best<- data.frame()

for (ploidy in unique(s$set)){
  # keep filter which best performance
  best<- rbind(best, s[s$F1_score == max(s$F1_score[s$set == ploidy]),])
}

write.table(best, paste0(dir,  "Performance_sim_errors_vc_filtered_best.tbl"), sep=",", quote=F, row.names=F, col.names=T)

#####################   Palettes    #########################
# palette <- c("NoFilter" ="#9e2142","QUAL1" ="#f8b95d", "QUAL30" ="#50917c", "QUAL400"  ="#61a9ed" )
palette<- c("100"="#084e3c", "50"="#56ae6c", "40"="royalblue", "30"="darkblue", "20"="#8a76bc", "10"="#e767c5", "7"="#9a1142", "5"="#FF7641", "2"="#ffaa2c")
palette <- c("100" = "#ffaa2c", "50" = "#FF7641", "40" = "#9a1142", "30" = "#e767c5", "20" = "#8a76bc", "10" = "darkblue", "7" = "royalblue", "5" = "#56ae6c", "2" = "#084e3c")

palette1<- c( "QUAL1"="#9a1142", "QUAL30"="#ffaa2c","QUAL400"="royalblue", "QUAL1000"="#488571")


palette2<-c("NoFilter_ploidy2"="#ECA7AE","NoFilter_ploidy5"="#E8436A","NoFilter_ploidy7"="#B13254","NoFilter_ploidy10"="#741434","NoFilter_ploidy50"="#510016", "NoFilter_ploidy100"="#3B0010",
            "QUAL1_ploidy2"="#ffc672","QUAL1_ploidy5"="#ffaa2c","QUAL1_ploidy7"="#f09000","QUAL1_ploidy10"="#bf7200","QUAL1_ploidy50"="#6D4100","QUAL1_ploidy100"="#4E2C00",
            "QUAL30_ploidy2"="#bae1d9","QUAL30_ploidy5"="#99e2ca","QUAL30_ploidy7"="#69bfa3","QUAL30_ploidy10"="#488571","QUAL30_ploidy50"="#184e3c","QUAL30_ploidy100"="#003e2a",
            "QUAL400_ploidy2"="#CBD8FF","QUAL400_ploidy5"="#8BBCF4","QUAL400_ploidy7"="#61A9ED","QUAL400_ploidy10"="#3C659C","QUAL400_ploidy50"="#2F3777","QUAL400_ploidy100"="#17224B")

palette2 <- c(  "NoFilter_ploidy2"="#ECA7AE", "NoFilter_ploidy5"="#E8436A", "NoFilter_ploidy7"="#B13254",  "NoFilter_ploidy10"="#741434", "NoFilter_ploidy20"="#5A101F", 
                "NoFilter_ploidy30"="#440C19","NoFilter_ploidy40"="#3D0916", "NoFilter_ploidy50"="#510016", "NoFilter_ploidy100"="#3B0010",
                
                "QUAL1_ploidy2"="#FFDDA2", "QUAL1_ploidy5"="#FFC9A1", "QUAL1_ploidy7"="#FFA577",   "QUAL1_ploidy10"="#FF8E5C", "QUAL1_ploidy20"="#FF7641",
                "QUAL1_ploidy30"="#FF5F29","QUAL1_ploidy40"="#FF4714", "QUAL1_ploidy50"="#FF3000", "QUAL1_ploidy100"="#CC2500",
                
                "QUAL30_ploidy2"="#ffc672", "QUAL30_ploidy5"="#ffaa2c", "QUAL30_ploidy7"="#f09000",   "QUAL30_ploidy10"="#bf7200", "QUAL30_ploidy20"="#a16000", 
                "QUAL30_ploidy30"="#864f00", "QUAL30_ploidy40"="#704100", "QUAL30_ploidy50"="#6D4100", "QUAL30_ploidy100"="#4E2C00",
              
                "QUAL400_ploidy2"="#CBD8FF", "QUAL400_ploidy5"="#8BBCF4", "QUAL400_ploidy7"="#61A9ED",  "QUAL400_ploidy10"="#3C659C", "QUAL400_ploidy20"="#30528C", 
                "QUAL400_ploidy30"="#28477B",  "QUAL400_ploidy40"="#1F3C6A", "QUAL400_ploidy50"="#2F3777", "QUAL400_ploidy100"="#17224B",
                
                "QUAL1000_ploidy2"="#bae1d9", "QUAL1000_ploidy5"="#99e2ca", "QUAL1000_ploidy7"="#69bfa3",  "QUAL1000_ploidy10"="#488571", "QUAL1000_ploidy20"="#377459", 
                "QUAL1000_ploidy30"="#2a6348",  "QUAL1000_ploidy40"="#20513b", "QUAL1000_ploidy50"="#184e3c", "QUAL1000_ploidy100"="#003e2a")
##############   PLOTTING precision recall  ########################
png(file=paste0(dirI, "VC_Precision_vs_Recall_filters.png"), res=350, width=1000, height=800);
ggplot(s, aes(x = recall, y = precision, color = factor(set), shape = factor(eval_level))) +
  geom_point(size=3, alpha=0.7) +
  labs(x = "Recall", y = "Precision", color = "Ploidy", shape = "Filter") +
  scale_shape_manual( values= c(18, 3, 2, 20, 0))+
  scale_color_manual(values=palette)+
  theme_minimal()+
  guides(color = guide_legend(override.aes = list(size=2))) +   
  guides(shape = guide_legend(override.aes = list(size=1.5)))+
  # geom_point(shape=3,size=2, alpha=0.6) +
  coord_cartesian(xlim =c(0.7, 1), ylim = c(0.7, 1))+
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
        legend.margin=margin(-5,0,0,0),
        legend.box.margin=margin(25,0,0,0))
dev.off();


##############   PLOTTING F1 score  ########################
png(file=paste0(dirI, "VC_F1_score_callable_filters.png"), res=180, width=1000, height=710);
ggplot(s, aes(x=F1_score, y = file , fill = factor(file))) +
  coord_cartesian(xlim=c(0.05, 1))+
  geom_bar(stat="identity")+
  # geom_text(aes(label=type), hjust=2,  color="white", fontface="bold")+
  scale_fill_manual(values=palette2)+
  guides(fill="none")+
  theme_minimal()+
  theme( axis.text.x = element_text(size=7,vjust=1, hjust=1),  axis.title.x = element_text(size=8),
         axis.title.y = element_blank(), axis.text.y = element_text(size=6))
dev.off();


png(file=paste0(dirI, "VC_F1_score_callable_filters_label.png"), res=180, width=1100, height=850);
ggplot(s, aes(x=F1_score, y = file , fill = factor(file))) +
  coord_cartesian(xlim=c(0.05, 1))+
  geom_bar(stat="identity")+
  geom_text(aes(label=round(F1_score, 3)), hjust=1.4,  size=2.5, color="white", fontface="bold")+
  scale_fill_manual(values=palette2)+
  guides(fill="none")+
  theme_minimal()+
  theme( axis.text.x = element_text(size=7,vjust=1, hjust=1),  axis.title.x = element_text(size=8),
         axis.title.y = element_blank(), axis.text.y = element_text(size=6))
dev.off();
##############   PLOTTING FBEST  ########################
png(file=paste0(dirI, "VC_Precision_vs_Recall_filters_zoom_best.png"), res=400, width=1000, height=800);
ggplot(best, aes(x = recall, y = precision, color = factor(set), shape = factor(eval_level))) +
  geom_point(size=2.5, alpha=0.7) +
  labs(x = "Recall", y = "Precision", color = "Ploidy", shape = "Filter") +
  scale_shape_manual( values= c(  3, 2, 20, 0))+
  scale_color_manual(values=palette)+
  theme_minimal()+
  guides(color = guide_legend(override.aes = list(size=2))) +   
  guides(shape = guide_legend(override.aes = list(size=1.5)))+
  # geom_point(shape=3,size=2, alpha=0.6) +
  coord_cartesian(xlim =c(0.7, 0.95), ylim = c(0.91, 0.94))+
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

png(file=paste0(dirI, "VC_F1_score_callable_filters_label_best.png"), res=180, width=1300, height=450);
ggplot(best, aes(x=F1_score, y = factor(set) , fill = factor(eval_level))) +
  coord_cartesian(xlim=c(0.6, 1))+
  geom_bar(stat="identity", alpha=0.9)+
  geom_text(aes(label=round(F1_score, 3)), hjust=1.4,  size=3, color="white", fontface="bold")+
  scale_fill_manual(values=palette1)+
  labs( y="Ploidy", fill="")+
  theme_minimal()+
  theme( axis.text.x = element_text(size=7,vjust=1, hjust=1),  axis.title = element_text(size=8),
         legend.key.width = unit(2, 'mm'),legend.margin=margin(-15,-4,-2,-22),    axis.text.y = element_text(size=7))
dev.off();
