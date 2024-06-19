library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)
####### GET DATA  ######
dir = "~/Desktop/ribosomal_RNAs/Alba/07_Performance/"
dirI<- paste0(dir, "IMAGES/")

ploidies<- c(2, 5, 7, 10, 50, 100)
# Get all tables
s1<-  read.table(paste0(dir, "Performance_perfect_reads.tbl"), sep=",", header = T)[,1:6]
s2<-  read.table(paste0(dir, "Performance_sim_errors_vc_all.tbl"),sep=",", header = T)[,1:6]
s3<-  read.table(paste0(dir, "Performance_sim_errors_vc_filtered_2sample.tbl"), sep=",",header = T)[,1:6]
s4<-  read.table(paste0(dir, "Performance_sim_errors_postjoin.tbl"), sep=",",header = T)[,1:6]
s5<-  read.table(paste0(dir, "Performance_sim_errors_joingen_filtered_best.tbl"), sep=",",header = T)[,1:6]

s1$eval_level<- "NoFilter"
s1$type<- "vc kmers"
s3$type<- "vc"
s4$type<- "genotype"
s5$type<- "genotype"

s<-  rbind(s1,s2,s3, s4, s5)

s$type <-factor( c("vc kmers 30x", "vc kmers 100x", "vc kmers 1000x", "vc ploidy 2",
                   "vc ploidy 5", "vc ploidy 7", "vc ploidy 10", "vc ploidy 20",
                   "vc ploidy 30", "vc ploidy 40", "vc ploidy 50", "vc ploidy 100", "vc ploidy 2",
                   "vc ploidy 5", "vc ploidy 7", "vc ploidy 10", "vc ploidy 20",
                   "vc ploidy 30", "vc ploidy 40", "vc ploidy 50", "vc ploidy 100",
                   "genotype ploidy 2", "genotype ploidy 5", "genotype ploidy 7", "genotype ploidy 10",
                   "genotype ploidy 50", "genotype ploidy 100", "genotype ploidy 2", "genotype ploidy 5",
                   "genotype ploidy 7", "genotype ploidy 10", "genotype ploidy 50", "genotype ploidy 100")      , 
                levels = c("vc kmers 30x", "vc kmers 100x", "vc kmers 1000x", "vc ploidy 2",
                           "vc ploidy 5", "vc ploidy 7", "vc ploidy 10", "vc ploidy 20",
                           "vc ploidy 30", "vc ploidy 40", "vc ploidy 50", "vc ploidy 100",
                           "genotype ploidy 2", "genotype ploidy 5", "genotype ploidy 7", "genotype ploidy 10",
                           "genotype ploidy 50", "genotype ploidy 100"))

s$Filter<- c( rep("Not filtered", 12), rep("Filtered", 9), rep("Not filtered", 6), rep("Filtered", 6))

s$Reads<- c(rep("Perfect reads", 3), rep("Sequencing errors", 30))

##############   PLOTTING precision recall  ########################
# palette <-  c( "vc kmers 30x" = "#9bb5ac", "vc kmers 100x" = "#50917c", "vc kmers 1000x"="#001e2d",
#                "vc 2"= "#ff9bc1", "vc 5"= "#ff2b83", "vc 7"="#9e2142","vc 10"="#590010",
#                "vc 2 filtered"= "#f8b95d", "vc 5 filtered"= "#ff850a", "vc 7 filtered"="#ad4700","vc 10 filtered"="#583540",
#                "join-genotype 2"="#7ac8e6", "join-genotype 5"="#61a9ed", "join-genotype 7"="#155bb0", "join-genotype 10"="#001f6a")
# 
# "NoFilter_p2"="#ECA7AE","NoFilter_p5"="#E8436A","NoFilter_p7"="#B13254","NoFilter_p10"="#741434","NoFilter_p50"="#510016", "NoFilter_p100"="#3B0010",
# "Filter1_p2"="#ffc672","Filter1_p5"="#ffaa2c","Filter1_p7"="#f09000","Filter1_p10"="#bf7200","Filter1_p50"="#6D4100","Filter1_p100"="#4E2C00",
# "Filter2_p2"="#bae1d4","Filter2_p5"="#99e2ca","Filter2_p7"="#69bfa3","Filter2_p10"="#488571","Filter2_p50"="#184e3c","Filter2_p100"="#003e2a",
# "Filter3_p2"="#CBD8FF","Filter3_p5"="#8BBCF4","Filter3_p7"="#61A9ED","Filter3_p10"="#3C659C","Filter3_p50"="#2F3777","Filter3_p100"="#17224B")

palette<- c("vc kmers 30x" = "#9bb5ac", "vc kmers 100x" = "#50917c", "vc kmers 1000x"="#001e2d",
            "vc ploidy 2"="#EaA7AE", "vc ploidy 5"="#FF8C8C", "vc ploidy 7"="#F05F5F", "vc ploidy 10"="#E03F3F","vc ploidy 20"="#D00000", 
            "vc ploidy 30"="#A20000", "vc ploidy 40"="#6A0000", "vc ploidy 50"="#3e0711", "vc ploidy 100"="#1B0010",
           "genotype ploidy 2"="#CBD8FF","genotype ploidy 5"="#8BBCF4","genotype ploidy 7"="#61A9ED","genotype ploidy 10"="#3C659C",  
           "genotype ploidy 50"="#2F3777",  "genotype ploidy 100"="#17224B")

png(file=paste0(dirI, "ALL_Precision_vs_Recall_callable.png"), res=350, width=1300, height=1000);
ggplot(s, aes(x = recall, y = precision, color = type, shape = Filter)) +
  geom_point(size=5, alpha=0.6) +
  labs(x = "Recall", y = "Precision", color = "Set of reads", shape = "Evaluation level") +
  scale_shape_manual( values= c(18, 20))+
  scale_color_manual(values=palette)+
  theme_minimal()+
  guides(color = guide_legend(override.aes = list(size=2.4))) +   
  guides(shape = guide_legend(override.aes = list(size=1.5)))+
  geom_point(shape=3,size=3, alpha=0.8) +
  coord_cartesian(xlim =c(0.55, 1), ylim = c(0.55, 1))+
  theme(legend.key.size = unit(0, 'lines'), 
        legend.key.height = unit(3, 'mm'), #change legend key height
        legend.key.width = unit(0, 'cm'),
        legend.spacing.x = unit(1, 'mm'),
        legend.spacing.y = unit(0, 'mm'),
        legend.text = element_text(size=6, margin = margin(l = 0.5, unit = "mm")), 
        legend.title = element_text(size=6), 
        axis.text.x = element_text(size=6),  
        axis.text.y = element_text(size=6),
        axis.title = element_text(size=9),
        axis.line = element_line(color = "black", size=0.2),
        legend.margin=margin(2,0,0,-5),
        legend.box.margin=margin(15,0,0,0))
dev.off();
##############   PLOTTING  F1 score ########################

s$typef1 <- factor(c("vc ploidy 5 kmers 30x", "vc ploidy 5 kmers 100x", "vc ploidy 5 kmers 1000x", 
                     "vc ploidy 2", "vc ploidy 5", "vc ploidy 7", "vc ploidy 10",   "vc ploidy 20", "vc ploidy 30", "vc ploidy 40", "vc ploidy 50",  "vc ploidy 100", 
                     "vc ploidy 2 2samples QUAL30", "vc ploidy 5 2samples QUAL30", "vc ploidy 7 2samples QUAL30", "vc ploidy 10 2samples QUAL1000", "vc ploidy 20 2samples QUAL400", "vc ploidy 30 2samples QUAL1000", 
                     "vc ploidy 40 2samples QUAL1000", "vc ploidy 50 2samples QUAL1000", "vc ploidy 100 2samples QUAL1000",  "genotype ploidy 2", "genotype ploidy 5", "genotype ploidy 7", 
                     "genotype ploidy 10", "genotype ploidy 50", "genotype ploidy 100",  "genotype ploidy 2 Filter2", "genotype ploidy 5 Filter2", 
                     "genotype ploidy 7 Filter2", "genotype ploidy 10 Filter2",   "genotype ploidy 50 Filter2", "genotype ploidy 100 Filter2"), 
                   levels= c("genotype ploidy 100 Filter2", "genotype ploidy 50 Filter2", "genotype ploidy 10 Filter2",  
                             "genotype ploidy 7 Filter2", "genotype ploidy 5 Filter2", "genotype ploidy 2 Filter2", 
                             "genotype ploidy 100", "genotype ploidy 50", "genotype ploidy 10", "genotype ploidy 7", "genotype ploidy 5", "genotype ploidy 2", 
                             "vc ploidy 100 2samples QUAL1000", "vc ploidy 50 2samples QUAL1000", "vc ploidy 40 2samples QUAL1000",  "vc ploidy 30 2samples QUAL1000", "vc ploidy 20 2samples QUAL400", "vc ploidy 10 2samples QUAL1000", 
                             "vc ploidy 7 2samples QUAL30", "vc ploidy 5 2samples QUAL30", "vc ploidy 2 2samples QUAL30", 
                             "vc ploidy 100", "vc ploidy 50", "vc ploidy 40", "vc ploidy 30", "vc ploidy 20", "vc ploidy 10", "vc ploidy 7", "vc ploidy 5", 
                             "vc ploidy 2", "vc ploidy 5 kmers 30x", "vc ploidy 5 kmers 100x", "vc ploidy 5 kmers 1000x"))

palette2 <- c("vc ploidy 5 kmers 30x"="#99d1bf", "vc ploidy 5 kmers 100x"= "#50917c", "vc ploidy 5 kmers 1000x"="#184e3c",
              "vc ploidy 2" = "#FFC672", "vc ploidy 5" = "#FFAA2C", "vc ploidy 7" = "#F09000", "vc ploidy 10" = "#D37F00", "vc ploidy 20" = "#BF7200", 
              "vc ploidy 30" = "#A15F00", "vc ploidy 40" = "#774f00", "vc ploidy 50" = "#674300", "vc ploidy 75" = "#563700", "vc ploidy 100" = "#463000",
              
              "vc ploidy 2 2samples QUAL30" = "#EaA7AE", "vc ploidy 5 2samples QUAL30" = "#E8436A", "vc ploidy 7 2samples QUAL30" = "#B13254", "vc ploidy 10 2samples QUAL1000" = "#941434", "vc ploidy 20 2samples QUAL400" = "#841434",
              "vc ploidy 30 2samples QUAL1000" = "#71101F", "vc ploidy 40 2samples QUAL1000" = "#540C19", "vc ploidy 50 2samples QUAL1000" = "#3B0010", "vc ploidy 100 2samples QUAL1000" = "#180006",
              
              "genotype ploidy 2" = "#9BBCF4", "genotype ploidy 5" = "#7ab7d9", "genotype ploidy 7" = "#61A9ED", "genotype ploidy 10" = "#3C659C", "genotype ploidy 50" = "#2F3777", "genotype ploidy 100" = "#17224B",
              "genotype ploidy 2 Filter2" = "#e3cfed", "genotype ploidy 5 Filter2" = "#d7aaff", "genotype ploidy 7 Filter2" = "#c187f4", "genotype ploidy 10 Filter2" = "#9b71c0", "genotype ploidy 50 Filter2" = "#5c1995",
              "genotype ploidy 100 Filter2" = "#3A1459")


png(file=paste0(dirI, "ALL_F1_score_callable_labeled.png"), res=200, width=1500, height=900);
ggplot(s, aes(x=F1_score, y = typef1 , fill = typef1)) +
  coord_cartesian(xlim=c(0.05, 0.96))+
  geom_bar(stat="identity")+
  geom_text(aes(label=round(F1_score, 4)), hjust=1.5, size=3,  color="white", fontface="bold")+
  scale_fill_manual(values=palette2)+
  guides(fill="none")+
  theme_minimal()+
  theme( axis.text.x = element_text(size=7,vjust=1, hjust=1), 
         axis.title.y = element_blank(), axis.title.x = element_text(size=8))
dev.off();

# png(file=paste0(dirI, "ALL_F1_score_callable.png"), res=180, width=1100, height=700);
# ggplot(s, aes(x=F1_score, y = typef1 , fill = typef1)) +
#   coord_cartesian(xlim=c(0.05, 1))+
#   geom_bar(stat="identity")+
#   # geom_text(aes(label=type), hjust=2,  color="white", fontface="bold")+
#   scale_fill_manual(values=palette2)+
#   guides(fill="none")+
#   theme_minimal()+
#   theme( axis.text.x = element_text(size=7,vjust=1, hjust=1), 
#          axis.title.y = element_blank(), axis.title.x = element_text(size=8))
# dev.off();

