# ## --------------------------------------  CHOOSE VC SET ---------------------------------   ############
# # GATK version either "4.5.0.0" or "4.3.0.0"
set="4.3.0.0"
# #set=  "4.5.0.0"
# ## STEP of the pipeline this is for gVCF files so either "VC" or"filtering"
step="VC"
step= "filtering"
# 
version="2sample_v"
# ## --------------------------------------  LOAD DATA ---------------------------------   ############
print("#############################################################################")

#####  Library #####  
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)
#####  Palette  #####  
palette <-  c("TP" ="#75bb6f","TP_FN_allele"="#aeb852","TP_pos_other"="#706f4d",  "FN_uncalled"="#ffc60e", "FN_ref_homozygous"= "#4c89ee","FN_alt_homozygous"="#71a5be",
              "FP"= "#bd0752", "TP_pos_FP"= "#ff4d7a", "TP_pos_FP_allele"="#ed8c96","TP_allele_FP_allele"="#ed8c96",  "err"= "#a6a6a6") 
## ------------------------------------ SET PARAMETERS --------------------------------- ############
dir <- "~/Desktop/ribosomal_RNAs/Alba/08_validation/"
dirf<- paste0(dir, "VC_out/")
dirI<- paste0(dir, "IMAGES/joingen_analysis/")
# dirI<- paste0(dir, "IMAGES/VC_analysis/")

## --------------------------------------  LOAD DATA ---------------------------------   ############
# d1<- read.table(paste0(dirf, "VCeval_prejoin_p50_VC.tbl"), sep="\t")
# d1$ploidy<- paste(d1$ploidy, "ploidy")
# d2<- read.table(paste0(dirf, "VCeval_prejoin_p50_filtering_2sample_v4.tbl"), sep="\t")
# d2$ploidy<- paste(d2$ploidy, "ploidy QUAL1000 2 samples")


d1<- read.table(paste0(dirf, "VCeval_postjoing_p7_joingen.tbl"), sep="\t")
d1$ploidy<- paste(d1$ploidy, "ploidy ")
d2<- read.table(paste0(dirf, "VCeval_postjoing_p7_filtering_v2.tbl"), sep="\t")
d2$ploidy<- paste(d2$ploidy, "ploidy Filter2")


d <- rbind(d1, d2)

d$eval_sample<-  factor(d$eval_sample, levels=c("FN_alt_homozygous", "FN_ref_homozygous",  "FN_uncalled","TP_FN_allele" ,"TP_pos_FN_allele", "FP","TP_pos_FP", "TP_allele_FP_allele", "TP"))

d_plot<- d[d$eval_sample!="FN_uncalled",]

########------------------------------- Plotting variables per position -----------------------------------########

png(file=paste0(dirI,"QUAL", "_distribution.png"), res=500, width=3000, height=1300);
ggplot(d_plot, aes(x=eval_sample, y=as.numeric(QUAL), color=eval_sample) )+
  # coord_cartesian(ylim = c(0,10000))+
  coord_cartesian(ylim = c(0,100000))+
  geom_hline(yintercept=1000, color="blue", linetype="dashed", linewidth=0.6)+
  geom_text(aes(label = "1000", x="FN_ref_homozygous", y =3500), color = "blue", size=3)+
  geom_jitter(alpha=0.3, size=1)+
  geom_boxplot(alpha=0.3, size=0.5, color="black", width = 0.4, outlier.shape = NA)+
  scale_color_manual(values=palette)+
  labs(y="quality", x="")+
  facet_grid(.~ploidy)+
  guides(color="none")+
  theme_classic() +
  theme(axis.text.x = element_text(size=8, angle=17, vjust=0.75, hjust=0.6),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"),
        panel.grid.minor = element_line(colour = "#e6e6e6"))
dev.off();


png(file=paste0(dirI,"QUAL", "30_distribution.png"), res=500, width=3000, height=1300);
ggplot(d_plot, aes(x=eval_sample, y=as.numeric(QUAL), color=eval_sample) )+
  # coord_cartesian(ylim = c(0,10000))+
  coord_cartesian(ylim = c(0,3000))+
  geom_hline(yintercept=30, color="blue", linetype="dashed", linewidth=0.6)+
  geom_text(aes(label = "30", x="FN_ref_homozygous", y =500), color = "blue", size=3)+
  geom_jitter(alpha=0.3, size=1)+
  geom_boxplot(alpha=0.3, size=0.5, color="black", width = 0.4, outlier.shape = NA)+
  scale_color_manual(values=palette)+
  labs(y="quality", x="")+
  facet_grid(.~ploidy)+
  guides(color="none")+
  theme_classic() +
  theme(axis.text.x = element_text(size=8, angle=17, vjust=0.75, hjust=0.6),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"),
        panel.grid.minor = element_line(colour = "#e6e6e6"))
dev.off();

# find possible filter know more about the test: https://gatk.broadinstitute.org/hc/en-us/articles/360042477772-BaseQualityRankSumTest
png(file=paste0(dirI,"BaseQRankSum", "_distribution.png"),  res=500, width=2500, height=2000);
ggplot(d_plot, aes(x=eval_sample, y=as.numeric(BaseQRankSum), color=eval_sample) )+
  geom_jitter(alpha=0.6, size=0.4, width = 0.4)+
  geom_boxplot(alpha=0.3, size=0.5, color="black", width = 0.4, outlier.shape = NA)+
  scale_color_manual(values=palette)+
  geom_hline(yintercept=-6, color="blue", linetype="dashed", linewidth=0.6)+
  geom_text(aes(label = "-6", x="FN_alt_homozygous", y = -5), color = "blue", size=3)+
  facet_grid(.~ploidy)+
  labs(y="BaseQRankSum", x="", color="")+
  guides(color = guide_legend(nrow = 1, override.aes = list(size=3)))+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.margin = margin(-15, -10, 0, 0),
        strip.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 6),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "#7d7d7d", size = 0.5))
dev.off();
# png(file=paste0(dirI,"ExcessHet", "_distribution.png"), res=500, width=7500, height=2000);
# ggplot(d_plot, aes(x=eval_sample, y=as.numeric(ExcessHet), color=eval_sample) )+
#   geom_jitter(alpha=0.3, size=1)+
#   scale_color_manual(values=palette)+
#   facet_grid(.~ploidy)+
#   guides(color="none")+
#   labs(y="ExcessHet", x="")+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle=45, vjust=0.5))
# dev.off();
png(file=paste0(dirI,"MQRankSum", "_distribution.png"), res=230, width=1500, height=500);
ggplot(d_plot, aes(x=eval_sample, y=as.numeric(MQRankSum), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  geom_boxplot(alpha=0.1, color="black",  outlier.shape = NA)+
  scale_color_manual(values=palette)+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="MQRankSum", x="", color="")+
  guides(color = guide_legend(nrow = 1, override.aes = list(size=3)))+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.margin = margin(-15, -10, 0, 0),
        strip.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 6),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "#7d7d7d", size = 0.5))
dev.off();

png(file=paste0(dirI,"ReadPosRankSum", "_4_distribution.png"), res=230, width=1500, height=500);
ggplot(d_plot, aes(x=eval_sample, y=as.numeric(ReadPosRankSum), color=eval_sample) )+
  geom_jitter(alpha=0.7, size=0.7)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.3, size=0.5, color="black", width = 0.4, outlier.shape = NA)+
  geom_hline(yintercept=4, color="blue", linetype="dashed", linewidth=0.6)+
  geom_text(aes(label = "4", x="FN_alt_homozygous", y=7), color = "blue", size=5)+
  facet_grid(.~ploidy)+
  guides(color = guide_legend(nrow = 1, override.aes = list(size=3)))+
  labs(y="ReadPosRankSum", x="", color="")+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.margin = margin(-15, -10, 0, 0),
        strip.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "#7d7d7d", size = 0.5))
dev.off();

png(file=paste0(dirI,"ReadPosRankSum", "_10_distribution.png"), res=230, width=1500, height=500);
ggplot(d_plot, aes(x=eval_sample, y=as.numeric(ReadPosRankSum), color=eval_sample) )+
  geom_jitter(alpha=0.7, size=0.7)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.3, size=0.5, color="black", width = 0.4, outlier.shape = NA)+
  geom_hline(yintercept=10, color="blue", linetype="dashed", linewidth=0.6)+
  geom_text(aes(label = "10", x="FN_alt_homozygous", y = 12), color = "blue", size=5)+
  facet_grid(.~ploidy)+
  guides(color = guide_legend(nrow = 1, override.aes = list(size=3)))+
  labs(y="ReadPosRankSum", x="", color="")+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.margin = margin(-15, -10, 0, 0),
        strip.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "#7d7d7d", size = 0.5))
dev.off();

s<- d[d$eval_sample != "FN_uncalled" & !is.na(d$SOR), ]

png(file=paste0(dirI,"SOR", "_distribution.png"), res=230, width=1500, height=500);
ggplot(s, aes(x=eval_sample, y=as.numeric(SOR), color=eval_sample) )+
  geom_jitter(alpha=0.6, size=0.7)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.3, size=0.5, color="black", width = 0.4, outlier.shape = NA)+
  facet_grid(.~ploidy)+
  guides(color = guide_legend(nrow = 1, override.aes = list(size=3)))+
  labs(y="SOR", x="", color="")+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.margin = margin(-15, -10, 0, 0),
        strip.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "#7d7d7d", size = 0.5))
dev.off();


png(file=paste0(dirI,"DP.1", "_distribution.png"), res=230, width=1500, height=500);
ggplot(d_plot, aes(x=eval_sample, y=as.numeric(DP.1), color=eval_sample) )+
  geom_jitter(alpha=0.6, size=0.7)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.3, size=0.5, color="black", width = 0.4, outlier.shape = NA)+
  facet_grid(.~ploidy)+
  guides(color = guide_legend(nrow = 1, override.aes = list(size=3)))+
  labs(y="DP.1", x="", color="")+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.margin = margin(-15, -10, 0, 0),
        strip.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "#7d7d7d", size = 0.5))
dev.off();


png(file=paste0(dirI,"GQ", "_distribution.png"), res=230, width=1500, height=500);
ggplot(d_plot, aes(x=eval_sample, y=GQ, color=eval_sample) )+
  geom_jitter(alpha=0.6, size=0.5)+
  geom_boxplot(alpha=0.1, color="black",  width=0.3, outlier.shape = NA)+
  scale_color_manual(values=palette)+
  facet_grid(.~ploidy)+
  guides(color = guide_legend(nrow = 1, override.aes = list(size=3)))+
  labs(y="GQ", x="", color="")+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.margin = margin(-15, -10, 0, 0),
        strip.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "#7d7d7d", size = 0.5))
dev.off();

png(file=paste0(dirI,"GQ", "_zoom.png"), res=130, width=700, height=300);
ggplot(d_plot, aes(x=eval_sample, y=GQ, color=eval_sample) )+
  geom_jitter(alpha=0.6, size=0.5)+
  geom_boxplot(aes(x=1.5), alpha=0.1, color="black",  width=1, outlier.shape = NA)+
  scale_color_manual(values=palette)+
  facet_grid(.~ploidy)+
  guides(color = "none")+
  labs(y="GQ", x="", color="")+
  theme_bw()+
  theme(
    # axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.margin = margin(-15, -10, 0, 0),
        strip.text = element_text(size = 10, face = "bold"),
        legend.text = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "#7d7d7d", size = 0.5))
dev.off();
print("All done!")
print("#############################################################################")
# #SB
