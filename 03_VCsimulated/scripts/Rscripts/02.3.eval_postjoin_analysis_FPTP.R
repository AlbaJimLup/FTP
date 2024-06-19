###### Libraries ######
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)
####  Palettes  #####
palette <-  c("TP" ="#75bb6f","TP_FN_allele"="#aeb852",    "TP_pos_other"="#706f4d", 
              "FN_uncalled"="#ffc60e", "FN_ref_homozygous"= "#4c89ee","FN_alt_homozygous"="#71a5be",
              "FP"= "#bd0752", "TP_pos_FP"= "#ff4d7a","TP_pos_FP_sample"= "#ad4752",  "TP_pos_FP_allele"="#ed8c96","TP_allele_FP_allele"="#ed8c96",  "err"= "#a6a6a6") 

palette2 <- c( "5_ETS" = "#2f4b5b", "18S" = "#6ad68e", "ITS1" = "#9dbccf", "5.8S" = "#90c230",  "ITS2" = "#4d8195",  "28S" = "#71d064",  "3_ETS" = "#375794")

palette3 <- c( "5_ETS" = "#2f4b5b", "18S" = "#ffe290",  "ITS1" = "#9dbccf","5.8S" = "#fdae61",  "ITS2" = "#4d8195",  "28S" = "#ff6e7d",  "3_ETS" = "#375794")

## --------------------------------------  LOAD DATA ---------------------------------   ############

## GATK version either "4.5.0.0" or "4.3.0.0"
# set="4.3.0.0" 
## STEP of the pipeline this is for gVCF files so either "VC" or"filtering"
# step="joingen"
# step= "filtering"
# version="v1"

##############------------------------ Data prep -------------------------################
#dir <-  "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/pre_join_gen/"
# dir <- "~/Desktop/ribosomal_RNAs/Alba/03_VCsimulated/v2/"
# dirf<- paste0(dir, "VC_out/Sample_level_eval/4.3.0.0/joingen/")
# dirt <- paste0(dir, "scripts/Rscripts/")
# # dirI<- paste0(dir, "IMAGES/02_POSTJOIN/4.3.0.0/", step, "/Analysis/")
# # # true_variants_sample <- read.table(paste0(dirt, "true_variants.tbl"), header=T)
# # # catalog <-  read.table(paste0(dirt, "variants_AF_Counts.tbl"))
# if(step=="joingen"){
#   dirI<- paste0(dir, "IMAGES/02_POSTJOIN/", set, "/joingen/Analysis/")
#   dfs <- read.table( paste0(dirf,"VCeval_postjoin_joingen.tbl"), sep="\t")
#   
# }else{
#   dirI<- paste0(dir, "IMAGES/02_POSTJOIN/", set, "/filtering_", version, "/Analysis/")
#   dfs <- read.table( paste0(dirf,"VCeval_postjoin_filtering_", version, ".tbl"), sep="\t")
# }

dirI<- paste0(dirI, "/Analysis/")
d<-dfs
####### Data frame to compare FP and TP called variants ##########
# tpfp<- d[d$eval_sample %in% c("TP", "FP", "TP_allele_FP_allele"), ]
tpfp<- d[d$eval_sample !="FN_uncalled", ]
tpfp$region<- factor(tpfp$region, levels =c("5_ETS", "18S", "ITS1", "5.8S", "ITS2",  "28S", "3_ETS") )
# 
###### Data frame to compare FN regions  ########
fn <-  d[d$eval_sample == "FN_uncalled", ]

# Get table
fnt<- as.data.frame(table(fn$region, fn$ploidy)) # get counts per regions
colnames(fnt)<- c("Region", "Ploidy", "Counts")
fnt$Region<-  factor(fnt$Region, levels =c("5_ETS", "18S", "ITS1", "5.8S", "ITS2",  "28S", "3_ETS"))
#lengths per region
fnt$length<- 0
fnt$length[fnt$Region =="5_ETS"]<- 3652
fnt$length[fnt$Region =="18S"]<- 1869
fnt$length[fnt$Region =="ITS1"]<- 1077
fnt$length[fnt$Region =="5.8S"]<- 157 
fnt$length[fnt$Region =="ITS2"]<- 1168
fnt$length[fnt$Region =="28S"]<- 5080
fnt$length[fnt$Region =="3_ETS"]<- 361 
# Normalize by length
fnt$norm <- fnt$Counts / fnt$length
##### Do it per sample #####
fn_sample<- as.data.frame(table(fn$region, fn$ploidy, fn$sample)) 
colnames(fn_sample)<- c("Region", "Ploidy", "Sample", "Counts")
fn_sample<-  fn_sample[fn_sample$Counts !=0,]
fn_sample$Region<-  factor(fn_sample$Region, levels =c("5_ETS", "18S", "ITS1", "5.8S", "ITS2",  "28S", "3_ETS") )
fn_sample$Sample<-  factor(fn_sample$Sample, levels =paste0("sample_", 1:100))
fn_sample$length<- 0
fn_sample$length[fn_sample$Region =="5_ETS"]<- 3652
fn_sample$length[fn_sample$Region =="18S"]<- 1869
fn_sample$length[fn_sample$Region =="ITS1"]<- 1077
fn_sample$length[fn_sample$Region =="5.8S"]<- 157 
fn_sample$length[fn_sample$Region =="ITS2"]<- 1168
fn_sample$length[fn_sample$Region =="28S"]<- 5080
fn_sample$length[fn_sample$Region =="3_ETS"]<- 361 
# Normalize by length
fn_sample$norm <- fn_sample$Counts / fn_sample$length

##############------------------------ Plotting FNs-------------------------################
png(file=paste0(dirI,"FN_variants_region.png"), res=500, width=2000, height=3000);
ggplot(fnt, aes(x=Region, y=Counts, fill=Region))+
  geom_bar(stat = "identity")+
  labs(y="Number of variants not called", x="region", fill="Region")+
  geom_text(aes(label=Counts),vjust = -0.2,  size=3 )+
  guides(fill="none")+
  scale_fill_manual(values=palette2)+
  facet_grid(Ploidy~.)+
  theme_classic()
dev.off();

png(file=paste0(dirI,"FN_normalized_variants_region.png"), res=450, width=2000, height=3000);
ggplot(fnt, aes(x=Region, y=norm, fill=Region))+
  geom_bar(stat = "identity")+
  geom_text(aes(label=Counts),vjust = 0.2,  size=4 )+
  labs(y="Normalized number of variants not called", x="Region", fill="Region")+
  guides(fill="none")+
  scale_fill_manual(values=palette2)+
  facet_grid(Ploidy~.)+
  theme_classic()
dev.off();

png(file=paste0(dirI,"FN_samples_region.png"), res=400, width=5000, height=3000);
ggplot(fn_sample, aes(x=Sample, y=Counts, fill=Region))+
  geom_bar(stat = "identity")+
  labs(y="Number of variants not called", x="region", fill="Region")+
  geom_text(aes(label=Counts),position = "stack", vjust=2, size=2, color="white" )+
  scale_fill_manual(values=palette3)+
  facet_grid(Ploidy~.)+
  theme_classic()+
  theme( legend.position = "bottom", axis.text.x = element_text(angle=60, vjust = 0.5, size=8))
dev.off();

png(file=paste0(dirI,"FN_samples_region_norm.png"),  res=400, width=5200, height=3000);
ggplot(fn_sample, aes(x=Sample, y=norm, fill=Region))+
  geom_bar(stat = "identity")+
  labs(y="Number of variants not called", x="region", fill="Region")+
  scale_fill_manual(values=palette3)+
  facet_grid(Ploidy~.)+
  theme_classic()+
  theme( legend.position = "bottom", axis.text.x = element_text(angle=60, vjust = 0.5))
dev.off();
##############------------------------ Plotting GQ -------------------------################
png(file=paste0(dirI,"00_GQ", "_boxplot_ploidy_distribution.png"), res=500, width=2000, height=1800);
ggplot(tpfp, aes(x = as.factor(ploidy), y = as.numeric(GQ))) +
  geom_jitter(aes(color = eval_sample), alpha = 0.3, size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
  scale_color_manual(values = palette) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "black", position = position_dodge(width = 0.75)) +
  stat_summary(  fun = function(x) quantile(x, probs = 0.25, na.rm = TRUE),
                  geom = "text",     aes(label = sprintf("%0.2f", ..y..)),
                  position = position_dodge(width = 0.75), vjust = -0.5, size = 3, color = "black" ) +
  guides(color = "none") +
  labs(y = "GQ", x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, size = 8, hjust = 1),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"),
        panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"00_GQ", "_boxplot_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(GQ), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black" )+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="GQ", x="")+
  theme_classic()+  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) , 
                     axis.title = element_text(size=8),
                     strip.text = element_text(size = 10, face = "bold"),
                     panel.grid.minor = element_line(colour = "#e6e6e6"),  
                     panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                     panel.grid.major = element_line(colour = "#e6e6e6"), 
                     strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"00_GQ", "_violin_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(GQ), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_violin(alpha=0.5,  color="black" )+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="GQ", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) , 
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"00_GQ", "_violin_zoom.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(GQ), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  coord_cartesian(yli=c(19, 50))+
  geom_hline(yintercept=40, color="blue", linetype="dashed", linewidth=0.6)+
  scale_color_manual(values=palette)+
  geom_violin(alpha=0.5,  color="black" )+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="GQ", x="")+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) ,
                     axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
                     strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();
##############------------------------ Plotting other Variants -------------------------################

png(file=paste0(dirI,"QUAL", "_distribution.png"), res=500, width=5500, height=1800);
ggplot(tpfp, aes(x=eval_sample, y=as.numeric(QUAL), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  labs(y="quality", x="", color="Variant")+
  guides(fill = guide_legend(nrow = 1))+
  facet_grid(.~ploidy)+
  guides(color="none")+
  theme_classic()+          
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) , 
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
  # theme(axis.text.x =  element_text(angle=25,margin = margin(t = 20) ))
dev.off();


png(file=paste0(dirI,"QUAL", "_30distribution.png"), res=500, width=5500, height=1800);
ggplot(tpfp, aes(x=eval_sample, y=as.numeric(QUAL), color=eval_sample) )+
  coord_cartesian(ylim = c(0,1000000))+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  geom_hline(yintercept=30, color="blue", linetype="dashed", linewidth=0.6)+
  labs(y="quality", x="")+
  facet_grid(.~ploidy)+
  guides(color="none")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) , 
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"AC", "_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(AC), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  scale_color_manual(values=palette)+  
  geom_text(aes(label = "100", x="FN_alt_homozygous", y = 200), color = "blue", size=3)+ 
  geom_hline(yintercept=100, color="blue", linetype="dashed", linewidth=0.5)+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="AC", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) ,  
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"AF", "_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(AF), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  geom_text(aes(label = "0.02", x="FN_alt_homozygous", y = 0.05), color = "blue", size=3)+ 
  geom_hline(yintercept=0.02, color="blue", linetype="dashed", linewidth=0.5)+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="AF", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) ,   
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"AN", "_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(AN), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="AN", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) , 
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"BaseQRankSum", "_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(BaseQRankSum), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  geom_text(aes(label = "-2", x="TP_pos_TP_allele_FP_allele", y = -3), color = "blue", size=3)+ 
  geom_hline(yintercept=-2, color="blue", linetype="dashed", linewidth=0.5)+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="BaseQRankSum", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) , 
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"DP", "_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(DP), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="DP", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) ,  
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"FS", "_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(FS), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="FS", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) ,  
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"MLEAC", "_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(MLEAC), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  geom_text(aes(label = "100", x="TP_allele_FP_allele", y = 200), color = "blue", size=3)+ 
  geom_hline(yintercept=100, color="blue", linetype="dashed", linewidth=0.5)+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="MLEAC", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) , 
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"MLEAF", "_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(MLEAF), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="MLEAF", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) ,  
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"MQ", "_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(MQ), color=eval_sample) )+
  geom_jitter(alpha=0.2, size=1)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="MQ", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) , 
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"MQRankSum", "_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(MQRankSum), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  geom_text(aes(label = "-2", x="TP_pos_TP_allele_FP_allele", y = -3), color = "blue", size=3)+ 
  geom_hline(yintercept=-2, color="blue", linetype="dashed", linewidth=0.5)+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="MQRankSum", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) ,  
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"QD", "_02_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(QD), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+ 
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  geom_hline(yintercept=2, color="blue", linetype="dashed", linewidth=0.6)+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="QD", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) , 
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"QD", "_01_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(QD), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+ 
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  geom_hline(yintercept=1, color="blue", linetype="dashed", linewidth=0.6)+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="QD", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) ,   
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"ReadPosRankSum", "_4_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(ReadPosRankSum), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  geom_hline(yintercept=4, color="blue", linetype="dashed", linewidth=0.6)+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="ReadPosRankSum", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) ,   
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"ReadPosRankSum", "_6_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(ReadPosRankSum), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  geom_hline(yintercept=6, color="blue", linetype="dashed", linewidth=0.6)+
  geom_text(aes(label = "6", x="TP", y = 6.3), color = "blue", size=3)+ 
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="ReadPosRankSum", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) , 
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();


png(file=paste0(dirI,"SOR", "3_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(SOR), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_hline(yintercept=3, color="blue", linetype="dashed", linewidth=0.6)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  geom_text(aes(label = "3", x="TP", y = 3.3), color = "blue", size=3)+ 
  # coord_cartesian(ylim = c(0.745,0.8))+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="SOR", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) ,    
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"SOR", "_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(SOR), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="SOR", x="")+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) ,   
                     axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
                     strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

png(file=paste0(dirI,"DP.1", "_distribution.png"), res=500, width=5000, height=1800);
ggplot(tpfp, 
       aes(x=eval_sample, y=as.numeric(DP.1), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="DP.1", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) ,   
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

##############------------------------ Plotting AD -------------------------################
ad <- tpfp %>%
      pivot_longer(cols = c("ad_ref", "ad_alt", "ad_alt2", "ad_alt3"),
                   names_to = c("allele"), 
                   values_to = "value",)           %>%
      drop_na()

ad<- ad[!ad$value=="", colnames(ad) %in% c("position", "region", "sample","ploidy", "eval_sample", "allele", "value")]

alt_ad<- ad[!ad$allele=="ad_ref",]

png(file=paste0(dirI,"AD", "_distribution.png"), res=500, width=5000, height=1800);
ggplot(alt_ad, 
       aes(x=eval_sample, y=as.numeric(value), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="alternative allele depth", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) , 
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();


png(file=paste0(dirI,"AD", "_2000_distribution.png"), res=500, width=5000, height=1800);
ggplot(alt_ad, 
       aes(x=eval_sample, y=as.numeric(value), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  scale_color_manual(values=palette)+
  geom_boxplot(alpha=0.5, outlier.shape = NA, color="black", width=0.4 )+
  coord_cartesian(ylim=c(0, 2000))+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="alternative allele depth", x="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25, size=8, hjust=1) ,  
        axis.title = element_text(size=8),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        strip.background = element_rect( color = "black"), panel.spacing = unit(1, "pt"))
dev.off();

# pca_result <- prcomp(tpfp[, ! colnames(tpfp) %in% c("ploidy", "region", "sample", "GT","AD", "PL", "encoded")], center = TRUE, scale. = TRUE)
# 
# library(ggplot2)
# 
# pca_df <- as.data.frame(pca_result$x) 
# pca_df$eval_sample <- eval_sample
# 
# ggplot(pca_df, aes(x = PC1, y = PC2, color = eval_sample)) +
#   geom_point() +
#   labs(x = "Principal Component 1", y = "Principal Component 2", color = "eval_sample") +
#   theme_classic()
