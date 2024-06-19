# ## --------------------------------------  CHOOSE VC SET ---------------------------------   ############
# # GATK version either "4.5.0.0" or "4.3.0.0"
# set="4.3.0.0"
# #set=  "4.5.0.0"
# ## STEP of the pipeline this is for gVCF files so either "VC" or"filtering"
# step="VC"
# # step= "filtering"
# 
#dir <- "~/Desktop/ribosomal_RNAs/Alba/03_VCsimulated/v2/"
# ## --------------------------------------  LOAD DATA ---------------------------------   ############
dir<- "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/"
# # Dataset containing all ploidies VCs
# # d<- read.table(paste0(dirf, "VCeval_prejoin_", step,".tbl"), sep="\t")
# # d<- dfs
# # dirI<- paste0(dirI,"Analysis/" )
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
#CHOOSE Join genotyping SET
args <- commandArgs(trailingOnly = TRUE)

#GATK version either "4.5.0.0" or "4.3.0.0"
set=args[1]
# set= "4.3.0.0"

# STEP of the pipeline this is for gVCF files so either "VC" or"filtering"
step=args[2]
# step= "VC"
# step= "filtering"

# version
version=args[3]
# version="v1"
# version="v1"
# version="v3"
# version="v4"
# version="5sample_v4"
# version="10sample_v4"
# version="5sample_v3"
# version="10sample_v3"
# version="5sample_v2"
# version="5sample_v3"
# version="10sample_v2"
# version="2sample_v4"
# version="2sample_v2"
# version="2sample_v3"
# #GATK version either "4.5.0.0" or "4.3.0.0"
## --------------------------------------  LOAD DATA ---------------------------------   ############

# Dataset containing all ploidies VCs
if(step=="filtering"){
  dirf<- paste0(dir, "VC_out/Sample_level_eval/", set, "/VC/")
  d<- read.table(  paste0(dirf, "All_VCeval_prejoin_filtering_", version, ".tbl"), sep="\t")
  if (version %in% c("v1", "v2", "v3", "v4")){
    dirI<- paste0(dir, "IMAGES/01_PREJOIN/", set, "/filtering_", version, "_all/")
  }else{
    dirI<- paste0(dir, "IMAGES/01_PREJOIN/", set, "/filtering_sample/filtering_", version, "_all/")
  }
}else{
  dirf<- paste0(dir, "VC_out/Sample_level_eval/", set, "/VC/")
  dirI<- paste0(dir, "IMAGES/01_PREJOIN/", set, "/VC_all/")
  d<- read.table(  paste0(dirf, "All_VCeval_prejoin_VC.tbl"), sep="\t")
}


print(paste("We got file", set, step))


# d<- dfs
dirI<- paste0(dirI, "Analysis/")
setwd(dirI)

# ploidies<-  c(2,5,7,10,20,30,40,50,100)
ploidies<- unique(d$ploidy)

print("#############################################################################")
print(paste("PLOTTING VC analysis for", set, step, ploidies))

########------------------------------- Get AD to frequency columns -----------------------------------########
d$eval_sample[d$eval_sample=="TP_FN_allele" ]<- "TP_pos_FN_allele"

d$eval_sample<-  factor(d$eval_sample, levels=c("FN_alt_homozygous", "FN_ref_homozygous", "FN_uncalled","TP_pos_FN_allele", "FP","TP_pos_FP", "TP_allele_FP_allele", "TP"))

nrows<- nrow(d)
d$AFs<- ""
numcols<- c()

# Get vector with AFs
for (i in 1: nrows){
  # GET only the alleles that were called in the genotype
  encoded <- strsplit(d$data[i], ":")[[1]][1]

  if (grepl("/", encoded)) {# using "/" delimiter
    alleles <- unlist(strsplit(encoded, "/"))

  } else if (grepl("\\|", encoded)) {# using "|" delimiter
    alleles <- unlist(strsplit(encoded, "\\|"))
  }
  AD <- as.numeric(unlist(strsplit(d$AD[i], ","))) # get all alleles
  # select only the alleles called
  ADs<- c()
  for (allele in unique(alleles)){
    ADs<-  c(ADs, AD[as.numeric(allele)+1])
  }
  # Generate AF vector
  d$AFs[i] <- paste(round(ADs/sum(ADs), 6),collapse = ",")
  if (length(ADs)>3) numcols<- c(numcols, length(ADs))
  print(paste0("Done ", i*100/nrows, " %"))
}
# Maximum number of alleles
max_a<- max(numcols)

# Generate the columns for AF alternative alleles
alleles_freqs<- data.frame(matrix(ncol=max_a, nrow=nrows))
colnames(alleles_freqs)<-  c("ref", paste0("alt", 1:(max_a-1)))

for (i in 1:nrows){
  AFs<- unlist(strsplit(d$AF[i], ",")) # get all alleles frequencies
  if (!NA %in% AFs ){
    for (col in 1:max_a){ # alternative alleles frequency columns
      alleles_freqs[i, col] <- as.numeric(AFs[col])
    }
  }
  print(paste0("Done ", i*100/nrows, " %"))
}

print("Done preping data")
#######------------------------------- PREP PLOTTING -----------------------------------########
# write.table(cbind(d, alleles_freqs), paste0(dirf, "VCeval_prejoin_alleles_", step, version,".tbl"), sep="\t")
# d<- read.table( paste0(dirf, "VCeval_prejoin_alleles_", step,".tbl"), sep="\t")
print("Done writing data frame")

## Add columns alternative alleles frequencies  in Long format
d_long<- pivot_longer(cbind(d, alleles_freqs), paste0("alt", 1:(max_a-1)), names_to = "allele", values_to = "AF")
## Remove empty rows
d_long<- d_long[!is.na(d_long$AF), ]

d_plot<- d[d$eval_sample!="FN_uncalled",]

print("Done preparing data")
########------------------------------- Plotting alternative alleles frequency -----------------------------------########
library(ggplot2)

p1 <- ggplot(d_long,  aes(x=eval_sample, y=as.numeric(AF), color=eval_sample) ) +
  geom_jitter(alpha=0.3, size=1) +
  geom_boxplot(alpha=0.1, color="black",  outlier.shape = NA) +
  scale_color_manual(values=palette) +
  facet_grid(.~ploidy) +
  guides(color="none") +
  labs(y="alternative allele depth", x="") +
  theme_classic() +
  theme(axis.text.x = element_text(size=8, angle=45, vjust=0.5), 
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"),
        panel.grid.minor = element_line(colour = "#e6e6e6"))

ggsave(filename=paste0(dirI,"AF_alt_boxplot.png"), plot=p1, dpi=500, width=8000/500, height=2200/500)

p2 <- ggplot(d_long,  aes(x=eval_sample, y=as.numeric(AF), color=eval_sample) ) +
  geom_jitter(alpha=0.3, size=1) +
  geom_boxplot(alpha=0.1, color="black",  outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 0.18)) +
  scale_color_manual(values=palette) +
  facet_grid(.~ploidy) +
  guides(color="none") +
  labs(y="alternative allele depth", x="") +
  theme_classic() +
  theme(axis.text.x = element_text(size=8, angle=45, vjust=0.5),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"),
        panel.grid.minor = element_line(colour = "#e6e6e6"))

ggsave(filename=paste0(dirI,"AF_alt_boxplot_zoom.png"), plot=p2, dpi=500, width=8000/500, height=2000/500)

p3 <- ggplot(d_long,  aes(x=eval_sample, y=as.numeric(ref), color=eval_sample) ) +
  geom_jitter(alpha=0.3, size=1) +
  geom_boxplot(alpha=0.1, color="black",  outlier.shape = NA) +
  scale_color_manual(values=palette) +
  facet_grid(.~ploidy) +
  guides(color="none") +
  labs(y="reference allele depth", x="") +
  theme_classic() +
  theme(axis.text.x = element_text(size=8, angle=45, vjust=0.5),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(colour = "#e6e6e6"),
        panel.grid.minor = element_line(colour = "#e6e6e6"))

ggsave(filename=paste0(dirI,"AF_ref_boxplot.png"), plot=p3, dpi=500, width=8000/500, height=2000/500)

p4 <- ggplot(d_long,  aes(x=eval_sample, y=as.numeric(ref), color=eval_sample) ) +
  geom_jitter(alpha=0.3, size=1) +
  geom_boxplot(alpha=0.1, color="black",  outlier.shape = NA) +
  geom_hline(yintercept=0.99, color="blue", linetype="dashed", linewidth=0.6) +
  geom_text(aes(label = "0.99", x="TP", y = 0.999), color = "blue", size=4) +
  scale_color_manual(values=palette) +
  coord_cartesian(ylim = c(0.9, 1)) +
  facet_grid(.~ploidy) +
  guides(color="none") +
  labs(y="reference allele depth", x="") +
  theme_classic() +
  theme(axis.text.x = element_text(size=8, angle=45, vjust=0.5), 
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.minor = element_line(colour = "#e6e6e6"))

ggsave(filename=paste0(dirI,"AF_ref_boxplot_zoom.png"), plot=p4, dpi=500, width=8000/500, height=2000/500)

########------------------------------- Plotting QUAL -----------------------------------########
p1 <- ggplot(d_plot, aes(x=eval_sample, y=as.numeric(QUAL), color=eval_sample)) +
  geom_jitter(alpha=0.3, size=1) +
  geom_boxplot(alpha=0.3, size=0.5, color="black", width=0.4, outlier.shape=NA) +
  scale_color_manual(values=palette) +
  labs(y="quality", x="") +
  facet_grid(.~ploidy) +
  guides(color="none") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust=0.5),  
        strip.text = element_text(size=10, face="bold"),
        panel.border = element_rect(color="black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour="#e6e6e6"))
ggsave(filename=paste0(dirI, "QUAL_distribution.png"), plot=p1, dpi=500, width=8500/500, height=2000/500, units="in")

p2 <- ggplot(d_plot, aes(x=eval_sample, y=as.numeric(QUAL), color=eval_sample)) +
  coord_cartesian(ylim=c(0, 100000)) +
  geom_hline(yintercept=1000, color="blue", linetype="dashed", linewidth=0.6) +
  geom_text(aes(label="1000", x="TP_pos_FN_allele", y=2000), color="blue", size=4) +
  geom_jitter(alpha=0.3, size=1) +
  geom_boxplot(alpha=0.3, size=0.5, color="black", width=0.4, outlier.shape=NA) +
  scale_color_manual(values=palette) +
  labs(y="quality", x="") +
  facet_grid(.~ploidy) +
  guides(color="none") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust=0.5), 
        strip.text = element_text(size=10, face="bold"),
        panel.border = element_rect(color="black", fill=NA, size=0.5),
        panel.grid.minor = element_line(colour="#e6e6e6"),
        panel.grid.major = element_line(colour="#e6e6e6"))
ggsave(filename=paste0(dirI, "QUAL_10000_distribution.png"), plot=p2, dpi=500, width=8500/500, height=2000/500, units="in")

p3 <- ggplot(d_plot, aes(x=eval_sample, y=as.numeric(QUAL), color=eval_sample)) +
  coord_cartesian(ylim=c(0, 2000)) +
  geom_jitter(alpha=0.8, size=0.4) +
  geom_boxplot(alpha=0.1, size=0.5, color="black", width=0.4, outlier.shape=NA) +
  scale_color_manual(values=palette) +
  geom_hline(yintercept=1000, color="blue", linetype="dashed", linewidth=0.6) +
  geom_text(data=subset(d_plot, d_plot$ploidy==2), aes(label="1000", x="FN_ref_homozygous", y=1150), color="blue", size=4) +
  geom_hline(yintercept=400, color="blue", linetype="dashed", linewidth=0.6) +
  geom_text(data=subset(d_plot, d_plot$ploidy==2), aes(label="400", x="FN_ref_homozygous", y=550), color="blue", size=4) +
  geom_hline(yintercept=30, color="blue", linetype="dashed", linewidth=0.6) +
  geom_text(data=subset(d_plot, d_plot$ploidy==2), aes(label="30", x="FN_ref_homozygous", y=100), color="blue", size=4) +
  labs(y="quality", x="") +
  facet_grid(.~ploidy) +
  guides(color="none") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust=0.5),
        strip.text = element_text(size=12, face="bold"),
        panel.grid.minor = element_line(colour="#e6e6e6"), 
        panel.border = element_rect(color="black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour="#e6e6e6"))
ggsave(filename=paste0(dirI, "QUAL_2000_distribution.png"), plot=p3, dpi=520, width=8500/520, height=2000/520, units="in")


p4 <- ggplot(d_plot, aes(x=eval_sample, y=as.numeric(QUAL), color=eval_sample)) +
  coord_cartesian(ylim=c(0, 900)) +
  geom_jitter(alpha=0.3, size=0.6) +
  geom_boxplot(alpha=0.3, size=0.5, color="black", width=0.4, outlier.shape=NA) +
  scale_color_manual(values=palette) +
  geom_hline(yintercept=400, color="blue", linetype="dashed", linewidth=0.6) +
  geom_text(aes(label="400", x="FN_alt_homozygous", y=450), color="blue", size=5) +
  geom_hline(yintercept=30, color="blue", linetype="dashed", linewidth=0.6) +
  geom_text(aes(label="30", x="FN_alt_homozygous", y=90), color="blue", size=5) +
  labs(y="quality", x="") +
  facet_grid(.~ploidy) +
  guides(color="none") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust=0.5),
        strip.text = element_text(size=10, face="bold"),
        panel.grid.minor = element_line(colour="#e6e6e6"), 
        panel.border = element_rect(color="black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour="#e6e6e6"))
ggsave(filename=paste0(dirI, "QUAL_400_distribution.png"), plot=p4, dpi=500, width=8500/500, height=2000/500, units="in")


p5 <- ggplot(d_plot, aes(x=eval_sample, y=as.numeric(QUAL), color=eval_sample)) +
  coord_cartesian(ylim=c(0, 100)) +
  geom_jitter(alpha=0.8, size=1) +
  geom_boxplot(alpha=0.3, size=0.5, color="black", width=0.4, outlier.shape=NA) +
  scale_color_manual(values=palette) +
  geom_hline(yintercept=30, color="blue", linetype="dashed", linewidth=0.6) +
  geom_text(aes(label="30", x="FN_alt_homozygous", y=40), color="blue", size=5) +
  labs(y="quality", x="") +
  facet_grid(.~ploidy) +
  guides(color="none") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust=0.5),
        strip.text = element_text(size=10, face="bold"),
        panel.grid.minor = element_line(colour="#e6e6e6"),  
        panel.border = element_rect(color="black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour="#e6e6e6"))
ggsave(filename=paste0(dirI, "QUAL_30_100_distribution.png"), plot=p5, dpi=500, width=8500/500, height=2000/500, units="in")

p6 <- ggplot(d_plot, aes(x=eval_sample, y=as.numeric(QUAL), color=eval_sample)) +
  coord_cartesian(ylim=c(0, 60)) +
  geom_jitter(alpha=0.8, size=1) +
  geom_boxplot(alpha=0.3, size=0.5, color="black", width=0.4, outlier.shape=NA) +
  scale_color_manual(values=palette) +
  geom_hline(yintercept=30, color="blue", linetype="dashed", linewidth=0.6) +
  geom_text(aes(label="30", x="FN_alt_homozygous", y=40), color="blue", size=5) +
  labs(y="quality", x="") +
  facet_grid(region~ploidy) +
  guides(color="none") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust=0.5, size=12),
        panel.grid.major = element_line(colour="#e6e6e6"), 
        panel.border = element_rect(color="black", fill=NA, size=0.5),
        strip.text = element_text(size=10, face="bold"),
        panel.grid.minor = element_line(colour="#e6e6e6"))
ggsave(filename=paste0(dirI, "QUAL_30_60_distribution.png"), plot=p6,dpi=500, width=8500/500, height=8000/500, units="in")

########------------------------------- Plotting other variables -----------------------------------########
p<- ggplot(d_plot, aes(x=eval_sample, y=as.numeric(BaseQRankSum), color=eval_sample)) +
  geom_jitter(alpha=0.6, size=0.4, width = 0.4) +
  geom_boxplot(alpha=0.3, size=0.5, color="black", width = 0.4, outlier.shape = NA) +
  scale_color_manual(values=palette) +
  facet_grid(.~ploidy) +
  labs(y="BaseQRankSum", x="", color="Sample level evaluation") +
  guides(color = guide_legend(nrow = 1, override.aes = list(size=3))) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.margin = margin(-22, -10, 0, 0),
        strip.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        panel.grid.major = element_line(colour = "#e6e6e6"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),
        strip.background = element_rect(fill = NA, color = "black"),
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))
ggsave(filename=paste0(dirI, "BaseQRankSum", "_distribution.png"), plot=p, dpi=500, width=12, height=2.5, units="in")

# ExcessHet plot
p<- ggplot(d_plot, aes(x=eval_sample, y=as.numeric(ExcessHet), color=eval_sample)) +
  geom_jitter(alpha=0.3,  size=0.5) +
  scale_color_manual(values=palette) +
  facet_grid(.~ploidy) +
  labs(y="ExcessHet", x="", color="Sample level evaluation") +
  guides(color = guide_legend(nrow = 1, override.aes = list(size=3))) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.margin = margin(-22, -10, 0, 0),
        strip.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        strip.background = element_rect(fill = NA, color = "black"),
        panel.grid.major = element_line(colour = "#e6e6e6"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "black",fill = NA,  size = 0.5))
ggsave(filename=paste0(dirI, "ExcessHet", "_distribution.png"),plot=p, dpi=500, width=12, height=2.5, units="in")

# MQRankSum plot
p<- ggplot(d_plot, aes(x=eval_sample, y=as.numeric(MQRankSum), color=eval_sample)) +
  geom_jitter(alpha=0.7,  size=0.5) +
  geom_boxplot(alpha=0.1, color="black",  outlier.shape = NA) +
  scale_color_manual(values=palette) +
  facet_grid(.~ploidy) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size=3))) +
  labs(y="MQRankSum", x="", color="Sample level evaluation") +
  theme_classic() +
  theme(axis.text.x = element_blank(), legend.position = "bottom", 
        panel.grid.major = element_line(colour = "#e6e6e6"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.margin = margin(-22, -10, 0, 0), axis.ticks = element_blank())
ggsave(filename=paste0(dirI, "MQRankSum", "_distribution.png"),plot=p, dpi=500, width=12, height=2.5, units="in")

# ReadPosRankSum (4) plot
p<- ggplot(d_plot, aes(x=eval_sample, y=as.numeric(ReadPosRankSum), color=eval_sample)) +
  geom_jitter(alpha=0.7,  size=0.5) +
  scale_color_manual(values=palette) +
  geom_boxplot(alpha=0.3, size=0.5, color="black", width = 0.4, outlier.shape = NA) +
  geom_hline(yintercept=4, color="blue", linetype="dashed", linewidth=0.6) +
  geom_text(aes(label = "4", x="FN_alt_homozygous", y=7), color = "blue", size=5) +
  facet_grid(.~ploidy) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size=3))) +
  labs(y="ReadPosRankSum", x="", color="Sample level evaluation") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.margin = margin(-22, -10, 0, 0),
        strip.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        strip.background = element_rect(fill = NA,  color = "black"),
        panel.grid.major = element_line(colour = "#e6e6e6"),
        panel.grid.minor = element_line(colour = "#e6e6e6"), 
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "black",fill = NA,  size = 0.5))
ggsave(filename=paste0(dirI, "ReadPosRankSum", "_4_distribution.png"), plot=p,dpi=500, width=12, height=2.5, units="in")

# ReadPosRankSum (10) plot
p<- ggplot(d_plot, aes(x=eval_sample, y=as.numeric(ReadPosRankSum), color=eval_sample)) +
  geom_jitter(alpha=0.7,  size=0.5) +
  scale_color_manual(values=palette) +
  geom_boxplot(alpha=0.3, size=0.5, color="black", width = 0.4, outlier.shape = NA) +
  geom_hline(yintercept=10, color="blue", linetype="dashed", linewidth=0.6) +
  geom_text(aes(label = "10", x="FN_alt_homozygous", y = 12), color = "blue", size=5) +
  facet_grid(.~ploidy) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size=3))) +
  labs(y="ReadPosRankSum", x="",color="Sample level evaluation") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.margin = margin(-22, -10, 0, 0),
        strip.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        strip.background = element_rect(fill = NA, color = "black"),
        panel.grid.major = element_line(colour = "#e6e6e6"),
        panel.grid.minor = element_line(colour = "#e6e6e6"), 
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "black",fill = NA,  size = 0.5))
ggsave(filename=paste0(dirI, "ReadPosRankSum", "_10_distribution.png"),plot=p, dpi=500, width=12, height=2.5, units="in")

# DP.1 plot
p<- ggplot(d_plot, aes(x=eval_sample, y=as.numeric(DP.1), color=eval_sample)) +
  geom_jitter(alpha=0.6,  size=0.5) +
  scale_color_manual(values=palette) +
  geom_boxplot(alpha=0.3, size=0.5, color="black", width = 0.4, outlier.shape = NA) +
  facet_grid(.~ploidy) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size=3))) +
  labs(y="DP.1", x="", color="Sample level evaluation") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.margin = margin(-22, -10, 0, 0),
        strip.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        strip.background = element_rect(fill = NA, color = "black"),
        panel.grid.major = element_line(colour = "#e6e6e6"), 
        panel.grid.minor = element_line(colour = "#e6e6e6"),
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5))
ggsave(filename=paste0(dirI, "DP.1", "_distribution.png"),plot=p, dpi=500, width=12, height=2.5, units="in")

# GQ plot
p<- ggplot(d_plot, aes(x=eval_sample, y=GQ, color=eval_sample)) +
  geom_jitter(alpha=0.6, size=0.5) +
  geom_boxplot(alpha=0.1, color="black",  width=0.3, outlier.shape = NA) +
  scale_color_manual(values=palette) +
  facet_grid(.~ploidy) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size=3))) +
  labs(y="GQ", x="", color="Sample level evaluation") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.margin = margin(-22, -10, 0, 0),
        strip.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        strip.background = element_rect(fill = NA, color = "black"),
        panel.grid.major = element_line(colour = "#e6e6e6"),
        panel.grid.minor = element_line(colour = "#e6e6e6"),
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "black",fill = NA,  size = 0.5))
ggsave(filename=paste0(dirI, "GQ", "_distribution.png"),plot=p,  dpi=500, width=12, height=2.5, units="in")

print("All done!")
print("#############################################################################")

########------------------------------- Varables that need to be split before plotting -----------------------------------########
# # "DP"
# #"MLEAC"
# # "MLEAF"
# # RAW_MQandDP
# #AD
# #PL
# #SB
