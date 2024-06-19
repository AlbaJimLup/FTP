library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)
## --------------------------------------  CHOOSE VC SET ---------------------------------   ############
# GATK version either "4.5.0.0" or "4.3.0.0"
set="4.3.0.0"
#set=  "4.5.0.0"
## STEP of the pipeline this is for gVCF files so either "VC" or"filtering"
step="VC"
step= "filtering"
version="v4"

dir <- "~/Desktop/ribosomal_RNAs/Alba/03_VCsimulated/v2/"

##  Palette  
palette <-  c("TP" ="#75bb6f","TP_FN_allele"="#aeb852", "TP_pos_other"="#706f4d", "FN_uncalled"="#ffc60e", "FN_ref_homozygous"= "#4c89ee","FN_alt_homozygous"="#71a5be",
              "FP"= "#bd0752", "TP_pos_FP"= "#ff4d7a", "TP_pos_FP_allele"="#ed8c96","TP_allele_FP_allele"="#ed8c96",  "err"= "#a6a6a6") 
#
## --------------------------------------  LOAD DATA ---------------------------------   ############
dirf<- paste0(dir, "VC_out/Sample_level_eval/4.3.0.0/VC/")
dirI<- paste0(dir, "IMAGES/01_PREJOIN/4.3.0.0/VC_all/")

##########################  VC ########################## 
# Dataset containing all ploidies VCs
d<- read.table(paste0(dirf, "VCeval_prejoin_", step,".tbl"), sep="\t")

compare<-data.frame()

FPa<- d[d$ploidy %in% c(5, 50) & d$eval_sample=="TP_allele_FP_allele", ]

for (i in 1:nrow(FPa)){
  row<- FPa[i , ]
  compare<- rbind(compare, d[d$ploidy==50 &  d$position==row$position & d$sample==row$sample, ])
}
compare<- rbind(compare, FPa)


comparet<- as.data.frame(table(compare$ploidy, compare$eval_sample))
names(comparet)<- c("ploidy", "eval_sample", "counts")
comparet<- comparet[comparet$counts!=0, ]


png(file=paste0(dirI,"5_vs_50.png"), res=180, width=900, height=1000);
ggplot(comparet, aes(x=as.factor(ploidy), y=counts, fill=eval_sample, label=counts)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=palette) +
  geom_text_repel(position="stack")+
  labs(x="Ploidy", y="Counts", fill="Evaluation") +
  theme_minimal()
dev.off();
############## VC 4##############
d<- read.table(paste0(dirf, "All_VCeval_prejoin_", step,"_", version, ".tbl"), sep="\t")

compare<-data.frame()

FPa<- d[d$ploidy %in% c(5, 50) & d$eval_sample=="TP_allele_FP_allele", ]

for (i in 1:nrow(FPa)){
  row<- FPa[i , ]
  compare<- rbind(compare, d[d$ploidy==50 &  d$position==row$position & d$sample==row$sample, ])
}
compare<- rbind(compare, FPa)


comparet<- as.data.frame(table(compare$ploidy, compare$eval_sample))
names(comparet)<- c("ploidy", "eval_sample", "counts")
comparet<- comparet[comparet$counts!=0, ]


png(file=paste0(dirI,"5_vs_50_v4.png"), res=180, width=900, height=1000);
ggplot(comparet, aes(x=as.factor(ploidy), y=counts, fill=eval_sample, label=counts)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=palette) +
  geom_text_repel(position="stack")+
  labs(x="Ploidy", y="Counts", fill="Evaluation") +
  theme_minimal()
dev.off();

## --------------------------------------  Ploidy 50 Samples DATA ---------------------------------   ############

# Dataset containing all ploidies VCs
dirI<-paste0( dirI, "Analysis/")
d1<- read.table(paste0(dirf, "VCeval_prejoin_p50_filtering_2sample_v4.tbl"), sep="\t")
d1$ploidy<- paste(d1$ploidy, "ploidy 5 samples")
d2<- read.table(paste0(dirf, "VCeval_prejoin_p50_filtering_15sample_v4.tbl"), sep="\t")
d2$ploidy<- paste(d2$ploidy, "ploidy 10 samples")

d <- rbind(d1, d2)

d$eval_sample<-  factor(d$eval_sample, levels=c("FN_alt_homozygous", "FN_ref_homozygous", "FN_uncalled","TP_FN_allele" ,"TP_pos_FN_allele", "FP","TP_pos_FP", "TP_allele_FP_allele", "TP"))

d_plot<- d[d$eval_sample!="FN_uncalled",]
########------------------------------- Get AD to frequency columns -----------------------------------########
# d$eval_sample[d$eval_sample=="TP_FN_allele" ]<- "TP_pos_FN_allele"

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
########------------------------------- PREP PLOTTING -----------------------------------########
# write.table(cbind(d, alleles_freqs), paste0(dirf, "VCeval_prejoin_alleles_", step,".tbl"), sep="\t")
# d<- read.table( paste0(dirf, "VCeval_prejoin_alleles_", step,".tbl"), sep="\t")

## Add columns alternative alleles frequencies  in Long format
d_long<- pivot_longer(cbind(d, alleles_freqs), paste0("alt", 1:(max_a-1)), names_to = "allele", values_to = "AF")
## Remove empty rows
d_long<- d_long[!is.na(d_long$AF), ]

print("Done preparing data")
########------------------------------- Plotting alternative alleles frequency -----------------------------------########

png(file=paste0(dirI,"AF_alt_boxplot.png"), res=500, width=3000, height=2500);
ggplot(d_long,  aes(x=eval_sample, y=as.numeric(AF), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  geom_boxplot(alpha=0.1, color="black",  outlier.shape = NA)+
  scale_color_manual(values=palette)+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="alternative allele depth", x="")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8, angle=45, vjust=0.5))
dev.off();

png(file=paste0(dirI,"AF_alt_boxplot_zoom.png"), res=500, width=3000, height=2000);
ggplot(d_long,  aes(x=eval_sample, y=as.numeric(AF), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  geom_boxplot(alpha=0.1, color="black",  outlier.shape = NA)+
  coord_cartesian(ylim = c(0, 0.18))+
  scale_color_manual(values=palette)+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="alternative allele depth", x="")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8, angle=45, vjust=0.5))
dev.off();

png(file=paste0(dirI,"AF_ref_boxplot.png"), res=500, width=3000, height=2500);
ggplot(d_long,  aes(x=eval_sample, y=as.numeric(ref), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  geom_boxplot(alpha=0.1, color="black",  outlier.shape = NA)+
  scale_color_manual(values=palette)+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="reference allele depth", x="")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8, angle=45, vjust=0.5))
dev.off();

png(file=paste0(dirI,"AF_ref_boxplot_zoom.png"), res=500, width=3000, height=2000);
ggplot(d_long,  aes(x=eval_sample, y=as.numeric(ref), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  geom_boxplot(alpha=0.1, color="black",  outlier.shape = NA)+
  geom_hline(yintercept=0.99, color="blue", linetype="dashed", linewidth=0.6)+
  geom_text(aes(label = "0.99", x="TP", y = 0.999), color = "blue", size=4)+
  scale_color_manual(values=palette)+
  coord_cartesian(ylim = c(0.9, 1))+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="reference allele depth", x="")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8, angle=45, vjust=0.5))
dev.off();

########------------------------------- Plotting variables per position -----------------------------------########
png(file=paste0(dirI,"QUAL", "_distribution.png"), res=500, width=2500, height=2500);
ggplot(d_plot, aes(x=eval_sample, y=as.numeric(QUAL), color=eval_sample) )+
  # coord_cartesian(ylim = c(0,10000))+
  geom_jitter(alpha=0.3, size=1)+
  geom_boxplot(alpha=0.3, size=0.5, color="black", width = 0.4, outlier.shape = NA)+
  scale_color_manual(values=palette)+
  labs(y="quality", x="")+
  facet_grid(.~ploidy)+
  guides(color="none")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, vjust=0.5))
dev.off();

png(file=paste0(dirI,"QUAL", "10000_distribution.png"), res=500, width=2500, height=2500);
ggplot(d_plot, aes(x=eval_sample, y=as.numeric(QUAL), color=eval_sample) )+
  # coord_cartesian(ylim = c(0,10000))+
  coord_cartesian(ylim = c(0,100000))+
  geom_hline(yintercept=1000, color="blue", linetype="dashed", linewidth=0.6)+
  geom_text(aes(label = "1000", x="TP_pos_FN_allele", y = 2000), color = "blue", size=4)+
  geom_jitter(alpha=0.3, size=1)+
  geom_boxplot(alpha=0.3, size=0.5, color="black", width = 0.4, outlier.shape = NA)+
  scale_color_manual(values=palette)+
  labs(y="quality", x="")+
  facet_grid(.~ploidy)+
  guides(color="none")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, vjust=0.5))
dev.off();

# find possible filter know more about the test: https://gatk.broadinstitute.org/hc/en-us/articles/360042477772-BaseQualityRankSumTest
png(file=paste0(dirI,"BaseQRankSum", "_distribution.png"), res=700, width=2500, height=2500);
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
        strip.background = element_rect(fill = "lightgray", color = "black"),
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
png(file=paste0(dirI,"MQRankSum", "_distribution.png"), res=500, width=2500, height=2500);
ggplot(d_plot, aes(x=eval_sample, y=as.numeric(MQRankSum), color=eval_sample) )+
  geom_jitter(alpha=0.3, size=1)+
  geom_boxplot(alpha=0.1, color="black",  outlier.shape = NA)+
  scale_color_manual(values=palette)+
  facet_grid(.~ploidy)+
  guides(color="none")+
  labs(y="MQRankSum", x="", color="")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, vjust=0.5))
dev.off();

png(file=paste0(dirI,"ReadPosRankSum", "_4_distribution.png"), res=500, width=2500, height=2500);
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
        strip.background = element_rect(fill = "lightgray", color = "black"),
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "#7d7d7d", size = 0.5))
dev.off();

png(file=paste0(dirI,"ReadPosRankSum", "_10_distribution.png"), res=500, width=2500, height=2500);
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
        strip.background = element_rect(fill = "lightgray", color = "black"),
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "#7d7d7d", size = 0.5))
dev.off();

png(file=paste0(dirI,"DP.1", "_distribution.png"), res=500, width=2500, height=2500);
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
        strip.background = element_rect(fill = "lightgray", color = "black"),
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "#7d7d7d", size = 0.5))
dev.off();


png(file=paste0(dirI,"GQ", "_distribution.png"), res=500, width=2500, height=2500);
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
        strip.background = element_rect(fill = "lightgray", color = "black"),
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "#7d7d7d", size = 0.5))
dev.off();

png(file=paste0(dirI,"GQ", "_zoom_distribution.png"), res=500, width=2500, height=1500);
ggplot(d_plot, aes(x=eval_sample, y=GQ, color=eval_sample) )+
  geom_jitter(alpha=0.6, size=0.5)+
  geom_boxplot(alpha=0.1, color="black",  width=0.3, outlier.shape = NA)+
  scale_color_manual(values=palette)+
  coord_cartesian(ylim=c(96, 99.5))+
  facet_grid(.~ploidy)+
  guides(color = guide_legend(nrow = 1, override.aes = list(size=3)))+
  labs(y="GQ", x="", color="")+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.margin = margin(-15, -10, 0, 0),
        strip.text = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        strip.background = element_rect(fill = "lightgray", color = "black"),
        strip.placement = "outside",
        panel.spacing = unit(0, "pt"),
        panel.border = element_rect(color = "#7d7d7d", size = 0.5))
dev.off();
print("All done!")
print("#############################################################################")