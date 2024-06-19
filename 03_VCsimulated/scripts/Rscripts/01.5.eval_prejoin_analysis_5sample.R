## --------------------------------------  CHOOSE VC SET ---------------------------------   ############
# args <- commandArgs(trailingOnly = TRUE)
# ## GATK version either "4.5.0.0" or "4.3.0.0"
# set=args[1]
# #set=  "4.3.0.0"
# ## STEP of the pipeline this is for gVCF files so either "VC" or"filtering"
# step=args[2]
# #step= "filtering"
# dir <-  "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/"

# GATK version either "4.5.0.0" or "4.3.0.0"
set="4.3.0.0"
#set=  "4.5.0.0"
## STEP of the pipeline this is for gVCF files so either "VC" or"filtering"
step="VC"
step= "filtering"
version="v4"

dir <- "~/Desktop/ribosomal_RNAs/Alba/03_VCsimulated/v2/"

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)


##  Palette  
palette <-  c("TP" ="#75bb6f","TP_FN_allele"="#aeb852", "TP_pos_other"="#706f4d", "FN_uncalled"="#ffc60e", "FN_ref_homozygous"= "#4c89ee","FN_alt_homozygous"="#71a5be",
              "FP"= "#bd0752", "TP_pos_FP"= "#ff4d7a", "TP_pos_FP_allele"="#ed8c96","TP_allele_FP_allele"="#ed8c96",  "err"= "#a6a6a6") 
#
## --------------------------------------  LOAD DATA ---------------------------------   ############
dirf<- paste0(dir, "VC_out/Sample_level_eval/4.3.0.0/VC/")
dirI<- paste0(dir, "IMAGES/01_PREJOIN/4.3.0.0/VC/")

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


