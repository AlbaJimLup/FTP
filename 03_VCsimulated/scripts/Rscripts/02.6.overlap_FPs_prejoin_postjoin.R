## --------------------------------------  LOAD DATA ---------------------------------   ############
#dir <-  "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/pre_join_gen/"
dir <- "~/Desktop/ribosomal_RNAs/Alba/03_VCsimulated/v2/"
dirf<- paste0(dir, "VC_out/Sample_level_eval/default_4.3.0.0/")
dirt <- paste0(dir, "scripts/Rscripts/")
dirI<- paste0(dir, "IMAGES/03_POSTJOIN/default_4.3.0.0/Analysis/")

prejoin <- read.table( paste0(dirf,"VCeval_prejoin.tbl"), sep="\t")
postjoin <- read.table( paste0(dirf,"VCeval_postjoin.tbl"), sep="\t")

true_variants_sample <- read.table(paste0(dirt, "true_variants.tbl"), header=T)
catalog <-  read.table(paste0(dirt, "variants_AF_Counts.tbl"))
###### Libraries ######
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)
####  Palettes  #####
palette <-  c("TP" ="#75bb6f","TP_FN_allele"="#aeb852",    "TP_pos_other"="#706f4d", 
              "FN_uncalled"="#ffc60e", "FN_ref_homozygous"= "#4c89ee","FN_alt_homozygous"="#71a5be",
              "FP"= "#bd0752", "TP_pos_FP"= "#ff4d7a", "TP_pos_FP_allele"="#ed8c96","TP_allele_FP_allele"="#ed8c96",  "err"= "#a6a6a6") 

palette2 <- c( "5_ETS" = "#2f4b5b",   "18S" = "#ffe290",  "ITS1" = "#9dbccf",
               "5.8S" = "#fdae61",  "ITS2" = "#4d8195",  "28S" = "#ff6e7d",  "3_ETS" = "#375794")
############## 

prejoin_q1<- prejoin[prejoin$QUAL < 1 & prejoin$eval_sample !="FN_uncalled", ]
postjoin<- postjoin[postjoin$eval_sample!="FN_uncalled", ]

t<-  as.data.frame(matrix(ncol=5, nrow = 400))
colnames(t)<- c("sample", "ploidy", "prejoin", "postjoin", "both")
t$sample<- rep(paste0("sample_", 1:100), 4)
t$ploidy<- c(rep(2, 100), rep(5, 100), rep(7, 100), rep(10, 100))

for (p in c(2, 5, 7, 10)){
  for (s in paste0("sample_", 1:100)){
    pre_q1pos<-  prejoin_q1$position[prejoin_q1$ploidy==p & prejoin_q1$sample==s]
    # Checking if they are coing it wrongly
    post_pos <- postjoin$position[postjoin$ploidy == p & postjoin$sample == s & postjoin$eval_sample != "TP"]
    # Only quality < 1 and wrongly called in prejoin
    t$prejoin[t$ploidy==p & t$sample==s] <- sum(setdiff(pre_q1pos, post_pos))
    # Only wrongly called in postjoin
    t$postjoin[t$ploidy==p & t$sample==s] <- sum(setdiff(post_pos, pre_q1pos))
    # Wrongly in both
    t$both<- sum(intersect(pre_q1pos, post_pos))
    
    print(paste("Done with sample", s))
  }
}
