#  Libraries 
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)
###### Functions######
get_metrics<- function(d){
  # Initialize recount table
  score <-  data.frame(matrix(nrow=0, ncol=4))
  colnames(score)<-  c("ploidy", "precision", "recall", "F1_score")
  
  for (ploidy in unique(d$ploidy)){
    FP <- 0
    TP <- 0
    FN <- 0
    rows <- d[d$ploidy==ploidy,]
    
    for (i in 1:nrow(rows)){ # Get estimates
      row<- rows[i,]
      if(row$eval == "FN_uncalled")      FN <- FN + as.numeric(row$counts)
      else if ( row$eval =="TP")         TP <- TP + as.numeric(row$counts)
      else    FP <-  FP+ as.numeric(row$counts)
    }
    # accuracy vs mistakes 
    precision <-   TP/(TP+FP) 
    # how good at finding
    recall <-  TP/(TP+FN)  
    
    score<- rbind(score, c(ploidy,  precision, recall,   2 * precision * recall / (precision +recall)))
  }
  return(score)
}
##########
step <- "filtering"

dir <-  "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/Sample_level_eval/4.3.0.0/VC/"
# dir <- paste0("~/Desktop/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/Sample_level_eval/4.3.0.0/VC/")
# dirf<- "~/Desktop/ribosomal_RNAs/Alba/07_Performance/"
dirf<-"/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/07_Performance/"
#dir <-  "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/pre_join_gen/"

for (v in 1:4){
  ## --------------------------------------  LOAD DATA ---------------------------------############
  vc<-  read.table( paste0(dir,"All_VCeval_prejoin_", step, "_v", v ,".tbl"), sep="\t")# after all analysis Rscript files
  
  ####### ---------------------------------- DATA PREP -----------------------------############
  vc<- vc[ vc$class =="callable", ]#  KEEP ONLY CALLABLE
  vc$eval_sample[ vc$eval_sample =="TP_pos_TP_allele_FP_allele"]<- "TP_allele_FP_allele"
  ####### ---------------------------------- COUNT EVALS -----------------------------############
  # Get total counts of each evaluation
  counts_t<- as.data.frame(table(vc$ploidy, vc$eval_sample))
  colnames(counts_t)<- c("ploidy", "eval", "counts")
  counts_t<- counts_t[counts_t$counts !=0, ]
  # # Get per sample
  # sample_t<-  as.data.frame.table(table(vc$ploidy, vc$eval_sample, vc$sample))
  # colnames(sample_t)<- c("ploidy", "eval", "sample","counts")
  # sample_t<- sample_t[sample_t$counts !=0, ]
  
  ####### ---------------------------------- 1. GET General STATS -----------------------------############
  
  # Fill table
  score <-  get_metrics(counts_t)
  colnames(score)<-  c("set", "precision",    "recall",  "F1_score")  
  score$type<- "vc"
  score$eval_level<- paste0("Filter", v)
  # Export performance table:
  write.table(score, paste0(dirf,  "Performance_sim_errors_vc_filtered_v", v, ".tbl"), sep=",", quote=F, row.names=F, col.names=T)
}


print("All done!")
# ####### ---------------------------------- 2. GET sample STATS -----------------------------############
# get_metrics<- function(d){
#   # Initialize recount table
#   score <-  data.frame(matrix(nrow=0, ncol=5))
#   colnames(score)<-  c("ploidy","sample", "precision", "recall", "F1_score")
#   
#   for (ploidy in unique(d$ploidy)){
#     FP <- 0
#     TP <- 0
#     FN <- 0
#     for (sample in paste0("sample_",1:100)){
#         rows <- d[d$ploidy==ploidy &d$sample==sample,]
#         
#         for (i in 1:nrow(rows)){ # Get estimates
#           row<- rows[i,]
#           if(row$eval == "FN_uncalled") {
#             FN <- FN + as.numeric(row$counts)
#           }else if ( row$eval =="TP"){
#             TP <- TP + as.numeric(row$counts)
#           }else{
#             FP <-  FP+ as.numeric(row$counts)
#           }
#         }
#         # accuracy vs mistakes 
#         precision <-   TP/(TP+FP) 
#         # how good at finding
#         recall <-  TP/(TP+FN)  
#         
#         score<- rbind(score,  c(ploidy,sample,  precision, recall,   2 * precision * recall / (precision +recall)))
#         #print(paste("done with", sample)
#     }
#   }
#   return(score)
# }
# 
# # Fill table
# score <-  get_metrics(sample_t)
# colnames(score)<-  c("set", "precision",    "recall",  "F1_score")  
# score$type<- "vc"
# score$eval_level<- "sample"
# # Export performance table:
# write.table(score, paste0(dir,  "Performance_sim_errors_vc_sample.tbl"), sep=",", quote=F, row.names=F, col.names=T)
# 


