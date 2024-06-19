#  Libraries 
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)
library(stringr)
###### Functions######
get_metrics<- function(d){
  # Initialize recount table
  score <-  data.frame(matrix(nrow=0, ncol=4))
  
  for (filter in unique(d$filter)){
    FP <- 0
    TP <- 0
    FN <- 0
    rows <- d[d$filter==filter,]
    
    for (i in 1:nrow(rows)){ # Get estimates
      row<- rows[i,]
      if(row$eval_sample == "FN_uncalled")      FN <- FN + as.numeric(row$count)
      else if ( row$eval_sample =="TP")         TP <- TP + as.numeric(row$count)
      else    FP <-  FP+ as.numeric(row$count)
    }
    # accuracy vs mistakes 
    precision <-   TP/(TP+FP) 
    # how good at finding
    recall <-  TP/(TP+FN)  
    
    score<- rbind(score, c(filter,  precision, recall,2 * precision * recall / (precision +recall)))
  } 
  colnames(score)<-  c("filter", "precision", "recall", "F1_score")
  return(score)
}
##########
step <- "filtering"

dirf<-"/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/07_Performance/"
dir <-  "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/Sample_level_eval/4.3.0.0/VC/"

# dir <- paste0("~/Desktop/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/Sample_level_eval/4.3.0.0/VC/")
# dirf<- "~/Desktop/ribosomal_RNAs/Alba/07_Performance/"

## --------------------------------------  LOAD DATA ---------------------------------############

# ~/Desktop/ribosomal_RNAs/Alba/03_VCsimulated/v2/scripts/Rscripts/01.3.eval_prejoin_plotting_sample.R
# write.table(dfs_c_counts_ploidy, paste0(dir, "VC_filtered_sample.tbl"),sep=",", quote=F, row.names=F, col.names=T)
vc<-  read.table( paste0(dir,"VC_filtered_sample.tbl"), sep=",", header =  T)# after all analysis Rscript files

####### ---------------------------------- 1. GET General STATS -----------------------------############
# Fill table
score <-  get_metrics(vc)
colnames(score)<-  c("set", "precision", "recall",  "F1_score") 

score$ploidy <- word(score$set, 1)
score$sample <- factor(word(score$set, 2))
score$qual <- factor(word(score$set, 3))

score$type<- "vc"
# Export performance table:
write.table(score, paste0(dirf,  "Performance_sim_errors_vc_filtered_sample.tbl"), sep=",", quote=F, row.names=F, col.names=T)

## --------------------------------------  LOAD DATA ---------------------------------############

# ~/Desktop/ribosomal_RNAs/Alba/03_VCsimulated/v2/scripts/Rscripts/01.3.eval_prejoin_plotting_sample.R
# write.table(dfs_c_counts_ploidy, paste0(dir, "VC_filtered_sample.tbl"),sep=",", quote=F, row.names=F, col.names=T)
vc<-  read.table( paste0(dir,"VC_filtered_2sample.tbl"), sep=",", header =  T)# after all analysis Rscript files

####### ---------------------------------- 1. GET General STATS -----------------------------############
# Fill table
score <-  get_metrics(vc)
colnames(score)<-  c("set", "precision", "recall",  "F1_score") 

# score$ploidy <- word(score$set, 1)
# score$sample <- factor(word(score$set, 2))
# score$qual <- factor(word(score$set, 3))

score$type<- "vc"
score$eval_level<- "Filter"
# Export performance table:
write.table(score, paste0(dirf,  "Performance_sim_errors_vc_filtered_2sample.tbl"), sep=",", quote=F, row.names=F, col.names=T)


print("All done!")
# ####### ---------------------------------- 2. GET sample STATS -----------------------------############
# get_metrics<- function(d){
#   # Initialize recount table
#   score <-  data.frame(matrix(nrow=0, ncol=5))
#   colnames(score)<-  c("filter","sample", "precision", "recall", "F1_score")
#   
#   for (filter in unique(d$filter)){
#     FP <- 0
#     TP <- 0
#     FN <- 0
#     for (sample in paste0("sample_",1:100)){
#         rows <- d[d$filter==filter &d$sample==sample,]
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
#         score<- rbind(score,  c(filter,sample,  precision, recall,   2 * precision * recall / (precision +recall)))
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


