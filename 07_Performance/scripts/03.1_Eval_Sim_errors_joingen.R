######### Directories ########
dirt ="~/Desktop/ribosomal_RNAs/Alba/07_Performance/"
dirt<-"/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/07_Performance/"
# dir ="~/Desktop/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/Sample_level_eval/4.3.0.0/joingen/"
dir <-  "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/Sample_level_eval/4.3.0.0/joingen/"
####### libraries ######
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)
####### Palette ######
palette <- c( "precision"="#3aa8ec",      "recall"="#fee08b",       "F1_score"="#66c2a5")
############----------------------------- GET Frequencies --------------------------- #############
## 1. POP level Load variants data
d<-  read.table(  paste0(dir,"VCeval_postjoin_joingen.tbl"), sep="\t") # after all analysis Rscript files
d<- d[d$class=="callable", ]
####### ---------------------------------- COUNT EVALS -----------------------------############
# Get total counts of each evaluation
counts_t<- as.data.frame(table(d$ploidy, d$eval_sample))
colnames(counts_t)<- c("ploidy", "eval", "counts")
counts_t<- counts_t[counts_t$counts !=0, ]
# Get per sample
# sample_t<-  as.data.frame.table(table(d$ploidy, d$eval_sample, d$sample))
# colnames(sample_t)<- c("ploidy", "eval", "sample","counts")
# sample_t<- sample_t[sample_t$counts !=0, ]

####### -------------------------------1. GET general STATS -----------------------------############
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

# Fill table
score <-  get_metrics(counts_t)
colnames(score)<-  c("set", "precision","recall",  "F1_score")  
score$type<- "join-genotyping"
score$eval_level<- "NoFilter"
# Export performance table:
write.table(score, paste0(dirt,  "Performance_sim_errors_postjoin.tbl"), sep=",", quote=F, row.names=F, col.names=T)


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
#       rows <- d[d$ploidy==ploidy &d$sample==sample,]
#       
#       for (i in 1:nrow(rows)){ # Get estimates
#         row<- rows[i,]
#         if(row$eval == "FN_uncalled")      FN <- FN + as.numeric(row$counts)
#         else if ( row$eval =="TP")         TP <- TP + as.numeric(row$counts)
#         else    FP <-  FP+ as.numeric(row$counts)
#       }
#       # accuracy vs mistakes 
#       precision <-   TP/(TP+FP) 
#       # how good at finding
#       recall <-  TP/(TP+FN)  
#       
#       score<- rbind(score,  c(ploidy,sample,  precision, recall,   2 * precision * recall / (precision +recall)))
#     }
#   }
#   return(score)
# }
# 
# # Fill table
# score <-  get_metrics(sample_t)
# colnames(score)<-  c("set", "precision",    "recall",  "F1_score")  
# score$type<- "post-join-gen"
# score$eval_level<- "sample"
# # Export performance table:
# write.table(score, paste0(dirt,  "Performance_sim_errors_d_sample.tbl"), sep=",", quote=F, row.names=F, col.names=T)
# 
# 
