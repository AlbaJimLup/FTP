####### GET DATA  ######
dir="~/Desktop/ribosomal_RNAs/Alba/01_VCkmers/"
#
s30 <- read.table(paste0(dir, "x30stats_evaluationGATKvsMSA.tbl"), 
                  header = TRUE, row.names = NULL, sep = ",")
s100 <- read.table(paste0(dir, "x100stats_evaluationGATKvsMSA.tbl"), 
                   header = TRUE, row.names = NULL, sep = ",")
s1000 <- read.table(paste0(dir, "x1000stats_evaluationGATKvsMSA.tbl"), 
                    header = TRUE, row.names = NULL, sep = ",")

dir = "~/Desktop/ribosomal_RNAs/Alba/07_Performance/"
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)

###### PREP DATA   ########

s30$coverage <- rep("30x", nrow(s30))
s100$coverage <- rep("100x", nrow(s100))
s1000$coverage <- rep("1000x", nrow(s1000))


counts_ev <-  rbind(s30, s100, s1000)

counts_ev$evaluation[counts_ev$evaluation =="same variant"] <-  "TP"
counts_ev$evaluation[counts_ev$evaluation =="missing variant"] <- "FN"
counts_ev$evaluation[counts_ev$evaluation =="missing allele"] <- "FN_allele" 
counts_ev$evaluation[counts_ev$evaluation =="different alleles"] <- "FP_allele"

counts_ev$evaluation <- factor(counts_ev$evaluation, 
                               levels = c( "FP_allele", "FP", "FN_allele","FN", "TP"))

## Select only callable positions
s30 <- counts_ev[counts_ev$coverage =="30x", ]
s30 <-  s30[s30$class=="callable",]
rownames(s30) <- s30$evaluation
s30 <- t(s30)
s30 <-  s30["counts", ]
s30 <- as.data.frame(t(s30))

s100  <- counts_ev[counts_ev$coverage =="100x", ]
s100 <-  s100[s100$class=="callable",]
rownames(s100) <- s100$evaluation
s100 <- t(s100)
s100 <-  s100["counts", ]
s100 <- as.data.frame(t(s100))


s1000 <- counts_ev[counts_ev$coverage =="1000x", ]
s1000 <- s1000[s1000$class=="callable",]
rownames(s1000) <- s1000$evaluation
s1000 <- t(s1000)
s1000 <-  s1000["counts", ]
s1000 <- as.data.frame(t(s1000))

########### GET SCORING ##########
getmetrics<- function(row){
  
  FP <-ifelse("FP" %in% names(row), as.numeric(row$FP), 0) 
  TP <- as.numeric(row$TP)
  FN <- as.numeric(row$FN) +  ifelse("FN_allele" %in% names(row), as.numeric(row$FN_allele), 0) 
  
  # accuracy vs mistakes 
  precision <-   TP/(TP+FP) 
  # how good at finding
  recall <-  TP/(TP+FN)  
  
  return(list(ploidy = row$ploidy,
              precision =  precision, 
              recall= recall,   
              #harmonic mean of precision and recall.
              F1_score= 2 * precision * recall / (precision +recall)))
}

kmers <-  rbind(s30, s100, s1000)
kmers$ploidy <-  c("30x", "100x", "1000x")

tkmers <-  data.frame(matrix(nrow=3, ncol=4))
colnames(tkmers)<-  c("depth", "precision", "recall", "F1_score")

#Fill tables
for (i in 1:nrow(tkmers)){ 
  tkmers[i,]<- getmetrics(kmers[i,])
}

colnames(tkmers)<-  c("set", "precision", "recall",  "F1_score")  
tkmers$type<- "kmers"
tkmers$eval_level<-"population"

# Export performance table:
write.table(tkmers, paste0(dir,  "Performance_perfect_reads.tbl"), sep = ",", quote = FALSE, row.names = FALSE)


