## ------------------------------------ SET PARAMETERS --------------------------------- ############
#### CHOOSE VC SET ######
args <- commandArgs(trailingOnly = TRUE)
#GATK version either "4.5.0.0" or "4.3.0.0"
set=args[1] 
# STEP of the pipeline this is for gVCF files so either "VC" or"filtering"
step=args[2]

if (step=="filtering"){
  version=args[3]
  ploydies<- unlist(strsplit(args[4: length(args)], " "))
  
}else{
  ploydies<- unlist(strsplit(args[3: length(args)], " "))
}


print("#############################################################################")
print(paste("Merging post-joingenotyping evaluation for", set, step, as.character(ploydies)))

## libraries 
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)

dir <-  "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/"
#dir <- "~/Desktop/ribosomal_RNAs/Alba/03_VCsimulated/v2/"
dirf<- paste0(dir, "VC_out/Sample_level_eval/", set, "/VC/")

#############---------------------------- GET TABLES  -------------------------##################
all<- data.frame()

if(step=="filtering"){
  for (ploidy in ploydies){
      print( paste0(dirf,  "VCeval_prejoin_p", ploidy, "_", step,"_", version,".tbl"))
      d<- read.table(paste0(dirf, "VCeval_prejoin_p", ploidy, "_", step, "_", version,".tbl"), sep="\t", header = T)
      # Save to have a single output file with all ploidies 
      all<- rbind(all, d)
  }
  write.table(all, paste0(dirf, "All_VCeval_prejoin_", step, "_", version, ".tbl"), sep="\t")
}else{
  for (ploidy in ploydies){
      print( paste0(dirf,  "VCeval_prejoin_p", ploidy, "_", step, ".tbl"))
      d<- read.table(paste0(dirf, "VCeval_prejoin_p", ploidy, "_", step, ".tbl"), sep="\t", header = T)
      # Save to have a single output file with all ploidies 
      all<- rbind(all, d)
  }
  write.table(all, paste0(dirf, "All_VCeval_prejoin_", step, ".tbl"), sep="\t")
}

print("Done!")
print("#############################################################################")