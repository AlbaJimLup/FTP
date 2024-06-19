## libraries 
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)

## ------------------------------------ SET PARAMETERS --------------------------------- ############
#### CHOOSE Join genotyping SET ######
args <- commandArgs(trailingOnly = TRUE)
#GATK version either "4.5.0.0" or "4.3.0.0"
set=args[1] 
# STEP of the pipeline this is for gVCF files so either "joingen" or"filtering"
step=args[2]

if (step=="filtering"){
  version=args[3]
  ploydies<- unlist(strsplit(args[4], " "))
  
}else{
  version=""
  ploydies<- unlist(strsplit(args[3], " "))
}

# ploydies<- c(2,5,7,10,20,50,100)

print("#############################################################################")
print(paste("Merging postjoingenotyping evaluation for", set, step, version, as.character(ploydies)))

dir <-  "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/"

#dir <- "~/Desktop/ribosomal_RNAs/Alba/03_VCsimulated/v2/"
dirf<- paste0(dir, "VC_out/Sample_level_eval/", set, "/joingen/")


all<- data.frame()
report <- data.frame() # keep counts of found variants per region

if (step=="filtering"){
  for (ploidy in ploydies){
    print( paste0(dirf,  "VCeval_postjoin_p", ploidy, "_", step, "_", version, ".tbl"))
    d<- read.table( paste0(dirf,  "VCeval_postjoin_p", ploidy, "_", step,"_", version, ".tbl"), sep="\t")
    # One dataframe all ploidies 
    all<-  rbind(all, d)
  }
  all <- all[, !colnames(all) %in% c("REF", "ALT","econded")]
  write.table(all, paste0(dirf,"VCeval_postjoin_", step,"_", version, ".tbl"), sep="\t")
    
}else{
  for (ploidy in ploydies){
    print( paste0(dirf,  "VCeval_postjoin_p", ploidy, "_", step, ".tbl"))
    d<- read.table( paste0(dirf,  "VCeval_postjoin_p", ploidy, "_", step,".tbl"), sep="\t")
    # One dataframe all ploidies 
    all<-  rbind(all, d)
  }
  all <- all[, !colnames(all) %in% c("REF", "ALT","econded")]
  write.table(all, paste0(dirf,"VCeval_postjoin_", step,".tbl"), sep="\t")
}

print("All done!")   
print("#############################################################################")
 

###########------------------------- Plot variants found per region  ------------------------------- #############
# colnames(report)<- c("ploidy", "region", "counts")
# 
# report$region<-  factor(report$region, levels =c("5_ETS", "18S", "ITS1", "5.8S", "ITS2",  "28S", "3_ETS") )
# report$ploidy<-  factor(report$ploidy, levels =c("2", "5", "7", "10") )
# report$counts<- as.numeric(report$counts)
# 
# 
# palette<- c( "#2f4b5b","#6ad68e", "#9dbccf", "#90c230", "#4d8195",  "#71d064",  "#375794")
# 
# 
# png(file=paste0(dirI,"Analysis/Called_variants_per_region.png"), res=500, width=2500, height=2500);
# ggplot(report, aes(x=ploidy, y=counts, fill=region))+
#   geom_bar(stat="identity")+
#   geom_text_repel(aes(label=region),force =0.001, position = position_stack(vjust = 0.5),
#                   fontface="bold", size=3.5, color="white", max.overlaps =50)+
#   scale_fill_manual(values=palette)+
#   labs(x="Ploidy", y="Number of variants called", fill="Region")+
#   theme_minimal()
# dev.off();
# 
# ggplot(report, aes(x=region, y=counts, fill=region))+
#   geom_bar(stat="identity")+
#   facet_grid(ploidy~.)+
#   scale_fill_manual(values=palette)
# 
# ####### INCLUDE TV FOR COMPARISON #######
# 
# tv<- as.data.frame(table(true_variants$region))
# names(tv)<- c("region", "counts")
# tv$set <- "true_variants"
# tv$region<- c("18S",   "28S", "3_ETS",  "5.8S", "5_ETS",  "ITS1",  "ITS2" )
# tv$region<-  factor(tv$region, levels =c("5_ETS", "18S", "ITS1", "5.8S", "ITS2",  "28S", "3_ETS") )
# 
# 
# colnames(report)<- c("set", "region", "counts")
# report$set<-  paste0("ploidy_", c(rep("2",7), rep("5",7), rep("7",7), rep("10",7)))
# report<- rbind(report, tv)
# report$set<-  factor(report$set, levels =c("true_variants","ploidy_2", "ploidy_5", "ploidy_7", "ploidy_10") )
# 
# 
# png(file=paste0(dirI,"Analysis/Called_tv_variants_per_region.png"), res=500, width=2500, height=2500);
# ggplot(report, aes(x=set, y=counts, fill=region))+
#   geom_bar(stat="identity")+
#   geom_text_repel(aes(label=region),force =0.001, position = position_stack(vjust = 0.5),
#                   fontface="bold", size=3.5, color="white", max.overlaps =50)+
#   scale_fill_manual(values=palette)+
#   labs(x="", y="Number of variants called", fill="Region")+
#   theme_minimal()
# dev.off();
