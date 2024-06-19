## libraries 
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)
#### CHOOSE VC SET ######
args <- commandArgs(trailingOnly = TRUE)
#GATK version either "4.5.0.0" or "4.3.0.0"
ploidy=args[1]
version=args[2]
ploidy=7
version="v2"
step="joingen"

samples<- c("SRR1997411",  "SRR3189741",  "SRR3189742",  "SRR3189743")

# set="4.3.0.0"
# ploydies<- c(2,5,7,10,20,50,100)
#########################-------------------------GET DIRECTORIES-------------------------############
dir <-  paste0("/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/08_validation/VC_out/" )

dir <- "~/Desktop/ribosomal_RNAs/Alba/08_validation/VC_out/"
dirf<- paste0(dir,  "ploidy_", ploidy)
dirout<- paste0(dir)

print("#############################################################################")
print(paste("Performing join genotyping filtering for", set," and ploidy"  ,ploidy))
print("#############################################################################")
###############----------------- GET file ---------------------###############
get_vc <-  function(p=ploidy) { # each ploidy 
  
  ts <- data.frame() # table of all regions to be returned for eahc ploidy
  variants <- data.frame(matrix(ncol=3, nrow=0)) #table with reported number of variants for each region and ploidy
  
  t_raw <- read.table(gzfile(paste0(dirf, "/WGS.", p ,".vcf.gz" )))  
  
  colnames(t_raw)<- c("CHROM",	"position",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT", samples)
  
  ## Create INFO columns
  info_cols<- c("AC","AF","AN","BaseQRankSum","DP","FS","MLEAC","MLEAF","MQ","MQRankSum","QD","ReadPosRankSum","SOR")
  
  for (cols in info_cols  ){
    t_raw[[cols]] <- gsub(paste0(".*", cols, "=([^;]+).*"), "\\1", t_raw$INFO)
  }
  # Initialize ploidy region and evaluation result
  t_raw$ploidy<- p
  t_raw$eval_sample<- ""
  
  ## Samples to long format
  t <- t_raw %>% pivot_longer(cols =c("SRR1997411",  "SRR3189741",  "SRR3189742",  "SRR3189743"), 
                              names_to = "sample", 
                              values_to = "GT_AD_DP_GQ_PL")
  ## Create GT:AD:DP:GQ:PL columns
  format_cols <- c("GT", "AD", "DP", "GQ", "PL")
  new_cols <- matrix(nrow = length(t$GT_AD_DP_GQ_PL), ncol = length(format_cols))
  colnames(new_cols) <- format_cols
  
  for (i in seq_along(t$GT_AD_DP_GQ_PL)) {
    row <- unlist(strsplit(t$GT_AD_DP_GQ_PL[i], ":"))
    new_cols[i, ] <- row[1:length(format_cols)]
  }
  ##  Add region to dataframe 
  t<- cbind(t, new_cols) # Add to t
  # Remove void or compacted info columns: 
  t<- t[, ! colnames(t) %in% c("CHROM", "ID",  "INFO", "FORMAT", "GT_AD_DP_GQ_PL")]
  
  print(paste0("Fetched with ploidy ", p, "   with ", length(t$position), " called variants"))

  return(as.data.frame(t))
}
################ ------------------------ Apply filter and export file-----------------------------#################33
print(paste("Fetching join genotyping vcf files for ploidy", ploidy, "and applying filter"))
## Get join gen vcf files
d<- as.data.frame(get_vc(ploidy))
print(paste0("Applying filters"))
# Keep ony those that passed the filter:
################################################ 
if (ploidy==2){
  GQ<- 20 # the GQ for ploidy 2 s 99 which removes all variants
}else{
  GQ<- quantile(as.numeric(d$GQ), probs = 0.25, na.rm = TRUE)
  GQ<-0
}
# How I got here:  https://docs.google.com/presentation/d/1yXJ4ixGVS7xoPpFzqItunfJVn5HqrQP7x5y5n1j9z0E/edit#slide=id.g2e03fd93fb6_0_199
if (version=="v1"){
  filtered <-  d[as.numeric(d$GQ)>GQ  , ]

}else if (version=="v2"){
    filtered <-  d[as.numeric(d$GQ)>GQ &
                     as.numeric(d$QD)>1 &
                     as.numeric(d$SOR)<3 &
                     as.numeric(d$QUAL)>30  , ]
}else if (version=="v3"){
  filtered <-  d[as.numeric(d$GQ)>(GQ)  &
                   as.numeric(d$QD)>1 &
                   as.numeric(d$SOR)<3 &
                   as.numeric(d$QUAL)>30 &
                    as.numeric(d$ReadPosRankSum)<6 , ]
}
################################################
print(paste("Done, now writting file into directory", dirf))
# write.table(d, paste0(dirout,  "/joingen_Rfiltering_p", ploidy, ".tbl"), sep="\t")
write.table(filtered, paste0(dirout,  "/joingen_Rfiltering_", version,"_p", ploidy, ".tbl"), sep="\t")

print(paste("Done with ", set, step," for ploidy ",  ploidy))
print(dim(d))
print(dim(filtered))
print("#######################################################################################")



ggplot(filtered, aes(x=1, y=as.numeric(GQ)))+
  geom_jitter()+
  geom_boxplot()
