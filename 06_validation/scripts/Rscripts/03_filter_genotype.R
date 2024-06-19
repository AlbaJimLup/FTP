## libraries 
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)
## ------------------------------------ SET PARAMETERS --------------------------------- ############
##CHOOSE VC SET 
dir <-  "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/08_validation/"

args <- commandArgs(trailingOnly = TRUE)
ploidy=args[1]
filter=args[2]
samples=strsplit(args[3], " ")

print(paste("Doing VC evaluation for", sample, ploidy,"and filtering", filter, "samples"))
print("#############################################################################")

# dir <- "~/Desktop/ribosomal_RNAs/Alba/08_validation/"
# ploidy="50"
# filter="4"
samples=unlist(strsplit("SRR1997411 SRR3189741 SRR3189742 SRR3189743", " "))
# Directory where files are
dirf<- paste0(dir, "VC_out/")
dirI<- paste0(dir, "IMAGES/")
## --------------------------------------  LOAD DATA --------------------------------- ############
# Get bed callable blacklisted regions  
bed <- read.table("~/Desktop/ribosomal_RNAs/Alba/data/47S_pre-rRNA.blacklisted_repetitive_seq.bed")
## --------------------------------------  FUNCTIONS   --------------------------------- ############
# Get and prep file
get_vc_filtering <-  function(sample) { # each ploidy we ran so far
  
    t <- read.table(gzfile(paste0(dirf, sample, "/", sample, ".ploidy_", ploidy,"_filtering.g.vcf.gz")))  
    
    t$sample<-sample
    t$ploidy<- ploidy
    
    # Remove "<NON_REF>" from all ALT rows
    t$V5 <- gsub(",<NON_REF>", "", t$V5)   
    # Select only called positions
    t <- t[t$V5 != "<NON_REF>", ] # non-called but possible non-reference alleles
    
    colnames(t)<- c("CHROM",	"position",	  "ID", 	"REF",  	"ALT",
                   "QUAL", "FILTER",	"INFO",	"FORMAT",	"data", "sample", "ploidy")
  
    #### Splice data #### 
    # Create INFO columns
    info_cols<- c( "BaseQRankSum","DP", "ExcessHet","MLEAC","MLEAF","MQRankSum","RAW_MQandDP","ReadPosRankSum")
      
    for (cols in info_cols  ){
        t[[cols]] <- gsub(paste0(".*", cols, "=([^;]+).*"), "\\1", t$INFO)
    }
    ## Create GT:AD:DP:GQ:PL columns
    format_cols <- c("GT", "AD", "DP", "GQ", "PL", "SB")
    new_cols <- matrix(nrow = length(t$FORMAT), ncol = length(format_cols))
    colnames(new_cols) <- format_cols
    
    for (i in seq_along(t$FORMAT)) {
      row <- unlist(strsplit(t$data[i], ":"))
      new_cols[i, ] <- row[1:length(format_cols)]
    }
    ##  Add region to dataframe 
    t<- cbind(t, new_cols) # Add to t
    # Remove void or compacted info columns: 
    print("We decoded additional annotations information info different columns")
    
  return(as.data.frame(t[, ! colnames(t) %in% c("CHROM", "INFO","FORMAT")]))
}
# Classify callable blacklisted  
get_black<-  function(df){
  
  for (pos in unique(df$position)){
    
    eval_black <-  unique(ifelse(pos > bed$V2 & pos <= bed$V3, "YES", "NO"))
    
    if(length(unique(eval_black)) == 1){
      df$class[df$position == pos] <- "callable"
    }
    else{
      df$class[df$position == pos] <- "blacklisted"
    }
  }
  print("We classified callable and blacklisted regions")
  return(df)
}
## Classify rDNA region 
add_regions_chrR <- function(d){
  positions<-as.numeric(unique(d$position))
  # <= only in IGS at start but all other should be > only but then variants are 0 bed based
  #so variant at the end of region be annotated to be in wrong region so we should use >=
  for (POS in positions){
    d$region[d$position==POS] <- ifelse(POS <= 9338, "IGS",
                                        ifelse(POS <= 12995, "5_ETS",
                                               ifelse(POS <= 14864, "18S",
                                                      ifelse(POS <= 15934, "ITS1",
                                                             ifelse(POS <= 16091, "5.8S",
                                                                    ifelse(POS <= 17258, "ITS2",
                                                                           ifelse(POS <= 22309, "28S",
                                                                                  ifelse(POS <= 22670, "3_ETS", "IGS"))))))))
  }
  print("We annotated chrR regions")
  return(d)
}

## --------------------------------  RUN FOR THE SAMPELS   --------------------------------- ############

ds<- data.frame()

for (sample in samples){
  d<- get_vc_filtering(sample)
  d<- get_black(d)
  d<- add_regions_chrR(d)
  
  ds<- rbind(ds, d)
}

report<- data.frame(matrix(ncol=4, nrow=7))
names(report)<- c("filter", "step","#var", "plot")
report$filter<- c("Variant calling", "QUAL1000","QUAL1000","QUAL1000",  "QUAL30",  "QUAL30","QUAL30")
report$step<- factor(c("Variant calling", "QUAL1000","2samples","4samples","QUAL30","2samples", "4samples"), 
                     levels=c("Variant calling","QUAL1000","QUAL30","2samples", "4samples"))

t<- as.data.frame(table(ds$position))
report$`#var`[1]<- dim(t)[1]
report$plot[1]<- dim(t)[1]

# Apply filter, remove those that don't pass the filter, 
# keep QUAL>1000 and positions present in at least 5 samples

qual1000<- ds[ ds$FILTER == "PASS" ,  ]
t<- as.data.frame(table(qual1000$position))
report$`#var`[2]<- dim(t)[1]
# check position frequency
report$`#var`[3]<- dim(t[t$Freq>1,])[1]
report$`#var`[4]<- dim(t[t$Freq>3,])[1]
report$plot[2]<-  dim(t)[1]-dim(t[t$Freq>1,])[1]
report$plot[3]<-  dim(t[t$Freq>1,])[1]-dim(t[t$Freq>3,])[1]
report$plot[4]<- dim(t[t$Freq>3,])[1]

qual30<-  ds[as.numeric(ds$QUAL)>30, ]
t<- as.data.frame(table(qual30$position))
report$`#var`[5]<-  dim(t)[1]
report$`#var`[6]<- dim(t[t$Freq>1,])[1]
report$`#var`[7]<- dim(t[t$Freq>3,])[1]

report$plot[5]<-  dim(t)[1]-dim(t[t$Freq>1,])[1]
report$plot[6]<-  dim(t[t$Freq>1,])[1]-dim(t[t$Freq>3,])[1]
report$plot[7]<- dim(t[t$Freq>3,])[1]
############################# Plot ###############################33

palette <-  c("Variant calling"= "#b1cae4","QUAL1000"="#9fb5cb" ,"QUAL30"="#82b5bf",
             "4samples"="#c0eca3", "2samples"="#8dbfaa") 

png(file=paste0(dirI, "Num_vairants_filter.png"), res=500, width=2100, height=2100);
ggplot(report, aes(y=plot, x=reorder(filter,-`#var`), fill=step))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=palette)+
  geom_text_repel(aes(label=`#var`), position="stack", vjust=2, force=0.005, fontface = "bold" )+
  # guides(fill="none")+
  labs(x="Filter", y="Number of reported variants", fill="Filter step")+
  theme_minimal()
dev.off();

