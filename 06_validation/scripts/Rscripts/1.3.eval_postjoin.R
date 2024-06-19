library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)

##  Palette  
palette <-  c("TP" ="#75bb6f","TP_pos_FN_allele"="#aeb852", "TP_pos_other"="#706f4d", "FN_uncalled"="#ffc60e", "FN_ref_homozygous"= "#4c89ee","FN_alt_homozygous"="#71a5be",
              "FP"= "#bd0752", "TP_pos_FP"= "#ff4d7a", "TP_pos_FP_allele"="#ed8c96","TP_allele_FP_allele"="#ed8c96",  "err"= "#a6a6a6") 
### ------------------------------------ SET PARAMETERS --------------------------------- ############
#CHOOSE VC SET 
args <- commandArgs(trailingOnly = TRUE)
#GATK version either "4.5.0.0" or "4.3.0.0"
set=args[1] 
# STEP of the pipeline this is for gVCF files so either "VC" or"filtering"
step=args[2]
ploidy=args[3]
version=args[4]

set="4.3.0.0"
step="joingen"
step="filtering"
ploidy=7
version="v2"

print(paste("Doing joingen evaluation for", set, step, ploidy, version))
print("#############################################################################")

## libraries 
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)

dir <-  "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/08_validation/"

dir <- "~/Desktop/ribosomal_RNAs/Alba/08_validation/"
dirf<- paste0(dir, "VC_out/")
dirI<- paste0(dir, "IMAGES/")
## --------------------------------------  LOAD DATA --------------------------------- ############
### Get true_variants ######
true_variants <- read.table(paste0(dir, "Rscripts/aligned_10.vcf"))
names(true_variants)<- c("chr", "position", "id", "ref", "alt")
true_variants$position<- true_variants$position+9338
not_FN<- c()
### Get bed callable blacklisted regions  ######
bed <- read.table(paste0("~/Desktop/ribosomal_RNAs/Alba/data/47S_pre-rRNA.blacklisted_repetitive_seq.bed"))
#############---------------------------- FUNCTIONS  -------------------------##################
## Get pre join gen g.vdf file output either HaplotypeCaller or VariantFilter ##
get_vc <-  function(p){ # each ploidy we ran so far
  if(step=="joingen"){
        t<-read.table(gzfile(paste0(dirf, "/joingen_Rfiltering_p", p ,".tbl" )))
  }else{
        t<- read.table(gzfile(paste0(dirf, "/joingen_Rfiltering_", version ,"_p", p ,".tbl" )))   
      }
                                   
  t<- as.data.frame(t)
  t$ploidy<- p
  t$eval_sample<- ""

  return(t)
}
#### Get GT column #### 
get_GT<- function(df){
  GTs <- c()
  for (i in 1:nrow(df)){
    row<- df[i,]
    encoded <- df$GT
    a <-  c(row$REF,  unlist(strsplit(row$ALT, ",")))
    
    if (grepl("/", encoded)) {# using "/" delimiter
      alleles <- unlist(strsplit(encoded, "/"))
      
    } else if (grepl("\\|", encoded)) {# using "|" delimiter
      alleles <- unlist(strsplit(encoded, "\\|"))
    }
        GT<- c()# Decode
    for (allele in alleles)  {
      GT<-  c(GT, a[as.numeric(allele)+1])
    }
    
    GTs<-c(GTs, paste(GT, collapse = "/"))# Collapse
  }
  print("Done fetching GT column")
  return(GTs)
}
#### VC evaluation sample level#### 
classify_variant_sample <- function(position, sample, df, true_variants) {
  # print(position)
  position<- as.numeric(position)
  # if the variant is not one of the true callable variants for this position
  if (!position %in% true_variants$position) {
    # called<- unique(unlist(strsplit(as.character(df[df$position == position & df$sample == sample, "GT"]), split = "[/|]")))
    # if ((position+1) %in% true_variants$position[true_variants$sample == sample]) { # shifted e.g.
    #   
    # }else{
      return("FP") # FP sample levelNon existent variant
    # }
  } else {
    if (df[df$position == position & df$sample == sample, "GT"] == ".") {
      return("FN_uncalled")
    } else {
      alleles <- unique(unlist(strsplit(as.character(df[df$position == position & df$sample == sample, "GT"]), split = "[/|]")))
      # if only one allele
      ref_allele <- unique(df[df$position == position & df$sample == sample, "REF"])
      alt_allele <- unique(unlist(strsplit(as.character(df[df$position == position & df$sample == sample, "ALT"]), split = ",")))
      
      if (length(alleles) == 1) {
        ifelse(alleles == ref_allele, return("FN_ref_homozygous"), return("FN_alt_homozygous"))

      } else {# if heterozygous check that this sample is supposed to have the variant or not
        true_alleles <- c(unique(true_variants[true_variants$position == position  , "ref"]),
                          sort(strsplit(as.character(true_variants[true_variants$position == position, "alt"]), ",")[[1]]))

        # gatk_alleles <- alleles
        gatk_alleles <- c(unique(as.character(df[df$position == position & df$sample == sample, "REF"])),
                          sort(strsplit(as.character(df[df$position == position & df$sample == sample, "ALT"]), ",")[[1]]))
        
        condition1 <- length(true_alleles) == length(gatk_alleles) # same number of alleles
        condition2 <- sum(gatk_alleles %in% true_alleles) == length(gatk_alleles) # all gatk alleles TP
         # all TP alt alleles are identified +  some FP alt alleles
        condition3 <- sum(true_alleles[true_alleles != ref_allele] %in%  gatk_alleles[gatk_alleles != ref_allele]) == length(true_alleles[true_alleles != ref_allele])
        # FP allele -> alt allele not in alt gatk alleles
        condition4 <- sum(true_alleles[true_alleles != ref_allele] %in%   gatk_alleles[gatk_alleles != ref_allele]) == 0
        
        if (condition1 & condition2) { # same number of alleles & same alleles
          return("TP")
        } else {
          if (!condition1 & condition2) { # Multi-allelic variants we are missing alleles  
                return("TP_pos_FN_allele") # all gatk alleles are true alleles but some true alleles were not identified 
          } else if (!condition1 & condition3) {
               return("TP_allele_FP_allele") # Found variant + FP some false allele were identified
          } else if (condition4) { # Diverging alleles even if correct position
                if (true_alleles[1] == gatk_alleles[1]) {
                  return("TP_pos_FP")
                }else {
                  if ("*" %in% alt_allele) { # Deletions may be wrongly called
                    ref <- substring(true_alleles[1], 1, nchar(true_alleles[1]) - 1)
                    # # 12469 GGCCCC  *,G,TGCCCC   0     sample_39    2  
                    # # 12469:GGCCCCG>G   sample_39    indel   blacklisted      5'ETS 12469 GGCCCCG  G
                    ifelse(ref == gatk_alleles[1], return("TP_allele_FP_allele"), return("TP_pos_FP"))
                  }else {
                    return("TP_pos_FP")
                  }
                }
          }else {
              return("TP_pos_FP") ## Same ref different alleles
            }
          }
        }
      }
    } 
}
fill_eval_sample <- function(df, true_variants){
  
  FNs <- data.frame(matrix(ncol=ncol(df), nrow=0)) # initialize missed variants
  samples <-c("SRR1997411", "SRR3189741", "SRR3189742", "SRR3189743") # num simulated samples
  
  for (sample in samples){ # each sample
    sample_pos <- unique(df$position[df$sample == sample])
    
    for (position in sample_pos){ #each sample's variant called
      df$eval_sample[df$position==position & df$sample==sample] <- classify_variant_sample(position = position, 
                                                                                           sample = sample, 
                                                                                           df = df, 
                                                                                           true_variants = true_variants )
    }
    ## SAMPLE FNs: evaluate which positions haven't been found for each of the samples  
    ##             that are in the simulated record. These are  GATK misses  
    FalseNegatives <- setdiff(true_variants$position,    sample_pos)
    FalseNegatives<- setdiff(FalseNegatives, not_FN)
    
    if (length(FalseNegatives) != 0){
      for (fn in FalseNegatives ){
        missed <- true_variants[true_variants$position == fn, ]
        newrow <-  c(fn, 
                     unique(missed$ref), 
                     paste(missed$alt, collapse = "/"), 
                     NA,
                     sample, 
                     unique(df$ploidy), 
                     "","", "", "", 
                     NA,NA,NA,NA,NA,NA,NA,NA,NA,
                     "FN_uncalled",sample,NA,NA,NA,NA,NA)
        FNs<- rbind(FNs, newrow)
      }
    }
  } 
  print("Done with evaluation of all samples")
  colnames(FNs) <-  colnames(df)
  return(rbind(df, FNs))
} # Evaluate for given GATK output all samples' variants found
## Classify callable blacklisted #### 
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
## Classify rDNA region #### 
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
#############---------------------------- GET TABLES  -------------------------##################
d<- get_vc(ploidy)
d<- as.data.frame(d)
# Get genotype column decoded e.g.  for ploidy 5: A/A/A/T/C
d$GT <- get_GT(d)  
## Get evaluation sample level for GATK ploidies
d <- fill_eval_sample(df=d, true_variants=true_variants)
# Classify callable and blacklisted positions
d <-  get_black(d)
# Annotate the region for the uncalled positions
d<- add_regions_chrR(d)
# Save evaluation into a tbl file 
print(paste("Writting into", dirf))

d$type<- ifelse(length(d$REF)==length(d$ALT), "SNP", "INDEL")

if(step=="filtering"){
  write.table(d, paste0(dirf, "VCeval_postjoing_p", ploidy, "_", step,"_", version,".tbl"), sep="\t")
}else{
  write.table(d, paste0(dirf, "VCeval_postjoing_p", ploidy, "_", step, ".tbl"), sep="\t")
}
print("Done!")
print("------------------------------------")
