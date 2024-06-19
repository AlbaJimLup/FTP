## ------------------------------------ SET PARAMETERS --------------------------------- ############
#### CHOOSE VC SET ######
args <- commandArgs(trailingOnly = TRUE)
#GATK version either "4.5.0.0" or "4.3.0.0"
set=args[1] 
# STEP of the pipeline this is for gVCF files so either "VC" or"filtering"
step=args[2]

ploidy=args[3]

version=args[4]

print(paste("Doing VC evaluation for", set, step, ploidy, version))
print("#############################################################################")

## libraries 
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)

root<- "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/"

# root<- "~/Desktop/ribosomal_RNAs/Alba/"
dir <- paste0(root,"03_VCsimulated/v2/")
dirf<- paste0(dir, "VC_out/Sample_level_eval/", set, "/VC/")
## --------------------------------------  LOAD DATA --------------------------------- ############
### Get true_variants ######
vcatalog <- read.table(paste0(dir, "scripts/Rscripts/variant_sharing.tbl"))
true_variants <- read.table(paste0(dir, "scripts/Rscripts/true_variants.tbl"))
not_FN<- c()
### Get bed callable blacklisted regions  ######
bed <- read.table(paste0(root, "/data/pre-rRNA_47S.included.bed"))
#############---------------------------- FUNCTIONS  -------------------------##################
## Get pre join gen g.vdf file output either HaplotypeCaller or VariantFilter ##
get_vc <-  function(p) { # each ploidy we ran so far
  
  ts <- data.frame()
  for (i in 1:100){ # each sample 
    for (region in c("5_ETS", "18S", "ITS1", "5.8S", "ITS2",  "28S", "3_ETS")){
      t <- read.table(gzfile(paste0(dir, "VC_out/hs1-rDNA_v1.0/simulated_WGS/rDNA_reads/", set,"/", step,"/sample_", i, "/ploidy_", p,
                                    "/sample_", i, ".ploidy_", p ,".", region,".g.vcf.gz" )))  
      
      t$sample<-paste0("sample_", i)
      
      t$ploidy<- p
      
      t$eval_sample<- ""
      
      t$region<- region
      
      # Remove "<NON_REF>" from all ALT rows
      t$V5 <- gsub(",<NON_REF>", "", t$V5)
      
      # Select only called positions
      ts <-  rbind(ts, t[t$V5 != "<NON_REF>", ]) # non-called but possible non-reference alleles
      #ts <-  rbind(ts, t)
      #print(paste0("Done with ", i))
    }
    
  }
  print(paste0("Done with prepering all samples for ploidy ", p))
  colnames(ts)<- c("CHROM",	"position",	  "ID", 	"REF",  	"ALT",
                   "QUAL", "FILTER",	"INFO",	"FORMAT",	"data", "sample", "ploidy", "eval_sample", "region")
  
  ts<- ts[ , c(2,4,5,6,11,12,13,14,9,10,8) ]
  return(ts)
}
#
get_vc_filtering <-  function(p, version="v1") { # each ploidy we ran so far
  
  ts <- data.frame()
  for (i in 1:100){ # each sample 
    for (region in c("5_ETS", "18S", "ITS1", "5.8S", "ITS2",  "28S", "3_ETS")){
      t <- read.table(gzfile(paste0(dir, "VC_out/hs1-rDNA_v1.0/simulated_WGS/rDNA_reads/", set,"/", step,"/sample_", i, "/ploidy_", p,
                                    "/sample_", i, ".ploidy_", p ,".", region,"_", step, ".g.vcf.gz" )))  
      
      t$sample<-paste0("sample_", i)
      t$ploidy<- p
      t$eval_sample<- ""
      t$region<- region
      
      # Remove "<NON_REF>" from all ALT rows
      t$V5 <- gsub(",<NON_REF>", "", t$V5)
      
      # Apply filter, remove those that are only filtered by QUAL<30 or not QUAL<1
      if (version=="v1"){ # keep QUAL>1
        t <- t[ t$V7 %in% c("PASS", "QUAL30", "QUAL30;QUAL400", "QUAL1000;QUAL30;QUAL400") ,  ]
      }
      else if(version=="v2"){ # keep QUAL>30
        t <- t[ t$V7 %in% c("PASS", "QUAL1000", "QUAL1000;QUAL400") ,  ]
      }
      else if(version=="v3"){ # keep QUAL>400
        t <- t[ t$V7 %in% c("PASS" , "QUAL1000"),  ]
      }
      else if(version=="v4"){ # keep QUAL>1000
        t <- t[ t$V7 == "PASS" ,  ]
      }
      # Select only called positions
      ts <-  rbind(ts, t[t$V5 != "<NON_REF>", ]) # non-called but possible non-reference alleles
      #ts <-  rbind(ts, t)
      #print(paste0("Done with ", i))
    }
  }
  print(paste0("Done with prepering all samples for ploidy ", p))
  colnames(ts)<- c("CHROM",	"position",	  "ID", 	"REF",  	"ALT",
                   "QUAL", "FILTER",	"INFO",	"FORMAT",	"data", "sample", "ploidy", "eval_sample", "region")
  
  ts<- ts[ , c(2,4,5,6,11,12,13,14,9,10,8) ]
  
  return(ts)
}
## Get GT column ##
get_GT<- function(df){
  GTs <- c()
  for (i in 1:nrow(df)){
    row<- df[i,]
    encoded <- strsplit(row$data, ":")[[1]][1]
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
classify_variant_sample <- function(position, sample, df, true_variants, vcatalog) {
  # print(position)
  position<- as.numeric(position)
  # if the variant is not one of the true callable variants for this position
  if (!position %in% true_variants$position[true_variants$sample == sample]) {
    
    called<- unique(unlist(strsplit(as.character(df[df$position == position & df$sample == sample, "GT"]), split = "[/|]")))
    if ((position+1) %in% true_variants$position[true_variants$sample == sample]) { # shifted e.g.
      # 12469:      GGCCCCG>G <- true_variants
      # 12468      GGGCCCC  G <- GATK
      # ref<- unique(as.character(df[df$position == position & df$sample == sample, "REF"]))
      # true_ref<- unique(true_variants[true_variants$position == position+1 & true_variants$sample == sample, "ref"])
      if (position==12468){
        if (length(setdiff(called, c("GGGCCCC", "G")))==0){
          return("TP")
          not_FN<- c(not_FN, position+1)
        }else{
          return("FP")
        }
      }else if(position==16595){
        if (length(setdiff(called, c("GCGTTTGT", "G")))==0){
          return("TP")
          not_FN<- c(not_FN, position+1)
        }else{
          return("FP")
        }
      }else if(position==13747){
        if (length(setdiff(called, c("GCTGA", "G")))==0){
          return("TP")
          not_FN<- c(not_FN, position+1)
        }else{
          return("FP")
        }
      }else if(position==13747){
        if (length(setdiff(called, c("GCTGA", "G")))==0){
          return("TP")
          not_FN<- c(not_FN, position+1)
        }else{
          return("FP")
        }
      }else if(position==22504){
        if (length(setdiff(called, c("CGCCG", "C")))==0){
          return("TP")
          not_FN<- c(not_FN, position+1)
        }else if (length(setdiff(called, c("CGCCG", "C")))==1){
          return("TP_allele_FP_allele")
          not_FN<- c(not_FN, position+1)
        }else{
          return("FP")
        }
      }else if(position==11544){
        if (length(setdiff(called, c("TCGGGTAC", "T")))==0){
          return("TP")
          not_FN<- c(not_FN, position+1)
        }else if (length(setdiff(called, c("CGCCG", "C")))==1){
          return("TP_allele_FP_allele")
          not_FN<- c(not_FN, position+1)
        }else{
          return("FP")
        }
      }else{
        return("FP")
      } 
    }else if ((position+2) %in% true_variants$position[true_variants$sample == sample]) { 
      if(position==14857){
        if (length(setdiff(called, c("GATC", "G")))==0){
          return("TP")
          not_FN<- c(not_FN, position+2)
        }else{
          return("FP")
        } 
      }else if(position==15998){
        if (length(setdiff(called, c("TAATGTG", "T")))==0){
          return("TP_allele_FP_allele")
          not_FN<- c(not_FN, position+2)
        }else{
          return("FP")
        }
      }else{
        return("FP")
      }
    }else{
      return("FP") # FP sample levelNon existent variant
    }
  } else {
    if (df[df$position == position & df$sample == sample, "GT"] == ".") {
      return("FN_uncalled")
    } else {
      alleles <- unique(unlist(strsplit(as.character(df[df$position == position & df$sample == sample, "GT"]), split = "[/|]")))
      # if only one allele
      ref_allele <- unique(df[df$position == position & df$sample == sample, "REF"])
      alt_allele <- unique(unlist(strsplit(as.character(df[df$position == position & df$sample == sample, "ALT"]), split = ",")))
      
      if (length(alleles) == 1) {
        if (position == 19715 & df[df$position == position & df$sample == sample, "ALT"] == "C") {
          return("TP")
        } else if (alleles == ref_allele) { # 1 allele & == REF
          return("FN_ref_homozygous")
        } else { # # 1 allele & == ALT
          return("FN_alt_homozygous")
        }
        # if heterozygous check that this sample is supposed to have the variant or not
      } else {
        true_alleles <- c(unique(true_variants[true_variants$position == position & true_variants$sample == sample, "ref"]),
                          sort(true_variants[true_variants$position == position & true_variants$sample == sample, "alt"]))
        
        gatk_alleles <- alleles
        # gatk_alleles <- c(unique(as.character(df[df$position == position & df$sample == sample, "REF"])),
        #                   sort(strsplit(as.character(df[df$position == position & df$sample == sample, "ALT"]), ",")[[1]]))
        # 
        condition1 <- length(true_alleles) == length(gatk_alleles) # same number of alleles
        condition2 <- sum(gatk_alleles %in% true_alleles) == length(gatk_alleles) # all gatk alleles TP
        
        condition3 <- sum(true_alleles[true_alleles != ref_allele] %in% # all TP alt alleles are identified +  some FP alt alleles
                            gatk_alleles[gatk_alleles != ref_allele]) == length(true_alleles[true_alleles != ref_allele])
        
        condition4 <- sum(true_alleles[true_alleles != ref_allele] %in%  # FP allele -> alt allele not in alt gatk alleles
                            gatk_alleles[gatk_alleles != ref_allele]) == 0
        
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
          }else if (position == 22309) {
            ifelse(sum(gatk_alleles %in% alt_allele) > 0,  return("TP"), return("TP_pos_FP"))
            
          }else {
            #### OTHER CASES ####
            if (true_alleles[1] != gatk_alleles[1]) {
              if ("*" %in% alt_allele) { # Deletions may be wrongly called
                # # 12469 GGCCCC  *,G,TGCCCC   0     sample_39    2  
                # # 12469:GGCCCCG>G   sample_39    indel   blacklisted      5'ETS 12469 GGCCCCG  G
                ref <- substring(true_alleles[1], 1, nchar(true_alleles[1]) - 1)
                
                if (ref == gatk_alleles[1]) {
                  return("TP_allele_FP_allele")
                }else if (position %in% c(18462,22505)){
                  return("TP_pos_FP")
                  # all ploidies call these positions the same way
                }else if (position %in% c(17347, 12399, 14552 )){ 
                  return("TP_allele_FP_allele")
                }else {
                  return("TP_pos_other")
                }
                
              }else if (position == 17340) {
                return("TP_pos_FP")
              }else if (position == 11830) {
                ifelse((ref_allele == "TGCG" & ("CGCG" %in% alt_allele)),  return("TP_allele_FP_allele"), return("TP_pos_FP"))
                
              }else if (position == 11651) {
                ifelse(ref_allele == "GC" & ("TC" %in% alt_allele), return("TP_allele_FP_allele"), return("TP_pos_FP"))
                
              }else if (position == 17347) {
                ifelse(ref_allele == "CG" & ("AC" == alt_allele), return("TP_allele_FP_allele"), return("TP_pos_FP"))
                
              }else {
                print(paste0("err: ", "diff_REF", " - ", position))
                return("TP_pos_other")
              }
            } else {
              return("TP_pos_FP") ## Same ref different alleles
            }
          }
        }
      }
    } 
  }
} # Comparison sample level 
#
fill_eval_sample <- function(df, true_variants, vcatalog){
  
  FNs <- data.frame(matrix(ncol=ncol(df), nrow=0)) # initialize missed variants
  samples <-paste0("sample_", 1:100) # num simulated samples
  
  for (sample in samples){ # each sample
    sample_pos <- unique(df$position[df$sample == sample])
    
    for (position in sample_pos){ #each sample's variant called
      df$eval_sample[df$position==position & df$sample==sample] <- classify_variant_sample(position = position, 
                                                                                           sample = sample, 
                                                                                           df = df, 
                                                                                           true_variants = true_variants, 
                                                                                           vcatalog = vcatalog)
    }
    ## SAMPLE FNs: evaluate which positions haven't been found for each of the samples  
    ##             that are in the simulated record. These are  GATK misses  
    FalseNegatives <- setdiff(true_variants$position[ true_variants$sample == sample],    sample_pos)
    
    FalseNegatives<- setdiff(FalseNegatives, not_FN)
    if (length(FalseNegatives) != 0){
      for (fn in FalseNegatives ){
        missed <- true_variants[true_variants$position == fn  & true_variants$sample==sample, ]
        newrow <-  c(fn, 
                     unique(missed$ref), 
                     paste(missed$alt, collapse = "/"), 
                     0,
                     sample, 
                     unique(df$ploidy), 
                     "FN_uncalled", 
                     "","", "", "", "")
        FNs<- rbind(FNs, newrow)
      }
    }
  } 
  print("Done with evaluation of all samples")
  colnames(FNs) <-  colnames(df)
  return(rbind(df, FNs))
} # Evaluate for given GATK output all samples' variants found
## Classify callable blacklisted  
get_black<-  function(df){
  
  for (pos in unique(df$position)){
    
    eval_black <-  unique(ifelse(pos > bed$V2 & pos <= bed$V3, "YES", "NO"))
    
    if(length(unique(eval_black)) == 1){
      df$class[df$position == pos] <- "blacklisted"
    }
    else{
      df$class[df$position == pos] <- "callable"
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
## Splice data 
get_annotations <-  function(p) { # each ploidy 
  ## Create INFO columns
  info_cols<- c( "BaseQRankSum","DP", "ExcessHet","MLEAC","MLEAF","MQRankSum","RAW_MQandDP","ReadPosRankSum")
  
  for (cols in info_cols  ){
    p[[cols]] <- gsub(paste0(".*", cols, "=([^;]+).*"), "\\1", p$INFO)
  }
  ## Create GT:AD:DP:GQ:PL columns
  format_cols <- c("GT", "AD", "DP", "GQ", "PL", "SB")
  new_cols <- matrix(nrow = length(p$FORMAT), ncol = length(format_cols))
  colnames(new_cols) <- format_cols
  
  for (i in seq_along(p$FORMAT)) {
    row <- unlist(strsplit(p$data[i], ":"))
    new_cols[i, ] <- row[1:length(format_cols)]
  }
  ##  Add region to dataframe 
  t<- cbind(p, new_cols) # Add to t
  # Remove void or compacted info columns: 
  print("We decoded additional annotations information info different columns")
  return(as.data.frame(t[, ! colnames(t) %in% c("INFO","FORMAT")]))
}
#### VC evaluation sample level #### 
# classify_variant_sample <- function(position, sample, df, true_variants, vcatalog) {
  # print(position)
#   # if the variant is not one of the true callable variants for this position
#   if (!position %in% true_variants$position[true_variants$sample == sample]) {
#     return("FP") # FP sample levelNon existent variant
#   } else {
#     if (df[df$position == position & df$sample == sample, "GT"] == ".") {
#       return("FN_uncalled")
#     } else {
#       alleles <- unique(unlist(strsplit(as.character(df[df$position == position & df$sample == sample, "GT"]), split = "[/|]")))
#       # if only one allele
#       ref_allele <- unique(df[df$position == position & df$sample == sample, "REF"])
#       alt_allele <- unique(unlist(strsplit(as.character(df[df$position == position & df$sample == sample, "ALT"]), split = ",")))
#       
#       if (length(alleles) == 1) {
#         if (position == 19715 & df[df$position == position & df$sample == sample, "ALT"] == "C") {
#           return("TP")
#         } else if (alleles == ref_allele) { # 1 allele & == REF
#           return("FN_ref_homozygous")
#         } else { # # 1 allele & == ALT
#           return("FN_alt_homozygous")
#         }
#         # if heterozygous check that this sample is supposed to have the variant or not
#       } else {
#         true_alleles <- c(unique(true_variants[true_variants$position == position & true_variants$sample == sample, "ref"]),
#                           sort(true_variants[true_variants$position == position & true_variants$sample == sample, "alt"]))
#         
#         gatk_alleles <- c(unique(as.character(df[df$position == position & df$sample == sample, "REF"])),
#                           sort(strsplit(as.character(df[df$position == position & df$sample == sample, "ALT"]), ",")[[1]]))
#         
#         condition1 <- length(true_alleles) == length(gatk_alleles) # same number of alleles
#         condition2 <- sum(gatk_alleles %in% true_alleles) == length(gatk_alleles) # all gatk alleles TP
#         
#         condition3 <- sum(true_alleles[true_alleles != ref_allele] %in% # all TP alt alleles are identified +  some FP alt alleles
#                             gatk_alleles[gatk_alleles != ref_allele]) == length(true_alleles[true_alleles != ref_allele])
#         
#         condition4 <- sum(true_alleles[true_alleles != ref_allele] %in%  # FP allele -> alt allele not in alt gatk alleles
#                             gatk_alleles[gatk_alleles != ref_allele]) == 0
#         
#         if (condition1 & condition2) { # same number of alleles & same alleles
#           return("TP")
#         } else {
#             if (!condition1 & condition2) { # Multi-allelic variants we are missing alleles  
#               return("TP_pos_FN_allele") # all gatk alleles are true alleles but some true alleles were not identified 
#             } else if (!condition1 & condition3) {
#               return("TP_allele_FP_allele") # Found variant + FP some false allele were identified
#             } else if (condition4) { # Diverging alleles even if correct position
#               if (true_alleles[1] == gatk_alleles[1]) {
#                 return("TP_pos_FP")
#               }else {
#                 if ("*" %in% alt_allele) { # Deletions may be wrongly called
#                   ref <- substring(true_alleles[1], 1, nchar(true_alleles[1]) - 1)
#                   # # 12469 GGCCCC  *,G,TGCCCC   0     sample_39    2  
#                   # # 12469:GGCCCCG>G   sample_39    indel   blacklisted      5'ETS 12469 GGCCCCG  G
#                   ifelse(ref == gatk_alleles[1], return("TP_allele_FP_allele"), return("TP_pos_FP"))
#               }else {
#                 return("TP_pos_FP")
#               }
#             }
#           }else if (position == 22309) {
#             ifelse(sum(gatk_alleles %in% alt_allele) > 0,  return("TP"), return("TP_pos_FP"))
#           
#           }else {
#             # OTHER CASES #
#             if (true_alleles[1] != gatk_alleles[1]) {
#               if ("*" %in% alt_allele) { # Deletions may be wrongly called
#                 # # 12469 GGCCCC  *,G,TGCCCC   0     sample_39    2  
#                 # # 12469:GGCCCCG>G   sample_39    indel   blacklisted      5'ETS 12469 GGCCCCG  G
#                 ref <- substring(true_alleles[1], 1, nchar(true_alleles[1]) - 1)
#                 
#                 if (ref == gatk_alleles[1]) {
#                   return("TP_allele_FP_allele")
#                 }else if (position %in% c(18462,22505)){
#                   return("TP_pos_FP")
#                 # all ploidies call these positions the same way
#                 }else if (position %in% c(17347, 12399, 14552 )){ 
#                   return("TP_allele_FP_allele")
#                 }else {
#                   return("TP_pos_other")
#                 }
#                 
#               }else if (position == 17340) {
#                 return("TP_pos_FP")
#               }else if (position == 11830) {
#                 ifelse((ref_allele == "TGCG" & ("CGCG" %in% alt_allele)),  return("TP_allele_FP_allele"), return("TP_pos_FP"))
#                 
#               }else if (position == 11651) {
#                 ifelse(ref_allele == "GC" & ("TC" %in% alt_allele), return("TP_allele_FP_allele"), return("TP_pos_FP"))
#               
#               }else if (position == 17347) {
#                 ifelse(ref_allele == "CG" & ("AC" == alt_allele), return("TP_allele_FP_allele"), return("TP_pos_FP"))
#                 
#               }else {
#                   print(paste0("err: ", "diff_REF", " - ", position))
#                   return("TP_pos_other")
#             }
#           } else {
#             return("TP_pos_FP") ## Same ref different alleles
#             }
#          }
#         }
#       }
#     } 
#   }
# } # Comparison sample level 
#############---------------------------- GET TABLES  -------------------------##################

if(step =="filtering"){ 
  d<- get_vc_filtering(ploidy, version)
}else {
  d<- get_vc(ploidy)
}
d<- as.data.frame(d)

version =paste0("2sample_", version)
# ## Finish to implement
# if (version =="2sample"){
t<-as.data.frame(table(d$position))
s2<- t$Var1[t$Freq>2]
### keep onlt those with more thank 2 samples with that position
print(paste("Remove variants present in less tha 2 samples, going from ", dim(d)[1] ))
d<- d[d$position %in% s2,  ]
print(paste("To ", dim(d)[1] ))
# }

# Get genotype column decoded e.g.  for ploidy 5: A/A/A/T/C
d$GT <- get_GT(d)  
## Get evaluation sample level for GATK ploidies
d <- fill_eval_sample(df=d, true_variants=true_variants, vcatalog)
# Classify callable and blacklisted positions
d <-  get_black(d)
# Annotate the region for the uncalled positions
d[d$eval_sample =="FN_uncalled",]<- add_regions_chrR(d[d$eval_sample =="FN_uncalled",])
# Decode other annotated information into different columns
d <- get_annotations(d)
# Save evaluation into a tbl file 
print(paste("Writting into", dirf))

#############---------------------------- SAVE TABLE  -------------------------##################
if(step=="filtering"){
  write.table(d, paste0(dirf, "VCeval_prejoin_p", ploidy, "_", step,"_", version,".tbl"), sep="\t")
}else{
  write.table(d, paste0(dirf, "VCeval_prejoin_p", ploidy, "_", step, ".tbl"), sep="\t")
}


print("Done!")
print("------------------------------------")
