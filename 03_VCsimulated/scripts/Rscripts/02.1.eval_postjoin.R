## libraries 
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)

#### CHOOSE VC SET ######
args <- commandArgs(trailingOnly = TRUE)
#GATK version either "4.5.0.0" or "4.3.0.0"
set=args[1] 
# STEP of the pipeline this is for gVCF files so either "VC" or"filtering"
step=args[2]
# ploidy gatk set
ploidy=args[3]
# version of the filter to be applyed
version=args[4]

# set="4.3.0.0"
# step="joingen"
# step="filtering"
# ploidy=7
# version="v3"
# ploydies<- c(2,5,7,10,20,50,100)

root<- "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/"
# root<- "~/Desktop/ribosomal_RNAs/Alba/"

dir <- paste0(root,"03_VCsimulated/v2/")
dirf<- paste0(dir, "VC_out/hs1-rDNA_v1.0/simulated_WGS/rDNA_reads/", set, "/", step , "/", "ploidy_", ploidy)
dirout<- paste0(dir, "VC_out/Sample_level_eval/", set, "/joingen/")

print("#############################################################################")
print(paste("Performing join genotyping  evaluation for", set, step, version ,ploidy))
print("#############################################################################")
## --------------------------------------  LOAD DATA --------------------------------- ############
#Get true_variants 
vcatalog <- read.table(paste0(dir, "scripts/Rscripts/variant_sharing.tbl"))
true_variants <- read.table(paste0(dir, "scripts/Rscripts/true_variants.tbl"))
not_FN<- c()
#Get bed callable blacklisted regions 
bed <- read.table(paste0(root, "/data/pre-rRNA_47S.included.bed"))
#
#############------------------------- FUNCTIONS-------------------------##################
# Get VC files
get_vc <-  function(p=ploidy) { # each ploidy 
  
  ts <- data.frame() # table of all regions to be returned for eahc ploidy
  variants <- data.frame(matrix(ncol=3, nrow=0)) #table with reported number of variants for each region and ploidy
  
  for (region in c("5_ETS", "18S", "ITS1", "5.8S", "ITS2",  "28S", "3_ETS")){
    
    t_raw <- read.table(gzfile(paste0(dirf, "/simulated_WGS.ploidy_", p ,".", region,".vcf.gz" )))  
    
    colnames(t_raw)<- c("CHROM",	"position",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT", paste0("sample_", sort(as.character(1:100))))
    
    ## Create INFO columns
    info_cols<- c("AC","AF","AN","BaseQRankSum","DP","FS","MLEAC","MLEAF","MQ","MQRankSum","QD","ReadPosRankSum","SOR")
    
    for (cols in info_cols  ){
      t_raw[[cols]] <- gsub(paste0(".*", cols, "=([^;]+).*"), "\\1", t_raw$INFO)
    }
    # Initialize ploidy region and evaluation result
    t_raw$ploidy<- p
    t_raw$eval_sample<- ""
    t_raw$region<- region
    
    ## Samples to long format
    t <- t_raw %>% pivot_longer(cols = "sample_1":"sample_99", 
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
    
    ts <-  rbind(ts, t)
    variants[nrow(variants)+1,]<- c(p, region, nrow(t))
    print(paste0("Fetched with ploidy ", p, "   region ", region, "   with ", length(t$position), " called variants"))
  }
  # colnames(ts)<- c("CHROM",	"position",	  "ID", 	"REF",  	"ALT",
  #                  "QUAL", "FILTER",	"INFO",	"FORMAT",	"data", "sample", "ploidy", "eval_sample", "region")
  # ts<- ts[ , c(2,4,5,6,11,12,13,9,10,8) ]
  # 
  return(list(data=as.data.frame(ts), report=variants))
}
get_vc_filtering <-  function(p=ploidy, version=v1) {
  d <- read.table(gzfile(paste0(dirf, "/joingen_Rfiltering_", version, "_p", ploidy, ".tbl"))) 
  print(paste0("Fetched filtered file for ploidy ", ploidy, " with ", length(d$position), " called variants"))
  return(d)
} 
# Get GT column
get_GT<- function(d){
  GTs <- c()
  for (i in 1:nrow(d)){
    row<- d[i,]
    encoded <- row$GT
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
  print("Resconstructed GT column")
  return(GTs)
}
# Comparison sample level 
classify_variant_sample <- function(position, sample, df, true_variants, vcatalog) {
  # print(position)
  called<- unique(unlist(strsplit(as.character(df[df$position == position & df$sample == sample, "GT"]), split = "[/|]")))
   
  # if the variant is not one of the true callable variants for this position
  if (!position %in% true_variants$position[true_variants$sample == sample]) {
    if (position %in% vcatalog$position){
      return("TP_pos_FP_sample")
    }else if ((position+1) %in% true_variants$position[true_variants$sample == sample]) { # shifted e.g.
      # 12469:      GGCCCCG>G <- true_variants
      # 12468      GGGCCCC  G <- GATK
      # ref<- unique(as.character(df[df$position == position & df$sample == sample, "REF"]))
      # true_ref<- unique(true_variants[true_variants$position == position+1 & true_variants$sample == sample, "ref"])
      if (position==12468){ #
        if (length(setdiff(called, c("GGGCCCC", "G")))==0){
          return("TP")
          not_FN<- c(not_FN, position+1)
        }else{
          return("FP")
        }
      }else if(position==17668){ #
        if (length(setdiff(called, c("AGGGCGT", "A")))==0){
          return("TP")
          not_FN<- c(not_FN, position+1)
        }else{
          return("FP")
        }
      }else if(position==16595){ #
        if (length(setdiff(called, c("GCGTTTGT", "G")))==0){
          return("TP")
          not_FN<- c(not_FN, position+1)
        }else{
          return("FP")
        }
      }else if(position==17128){ #
        if (length(setdiff(called, c("C", "G")))==0){
          return("TP_pos_FP")
          not_FN<- c(not_FN, position+1)
        }else{
          return("FP")
        }
      }else if(position==13747){ #
        if (length(setdiff(called, c("GCTGA", "G")))==0){
          return("TP")
          not_FN<- c(not_FN, position+1)
        }else{
          return("FP")
        }
      }else if(position==22504){#
        if (length(setdiff(called, c("CGCCG", "C")))==0){
          return("TP")
          not_FN<- c(not_FN, position+1)
        }else{
          return("FP")
        }
      }else if(position==11544){ #
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
      if(position==15998){ #
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
    # print(paste(df[df$position == position & df$sample == sample, "GT"], position))
    if (df[df$position == position & df$sample == sample, "GT"]==".") {
      return("FN_uncalled")
    } else {
      alleles <- unique(unlist(strsplit(as.character(df[df$position == position & df$sample == sample, "GT"]), split = "/")))
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
        gatk_alleles<- alleles
        # gatk_alleles <- c(unique(as.character(df[df$position == position & df$sample == sample, "REF"])),
        #                   sort(strsplit(as.character(df[df$position == position & df$sample == sample, "ALT"]), ",")[[1]]))
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
                    # # 12469 GGCCCC  *,G,TGCCCC   0     sample_39    2  
                    # # 12469:GGCCCCG>G   sample_39    indel   blacklisted      5'ETS 12469 GGCCCCG  G
                    ref <- substring(true_alleles[1], 1, nchar(true_alleles[1]) - 1)
                    
                    if (ref == gatk_alleles[1]) {
                      return("TP_allele_FP_allele")
                    } else {
                      return("TP_pos_FP")
                    }
                  } else {
                    print(paste0("Check out position ", position, "   for sample ",sample))
                    return("TP_pos_FP")
                  }
                }
          } else if (position == 22309) {
            ifelse(sum(gatk_alleles %in% alt_allele) > 0,     return("TP"), return("FP"))
          } else {
            ### OTHER CASES ###
            if (true_alleles[1] != gatk_alleles[1]) {
              if ("*" %in% alt_allele) { # Deletions may be wrongly called
                # # 12469 GGCCCC  *,G,TGCCCC   0     sample_39    2  
                # # 12469:GGCCCCG>G   sample_39    indel   blacklisted      5'ETS 12469 GGCCCCG  G
                ref <- substring(true_alleles[1], 1, nchar(true_alleles[1]) - 1)
                
                if (ref == gatk_alleles[1]) {
                  return("TP_allele_FP_allele")
                } else {
                  return("TP_pos_other")
                }
              } else if (position == 17340) {
                return("TP_pos_FP")
              } else if (position == 17347) {
                return("TP_allele_FP_allele")
              } else if (position == 11830) {
                if (ref_allele == "TGCG" & ("CGCG" %in% alt_allele)) return("TP_allele_FP_allele")
                else {
                  return("TP_pos_FP")
                }
              } else if (position == 11651) {
                if (ref_allele == "GC" & ("TC" %in% alt_allele)) return("TP_allele_FP_allele")
                else {
                  return("TP_pos_FP")
                }
              } else {
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
# Evaluate for given GATK output all samples' variants found
fill_eval_sample <- function(df, true_variants, vcatalog){
  
  
  FNs <- data.frame() # initialize missed variants
  samples <-paste0("sample_", 1:100) # num simulated samples
  
  for (sample in samples){ # each sample
    sample_pos <- unique(df$position[df$sample == sample])
    
    ## Evaluate positions GATK called
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
        
        newrow <-  c(
          fn,  # this is POS
          unique(missed$ref),  # this is REF
          paste(missed$alt, collapse = ","),  # this is ALT
          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA , NA ,  # these are QUAL, AC, AF, AN, BaseQRankSum, DP, FS, MLEAC, MLEAF, MQ, MQRankSum, QD, ReadPosRankSum, SOR
          unique(df$ploidy), # this is ploidy
          "FN_uncalled",  # this is eval_sample
          "unknown",  # this is region
          sample,  # this is sample
          NA, NA, NA, NA, NA, "" ) # these are GT, AD, DP.1, GQ, PL
        
        FNs<- rbind(FNs, newrow)
      }
    }
    colnames(FNs) <-  colnames(df)
  } 
  print("Evaluatied all samples")
  return(rbind(df, FNs))
}
# Classify callable blacklisted regions 
get_black<-  function(df){
  
  df$class<- "" #initialize
  
  for (pos in unique(df$position)){
    
    eval_black <-  unique(ifelse(pos > bed$V2 & pos <= bed$V3, "YES", "NO"))
    
    if(length(unique(eval_black)) == 1){ #If all "NO"
      df$class[df$position == pos] <- "blacklisted"
    }
    else{ # IF one "YES"
      df$class[df$position== pos] <- "callable"
    }
  }
  print("We classified callable and blacklisted regions")
  return(df)
}
## Classify region
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
## Get ADs columns
get_prep<- function(d){ 
  # keep callable
  # d<-  d[d$class=="callable", ] 
  ## Split allele depth ##
  d$ad_ref<- ""
  d$ad_alt<- ""
  d$ad_alt2<- ""
  d$ad_alt3<- ""
  # Fill the AD columns  
  for (position in unique(d$position)){
    # each sample has the same ad values
    for (sample in d$sample[d$position == position]){# the number of reads that support each of the reported alleles.
      #print(paste0("Trying with postion ", position, " and ", sample))
      if (d$eval_sample[d$position==position & d$sample==sample ] != "FN_uncalled"){
        # Get the alleles present in the sample
        gts<-  as.character(d$encoded[d$position==position & d$sample==sample ])
        
        # Get which alleles are present in the called genotypes
        gts <- as.numeric(unique(unlist(ifelse( grepl("/", gts),  strsplit(gts, "/"), strsplit(gts, "\\|") ))))
        
        # Get allele depths 
        depths <-  as.numeric(unlist(strsplit(d$AD[d$position==position & d$sample==sample ], ",")))
        
        ref <-  depths[gts[1]+1]
        alt <-  depths[gts[2]+1]
        
        nalleles<- length(depths)
        
        if(nalleles==2){
          d[d$position==position & d$sample==sample, c("ad_ref", "ad_alt")] <-  c(ref, alt)
          
        }else if(nalleles==3){
          d[d$position==position & d$sample==sample, 
            c("ad_ref", "ad_alt", "ad_alt2")] <-  c(ref, alt, depths[gts[3]+1])
          
        }else if(nalleles==4){
          #print(paste0("There are 4 alleles ", position))
          d[d$position==position & d$sample==sample, 
            c("ad_ref", "ad_alt", "ad_alt2", "ad_alt3")] <- c(ref, alt, depths[gts[3]+1], depths[gts[4]+1])
          
        }else{
          print(c(depths, gts))
          print(paste0("More than 4 alleles in position ", position))
        }
      }
    }
  }
  return(d)
}

#############------------------------- Perform evaluation -------------------------##################
## Get join gen vcf files
if(step=="joingen"){
  d<- get_vc(ploidy)
  d<- as.data.frame(d$data)
  # report<-rbind(report , d$report)
}else{
  d<- get_vc_filtering(ploidy, version)
  d<- as.data.frame(d)
}
## Get GT column 
d$encoded<- d$GT#keep encoded to get allele frequency in analysis
d$GT <- get_GT(d)
## VC evaluation sample level, get evaluation  current ploidy
d<- fill_eval_sample(df=d, true_variants=true_variants, vcatalog)
## Classify callable blacklisted
d <-  get_black(d)
## Classify region
d[d$eval_sample =="FN_uncalled",]<- add_regions_chrR(d[d$eval_sample =="FN_uncalled",])
##  Get ADs columns
d<-  get_prep(d)

print(paste("Writting file into", dirout))

if (step=="filtering"){
  write.table(d, paste0(dirout,  "VCeval_postjoin_p", ploidy, "_", step,"_", version, ".tbl"), sep="\t")
}else{
  write.table(d, paste0(dirout,  "VCeval_postjoin_p", ploidy, "_", step,".tbl"), sep="\t")
}

print(paste("Done with ", set, step, version, " of ploidy ", ploidy))

# ###########------------------------- Plot variants found per region  ------------------------------- #############
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



