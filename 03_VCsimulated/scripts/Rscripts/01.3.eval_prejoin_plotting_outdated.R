#####  Library #####  
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)
#####  Palette  #####  
palette <-  c("TP" ="#75bb6f","TP_FN_allele"="#aeb852","TP_pos_other"="#706f4d",  "FN_uncalled"="#ffc60e", "FN_ref_homozygous"= "#4c89ee","FN_alt_homozygous"="#71a5be",
              "FP"= "#bd0752", "TP_pos_FP"= "#ff4d7a", "TP_pos_FP_allele"="#ed8c96","TP_allele_FP_allele"="#ed8c96",  "err"= "#a6a6a6") 
#####  Functions  #####  
# d[0<as.numeric(d$QUAL) & as.numeric(d$QUAL)<1000 ,]
get_counts_per_ploidy <- function(dfs){
  
  counts<- as.data.frame(table(dfs$ploidy, dfs$eval_sample))
  names(counts)<- c("ploidy", "eval_sample", "count")
  
  ploydies<- unique(counts$ploidy)
  # Get as frequencies as well
  ## remove 0 values so they dont appear in plot:
  counts <- counts[counts$count != 0,]
  
  counts$freq <- rep(0, nrow(counts))
  
  for (p in ploydies) {
    sum_val <- sum(counts$count[counts$ploidy == p])
    for (e in counts$eval_sample[counts$ploidy == p]) {
      counts$freq[counts$ploidy == p & counts$eval_sample == e] <- round(counts$count[counts$ploidy == p & counts$eval == e] / sum_val, 8)
    }
  }
  # Order the evaluation so they are descending in plot counts
  counts<-   counts %>%  arrange(count)
  
  counts$eval_sample <- factor(counts$eval_sample,  levels = c( "FP", "TP_pos_FP", "TP_pos_other", "TP_allele_FP_allele", "TP_FN_allele", 
                                                                "FN_alt_homozygous", "FN_ref_homozygous","FN_uncalled","TP" ))
  return(counts)
}
get_counts_per_sample <- function(dfs){
  
  df_counts<- as.data.frame(table(dfs$ploidy, dfs$sample, dfs$eval_sample))
  names(df_counts)<- c("ploidy","sample", "eval_sample", "count")
  # Get frequencies
  ploydies<- unique(df_counts$ploidy)
  eval<- unique(df_counts$eval_sample)
  samples<- paste0("sample_", 1:100)
  df_counts$freq <- 0
  
  for (p in ploydies) {
    for (sample in samples){
      sum_val <- sum(df_counts$count[df_counts$ploidy == p & df_counts$sample == sample])
      for (e in eval) {
        freq <- round(df_counts$count[df_counts$ploidy == p & df_counts$sample == sample & df_counts$eval_sample == e] / sum_val, 8)
        df_counts$freq[df_counts$ploidy == p & df_counts$sample == sample & df_counts$eval_sample == e] <- freq
      }
    }
  }
  # Order the evaluation types in descending order of counts
  df_counts <- df_counts %>% arrange(desc(count))
  
  # Remove rows with count = 0
  df_counts <- df_counts[df_counts$count != 0,]
  
  df_counts$eval_sample<- factor(df_counts$eval_sample,  levels =  c( "FP", "TP_pos_FP", "TP_pos_other", "TP_allele_FP_allele", "TP_FN_allele", 
                                                                      "FN_alt_homozygous", "FN_ref_homozygous","FN_uncalled","TP" ))
  return(df_counts)
}
## ------------------------------------ SET PARAMETERS --------------------------------- ############
# #### CHOOSE Join genotyping SET ######
args <- commandArgs(trailingOnly = TRUE)

#GATK version either "4.5.0.0" or "4.3.0.0"
set=args[1]
# set= "4.3.0.0"

# STEP of the pipeline this is for gVCF files so either "VC" or"filtering"
step=args[2]
# step= "VC"
# step= "filtering"

# version
version=args[3]
# version="v1"
# version="v1"
# version="v3"
# version="v4"
# version="5sample_v4"
# version="10sample_v4"
# version="5sample_v3"
# version="10sample_v3"
# version="5sample_v2"
# version="5sample_v3"
# version="10sample_v2"
# version="2sample_v4"
# version="2sample_v2"
# version="2sample_v3"
# #GATK version either "4.5.0.0" or "4.3.0.0"
## --------------------------------------  LOAD DATA ---------------------------------   ############

dir <-  "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/"

# dir <- "~/Desktop/ribosomal_RNAs/Alba/03_VCsimulated/v2/"

# Dataset containing all ploidies VCs
if(step=="filtering"){
  dirf<- paste0(dir, "VC_out/Sample_level_eval/4.3.0.0/VC/")
  dfs<- read.table(  paste0(dirf, "All_VCeval_prejoin_filtering_", version, ".tbl"), sep="\t")
  if (version %in% c("v1", "v2", "v3", "v4")){
        dirI<- paste0(dir, "IMAGES/01_PREJOIN/4.3.0.0/filtering_", version, "_all/")
  }else{
         dirI<- paste0(dir, "IMAGES/01_PREJOIN/4.3.0.0/filtering_sample/filtering_", version, "_all/")
  }
}else{
  dirf<- paste0(dir, "VC_out/Sample_level_eval/4.3.0.0/VC/")
  dirI<- paste0(dir, "IMAGES/01_PREJOIN/4.3.0.0/VC_all/")
  dfs<- read.table(  paste0(dirf, "All_VCeval_prejoin_VC.tbl"), sep="\t")
}

setwd(dirI)

# ploidies<-  c(2,5,7,10,20,30,40,50,100)
ploidies<- unique(dfs$ploidy)

print("#############################################################################")
print(paste("PLOTTING VC evaluation for", set, step, ploidies))


######### --------------------------  Preparing for plotting  ----------------------------- #############

dfs$eval_sample[dfs$eval_sample =="TP_pos_FN_allele"] <- "TP_FN_allele"

dfs_c <- dfs[dfs$class =="callable",]
# Get statistics ploidy level
dfs_c_counts_ploidy <- get_counts_per_ploidy(dfs_c)

# Per sample 
dfs_c_counts_sample <- get_counts_per_sample(dfs_c)

print("Done preparing data")

########## ---------------------------- Plotting only callable ----------------------------------------------------

# Per ploidy
png(file=paste0(dirI, "VCevaluation_callable_ploidy_prejoingen.png"), res=500, width=3000, height=1600);
ggplot(dfs_c_counts_ploidy, aes(x = factor(ploidy, levels = ploidies), y = freq, fill = eval_sample)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim=c(0.035, 0.97))+
  scale_fill_manual(values = palette)+ 
  geom_text_repel(aes(label = reorder(count, count)),  position = position_stack(vjust = 0.5), force = .005,size = 3, fontface = "bold") +
  labs(x = "Ploidy", y = "Frequency", fill="")+
  theme_classic()+
  theme(legend.margin=margin(0,0,0,-11), legend.key.size = unit(3.4, "mm" ), legend.text = element_text(size=8))
dev.off()

png(file=paste0(dirI, "VCevaluation_callable_samples_prejoingen.png"), res=300, width=3000, height=4000);
ggplot(dfs_c_counts_sample, aes(x = factor(sample, levels = paste0("sample_", 1:100)), y = count, fill = reorder(eval_sample, count))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette)+ 
  labs( x = "", y = "Number of variants", fill="")+
  guides(fill = guide_legend(nrow = 1, override.aes = list(size=3)))+
  facet_grid(ploidy~.)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 65,   size = 7.5, margin = margin(t = 16, r = 0, b = 0, l = 0)), 
        legend.margin=margin(-40,-10,0,0),     legend.text = element_text(size=9),
        legend.key.size = unit(2, "mm" ),   legend.position = "bottom")
dev.off()

# PER REGION 
t<-  as.data.frame(table(dfs_c$ploidy, dfs_c$region,  dfs_c$eval_sample))
colnames(t)<- c("ploidy", "region", "eval", "counts")
t<-  t[t$counts!=0,]
t$region<-  factor(t$region, c("5_ETS", "18S", "ITS1", "5.8S", "ITS2",  "28S", "3_ETS"))

png(file=paste0(dirI, "VCevaluation_region_prejoingen.png"), res=300, width=6000, height=1300);
ggplot(t, aes(x = region, y = counts, fill = reorder(eval, counts))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette)+ 
  guides(fill = guide_legend(nrow = 1, override.aes = list(size=1))) +
  facet_grid(.~ploidy)+
  geom_text_repel(aes(label = reorder(counts, counts)),  position = position_stack(vjust = 0.5), force = .004,size = 3, fontface = "bold") +
  labs( x = "rDNA region", y = "Number of variants", fill="")+
  theme_bw() +
  theme(legend.position = "bottom",    legend.text = element_text(size = 7),  
        legend.key.size = unit(4, "mm"),      legend.margin = margin(0, 10, 0, -30),     legend.box = "horizontal")
dev.off()

##### ----------------------------CALLABLE VS BLACKLISTED---------------------------- #######

dfs_b <- dfs[!dfs$class =="callable",]
# Get statistics ploidy level
dfs_b_counts_ploidy <- get_counts_per_ploidy(dfs_b)
# Get statistics sample level
dfs_b_counts_sample <- get_counts_per_sample(dfs_b)

dfs_c_counts_ploidy$region <- "callable"
dfs_b_counts_ploidy$region <- "blacklisted"

dfs_bc_ploidy<- rbind(dfs_c_counts_ploidy, dfs_b_counts_ploidy)

png(file=paste0(dirI, "VC_b_vs_c_prejoingen.png"), res=400, width=3800, height=1600);
ggplot(dfs_bc_ploidy, aes(x = factor(ploidy, levels = ploidies), y = freq, fill = reorder(eval_sample, count))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette)+ 
  geom_text_repel(aes(label = reorder(count, count)),  position = position_stack(vjust = 0.5), force = .01, size = 3, fontface = "bold") +
  labs(title ="", x = "Ploidy", y = "Frequency", fill="")+
  facet_grid(.~factor(region, levels = c("callable", "blacklisted")))+
  theme_classic()+
  theme(axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        legend.text = element_text(size=9), legend.position = "bottom", 
        legend.key.size = unit(3, "mm" ),
        legend.margin=margin(-10,0,0,-10))
dev.off()

print("All done!")
print("#############################################################################")
# png(file=paste0(dirI, "VCevaluation_blacklisted_samples_prejoingen.png"), res=200, width=2300, height=1300);
# ggplot(dfs_b_counts_sample, aes(x = factor(sample, levels = paste0("sample_", 1:100)), y = count, fill = reorder(eval_sample, count))) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = palette)+ 
#   labs(title = "Sample level evaluation blacklisted positions", x = "Samples", y = "Counts", fill="")+
#   guides(fill = guide_legend(nrow = 1, override.aes = list(size=3)))+
#   facet_grid(ploidy~.)+
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8), 
#         legend.text = element_text(size=9),legend.key.size = unit(4, "mm" ),  legend.margin=margin(-5,0,0,0), legend.position = "bottom")
# dev.off() 

###########################################################################################


###########################################################################################
##### ---------------------------- PLOT CALLABLE ploidy 100 --------------------------- #######
############################################################################################################
# d100<-  read.table(paste0(dirf, "VCeval_prejoin_p100_VC.tbl"), sep="\t")
# d100<-  read.table(paste0(dirf, "VCeval_prejoin_p100_.tbl"), sep="\t")
# 
# d100$eval_sample[d100$eval_sample =="TP_pos_FN_allele"] <- "TP_FN_allele"
# d <- rbind(dfs_c, d100[d100$class =="callable",])
# 
# dq1<- rbind(d[d$QUAL>1 , ], d[d$eval_sample=="FN_uncalled", ])
# # Get statistics ploidy level
# d_counts_ploidy <- get_counts_per_ploidy(d)
# dq1_counts_ploidy <- get_counts_per_ploidy(dq1)
# dq1_counts_ploidy$ploidy<- "100_QUAL1"
# 
# d100_counts_ploidy$eval_sample<- factor(d100_counts_ploidy$eval_sample, levels = c("TP","FP","FN_uncalled","TP_pos_FP",
#                                         "TP_allele_FP_allele","FN_alt_homozygous", "FN_ref_homozygous","TP_FN_allele",  "TP_pos_other"))
# d100q1_counts_ploidy$eval_sample<- factor(d100q1_counts_ploidy$eval_sample, levels = c("TP","FP","FN_uncalled","TP_pos_FP",
#                                        "TP_allele_FP_allele","FN_alt_homozygous", "FN_ref_homozygous","TP_FN_allele",  "TP_pos_other"))
# 
# png(file=paste0(dirI, "VCevaluation_callable_ploidy100_prejoingen.png"), res=300, width=1500, height=1500);
# ggplot(rbind(d100_counts_ploidy, d100q1_counts_ploidy), aes(x = ploidy, y = freq, fill = reorder(eval_sample, count))) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = palette)+ 
#   geom_text_repel(aes(label = reorder(count, count)),  position = position_stack(vjust = 0.5), force = .005,size = 3.3, fontface = "bold") +
#   labs(x = "Ploidy", y = "Frequency", fill="")+
#   theme_minimal()+
#   theme( legend.text = element_text(size=7))
# dev.off()
# 
# d100_counts_sample <- get_counts_per_sample(d100)
# d100_counts_sample$eval_sample<- factor(d100_counts_sample$eval_sample, levels = c("TP", "FP", "FN_uncalled","TP_pos_FP","TP_allele_FP_allele", 
#                                                                                      "FN_alt_homozygous", "FN_ref_homozygous", "TP_FN_allele",  "TP_pos_other"))
# 
# png(file=paste0(dirI, "VCevaluation_callable_ploidy100_samples_prejoingen.png"), res=200, width=2300, height=1000);
# ggplot(d100_counts_sample, aes(x = factor(sample, levels = paste0("sample_", 1:100)), y = count, fill = reorder(eval_sample, count))) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = palette)+ 
#   labs(title = "Sample level evaluation callable positions", x = "", y = "Counts", fill="")+
#   guides(fill = guide_legend(nrow = 1, override.aes = list(size=3)))+
#   facet_grid(ploidy~.)+
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 65,   size = 8, margin = margin(t = 13, r = 0, b = 0, l = 0)), 
#         legend.margin=margin(-40,-10,0,0),
#         legend.text = element_text(size=9),
#         legend.key.size = unit(3, "mm" ),
#         legend.position = "bottom")
# dev.off()


