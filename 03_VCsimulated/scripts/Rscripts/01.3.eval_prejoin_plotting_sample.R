#####  Library #####  
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)
library(stringr)
#####  Palette  #####  
palette <-  c("TP" ="#75bb6f","TP_FN_allele"="#ccB852", "TP_pos_FN_allele"="#d1cf6b" ,"TP_pos_other"="#706f4d",  
              "FN_uncalled"="#ffc60e", "FN_ref_homozygous"= "#4c89ee","FN_alt_homozygous"="#71a5be",
              "FP"= "#bd0752", "TP_pos_FP"= "#ff4d7a", "TP_pos_FP_allele"="#ed8c96","TP_allele_FP_allele"="#ed8c96",  "err"= "#a6a6a6") 
######  Functions  #####  
# d[0<as.numeric(d$QUAL) & as.numeric(d$QUAL)<1000 ,]
get_counts_per_ploidy <- function(dfs){
  
  counts<- as.data.frame(table(dfs$filter ,dfs$eval_sample))
  names(counts)<- c("filter", "eval_sample", "count")
  
  ploydies<- unique(counts$filter)
  # Get as frequencies as well
  ## remove 0 values so they dont appear in plot:
  counts <- counts[counts$count != 0,]
  
  counts$freq <- rep(0, nrow(counts))
  
  for (p in ploydies) {
    sum_val <- sum(counts$count[counts$filter == p])
    for (e in counts$eval_sample[counts$filter == p]) {
      counts$freq[counts$filter == p & counts$eval_sample == e] <- round(counts$count[counts$filter == p & counts$eval == e] / sum_val, 8)
    }
  }
  # Order the evaluation so they are descending in plot counts
  counts<-   counts %>%  arrange(count)
  
  counts$eval_sample <- factor(counts$eval_sample,  levels = c("FN_alt_homozygous", "FN_ref_homozygous","TP_FN_allele","TP_pos_FN_allele", "FN_uncalled",
                                                               "FP", "TP_pos_FP", "TP_pos_other", "TP_allele_FP_allele", "TP" ))
  return(counts)
}
## ------------------------------------ SET PARAMETERS --------------------------------- ############
#GATK version either "4.5.0.0" or "4.3.0.0"
set= "4.3.0.0"
# STEP of the pipeline this is for gVCF files so either "VC" or"filtering"
step= "filtering"
dir<-  "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/"
# dir <- "~/Desktop/ribosomal_RNAs/Alba/03_VCsimulated/v2/"
dirf<- paste0(dir, "VC_out/Sample_level_eval/4.3.0.0/VC/")
dirI<- paste0(dir, "IMAGES/01_PREJOIN/4.3.0.0/filtering_sample/")
#
print("#############################################################################")
print(paste("PLOTTING VC evaluation for minimum common samples filter", set, step))
# ## --------------------------------------  LOAD DATA ---------------------------------   ############
# dfs<- data.frame()
# 
# for (ploidy in c(50, 100)){
#   for (version in 2:4){
#     for(filter in c(2,5,10,15)){
#       d<- read.table(paste0(dirf, "VCeval_prejoin_p", ploidy, "_filtering_", filter, "sample_v", version, ".tbl"), header = T)
#       d$filter<- paste0( ploidy," ",filter, "samples ", ifelse( version==2, "QUAL30", ifelse(version==3, "QUAL400", "QUAL1000")))
#       print(paste0( ploidy," ",filter, "samples ", ifelse( version==2, "QUAL30", ifelse(version==3, "QUAL400", "QUAL1000"))))
#       dfs<- rbind(dfs, d)
#     }
#   }
# }
# 
# ######### --------------------------  Preparing for plotting  ----------------------------- #############
# dfs_c <- dfs[dfs$class =="callable",]
# levels<-c("50 2samples QUAL30",  "50 2samples QUAL400",  "50 2samples QUAL1000",
#           "50 5samples QUAL30",  "50 5samples QUAL400",  "50 5samples QUAL1000",
#           "50 10samples QUAL30", "50 10samples QUAL400", "50 10samples QUAL1000",
#           "50 15samples QUAL30",  "50 15samples QUAL400",  "50 15samples QUAL1000",
#           
#           "100 2samples QUAL30",  "100 2samples QUAL400",  "100 2samples QUAL1000",
#           "100 5samples QUAL30",  "100 5samples QUAL400",  "100 5samples QUAL1000",
#           "100 10samples QUAL30", "100 10samples QUAL400", "100 10samples QUAL1000",
#           "100 15samples QUAL30",  "100 15samples QUAL400",  "100 15samples QUAL1000")
# # Get statistics ploidy level
# dfs_c_counts_ploidy <- get_counts_per_ploidy(dfs_c)
# 
# dfs_c_counts_ploidy$eval_sample<- factor(dfs_c_counts_ploidy$eval_sample)
# # levels(dfs_c_counts_ploidy$eval_sample)= c("TP_pos_FP","TP_pos_FN_allele","FP","FN_uncalled","TP" )
# # Extract ploidy from the filter column
# dfs_c_counts_ploidy$ploidy <- word(dfs_c_counts_ploidy$filter, 1)
# dfs_c_counts_ploidy$num_samples <- factor(word(dfs_c_counts_ploidy$filter, 2))
# # levels(dfs_c_counts_ploidy$num_samples ) =  c("2samples", "5samples", "10samples", "15samples"))
# dfs_c_counts_ploidy$qual <- factor(word(dfs_c_counts_ploidy$filter, 3))
# # , levels =  c("QUAL400","QUAL30",  "QUAL1000"))
# write.table(dfs_c_counts_ploidy, paste0(dirf, "VC_filtered_sample.tbl"),sep=",", quote=F, row.names=F, col.names=T)
# 
# print("Done preparing data")
# ########## ---------------------------- Plotting only callable ----------------------------------------------------
# 
# p<- ggplot(dfs_c_counts_ploidy, aes(x = num_samples, y = freq, label =count,fill =factor(eval_sample, levels =c("TP_pos_FP","TP_pos_FN_allele","FP","FN_uncalled","TP" )))) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = palette)+
#   facet_grid(ploidy~qual, scales = "free", space = "free")+
#   geom_text_repel(data = subset(dfs_c_counts_ploidy, freq > 0.4),   position = position_stack(vjust = 0.5), force = 0,size= 3, fontface = "bold") +
#   geom_text_repel(data = subset(subset(dfs_c_counts_ploidy, freq > 0.05), freq < 0.4),  
#                   nudge_y = 0.8, force = .01,size = 3, fontface = "bold",  segment.color = NA) +
#   geom_text_repel(data = subset(dfs_c_counts_ploidy, freq < 0.05), 
#                   nudge_y = 1.1, force = .005,size = 3, fontface = "bold",  segment.color = NA) +
#   labs(x = "Number of samples for a position to be considered", y = "Frequency", fill="")+
#   theme_classic()+
#   theme( axis.text.x = element_text(angle=45 , hjust=1.1), legend.position = "bottom",
#         legend.margin=margin(0,0,0,-11), legend.key.size = unit(3.4, "mm" ), legend.text = element_text(size=8))
# ggsave(filename=paste0(dirI,"VCevaluation_callable_filter_ploidy_labels.png"), plot=p, dpi=400, width=2300/400, height=2500/400)
# 
# 
# ggplot(dfs_c_counts_ploidy, aes(x = factor(num_samples, levels = c("2samples", "5samples", "10samples", "15samples")), y = freq, fill =factor(eval_sample, levels =c("TP_pos_FP","TP_pos_FN_allele","FP","FN_uncalled","TP" ))))+
#   geom_bar(stat = "identity", position = "stack") +
#   scale_fill_manual(values = palette) +
#   facet_grid(ploidy ~ factor(qual, levels =c("QUAL30", "QUAL400", "QUAL1000") ), scales = "free", space = "free") +
# 
#   labs(x = "Ploidy", y = "Frequency", fill = "") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1.1), 
#         legend.position = "bottom",
#         legend.margin = margin(0, 0, 0, -11), 
#         legend.key.size = unit(3.4, "mm"), 
#         legend.text = element_text(size = 8))
# ggsave(filename=paste0(dirI,"VCevaluation_callable_filter_ploidy.png"), plot=p, dpi=400, width=2500/400, height=2000/400)
# 
# 
# 
# p<- ggplot(dfs_c_counts_ploidy, aes(x =factor( num_samples, levels = c("2samples", "5samples", "10samples", "15samples")), y = count, fill = factor(eval_sample, levels = c("TP_pos_FP","TP_pos_FN_allele","FP","FN_uncalled","TP")))) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = palette) +
#   scale_y_continuous(breaks = c(0, 12500, 25000)) +
#   facet_grid(ploidy ~ factor(qual, levels = c("QUAL30", "QUAL400", "QUAL1000")), scales = "free", space = "free") +
#   geom_text_repel(data = subset(dfs_c_counts_ploidy, eval_sample == "FP"), 
#                   aes(label = reorder(count, -count)), 
#                   nudge_y = 20000, 
#                   force = 0, 
#                   size = 3, 
#                   fontface = "bold", segment.color = NA)+
#   labs(x = "Filter", y = "Number of Variant Positions", fill = "") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1.1), axis.text.y = element_text(size = 7), 
#         legend.position = "bottom",
#         legend.margin = margin(-10, 0, 0, -25), 
#         legend.key.size = unit(3.4, "mm"), 
# 
#         legend.text = element_text(size = 8))
# 
# ggsave(filename=paste0(dirI,"VCevaluation_callable_filter_ploidy_counts.png"), plot=p, dpi=450, width=2000/450, height=2000/450)
# ggsave(filename= "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/IMAGES TFG/VCevaluation_callable_filter_ploidy_counts.png" , 
#        plot=p, dpi=450, width=2000/500, height=2000/500)
# 
# 
# 
# # PER REGION 
# t<-  as.data.frame(table(dfs_c$filter, dfs_c$region,  dfs_c$eval_sample))
# colnames(t)<- c("filter", "region", "eval", "counts")
# t<-  t[t$counts!=0,]
# t$region<-  factor(t$region, c("5_ETS", "18S", "ITS1", "5.8S", "ITS2",  "28S", "3_ETS"))
# 
# ggplot(t, aes(x = region, y = counts, fill = reorder(eval, counts))) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = palette)+ 
#   guides(fill = guide_legend(nrow = 1, override.aes = list(size=1))) +
#   facet_grid(.~filter)+
#   geom_text_repel(aes(label = reorder(counts, counts)),  position = position_stack(vjust = 0.5), force = .001,size = 3, fontface = "bold") +
#   labs( x = "rDNA region", y = "Number of variants", fill="")+
#   theme_bw() +
#   theme(legend.position = "bottom",    legend.text = element_text(size = 7),  
#         legend.key.size = unit(4, "mm"),      legend.margin = margin(0, 10, 0, -30),     legend.box = "horizontal")
# 
# ggsave(filename=paste0(dirI,"VCevaluation_region_prejoingen.png"), plot=p, dpi=300, width=8000/300, height=1000/300)
# 
# ############################     LET'S GET BEST PERFORMING SETS   ####################################
# dfs<- data.frame()
# 
# for (ploidy in c(50, 100)){
#   for (version in 2:4){
#     for(filter in c(2,5,10,15)){
#       d<- read.table(paste0(dirf, "VCeval_prejoin_p", ploidy, "_filtering_", filter, "sample_v", version, ".tbl"), header = T)
#       d$filter<- paste0( ploidy," ",filter, "samples ", ifelse( version==2, "QUAL30", ifelse(version==3, "QUAL400", "QUAL1000")))
#       print(paste0( ploidy," ",filter, "samples ", ifelse( version==2, "QUAL30", ifelse(version==3, "QUAL400", "QUAL1000"))))
#       dfs<- rbind(dfs, d)
#     }
#   }
# }
# 
# ######### --------------------------  Preparing for plotting  ----------------------------- #############
dfs<- data.frame()

d<- read.table(paste0(dirf, "VCeval_prejoin_p", 2, "_filtering_2sample_v2.tbl"), header = T)
d$filter<-  "2 2samples QUAL30"
d <- d[d$class =="callable",]
dfs<- rbind(dfs, d)
d<- read.table(paste0(dirf, "VCeval_prejoin_p", 5, "_filtering_2sample_v2.tbl"), header = T)
d$filter<-  "5 2samples QUAL30"
d <- d[d$class =="callable",]
dfs<- rbind(dfs, d)
d<- read.table(paste0(dirf, "VCeval_prejoin_p", 7, "_filtering_2sample_v2.tbl"), header = T)
d$filter<-  "7 2samples QUAL30"
d <- d[d$class =="callable",]
dfs<- rbind(dfs, d)
d<- read.table(paste0(dirf, "VCeval_prejoin_p", 10, "_filtering_2sample_v4.tbl"), header = T)
d$filter<-  "10 2samples QUAL1000"
d <- d[d$class =="callable",]
dfs<- rbind(dfs, d)
d<- read.table(paste0(dirf, "VCeval_prejoin_p", 20, "_filtering_2sample_v3.tbl"), header = T)
d$filter<-  "20 2samples QUAL400"
d <- d[d$class =="callable",]
dfs<- rbind(dfs, d)
d<- read.table(paste0(dirf, "VCeval_prejoin_p", 30, "_filtering_2sample_v4.tbl"), header = T)
d$filter<-  "30 2samples QUAL1000"
d <- d[d$class =="callable",]
dfs<- rbind(dfs, d)
d<- read.table(paste0(dirf, "VCeval_prejoin_p", 40, "_filtering_2sample_v4.tbl"), header = T)
d$filter<-  "40 2samples QUAL1000"
d <- d[d$class =="callable",]
dfs<- rbind(dfs, d)
d<- read.table(paste0(dirf, "VCeval_prejoin_p", 50, "_filtering_2sample_v4.tbl"), header = T)
d$filter<-  "50 2samples QUAL1000"
d <- d[d$class =="callable",]
dfs<- rbind(dfs, d)
d<- read.table(paste0(dirf, "VCeval_prejoin_p", 100, "_filtering_2sample_v4.tbl"), header = T)
d$filter<-  "100 2samples QUAL1000"
d <- d[d$class =="callable",]
dfs<- rbind(dfs, d)


levels<-c("2 2samples QUAL30",   "5 2samples QUAL30",      "7 2samples QUAL30",  
          "10 2samples QUAL1000", "20 2samples QUAL400",  "30 2samples QUAL1000",
          "40 2samples QUAL1000", "50 2samples QUAL1000",  "100 2samples QUAL1000")

# Get statistics ploidy level
dfs_c<- dfs[dfs$class =="callable",]
dfs_c_counts_ploidy <- get_counts_per_ploidy(dfs_c)

dfs_c_counts_ploidy$ploidy <- word(dfs_c_counts_ploidy$filter, 1)
dfs_c_counts_ploidy$num_samples <- factor(word(dfs_c_counts_ploidy$filter, 2))
dfs_c_counts_ploidy$qual <- factor(word(dfs_c_counts_ploidy$filter, 3))

write.table(dfs_c_counts_ploidy, paste0(dirf, "VC_filtered_2sample.tbl"),sep=",", quote=F, row.names=F, col.names=T)


print("written in ")
print(paste0(dirf, "VC_filtered_2sample.tbl"))
print("####################################")
###############
