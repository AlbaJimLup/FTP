library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)

##  Palette  
palette <-  c("TP" ="#75bb6f","TP_pos_FN_allele"="#aeb852", "TP_pos_other"="#706f4d", "FN_uncalled"="#ffc60e", 
              "FN_ref_homozygous"= "#4c89ee","FN_alt_homozygous"="#71a5be",
              "FP"= "#bd0752", "TP_pos_FP"= "#ff4d7a", "TP_pos_FP_allele"="#ed8c96","TP_allele_FP_allele"="#ed8c96",  "err"= "#a6a6a6") 
### ------------------------------------ SET PARAMETERS --------------------------------- ############
#### CHOOSE VC SET ######
args <- commandArgs(trailingOnly = TRUE)
#GATK version either "4.5.0.0" or "4.3.0.0"
set=args[1] 
# STEP of the pipeline this is for gVCF files so either "VC" or"filtering"
step=args[2]
ploidy=args[3]
version=args[4]

step="filtering"
ploidy=7
version="v2"
# version =paste0("2sample_", version)

print(paste("Doing VC evaluation for", set, step, ploidy, version))
print("#############################################################################")

dir <- "~/Desktop/ribosomal_RNAs/Alba/08_validation/"
dirf<- paste0(dir, "VC_out/")
## --------------------------------------  LOAD DATA --------------------------------- ############
#Get true_variants 
# true_variants <- read.table(paste0(dir, "Rscripts/aligned_10.vcf"))
# names(true_variants)<- c("chr", "position", "id", "ref", "alt")
# true_variants$position<- true_variants$position+9338
# bed <- read.table(paste0("~/Desktop/ribosomal_RNAs/Alba/data/47S_pre-rRNA.blacklisted_repetitive_seq.bed"))
# # 
# d<- read.table(paste0(dirf,"VCeval_postjoing_p7_joingen.tbl"), sep="\t")
# dirI<- paste0(dir, "IMAGES/joingen7/")
d<- read.table(paste0(dirf,"VCeval_postjoing_p7_filtering_v2.tbl"), sep="\t")
dirI<- paste0(dir, "IMAGES/joingen7_v2/")

d<-d[d$class=="callable",]

## -------------------------------------- PLOT RESULTS  ---------------------------------   ############
t<- as.data.frame(table(d$position, d$eval_sample))
names(t)<- c("position", "eval", "counts")
t<- t[t$counts!=0,]

t$type<- ""

t$eval<- factor(t$eval, levels=c("FN_uncalled", "FN_alt_homozygous", "TP_pos_FN_allele", "TP", "TP_pos_FP", "FP" ))

png(file=paste0(dirI,"00_Joingen_Positions_eval_boxplot_ploidy7.png"), res=270, width=1100, height=1000);
ggplot(t, aes(x=as.factor(eval), y=counts, fill=eval, color=eval))+
  geom_jitter( width=0.35, size=1.5, alpha=0.9)+
  geom_boxplot(alpha=0.02, color="black", width=0.9,outlier.shape = NA) +
  scale_fill_manual(values=palette) +
  scale_color_manual(values=palette) +
  geom_hline(yintercept=2, color="blue", linetype="dashed", linewidth=0.6)+
  geom_text(aes(label = "2", x="TP_pos_FN_allele", y = 2.15), color = "blue", size=4)+
  guides(fill="none", color="none")+
  # geom_text_repel(position="stack")+
  labs(x="rDNA Region", y="Position evaluation counts", fill="Evaluation") +
  theme_minimal()+
  theme( axis.title = element_text(color="black", size=12), axis.text.x = element_text(color="black", angle=25,size=10, vjust=0.8))
dev.off();

## -------------------------------------- PLOT RESULTS  ---------------------------------   ############
# Select only callable
dfs_c <- d
#
get_counts_per_ploidy <- function(dfs){
  
  counts<- as.data.frame(table(dfs$eval_sample))
  names(counts)<- c("eval_sample", "count")
  
  ploydies<- unique(counts$ploidy)
  # Get as frequencies as well
  ## remove 0 values so they dont appear in plot:
  counts <- counts[counts$count != 0,]
  
  counts$freq <- rep(0, nrow(counts))
  
  sum_val <- sum(counts$count )
  for (e in counts$eval_sample ) {
    counts$freq[ counts$eval_sample == e] <- round(counts$count[ counts$eval == e] / sum_val, 8)
  }
  
  # Order the evaluation so they are descending in plot counts
  counts<-   counts %>%  arrange(count)
  
  return(counts)
}
get_counts_per_sample <- function(dfs){
  
  df_counts<- as.data.frame(table( dfs$sample, dfs$eval_sample))
  names(df_counts)<- c( "sample", "eval_sample", "count")
  # Get frequencies
  eval<- unique(df_counts$eval_sample)
  samples<- c("SRR1997411", "SRR3189741", "SRR3189742", "SRR3189743")
  df_counts$freq <- 0
  
  for (sample in samples){
    sum_val <- sum(df_counts$count[  df_counts$sample == sample])
    for (e in eval) {
      freq <- round(df_counts$count[  df_counts$sample == sample & df_counts$eval_sample == e] / sum_val, 8)
      df_counts$freq[ df_counts$sample == sample & df_counts$eval_sample == e] <- freq
      }
    }
  
  # Order the evaluation types in descending order of counts
  df_counts <- df_counts %>% arrange(desc(count))
  
  # Remove rows with count = 0
  df_counts <- df_counts[df_counts$count != 0,]
  
  return(df_counts)
}
#
########## ---------------------------- PLOT DIST VARIANTS evalutions-----------------------------------#################  
 
positions_t<- as.data.frame.table(table(dfs_c$position, dfs_c$eval_sample, dfs_c$ploidy, dfs_c$region))
names(positions_t)<- c("pos", "eval",  "ploidy","region","counts")
positions_t<- positions_t[positions_t$counts!=0,]
positions_t$region<- factor(positions_t$region, levels =c("5_ETS", "18S", "ITS1", "5.8S", "ITS2",  "28S", "3_ETS") )
write.table(positions_t, paste0(dirf, "Positions_per_evaluation.tbl"))

get_plot<- function(positions, title, size){
  png(file=paste0(dirI, title ,"_Joingen_evaluation_callable_ploidy7.png"), res=300, width=4500, height=1200);
  print(ggplot(positions, aes(x=pos,y=counts, fill=eval))+
          geom_bar(stat="identity")+
          # coord_cartesian(xlim = c(min(positions_t$pos)+12,max(positions_t$pos)-12))+
          facet_grid(ploidy~region, scales="free_x",  space="free_x")+
          theme_bw()+
          labs(x="variant positions", y="variant counts", fill="")+
          guides(fill = guide_legend(nrow = 1))+
          scale_fill_manual(values = palette)+
          theme(legend.position = "bottom", legend.direction = "horizontal",
                legend.text = element_text(size = 8), legend.key.size = unit(3, "mm" ),
                legend.margin=margin(-5,-10,0,0),
                axis.text.x = element_text(angle=45, size=size), strip.text = element_text(size = 7, face = "bold"),  
                axis.title = element_text(size=8),
                strip.background = element_rect(fill = "lightgray", color = "black"), panel.spacing = unit(1, "pt")))
  dev.off();
  return("Done!")
}

get_plot(positions_t, "Positions", 3.6)
positions_t$pos<- as.factor(positions_t$pos)

#Select those positions that are always FP or FN
positions_FPFN <- data.frame()
for (p in ploidy){
  subset<- positions_t[positions_t$ploidy == p, ]
  for (pos in subset$pos){
    eval <-  subset$eval[subset$pos == pos]
    if (length(eval) ==1){
      positions_FPFN<- rbind(positions_FPFN, subset[subset$pos == pos, ])
    }        
  }
}
get_plot(positions_FPFN, "Positions_FP_FN", 5)
write.table(positions_FPFN, paste0(dirf, "FPFN_Positions_per_evaluation_Joingen.tbl"))

##### Get statistics ploidy level ######
dfs_c_counts_ploidy <- get_counts_per_ploidy(dfs_c)
dfs_c_counts_ploidy$eval_sample<- factor(dfs_c_counts_ploidy$eval_sample, levels = c("TP_pos_FP","TP_allele_FP_allele",  "FP",  "TP_pos_FP_sample","TP_pos_other", 
                                                                                     "FN_alt_homozygous", "FN_ref_homozygous", "FN_uncalled",  "TP_pos_FN_allele","TP" ))
dfs_c_counts_ploidy$position<- dfs_c_counts_ploidy$freq

png(file=paste0(dirI, "00_Joingen_evaluation_callable_ploidy7.png"), res=500, width=1500, height=1500);
ggplot(dfs_c_counts_ploidy, aes(x =  ploidy,   y = freq, fill = eval_sample)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim=c(0.035, 0.97))+
  scale_fill_manual(values = palette)+ 
  geom_text_repel(aes(label = reorder(count, count)),  position = position_stack(vjust = 0.5), force = .005,size = 3.2, fontface = "bold") +
  labs(x = "Ploidy", y = "Frequency", fill="")+
  theme_classic()+
  theme(legend.margin=margin(0,0,0,-11), legend.key.size = unit(3.4, "mm" ), legend.text = element_text(size=8))
dev.off()

############## MOST IMPORTANT PLOT ###########
# Get statistics sample level
dfs_c_counts_sample <- get_counts_per_sample(dfs_c)

dfs_c_counts_ploidy$eval_sample<- factor(dfs_c_counts_ploidy$eval_sample, levels = c("TP_pos_FP","TP_allele_FP_allele",  "FP",  "TP_pos_FP_sample","TP_pos_other",
                                                                                     "FN_alt_homozygous", "FN_ref_homozygous", "FN_uncalled",  "TP_FN_allele", "TP" ))


png(file=paste0(dirI, "00_Joingen_evaluation_callable_samples_ploidy7.png"), res=300, width=1300, height=1000);
ggplot(dfs_c_counts_sample, aes(x = factor(sample, levels =c("SRR1997411",  "SRR3189741",  "SRR3189742",  "SRR3189743")), y = count, fill = reorder(eval_sample, count))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette)+ 
  labs( x = "", y = "Counts", fill="")+
  geom_text_repel(aes(label = reorder(count, count)),  position = position_stack(vjust = 0.5), force = 0.0001,size =4, fontface = "bold") +
  # guides(fill = guide_legend(nrow = 1, override.aes = list(size=3)))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=15, size = 8,color="black", margin = margin(t = 13, r = 0, b = 0, l = 0)), 
        legend.margin=margin(-40,-2,0,-10),
        legend.text = element_text(size=9),
        legend.key.size = unit(3, "mm" ),
        legend.position = "right")
dev.off()
