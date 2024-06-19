library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)
#####
set="4.3.0.0"
# step="VC"
step="filtering"
# version="v3"
# version="v2"
# version="v1"
ploidies=c(2,5,7,10,50,100)
## --------------------------------------  LOAD DATA ---------------------------------   ############
#dir <-  "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/pre_join_gen/"
dir <- "~/Desktop/ribosomal_RNAs/Alba/03_VCsimulated/v2/"
dirf<- paste0(dir, "VC_out/Sample_level_eval/", set,"/VC/")

if(step=="VC"){
  dirI<- paste0(dir, "IMAGES/01_PREJOIN/", set, "/VC/")
  dfs <- read.table( paste0(dirf,"VCeval_postjoin_VC.tbl"), sep="\t")
  
}else{
  dirI<- paste0(dir, "IMAGES/01_PREJOIN/", set, "/filtering_", version, "/")
  dfs <- read.table( paste0(dirf,"VCeval_postjoin_filtering_", version, ".tbl"), sep="\t")
}
##  Palette  
palette <-  c("TP" ="#75bb6f","TP_FN_allele"="#aeb852",    "TP_pos_other"="#706f4d", 
              "FN_uncalled"="#ffc60e", "FN_ref_homozygous"= "#4c89ee","FN_alt_homozygous"="#71a5be",
              "FP"= "#bd0752", "TP_pos_FP"= "#ff4d7a", "TP_pos_FP_sample"= "#ad4752", "TP_pos_FP_allele"="#ed8c96","TP_allele_FP_allele"="#ed8c96",  "err"= "#a6a6a6") 
#
######### --------------------------  Preparing for plotting  ----------------------------- #############

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
  
  return(df_counts)
}

# dfs$eval_sample[dfs$eval_sample =="TP_FP_allele"] <-  "TP_pos_FP_allele"
dfs$eval_sample[dfs$eval_sample =="TP_pos_FN_allele"] <- "TP_FN_allele"
dfs$eval_sample[dfs$eval_sample =="FP_TP_pos"] <- "TP_pos_FP"
dfs$eval_sample[dfs$eval_sample =="TP_pos_TP_allele_FP_allele"] <- "TP_allele_FP_allele"

# Get statistics ploidy level
dfs_counts_ploidy <- get_counts_per_ploidy(dfs)
# Get statistics sample level
dfs_counts_sample <- get_counts_per_sample(dfs)

########## ---------------------------- Select only callable ----------------------------------------------------
dfs_c <- dfs[dfs$class =="callable",]
################# PLOT DIST VARIANTS evalutions ##################
positions_t<- as.data.frame.table(table(dfs_c$position, dfs_c$eval_sample, dfs_c$ploidy, dfs_c$region))
names(positions_t)<- c("pos", "eval",  "ploidy","region","counts")
positions_t<- positions_t[positions_t$counts!=0,]
positions_t$region<- factor(positions_t$region, levels =c("5_ETS", "18S", "ITS1", "5.8S", "ITS2",  "28S", "3_ETS") )
write.table(positions_t, paste0(dirf, "Positions_per_evaluation.tbl"))

get_plot<- function(positions, title, size){
  png(file=paste0(dirI, title ,"_VCevaluation_callable_ploidy.png"), res=300, width=4500, height=1200);
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
for (p in ploidies){
  subset<- positions_t[positions_t$ploidy == p, ]
  for (pos in subset$pos){
    eval <-  subset$eval[subset$pos == pos]
    if (length(eval) ==1){
      positions_FPFN<- rbind(positions_FPFN, subset[subset$pos == pos, ])
    }        
  }
}
get_plot(positions_FPFN, "Positions_FP_FN", 5)
write.table(positions_FPFN, paste0(dirf, "FPFN_Positions_per_evaluation.tbl"))

##### Get statistics ploidy level ######
dfs_c_counts_ploidy <- get_counts_per_ploidy(dfs_c)
dfs_c_counts_ploidy$eval_sample<- factor(dfs_c_counts_ploidy$eval_sample, levels = c("TP_pos_FP","TP_allele_FP_allele",  "FP",  "TP_pos_FP_sample","TP_pos_other",  "FN_alt_homozygous", "FN_ref_homozygous",  "TP_FN_allele", "FN_uncalled","TP" ))
dfs_c_counts_ploidy$position<- dfs_c_counts_ploidy$freq


if (step=="filtering"){
  for (ploidy in unique(dfs_c_counts_ploidy$ploidy)){
    for (freq in dfs_c_counts_ploidy$freq[!dfs_c_counts_ploidy$eval_sample %in% c("TP", "FP","FN_uncalled") & dfs_c_counts_ploidy$ploidy==ploidy]){
      dfs_c_counts_ploidy$position[ dfs_c_counts_ploidy$freq==freq & dfs_c_counts_ploidy$ploidy==ploidy] <- sum(freq, dfs_c_counts_ploidy$freq[dfs_c_counts_ploidy$eval_sample=="TP"& dfs_c_counts_ploidy$ploidy==ploidy],
                                                                                                                dfs_c_counts_ploidy$freq[dfs_c_counts_ploidy$eval_sample=="FN_uncalled"& dfs_c_counts_ploidy$ploidy==ploidy])
    }
  }
  position=dfs_c_counts_ploidy$position   
  plot1 <- ggplot(dfs_c_counts_ploidy, aes(x = factor(ploidy, levels = ploidies), y = freq, fill = eval_sample )) +
    geom_bar(stat = "identity") +
    coord_cartesian(ylim=c(0.035, 1))+
    scale_fill_manual(values = palette)+
    geom_text_repel(data=subset(dfs_c_counts_ploidy,as.numeric(dfs_c_counts_ploidy$count)>400), aes(label = reorder(count, count)),  
                    alpha=0.8, position = position_stack(vjust = 0.6), force = .001, size = 3, fontface = "bold") +
    geom_text_repel(data=subset(dfs_c_counts_ploidy, as.numeric(dfs_c_counts_ploidy$count)<400 & as.numeric(dfs_c_counts_ploidy$count)>5),  aes(label = reorder(count, count), y =position),vjust = 0.9,hjust=2, 
                    force = .0015,  arrow=NULL, size = 1.8, color="black" , min.segment.length = 1, max.overlaps = 12) +
    labs(x = "Ploidy", y = "Frequency", fill="")+
    theme_classic()+
    theme(legend.margin=margin(0,0,0,-11), legend.key.size = unit(3.4, "mm" ), legend.text = element_text(size=8))
}else{
  for (ploidy in unique(dfs_c_counts_ploidy$ploidy)){
    for (freq in dfs_c_counts_ploidy$freq[!dfs_c_counts_ploidy$eval_sample %in% c("TP", "FP", "TP_pos_FP_sample","FN_uncalled") & dfs_c_counts_ploidy$ploidy==ploidy]){
      dfs_c_counts_ploidy$position[ dfs_c_counts_ploidy$freq==freq & dfs_c_counts_ploidy$ploidy==ploidy] <- sum(freq, dfs_c_counts_ploidy$freq[dfs_c_counts_ploidy$eval_sample=="TP"& dfs_c_counts_ploidy$ploidy==ploidy],
                                                                                                                dfs_c_counts_ploidy$freq[dfs_c_counts_ploidy$eval_sample=="FN_uncalled"& dfs_c_counts_ploidy$ploidy==ploidy])
    }
  }
  position=dfs_c_counts_ploidy$position 
  plot1 <- ggplot(dfs_c_counts_ploidy, aes(x = factor(ploidy, levels = ploidies), y = freq, fill = eval_sample )) +
    geom_bar(stat = "identity") +
    coord_cartesian(ylim=c(0.035, 0.97))+
    scale_fill_manual(values = palette)+
    geom_text_repel(data=subset(dfs_c_counts_ploidy, as.numeric(dfs_c_counts_ploidy$count)>1000), aes(label = reorder(count, count)),  
                    alpha=0.8, position = position_stack(vjust = 0.5), force = .005, size = 3, fontface = "bold") +
    geom_text_repel(data=subset(dfs_c_counts_ploidy, as.numeric(dfs_c_counts_ploidy$count)<1000& as.numeric(dfs_c_counts_ploidy$count)>5),  aes(label = reorder(count, count), y =position),vjust = 0.5,hjust=1.5, 
                    force = .0015,  direction="x", arrow=NULL, size = 2, color="black" , min.segment.length = 1) +
    labs(x = "Ploidy", y = "Frequency", fill="")+
    theme_classic()+
    theme(legend.margin=margin(0,0,0,-11), legend.key.size = unit(3.4, "mm" ), legend.text = element_text(size=8))
}
ggsave(filename = paste0(dirI, "00_VCevaluation_callable_ploidy_VC.png"), plot = plot1, dpi = 500, width = 2500/500, height = 1500/500)


png(file=paste0(dirI, "00_VCevaluation_callable_ploidy.png"), res=500, width=2500, height=1500);
ggplot(dfs_c_counts_ploidy, aes(x = factor(ploidy, levels = ploidies), y = freq, fill = eval_sample)) +
  
  geom_bar(stat = "identity") +
  coord_cartesian(ylim=c(0.035, 0.97))+
  scale_fill_manual(values = palette)+ 
  geom_text_repel(aes(label = reorder(count, count)),  position = position_stack(vjust = 0.5), force = .005,size = 3.2, fontface = "bold") +
  labs(x = "Ploidy", y = "Frequency", fill="")+
  theme_classic()+
  theme(legend.margin=margin(0,0,0,-11), legend.key.size = unit(3.4, "mm" ), legend.text = element_text(size=8))
dev.off()

png(file=paste0(dirI, "VCevaluation_callable_ploidy_preVC.png"), res=500, width=3000, height=1600);
ggplot(dfs_c_counts_ploidy, aes(x = factor(ploidy, levels = ploidies), y = freq, fill = eval_sample)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim=c(0.035, 0.97))+
  scale_fill_manual(values = palette)+ 
  geom_text_repel(aes(label = reorder(count, count)),  position = position_stack(vjust = 0.5), force = .005,size = 3, fontface = "bold") +
  labs(x = "Ploidy", y = "Frequency", fill="")+
  theme_classic()+
  theme(legend.margin=margin(0,0,0,-11), legend.key.size = unit(3.4, "mm" ), legend.text = element_text(size=8))
dev.off()

############## MOST IMPORTANT PLOT ###########
# Get statistics sample level
dfs_c_counts_sample <- get_counts_per_sample(dfs_c)

dfs_c_counts_ploidy$eval_sample<- factor(dfs_c_counts_ploidy$eval_sample, levels = c("TP_pos_FP","TP_allele_FP_allele",  "FP",  "TP_pos_FP_sample","TP_pos_other",
                                                                                     "FN_alt_homozygous", "FN_ref_homozygous",
                                                                                     "TP_FN_allele", "FN_uncalled","TP" ))

write.table(dfs_c_counts_sample, paste0(dirf, "callable_eval_counts_sample_evaluation.tbl"))

png(file=paste0(dirI, "VCevaluation_callable_samples.png"), res=200, width=2300, height=2000);
ggplot(dfs_c_counts_sample, aes(x = factor(sample, levels = paste0("sample_", 1:100)), y = count, fill = reorder(eval_sample, count))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette)+ 
  labs(title = "Sample level evaluation callable positions after join genotyping", x = "", y = "Counts", fill="")+
  guides(fill = guide_legend(nrow = 1, override.aes = list(size=3)))+
  facet_grid(ploidy~.)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 65,   size = 8, margin = margin(t = 13, r = 0, b = 0, l = 0)), 
        legend.margin=margin(-40,-10,0,0),
        legend.text = element_text(size=9),
        legend.key.size = unit(3, "mm" ),
        legend.position = "bottom")
dev.off()

######### ---------------------------- OTHERS ------------------------------- #############
# whole set c and b 
png(file=paste0(dirI, "VCevaluation_ploidy.png"), res=400, width=2000, height=1500);
ggplot(dfs_counts_ploidy, aes(x = factor(ploidy, levels = ploidies), y = freq, fill = reorder(eval_sample, count))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette)+ 
  geom_text_repel(aes(label = reorder(count, count)),  position = position_stack(vjust = 0.5), force = .003,size = 3, fontface = "bold") +
  labs( x = "Ploidy", y = "Frequency", fill="")+
  theme_minimal()+
  theme( legend.text = element_text(size=7),  legend.key.size = unit(4, "mm" ),  legend.margin=margin(0,10,0,-10))
dev.off()

png(file=paste0(dirI, "VCevaluation_samples.png"), res=200, width=2000, height=1500);
ggplot(dfs_counts_sample, aes(x = factor(sample, levels = paste0("sample_", 1:100)), y = count, fill = reorder(eval_sample, count))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette)+ 
  labs( x = "Samples", y = "Counts", fill="")+
  guides(fill = guide_legend(nrow = 1, override.aes = list(size=3)))+
  facet_grid(ploidy~.)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 65, vjust = 0.5, hjust=0.3, size = 8),
        legend.text = element_text(size=9),legend.key.size = unit(4, "mm" ),  legend.margin=margin(-5,0,0,0), legend.position = "bottom")
dev.off()

## CALLABLE VS BLACKLISTED 
dfs_b <- dfs[!dfs$class =="callable",]
# Get statistics ploidy level
dfs_b_counts_ploidy <- get_counts_per_ploidy(dfs_b)
# Get statistics sample level
dfs_b_counts_sample <- get_counts_per_sample(dfs_b)

dfs_c_counts_ploidy$region <- "callable"
dfs_b_counts_ploidy$region <- "blacklisted"


dfs_b_counts_ploidy$position<- dfs_b_counts_ploidy$freq

dfs_bc_ploidy<- rbind(dfs_c_counts_ploidy, dfs_b_counts_ploidy)


png(file=paste0(dirI, "VCpopulation_b_vs_c.png"), res=200, width=1300, height=1200);
ggplot(dfs_bc_ploidy, aes(x = factor(ploidy, levels = ploidies), y = freq, fill = reorder(eval_sample, count))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette)+ 
  geom_text_repel(aes(label = reorder(count, count)),  position = position_stack(vjust = 0.5), force = .0005, size = 3, fontface = "bold") +
  labs(title ="", x = "Ploidy", y = "Frequency", fill="")+
  facet_grid(.~factor(region, levels = c("callable", "blacklisted")))+
  theme_classic()+
  theme(axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        legend.text = element_text(size=9), legend.position = "bottom", 
        legend.key.size = unit(3, "mm" ),
        legend.margin=margin(-10,20,0,-10))
dev.off()

# dfs_c_counts_sample$region <- "callable"
# dfs_b_counts_sample$region <- "blacklisted"
# 
# dfs_bc_sample<- rbind(dfs_c_counts_sample, dfs_b_counts_sample)
png(file=paste0(dirI, "VCevaluation_blacklisted_samples.png"), res=200, width=2300, height=1300);
ggplot(dfs_b_counts_sample, aes(x = factor(sample, levels = paste0("sample_", 1:100)), y = count, fill = reorder(eval_sample, count))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette)+ 
  labs(title = "Sample level evaluation blacklisted positions", x = "Samples", y = "Counts", fill="")+
  facet_grid(ploidy~.)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8), 
        legend.text = element_text(size=9),legend.key.size = unit(4, "mm" ),  legend.margin=margin(-5,0,0,0), legend.position = "bottom")
dev.off()




