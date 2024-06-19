#  TP_allele_FP_allele : variant postions where an allele is correctly identified but 
# 
#dir <-  "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/pre_join_gen/"
dir <- "~/Desktop/ribosomal_RNAs/Alba/03_VCsimulated/v2/"
dirf<- paste0(dir, "VC_out/Sample_level_eval/default_4.3.0.0/")
dirt <- paste0(dir, "scripts/Rscripts/")
dirI<- paste0(dir, "IMAGES/03_POSTJOIN/AnalysisFN/")
###### Libraries ######
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)
####  Palettes  #####
palette <-  c("TP" ="#75bb6f","TP_FN_allele"="#aeb852",    "TP_pos_other"="#706f4d", 
              "FN_uncalled"="#ffc60e", "FN_ref_homozygous"= "#4c89ee","FN_alt_homozygous"="#71a5be",
              "FP"= "#bd0752", "TP_pos_FP"= "#ff4d7a", "TP_pos_FP_sample"= "#ad4752",  "TP_pos_FP_allele"="#ed8c96","TP_allele_FP_allele"="#ed8c96",  "err"= "#a6a6a6") 
#
palette2 <- c( "5_ETS" = "#2f4b5b",   "18S" = "#6ad68e",  "ITS1" = "#9dbccf",
               "5.8S" = "#90c230",  "ITS2" = "#4d8195",  "28S" = "#71d064",  "3_ETS" = "#375794")

palette3 <- c( "5_ETS" = "#2f4b5b",   "18S" = "#ffe290",  "ITS1" = "#9dbccf",
               "5.8S" = "#fdae61",  "ITS2" = "#4d8195",  "28S" = "#ff6e7d",  "3_ETS" = "#375794")

##############------------------------ LOAD DATA   ----------------------------##################
prejoin<-  read.table( paste0(dirf,"VCeval_prejoin_VC.tbl"), sep="\t")
postjoin<-  read.table(  paste0(dirf,"VCeval_postjoin_filtering_v3.tbl"), sep="\t")
######### --------------------------------------  PREP DATA ---------------------------------############

#keep callable prejoin
prejoin<-  prejoin[prejoin$class == "callable", ]
prejoin$eval_sample[prejoin$eval_sample=="TP_pos_TP_allele_FP_allele"]<- "TP_allele_FP_allele"
# postjoin are already only callable

# label regions
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
                                                                                  ifelse(POS <= 22670, "3_ETS", "IGS"))))))))}
  
  d$region<- factor(d$region, levels =c("5_ETS", "18S", "ITS1", "5.8S", "ITS2",  "28S", "3_ETS") )
  
  return(d)
}
prejoin<- add_regions_chrR(prejoin)
postjoin<- add_regions_chrR(postjoin)
# 
dim(prejoin)
# 103422     25
dim(postjoin)
#[1] 109895     29

##############------------- Evaluate FNs postjoin and TP_allele_FP_allele prejoin  -------------------################

# Lets select the positions that after vc are evaluated as TP_allele_FP_allele and see what happens after join genotyping
pre_FPTP<- prejoin[prejoin$eval_sample =="TP_allele_FP_allele",]
pre_FPTP$eval_post<-  ""
# Lets see what positions are becoming FN after the joingenotyping
post_FN<- postjoin[postjoin$eval_sample =="FN_uncalled",  c("sample", "position", "ploidy", "eval_sample")]
pre_post_FN<- data.frame()

# RETRIEVE the other evaluation in both cases
for(p in c(2,5,7,10)){
  for (s in paste0("sample_", 1:100)){
    # Get what they become after join genotyping
    for (pos in unique(pre_FPTP$position[pre_FPTP$sample==s & pre_FPTP$ploidy==p])){

        pre_FPTP$eval_post[pre_FPTP$sample==s &
                             pre_FPTP$ploidy==p & 
                             pre_FPTP$position == pos]<- postjoin$eval_sample[postjoin$sample==s &
                                                                               postjoin$ploidy==p & 
                                                                               postjoin$position == pos]
    }
    # Get what they were beofore join genotyping
    for (pos in unique(post_FN$position[post_FN$sample==s & post_FN$ploidy==p])){
      pre_post_FN<- rbind(pre_post_FN, 
                          prejoin[prejoin$sample==s & prejoin$ploidy==p & prejoin$position == pos,])
      
    }
  print(paste(p,s))
  }
}


##############----------------------- PREP TABLES  ----------------------------##################

FP_ploidy<-  as.data.frame.table(table(pre_FPTP$ploidy, pre_FPTP$eval_sample, pre_FPTP$eval_post))
colnames(FP_ploidy)<- c("ploidy", "eval_pre", "eval_post", "counts")
FP_sample<-  as.data.frame.table(table(pre_FPTP$ploidy, pre_FPTP$eval_sample, pre_FPTP$eval_post, pre_FPTP$sample))
colnames(FP_sample)<- c("ploidy", "eval_pre", "eval_post", "sample", "counts")
FP_region<-  as.data.frame.table(table(pre_FPTP$ploidy, pre_FPTP$eval_sample,  pre_FPTP$eval_post, pre_FPTP$region))
colnames(FP_region)<- c("ploidy", "eval_pre", "eval_post", "region","counts")

pre_post_FN$eval_post<- "FN_uncalled"
FN_ploidy<-  as.data.frame.table(table(pre_post_FN$ploidy, pre_post_FN$eval_post, pre_post_FN$eval_sample))
colnames(FN_ploidy)<- c("ploidy", "eval_post", "eval_pre", "counts")
FN_sample<-  as.data.frame.table(table(pre_post_FN$ploidy, pre_post_FN$eval_post, pre_post_FN$eval_sample, pre_post_FN$sample))
colnames(FN_sample)<- c("ploidy", "eval_post", "eval_pre", "sample","counts")
FN_region<-  as.data.frame.table(table(pre_post_FN$ploidy, pre_post_FN$eval_post, pre_post_FN$eval_sample, pre_post_FN$region))
colnames(FN_region)<- c("ploidy", "eval_post", "eval_pre", "region","counts") 

FP_ploidy <- FP_ploidy[FP_ploidy$counts != 0, ]
FP_sample <- FP_sample[FP_sample$counts != 0, ]
FP_region <- FP_region[FP_region$counts != 0, ]
FN_ploidy <- FN_ploidy[FN_ploidy$counts != 0, ]
FN_sample <- FN_sample[FN_sample$counts != 0, ]
FN_region <- FN_region[FN_region$counts != 0, ]

########-------------------- TP_allele_FP_allele  -------------------------################
png(file=paste0(dirI,"prejoin_TP_FP_alleles_ploidy.png"), res=600, width=2500, height=2500);
ggplot(FP_ploidy, aes(x=ploidy, y=counts, fill=eval_post))+
  geom_bar(stat="identity")+
  geom_text_repel(aes(label=counts),force =0.001, position = position_stack(vjust = 0.5),
                  fontface="bold", size=3.5, color="black", max.overlaps = 50)+
  scale_fill_manual(values=palette)+
  labs(x="Ploidy", y="Number of variants", fill="Region")+
  theme_minimal()+
  theme(legend.position = "bottom", title = element_text(size = 9))
dev.off();

png(file=paste0(dirI,"prejoin_TP_FP_alleles_sample.png"), res=120, width=2000, height=1000);
ggplot(FP_sample, aes(x=factor(sample, levels = paste0("sample_", 1:100)), y=counts, fill=eval_post))+
  geom_bar(stat="identity")+
  geom_text_repel(aes(label=counts),force =0.001, position = position_stack(vjust = 0.5),
                  fontface="bold", size=3, color="black", max.overlaps = 50)+
  scale_fill_manual(values=palette)+
  facet_grid(ploidy~.)+
  labs( y="Number of variants",x="", fill="Region")+
  theme_bw()+
  theme(legend.position = "bottom", title = element_text(size = 9), axis.text.x = element_text(angle=45,  vjust = 0.5))
dev.off();

png(file=paste0(dirI,"prejoin_TP_FP_alleles_region.png"), res=300, width=1500, height=2000);
ggplot(FP_region, aes(x=region, y=counts, fill=eval_post))+
  geom_bar(stat="identity")+
  geom_text_repel(aes(label=counts),force =0.001, position = position_stack(vjust = 0.5),
                  fontface="bold", size=3, color="black", max.overlaps = 50)+
  scale_fill_manual(values=palette)+
  facet_grid(ploidy~.)+
  labs(y="Number of variants", fill="Region")+
  theme_bw()+
  theme(title = element_text(size = 9), legend.text = element_text(size=7.5,  margin = margin(r=0, unit = "pt")),
        legend.box.spacing = unit(0.2, "pt"))
dev.off();

########----------------------------------- FNs ----------------------------------################
png(file=paste0(dirI,"postjoin_FN_ploidy.png"), res=600, width=3000, height=2000);
ggplot(FN_ploidy, aes(x=ploidy, y=counts, fill=eval_pre))+
  geom_bar(stat="identity")+
  geom_text_repel(aes(label=counts),force =0.001, position = position_stack(vjust = 0.5),
                  fontface="bold", size=3.5, color="black", max.overlaps = 50)+
  scale_fill_manual(values=palette)+
  labs(x="Ploidy", y="Number of variants", fill="Region")+
  theme_minimal()+
  theme(title = element_text(size = 9),
        legend.text = element_text(size=10),legend.key.size = unit(0.5, 'cm'))
dev.off();

png(file=paste0(dirI,"postjoin_FN_sample.png"), res=120, width=2000, height=1000);
ggplot(FN_sample, aes(x=factor(sample, levels = paste0("sample_", 1:100)), y=counts, fill=eval_pre))+
  geom_bar(stat="identity")+
  geom_text_repel(aes(label=counts),force =0.001, position = position_stack(vjust = 0.5),
                  size=3, color="black", max.overlaps = 50)+
  scale_fill_manual(values=palette)+
  facet_grid(ploidy~.)+
  labs( y="Number of variants",x="", fill="Region")+
  theme_bw()+
  theme(legend.position = "bottom", title = element_text(size = 9), axis.text.x = element_text(angle=45,  vjust = 0.5))
dev.off();

png(file=paste0(dirI,"postjoin_FN_region.png"), res=300, width=1500, height=2400);
ggplot(FN_region, aes(x=region, y=counts, fill=eval_pre))+
  geom_bar(stat="identity")+
  geom_text_repel(aes(label=counts),force =0.001, position = position_stack(vjust = 0.5),
                   size=3, color="black", max.overlaps = 50)+
  scale_fill_manual(values=palette)+
  facet_grid(ploidy~.)+
  labs(y="Number of variants", fill="Region")+
  theme_bw()+
  theme(title = element_text(size = 9), legend.text = element_text(size=7.5,  margin = margin(r=0, unit = "pt")),
        legend.box.spacing = unit(0.2, "pt"))
dev.off();