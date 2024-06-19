## --------------------------------------  CHOOSE VC SET ---------------------------------   ############
# GATK version either "4.5.0.0" or "4.3.0.0"
set="4.3.0.0"
#set=  "4.5.0.0"
## STEP of the pipeline this is for gVCF files so either "VC" or"filtering"
step="VC"
step= "filtering"
version="v4"

dir <- "~/Desktop/ribosomal_RNAs/Alba/03_VCsimulated/v2/"

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)


##  Palette  
palette <-  c("TP" ="#75bb6f","TP_pos_FN_allele"="#aeb852", "TP_pos_other"="#706f4d", "FN_uncalled"="#ffc60e", "FN_ref_homozygous"= "#4c89ee","FN_alt_homozygous"="#71a5be",
              "FP"= "#bd0752",  "TP_pos_FP"= "#ff4d7a","TP_pos_FP_sample"= "#ad4752", "TP_pos_FP_allele"="#ed8c96","TP_allele_FP_allele"="#ed8c96",  "err"= "#a6a6a6") 
#
## --------------------------------------  LOAD DATA ---------------------------------   ############
# dirf<- paste0(dir, "VC_out/Sample_level_eval/4.3.0.0/VC/")
dirf<- paste0(dir, "VC_out/Sample_level_eval/4.3.0.0/VC/outdated/")
d<- read.table(paste0(dirf, "VCeval_postjoing_p7_filtering_v4.tbl"), sep="\t")
true_variants <- read.table(paste0(dir, "scripts/Rscripts/true_variants.tbl"))

a<-d
d<-d[d$class=="callable",]
dirI<- paste0(dir, "IMAGES/02_POSTJOIN/4.3.0.0/joingen/")
## --------------------------------------   ---------------------------------   ############
t<- as.data.frame(table(d$position, d$eval_sample))
names(t)<- c("position", "eval", "counts")
t<- t[t$counts!=0,]

t$type<- ""

for (i in 1:nrow(true_variants)){
  t$type[t$position==true_variants$position[i]]<- true_variants$type[i]
}

t$eval<- factor(t$eval, levels = c("TP_pos_FP","TP_allele_FP_allele",  "FP",  "TP_pos_FP_sample","TP_pos_other",
                                   "FN_alt_homozygous", "FN_ref_homozygous",
                                   "TP_FN_allele", "FN_uncalled","TP" ))

png(file=paste0(dirI,"Positions_eval_boxplot_ploidy7.png"), res=270, width=1300, height=1000);
ggplot(t, aes(x=as.factor(eval), y=counts, fill=eval, color=eval))+
  geom_jitter( width=0.35, size=1.5, alpha=0.9)+
  geom_boxplot(alpha=0.02, color="black", width=0.9,outlier.shape = NA) +
  scale_fill_manual(values=palette) +
  scale_color_manual(values=palette) +
  geom_hline(yintercept=15, color="blue", linetype="dashed", linewidth=0.6)+
  geom_text(aes(label = "15", x="FN_uncalled", y = 18), color = "blue", size=4)+
  geom_hline(yintercept=2, color="blue", linetype="dashed", linewidth=0.6)+
  geom_text(aes(label = "2", x="FN_uncalled", y = 6), color = "blue", size=4)+
  guides(fill="none", color="none")+
  # geom_text_repel(position="stack")+
  labs(x="Evaluation", y="Position evaluation counts", fill="Evaluation") +
  theme_minimal()+
  theme( axis.title = element_text(color="black", size=12), 
         axis.text.x = element_text(color="black", angle=25, size=10, vjust=0.8, hjust=0.7))
dev.off();


for (i in 1:nrow(d)){
  d$type[d$position==true_variants$position[i]]<- true_variants$type[i]
}
t<- as.data.frame(table(d$position, d$eval_sample, d$type))
names(t)<- c("position", "eval", "type","counts")
t<- t[t$counts!=0,]

t$eval<- factor(t$eval,levels = c("TP_pos_FP","TP_allele_FP_allele",  "FP",  "TP_pos_FP_sample","TP_pos_other",
                                  "FN_alt_homozygous", "FN_ref_homozygous",
                                  "TP_FN_allele", "FN_uncalled","TP" ))


png(file=paste0(dirI,"Positions_Type_ploidy7_v4.png"), res=300, width=1900, height=1000);
ggplot(t, aes(x=as.factor(eval), y=counts, fill=eval, color=eval))+
  geom_jitter( width=0.35, size=1.5, alpha=0.9)+
  geom_boxplot(alpha=0.02, color="black", width=0.9,outlier.shape = NA) +
  scale_fill_manual(values=palette) +
  scale_color_manual(values=palette) +
  facet_grid(.~type,  scales="free_x",space="free")+
  guides(fill="none", color="none")+
  # geom_text_repel(position="stack")+
  labs(x="rDNA Region", y="Position evaluation counts", fill="Evaluation") +
  theme_minimal()+
  theme( axis.title = element_text(color="black", size=12), axis.text.x = element_text(color="black", angle=25,size=10, vjust=0.8))
dev.off();

t<- as.data.frame(table( d$eval_sample, d$type))
names(t)<- c(  "eval", "type","counts")
t<- t[t$counts!=0,]

t$eval<- factor(t$eval,levels = c("TP_pos_FP","TP_allele_FP_allele",  "FP",  "TP_pos_FP_sample","TP_pos_other",
                                  "FN_alt_homozygous", "FN_ref_homozygous",
                                  "TP_FN_allele", "FN_uncalled","TP" ))

t$freq<- 0
t$freq[t$type=="SNP"]<- t$counts[t$type=="SNP"]/sum(t$counts[t$type=="SNP"])
t$freq[t$type=="indel"]<- t$counts[t$type=="indel"]/sum(t$counts[t$type=="indel"])
t$freq[t$type=="mixed"]<- t$counts[t$type=="mixed"]/sum(t$counts[t$type=="mixed"])

png(file=paste0(dirI,"Positions_Type_bar_ploidy7_v4.png"), res=300, width=1000, height=1000);
ggplot(t, aes(x=as.factor(type), y=counts, fill=eval, label=counts ))+
  geom_bar( stat="identity",   size=1.5, alpha=0.9)+
  geom_text_repel(position = position_stack(vjust = 0.5), force = 0.001 ) +
  scale_fill_manual(values=palette) +
  # facet_grid(.~type,  scales="free_x",space="free")+
  guides(fill="none", color="none")+
  labs(x="rDNA Region", y="Position evaluation counts", fill="Evaluation") +
  theme_minimal()+
  theme( axis.title = element_text(color="black", size=12), axis.text.x = element_text(color="black", angle=25,size=10, vjust=0.8))
dev.off();


png(file=paste0(dirI,"Positions_Type_norm_bar_ploidy7_v4.png"), res=300, width=1000, height=1000);
ggplot(t, aes(x=as.factor(type), y=freq, fill=eval, label=counts ))+
  geom_bar( stat="identity", width=0.8, size=1.5, alpha=0.9)+
  geom_text_repel(position = position_stack(vjust = 0.5), force = 0.001) +
  scale_fill_manual(values=palette) +
  # facet_grid(.~type,  scales="free_x",space="free")+
  guides(fill="none", color="none")+
  labs(x="rDNA Region", y="Position evaluation counts", fill="Evaluation") +
  theme_minimal()+
  theme( axis.title = element_text(color="black", size=12), axis.text.x = element_text(color="black", size=10, vjust=0.8))
dev.off();