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
              "FP"= "#bd0752", "TP_pos_FP"= "#ff4d7a", "TP_pos_FP_allele"="#ed8c96","TP_allele_FP_allele"="#ed8c96",  "err"= "#a6a6a6") 
#
## --------------------------------------  LOAD DATA ---------------------------------   ############
# dirf<- paste0(dir, "VC_out/Sample_level_eval/4.3.0.0/VC/")
dirf<- paste0(dir, "VC_out/Sample_level_eval/4.3.0.0/VC/")
dirI<- paste0(dir, "IMAGES/01_PREJOIN/4.3.0.0/final/2sample/")
d<- read.table(paste0(dirf, "VCeval_prejoin_p50_filtering_2sample_v4.tbl"), sep="\t")
true_variants <- read.table(paste0(dir, "scripts/Rscripts/true_variants.tbl"))

a<-d
d<-d[d$class=="callable",]

## --------------------------------------   ---------------------------------   ############
t<- as.data.frame(table(d$position, d$eval_sample))
names(t)<- c("position", "eval", "counts")
t<- t[t$counts!=0,]

t$type<- ""

for (i in 1:nrow(true_variants)){
  t$type[t$position==true_variants$position[i]]<- true_variants$type[i]
}

t$eval<- factor(t$eval, levels=c("FN_uncalled", "TP_pos_FN_allele","TP_allele_FP_allele","TP", "TP_pos_FP", "FP" ))

png(file=paste0(dirI,"Positions_eval_boxplot_ploidy50_v4.png"), res=270, width=1000, height=1000);
ggplot(t, aes(x=as.factor(eval), y=counts, fill=eval, color=eval))+
  geom_jitter( width=0.35, size=1.5, alpha=0.9)+
  geom_boxplot(alpha=0.02, color="black", width=0.9,outlier.shape = NA) +
  scale_fill_manual(values=palette) +
  scale_color_manual(values=palette) +
  geom_hline(yintercept=15, color="blue", linetype="dashed", linewidth=0.6)+
  geom_text(aes(label = "15", x="TP_pos_FN_allele", y = 18), color = "blue", size=4)+
  geom_hline(yintercept=2, color="blue", linetype="dashed", linewidth=0.6)+
  geom_text(aes(label = "2", x="TP_pos_FN_allele", y = 6), color = "blue", size=4)+
  guides(fill="none", color="none")+
  # geom_text_repel(position="stack")+
  labs(x="rDNA Region", y="Position evaluation counts", fill="Evaluation") +
  theme_minimal()+
  theme( axis.title = element_text(color="black", size=12), axis.text.x = element_text(color="black", angle=25,size=10, vjust=0.8))
dev.off();

png(file=paste0(dirI,"Positions_eval_boxplot_ploidy50_v4_labels.png"), res=470, width=2000, height=2000);
ggplot(t, aes(x=as.factor(eval), y=counts, fill=eval, color=eval))+
  geom_jitter( width=0.1, size=1.5, alpha=0.9)+
  geom_boxplot(alpha=0.02, color="black", width=0.7,outlier.shape = NA) +
  scale_fill_manual(values=palette) +
  scale_color_manual(values=palette) +
  # geom_hline(yintercept=15, color="blue", linetype="dashed", linewidth=0.6)+
  # geom_text(aes(label = "15", x="TP_pos_FN_allele", y = 18), color = "blue", size=4)+
  # geom_hline(yintercept=2, color="blue", linetype="dashed", linewidth=0.6)+
  # geom_text(aes(label = "2", x="TP_pos_FN_allele", y = 6), color = "blue", size=4)+
  guides(fill="none", color="none")+
  geom_text_repel(data=subset(t, t$counts>55 & t$eval %in% c("FP", "FN_uncalled")), aes(x=eval, y=counts,label=position), 
                  max.overlaps = 15, color="black", 
                  position="identity", size=3)+
  labs(x="rDNA Region", y="Position evaluation counts", fill="Evaluation") +
  theme_minimal()+
  theme( axis.title = element_text(color="black", size=12), axis.text.x = element_text(color="black", angle=25,size=10, vjust=0.8))
dev.off();
## -------------------------------------- min dist  ---------------------------------   ############


fn <- data.frame(position = sort(as.numeric(as.character(t$position[t$eval == "FN_uncalled" & t$counts>50 ]))))
fp <- data.frame(position = sort(as.numeric(as.character(t$position[t$eval == "FP"  & t$counts>50 ]))))

fn <- data.frame(position = sort(as.numeric(as.character(t$position[t$eval == "FN_uncalled"  ]))))
fp <- data.frame(position = sort(as.numeric(as.character(t$position[t$eval == "FP"  ]))))
tp <- data.frame(position = sort(as.numeric(as.character(t$position[t$eval == "TP"  ]))))

fn$min_len[1] <- fn$position[2] - fn$position[1]
fn$min_len[nrow(fn)] <- fn$position[nrow(fn)] - fn$position[nrow(fn)-1]

for (i in 2:(nrow(fn)-1)){
  prev<- fn$position[i-1]
  pos<- fn$position[i]
  nex <- fn$position[i+1] 
  fn$min_len[i]<- ifelse( (pos-prev) > nex-pos,   nex-pos, pos-prev )
}

fp$min_len[1] <- fp$position[2] - fp$position[1]
fp$min_len[nrow(fp)] <- fp$position[nrow(fp)] - fp$position[nrow(fp)-1]
for (i in 2:(nrow(fp)-1)){
  prev<- fp$position[i-1]
  pos<- fp$position[i]
  nex <- fp$position[i+1] 
  fp$min_len[i]<- ifelse( (pos- prev) > nex-pos, nex-pos, pos-prev )
}

tp$min_len[1] <- tp$position[2] - tp$position[1]
tp$min_len[nrow(tp)] <- tp$position[nrow(tp)] - tp$position[nrow(tp)-1]
for (i in 2:(nrow(tp)-1)){
  prev<- tp$position[i-1]
  pos<- tp$position[i]
  nex <- tp$position[i+1] 
  tp$min_len[i]<- ifelse( (pos- prev) > nex-pos, nex-pos, pos-prev )
}

fn$eval<- "FN_uncalled"
fp$eval<- "FP"
tp$eval<- "TP"

a<- rbind(fn, fp, tp)

png(file=paste0(dirI,"Positions_length_ploidy50_v4.png"), res=400, width=1000, height=3000);
ggplot(a, aes(x=eval, y=min_len))+
  geom_jitter(aes(color=eval), width=0.35, size=1.5, alpha=0.9)+
  geom_boxplot(alpha=0.02, color="black", width=0.9,outlier.shape = NA) +
  scale_fill_manual(values=palette) +
  scale_color_manual(values=palette) +
  labs(y="Distance to closest variant", x="Evaluation") +
  guides(color="none", fill="none") +
  coord_cartesian(ylim=c(20, 780))+
  theme_minimal()
dev.off();

## --------------------------------------  fp and fn superpose  ---------------------------------   ############

a<- t[t$counts>50, ]
a$position <- as.numeric(as.character(a$position))

png(file=paste0(dirI, "Positions_fpfn_ploidy50_v4.png"), res=400, width=4300, height=1500)
ggplot(a, aes(x=position, y=factor(eval, levels=c("FN_uncalled", "TP_pos_FN_allele","TP_allele_FP_allele", "TP", "TP_pos_FP", "FP")), color=eval)) +
  geom_point(size=1.5, alpha=0.9) +
  scale_color_manual(values=palette) +
  scale_x_continuous(breaks=seq(9950, max(a$position), by=1000), limits=c(min(a$position), max(a$position))) +
  labs(y="Evaluation", x="Position") +
  guides(color="none") +
  geom_text_repel(data=subset(a, eval %in% c("FP", "FN_uncalled")),
                  aes(x=position, y=eval, label=position),
                  angle=0,   direction="y",   max.overlaps=15,   color="black",  
                  position="identity",  size=4) +
  coord_cartesian(ylim=c(0, 6.3)) +
  theme_minimal() +
  theme(axis.text.y=element_text(color="black", size=10))
dev.off()


## --------------------------------------   ---------------------------------   ############

for (i in 1:nrow(d)){
  d$type[d$position==true_variants$position[i]]<- true_variants$type[i]
}
t<- as.data.frame(table(d$position, d$eval_sample, d$type))
names(t)<- c("position", "eval", "type","counts")
t<- t[t$counts!=0,]

t$eval<- factor(t$eval, levels=c("FN_uncalled", "TP_pos_FN_allele","TP_allele_FP_allele","TP", "TP_pos_FP", "FP" ))

png(file=paste0(dirI,"Positions_Type_ploidy50_v4.png"), res=300, width=1900, height=1000);
ggplot(t, aes(x=as.factor(eval), y=counts, fill=eval, color=eval))+
  geom_jitter( width=0.35, size=1.5, alpha=0.9)+
  geom_boxplot(alpha=0.02, color="black", width=0.9,outlier.shape = NA) +
  scale_fill_manual(values=palette) +
  scale_color_manual(values=palette) +
  facet_grid(.~type,  scales="free_x",space="free")+
  guides(fill="none", color="none")+
  # geom_text_repel(position="stack")+
  labs(x="rDNA Region", y="Position evaluation counts", fill="Evaluation") +
  theme_classic()+
  theme( axis.title = element_text(color="black", size=12), 
         axis.text.x = element_text(color="black", angle=25,size=10, vjust=0.6))
dev.off();

t<- as.data.frame(table( d$eval_sample, d$type))
names(t)<- c(  "eval", "type","counts")
t<- t[t$counts!=0,]

t$eval<- factor(t$eval, levels=c("TP_pos_FP","FN_uncalled","TP_pos_FN_allele","TP_allele_FP_allele",  "TP", "FP" ))

t$freq<- 0
t$freq[t$type=="SNP"]<- t$counts[t$type=="SNP"]/sum(t$counts[t$type=="SNP"])
t$freq[t$type=="indel"]<- t$counts[t$type=="indel"]/sum(t$counts[t$type=="indel"])
t$freq[t$type=="mixed"]<- t$counts[t$type=="mixed"]/sum(t$counts[t$type=="mixed"])

png(file=paste0(dirI,"Positions_Type_bar_ploidy50_v4.png"), res=300, width=1000, height=1000);
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


png(file=paste0(dirI,"Positions_Type_norm_bar_ploidy50_v4.png"), res=300, width=1500, height=1200);
ggplot(t, aes(x=as.factor(type), y=freq, fill=eval, label=counts ))+
geom_bar( stat="identity", width=0.8, size=1.5, alpha=0.9)+
geom_text_repel(position = position_stack(vjust = 0.5), force = 0.001) +
scale_fill_manual(values=palette) +
# facet_grid(.~type,  scales="free_x",space="free")+
guides( color="none")+
labs(x="rDNA Region", y="Evaluation counts per position", fill="Evaluation") +
theme_minimal()+
theme( axis.title = element_text(color="black", size=12),
       legend.margin=margin(0,-1,-15,-15),legend.key.size = unit(3.4, "mm" ),
       axis.text.x = element_text(color="black", size=10, vjust=0.8))
dev.off();