#  TP_allele_FP_allele : variant postions where an allele is correctly identified but 
# 
#dir <-  "/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/
dir <- "~/Desktop/ribosomal_RNAs/Alba/"
#
dirf<- paste0(dir, "03_VCsimulated/v2/VC_out/Sample_level_eval/4.3.0.0/VC/")
dirt <- paste0(dir, "03_VCsimulated/v2/scripts/Rscripts/")
dirI<- paste0(dir, "03_VCsimulated/v2/IMAGES/03_POSTJOIN/AnalysisGC/")
# 
dirBED<- paste0(dir, "data/hs_GC/")
###### Libraries ######
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggrepel)

####  Palettes  #####
palette <-  c("#9ecccf", "#4495b4","#2f4b5b", "black" )
palette2<- c( "#2f4b5b","#6ad68e", "#9dbccf", "#90c230", "#4d8195",  "#71d064",  "#375794")
palette3 <-  c("TP" ="#75bb6f","TP_FN_allele"="#aeb852","TP_pos_other"="#706f4d", "FN_uncalled"="#ffc60e", "FN_ref_homozygous"= "#4c89ee","FN_alt_homozygous"="#71a5be",
              "FP"= "#bd0752", "TP_pos_FP"= "#ff4d7a", "TP_pos_FP_allele"="#ed8c96","TP_allele_FP_allele"="#ed8c96",  "err"= "#a6a6a6", 
             "20"= "#9ecccf", "50"="#4495b4", "70"="#2f4b5b", "100"="black" ) 
#
####### Functions ##########

add_regions_chrR <- function(d){
  d$region<- ""
  # <= only in IGS at start but all other should be > only but then variants are 0 bed based
  #so variant at the end of region be annotated to be in wrong region so we should use >=
  for (POS in unique(d$middle)){
    d$region[d$middle==POS] <- ifelse(POS <= 9338, "IGS1",
                                        ifelse(POS <= 12995, "5_ETS",
                                               ifelse(POS <= 14864, "18S",
                                                      ifelse(POS <= 15934, "ITS1",
                                                             ifelse(POS <= 16091, "5.8S",
                                                                    ifelse(POS <= 17258, "ITS2",
                                                                           ifelse(POS <= 22309, "28S",
                                                                                  ifelse(POS <= 22670, "3_ETS", "IGS2"))))))))
  }
  return(d$region)
}

prep_GC_bed<- function(ws, window, step){
  names(ws)<- c("chr", "start", "end", "gc")
  # 1st position of the chrR we are accointing for is the 1st position of 5' ETS  
  ws$start<-as.numeric(ws$start)
  ws$end<- as.numeric(ws$end)
  ws$middle<- ws$start + (ws$end - ws$start) / 2 # get middle point
  ws$region<- add_regions_chrR(ws) 
  
  ws$window<- window
  ws$step<- step 
  
  # ws<-  ws[ws$gc != 0.00000, ]
  
  return(ws)
}

plot_gc<- function(ws, name, width, height, by){
  
  png(file=paste0(dirI, name ,"_GC.png"), res=300, width=width, height=height);
  print(ggplot(ws, aes(x = middle, y = gc, color=window)) +
    # geom_segment(aes(x = start, xend = end, y = gc, yend = gc )) +   # segments of gc
    geom_line() +   # join them at the middle
    scale_color_manual(values= palette )+
    labs(x = "Position", y = "GC Content") +  
    facet_grid(.~region, scales="free_x",  space="free_x")+
    scale_x_continuous(breaks = seq(min(ws$middle), max(ws$middle), by = by))+
    theme_bw()  +
    theme(  axis.text.x = element_text(size = 5),
            axis.title = element_text(size = 8),
            legend.title = element_text(size = 8),
            legend.position = "bottom",
            legend.margin = margin(-5, -10, 0, 0),
            strip.text = element_text(size = 7, face = "bold"),
            legend.text = element_text(size = 7),
            strip.background = element_rect(fill = "lightgray", color = "black"),
            strip.placement = "outside",
            panel.spacing = unit(0, "pt"),
            panel.border = element_rect(color = "#7d7d7d", size = 0.5)))
  dev.off();
  return(paste0( name, " done!"))
}

get_gc<- function(vc, ws, margin){
  vc$gc <- 0
  for (position in unique(vc$pos)){
    # Get all gcs values of  windows having our position 
    gcs<-   ws$gc[ws$start + margin <= position & ws$end-margin >= position]
    if (length(gcs)<1){
      print("margin is too large not fitting in any window")
    }
    # assign the mean 
    vc$gc[vc$pos == position]  <- round(mean(gcs), 5)
  }
  return(vc)
}

plot_gc_variants<- function(v, window, step ){
  png(file=paste0(dirI, "eval_pos_w", window, "_s", step, "_.png"), res=180, width=2000, height=500);
  print(ggplot(v, aes(x = eval, y = gc, color=eval)) +
    geom_jitter(alpha=0.9, size=1, width = 0.35) +
    facet_grid(.~ploidy)+
    labs(title= paste0("Using window ", window, " and step ", step) ,x = "Evaluation Category", y = "GC Content") +
    guides(color="none")+
    scale_color_manual(values = palette3)+   
    theme_bw()  +
    theme(  axis.text.x = element_text(size = 7),
            axis.title = element_text(size = 9),
            strip.text = element_text(size = 7, face = "bold"),
            strip.background = element_rect(fill = "lightgray", color = "black"),
            strip.placement = "outside",
            panel.spacing = unit(0, "pt"),
            panel.border = element_rect(color = "#7d7d7d", size = 0.5)))
  dev.off();
  return(0)
}

plot_v_gc <-  function(ws, subwc, name, width, height, by){
    png(file=paste0(dirI,"GC_variants_", name ,".png"), res=300, width=width, height=height);
    print(ggplot(ws) +
      # geom_segment(aes(x = start, xend = end, y = gc, yend = gc )) +   # segments of gc
      geom_line(aes(x = middle, y = gc, color=as.factor(window))) +   # join them at the middle
      geom_point(data = subset(subwc,  eval != ""), aes(x = middle, y = gc,  color = as.factor(eval)),
                 size = 2, alpha=0.6) + 
      coord_cartesian(xlim = c(min(ws$start)+20,max(ws$end)-20))+
      scale_x_continuous(breaks = seq(min(subwc$middle), max(subwc$middle), by = as.numeric(by)))+
      scale_color_manual(values = palette3 )+
      labs(x = "Position", y = "GC Content", color ="") +
      guides(colour = guide_legend(nrow = 1))+
      facet_grid(.~region, scales="free_x",  space="free_x")+
      theme_classic()  +
      theme(  axis.text.x = element_text(size = 5),
              axis.title = element_text(size = 8),
              legend.title = element_text(size = 8),
              legend.position = "bottom",
              legend.margin = margin(-5, -10, 0, 0),
              strip.text = element_text(size = 7, face = "bold"),
              legend.text = element_text(size = 7),
              strip.background = element_rect(fill = "white", color = "black"),
              strip.placement = "outside",
              panel.spacing = unit(0, "pt"), 
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      panel.grid.major = element_line(colour = "#e6e6e6"),
      panel.grid.minor = element_line(colour = "#e6e6e6")))
    dev.off();
    return(paste0( name, " done!"))
}
#
## --------------------------------------  LOAD DATA ---------------------------------   ############
### Get variant calling dataframe
dirf<- paste0(dir, "03_VCsimulated/v2/VC_out/Sample_level_eval/4.3.0.0/VC/")
dirI<- paste0(dir, "03_VCsimulated/v2/IMAGES/01_PREJOIN/4.3.0.0/final/")
vc<- read.table(paste0(dirf, "VCeval_prejoin_p50_filtering_2sample_v4.tbl"), sep="\t")
# 
# vc <- read.table( paste0(dirf,"VCeval_postjoin.tbl"), sep="\t")
vc <- vc[vc$class =="callable",]
version<- "2sample_v4"
variants_t<- read.table(paste0(dirf, "Positions_per_evaluation",version ,".tbl"))
FPFNpos<- read.table(paste0(dirf, "FPFN_Positions_per_evaluation",version ,".tbl"))
# ### Get variants per evaluation table from 02.2.eval_postjoin_plotting.R
# variants_t<- read.table(paste0(dirf, "Positions_per_evaluation.tbl"))
# FPFNpos<- read.table(paste0(dirf, "FPFN_Positions_per_evaluation.tbl"))
### Table position evaluated samples
# sample_t<- read.table(paste0(dirf, "callable_eval_counts_sample_evaluation.tbl"))
# GET BED FILES #
w20s5<- as.data.frame(read.table(  paste0(dirBED,"hs1-rDNA_v1.0.chrR.20window_5sliding.GC_content.bedgraph"), sep="\t"))
w70s5<-  as.data.frame(read.table(  paste0(dirBED,"hs1-rDNA_v1.0.chrR.70window_5sliding.GC_content.bedgraph"), sep="\t"))
w50s1<-  as.data.frame(read.table(  paste0(dirBED,"hs1-rDNA_v1.0.chrR.50window_1sliding.GC_content.bedgraph"), sep="\t"))
w100s50<- as.data.frame(read.table(  paste0(dirBED,"hs1-rDNA_v1.0.chrR.100window_50sliding.GC_content.bedgraph"), sep="\t"))
## --------------------------------------  PREP DATA ---------------------------------   ############
w20s5<-prep_GC_bed(w20s5, 20, 5) # step 
w70s5<-prep_GC_bed(w70s5, 70, 5)# step 5
w50s1<-prep_GC_bed(w50s1, 50, 1)# step 1
w100s50<-prep_GC_bed(w100s50, 100, 50)# step 50

ws<- rbind(w20s5, w50s1, w70s5, w100s50)
ws$window<- as.factor(ws$window)
ws$region<- factor(ws$region, levels =c("IGS1", "5_ETS", "18S", "ITS1", "5.8S", "ITS2",  "28S", "3_ETS", "IGS2"))
# lets compare the ones that are correctly called TP with variant possitions always uncalled and some non-variant positions 
# are always being wrongly calling in all samples
v <- rbind(FPFNpos, variants_t[variants_t$eval =="TP",])
## ---------------------------------- PLOT GC content ---------------------------------   ############
# plot gc of all chrR
plot_gc(ws, "00_full_chrR", 5000, 1000, 5000)
# plot gc of all eDNA without IGS
plot_gc(ws[!ws$region %in% c("IGS1", "IGS2"),], "01_no_IGS", 5000, 1000, 1000)
# plot gc of  ITSs
plot_gc(ws[ws$region == "ITS1",], "ITS1", 4000, 1000,  500)
plot_gc(ws[ws$region == "ITS2",], "ITS2", 4000, 1000, 500)
# plot gc of  ITSs
plot_gc(ws[ws$region == "5_ETS",], "5_ETS",5000, 1000, 1000)
plot_gc(ws[ws$region == "3_ETS",], "3_ETS", 4000, 1000, 100)
# plot genes
plot_gc(ws[ws$region == "18S",], "18S", 4000, 1000, 500)
plot_gc(ws[ws$region == "5.8S",], "5.8S", 3000, 1000, 50)
plot_gc(ws[ws$region == "28S",], "28S", 5000, 1000, 1000)
## ---------------------------------- PLOT DIST VARIANTS GC ---------------------------------   ############
# Positions that where always FN or FP vs pos that got TPs
v<- get_gc(v, w20s5, 3)
plot_gc_variants(v, 20, 5) #plot

v<- get_gc(v, w70s5, 3)
plot_gc_variants(v, 70, 5)

v<- get_gc(v, w50s1, 3)
plot_gc_variants(v, 50, 1)

v<- get_gc(v, w100s50, 3)
plot_gc_variants(v, 100, 50)
## -------------------------------  DIST VARIANTS in GC plot---------------------------------   ############
# w20s5<- as.data.frame(read.table(  paste0(dirBED,"hs1-rDNA_v1.0.chrR.20window_5sliding.GC_content.bedgraph"), sep="\t"))
# w20s5<-prep_GC_bed(w20s5, 20, 5) # step 
ploidy<- "50"
# cut out IGS as we have them masked there will be no positions there
w20s5<- w20s5[w20s5$start>9338  & w20s5$end < 22670 ,]
w20s5$ploidy <-  ploidy
w20s5$eval <- ""
# fill w20s5 with the variants positions there 
for (position in v$pos[v$ploidy==ploidy]){
  eval<-   unique(v$eval[v$ploidy==ploidy & v$pos==position])
  wc_eval<- unique(w20s5$eval[w20s5$middle-3 <= position & w20s5$middle+3 >= position])
  # position fits only in one window reducing window from 20 to 14
  if (length(wc_eval)==1){ 
    # if its not empty write evaluation
    if (wc_eval==""){
      w20s5$eval[w20s5$middle-3 <= position & w20s5$middle+3 >= position] <-  eval
    }else{
      print(paste0("Superposing positions 1 ", position," ",  wc_eval))
      print(w20s5[w20s5$middle-3 <= position & w20s5$middle+3 >= position, ])
      row <-  w20s5[w20s5$middle-3 <= position & w20s5$middle+3 >= position, ]
      row$eval<- eval
      w20s5<- rbind(w20s5, row)
    }
  }else{ # in two of the windows
    if( "" %in% wc_eval){ # if one of the two is empty
       w20s5$eval[w20s5$middle-3 <= position & w20s5$middle+3 >= position & w20s5$eval == ""] <- eval
   }else{
      print(paste0("Superposing positions 2 ", position))
      print( w20s5[w20s5$middle-3 <= position & w20s5$middle+3 >= position, ])
      row <-  w20s5[w20s5$middle-3 <= position & w20s5$middle+3 >= position, ][1,]
      row$eval<- eval
      w20s5<- rbind(w20s5, row)
  }
 }
}
# cut out IGS as we have them masked there will be no positions there
ws<- ws[ws$start>9338  & ws$end < 22670 ,]
## ------------------------------- PLOT DIST VARIANTS in GC plot---------------------------------   ############
plot_v_gc(ws[!ws$region %in% c("IGS1", "IGS2"),], w20s5, "noIGS_chrR", 5000, 1000, 1000)

plot_v_gc(ws[ws$region == "ITS1",],  w20s5[w20s5$region == "ITS1",], "ITS1", 4000, 1000,  500)
plot_v_gc(ws[ws$region == "ITS2",],  w20s5[w20s5$region == "ITS2",], "ITS2", 4000, 1000, 500)
# plot gc of  ITSs
plot_v_gc(ws[ws$region == "5_ETS",], w20s5[w20s5$region == "5_ETS",], "5_ETS",5000, 1000, 1000)
plot_v_gc(ws[ws$region == "3_ETS",],  w20s5[w20s5$region == "3_ETS",], "3_ETS", 4000, 1000, 100)
# plot genes
plot_v_gc(ws[ws$region == "18S",],  w20s5[w20s5$region == "18S",], "18S", 4000, 1000, 500)
plot_v_gc(ws[ws$region == "5.8S",],  w20s5[w20s5$region == "5.8S",], "5.8S", 3000, 1000, 50)
plot_v_gc(ws[ws$region == "28S",], w20s5[w20s5$region == "28S",],  "28S", 5000, 1000, 1000)
