# libraries ----
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(cowplot)
library(reshape2)
library(ggrepel)
################# GET AND PREP FILES #################
dir="~/Desktop/ribosomal_RNAs/Alba/"
dirI=paste0(dir, "05_rDNA/IMAGES/")

d<- read.table(paste0(dir, "data/hg38_rDNA-like_regions.bed"))
names(d)<- c("chr", "start", "end")

################# GET LENGTHS  #################

d$length<- d$end-d$start
# sort by length
d <- d %>%  arrange(desc(length))
d$seq<- 1:nrow(d)
# get mean
mean_length <- mean(d$length)

# 
# sum(d$length)
# [1] 880043
#################  PLOT LONG FORMAT  ##############################
palette<- c("#75bb6f","#bd0752")

png(file=paste0(dirI,"rDNA-like_regions_lengths_barplot.png"), res=300, width=1000, height=1000);
ggplot(d, aes(x = 1:473, y = length)) +
  geom_bar(stat = "identity", fill = "#8bb5a6", color = "#8bb5a6") +
  coord_cartesian(xlim = c(1, nrow(d)), ylim = c(10, max(d$length))) +
  scale_x_continuous(breaks = seq(0, nrow(d)+2000, by = 150), expand = c(0.005, 0)) +
  scale_y_continuous(breaks = seq(0,40000000, by = 50000), expand = c(0.005, 0)) +
  geom_hline(yintercept = mean_length, color = "black", size = 0.5, linetype = "dashed") +
  labs(x = "rDNA-like regions", y = "Length" ) +
  annotate("text", x = 300, y = 9500, label = paste(round(mean_length, 2), "bp"), color = "black", size = 3) +
  theme_minimal() +
  theme( legend.margin = margin(0, 0, 0, 0),   axis.title = element_text(size=8), axis.text = element_text(size = 7)  )
dev.off();

# d<- d[2:nrow(d),]

png(file=paste0(dirI,"rDNA-like_regions_lengths_boxplot.png"), res=300, width=500, height=1200);
ggplot(d, aes(x = 1, y = length)) +
  geom_jitter(color = "#6bb5a6", width = 0.2,  alpha=0.4, size=1.5) +
  scale_y_continuous(breaks = seq(0,3411000, by = 20000), expand = c(0.005, 0)) +
  geom_boxplot(outlier.alpha = 0) +
  labs(x = "", y = "Length rDNA-like regions" ) +
  theme_minimal() +
  theme( legend.margin = margin(0, 0, 0, -11), axis.text.x = element_blank(),  axis.title = element_text(size=8), axis.text = element_text(size = 7)  )
dev.off();


png(file=paste0(dirI,"rDNA-like_regions_lengths_boxplot_noout.png"), res=300, width=500, height=1200);
ggplot(d, aes(x = 1, y = length)) +
  coord_cartesian( ylim = c(10, 2045)) +
  scale_y_continuous(breaks = seq(0,3000, by = 200), expand = c(0.005, 0)) +
  geom_jitter(color = "#6bb5a6", width = 0.3,  alpha=0.4, size=1.5) +
  geom_boxplot(outlier.alpha = 0, alpha=0) +
  labs(x = "", y = "Length rDNA-like regions" ) +
  theme_minimal() +
  theme( legend.margin = margin(0,-10,0 , -11), axis.text.x = element_blank(),  axis.title = element_text(size=8), axis.text = element_text(size = 7)  )
dev.off();


