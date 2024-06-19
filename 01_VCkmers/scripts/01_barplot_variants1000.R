############----------------------------- DATA prep ----------------------------#############
# Directory
dir <- "~/Desktop/ribosomal_RNAs/Alba/01_VCkmers/"
dirI<-  paste0(dir, "IMAGES/")
# libraries 
library(readr)
library(forcats)
library(ggplot2)
library(dplyr)
# Get table
s <- read_csv(paste0(dir, "CHM13_variants.msa_vs_diploid_GATK - 1000x.csv"), show_col_types = FALSE)
s <-  as.data.frame(s)
##########----------------------------- DATA prep ----------------------------#############

# Setting column names
colnames(s) <-  c(s[1,1:10], "class2","region2" ,"evaluation", "annotations", s[1,15])
s <- s[-1,]
# Removing two last rows
s <- s[-nrow(s), ]
s <- s[-nrow(s), ]
#Change name "extra variant"
s$evaluation[s$evaluation =="extra variant" ] <- "FP"

# Diverging alleles postions table:
diverging <- data.frame(  position = character(),
                          GATKref = character(),
                          GATKalt = character(),
                          MSAref = character(),
                          MSAalt = character())
for (row in 1:nrow(s)){
  if(is.na(s$class[row])){
    s$class[row] <-  s$class2[row]
  }
  if(s$evaluation[row] != "same variant"){
    diverging <- rbind(diverging, 
                       data.frame(  position = ifelse(is.na(s$POS[row]),  s$Pos[row], s$POS[row]),
                                    GATKref = s$REF[row],
                                    GATKalt = s$ALT[row],
                                    MSAref = s$Ref[row],
                                    MSAalt = s$alt[row]))} 
}
# Export alleles postions table:
write.table(diverging, paste0(dir, "x1000diverging_positions_table.tbl"), sep = "\t", quote = FALSE, row.names = FALSE)
## We could remove class 2 now
# s <- s[,-11]

##########------------------------ CLASSIFY CALLABLE/BLACKLISTED  -----------------------#############

bed <- read.table(paste0(dir, "47S_pre-rRNA.blacklisted_repetitive_seq.bed"))

get_black<-  function(df){
  
  for (pos in unique(df$POS)){
    
    eval_black <-  unique(ifelse(pos > bed$V2 & pos <= bed$V3, "YES", "NO"))
    
    if( "YES" %in% eval_black ){
      df$class[df$POS == pos] <- "blacklisted"
    }
    else{
      df$class[df$POS == pos] <- "callable"
    }
  }
  return(df)
}
#

s<- get_black(s)

##########---------------------------- Plotting variants  -----------------------#############

# Plotting variants evaluation:
u_eval <-  unique(s$evaluation)

counts_ev <-  data.frame("evaluation" =  c(u_eval, u_eval),
                         "counts" = rep(0, length(u_eval)*2),
                         "freqs" = rep(0, length(u_eval)*2),
                         "class" = c(rep("callable", length(u_eval)), 
                                     rep("blacklisted", length(u_eval))))

callable <- sum((s$class == "callable"))
blacklisted <- sum((s$class == "blacklisted"))

for (row in 1:length(u_eval)){
  counts_ev$counts[row] = sum(s$evaluation == u_eval[row] & (s$class == "callable"))
  counts_ev$freqs[row] = sum(s$evaluation == u_eval[row] & (s$class == "callable"))/callable
  counts_ev$counts[row+5] = sum(s$evaluation == u_eval[row] & (s$class == "blacklisted"))
  counts_ev$freqs[row+5] = sum(s$evaluation == u_eval[row] & (s$class == "blacklisted"))/blacklisted
}


palette <- setNames(c( "#5bd18e",  "#ffc60e", "#ecd189",  "#ff6072","#bd0752"),
                    c("same variant","missing variant", "missing allele",
                      "different alleles",  "FP" ))


counts_ev <- counts_ev %>%
  filter(counts != 0) 
# Export counts table:
write.table(counts_ev, paste0(dir,  "x1000stats_evaluationGATKvsMSA.tbl"), sep = ",", quote = FALSE, row.names = FALSE)

png(file=paste0(dirI,"x1000GATKvsMSAalleles.png"), res=300, width=2300, height=1300);
ggplot(counts_ev, aes(x = fct_reorder(class, freqs), y =freqs)) +
  geom_bar(stat = "identity", aes(fill = evaluation))+
  geom_text(aes(label = c("39","1","1","56","45","8","17","11")), 
            position = position_stack(vjust = 0.3), size = 3.5,fontface = "bold") +
  theme_minimal()+
  labs(title ="GATK compared to MSA Allele Variants",
       subtitle = "Kmers of depth 1000x",
       y="Frequency")+
  scale_fill_manual(values=palette, name = "")+
  theme(axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        legend.position = "right")+
  labs(fill ="")
dev.off()

