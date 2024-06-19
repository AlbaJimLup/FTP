#!/bin/bash

#SBATCH --job-name=merge
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/scripts/Rscripts/
#SBATCH --output=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/out/Analysis/merge.%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/err/Analysis/merge.%A.err
#SBATCH --cpus-per-task=8
#SBATCH --qos=gp_bscls
#SBATCH --account=bsc83
#SBATCH --time=03:00:00

module load R/4.3.2

#Rscript 01.2.eval_prejoin_merge.R "4.3.0.0" "VC" "2 5 7 10 20 30 40 50 100"

#Rscript 01.2.eval_prejoin_merge.R "4.3.0.0" "filtering" "v1"  "2 5 7 10 20 30 40 50 100"  
#Rscript 01.2.eval_prejoin_merge.R "4.3.0.0" "filtering" "v2" "2 5 7 10 20 30 40 50 100"  
#Rscript 01.2.eval_prejoin_merge.R "4.3.0.0" "filtering" "v3" "2 5 7 10 20 30 40 50 100"  
#Rscript 01.2.eval_prejoin_merge.R "4.3.0.0" "filtering" "v4" "2 5 7 10 20 30 40 50 100"  

#Rscript "01.2.eval_prejoin_merge.R" "4.3.0.0" "filtering" "2sample_v2" "2 5 7 10 20 30 40 50 100"
#Rscript "01.2.eval_prejoin_merge.R" "4.3.0.0" "filtering" "2sample_v3" "2 5 7 10 20 30 40 50 100"
#Rscript "01.2.eval_prejoin_merge.R" "4.3.0.0" "filtering" "2sample_v4" "2 5 7 10 20 30 40 50 100"

#Rscript 01.2.eval_prejoin_merge.R "4.3.0.0" "filtering" "5sample_v2" "2 5 7 10 20 30 40 50 100" 
#Rscript 01.2.eval_prejoin_merge.R "4.3.0.0" "filtering" "5sample_v3" "2 5 7 10 20 30 40 50 100" 
#Rscript 01.2.eval_prejoin_merge.R "4.3.0.0" "filtering" "5sample_v4" "2 5 7 10 20 30 40 50 100"  

#Rscript 01.2.eval_prejoin_merge.R "4.3.0.0" "filtering" "10sample_v2" "2 5 7 10 20 30 40 50 100"
#Rscript 01.2.eval_prejoin_merge.R "4.3.0.0" "filtering" "10sample_v3" "2 5 7 10 20 30 40 50 100"
#Rscript 01.2.eval_prejoin_merge.R "4.3.0.0" "filtering" "10sample_v4" "2 5 7 10 20 30 40 50 100"  

#Rscript 01.2.eval_prejoin_merge.R "4.3.0.0" "filtering" "15sample_v2" "2 5 7 10 20 30 40 50 100"
#Rscript 01.2.eval_prejoin_merge.R "4.3.0.0" "filtering" "15sample_v3" "2 5 7 10 20 30 40 50 100"
#Rscript 01.2.eval_prejoin_merge.R "4.3.0.0" "filtering" "15sample_v4" "2 5 7 10 20 30 40 50 100"


Rscript 02.1.eval_postjoin_merge.R "4.3.0.0" "joingen" "2 5 7 10 50 100"
#Rscript 02.1.eval_postjoin_merge.R "4.3.0.0" "filtering" "v1" "2 5 7 10 50 100"
#Rscript 02.1.eval_postjoin_merge.R "4.3.0.0" "filtering" "v2" "2 5 7 10 50 100"
#Rscript 02.1.eval_postjoin_merge.R "4.3.0.0" "filtering" "v3" "2 5 7 10 50 100"



