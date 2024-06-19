#!/bin/bash

#SBATCH --job-name=PlotVC_Analysis
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/scripts/Rscripts/
#SBATCH --output=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/out/Analysis/VC_gVCF.%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/err/Analysis/VC_gVCF.%A.err
#SBATCH --cpus-per-task=8
#SBATCH --qos=gp_bscls
#SBATCH --account=bsc83
#SBATCH --time=02:00:00


module load R/4.3.2

Rscript 01.1.eval_prejoin.R "4.3.0.0" "VC"
#Rscript 01.2.eval_prejoin_plotting.R "4.3.0.0" "VC"
#Rscript 01.3.eval_prejoin_analysis.R "4.3.0.0" "VC"
