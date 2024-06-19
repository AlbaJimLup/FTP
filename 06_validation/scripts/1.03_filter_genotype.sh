#!/bin/bash

#SBATCH --job-name=FEvaluating
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/scripts/Rscripts/
#SBATCH --output=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/out/Analysis/Fevaluation.%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/err/Analysis/Fevaluation.%A.err
#SBATCH --cpus-per-task=8
#SBATCH --qos=gp_bscls
#SBATCH --account=bsc83
#SBATCH --time=01:00:00

file="03_filter_genotype.R"
sample=$1
filter_sample=$2
ploidy=$3

module load R/4.3.2

echo "TRYING" ${file} ${ploidy} ${filter_sample} ${sample}

Rscript ${file} ${ploidy} ${filter_sample} ${samples} 





