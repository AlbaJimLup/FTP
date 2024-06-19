#!/bin/bash

#SBATCH --job-name=04.Filtering
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/scripts/
#SBATCH --output=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/out/04_Rfilter/Filtering.%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/err/04_Rfilter/Filtering.%A.err
#SBATCH --cpus-per-task=8
#SBATCH --qos=gp_bscls
#SBATCH --account=bsc83
#SBATCH --time=01:30:00

file=$1
version=$2
ploidy=$3
v=$4

module load R/4.3.2

Rscript ${file} ${version} ${ploidy} ${v} 


