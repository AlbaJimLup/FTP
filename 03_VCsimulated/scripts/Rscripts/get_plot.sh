#!/bin/bash

#SBATCH --job-name=Plotting
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/scripts/Rscripts/
#SBATCH --output=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/out/Analysis/plotting.%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/err/Analysis/plotting.%A.err
#SBATCH --cpus-per-task=48
#SBATCH --qos=gp_bscls
#SBATCH --account=bsc83
#SBATCH --time=01:00:00

file=$1
version=$2
step=$3
filter=$4

module load R/4.3.2

echo "TRYING" ${file} ${version} ${step} ${filter}

Rscript ${file} ${version} ${step} ${filter} 




