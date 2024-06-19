#!/bin/bash

#SBATCH --job-name=Tables
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/07_Performance/
#SBATCH --output=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/07_Performance/out/tables.%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/07_Performance/err/tables.%A.err
#SBATCH --cpus-per-task=48
#SBATCH --qos=gp_bscls
#SBATCH --account=bsc83
#SBATCH --time=01:00:00

file=$1

module load R/4.3.2

echo "TRYING" ${file} 

Rscript ${file} 

