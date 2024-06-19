#!/bin/bash

#SBATCH --job-name=45S_clustal
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/05_rDNA
#SBATCH --output=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/05_rDNA/MSA/out/all_45S_rDNA_clustal.%A_%a.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/05_rDNA/MSA/err/all_45S_rDNA_clustal.%A_%a.err
#SBATCH --account=bsc83
#SBATCH --qos=gp_bscls
#SBATCH --cpus-per-task=48
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --array=25

module load clustal/1.2.4

i=${SLURM_ARRAY_TASK_ID}
outpath=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/05_rDNA/MSA/All_it${i}
mkdir -p ${outpath}

clustalo -i rDNAs_45S.fa  -o ${outpath}/All_45S.aln --distmat-out=${outpath}/All_45S.dist --guidetree-out=${outpath}/All_45S.nwk --percent-id --full  --iterations ${i}

