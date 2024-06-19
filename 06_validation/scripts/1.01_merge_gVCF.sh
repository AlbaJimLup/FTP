#!/bin/bash

#SBATCH --job-name=01.merge_gVCFs
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/08_validation/
#SBATCH --output=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/08_validation/out/01.merge_gVCFs.%A_%a.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/08_validation/err/01.merge_gVCFs.%A_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --qos=acc_debug 
#SBATCH --account=bsc83
#SBATCH --time=00:30:00

# load modules
#module load samtools/1.16.1 bwa/0.7.17 # MN
#module load java gatk/4.3.0.0 # nord3v2
module load bcftools

# sample (library)
version=$1
seq_type=$2
sample_id=$3

echo ${sample_id}

ploidy_value=$4
ploidy="ploidy_"${ploidy_value}

# paths
inpath=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/08_validation/VC_out/${sample_id}
outpath=${inpath}

file_basename="${sample_id}.${ploidy}"

# merge VCF files
bcftools concat \
        ${outpath}/${file_basename}.5_ETS.g.vcf.gz \
        ${outpath}/${file_basename}.18S.g.vcf.gz \
        ${outpath}/${file_basename}.ITS1.g.vcf.gz \
        ${outpath}/${file_basename}.5.8S.g.vcf.gz \
        ${outpath}/${file_basename}.ITS2.g.vcf.gz \
        ${outpath}/${file_basename}.28S.g.vcf.gz \
        ${outpath}/${file_basename}.3_ETS.g.vcf.gz \
        -o ${outpath}/${file_basename}.g.vcf.gz
bcftools index -t ${outpath}/${file_basename}.g.vcf.gz
