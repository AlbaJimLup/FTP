#!/bin/bash

#SBATCH --job-name=05_07.merge_VCFs
#SBATCH --chdir=/gpfs/projects/bsc83/bMN4/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/scripts/
#SBATCH --output=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/out/05_merge_out/05_07.merge_VCFs.%A_%a.out
#SBATCH --error=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/err/05_merge_err/05_07.merge_VCFs.%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --qos=acc_debug
#SBATCH --account=bsc83
#SBATCH --time=00:15:00

# load modules
module load samtools/1.16.1 bwa/0.7.17 # MN
#module load java gatk/4.3.0.0 # nord3v2
module load bcftools

# sample (library)
#read_length=$1
#coverage=$2
version=$1
seq_type=$2
ploidy_value=$3
ploidy="ploidy_"${ploidy_value}

# paths
#inpath=/gpfs/projects/bsc83/bMN4/bsc83/Projects/ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/T2T/03_variants/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${ploidy}
#outpath=/gpfs/projects/bsc83/bMN4/bsc83/Projects/ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/T2T/03_variants/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${ploidy}/
inpath=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${ploidy}
outpath=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${ploidy}
file_basename=simulated_WGS

# merge VCF files
bcftools concat \
        ${outpath}/${file_basename}.${ploidy}.5_ETS.vcf.gz \
        ${outpath}/${file_basename}.${ploidy}.18S.vcf.gz \
        ${outpath}/${file_basename}.${ploidy}.ITS1.vcf.gz \
        ${outpath}/${file_basename}.${ploidy}.5.8S.vcf.gz \
        ${outpath}/${file_basename}.${ploidy}.ITS2.vcf.gz \
        ${outpath}/${file_basename}.${ploidy}.28S.vcf.gz \
        ${outpath}/${file_basename}.${ploidy}.3_ETS.vcf.gz \
        -o ${outpath}/${file_basename}.${ploidy}.vcf.gz
bcftools index -t ${outpath}/${file_basename}.${ploidy}.vcf.gz


