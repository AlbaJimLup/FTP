#!/bin/bash


#SBATCH --job-name=05_08.merge_gVCFs
#SBATCH --chdir=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/scripts/
#SBATCH --output=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/out/05_merge_out/05_08.merge_gVCFs.%A_%a.out
#SBATCH --error=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/err/05_merge_err/05_08.merge_gVCFs.%A_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --qos=acc_debug
#SBATCH --account=bsc83
#SBATCH --time=00:30:00
#SBATCH --array=1-100

# load modules
#module load samtools/1.16.1 bwa/0.7.17 # MN
#module load java gatk/4.3.0.0 # nord3v2
module load bcftools

# sample (library)
#read_length=$1
#coverage=$2
version=$1
seq_type=$2
ploidy_value=$3
ploidy="ploidy_"${ploidy_value}
sample_id="sample_"${SLURM_ARRAY_TASK_ID}

# paths
inpath=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${sample_id}/${ploidy}/
outpath=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${sample_id}/${ploidy}/
file_basename="sample_"${SLURM_ARRAY_TASK_ID}

# merge VCF files
bcftools concat \
        ${outpath}/${file_basename}.${ploidy}.5_ETS.g.vcf.gz \
        ${outpath}/${file_basename}.${ploidy}.18S.g.vcf.gz \
        ${outpath}/${file_basename}.${ploidy}.ITS1.g.vcf.gz \
        ${outpath}/${file_basename}.${ploidy}.5.8S.g.vcf.gz \
        ${outpath}/${file_basename}.${ploidy}.ITS2.g.vcf.gz \
        ${outpath}/${file_basename}.${ploidy}.28S.g.vcf.gz \
        ${outpath}/${file_basename}.${ploidy}.3_ETS.g.vcf.gz \
        -o ${outpath}/${file_basename}.${ploidy}.g.vcf.gz
bcftools index -t ${outpath}/${file_basename}.${ploidy}.g.vcf.gz


