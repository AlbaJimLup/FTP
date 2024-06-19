#!/bin/bash

#SBATCH --job-name=VC.gVCF
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/08_validation/
#SBATCH --output=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/08_validation/out/00VC.gVCF.%A_%a.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/08_validation/err/00VC_.gVCF.%A_%a.err
#SBATCH --cpus-per-task=48
#SBATCH --qos=acc_bscls  #acc_debug
#SBATCH --account=bsc83
#SBATCH --time=00:50:00

## GATk version
gatk=4.3.0.0
#gatk=4.5.0.0

# load modules
module load  java-jdk/8u131  gatk/${gatk} #MN5
#module load java/8u201 gatk/4.3.0.0 # MN
#module load java gatk/4.3.0.0 # Nord3v2

# sample (library)
sample_id=$4
version=$1
seq_type=$2
region=$3

#WE WILL USE PLOIDY 50 100 and 7#
ploidy_value=$5
ploidy="ploidy_"${ploidy_value}

# paths
inpath=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Jose/04_Pipeline/results/
outpath=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/08_validation/VC_out/${sample_id}
mkdir -p ${outpath}

# infiles
reference=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/data/hs1-rDNA_${version}.fa

To change
bam=${inpath}/${sample_id}.sorted.chrR.f2F2308q20.wo_XA.bam


# Single-sample gVCF calling
gatk --java-options "-Xmx115g -Xms100g" HaplotypeCaller \
-I ${bam} \
-R ${reference} \
-O ${outpath}/${sample_id}.${ploidy}.${region}.g.vcf.gz \
-ERC BP_RESOLUTION \
--sample-ploidy ${ploidy_value} \
--intervals /gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/data/pre-rRNA_47S.included_intervals.${region}.bed \
--max-reads-per-alignment-start 0 \
--dont-use-soft-clipped-bases true \
--native-pair-hmm-threads 16 \
--max-alternate-alleles 6 \
--max-genotype-count 1024 \
#--max-alternate-alleles 6 \ # ploidy 10
#--max-genotype-count 1024 \ # ploidy 10
#--minimum-mapping-quality 30 \
#--mapping-quality-threshold-for-genotyping 30 \
#--min-base-quality-score 15 \

echo "Done with  ${ploidy} ${region} VC"


