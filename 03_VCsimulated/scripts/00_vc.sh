#!/bin/bash

#SBATCH --job-name=VC.gVCF
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/scripts/
#SBATCH --output=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/out/00_vc_v2_out/00VC.gVCF.%A_%a.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/err/00_vc_v2_err/00VC_.gVCF.%A_%a.err
#SBATCH --cpus-per-task=48
#SBATCH --qos=acc_bscls  #acc_debug
#SBATCH --account=bsc83
#SBATCH --array=1-100  
#SBATCH --time=00:50:00


## GATk version
gatk=4.3.0.0
#gatk=4.5.0.0

# load modules
module load  java-jdk/8u131  gatk/${gatk} #MN5
#module load java/8u201 gatk/4.3.0.0 # MN
#module load java gatk/4.3.0.0 # Nord3v2


# sample (library)
sample_id="sample_"${SLURM_ARRAY_TASK_ID}
read_length=151
coverage=$(sed -n "${SLURM_ARRAY_TASK_ID}p" simulated_WGS.tsv | cut -f2)
version=$1
seq_type=$2
region=$3
ploidy_value=$4
ploidy="ploidy_"${ploidy_value}
# directory for current filters


# paths
inpath=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/T2T/02_mapping/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${sample_id}/
outpath=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${gatk}/VC/${sample_id}/${ploidy}/
mkdir -p ${outpath}

# infiles
reference=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/data/hs1-rDNA_${version}.fa
bam=${inpath}/${sample_id}.${read_length}PE_${coverage}x.sorted.chrR.f2F2308q20.wo_XA.bam

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
--max-genotype-count 1024 
#--max-alternate-alleles 100 \
#--max-genotype-count 100000 
#--max-alternate-alleles 6 \ # ploidy 10
#--max-genotype-count 1024 \ # ploidy 10
#--minimum-mapping-quality 30 \
#--mapping-quality-threshold-for-genotyping 30 \
#--min-base-quality-score 15 \

echo "Done with  ${ploidy} ${region} VC"


