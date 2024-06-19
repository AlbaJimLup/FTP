#!/bin/bash

#SBATCH --job-name=02filter.gVCF
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/08_validation/
#SBATCH --output=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/08_validation/out/02filter.%A_%a.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/08_validation/err/02filter.gVCF.%A_%a.err
#SBATCH --cpus-per-task=48
#SBATCH --qos=acc_bscls  #acc_debug
#SBATCH --account=bsc83
#SBATCH --time=01:15:00

## GATk version
gatk=4.3.0.0
#gatk=4.5.0.0

# load modules
module load  java-jdk/8u131  gatk/${gatk} #MN5
#module load java/8u201 gatk/4.3.0.0 # MN
#module load java gatk/4.3.0.0 # Nord3v2


# sample (library)
version=$1
seq_type=$2
sample_id=$3

ploidy_value=$4
ploidy="ploidy_"${ploidy_value}
# directory for current filters

# paths
inpath=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/08_validation/VC_out/${sample_id}
outpath=${inpath}
mkdir -p ${outpath}

reference=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/data/hs1-rDNA_${version}.fa
inputVCF=${inpath}/${sample_id}.${ploidy}.g.vcf.gz

# Here we annotate the VCF files' FILTER column with the filters applied based on 
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
# we should distiguish between SNPs and INDELS we'll first try SNPs recommendations

gatk --java-options "-Xmx3g -Xms3g" VariantFiltration \
-R ${reference} \
-V ${inputVCF} \
-O ${outpath}/${sample_id}.${ploidy}_filtering.g.vcf.gz \
--filter-expression "QUAL < 30.0"  \
--filter-name "QUAL30" \
--filter-expression "QUAL < 400.0"  \
--filter-name "QUAL400" \
--filter-expression "QUAL < 1000.0"  \
--filter-name "QUAL1000" 


## Probably we'll have to remove variables we don't have:  QD SOR FS MQ  
#--filter-expression "QD < 2.0"  \
#--filter-name "QD2" \
#--filter-expression "QUAL < 30.0"  \
#--filter-name "QUAL30"   \   
#--filter-expression "SOR > 3.0"  \
#--filter-name "SOR3" \
#--filter-expression "FS > 60.0"  \
#--filter-name "FS60" \
#--filter-expression "MQ < 40.0"  \
#--filter-name "MQ40" \
#--filter-expression "MQRankSum < -12.5"  \
#--filter-name "MQRankSum-12.5" \
#--filter-expression "ReadPosRankSum < -8.0" \
#--filter-name "ReadPosRankSum-8" \
## FOR INDELS
#--filter-expression "FS > 200.0"  \
#--filter-name "FS200" \
#--filter-expression "ReadPosRankSum < -20.0" \
#--filter-name "ReadPosRankSum-20" \
# Add this based on what we see in https://docs.google.com/presentation/d/1VMbayj-oQh6zih3BFuiWD7u6MiBZrtgnF6NfkJoI3v8/edit#slide=id.g2dbd7327e25_0_27
# we need to check if QD<2 removes them or we need to add them or we need to add it 
#--filter-expression "BaseQRankSum < -6" \
#--filter-name "BaseQRankSum-6"

echo "Done with ${sample_id}  ${ploidy} filtering with filtering"

