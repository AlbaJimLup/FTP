#!/bin/bash

#SBATCH --job-name=03.genotype_gVCF
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/scripts/
#SBATCH --output=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/out/03_joingen_out/03.genotype_gVCF.%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/err/03_joingen_err/03.genotype_gVCF.%A.err
#SBATCH --cpus-per-task=32
#SBATCH --qos=gp_bscls
#SBATCH --account=bsc83
#SBATCH --time=01:00:00

## GATk version
gatk=4.3.0.0
#gatk=4.5.0.0

# load modules
module load  java-jdk/8u131  gatk/${gatk} #MN5
#module load java/8u201 gatk/4.3.0.0 # MN
#module load java gatk/4.3.0.0 # nord3v2

# sample (library)
version=$1
seq_type=$2
region=$3
ploidy_value=$4
ploidy="ploidy_"${ploidy_value}

# paths
#inpath=/gpfs/projects/bsc83/PMN4/bsc83/rojects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_${version}/simulated_WGS/${seq_type}/
#outpath=/gpfs/projects/bsc83/PMN4/bsc83/rojects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${ploidy}
inpath=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${gatk}/VC/
outpath=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${gatk}/joingen/${ploidy}

#mkdir -p ${outpath}

# infiles
reference=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/data/hs1-rDNA_${version}.fa

echo "Trying with ploidy  ${ploidy}  and region  ${region}"

# Perform joint genotyping using all the samples from the previously created DB
# Check if ploidy equals to 100 (special settings) 
database=${outpath}/simulated_WGS.${ploidy}.${region}.database

gatk --java-options "-Xmx16g -Xms14g" GenotypeGVCFs \
-R ${reference} \
-V gendb://${database} \
-O ${outpath}/simulated_WGS.${ploidy}.${region}.vcf.gz \
--sample-ploidy ${ploidy_value} \
--max-genotype-count 2999999999
		
echo "Done!" 
