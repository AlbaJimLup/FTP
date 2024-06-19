#!/bin/bash

#SBATCH --job-name=02.create_db
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/scripts/
#SBATCH --output=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/out/02_DB_out/02.create_db.%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/err/02_DB_err/02.create_db.%A.err
#SBATCH --cpus-per-task=16
#SBATCH --qos=gp_bscls # ploidy 100 regions 28S and 5_ETS  for others #SBATCH --qos=acc_bscls 
#SBATCH --account=bsc83
#SBATCH --time=06:00:00


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
#inpath=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_${version}/simulated_WGS/${seq_type}/
#outpath=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${ploidy}
inpath=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${gatk}/VC/
outpath=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${gatk}/joingen/${ploidy}
mkdir -p ${outpath}

# infiles
reference=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/data/hs1-rDNA_${version}.fa

# GenomicsDBImport https://gatk.broadinstitute.org/hc/en-us/articles/5358869876891-GenomicsDBImport
# Import single-sample GVCFs into GenomicsDB before joint genotyping
# creating sample map file
for i in `find ${inpath}/ -name *.${ploidy}.${region}.g.vcf.gz`;do file=`basename $i`; sample="$(cut -d'.' -f1 <<<"$file")" ;echo -e ${sample}"\t"${i} >> ${outpath}/simulated_WGS.${ploidy}.${region}.g.vcf.map;done
# Provide sample map with sample IDs and paths to each g.vcf to consolidate all samples before performing joint genotyping
# NOTE: for big ploidy values (such as 100), you might have to increase the default buffer size to solve the issue. 
gatk --java-options "-Xmx16g -Xms14g"  GenomicsDBImport \
--sample-name-map ${outpath}/simulated_WGS.${ploidy}.${region}.g.vcf.map \
--genomicsdb-workspace-path ${outpath}/simulated_WGS.${ploidy}.${region}.database \
--intervals /gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/data/pre-rRNA_47S.included_intervals.${region}.bed \
--reader-threads 4 \
--batch-size 50 \
--genomicsdb-vcf-buffer-size 16384000 

