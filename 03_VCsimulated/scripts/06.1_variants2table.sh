#!/bin/bash

#SBATCH --job-name=05_09.variants2table
#SBATCH --chdir=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/scripts/
#SBATCH --output=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/out/06_variants2table_out/05_09.variants2table.%A_%a.out
#SBATCH --error=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/err/06_variants2table_err/05_09.variants2table.%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --qos=acc_debug
#SBATCH --account=bsc83
#SBATCH --time=01:00:00

# load modules
module load samtools/1.19.2
#module load samtools/1.16.1 bwa/0.7.17 # MN
#module load java gatk/4.3.0.0 # nord3v2
#module load java gatk/4.3.0.0 # nord3v2

# sample (library)
#read_length=$1
#coverage=$2
version=$1
seq_type=$2
ploidy_value=$3
ploidy="ploidy_"${ploidy_value}


# paths
#inpath=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/T2T/03_variants/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${ploidy}
#outpath=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/T2T/03_variants/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${ploidy}/
inpath=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${ploidy}
outpath=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${ploidy}
file_basename=simulated_WGS.${ploidy}

# infiles
reference=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Alba/data/hs1-rDNA_${version}.fa

# variants 2 table
gatk --java-options "-Xmx4g -Xms2g" VariantsToTable \
        -V ${outpath}/${file_basename}.vcf.gz \
        -F CHROM -F POS -F TYPE -F REF -F ALT -F QUAL -F QD -F FS -F SOR -F MQ -F MQRankSum -F ReadPosRankSum \
        -GF GT -GF AD -GF DP -GF GQ \
        -O ${outpath}/${file_basename}.variants.GQ.tab
