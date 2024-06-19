#!/bin/bash

#SBATCH --job-name=05_10.gVCF2table
#SBATCH --chdir=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/T2T/scripts/
#SBATCH --output=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/T2T/03_variants/out/06_simulated_WGS/05_10.gVCF2table.%A_%a.out
#SBATCH --error=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/T2T/03_variants/err/06_simulated_WGS/05_10.gVCF2table.%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --qos=acc_debug
#SBATCH --account=bsc83
#SBATCH --time=01:00:00
#SBATCH --array=1-100

# load modules
#module load samtools/1.16.1 bwa/0.7.17 # MN
#module load java gatk/4.3.0.0 # nord3v2
module load java gatk/4.3.0.0 # nord3v2

# sample (library)
sample_id="sample_"${SLURM_ARRAY_TASK_ID}
#read_length=$1
#coverage=$2
version=$1
seq_type=$2
ploidy_value=$3
ploidy="ploidy_"${ploidy_value}

# paths
inpath=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/T2T/03_variants/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${sample_id}/${ploidy}
outpath=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/T2T/03_variants/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${sample_id}/${ploidy}/
file_basename=${sample_id}.${ploidy}

# infiles
reference=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/data/Human_hs1-rDNA_genome_${version}/hs1-rDNA_${version}.fa

# variants 2 table
gatk --java-options "-Xmx4g -Xms2g" VariantsToTable \
        -V ${outpath}/${file_basename}.g.vcf.gz \
        -F CHROM -F POS -F TYPE \
	-GF GT -GF AD -GF DP -GF GQ \
        -O ${outpath}/${file_basename}.g.variants_depth.tab


