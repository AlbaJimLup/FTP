#!/bin/bash

#SBATCH --job-name=VC.gVCF
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/scripts/
#SBATCH --output=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/out/00_vc_v2_out/00VC.gVCF.%A_%a.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/err/00_vc_v2_err/00VC_.gVCF.%A_%a.err
#SBATCH --cpus-per-task=112
#SBATCH --qos=gp_bscls  #gp_debug
#SBATCH --account=bsc83
#SBATCH --array=1-100  
#SBATCH --time=24:00:00

## GATk version
#gatk=4.3.0.0
gatk=4.5.0.0

# load modules
module load  java-openjdk  gatk/${gatk} #MN5
#module load  java-jdk/8u131  gatk/${gatk} #MN5
#module load java/8u201 gatk/4.3.0.0 # MN
#module load java gatk/4.3.0.0 # Nord3v2


# sample (library)
sample_id="sample_"${SLURM_ARRAY_TASK_ID}
read_length=151
coverage=$(sed -n "${SLURM_ARRAY_TASK_ID}p" simulated_WGS.tsv | cut -f2)
version=$1
seq_type=$2
ploidy_value=$3
ploidy="ploidy_"${ploidy_value}
maxgen=$4

echo "Generating VC for " ${ploidy} " with max value" ${maxgen}


# directory for current filters


# paths
inpath=/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_RNAs/Raquel/rDNA-Mapping-Genomes/T2T/02_mapping/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${sample_id}/
outpath=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${gatk}/VC/${sample_id}/${ploidy}/
mkdir -p ${outpath}

outpath_filtering=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_${version}/simulated_WGS/${seq_type}/${gatk}/filtering/${sample_id}/${ploidy}/
mkdir -p ${outpath_filtering}
	
# infiles
reference=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/data/hs1-rDNA_${version}.fa
bam=${inpath}/${sample_id}.${read_length}PE_${coverage}x.sorted.chrR.f2F2308q20.wo_XA.bam

# Single-sample gVCF calling
#for region in ITS1 ITS2 18S 5.8S 28S 5_ETS 3_ETS;do
for region in 5.8S 28S 5_ETS 3_ETS 18S;do
#	gatk --java-options "-Xmx115g -Xms100g" HaplotypeCaller \
#	-I ${bam} \
#	-R ${reference} \
#	-O ${outpath}/${sample_id}.${ploidy}.${region}.g.vcf.gz \
#	-ERC BP_RESOLUTION \
#	--sample-ploidy ${ploidy_value} \
#	--intervals /gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/data/pre-rRNA_47S.included_intervals.${region}.bed \
#	--max-reads-per-alignment-start 0 \
#	--dont-use-soft-clipped-bases true \
#	--native-pair-hmm-threads 48 \
#	--max-alternate-alleles 6 \
#	--max-genotype-count ${maxgen}
	
#	echo "Done with  ${ploidy} ${region} VC"
	
	inputVCF=${outpath}/${sample_id}.${ploidy}.${region}.g.vcf.gz

	gatk --java-options "-Xmx3g -Xms3g" VariantFiltration \
	-R ${reference} \
	-V ${inputVCF} \
	-O ${outpath_filtering}/${sample_id}.${ploidy}.${region}_filtering.g.vcf.gz \
	--filter-expression "QUAL < 30.0"  \
	--filter-name "QUAL30"   \
	--filter-expression "QUAL < 400.0"  \
	--filter-name "QUAL400"   \
	--filter-expression "QUAL < 1000.0"  \
	--filter-name "QUAL1000" 
	--filter-expression "QUAL < 1"  \
	--filter-name "QUAL1"   
done




