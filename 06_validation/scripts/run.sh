#!/bin/bash

regions="5_ETS 18S ITS1 5.8S ITS2 28S 3_ETS"
samples="SRR1997411 SRR3189741 SRR3189742 SRR3189743"

ploidy=150

for sample in $samples ;   do
	### 0. VC
#        for region in $regions ;do # create DB
#                sbatch 00_vc.sh v1.0 rDNA_reads ${region} ${sample}  ${ploidy};done
	
	### 1. MERGE gVCF FILES	
#        sbatch 1.01_merge_gVCF.sh v1.0 rDNA_reads ${sample} ${ploidy};

	### 2. Filter
#        sbatch 1.02_filter.sh v1.0 rDNA_reads ${sample} ${ploidy};

	### 3. Manual "Join-genotype" and filtering of uncalled positions
	sbatch 1.03_filter_genotype.sh ${ploidy} "4" ${samples};

done


### VERSION 2 OF THE PIPELINE ####	
### 2.  CREATE BD
#for region in $regions ;do # create DB
#	sbatch 2.01_DB.sh v1.0 rDNA_reads ${region} ${ploidy};done

        ### 3.  JOIN-GENOTYPING
#       for region in $regions ;do # perform join-genotyping
#              sbatch 03_joingen.sh v1.0 rDNA_reads ${region} ${ploidy};done

        ### 4.  FILTERING
#       #GATK'S VARIANT FILTERING DOES NOT WORK #
#       for region in $regions ;do # perform join-genotyping
#             sbatch 04_VariantFiltration.sh v1.0 rDNA_reads ${region} ${ploidy};done
#       #We will use manual filter:
#        outpath=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_v1.0/simulated_WGS/rDNA_reads/4.3.0.0/filtering/ploidy_${ploidy}/
#        mkdir -p ${outpath}
#      sbatch 04_Filtering.sh 04_Filtering.R "4.3.0.0" ${ploidy} "v1"
#       sbatch 04_Filtering.sh 04_Filtering.R "4.3.0.0" ${ploidy} "v2"
#       sbatch 04_Filtering.sh 04_Filtering.R "4.3.0.0" ${ploidy} "v3"

