#!/bin/bash

regions="ITS1" # "5_ETS 18S ITS1 5.8S ITS2 28S 3_ETS"
ploidies="50 100" #"2 5 7 10 20 30 40 50 100" # 2 5 7 10 20 50 100"

for ploidy in $ploidies ;   do
	### 0. VC
        for region in $regions ;do # create DB
                sbatch 00_vc.sh v1.0 rDNA_reads ${region} ${ploidy};done
	
	### 1. Filter
#        for region in $regions ;do # create DB
#                sbatch 01_filter.sh v1.0 rDNA_reads ${region} ${ploidy};done

	### 2.  CREATE BD
#        for region in $regions ;do # create DB
#                sbatch 02_DB.sh v1.0 rDNA_reads ${region} ${ploidy};done

        ### 3.  JOIN-GENOTYPING
#       for region in $regions ;do # perform join-genotyping
#              sbatch 03_joingen.sh v1.0 rDNA_reads ${region} ${ploidy};done

	### 4.  FILTERING
#	#GATK'S VARIANT FILTERING DOES NOT WORK #
#       for region in $regions ;do # perform join-genotyping
#             sbatch 04_VariantFiltration.sh v1.0 rDNA_reads ${region} ${ploidy};done
#	#We will use manual filter:
#        outpath=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/03_VCsimulated/v2/VC_out/hs1-rDNA_v1.0/simulated_WGS/rDNA_reads/4.3.0.0/filtering/ploidy_${ploidy}/
#        mkdir -p ${outpath}
#      sbatch 04_Filtering.sh 04_Filtering.R "4.3.0.0" ${ploidy} "v1"
#	sbatch 04_Filtering.sh 04_Filtering.R "4.3.0.0" ${ploidy} "v2"
#	sbatch 04_Filtering.sh 04_Filtering.R "4.3.0.0" ${ploidy} "v3"
done
