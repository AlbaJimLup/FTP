#!/bin/bash

#ploidies="2 5 7 10 20 30 40 50 100" # 2 5 7 10 20 50 100"

ploidies="10"
maxgen=1024

for ploidy in $ploidies ;   do
	### 0. VC (100 jobs per region per ploidy) + filtering
#	maxgen=$(python3 -c "import math; num_alleles=6; ploidy=$ploidy; maxgen=math.factorial(num_alleles + ploidy - 1) // (math.factorial(num_alleles - 1) * math.factorial(ploidy)); print(maxgen)")
        echo ${maxgen}
    	sbatch 00_vc_v3.sh v1.0 rDNA_reads ${ploidy} ${maxgen}
	#Note on this maxgen we cannot do joint genotyping as it only works with their defaults

	### 2.  CREATE BD (100 jobs per region per ploidy)
#        for region in $regions ;do 
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
