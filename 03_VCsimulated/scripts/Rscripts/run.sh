#!/bin/bash

ploidies="2 5 7 10 20 30 40 50 100"
ploidies="2 5 7 10 50 100"

for ploidy in $ploidies ;   do
	sbatch get_eval.sh 02.1.eval_postjoin.R "4.3.0.0" "joingen" ${ploidy} ""
	sbatch get_eval.sh 02.1.eval_postjoin.R "4.3.0.0" "filtering" ${ploidy} "v1" 
	sbatch get_eval.sh 02.1.eval_postjoin.R "4.3.0.0" "filtering" ${ploidy} "v2" 
	sbatch get_eval.sh 02.1.eval_postjoin.R "4.3.0.0" "filtering" ${ploidy} "v3" 

#	sbatch get_eval.sh 01.1.eval_prejoin.R "4.3.0.0" "VC" ${ploidy}
#	sbatch get_eval.sh 01.1.eval_prejoin.R "4.3.0.0" "filtering" ${ploidy} "v1"
#	sbatch get_eval.sh 01.1.eval_prejoin.R "4.3.0.0" "filtering" ${ploidy} "v2"
#	sbatch get_eval.sh 01.1.eval_prejoin.R "4.3.0.0" "filtering" ${ploidy} "v3"
#       sbatch get_eval.sh 01.1.eval_prejoin.R "4.3.0.0" "filtering" ${ploidy} "v4"

#	sbatch get_eval.sh "01.1.eval_prejoin_5sample.R" "4.3.0.0" "filtering" ${ploidy} "v2"
#	sbatch get_eval.sh "01.1.eval_prejoin_5sample.R" "4.3.0.0" "filtering" ${ploidy} "v3"
#	sbatch get_eval.sh "01.1.eval_prejoin_5sample.R" "4.3.0.0" "filtering" ${ploidy} "v4"
done
