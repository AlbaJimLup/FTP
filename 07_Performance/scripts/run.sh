#!/bin/bash

sbatch get_tables.sh 02.1_Eval_Sim_errors_VC.R
sbatch get_tables.sh 02.2_Eval_Sim_errors_VC_filtered.R
sbatch get_tables.sh 02.3_Eval_Sim_errors_VC_filtered_sample.R
sbatch get_tables.sh 03.1_Eval_Sim_errors_joingen.R
sbatch get_tables.sh 03.2_Eval_Sim_errors_joingen_filtered.R


