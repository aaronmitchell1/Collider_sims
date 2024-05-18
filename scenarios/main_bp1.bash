#!/bin/bash
#SBATCH --job-name=collider_sims_corr
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --time=47:59:59
#SBATCH --mail-user=aaron.mitchell@bristol.ac.uk
#SBATCH --nodes=4
#SBATCH --account=sscm013902
#SBATCH --mem=16G
module load apps/jags/4.3.0
module load lang/r/4.3.0-gcc
Rscript sim_HPC_increased_corr_cont_prog.R