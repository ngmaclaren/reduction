#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=52
#SBATCH --mem=100000
#SBATCH --time=01:30:00

#SBATCH --mail-user=neil.g.maclaren@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=fMRI-fix
#SBATCH --output=./output/fMRI-fix.out
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc openmpi r

# Rscript analyze-fMRI.R
Rscript add-degree-preserving.R
