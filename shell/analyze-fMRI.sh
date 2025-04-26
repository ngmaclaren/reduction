#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1 # 52
#SBATCH --mem=10000 # 0
#SBATCH --time=06:00:00

#SBATCH --mail-user=neil.g.maclaren@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=fMRI
#SBATCH --output=./output/fMRI.out
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc openmpi r

Rscript analyze-fMRI.R
