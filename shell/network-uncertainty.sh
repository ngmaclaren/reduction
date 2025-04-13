#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=52
#SBATCH --mem=100000
#SBATCH --job-name=nu-email
#SBATCH --output=./output/nu-email.out
#SBATCH --mail-user=neil.g.maclaren@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc openmpi r

Rscript ../analysis/network-uncertainty.R
